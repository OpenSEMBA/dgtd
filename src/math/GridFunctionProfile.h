#pragma once

#include "math/Function.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "general/kdtree.hpp"

namespace maxwell {

/// Pointwise evaluation of a scalar H1 GridFunction stored on an external
/// MFEM mesh (plasma_kit archive: ``mesh.mesh`` + ``omega_p.gf``/``nu_c.gf``).
///
/// The Lorentz machinery of both evolution operators evaluates spatial rate
/// profiles at element quadrature points and pre-multiplies the resulting
/// mass matrix by the element inverse mass.  Returning the H1 value at those
/// points therefore performs the L2 projection onto the DG space of the
/// simulation, restricted to the elements carrying the Lorentz material tag.
///
/// For a 2D source mesh the query point (x, y, z) is reduced to (x, r) with
/// r = sqrt(y^2 + z^2), which is exact for axisymmetric data.  A 3D source
/// mesh is located directly.  ``scale`` converts the file value to solver
/// units (1 / c0 for rad/s files); ``outside_value`` is returned when the
/// point is not found in the source mesh (default 0).
class GridFunctionProfile : public Function {
public:
	GridFunctionProfile(
		const std::string& mesh_path,
		const std::string& gf_path,
		bool axisymmetric,
		double scale,
		double outside_value = 0.0,
		bool fallback_to_nearest = false)
		: data_{ load(mesh_path, gf_path, axisymmetric, scale, outside_value,
		              fallback_to_nearest) }
	{
	}

	int dimension() const override
	{
		return data_->reduce_axisymmetric ? 3 : data_->dimension;
	}

	std::unique_ptr<Function> clone() const override
	{
		return std::make_unique<GridFunctionProfile>(*this);
	}

	double eval(const Position& pos) const override
	{
		auto& d = *data_;
		mfem::Vector query(d.query_dimension);
		query = 0.0;
		if (d.reduce_axisymmetric && pos.Size() >= 3) {
			query(0) = pos(0);
			query(1) = std::hypot(pos(1), pos(2));
		}
		else {
			const int n = std::min(d.query_dimension, pos.Size());
			for (int i = 0; i < n; ++i) { query(i) = pos(i); }
		}

		int element = -1;
		mfem::IntegrationPoint ip;
		if (!d.locate(query, element, ip)) {
			if (d.fallback_to_nearest) {
				const int nearest = d.kdtree->FindClosestPoint(query.GetData());
				if (nearest >= 0 && nearest < d.mesh.GetNE()) {
					const mfem::IntegrationPoint& center =
						mfem::Geometries.GetCenter(
							d.mesh.GetElementBaseGeometry(nearest));
					return d.scale * d.gf.GetValue(nearest, center);
				}
			}
			return d.scale * d.outside_value;
		}
		return d.scale * d.gf.GetValue(element, ip);
	}

	/// Maximum and minimum nodal values of the source grid function, in file
	/// units (before ``scale``).
	double maxValue() const { return data_->max_value; }
	double minValue() const { return data_->min_value; }

	const std::string& meshPath() const { return data_->mesh_path; }
	const std::string& fieldPath() const { return data_->gf_path; }
	bool axisymmetric() const { return data_->reduce_axisymmetric; }

private:
	struct Data {
		mfem::Mesh mesh;
		mfem::GridFunction gf;
		std::unique_ptr<mfem::KDTreeBase<int, mfem::real_t>> kdtree;
		std::vector<std::vector<int>> vertex_to_elements;
		std::string mesh_path;
		std::string gf_path;
		bool reduce_axisymmetric = false;
		int query_dimension = 3;
		int dimension = 3;
		double scale = 1.0;
		double outside_value = 0.0;
		bool fallback_to_nearest = false;
		double max_value = 0.0;
		double min_value = 0.0;

		Data(
			const std::string& mesh_file,
			const std::string& gf_file,
			bool axisymmetric,
			double value_scale,
			double outside,
			bool nearest_fallback)
			: mesh{mfem::Mesh::LoadFromFile(mesh_file, 1, 0, true)}
			, gf{loadGridFunction(mesh, gf_file, mesh_file)}
			, mesh_path{mesh_file}
			, gf_path{gf_file}
			, scale{value_scale}
			, outside_value{outside}
			, fallback_to_nearest{nearest_fallback}
		{
			const int sdim = mesh.SpaceDimension();
			reduce_axisymmetric = axisymmetric && mesh.Dimension() == 2;
			query_dimension = sdim;
			dimension = sdim;
			buildLocator();
			computeRange();
		}

		static mfem::GridFunction loadGridFunction(
			mfem::Mesh& m, const std::string& path, const std::string& mesh_path)
		{
			std::ifstream input(path);
			if (!input) {
				throw std::runtime_error(
					"GridFunctionProfile: cannot open '" + path + "'.");
			}
			try {
				return mfem::GridFunction(&m, input);
			}
			catch (const std::exception& e) {
				throw std::runtime_error(
					"GridFunctionProfile: could not read '" + path +
					"' on mesh '" + mesh_path + "': " + e.what());
			}
		}

		void buildLocator()
		{
			const int ne = mesh.GetNE();
			if (ne == 0) {
				throw std::runtime_error(
					"GridFunctionProfile: mesh '" + mesh_path + "' has no elements.");
			}
			const int sdim = mesh.SpaceDimension();
			mfem::Vector center(sdim);
			if (sdim == 2) {
				kdtree.reset(new mfem::KDTree2D);
			}
			else if (sdim == 3) {
				kdtree.reset(new mfem::KDTree3D);
			}
			else {
				throw std::runtime_error(
					"GridFunctionProfile: only 2D and 3D source meshes are supported.");
			}
			for (int e = 0; e < ne; ++e) {
				mfem::ElementTransformation* tr = mesh.GetElementTransformation(e);
				const mfem::IntegrationPoint& ref =
					mfem::Geometries.GetCenter(mesh.GetElementBaseGeometry(e));
				tr->Transform(ref, center);
				kdtree->AddPoint(center.GetData(), e);
			}
			kdtree->Sort();

			mfem::Array<int> vertices;
			vertex_to_elements.assign(mesh.GetNV(), {});
			for (int e = 0; e < ne; ++e) {
				mesh.GetElementVertices(e, vertices);
				for (int i = 0; i < vertices.Size(); ++i) {
					vertex_to_elements[vertices[i]].push_back(e);
				}
			}
		}

		void computeRange()
		{
			if (gf.Size() == 0) {
				max_value = min_value = 0.0;
				return;
			}
			max_value = min_value = gf(0);
			for (int i = 1; i < gf.Size(); ++i) {
				max_value = std::max(max_value, gf(i));
				min_value = std::min(min_value, gf(i));
			}
		}

		bool tryElement(const mfem::Vector& x, int element, mfem::IntegrationPoint& ip)
		{
			mfem::InverseElementTransformation inv;
			inv.SetTransformation(*mesh.GetElementTransformation(element));
			const int result = inv.Transform(x, ip);
			return result == mfem::InverseElementTransformation::Inside;
		}

		bool locate(const mfem::Vector& x, int& element, mfem::IntegrationPoint& ip)
		{
			const int nearest = kdtree->FindClosestPoint(x.GetData());
			if (nearest >= 0 && nearest < mesh.GetNE() && tryElement(x, nearest, ip)) {
				element = nearest;
				return true;
			}

			// Walk the vertex neighbourhood (two rings) of the closest element.
			std::vector<int> frontier{ nearest };
			std::vector<char> visited(mesh.GetNE(), 0);
			if (nearest >= 0 && nearest < mesh.GetNE()) { visited[nearest] = 1; }
			mfem::Array<int> vertices;
			for (int ring = 0; ring < 2 && !frontier.empty(); ++ring) {
				std::vector<int> next;
				for (int e : frontier) {
					if (e < 0 || e >= mesh.GetNE()) { continue; }
					if (e != nearest && tryElement(x, e, ip)) {
						element = e;
						return true;
					}
					mesh.GetElementVertices(e, vertices);
					for (int i = 0; i < vertices.Size(); ++i) {
						for (int candidate : vertex_to_elements[vertices[i]]) {
							if (!visited[candidate]) {
								visited[candidate] = 1;
								next.push_back(candidate);
							}
						}
					}
				}
				frontier.swap(next);
			}
			return false;
		}
	};

	std::shared_ptr<Data> data_;

	static std::shared_ptr<Data> load(
		const std::string& mesh_path,
		const std::string& gf_path,
		bool axisymmetric,
		double scale,
		double outside_value,
		bool fallback_to_nearest)
	{
		namespace fs = std::filesystem;
		if (!fs::exists(mesh_path)) {
			throw std::runtime_error(
				"GridFunctionProfile: source mesh '" + mesh_path + "' does not exist.");
		}
		if (!fs::exists(gf_path)) {
			throw std::runtime_error(
				"GridFunctionProfile: grid function '" + gf_path + "' does not exist.");
		}
		return std::make_shared<Data>(
			mesh_path, gf_path, axisymmetric, scale, outside_value,
			fallback_to_nearest);
	}
};

} // namespace maxwell
