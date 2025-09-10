#pragma once

#include <mfem.hpp>
#include <math.h>
#include <filesystem>
#include <string>
#include "components/Types.h"
#include "math/Function.h"

namespace maxwell {

using namespace mfem;
using Position = mfem::Vector;

enum FunctionType 
{
    TimeResonantSin = 0
};

class TimeFunction {
public:

	virtual ~TimeFunction() = default;

	virtual std::unique_ptr<TimeFunction> clone() const = 0;

	virtual double eval(const Position&, const Time&) const = 0;

};

class TimeResonantSinusoidalMode : public TimeFunction
{
public:
	TimeResonantSinusoidalMode(const std::vector<std::size_t>& modes, const std::vector<double>& box_size)
    {
        int dim = modes.size();
        if(modes_.size() != box_size_.size()){
            modes_.size() > box_size_.size() ? dim = modes_.size() : dim = box_size_.size();
        }
        modes_.resize(dim);
        box_size_.resize(dim);
        for (auto d = 0; d < dim; d++){
            modes_[d] = modes[d];
            box_size_[d] = box_size[d];
        }
    }

	std::unique_ptr<TimeFunction> clone() const {
		return std::make_unique<TimeResonantSinusoidalMode>(*this);
	}

	double eval(const Position& pos, const Time& t) const
	{
		double w = M_PI;
        for (auto d = 0; d < modes_.size(); d++){
            w *= std::pow(modes_[d] / box_size_[d], 2);
        }

        double res = std::cos(w * t);
        for (auto d = 0; d < modes_.size(); d++){
            res *= std::sin(modes_[d] * M_PI * pos[d] / box_size_[d]);
        }

        return res;
	}

private:
	std::vector<std::size_t> modes_;
    std::vector<double> box_size_;
};

double evalTimeResonantSinMode(const Position& pos, const Time& time, const std::vector<int>& modes, const std::vector<double>& box_size)
{
    return 0.0;
}

Mesh loadMeshFromFile(const std::string& mesh_path);
GridFunction loadGridFunctionFromFile(const std::string& file_path, Mesh& mesh);
int getRankAmount(const std::string& data_path);

class L2SimDataCalculator
{
    L2SimDataCalculator(const std::string& data_path, const FunctionType function_type);

private:

    void assignFunctionType(const FunctionType ft);
    void loadMeshes(const std::string& data_path);

    std::vector<Mesh> meshes_;
    std::vector<FiniteElementSpace> feses_;
    std::unique_ptr<TimeFunction> function_;
};


}

