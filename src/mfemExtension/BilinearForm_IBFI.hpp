#pragma once

#include <mfem.hpp>

namespace maxwell {
namespace mfemExtension{

using namespace mfem;
class BilinearFormIBFI : public BilinearForm
{
public:

	/// Creates bilinear form associated with FE space @a *f.
	/** The pointer @a f is not owned by the newly constructed object. */
	BilinearFormIBFI(FiniteElementSpace* f);

	/// Access all the integrators added with AddInteriorBoundaryFaceIntegrator().
	Array<BilinearFormIntegrator*>* GetIBFI() { return &internal_boundary_face_integs; }

	/** @brief Access all boundary markers added with AddInteriorBoundaryFaceIntegrator().*/
	Array<Array<int>*>* GetIBFI_Marker() { return &internal_boundary_face_integs_marker; }
	
	/// Adds new Boundary Integrator restricted to certain Interior faces specified by
	/// the @a int_attr_marker.
	void AddInteriorBoundaryFaceIntegrator(BilinearFormIntegrator* bfi,
		Array<int>& internal_bdr_attr_marker);

	/// Adds new Boundary Integrator restricted to certain Interior faces.
	void AddInternalBoundaryFaceIntegrator(BilinearFormIntegrator* bfi);

	/// Assembles the form i.e. sums over all domain/bdr integrators.
	void Assemble(int skip_zeros = 1);

protected:

	/// Interior Boundary integrators.
	Array<BilinearFormIntegrator*> internal_boundary_face_integs;
	Array<Array<int>*> internal_boundary_face_integs_marker; ///< Entries are not owned.


private:

	/// Copy construction is not supported; body is undefined.
	BilinearFormIBFI(const BilinearFormIBFI&);

	/// Copy assignment is not supported; body is undefined.
	BilinearFormIBFI& operator=(const BilinearFormIBFI&);

};

}
}