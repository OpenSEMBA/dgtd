#pragma once

#include <mfem.hpp>

namespace maxwell {
namespace mfemExtension{

using namespace mfem;
class BilinearFormTF : public BilinearForm
{
public:

	/// Creates bilinear form associated with FE space @a *f.
	/** The pointer @a f is not owned by the newly constructed object. */
	BilinearFormTF(FiniteElementSpace* f);

	/// Access all the integrators added with AddInteriorBoundaryFaceIntegrator().
	Array<BilinearFormIntegrator*>* GetIBFI() { return &interior_boundary_face_integs; }

	/** @brief Access all boundary markers added with AddInteriorBoundaryFaceIntegrator().*/
	Array<Array<int>*>* GetIBFI_Marker() { return &interior_boundary_face_integs_marker; }
	
	/// Adds new Boundary Integrator restricted to certain Interior faces specified by
	/// the @a int_attr_marker.
	void AddInteriorBoundaryFaceIntegrator(BilinearFormIntegrator* bfi,
		Array<int>& int_bdr_marker);

	/// Assembles the form i.e. sums over all domain/bdr integrators.
	void Assemble(int skip_zeros = 1);

protected:

	/// Interior Boundary integrators.
	Array<BilinearFormIntegrator*> interior_boundary_face_integs;
	Array<Array<int>*> interior_boundary_face_integs_marker; ///< Entries are not owned.


private:

	/// Copy construction is not supported; body is undefined.
	BilinearFormTF(const BilinearFormTF&);

	/// Copy assignment is not supported; body is undefined.
	BilinearFormTF& operator=(const BilinearFormTF&);

};

}
}