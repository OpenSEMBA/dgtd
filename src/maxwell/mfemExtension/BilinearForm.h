#pragma once

#include <mfem.hpp>

namespace maxwell {
namespace mfemExtension{

using namespace mfem;
class BilinearFormTF : public BilinearForm
{
public:


	BilinearFormTF(FiniteElementSpace* f);

	Array<BilinearFormIntegrator*>* GetIBFI() { return &interior_boundary_face_integs; }
	Array<Array<int>*>* GetIBFI_Marker() { return &interior_boundary_face_integs_marker; }

	void AddInteriorBoundaryFaceIntegrator(BilinearFormIntegrator* bfi,
		Array<int>& int_bdr_marker);

	void Assemble(int skip_zeros = 1);

protected:

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