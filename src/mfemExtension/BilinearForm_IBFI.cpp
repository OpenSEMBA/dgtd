#include "BilinearForm_IBFI.hpp"
#include <cstddef>

namespace maxwell {
namespace mfemExtension {

using namespace mfem;

BilinearFormIBFI::BilinearFormIBFI(FiniteElementSpace* f) : BilinearForm(f)
{
}

void BilinearFormIBFI::AddInternalBoundaryFaceIntegrator(BilinearFormIntegrator
    * bfi)
{
    internal_boundary_face_integs.Append(bfi);
    // nullptr -> all attributes are active
    internal_boundary_face_integs_marker.Append(nullptr);
}

void BilinearFormIBFI::AddInteriorBoundaryFaceIntegrator(BilinearFormIntegrator* bfi,
     Array<int>& internal_bdr_attr_marker)
{
    internal_boundary_face_integs.Append(bfi);
    internal_boundary_face_integs_marker.Append(&internal_bdr_attr_marker);
}

void BilinearFormIBFI::Assemble(int skip_zeros) {
    
    Mesh* mesh = fes->GetMesh();
    
    BilinearForm::Assemble(skip_zeros);

    if (internal_boundary_face_integs.Size())
    {
        // Which internal boundary attributes need to be processed?
        Array<int> bdr_attr_marker(mesh->bdr_attributes.Size() ?
            mesh->bdr_attributes.Max() : 0);
        bdr_attr_marker = 0;
        for (int k = 0; k < internal_boundary_face_integs.Size(); k++)
        {
            if (internal_boundary_face_integs_marker[k] == NULL)
            {
                bdr_attr_marker = 1;
                break;
            }
            auto& bdr_marker = *internal_boundary_face_integs_marker[k];
            MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
                "invalid boundary marker for internal boundary face "
                "integrator #" << k << ", counting from zero");
            for (int i = 0; i < bdr_attr_marker.Size(); i++)
            {
                bdr_attr_marker[i] |= bdr_marker[i];
            }
        }

        Array<int> vdofs2;
        for (int i = 0; i < mesh->GetNBE(); i++)
        {
            const int bdr_attr = mesh->GetBdrAttribute(i);
            if (bdr_attr_marker[bdr_attr - 1] == 0) { continue; }

            auto* tr = mesh->GetInternalBdrFaceTransformations(i);
            if (tr != nullptr)
            {
                fes->GetElementVDofs(tr->Elem1No, vdofs);
                fes->GetElementVDofs(tr->Elem2No, vdofs2);
                vdofs.Append(vdofs2);
                const auto* fe1 = fes->GetFE(tr->Elem1No);
                const auto* fe2 = fes->GetFE(tr->Elem2No);
                for (int k = 0; k < internal_boundary_face_integs.Size(); k++)
                {
                    if (internal_boundary_face_integs_marker[k] &&
                        (*internal_boundary_face_integs_marker[k])[bdr_attr - 1] == 0)
                    {
                        continue;
                    }

                    internal_boundary_face_integs[k]->AssembleFaceMatrix(
                        *fe1, *fe2, *tr, elemmat);
                    mat->AddSubMatrix(vdofs, vdofs, elemmat, skip_zeros);
                }
            }
        }
    }

#ifdef MFEM_USE_LEGACY_OPENMP
    if (free_element_matrices)
    {
        FreeElementMatrices();
    }
#endif
}

}
}