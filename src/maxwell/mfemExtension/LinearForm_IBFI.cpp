#include "LinearForm_IBFI.hpp"

namespace maxwell {
namespace mfemExtension {

using namespace mfem;

LinearFormIBFI::LinearFormIBFI(FiniteElementSpace* f) : LinearForm(f)
{
}

void LinearFormIBFI::AddInteriorBoundaryFaceIntegrator(LinearFormIntegrator* lfi,
    Array<int>& int_bdr_marker)
{
    interior_boundary_face_integs.Append(lfi);
    interior_boundary_face_integs_marker.Append(&int_bdr_marker);
}

void LinearFormIBFI::Assemble() {

    Mesh* mesh = fes->GetMesh();

    LinearForm::Assemble();

    if (interior_boundary_face_integs.Size())
    {
        FaceElementTransformations* tr;
        Mesh* mesh = fes->GetMesh();

        Array<int> vdofs, vdofs2;
        const FiniteElement* fe1, * fe2;
        Vector elemvect, elem2vect;

        // Which interior boundary attributes need to be processed?
        Array<int> int_bdr_attr_marker(mesh->bdr_attributes.Size() ?
            mesh->bdr_attributes.Max() : 0);
        int_bdr_attr_marker = 0;
        for (int k = 0; k < interior_boundary_face_integs.Size(); k++)
        {
            if (interior_boundary_face_integs_marker[k] == NULL)
            {
                int_bdr_attr_marker = 1;
                break;
            }
            Array<int>& int_bdr_marker = *interior_boundary_face_integs_marker[k];
            MFEM_ASSERT(int_bdr_marker.Size() == int_bdr_attr_marker.Size(),
                "invalid boundary marker for interior boundary face integrator #"
                << k << ", counting from zero");
            for (int i = 0; i < int_bdr_attr_marker.Size(); i++)
            {
                int_bdr_attr_marker[i] |= int_bdr_marker[i];
            }
        }

        for (int i = 0; i < mesh->GetNBE(); i++)
        {
            const int bdr_attr = mesh->GetBdrAttribute(i);
            if (int_bdr_attr_marker[bdr_attr - 1] == 0) { continue; }

            tr = mesh->GetInteriorFaceTransformations(mesh->GetBdrFace(i));
            if (tr != NULL)
            {
                fes->GetElementVDofs(tr->Elem1No, vdofs);
                fes->GetElementVDofs(tr->Elem2No, vdofs2);
                fe1 = fes->GetFE(tr->Elem1No);
                fe2 = fes->GetFE(tr->Elem2No);
                for (int k = 0; k < interior_boundary_face_integs.Size(); k++)
                {
                    if (interior_boundary_face_integs_marker[k] &&
                        (*interior_boundary_face_integs_marker[k])[bdr_attr - 1] == 0)
                    {
                        continue;
                    }

                    interior_boundary_face_integs[k]->AssembleRHSElementVect(*fe1, *tr, elemvect);
                    AddElementVector(vdofs,  elemvect);

                    interior_boundary_face_integs[k]->AssembleRHSElementVect(*fe2, *tr, elem2vect);
                    AddElementVector(vdofs2, elem2vect);
                }
            }
        }
    }
}

}
}