#include "BilinearForm.h"

namespace maxwell {
namespace mfemExtension {

using namespace mfem;

BilinearFormTF::BilinearFormTF(FiniteElementSpace* f) : BilinearForm(f)
{
}

void BilinearFormTF::Assemble(int skip_zeros) {
    if (ext)
    {
        ext->Assemble();
        return;
    }

    ElementTransformation* eltrans;
    DofTransformation* doftrans;
    Mesh* mesh = fes->GetMesh();
    DenseMatrix elmat, * elmat_p;

    if (mat == NULL)
    {
        AllocMat();
    }

#ifdef MFEM_USE_LEGACY_OPENMP
    int free_element_matrices = 0;
    if (!element_matrices)
    {
        ComputeElementMatrices();
        free_element_matrices = 1;
    }
#endif

    if (domain_integs.Size())
    {
        for (int k = 0; k < domain_integs.Size(); k++)
        {
            if (domain_integs_marker[k] != NULL)
            {
                MFEM_VERIFY(domain_integs_marker[k]->Size() ==
                    (mesh->attributes.Size() ? mesh->attributes.Max() : 0),
                    "invalid element marker for domain integrator #"
                    << k << ", counting from zero");
            }
        }

        for (int i = 0; i < fes->GetNE(); i++)
        {
            int elem_attr = fes->GetMesh()->GetAttribute(i);
            doftrans = fes->GetElementVDofs(i, vdofs);
            if (element_matrices)
            {
                elmat_p = &(*element_matrices)(i);
            }
            else
            {
                elmat.SetSize(0);
                for (int k = 0; k < domain_integs.Size(); k++)
                {
                    if (domain_integs_marker[k] == NULL ||
                        (*(domain_integs_marker[k]))[elem_attr - 1] == 1)
                    {
                        const FiniteElement& fe = *fes->GetFE(i);
                        eltrans = fes->GetElementTransformation(i);
                        domain_integs[k]->AssembleElementMatrix(fe, *eltrans, elemmat);
                        if (elmat.Size() == 0)
                        {
                            elmat = elemmat;
                        }
                        else
                        {
                            elmat += elemmat;
                        }
                    }
                }
                if (elmat.Size() == 0)
                {
                    continue;
                }
                else
                {
                    elmat_p = &elmat;
                }
                if (doftrans)
                {
                    doftrans->TransformDual(elmat);
                }
                elmat_p = &elmat;
            }
            if (static_cond)
            {
                static_cond->AssembleMatrix(i, *elmat_p);
            }
            else
            {
                mat->AddSubMatrix(vdofs, vdofs, *elmat_p, skip_zeros);
                if (hybridization)
                {
                    hybridization->AssembleMatrix(i, *elmat_p);
                }
            }
        }
    }

    if (boundary_integs.Size())
    {
        // Which boundary attributes need to be processed?
        Array<int> bdr_attr_marker(mesh->bdr_attributes.Size() ?
            mesh->bdr_attributes.Max() : 0);
        bdr_attr_marker = 0;
        for (int k = 0; k < boundary_integs.Size(); k++)
        {
            if (boundary_integs_marker[k] == NULL)
            {
                bdr_attr_marker = 1;
                break;
            }
            Array<int>& bdr_marker = *boundary_integs_marker[k];
            MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
                "invalid boundary marker for boundary integrator #"
                << k << ", counting from zero");
            for (int i = 0; i < bdr_attr_marker.Size(); i++)
            {
                bdr_attr_marker[i] |= bdr_marker[i];
            }
        }

        for (int i = 0; i < fes->GetNBE(); i++)
        {
            const int bdr_attr = mesh->GetBdrAttribute(i);
            if (bdr_attr_marker[bdr_attr - 1] == 0) { continue; }

            const FiniteElement& be = *fes->GetBE(i);
            doftrans = fes->GetBdrElementVDofs(i, vdofs);
            eltrans = fes->GetBdrElementTransformation(i);
            int k = 0;
            for (; k < boundary_integs.Size(); k++)
            {
                if (boundary_integs_marker[k] &&
                    (*boundary_integs_marker[k])[bdr_attr - 1] == 0) {
                    continue;
                }

                boundary_integs[k]->AssembleElementMatrix(be, *eltrans, elmat);
                k++;
                break;
            }
            for (; k < boundary_integs.Size(); k++)
            {
                if (boundary_integs_marker[k] &&
                    (*boundary_integs_marker[k])[bdr_attr - 1] == 0) {
                    continue;
                }

                boundary_integs[k]->AssembleElementMatrix(be, *eltrans, elemmat);
                elmat += elemmat;
            }
            if (doftrans)
            {
                doftrans->TransformDual(elmat);
            }
            elmat_p = &elmat;
            if (!static_cond)
            {
                mat->AddSubMatrix(vdofs, vdofs, *elmat_p, skip_zeros);
                if (hybridization)
                {
                    hybridization->AssembleBdrMatrix(i, *elmat_p);
                }
            }
            else
            {
                static_cond->AssembleBdrMatrix(i, *elmat_p);
            }
        }
    }

    if (interior_face_integs.Size())
    {
        FaceElementTransformations* tr;
        Array<int> vdofs2;

        int nfaces = mesh->GetNumFaces();
        for (int i = 0; i < nfaces; i++)
        {
            tr = mesh->GetInteriorFaceTransformations(i);
            if (tr != NULL)
            {
                fes->GetElementVDofs(tr->Elem1No, vdofs);
                fes->GetElementVDofs(tr->Elem2No, vdofs2);
                vdofs.Append(vdofs2);
                for (int k = 0; k < interior_face_integs.Size(); k++)
                {
                    interior_face_integs[k]->
                        AssembleFaceMatrix(*fes->GetFE(tr->Elem1No),
                            *fes->GetFE(tr->Elem2No),
                            *tr, elemmat);
                    mat->AddSubMatrix(vdofs, vdofs, elemmat, skip_zeros);
                }
            }
        }
    }

    if (boundary_face_integs.Size())
    {
        FaceElementTransformations* tr;
        Array<int> vdofs2;
        const FiniteElement* fe1, * fe2;

        // Which boundary attributes need to be processed?
        Array<int> bdr_attr_marker(mesh->bdr_attributes.Size() ?
            mesh->bdr_attributes.Max() : 0);
        bdr_attr_marker = 0;
        for (int k = 0; k < boundary_face_integs.Size(); k++)
        {
            if (boundary_face_integs_marker[k] == NULL)
            {
                bdr_attr_marker = 1;
                break;
            }
            Array<int>& bdr_marker = *boundary_face_integs_marker[k];
            MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
                "invalid boundary marker for boundary face integrator #"
                << k << ", counting from zero");
            for (int i = 0; i < bdr_attr_marker.Size(); i++)
            {
                bdr_attr_marker[i] |= bdr_marker[i];
            }
        }

        for (int i = 0; i < fes->GetNBE(); i++)
        {
            const int bdr_attr = mesh->GetBdrAttribute(i);
            if (bdr_attr_marker[bdr_attr - 1] == 0) { continue; }

            tr = mesh->GetInteriorFaceTransformations(i);
            if (tr != NULL)
            {
                fes->GetElementVDofs(tr->Elem1No, vdofs);
                fe1 = fes->GetFE(tr->Elem1No);
                // The fe2 object is really a dummy and not used on the boundaries,
                // but we can't dereference a NULL pointer, and we don't want to
                // actually make a fake element.
                fe2 = fe1;
                for (int k = 0; k < boundary_face_integs.Size(); k++)
                {
                    if (boundary_face_integs_marker[k] &&
                        (*boundary_face_integs_marker[k])[bdr_attr - 1] == 0)
                    {
                        continue;
                    }

                    boundary_face_integs[k]->AssembleFaceMatrix(*fe1, *fe2, *tr,
                        elemmat);
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