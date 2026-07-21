#include "RCSSurfaceExporter.h"

#include "general/forall.hpp"

namespace maxwell {

void rcs_gather_surface_fields_gpu(const mfem::Array<int>& parent_dof_ids,
                                   const mfem::Vector& parent_dof_signs,
                                   const Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction>& globalFields,
                                   Fields<mfem::FiniteElementSpace, mfem::GridFunction>& surfaceFields,
                                   int num_dofs)
{
    if (num_dofs <= 0) {
        return;
    }

    const int* parent_ids = parent_dof_ids.Read();
    const double* signs = parent_dof_signs.Read();

    const double* ex = globalFields.get(E, X).Read();
    const double* ey = globalFields.get(E, Y).Read();
    const double* ez = globalFields.get(E, Z).Read();
    const double* hx = globalFields.get(H, X).Read();
    const double* hy = globalFields.get(H, Y).Read();
    const double* hz = globalFields.get(H, Z).Read();

    double* sex = surfaceFields.get(E, X).Write();
    double* sey = surfaceFields.get(E, Y).Write();
    double* sez = surfaceFields.get(E, Z).Write();
    double* shx = surfaceFields.get(H, X).Write();
    double* shy = surfaceFields.get(H, Y).Write();
    double* shz = surfaceFields.get(H, Z).Write();

    mfem::forall(num_dofs, [=] MFEM_DEVICE(int i) {
        const int p = parent_ids[i];
        const double s = signs[i];
        sex[i] = s * ex[p];
        sey[i] = s * ey[p];
        sez[i] = s * ez[p];
        shx[i] = s * hx[p];
        shy[i] = s * hy[p];
        shz[i] = s * hz[p];
    });
}

} // namespace maxwell
