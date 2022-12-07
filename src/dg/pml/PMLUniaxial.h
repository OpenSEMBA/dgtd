#pragma once

#include "PML.h"

template<int D>
class DGPMLUniaxial : public PML {
public:
    DGPMLUniaxial(
            const PMVolumePML& mat,
            const CellGroup& cells,
            const bool useConductivity,
            const Math::Real conductivity);
    virtual ~DGPMLUniaxial() = default;
    void addRHSToRes(
            const size_t e1, const size_t e2,
            const Math::Real rka, const Math::Real dt);
    void updateWithRes(
            const size_t e1,
            const size_t e2,
            const Math::Real rkb);
    void computeRHSElectric(
            FieldR3& rhsE,
            const FieldR3& E,
            const size_t e1, const size_t e2) const;
    void computeRHSMagnetic(
            FieldR3& rhsH,
            const FieldR3& H,
            const size_t e1, const size_t e2) const;
    void computeRHSElectricPolarizationCurrents(
            const FieldR3& E,
            const size_t e1, const size_t e2);
    void computeRHSMagneticPolarizationCurrents(
            const FieldR3& H,
            const size_t e1, const size_t e2);
protected:
    FieldR1 J, M, resJ, resM, rhsJ, rhsM;
private:
    bool check() const;
    static const CartesianAxis dir1 = CartesianAxis(D);
    static const CartesianAxis dir2 = CartesianAxis((D + 1)%3);
    static const CartesianAxis dir3 = CartesianAxis((D + 2)%3);
};


template<Math::Int D>
DGPMLUniaxial<D>::DGPMLUniaxial(
    const PMVolumePML& mat,
    const CellGroup& cells,
    const bool useConductivity,
    const Math::Real conductivity) :
    DGPML(mat) {
    J.set(dof, 0.0);
    resJ.set(dof, 0.0);
    rhsJ.set(dof, 0.0);
    M.set(dof, 0.0);
    resM.set(dof, 0.0);
    rhsM.set(dof, 0.0);
    assert(check());
}

template<Math::Int D>
void DGPMLUniaxial<D>::addRHSToRes(
    const size_t e1, const size_t e2,
    const Math::Real rka, const Math::Real dt) {
    //    size_t i,e;
    //#ifdef SOLVER_USE_OPENMP
    //#pragma omp parallel for private(i,e)
    //#endif
    //    for (i = 0; i < dof; i++) {
    //        e = elem[(i / np) % nElem];
    //        if (e1 <= e && e < e2) {
    //            resJ[i] *= rka;
    //            resJ[i] += rhsJ[i] * dt;
    //            resM[i] *= rka;
    //            resM[i] += rhsM[i] * dt;
    //        }
    //    }
}

template<Math::Int D>
void DGPMLUniaxial<D>::computeRHSElectric(
    FieldR3& rhsE, const FieldR3& E,
    const size_t e1, const size_t e2) const {
    //    if (useConstantConductivity) {
    //        size_t i, j, e, n;
    //#ifdef SOLVER_USE_OPENMP
    //#pragma omp parallel for private(i,j,e,n)
    //#endif
    //        for (i = 0; i < dof; i++) {
    //            e = elem[(i / np) % nElem];
    //            if (e1 <= e && e < e2) {
    //                n = i % np;
    //                j = e * np + n ;
    //                rhsE(dir1)[j] += (Constants::eps0*sig) * E(dir1)[j] - Constants::eps0 * J[i];
    //                rhsE(dir2)[j] -= (Constants::eps0*sig) * E(dir2)[j];
    //                rhsE(dir3)[j] -= (Constants::eps0*sig) * E(dir3)[j];
    //            }
    //        }
    //    } else {
    //        size_t i,j,e;
    //#ifdef SOLVER_USE_OPENMP
    //#pragma omp parallel for private(i,j,e)
    //#endif
    //        for (e = 0; e < nElem; e++) {
    //            if (e1 <= elem[e] && elem[e] < e2) {
    //                i = e * np;
    //                j = elem[e] * np;
    //                //rhsE1[j] += (Constants::eps0*sig1) * E1[j] - Constants::eps0 * J[i];
    //                add_am_v_prod<Math::Real,np,np>(&rhsE.set(dir1)[j], sig1[e], &E(dir1)[j], Constants::eps0);
    //                sub_a_v_prod<Math::Real,np>(&rhsE(dir1)[j], &J[i], Constants::eps0);
    //                //rhsE2[j] -= (Constants::eps0*sig1) * E2[j];
    //                sub_am_v_prod<Math::Real,np,np>(&rhsE(dir2)[j], sig1[e], &E(dir2)[j], Constants::eps0);
    //                //rhsE3[j] -= (Constants::eps0*sig1) * E3[j];
    //                sub_am_v_prod<Math::Real,np,np>(&rhsE(dir3)[j], sig1[e], &E(dir3)[j], Constants::eps0);
    //            }
    //        }
    //    }
}

template<Math::Int D>
void DGPMLUniaxial<D>::computeRHSMagnetic(
    FieldR3& rhsH, const FieldR3& H,
    const size_t e1, const size_t e2) const {
    //    if (useConstantConductivity) {
    //        size_t i, j, e, n;
    //#ifdef SOLVER_USE_OPENMP
    //#pragma omp parallel for private(i,j,e,n)
    //#endif
    //        for (i = 0; i < dof; i++) {
    //            e = elem[(i / np) % nElem];
    //            if (e1 <= e && e < e2) {
    //                n = i % np;
    //                j = e * np + n ;
    //                rhsH(dir1)[j] += (Constants::mu0*sig) * H(dir1)[j] - Constants::mu0 * M[i];
    //                rhsH(dir2)[j] -= (Constants::mu0*sig) * H(dir2)[j];
    //                rhsH(dir3)[j] -= (Constants::mu0*sig) * H(dir3)[j];
    //            }
    //        }
    //    } else {
    //        size_t i, j, e;
    //#ifdef SOLVER_USE_OPENMP
    //#pragma omp parallel for private(i,j,e)
    //#endif
    //        for (e = 0; e < nElem; e++) {
    //            if (e1 <= elem[e] && elem[e] < e2) {
    //                i = e * np;
    //                j = elem[e] * np;
    //                //rhsH1[j] += (Constants::mu0*sigma1) * H1[j] - Constants::mu0 * M[i];
    //                add_am_v_prod<Math::Real,np,np>(&rhsH(dir1)[j], sig1[e], &H(dir1)[j], Constants::mu0);
    //                sub_a_v_prod<Math::Real,np>(&rhsH(dir1)[j], &M[i], Constants::mu0);
    //                //rhsH2[j] -= (Constants::mu0*sigma1) * H2[j];
    //                sub_am_v_prod<Math::Real,np,np>(&rhsH(dir2)[j], sig1[e], &H(dir2)[j], Constants::mu0);
    //                //rhsH3[j] -= (Constants::mu0*sigma1) * H3[j];
    //                sub_am_v_prod<Math::Real,np,np>(&rhsH(dir3)[j], sig1[e], &H(dir3)[j], Constants::mu0);
    //            }
    //        }
    //    }
}

template<Math::Int D>
void DGPMLUniaxial<D>::computeRHSElectricPolarizationCurrents(
    const FieldR3& E,
    const size_t e1, const size_t e2) {
    //    if (useConstantConductivity) {
    //        size_t i, j, e, n;
    //#ifdef SOLVER_USE_OPENMP
    //#pragma omp parallel for private(i,j,e,n)
    //#endif
    //        for (i = 0; i < dof; i++) {
    //            e = elem[(i / np) % nElem];
    //            if (e1 <= e && e < e2) {
    //                n = i % np;
    //                j = e * np + n ;
    //                rhsJ[i] = E(dir1)[j] * (sig*sig) - sig * J[i];
    //            }
    //        }
    //    } else {
    //        size_t i,j,e;
    //#ifdef SOLVER_USE_OPENMP
    //#pragma omp parallel for private(i,j,e)
    //#endif
    //        for (e = 0; e < nElem; e++) {
    //            if (e1 <= elem[e] && elem[e] < e2) {
    //                i = e * np;
    //                j = elem[e] * np;
    //                //rhsJ[i] = E1[j] * (sig11) - sig1 * J[i];
    //                m_v_prod<Math::Real,np,np>(&rhsJ[i], sig11[e], &E(dir1)[j]);
    //                sub_m_v_prod<Math::Real,np,np>(&rhsJ[i], sig1[e], &J[i]);
    //            }
    //        }
    //    }
}

template<Math::Int D>
void DGPMLUniaxial<D>::computeRHSMagneticPolarizationCurrents(
    const FieldR3& H,
    const size_t e1, const size_t e2) {
    //    if (useConstantConductivity) {
    //        size_t i, j, e, n;
    //#ifdef SOLVER_USE_OPENMP
    //#pragma omp parallel for private(i,j,e,n)
    //#endif
    //        for (i = 0; i < dof; i++) {
    //            e = elem[(i / np) % nElem];
    //            if (e1 <= e && e < e2) {
    //                n = i % np;
    //                j = e * np + n ;
    //                rhsM[i] = H(dir1)[j] * (sig*sig) - sig * M[i];
    //            }
    //        }
    //    } else {
    //        size_t i, j, e;
    //#ifdef SOLVER_USE_OPENMP
    //#pragma omp parallel for private(i,j,e)
    //#endif
    //        for (e = 0; e < nElem; e++) {
    //            if (e1 <= elem[e] && elem[e] < e2) {
    //                i = e * np;
    //                j = elem[e] * np;
    //                //rhsM[i] = H1[j] * (sigma1*sigma1) - sigma1 * M[i];
    //                m_v_prod<Math::Real,np,np>(&rhsM[i], sig11[e], &H(dir1)[j]);
    //                sub_m_v_prod<Math::Real,np,np>(&rhsM[i], sig1[e], &M[i]);
    //            }
    //        }
    //    }
}

template<Math::Int D>
bool DGPMLUniaxial<D>::check() const {
    bool sigInitialized = true;
    if (!useConstantConductivity) {
        sigInitialized &= (sig1 != NULL);
        sigInitialized &= (sig2 == NULL);
        sigInitialized &= (sig3 == NULL);
        sigInitialized &= (sig11 != NULL);
        sigInitialized &= (sig22 == NULL);
        sigInitialized &= (sig33 == NULL);
        sigInitialized &= (sig12 == NULL);
        sigInitialized &= (sig23 == NULL);
        sigInitialized &= (sig31 == NULL);
    }
    return sigInitialized;
}

template<Math::Int D>
void DGPMLUniaxial<D>::updateWithRes(
    const size_t e1,
    const size_t e2,
    const Math::Real rkb) {
    size_t i, e;
#ifdef SOLVER_USE_OPENMP
#pragma omp parallel for private(i, e)
#endif
    for (i = 0; i < dof; i++) {
        e = elem[(i / np) % nElem];
        if (e1 <= e && e < e2) {
            J[i] += resJ[i] * rkb;
            M[i] += resM[i] * rkb;
        }
    }
}

typedef DGPMLUniaxial<x> DGPMLx;
typedef DGPMLUniaxial<y> DGPMLy;
typedef DGPMLUniaxial<z> DGPMLz;
