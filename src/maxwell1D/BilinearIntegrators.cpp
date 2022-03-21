#include "BilinearIntegrators.h"

namespace Maxwell1D {

/*########################### EA START ###########################*/

static void EADGTraceAssemble1DInt(const int NF,
    const Array<double>& basis,
    const Vector& padata,
    Vector& eadata_int,
    Vector& eadata_ext,
    const bool add)
{
    auto D = Reshape(padata.Read(), 2, 2, NF);
    auto A_int = Reshape(eadata_int.ReadWrite(), 2, NF);
    auto A_ext = Reshape(eadata_ext.ReadWrite(), 2, NF);
    MFEM_FORALL(f, NF,
        {
            double val_int0, val_int1, val_ext01, val_ext10;
            val_int0 = D(0, 0, f);
            val_ext10 = D(1, 0, f);
            val_ext01 = D(0, 1, f);
            val_int1 = D(1, 1, f);
            if (add)
            {
                A_int(0, f) += val_int0;
                A_int(1, f) += val_int1;
                A_ext(0, f) += val_ext01;
                A_ext(1, f) += val_ext10;
            }
            else
            {
                A_int(0, f) = val_int0;
                A_int(1, f) = val_int1;
                A_ext(0, f) = val_ext01;
                A_ext(1, f) = val_ext10;
            }
        });
}

static void EADGTraceAssemble1DBdr(const int NF,
    const Array<double>& basis,
    const Vector& padata,
    Vector& eadata_bdr,
    const bool add)
{
    auto D = Reshape(padata.Read(), 2, 2, NF);
    auto A_bdr = Reshape(eadata_bdr.ReadWrite(), NF);
    MFEM_FORALL(f, NF,
        {
            if (add)
            {
                A_bdr(f) += D(0, 0, f);
            }
            else
            {
                A_bdr(f) = D(0, 0, f);
            }
        });
}

template<int T_D1D = 0, int T_Q1D = 0>
static void EADGTraceAssemble2DInt(const int NF,
    const Array<double>& basis,
    const Vector& padata,
    Vector& eadata_int,
    Vector& eadata_ext,
    const bool add,
    const int d1d = 0,
    const int q1d = 0)
{
    const int D1D = T_D1D ? T_D1D : d1d;
    const int Q1D = T_Q1D ? T_Q1D : q1d;
    MFEM_VERIFY(D1D <= MAX_D1D, "");
    MFEM_VERIFY(Q1D <= MAX_Q1D, "");
    auto B = Reshape(basis.Read(), Q1D, D1D);
    auto D = Reshape(padata.Read(), Q1D, 2, 2, NF);
    auto A_int = Reshape(eadata_int.ReadWrite(), D1D, D1D, 2, NF);
    auto A_ext = Reshape(eadata_ext.ReadWrite(), D1D, D1D, 2, NF);
    MFEM_FORALL_3D(f, NF, D1D, D1D, 1,
        {
            const int D1D = T_D1D ? T_D1D : d1d;
            const int Q1D = T_Q1D ? T_Q1D : q1d;
            MFEM_FOREACH_THREAD(i1,x,D1D)
            {
                MFEM_FOREACH_THREAD(j1,y,D1D)
                {
                    double val_int0 = 0.0;
                    double val_int1 = 0.0;
                    double val_ext01 = 0.0;
                    double val_ext10 = 0.0;
                    for (int k1 = 0; k1 < Q1D; ++k1)
                    {
                    val_int0 += B(k1,i1) * B(k1,j1) * D(k1, 0, 0, f);
                    val_ext01 += B(k1,i1) * B(k1,j1) * D(k1, 0, 1, f);
                    val_ext10 += B(k1,i1) * B(k1,j1) * D(k1, 1, 0, f);
                    val_int1 += B(k1,i1) * B(k1,j1) * D(k1, 1, 1, f);
                    }
                    if (add)
                    {
                    A_int(i1, j1, 0, f) += val_int0;
                    A_int(i1, j1, 1, f) += val_int1;
                    A_ext(i1, j1, 0, f) += val_ext01;
                    A_ext(i1, j1, 1, f) += val_ext10;
                    }
                    else
                    {
                    A_int(i1, j1, 0, f) = val_int0;
                    A_int(i1, j1, 1, f) = val_int1;
                    A_ext(i1, j1, 0, f) = val_ext01;
                    A_ext(i1, j1, 1, f) = val_ext10;
                    }
                }
            }
        });
}

template<int T_D1D = 0, int T_Q1D = 0>
static void EADGTraceAssemble2DBdr(const int NF,
    const Array<double>& basis,
    const Vector& padata,
    Vector& eadata_bdr,
    const bool add,
    const int d1d = 0,
    const int q1d = 0)
{
    const int D1D = T_D1D ? T_D1D : d1d;
    const int Q1D = T_Q1D ? T_Q1D : q1d;
    MFEM_VERIFY(D1D <= MAX_D1D, "");
    MFEM_VERIFY(Q1D <= MAX_Q1D, "");
    auto B = Reshape(basis.Read(), Q1D, D1D);
    auto D = Reshape(padata.Read(), Q1D, 2, 2, NF);
    auto A_bdr = Reshape(eadata_bdr.ReadWrite(), D1D, D1D, NF);
    MFEM_FORALL_3D(f, NF, D1D, D1D, 1,
        {
            const int D1D = T_D1D ? T_D1D : d1d;
            const int Q1D = T_Q1D ? T_Q1D : q1d;
            MFEM_FOREACH_THREAD(i1,x,D1D)
            {
                MFEM_FOREACH_THREAD(j1,y,D1D)
                {
                    double val_bdr = 0.0;
                    for (int k1 = 0; k1 < Q1D; ++k1)
                    {
                    val_bdr += B(k1,i1) * B(k1,j1) * D(k1, 0, 0, f);
                    }
                    if (add)
                    {
                    A_bdr(i1, j1, f) += val_bdr;
                    }
                    else
                    {
                    A_bdr(i1, j1, f) = val_bdr;
                    }
                }
            }
        });
}

template<int T_D1D = 0, int T_Q1D = 0>
static void EADGTraceAssemble3DInt(const int NF,
    const Array<double>& basis,
    const Vector& padata,
    Vector& eadata_int,
    Vector& eadata_ext,
    const bool add,
    const int d1d = 0,
    const int q1d = 0)
{
    const int D1D = T_D1D ? T_D1D : d1d;
    const int Q1D = T_Q1D ? T_Q1D : q1d;
    MFEM_VERIFY(D1D <= MAX_D1D, "");
    MFEM_VERIFY(Q1D <= MAX_Q1D, "");
    auto B = Reshape(basis.Read(), Q1D, D1D);
    auto D = Reshape(padata.Read(), Q1D, Q1D, 2, 2, NF);
    auto A_int = Reshape(eadata_int.ReadWrite(), D1D, D1D, D1D, D1D, 2, NF);
    auto A_ext = Reshape(eadata_ext.ReadWrite(), D1D, D1D, D1D, D1D, 2, NF);
    MFEM_FORALL_3D(f, NF, D1D, D1D, 1,
        {
            const int D1D = T_D1D ? T_D1D : d1d;
            const int Q1D = T_Q1D ? T_Q1D : q1d;
            constexpr int MD1 = T_D1D ? T_D1D : MAX_D1D;
            constexpr int MQ1 = T_Q1D ? T_Q1D : MAX_Q1D;
            double r_B[MQ1][MD1];
            for (int d = 0; d < D1D; d++)
            {
                for (int q = 0; q < Q1D; q++)
                {
                    r_B[q][d] = B(q,d);
                }
            }
            MFEM_SHARED double s_D[MQ1][MQ1][2][2];
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    MFEM_FOREACH_THREAD(k1,x,Q1D)
                    {
                    MFEM_FOREACH_THREAD(k2,y,Q1D)
                    {
                        s_D[k1][k2][i][j] = D(k1,k2,i,j,f);
                    }
                    }
                }
            }
            MFEM_SYNC_THREAD;
            MFEM_FOREACH_THREAD(i1,x,D1D)
            {
                MFEM_FOREACH_THREAD(i2,y,D1D)
                {
                    for (int j1 = 0; j1 < D1D; ++j1)
                    {
                    for (int j2 = 0; j2 < D1D; ++j2)
                    {
                        double val_int0 = 0.0;
                        double val_int1 = 0.0;
                        double val_ext01 = 0.0;
                        double val_ext10 = 0.0;
                        for (int k1 = 0; k1 < Q1D; ++k1)
                        {
                            for (int k2 = 0; k2 < Q1D; ++k2)
                            {
                                val_int0 += r_B[k1][i1] * r_B[k1][j1]
                                            * r_B[k2][i2] * r_B[k2][j2]
                                            * s_D[k1][k2][0][0];
                                val_int1 += r_B[k1][i1] * r_B[k1][j1]
                                            * r_B[k2][i2] * r_B[k2][j2]
                                            * s_D[k1][k2][1][1];
                                val_ext01 += r_B[k1][i1] * r_B[k1][j1]
                                            * r_B[k2][i2] * r_B[k2][j2]
                                            * s_D[k1][k2][0][1];
                                val_ext10 += r_B[k1][i1] * r_B[k1][j1]
                                            * r_B[k2][i2] * r_B[k2][j2]
                                            * s_D[k1][k2][1][0];
                            }
                        }
                        if (add)
                        {
                            A_int(i1, i2, j1, j2, 0, f) += val_int0;
                            A_int(i1, i2, j1, j2, 1, f) += val_int1;
                            A_ext(i1, i2, j1, j2, 0, f) += val_ext01;
                            A_ext(i1, i2, j1, j2, 1, f) += val_ext10;
                        }
                        else
                        {
                            A_int(i1, i2, j1, j2, 0, f) = val_int0;
                            A_int(i1, i2, j1, j2, 1, f) = val_int1;
                            A_ext(i1, i2, j1, j2, 0, f) = val_ext01;
                            A_ext(i1, i2, j1, j2, 1, f) = val_ext10;
                        }
                    }
                    }
                }
            }
        });
}

template<int T_D1D = 0, int T_Q1D = 0>
static void EADGTraceAssemble3DBdr(const int NF,
    const Array<double>& basis,
    const Vector& padata,
    Vector& eadata_bdr,
    const bool add,
    const int d1d = 0,
    const int q1d = 0)
{
    const int D1D = T_D1D ? T_D1D : d1d;
    const int Q1D = T_Q1D ? T_Q1D : q1d;
    MFEM_VERIFY(D1D <= MAX_D1D, "");
    MFEM_VERIFY(Q1D <= MAX_Q1D, "");
    auto B = Reshape(basis.Read(), Q1D, D1D);
    auto D = Reshape(padata.Read(), Q1D, Q1D, 2, 2, NF);
    auto A_bdr = Reshape(eadata_bdr.ReadWrite(), D1D, D1D, D1D, D1D, NF);
    MFEM_FORALL_3D(f, NF, D1D, D1D, 1,
        {
            const int D1D = T_D1D ? T_D1D : d1d;
            const int Q1D = T_Q1D ? T_Q1D : q1d;
            constexpr int MD1 = T_D1D ? T_D1D : MAX_D1D;
            constexpr int MQ1 = T_Q1D ? T_Q1D : MAX_Q1D;
            double r_B[MQ1][MD1];
            for (int d = 0; d < D1D; d++)
            {
                for (int q = 0; q < Q1D; q++)
                {
                    r_B[q][d] = B(q,d);
                }
            }
            MFEM_SHARED double s_D[MQ1][MQ1][2][2];
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    MFEM_FOREACH_THREAD(k1,x,Q1D)
                    {
                    MFEM_FOREACH_THREAD(k2,y,Q1D)
                    {
                        s_D[k1][k2][i][j] = D(k1,k2,i,j,f);
                    }
                    }
                }
            }
            MFEM_SYNC_THREAD;
            MFEM_FOREACH_THREAD(i1,x,D1D)
            {
                MFEM_FOREACH_THREAD(i2,y,D1D)
                {
                    for (int j1 = 0; j1 < D1D; ++j1)
                    {
                    for (int j2 = 0; j2 < D1D; ++j2)
                    {
                        double val_bdr = 0.0;
                        for (int k1 = 0; k1 < Q1D; ++k1)
                        {
                            for (int k2 = 0; k2 < Q1D; ++k2)
                            {
                                val_bdr += r_B[k1][i1] * r_B[k1][j1]
                                        * r_B[k2][i2] * r_B[k2][j2]
                                        * s_D[k1][k2][0][0];
                            }
                        }
                        if (add)
                        {
                            A_bdr(i1, i2, j1, j2, f) += val_bdr;
                        }
                        else
                        {
                            A_bdr(i1, i2, j1, j2, f) = val_bdr;
                        }
                    }
                    }
                }
            }
        });
}

/*########################### PA START ###########################*/

// PA DG Trace Integrator
static void PADGTraceSetup2D(const int Q1D,
    const int NF,
    const Array<double>& w,
    const Vector& det,
    const Vector& nor,
    const Vector& rho,
    const Vector& vel,
    const double alpha,
    const double beta,
    const double gamma,
    Vector& op)
{
    const int VDIM = 2;

    auto d = Reshape(det.Read(), Q1D, NF);
    auto n = Reshape(nor.Read(), Q1D, VDIM, NF);
    const bool const_r = rho.Size() == 1;
    auto R =
        const_r ? Reshape(rho.Read(), 1, 1) : Reshape(rho.Read(), Q1D, NF);
    const bool const_v = vel.Size() == 2;
    auto V =
        const_v ? Reshape(vel.Read(), 2, 1, 1) : Reshape(vel.Read(), 2, Q1D, NF);
    auto W = w.Read();
    auto qd = Reshape(op.Write(), Q1D, 2, 2, NF);

    MFEM_FORALL(f, NF, // can be optimized with Q1D thread for NF blocks
        {
            for (int q = 0; q < Q1D; ++q)
            {
                const double r = const_r ? R(0,0) : R(q,f);
                const double v0 = const_v ? V(0,0,0) : V(0,q,f);
                const double v1 = const_v ? V(1,0,0) : V(1,q,f);
                const double dot = n(q,0,f) * v0 + n(q,1,f) * v1;
                const double negdot = n(q, 0, f) * v0 - n(q, 1, f) * v1;
                const double abs = dot > 0.0 ? dot : -dot;
                const double w = W[q] * r * d(q,f);
                qd(q,0,0,f) = w * (alpha / 2 * dot + gamma * negdot);
                qd(q,1,0,f) = w * (alpha / 2 * dot - gamma * negdot);
                qd(q,0,1,f) = w * (-alpha / 2 * dot - gamma * negdot);
                qd(q,1,1,f) = w * (-alpha / 2 * dot + gamma * negdot);
            }
        });
}

static void PADGTraceSetup3D(const int Q1D,
    const int NF,
    const Array<double>& w,
    const Vector& det,
    const Vector& nor,
    const Vector& rho,
    const Vector& vel,
    const double alpha,
    const double beta,
    const double gamma,
    Vector& op)
{
    const int VDIM = 3;

    auto d = Reshape(det.Read(), Q1D, Q1D, NF);
    auto n = Reshape(nor.Read(), Q1D, Q1D, VDIM, NF);
    const bool const_r = rho.Size() == 1;
    auto R =
        const_r ? Reshape(rho.Read(), 1, 1, 1) : Reshape(rho.Read(), Q1D, Q1D, NF);
    const bool const_v = vel.Size() == 3;
    auto V =
        const_v ? Reshape(vel.Read(), 3, 1, 1, 1) : Reshape(vel.Read(), 3, Q1D, Q1D, NF);
    auto W = w.Read();
    auto qd = Reshape(op.Write(), Q1D, Q1D, 2, 2, NF);

    MFEM_FORALL(f, NF, // can be optimized with Q1D*Q1D threads for NF blocks
        {
            for (int q1 = 0; q1 < Q1D; ++q1)
            {
                for (int q2 = 0; q2 < Q1D; ++q2)
                {
                    const double r = const_r ? R(0,0,0) : R(q1,q2,f);
                    const double v0 = const_v ? V(0,0,0,0) : V(0,q1,q2,f);
                    const double v1 = const_v ? V(1,0,0,0) : V(1,q1,q2,f);
                    const double v2 = const_v ? V(2,0,0,0) : V(2,q1,q2,f);
                    const double dot = n(q1,q2,0,f) * v0 + n(q1,q2,1,f) * v1 +
                        /* */              n(q1,q2,2,f) * v2;
                        const double abs = dot > 0.0 ? dot : -dot;
                        const double w = W[q1 + q2 * Q1D] * r * d(q1,q2,f);
                        qd(q1,q2,0,0,f) = w * (alpha / 2 * dot + gamma * abs);
                        qd(q1,q2,1,0,f) = w * (alpha / 2 * dot - gamma * abs);
                        qd(q1,q2,0,1,f) = w * (-alpha / 2 * dot - gamma * abs);
                        qd(q1,q2,1,1,f) = w * (-alpha / 2 * dot + gamma * abs);
                    }
                }
        });
}

static void PADGTraceSetup(const int dim,
    const int D1D,
    const int Q1D,
    const int NF,
    const Array<double>& W,
    const Vector& det,
    const Vector& nor,
    const Vector& rho,
    const Vector& u,
    const double alpha,
    const double beta,
    const double gamma,
    Vector& op)
{
    if (dim == 1) { MFEM_ABORT("dim==1 not supported in PADGTraceSetup"); }
    if (dim == 2)
    {
        PADGTraceSetup2D(Q1D, NF, W, det, nor, rho, u, alpha, beta, gamma, op);
    }
    if (dim == 3)
    {
        PADGTraceSetup3D(Q1D, NF, W, det, nor, rho, u, alpha, beta, gamma, op);
    }
}

// PA DGTrace Apply 2D kernel for Gauss-Lobatto/Bernstein
template<int T_D1D = 0, int T_Q1D = 0> static
    void PADGTraceApply2D(const int NF,
        const Array<double>& b,
        const Array<double>& bt,
        const Vector& op_,
        const Vector& x_,
        Vector& y_,
        const int d1d = 0,
        const int q1d = 0)
{
    const int VDIM = 1;
    const int D1D = T_D1D ? T_D1D : d1d;
    const int Q1D = T_Q1D ? T_Q1D : q1d;
    MFEM_VERIFY(D1D <= MAX_D1D, "");
    MFEM_VERIFY(Q1D <= MAX_Q1D, "");
    auto B = Reshape(b.Read(), Q1D, D1D);
    auto Bt = Reshape(bt.Read(), D1D, Q1D);
    auto op = Reshape(op_.Read(), Q1D, 2, 2, NF);
    auto x = Reshape(x_.Read(), D1D, VDIM, 2, NF);
    auto y = Reshape(y_.ReadWrite(), D1D, VDIM, 2, NF);

    MFEM_FORALL(f, NF,
        {
            const int VDIM = 1;
            const int D1D = T_D1D ? T_D1D : d1d;
            const int Q1D = T_Q1D ? T_Q1D : q1d;
            // the following variables are evaluated at compile time
            constexpr int max_D1D = T_D1D ? T_D1D : MAX_D1D;
            constexpr int max_Q1D = T_Q1D ? T_Q1D : MAX_Q1D;
            double u0[max_D1D][VDIM];
            double u1[max_D1D][VDIM];
            for (int d = 0; d < D1D; d++)
            {
                for (int c = 0; c < VDIM; c++)
                {
                    u0[d][c] = x(d,c,0,f);
                    u1[d][c] = x(d,c,1,f);
                }
            }
            double Bu0[max_Q1D][VDIM];
            double Bu1[max_Q1D][VDIM];
            for (int q = 0; q < Q1D; ++q)
            {
                for (int c = 0; c < VDIM; c++)
                {
                    Bu0[q][c] = 0.0;
                    Bu1[q][c] = 0.0;
                }
                for (int d = 0; d < D1D; ++d)
                {
                    const double b = B(q,d);
                    for (int c = 0; c < VDIM; c++)
                    {
                    Bu0[q][c] += b * u0[d][c];
                    Bu1[q][c] += b * u1[d][c];
                    }
                }
            }
            double DBu[max_Q1D][VDIM];
            for (int q = 0; q < Q1D; ++q)
            {
                for (int c = 0; c < VDIM; c++)
                {
                    DBu[q][c] = op(q,0,0,f) * Bu0[q][c] + op(q,1,0,f) * Bu1[q][c];
                }
            }
            double BDBu[max_D1D][VDIM];
            for (int d = 0; d < D1D; ++d)
            {
                for (int c = 0; c < VDIM; c++)
                {
                    BDBu[d][c] = 0.0;
                }
                for (int q = 0; q < Q1D; ++q)
                {
                    const double b = Bt(d,q);
                    for (int c = 0; c < VDIM; c++)
                    {
                    BDBu[d][c] += b * DBu[q][c];
                    }
                }
                for (int c = 0; c < VDIM; c++)
                {
                    y(d,c,0,f) += BDBu[d][c];
                    y(d,c,1,f) += -BDBu[d][c];
                }
            }
        });
}

// PA DGTrace Apply 3D kernel for Gauss-Lobatto/Bernstein
template<int T_D1D = 0, int T_Q1D = 0> static
    void PADGTraceApply3D(const int NF,
        const Array<double>& b,
        const Array<double>& bt,
        const Vector& op_,
        const Vector& x_,
        Vector& y_,
        const int d1d = 0,
        const int q1d = 0)
{
    const int VDIM = 1;
    const int D1D = T_D1D ? T_D1D : d1d;
    const int Q1D = T_Q1D ? T_Q1D : q1d;
    MFEM_VERIFY(D1D <= MAX_D1D, "");
    MFEM_VERIFY(Q1D <= MAX_Q1D, "");
    auto B = Reshape(b.Read(), Q1D, D1D);
    auto Bt = Reshape(bt.Read(), D1D, Q1D);
    auto op = Reshape(op_.Read(), Q1D, Q1D, 2, 2, NF);
    auto x = Reshape(x_.Read(), D1D, D1D, VDIM, 2, NF);
    auto y = Reshape(y_.ReadWrite(), D1D, D1D, VDIM, 2, NF);

    MFEM_FORALL(f, NF,
        {
            const int VDIM = 1;
            const int D1D = T_D1D ? T_D1D : d1d;
            const int Q1D = T_Q1D ? T_Q1D : q1d;
            // the following variables are evaluated at compile time
            constexpr int max_D1D = T_D1D ? T_D1D : MAX_D1D;
            constexpr int max_Q1D = T_Q1D ? T_Q1D : MAX_Q1D;
            double u0[max_D1D][max_D1D][VDIM];
            double u1[max_D1D][max_D1D][VDIM];
            for (int d1 = 0; d1 < D1D; d1++)
            {
                for (int d2 = 0; d2 < D1D; d2++)
                {
                    for (int c = 0; c < VDIM; c++)
                    {
                    u0[d1][d2][c] = x(d1,d2,c,0,f);
                    u1[d1][d2][c] = x(d1,d2,c,1,f);
                    }
                }
            }
            double Bu0[max_Q1D][max_D1D][VDIM];
            double Bu1[max_Q1D][max_D1D][VDIM];
            for (int q = 0; q < Q1D; ++q)
            {
                for (int d2 = 0; d2 < D1D; d2++)
                {
                    for (int c = 0; c < VDIM; c++)
                    {
                    Bu0[q][d2][c] = 0.0;
                    Bu1[q][d2][c] = 0.0;
                    }
                    for (int d1 = 0; d1 < D1D; ++d1)
                    {
                    const double b = B(q,d1);
                    for (int c = 0; c < VDIM; c++)
                    {
                        Bu0[q][d2][c] += b * u0[d1][d2][c];
                        Bu1[q][d2][c] += b * u1[d1][d2][c];
                    }
                    }
                }
            }
            double BBu0[max_Q1D][max_Q1D][VDIM];
            double BBu1[max_Q1D][max_Q1D][VDIM];
            for (int q1 = 0; q1 < Q1D; ++q1)
            {
                for (int q2 = 0; q2 < Q1D; q2++)
                {
                    for (int c = 0; c < VDIM; c++)
                    {
                    BBu0[q1][q2][c] = 0.0;
                    BBu1[q1][q2][c] = 0.0;
                    }
                    for (int d2 = 0; d2 < D1D; ++d2)
                    {
                    const double b = B(q2,d2);
                    for (int c = 0; c < VDIM; c++)
                    {
                        BBu0[q1][q2][c] += b * Bu0[q1][d2][c];
                        BBu1[q1][q2][c] += b * Bu1[q1][d2][c];
                    }
                    }
                }
            }
            double DBBu[max_Q1D][max_Q1D][VDIM];
            for (int q1 = 0; q1 < Q1D; ++q1)
            {
                for (int q2 = 0; q2 < Q1D; q2++)
                {
                    for (int c = 0; c < VDIM; c++)
                    {
                    DBBu[q1][q2][c] = op(q1,q2,0,0,f) * BBu0[q1][q2][c] +
                                        op(q1,q2,1,0,f) * BBu1[q1][q2][c];
                    }
                }
            }
            double BDBBu[max_Q1D][max_D1D][VDIM];
            for (int q1 = 0; q1 < Q1D; ++q1)
            {
                for (int d2 = 0; d2 < D1D; d2++)
                {
                    for (int c = 0; c < VDIM; c++)
                    {
                    BDBBu[q1][d2][c] = 0.0;
                    }
                    for (int q2 = 0; q2 < Q1D; ++q2)
                    {
                    const double b = Bt(d2,q2);
                    for (int c = 0; c < VDIM; c++)
                    {
                        BDBBu[q1][d2][c] += b * DBBu[q1][q2][c];
                    }
                    }
                }
            }
            double BBDBBu[max_D1D][max_D1D][VDIM];
            for (int d1 = 0; d1 < D1D; ++d1)
            {
                for (int d2 = 0; d2 < D1D; d2++)
                {
                    for (int c = 0; c < VDIM; c++)
                    {
                    BBDBBu[d1][d2][c] = 0.0;
                    }
                    for (int q1 = 0; q1 < Q1D; ++q1)
                    {
                    const double b = Bt(d1,q1);
                    for (int c = 0; c < VDIM; c++)
                    {
                        BBDBBu[d1][d2][c] += b * BDBBu[q1][d2][c];
                    }
                    }
                    for (int c = 0; c < VDIM; c++)
                    {
                    y(d1,d2,c,0,f) += BBDBBu[d1][d2][c];
                    y(d1,d2,c,1,f) += -BBDBBu[d1][d2][c];
                    }
                }
            }
        });
}

// Optimized PA DGTrace Apply 3D kernel for Gauss-Lobatto/Bernstein
template<int T_D1D = 0, int T_Q1D = 0, int T_NBZ = 0> static
    void SmemPADGTraceApply3D(const int NF,
        const Array<double>& b,
        const Array<double>& bt,
        const Vector& op_,
        const Vector& x_,
        Vector& y_,
        const int d1d = 0,
        const int q1d = 0)
{
    const int D1D = T_D1D ? T_D1D : d1d;
    const int Q1D = T_Q1D ? T_Q1D : q1d;
    constexpr int NBZ = T_NBZ ? T_NBZ : 1;
    MFEM_VERIFY(D1D <= MAX_D1D, "");
    MFEM_VERIFY(Q1D <= MAX_Q1D, "");
    auto B = Reshape(b.Read(), Q1D, D1D);
    auto Bt = Reshape(bt.Read(), D1D, Q1D);
    auto op = Reshape(op_.Read(), Q1D, Q1D, 2, 2, NF);
    auto x = Reshape(x_.Read(), D1D, D1D, 2, NF);
    auto y = Reshape(y_.ReadWrite(), D1D, D1D, 2, NF);

    MFEM_FORALL_2D(f, NF, Q1D, Q1D, NBZ,
        {
            const int tidz = MFEM_THREAD_ID(z);
            const int D1D = T_D1D ? T_D1D : d1d;
            const int Q1D = T_Q1D ? T_Q1D : q1d;
            // the following variables are evaluated at compile time
            constexpr int NBZ = T_NBZ ? T_NBZ : 1;
            constexpr int max_D1D = T_D1D ? T_D1D : MAX_D1D;
            constexpr int max_Q1D = T_Q1D ? T_Q1D : MAX_Q1D;
            MFEM_SHARED double u0[NBZ][max_D1D][max_D1D];
            MFEM_SHARED double u1[NBZ][max_D1D][max_D1D];
            MFEM_FOREACH_THREAD(d1,x,D1D)
            {
                MFEM_FOREACH_THREAD(d2,y,D1D)
                {
                    u0[tidz][d1][d2] = x(d1,d2,0,f + tidz);
                    u1[tidz][d1][d2] = x(d1,d2,1,f + tidz);
                }
            }
            MFEM_SYNC_THREAD;
            MFEM_SHARED double Bu0[NBZ][max_Q1D][max_D1D];
            MFEM_SHARED double Bu1[NBZ][max_Q1D][max_D1D];
            MFEM_FOREACH_THREAD(q1,x,Q1D)
            {
                MFEM_FOREACH_THREAD(d2,y,D1D)
                {
                    double Bu0_ = 0.0;
                    double Bu1_ = 0.0;
                    for (int d1 = 0; d1 < D1D; ++d1)
                    {
                    const double b = B(q1,d1);
                    Bu0_ += b * u0[tidz][d1][d2];
                    Bu1_ += b * u1[tidz][d1][d2];
                    }
                    Bu0[tidz][q1][d2] = Bu0_;
                    Bu1[tidz][q1][d2] = Bu1_;
                }
            }
            MFEM_SYNC_THREAD;
            MFEM_SHARED double BBu0[NBZ][max_Q1D][max_Q1D];
            MFEM_SHARED double BBu1[NBZ][max_Q1D][max_Q1D];
            MFEM_FOREACH_THREAD(q1,x,Q1D)
            {
                MFEM_FOREACH_THREAD(q2,y,Q1D)
                {
                    double BBu0_ = 0.0;
                    double BBu1_ = 0.0;
                    for (int d2 = 0; d2 < D1D; ++d2)
                    {
                    const double b = B(q2,d2);
                    BBu0_ += b * Bu0[tidz][q1][d2];
                    BBu1_ += b * Bu1[tidz][q1][d2];
                    }
                    BBu0[tidz][q1][q2] = BBu0_;
                    BBu1[tidz][q1][q2] = BBu1_;
                }
            }
            MFEM_SYNC_THREAD;
            MFEM_SHARED double DBBu[NBZ][max_Q1D][max_Q1D];
            MFEM_FOREACH_THREAD(q1,x,Q1D)
            {
                MFEM_FOREACH_THREAD(q2,y,Q1D)
                {
                    DBBu[tidz][q1][q2] = op(q1,q2,0,0,f + tidz) * BBu0[tidz][q1][q2] +
                                        op(q1,q2,1,0,f + tidz) * BBu1[tidz][q1][q2];
                }
            }
            MFEM_SYNC_THREAD;
            MFEM_SHARED double BDBBu[NBZ][max_Q1D][max_D1D];
            MFEM_FOREACH_THREAD(q1,x,Q1D)
            {
                MFEM_FOREACH_THREAD(d2,y,D1D)
                {
                    double BDBBu_ = 0.0;
                    for (int q2 = 0; q2 < Q1D; ++q2)
                    {
                    const double b = Bt(d2,q2);
                    BDBBu_ += b * DBBu[tidz][q1][q2];
                    }
                    BDBBu[tidz][q1][d2] = BDBBu_;
                }
            }
            MFEM_SYNC_THREAD;
            MFEM_FOREACH_THREAD(d1,x,D1D)
            {
                MFEM_FOREACH_THREAD(d2,y,D1D)
                {
                    double BBDBBu_ = 0.0;
                    for (int q1 = 0; q1 < Q1D; ++q1)
                    {
                    const double b = Bt(d1,q1);
                    BBDBBu_ += b * BDBBu[tidz][q1][d2];
                    }
                    y(d1,d2,0,f + tidz) += BBDBBu_;
                    y(d1,d2,1,f + tidz) += -BBDBBu_;
                }
            }
        });
}

static void PADGTraceApply(const int dim,
    const int D1D,
    const int Q1D,
    const int NF,
    const Array<double>& B,
    const Array<double>& Bt,
    const Vector& op,
    const Vector& x,
    Vector& y)
{
    if (dim == 2)
    {
        switch ((D1D << 4) | Q1D)
        {
        case 0x22: return PADGTraceApply2D<2, 2>(NF, B, Bt, op, x, y);
        case 0x33: return PADGTraceApply2D<3, 3>(NF, B, Bt, op, x, y);
        case 0x44: return PADGTraceApply2D<4, 4>(NF, B, Bt, op, x, y);
        case 0x55: return PADGTraceApply2D<5, 5>(NF, B, Bt, op, x, y);
        case 0x66: return PADGTraceApply2D<6, 6>(NF, B, Bt, op, x, y);
        case 0x77: return PADGTraceApply2D<7, 7>(NF, B, Bt, op, x, y);
        case 0x88: return PADGTraceApply2D<8, 8>(NF, B, Bt, op, x, y);
        case 0x99: return PADGTraceApply2D<9, 9>(NF, B, Bt, op, x, y);
        default:   return PADGTraceApply2D(NF, B, Bt, op, x, y, D1D, Q1D);
        }
    }
    else if (dim == 3)
    {
        switch ((D1D << 4) | Q1D)
        {
        case 0x23: return SmemPADGTraceApply3D<2, 3, 1>(NF, B, Bt, op, x, y);
        case 0x34: return SmemPADGTraceApply3D<3, 4, 2>(NF, B, Bt, op, x, y);
        case 0x45: return SmemPADGTraceApply3D<4, 5, 2>(NF, B, Bt, op, x, y);
        case 0x56: return SmemPADGTraceApply3D<5, 6, 1>(NF, B, Bt, op, x, y);
        case 0x67: return SmemPADGTraceApply3D<6, 7, 1>(NF, B, Bt, op, x, y);
        case 0x78: return SmemPADGTraceApply3D<7, 8, 1>(NF, B, Bt, op, x, y);
        case 0x89: return SmemPADGTraceApply3D<8, 9, 1>(NF, B, Bt, op, x, y);
        default:   return PADGTraceApply3D(NF, B, Bt, op, x, y, D1D, Q1D);
        }
    }
    MFEM_ABORT("Unknown kernel.");
}

// PA DGTrace Apply 2D kernel for Gauss-Lobatto/Bernstein
template<int T_D1D = 0, int T_Q1D = 0> static
    void PADGTraceApplyTranspose2D(const int NF,
        const Array<double>& b,
        const Array<double>& bt,
        const Vector& op_,
        const Vector& x_,
        Vector& y_,
        const int d1d = 0,
        const int q1d = 0)
{
    const int VDIM = 1;
    const int D1D = T_D1D ? T_D1D : d1d;
    const int Q1D = T_Q1D ? T_Q1D : q1d;
    MFEM_VERIFY(D1D <= MAX_D1D, "");
    MFEM_VERIFY(Q1D <= MAX_Q1D, "");
    auto B = Reshape(b.Read(), Q1D, D1D);
    auto Bt = Reshape(bt.Read(), D1D, Q1D);
    auto op = Reshape(op_.Read(), Q1D, 2, 2, NF);
    auto x = Reshape(x_.Read(), D1D, VDIM, 2, NF);
    auto y = Reshape(y_.ReadWrite(), D1D, VDIM, 2, NF);

    MFEM_FORALL(f, NF,
        {
            const int VDIM = 1;
            const int D1D = T_D1D ? T_D1D : d1d;
            const int Q1D = T_Q1D ? T_Q1D : q1d;
            // the following variables are evaluated at compile time
            constexpr int max_D1D = T_D1D ? T_D1D : MAX_D1D;
            constexpr int max_Q1D = T_Q1D ? T_Q1D : MAX_Q1D;
            double u0[max_D1D][VDIM];
            double u1[max_D1D][VDIM];
            for (int d = 0; d < D1D; d++)
            {
                for (int c = 0; c < VDIM; c++)
                {
                    u0[d][c] = x(d,c,0,f);
                    u1[d][c] = x(d,c,1,f);
                }
            }
            double Bu0[max_Q1D][VDIM];
            double Bu1[max_Q1D][VDIM];
            for (int q = 0; q < Q1D; ++q)
            {
                for (int c = 0; c < VDIM; c++)
                {
                    Bu0[q][c] = 0.0;
                    Bu1[q][c] = 0.0;
                }
                for (int d = 0; d < D1D; ++d)
                {
                    const double b = B(q,d);
                    for (int c = 0; c < VDIM; c++)
                    {
                    Bu0[q][c] += b * u0[d][c];
                    Bu1[q][c] += b * u1[d][c];
                    }
                }
            }
            double DBu0[max_Q1D][VDIM];
            double DBu1[max_Q1D][VDIM];
            for (int q = 0; q < Q1D; ++q)
            {
                for (int c = 0; c < VDIM; c++)
                {
                    DBu0[q][c] = op(q,0,0,f) * Bu0[q][c] + op(q,0,1,f) * Bu1[q][c];
                    DBu1[q][c] = op(q,1,0,f) * Bu0[q][c] + op(q,1,1,f) * Bu1[q][c];
                }
            }
            double BDBu0[max_D1D][VDIM];
            double BDBu1[max_D1D][VDIM];
            for (int d = 0; d < D1D; ++d)
            {
                for (int c = 0; c < VDIM; c++)
                {
                    BDBu0[d][c] = 0.0;
                    BDBu1[d][c] = 0.0;
                }
                for (int q = 0; q < Q1D; ++q)
                {
                    const double b = Bt(d,q);
                    for (int c = 0; c < VDIM; c++)
                    {
                    BDBu0[d][c] += b * DBu0[q][c];
                    BDBu1[d][c] += b * DBu1[q][c];
                    }
                }
                for (int c = 0; c < VDIM; c++)
                {
                    y(d,c,0,f) += BDBu0[d][c];
                    y(d,c,1,f) += BDBu1[d][c];
                }
            }
        });
}

// PA DGTrace Apply Transpose 3D kernel for Gauss-Lobatto/Bernstein
template<int T_D1D = 0, int T_Q1D = 0> static
    void PADGTraceApplyTranspose3D(const int NF,
        const Array<double>& b,
        const Array<double>& bt,
        const Vector& op_,
        const Vector& x_,
        Vector& y_,
        const int d1d = 0,
        const int q1d = 0)
{
    const int VDIM = 1;
    const int D1D = T_D1D ? T_D1D : d1d;
    const int Q1D = T_Q1D ? T_Q1D : q1d;
    MFEM_VERIFY(D1D <= MAX_D1D, "");
    MFEM_VERIFY(Q1D <= MAX_Q1D, "");
    auto B = Reshape(b.Read(), Q1D, D1D);
    auto Bt = Reshape(bt.Read(), D1D, Q1D);
    auto op = Reshape(op_.Read(), Q1D, Q1D, 2, 2, NF);
    auto x = Reshape(x_.Read(), D1D, D1D, VDIM, 2, NF);
    auto y = Reshape(y_.ReadWrite(), D1D, D1D, VDIM, 2, NF);

    MFEM_FORALL(f, NF,
        {
            const int VDIM = 1;
            const int D1D = T_D1D ? T_D1D : d1d;
            const int Q1D = T_Q1D ? T_Q1D : q1d;
            // the following variables are evaluated at compile time
            constexpr int max_D1D = T_D1D ? T_D1D : MAX_D1D;
            constexpr int max_Q1D = T_Q1D ? T_Q1D : MAX_Q1D;
            double u0[max_D1D][max_D1D][VDIM];
            double u1[max_D1D][max_D1D][VDIM];
            for (int d1 = 0; d1 < D1D; d1++)
            {
                for (int d2 = 0; d2 < D1D; d2++)
                {
                    for (int c = 0; c < VDIM; c++)
                    {
                    u0[d1][d2][c] = x(d1,d2,c,0,f);
                    u1[d1][d2][c] = x(d1,d2,c,1,f);
                    }
                }
            }
            double Bu0[max_Q1D][max_D1D][VDIM];
            double Bu1[max_Q1D][max_D1D][VDIM];
            for (int q1 = 0; q1 < Q1D; ++q1)
            {
                for (int d2 = 0; d2 < D1D; ++d2)
                {
                    for (int c = 0; c < VDIM; c++)
                    {
                    Bu0[q1][d2][c] = 0.0;
                    Bu1[q1][d2][c] = 0.0;
                    }
                    for (int d1 = 0; d1 < D1D; ++d1)
                    {
                    const double b = B(q1,d1);
                    for (int c = 0; c < VDIM; c++)
                    {
                        Bu0[q1][d2][c] += b * u0[d1][d2][c];
                        Bu1[q1][d2][c] += b * u1[d1][d2][c];
                    }
                    }
                }
            }
            double BBu0[max_Q1D][max_Q1D][VDIM];
            double BBu1[max_Q1D][max_Q1D][VDIM];
            for (int q1 = 0; q1 < Q1D; ++q1)
            {
                for (int q2 = 0; q2 < Q1D; ++q2)
                {
                    for (int c = 0; c < VDIM; c++)
                    {
                    BBu0[q1][q2][c] = 0.0;
                    BBu1[q1][q2][c] = 0.0;
                    }
                    for (int d2 = 0; d2 < D1D; ++d2)
                    {
                    const double b = B(q2,d2);
                    for (int c = 0; c < VDIM; c++)
                    {
                        BBu0[q1][q2][c] += b * Bu0[q1][d2][c];
                        BBu1[q1][q2][c] += b * Bu1[q1][d2][c];
                    }
                    }
                }
            }
            double DBu0[max_Q1D][max_Q1D][VDIM];
            double DBu1[max_Q1D][max_Q1D][VDIM];
            for (int q1 = 0; q1 < Q1D; ++q1)
            {
                for (int q2 = 0; q2 < Q1D; ++q2)
                {
                    const double D00 = op(q1,q2,0,0,f);
                    const double D01 = op(q1,q2,0,1,f);
                    const double D10 = op(q1,q2,1,0,f);
                    const double D11 = op(q1,q2,1,1,f);
                    for (int c = 0; c < VDIM; c++)
                    {
                    DBu0[q1][q2][c] = D00 * BBu0[q1][q2][c] + D01 * BBu1[q1][q2][c];
                    DBu1[q1][q2][c] = D10 * BBu0[q1][q2][c] + D11 * BBu1[q1][q2][c];
                    }
                }
            }
            double BDBu0[max_D1D][max_Q1D][VDIM];
            double BDBu1[max_D1D][max_Q1D][VDIM];
            for (int d1 = 0; d1 < D1D; ++d1)
            {
                for (int q2 = 0; q2 < Q1D; ++q2)
                {
                    for (int c = 0; c < VDIM; c++)
                    {
                    BDBu0[d1][q2][c] = 0.0;
                    BDBu1[d1][q2][c] = 0.0;
                    }
                    for (int q1 = 0; q1 < Q1D; ++q1)
                    {
                    const double b = Bt(d1,q1);
                    for (int c = 0; c < VDIM; c++)
                    {
                        BDBu0[d1][q2][c] += b * DBu0[q1][q2][c];
                        BDBu1[d1][q2][c] += b * DBu1[q1][q2][c];
                    }
                    }
                }
            }
            double BBDBu0[max_D1D][max_D1D][VDIM];
            double BBDBu1[max_D1D][max_D1D][VDIM];
            for (int d1 = 0; d1 < D1D; ++d1)
            {
                for (int d2 = 0; d2 < D1D; ++d2)
                {
                    for (int c = 0; c < VDIM; c++)
                    {
                    BBDBu0[d1][d2][c] = 0.0;
                    BBDBu1[d1][d2][c] = 0.0;
                    }
                    for (int q2 = 0; q2 < Q1D; ++q2)
                    {
                    const double b = Bt(d2,q2);
                    for (int c = 0; c < VDIM; c++)
                    {
                        BBDBu0[d1][d2][c] += b * BDBu0[d1][q2][c];
                        BBDBu1[d1][d2][c] += b * BDBu1[d1][q2][c];
                    }
                    }
                    for (int c = 0; c < VDIM; c++)
                    {
                    y(d1,d2,c,0,f) += BBDBu0[d1][d2][c];
                    y(d1,d2,c,1,f) += BBDBu1[d1][d2][c];
                    }
                }
            }
        });
}

// Optimized PA DGTrace Apply Transpose 3D kernel for Gauss-Lobatto/Bernstein
template<int T_D1D = 0, int T_Q1D = 0, int T_NBZ = 0> static
    void SmemPADGTraceApplyTranspose3D(const int NF,
        const Array<double>& b,
        const Array<double>& bt,
        const Vector& op_,
        const Vector& x_,
        Vector& y_,
        const int d1d = 0,
        const int q1d = 0)
{
    const int D1D = T_D1D ? T_D1D : d1d;
    const int Q1D = T_Q1D ? T_Q1D : q1d;
    constexpr int NBZ = T_NBZ ? T_NBZ : 1;
    MFEM_VERIFY(D1D <= MAX_D1D, "");
    MFEM_VERIFY(Q1D <= MAX_Q1D, "");
    auto B = Reshape(b.Read(), Q1D, D1D);
    auto Bt = Reshape(bt.Read(), D1D, Q1D);
    auto op = Reshape(op_.Read(), Q1D, Q1D, 2, 2, NF);
    auto x = Reshape(x_.Read(), D1D, D1D, 2, NF);
    auto y = Reshape(y_.ReadWrite(), D1D, D1D, 2, NF);

    MFEM_FORALL_2D(f, NF, Q1D, Q1D, NBZ,
        {
            const int tidz = MFEM_THREAD_ID(z);
            const int D1D = T_D1D ? T_D1D : d1d;
            const int Q1D = T_Q1D ? T_Q1D : q1d;
            // the following variables are evaluated at compile time
            constexpr int NBZ = T_NBZ ? T_NBZ : 1;
            constexpr int max_D1D = T_D1D ? T_D1D : MAX_D1D;
            constexpr int max_Q1D = T_Q1D ? T_Q1D : MAX_Q1D;
            MFEM_SHARED double u0[NBZ][max_D1D][max_D1D];
            MFEM_SHARED double u1[NBZ][max_D1D][max_D1D];
            MFEM_FOREACH_THREAD(d1,x,D1D)
            {
                MFEM_FOREACH_THREAD(d2,y,D1D)
                {
                    u0[tidz][d1][d2] = x(d1,d2,0,f + tidz);
                    u1[tidz][d1][d2] = x(d1,d2,1,f + tidz);
                }
            }
            MFEM_SYNC_THREAD;
            MFEM_SHARED double Bu0[NBZ][max_Q1D][max_D1D];
            MFEM_SHARED double Bu1[NBZ][max_Q1D][max_D1D];
            MFEM_FOREACH_THREAD(q1,x,Q1D)
            {
                MFEM_FOREACH_THREAD(d2,y,D1D)
                {
                    double Bu0_ = 0.0;
                    double Bu1_ = 0.0;
                    for (int d1 = 0; d1 < D1D; ++d1)
                    {
                    const double b = B(q1,d1);
                    Bu0_ += b * u0[tidz][d1][d2];
                    Bu1_ += b * u1[tidz][d1][d2];
                    }
                    Bu0[tidz][q1][d2] = Bu0_;
                    Bu1[tidz][q1][d2] = Bu1_;
                }
            }
            MFEM_SYNC_THREAD;
            MFEM_SHARED double BBu0[NBZ][max_Q1D][max_Q1D];
            MFEM_SHARED double BBu1[NBZ][max_Q1D][max_Q1D];
            MFEM_FOREACH_THREAD(q1,x,Q1D)
            {
                MFEM_FOREACH_THREAD(q2,y,Q1D)
                {
                    double BBu0_ = 0.0;
                    double BBu1_ = 0.0;
                    for (int d2 = 0; d2 < D1D; ++d2)
                    {
                    const double b = B(q2,d2);
                    BBu0_ += b * Bu0[tidz][q1][d2];
                    BBu1_ += b * Bu1[tidz][q1][d2];
                    }
                    BBu0[tidz][q1][q2] = BBu0_;
                    BBu1[tidz][q1][q2] = BBu1_;
                }
            }
            MFEM_SYNC_THREAD;
            MFEM_SHARED double DBBu0[NBZ][max_Q1D][max_Q1D];
            MFEM_SHARED double DBBu1[NBZ][max_Q1D][max_Q1D];
            MFEM_FOREACH_THREAD(q1,x,Q1D)
            {
                MFEM_FOREACH_THREAD(q2,y,Q1D)
                {
                    const double D00 = op(q1,q2,0,0,f + tidz);
                    const double D01 = op(q1,q2,0,1,f + tidz);
                    const double D10 = op(q1,q2,1,0,f + tidz);
                    const double D11 = op(q1,q2,1,1,f + tidz);
                    const double u0 = BBu0[tidz][q1][q2];
                    const double u1 = BBu1[tidz][q1][q2];
                    DBBu0[tidz][q1][q2] = D00 * u0 + D01 * u1;
                    DBBu1[tidz][q1][q2] = D10 * u0 + D11 * u1;
                }
            }
            MFEM_SYNC_THREAD;
            MFEM_SHARED double BDBBu0[NBZ][max_Q1D][max_D1D];
            MFEM_SHARED double BDBBu1[NBZ][max_Q1D][max_D1D];
            MFEM_FOREACH_THREAD(q1,x,Q1D)
            {
                MFEM_FOREACH_THREAD(d2,y,D1D)
                {
                    double BDBBu0_ = 0.0;
                    double BDBBu1_ = 0.0;
                    for (int q2 = 0; q2 < Q1D; ++q2)
                    {
                    const double b = Bt(d2,q2);
                    BDBBu0_ += b * DBBu0[tidz][q1][q2];
                    BDBBu1_ += b * DBBu1[tidz][q1][q2];
                    }
                    BDBBu0[tidz][q1][d2] = BDBBu0_;
                    BDBBu1[tidz][q1][d2] = BDBBu1_;
                }
            }
            MFEM_SYNC_THREAD;
            MFEM_FOREACH_THREAD(d1,x,D1D)
            {
                MFEM_FOREACH_THREAD(d2,y,D1D)
                {
                    double BBDBBu0_ = 0.0;
                    double BBDBBu1_ = 0.0;
                    for (int q1 = 0; q1 < Q1D; ++q1)
                    {
                    const double b = Bt(d1,q1);
                    BBDBBu0_ += b * BDBBu0[tidz][q1][d2];
                    BBDBBu1_ += b * BDBBu1[tidz][q1][d2];
                    }
                    y(d1,d2,0,f + tidz) += BBDBBu0_;
                    y(d1,d2,1,f + tidz) += BBDBBu1_;
                }
            }
        });
}

static void PADGTraceApplyTranspose(const int dim,
    const int D1D,
    const int Q1D,
    const int NF,
    const Array<double>& B,
    const Array<double>& Bt,
    const Vector& op,
    const Vector& x,
    Vector& y)
{
    if (dim == 2)
    {
        switch ((D1D << 4) | Q1D)
        {
        case 0x22: return PADGTraceApplyTranspose2D<2, 2>(NF, B, Bt, op, x, y);
        case 0x33: return PADGTraceApplyTranspose2D<3, 3>(NF, B, Bt, op, x, y);
        case 0x44: return PADGTraceApplyTranspose2D<4, 4>(NF, B, Bt, op, x, y);
        case 0x55: return PADGTraceApplyTranspose2D<5, 5>(NF, B, Bt, op, x, y);
        case 0x66: return PADGTraceApplyTranspose2D<6, 6>(NF, B, Bt, op, x, y);
        case 0x77: return PADGTraceApplyTranspose2D<7, 7>(NF, B, Bt, op, x, y);
        case 0x88: return PADGTraceApplyTranspose2D<8, 8>(NF, B, Bt, op, x, y);
        case 0x99: return PADGTraceApplyTranspose2D<9, 9>(NF, B, Bt, op, x, y);
        default: return PADGTraceApplyTranspose2D(NF, B, Bt, op, x, y, D1D, Q1D);
        }
    }
    else if (dim == 3)
    {
        switch ((D1D << 4) | Q1D)
        {
        case 0x23: return SmemPADGTraceApplyTranspose3D<2, 3>(NF, B, Bt, op, x, y);
        case 0x34: return SmemPADGTraceApplyTranspose3D<3, 4>(NF, B, Bt, op, x, y);
        case 0x45: return SmemPADGTraceApplyTranspose3D<4, 5>(NF, B, Bt, op, x, y);
        case 0x56: return SmemPADGTraceApplyTranspose3D<5, 6>(NF, B, Bt, op, x, y);
        case 0x67: return SmemPADGTraceApplyTranspose3D<6, 7>(NF, B, Bt, op, x, y);
        case 0x78: return SmemPADGTraceApplyTranspose3D<7, 8>(NF, B, Bt, op, x, y);
        case 0x89: return SmemPADGTraceApplyTranspose3D<8, 9>(NF, B, Bt, op, x, y);
        default: return PADGTraceApplyTranspose3D(NF, B, Bt, op, x, y, D1D, Q1D);
        }
    }
    MFEM_ABORT("Unknown kernel.");
}

/*########################## MDG START  ##########################*/
//Has alpha (Done?)
void MaxwellDGTraceIntegrator::AssembleFaceMatrix(const FiniteElement& el1,
    const FiniteElement& el2,
    FaceElementTransformations& Trans,
    DenseMatrix& elmat)
{
    int dim, ndof1, ndof2;

    double un, a, b, g, w;

    dim = el1.GetDim();
    ndof1 = el1.GetDof();
    Vector vu(dim), nor(dim);

    if (Trans.Elem2No >= 0)
    {
        ndof2 = el2.GetDof();
    }
    else
    {
        ndof2 = 0;
    }

    shape1_.SetSize(ndof1);
    shape2_.SetSize(ndof2);
    elmat.SetSize(ndof1 + ndof2);
    elmat = 0.0;

    const IntegrationRule* ir = IntRule;
    if (ir == NULL)
    {
        int order;
        // Assuming order(u)==order(mesh)
        if (Trans.Elem2No >= 0)
            order = (std::min(Trans.Elem1->OrderW(), Trans.Elem2->OrderW()) +
                2 * std::max(el1.GetOrder(), el2.GetOrder()));
        else
        {
            order = Trans.Elem1->OrderW() + 2 * el1.GetOrder();
        }
        if (el1.Space() == FunctionSpace::Pk)
        {
            order++;
        }
        ir = &IntRules.Get(Trans.GetGeometryType(), order);
    }

    for (int p = 0; p < ir->GetNPoints(); p++)
    {
        const IntegrationPoint& ip = ir->IntPoint(p);

        // Set the integration point in the face and the neighboring elements
        Trans.SetAllIntPoints(&ip);

        // Access the neighboring elements' integration points
        // Note: eip2 will only contain valid data if Elem2 exists
        const IntegrationPoint& eip1 = Trans.GetElement1IntPoint();
        const IntegrationPoint& eip2 = Trans.GetElement2IntPoint();

        el1.CalcShape(eip1, shape1_);

        u->Eval(vu, *Trans.Elem1, eip1);

        if (dim == 1)
        {
            nor(0) = 2 * eip1.x - 1.0;
        }
        else
        {
            CalcOrtho(Trans.Jacobian(), nor);
        }

        un = vu * nor;
        a = 0.5 * alpha * un;
        //b = beta * fabs(un);
        g = gamma * fabs(un);
        // note: if |alpha/2|==|beta| then |a|==|b|, i.e. (a==b) or (a==-b)
        //       and therefore two blocks in the element matrix contribution
        //       (from the current quadrature point) are 0

        if (rho)
        {
            double rho_p;
            if (un >= 0.0 && ndof2)
            {
                rho_p = rho->Eval(*Trans.Elem2, eip2);
            }
            else
            {
                rho_p = rho->Eval(*Trans.Elem1, eip1);
            }
            a *= rho_p;
            //b *= rho_p;
            g *= rho_p;
        }

        w = ip.weight * (a + g);
        if (w != 0.0)
        {
            for (int i = 0; i < ndof1; i++)
                for (int j = 0; j < ndof1; j++)
                {
                    elmat(i, j) += w * shape1_(i) * shape1_(j);
                }
        }

        if (ndof2)
        {
            el2.CalcShape(eip2, shape2_);

            if (w != 0.0)
                for (int i = 0; i < ndof2; i++)
                    for (int j = 0; j < ndof1; j++)
                    {
                        elmat(ndof1 + i, j) -= w * shape2_(i) * shape1_(j);
                    }

            w = ip.weight * (g - a);
            if (w != 0.0)
            {
                for (int i = 0; i < ndof2; i++)
                    for (int j = 0; j < ndof2; j++)
                    {
                        elmat(ndof1 + i, ndof1 + j) += w * shape2_(i) * shape2_(j);
                    }

                for (int i = 0; i < ndof1; i++)
                    for (int j = 0; j < ndof2; j++)
                    {
                        elmat(i, ndof1 + j) -= w * shape1_(i) * shape2_(j);
                    }
            }
        }
    }
}

//Has alpha through SetupPA
void MaxwellDGTraceIntegrator::AssemblePAInteriorFaces(const FiniteElementSpace& fes)
{
	SetupPA(fes, FaceType::Interior);
}

//Has alpha through SetupPA
void MaxwellDGTraceIntegrator::AssemblePABoundaryFaces(const FiniteElementSpace& fes)
{
	SetupPA(fes, FaceType::Boundary);
}

//Does not have alpha
void MaxwellDGTraceIntegrator::AddMultTransposePA(const Vector& x, Vector& y) const
{
    PADGTraceApplyTranspose(dim, dofs1D, quad1D, nf,
        maps->B, maps->Bt,
        pa_data, x, y);
}

//Does not have alpha
void MaxwellDGTraceIntegrator::AddMultPA(const Vector& x, Vector& y) const
{
    PADGTraceApply(dim, dofs1D, quad1D, nf,
        maps->B, maps->Bt,
        pa_data, x, y);
}

//Does not have alpha
void MaxwellDGTraceIntegrator::AssembleEAInteriorFaces(const FiniteElementSpace& fes,
    Vector& ea_data_int,
    Vector& ea_data_ext,
    const bool add)
{
    SetupPA(fes, FaceType::Interior);
    nf = fes.GetNFbyType(FaceType::Interior);
    if (nf == 0) { return; }
    const Array<double>& B = maps->B;
    if (dim == 1)
    {
        return EADGTraceAssemble1DInt(nf, B, pa_data, ea_data_int, ea_data_ext, add);
    }
    else if (dim == 2)
    {
        switch ((dofs1D << 4) | quad1D)
        {
        case 0x22:
            return EADGTraceAssemble2DInt<2, 2>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        case 0x33:
            return EADGTraceAssemble2DInt<3, 3>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        case 0x44:
            return EADGTraceAssemble2DInt<4, 4>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        case 0x55:
            return EADGTraceAssemble2DInt<5, 5>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        case 0x66:
            return EADGTraceAssemble2DInt<6, 6>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        case 0x77:
            return EADGTraceAssemble2DInt<7, 7>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        case 0x88:
            return EADGTraceAssemble2DInt<8, 8>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        case 0x99:
            return EADGTraceAssemble2DInt<9, 9>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        default:
            return EADGTraceAssemble2DInt(nf, B, pa_data, ea_data_int,
                ea_data_ext, add, dofs1D, quad1D);
        }
    }
    else if (dim == 3)
    {
        switch ((dofs1D << 4) | quad1D)
        {
        case 0x23:
            return EADGTraceAssemble3DInt<2, 3>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        case 0x34:
            return EADGTraceAssemble3DInt<3, 4>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        case 0x45:
            return EADGTraceAssemble3DInt<4, 5>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        case 0x56:
            return EADGTraceAssemble3DInt<5, 6>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        case 0x67:
            return EADGTraceAssemble3DInt<6, 7>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        case 0x78:
            return EADGTraceAssemble3DInt<7, 8>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        case 0x89:
            return EADGTraceAssemble3DInt<8, 9>(nf, B, pa_data, ea_data_int,
                ea_data_ext, add);
        default:
            return EADGTraceAssemble3DInt(nf, B, pa_data, ea_data_int,
                ea_data_ext, add, dofs1D, quad1D);
        }
    }
    MFEM_ABORT("Unknown kernel.");
}

//Does not have alpha
void MaxwellDGTraceIntegrator::AssembleEABoundaryFaces(const FiniteElementSpace& fes,
    Vector& ea_data_bdr,
    const bool add)
{
    SetupPA(fes, FaceType::Boundary);
    nf = fes.GetNFbyType(FaceType::Boundary);
    if (nf == 0) { return; }
    const Array<double>& B = maps->B;
    if (dim == 1)
    {
        return EADGTraceAssemble1DBdr(nf, B, pa_data, ea_data_bdr, add);
    }
    else if (dim == 2)
    {
        switch ((dofs1D << 4) | quad1D)
        {
        case 0x22: return EADGTraceAssemble2DBdr<2, 2>(nf, B, pa_data, ea_data_bdr, add);
        case 0x33: return EADGTraceAssemble2DBdr<3, 3>(nf, B, pa_data, ea_data_bdr, add);
        case 0x44: return EADGTraceAssemble2DBdr<4, 4>(nf, B, pa_data, ea_data_bdr, add);
        case 0x55: return EADGTraceAssemble2DBdr<5, 5>(nf, B, pa_data, ea_data_bdr, add);
        case 0x66: return EADGTraceAssemble2DBdr<6, 6>(nf, B, pa_data, ea_data_bdr, add);
        case 0x77: return EADGTraceAssemble2DBdr<7, 7>(nf, B, pa_data, ea_data_bdr, add);
        case 0x88: return EADGTraceAssemble2DBdr<8, 8>(nf, B, pa_data, ea_data_bdr, add);
        case 0x99: return EADGTraceAssemble2DBdr<9, 9>(nf, B, pa_data, ea_data_bdr, add);
        default:
            return EADGTraceAssemble2DBdr(nf, B, pa_data, ea_data_bdr, add, dofs1D, quad1D);
        }
    }
    else if (dim == 3)
    {
        switch ((dofs1D << 4) | quad1D)
        {
        case 0x23: return EADGTraceAssemble3DBdr<2, 3>(nf, B, pa_data, ea_data_bdr, add);
        case 0x34: return EADGTraceAssemble3DBdr<3, 4>(nf, B, pa_data, ea_data_bdr, add);
        case 0x45: return EADGTraceAssemble3DBdr<4, 5>(nf, B, pa_data, ea_data_bdr, add);
        case 0x56: return EADGTraceAssemble3DBdr<5, 6>(nf, B, pa_data, ea_data_bdr, add);
        case 0x67: return EADGTraceAssemble3DBdr<6, 7>(nf, B, pa_data, ea_data_bdr, add);
        case 0x78: return EADGTraceAssemble3DBdr<7, 8>(nf, B, pa_data, ea_data_bdr, add);
        case 0x89: return EADGTraceAssemble3DBdr<8, 9>(nf, B, pa_data, ea_data_bdr, add);
        default:
            return EADGTraceAssemble3DBdr(nf, B, pa_data, ea_data_bdr, add, dofs1D, quad1D);
        }
    }
    MFEM_ABORT("Unknown kernel.");
}

//Has alpha (Done?)
void MaxwellDGTraceIntegrator::SetupPA(const FiniteElementSpace& fes, FaceType type)
{
    nf = fes.GetNFbyType(type);
    if (nf == 0) { return; }
    // Assumes tensor-product elements
    Mesh* mesh = fes.GetMesh();
    const FiniteElement& el =
        *fes.GetTraceElement(0, fes.GetMesh()->GetFaceBaseGeometry(0));
    FaceElementTransformations& T =
        *fes.GetMesh()->GetFaceElementTransformations(0);
    const IntegrationRule* ir = IntRule ?
        IntRule :
        &GetRule(el.GetGeomType(), el.GetOrder(), T);
    const int symmDims = 4;
    const int nq = ir->GetNPoints();
    dim = mesh->Dimension();
    geom = mesh->GetFaceGeometricFactors(
        *ir,
        FaceGeometricFactors::DETERMINANTS |
        FaceGeometricFactors::NORMALS, type);
    maps = &el.GetDofToQuad(*ir, DofToQuad::TENSOR);
    dofs1D = maps->ndof;
    quad1D = maps->nqpt;
    pa_data.SetSize(symmDims * nq * nf, Device::GetDeviceMemoryType());
    Vector vel;
    if (VectorConstantCoefficient* c_u = dynamic_cast<VectorConstantCoefficient*>
        (u))
    {
        vel = c_u->GetVec();
    }
    else if (VectorQuadratureFunctionCoefficient* c_u =
        dynamic_cast<VectorQuadratureFunctionCoefficient*>(u))
    {
        // Assumed to be in lexicographical ordering
        const QuadratureFunction& qFun = c_u->GetQuadFunction();
        MFEM_VERIFY(qFun.Size() == dim * nq * nf,
            "Incompatible QuadratureFunction dimension \n");

        MFEM_VERIFY(ir == &qFun.GetSpace()->GetElementIntRule(0),
            "IntegrationRule used within integrator and in"
            " QuadratureFunction appear to be different");
        qFun.Read();
        vel.MakeRef(const_cast<QuadratureFunction&>(qFun), 0);
    }
    else
    {
        vel.SetSize(dim * nq * nf);
        auto C = Reshape(vel.HostWrite(), dim, nq, nf);
        Vector Vq(dim);
        int f_ind = 0;
        for (int f = 0; f < fes.GetNF(); ++f)
        {
            int e1, e2;
            int inf1, inf2;
            fes.GetMesh()->GetFaceElements(f, &e1, &e2);
            fes.GetMesh()->GetFaceInfos(f, &inf1, &inf2);
            int face_id = inf1 / 64;
            if ((type == FaceType::Interior && (e2 >= 0 || (e2 < 0 && inf2 >= 0))) ||
                (type == FaceType::Boundary && e2 < 0 && inf2 < 0))
            {
                FaceElementTransformations& T =
                    *fes.GetMesh()->GetFaceElementTransformations(f);
                for (int q = 0; q < nq; ++q)
                {
                    // Convert to lexicographic ordering
                    int iq = ToLexOrdering(dim, face_id, quad1D, q);
                    T.SetAllIntPoints(&ir->IntPoint(q));
                    const IntegrationPoint& eip1 = T.GetElement1IntPoint();
                    u->Eval(Vq, *T.Elem1, eip1);
                    for (int i = 0; i < dim; ++i)
                    {
                        C(i, iq, f_ind) = Vq(i);
                    }
                }
                f_ind++;
            }
        }
        MFEM_VERIFY(f_ind == nf, "Incorrect number of faces.");
    }
    Vector r;
    if (rho == nullptr)
    {
        r.SetSize(1);
        r(0) = 1.0;
    }
    else if (ConstantCoefficient* c_rho = dynamic_cast<ConstantCoefficient*>(rho))
    {
        r.SetSize(1);
        r(0) = c_rho->constant;
    }
    else if (QuadratureFunctionCoefficient* c_rho =
        dynamic_cast<QuadratureFunctionCoefficient*>(rho))
    {
        const QuadratureFunction& qFun = c_rho->GetQuadFunction();
        MFEM_VERIFY(qFun.Size() == nq * nf,
            "Incompatible QuadratureFunction dimension \n");

        MFEM_VERIFY(ir == &qFun.GetSpace()->GetElementIntRule(0),
            "IntegrationRule used within integrator and in"
            " QuadratureFunction appear to be different");
        qFun.Read();
        r.MakeRef(const_cast<QuadratureFunction&>(qFun), 0);
    }
    else
    {
        r.SetSize(nq * nf);
        auto C_vel = Reshape(vel.HostRead(), dim, nq, nf);
        auto n = Reshape(geom->normal.HostRead(), nq, dim, nf);
        auto C = Reshape(r.HostWrite(), nq, nf);
        int f_ind = 0;
        for (int f = 0; f < fes.GetNF(); ++f)
        {
            int e1, e2;
            int inf1, inf2;
            fes.GetMesh()->GetFaceElements(f, &e1, &e2);
            fes.GetMesh()->GetFaceInfos(f, &inf1, &inf2);
            int face_id = inf1 / 64;
            if ((type == FaceType::Interior && (e2 >= 0 || (e2 < 0 && inf2 >= 0))) ||
                (type == FaceType::Boundary && e2 < 0 && inf2 < 0))
            {
                FaceElementTransformations& T =
                    *fes.GetMesh()->GetFaceElementTransformations(f);
                for (int q = 0; q < nq; ++q)
                {
                    // Convert to lexicographic ordering
                    int iq = ToLexOrdering(dim, face_id, quad1D, q);

                    T.SetAllIntPoints(&ir->IntPoint(q));
                    const IntegrationPoint& eip1 = T.GetElement1IntPoint();
                    const IntegrationPoint& eip2 = T.GetElement2IntPoint();
                    double r;

                    if (inf2 < 0)
                    {
                        r = rho->Eval(*T.Elem1, eip1);
                    }
                    else
                    {
                        double udotn = 0.0;
                        for (int d = 0; d < dim; ++d)
                        {
                            udotn += C_vel(d, iq, f_ind) * n(iq, d, f_ind);
                        }
                        if (udotn >= 0.0) { r = rho->Eval(*T.Elem2, eip2); }
                        else { r = rho->Eval(*T.Elem1, eip1); }
                    }
                    C(iq, f_ind) = r;
                }
                f_ind++;
            }
        }
        MFEM_VERIFY(f_ind == nf, "Incorrect number of faces.");
    }
    PADGTraceSetup(dim, dofs1D, quad1D, nf, ir->GetWeights(),
    geom->detJ, geom->normal, r, vel,
    alpha, beta, gamma, pa_data);
}

//Does not have alpha
const IntegrationRule& MaxwellDGTraceIntegrator::GetRule(
    Geometry::Type geom, int order, FaceElementTransformations& T)
{
    int int_order = T.Elem1->OrderW() + 2 * order;
    return IntRules.Get(geom, int_order);
}

}


