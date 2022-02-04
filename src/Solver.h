#pragma once

#include "mfem.hpp"

namespace Maxwell {

typedef std::size_t Direction;

const Direction X = 0;
const Direction Y = 1;

//mfem::ODESolver* ode_solver = new mfem::RK4Solver;

class Solver {
public:
    typedef double ElectricField;
    typedef mfem::Vector Position;

    struct Options {
        int order = 3;
        double t_final = 10.0;
        double dt = 0.01;
        int vis_steps = 5;
        int precision = 8;
    };

	Solver(const Options&, const mfem::Mesh&);
    
    void setInitialElectricField(std::function<ElectricField(const Position&)>);
    
    mfem::Mesh& getMesh() { return mesh_;  }
    
    ElectricField getElectricFieldAtPosition(const Position&) const;

    void run();

private:
    
    Options opts_;

    std::unique_ptr<mfem::DG_FECollection> fec_;
    std::unique_ptr<mfem::FiniteElementSpace> fes_;

    mfem::Mesh mesh_;

    std::unique_ptr<mfem::LinearForm> inflowForm_;

    std::unique_ptr<mfem::BilinearForm> MInv_;
    std::unique_ptr<mfem::BilinearForm> Kx_;
    std::unique_ptr<mfem::BilinearForm> Ky_;

    mfem::GridFunction ez_, hx_, hy_;

    std::unique_ptr<mfem::ParaViewDataCollection> pd_;
    
    void checkOptionsAreValid(const Options&, const mfem::Mesh&);
    std::unique_ptr<mfem::LinearForm> buildInflowForm() const;
    std::unique_ptr<mfem::BilinearForm> buildMassMatrix() const;
    std::unique_ptr<mfem::BilinearForm> buildDerivativeOperator(const Direction&) const;

    void initializeParaviewData();
};

class DG_Solver : public mfem::Solver
{
private:
    mfem::SparseMatrix& M, & K, A;
    mfem::GMRESSolver linear_solver;
    mfem::BlockILU prec;
    double dt;
public:
    DG_Solver(mfem::SparseMatrix& M_, mfem::SparseMatrix& K_, const mfem::FiniteElementSpace& fes)
        : M(M_),
        K(K_),
        prec(fes.GetFE(0)->GetDof(),
            mfem::BlockILU::Reordering::MINIMUM_DISCARDED_FILL),
        dt(-1.0)
    {
        linear_solver.iterative_mode = false;
        linear_solver.SetRelTol(1e-9);
        linear_solver.SetAbsTol(0.0);
        linear_solver.SetMaxIter(100);
        linear_solver.SetPrintLevel(0);
        linear_solver.SetPreconditioner(prec);
    }

    void SetTimeStep(double dt_)
    {
        if (dt_ != dt)
        {
            dt = dt_;
            // Form operator A = M - dt*K
            A = K;
            A *= -dt;
            A += M;

            // this will also call SetOperator on the preconditioner
            linear_solver.SetOperator(A);
        }
    }

    void SetOperator(const mfem::Operator& op)
    {
        linear_solver.SetOperator(op);
    }

    virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const
    {
        linear_solver.Mult(x, y);
    }
};

class FE_Evolution : public mfem::TimeDependentOperator
{
private:
    mfem::BilinearForm& M, & K;
    const mfem::Vector& b;
    mfem::Solver* M_prec;
    mfem::CGSolver M_solver;
    DG_Solver* dg_solver;

    mutable mfem::Vector z;

public:
    FE_Evolution(mfem::BilinearForm& M_, mfem::BilinearForm& K_, const mfem::Vector& b_);

    virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const;
    virtual void ImplicitSolve(const double timestep, const mfem::Vector& x, mfem::Vector& k);

    virtual ~FE_Evolution();
};

}

