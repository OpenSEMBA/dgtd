#include "Solver1D.h"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace mfem;




namespace Maxwell {

	double inflow_function(const Vector& x)
	{
		return 0.0;
	}

	void velocity_function(const Vector& x, Vector& v)
	{
		int dim = x.Size();

		// map to the reference [-1,1] domain
		Vector pos(dim);
		double normalizedPos;
		double center = (0.0 + 1.0) * 0.5;
		normalizedPos = 2 * (pos[0] - center) / (1.0 - 0.0);

		v(0) = 1.0;
	}

	
	class DG_Solver : public Solver
	{
	private:
		SparseMatrix& M, & K, A;
		GMRESSolver linear_solver;
		BlockILU prec;
		double dt;
	public:
		DG_Solver(SparseMatrix& M_, SparseMatrix& K_, const FiniteElementSpace& fes)
			: M(M_),
			K(K_),
			prec(fes.GetFE(0)->GetDof(),
				BlockILU::Reordering::MINIMUM_DISCARDED_FILL),
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

		void SetOperator(const Operator& op)
		{
			linear_solver.SetOperator(op);
		}

		virtual void Mult(const Vector& x, Vector& y) const
		{
			linear_solver.Mult(x, y);
		}
	};

	class FE_Evolution : public TimeDependentOperator
	{
	private:
		BilinearForm& M, & K;
		const Vector& b;
		Solver* M_prec;
		CGSolver M_solver;
		DG_Solver* dg_solver;

		mutable Vector z;

	public:
		FE_Evolution(BilinearForm& M_, BilinearForm& K_, const Vector& b_);

		virtual void Mult(const Vector& x, Vector& y) const;
		virtual void ImplicitSolve(const double dt, const Vector& x, Vector& k);

		virtual ~FE_Evolution();
	};

	// Implementation of class FE_Evolution
	FE_Evolution::FE_Evolution(BilinearForm& M_, BilinearForm& K_, const Vector& b_)
		: TimeDependentOperator(M_.Height()), M(M_), K(K_), b(b_), z(M_.Height())
	{
		Array<int> ess_tdof_list;
		if (M.GetAssemblyLevel() == AssemblyLevel::LEGACY)
		{
			M_prec = new DSmoother(M.SpMat());
			M_solver.SetOperator(M.SpMat());
			dg_solver = new DG_Solver(M.SpMat(), K.SpMat(), *M.FESpace());
		}
		else
		{
			M_prec = new OperatorJacobiSmoother(M, ess_tdof_list);
			M_solver.SetOperator(M);
			dg_solver = NULL;
		}
		M_solver.SetPreconditioner(*M_prec);
		M_solver.iterative_mode = false;
		M_solver.SetRelTol(1e-9);
		M_solver.SetAbsTol(0.0);
		M_solver.SetMaxIter(100);
		M_solver.SetPrintLevel(0);
	}

	void FE_Evolution::Mult(const Vector& x, Vector& y) const
	{
		// y = M^{-1} (K x + b)
		K.Mult(x, z);
		z += b;
		M_solver.Mult(z, y);
	}

	void FE_Evolution::ImplicitSolve(const double dt, const Vector& x, Vector& k)
	{
		MFEM_VERIFY(dg_solver != NULL,
			"Implicit time integration is not supported with partial assembly");
		K.Mult(x, z);
		z += b;
		dg_solver->SetTimeStep(dt);
		dg_solver->Mult(z, k);
	}

	FE_Evolution::~FE_Evolution()
	{
		delete M_prec;
		delete dg_solver;
	}

	
	Solver1D::Solver1D(const Options& opts, const Mesh& mesh)
	{
		checkOptionsAreValid(opts, mesh);

		mesh_ = mfem::Mesh(mesh, true);
		opts_ = opts;
		int dim = mesh.Dimension();

		fec_ = std::make_unique<DG_FECollection>(opts_.order, mesh_.Dimension(), BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());

		boundaryTDoF_ = buildEssentialTrueDOF();

		M_ = buildMassMatrix();
		MInv_ = buildInverseMassMatrix();
		Kx_ = buildDerivativeAndFluxOperator(X);

		Ez_.SetSpace(fes_.get());
		Ez_.ProjectCoefficient(ConstantCoefficient(0.0));
		Ez_.ProjectBdrCoefficient(ConstantCoefficient(0.0), boundaryTDoF_);

		Hy_.SetSpace(fes_.get());
		Hy_.ProjectCoefficient(ConstantCoefficient(0.0));

		VectorFunctionCoefficient velocity(dim, velocity_function);
		FunctionCoefficient inflow(inflow_function);
		B_.Update(fes_.get());
		B_.AddBdrFaceIntegrator(
			new BoundaryFlowIntegrator(inflow, velocity, 1.0));
		B_.Assemble();

		initializeParaviewData();
	}

	void Solver1D::checkOptionsAreValid(const Options& opts, const Mesh& mesh)
	{
		if (mesh.Dimension() != 1) {
			throw std::exception("Incorrect Dimension for mesh");
		}
		if ((opts.order < 0) ||
			(opts.t_final < 0) ||
			(opts.dt < 0) ||
			(opts.vis_steps < 1) ||
			(opts.precision < 1)) {
			throw std::exception("Incorrect parameters in Options");
		}
	}

	mfem::Array<int> Solver1D::buildEssentialTrueDOF()
	{
		Array<int> ess_tdof_list;
		if (mesh_.bdr_attributes.Size())
		{
			Array<int> ess_bdr(mesh_.bdr_attributes.Max());
			ess_bdr = 1;
			fes_.get()->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
		}
		return ess_tdof_list;
	}

	std::unique_ptr<mfem::BilinearForm> Solver1D::buildMassMatrix() const
	{
		auto M = std::make_unique<BilinearForm>(fes_.get());
		M->AddDomainIntegrator((new MassIntegrator));
		M->Assemble();
		M->Finalize();
		return M;
	}

	std::unique_ptr<mfem::BilinearForm> Solver1D::buildInverseMassMatrix() const
	{
		auto MInv = std::make_unique<BilinearForm>(fes_.get());
		MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
		MInv->Assemble();
		MInv->Finalize();
		return MInv;
	}

	std::unique_ptr<mfem::BilinearForm> Solver1D::buildDerivativeAndFluxOperator(const Direction& d) const
	{

		assert(d == X, "Incorrect argument for direction.");

		ConstantCoefficient one(1.0);

		auto K = std::make_unique<BilinearForm>(fes_.get());
		K->AddDomainIntegrator(
			new TransposeIntegrator(
				new DerivativeIntegrator(one, X)));

		std::vector<VectorConstantCoefficient> n = {
			VectorConstantCoefficient(Vector({1.0})),
		};

		K->AddInteriorFaceIntegrator(
			new DGTraceIntegrator(n[0], -1.0, 0.0));
		K->AddBdrFaceIntegrator(
			new DGTraceIntegrator(n[0], -1.0, 0.0));

		K->Assemble();
		K->Finalize();

		return K;
	}

	void Solver1D::setInitialElectricField(std::function<ElectricField(const Position&)> f)
	{
		Ez_.ProjectCoefficient(FunctionCoefficient(f));
	}

	void Solver1D::initializeParaviewData()
	{
		pd_ = NULL;
		pd_ = std::make_unique<ParaViewDataCollection>("MaxwellView1D", &mesh_);
		pd_->SetPrefixPath("ParaView");
		pd_->RegisterField("Ez", &Ez_);
		pd_->RegisterField("Hy", &Hy_);
		pd_->SetLevelsOfDetail(opts_.order);
		pd_->SetDataFormat(VTKFormat::BINARY);
		opts_.order > 0 ? pd_->SetHighOrderOutput(true) : pd_->SetHighOrderOutput(false);
	}

	void Solver1D::run()
	{

		double time = 0.0;

		Vector aux(fes_->GetVSize());
		Vector ezNew(fes_->GetVSize());
		Vector hyNew(fes_->GetVSize());

		pd_->SetCycle(0);
		pd_->SetTime(0.0);
		pd_->Save();

		bool done = false;
		for (int cycle = 0; !done;)
		{

			// Update E.
			Kx_->Mult(Hy_, aux);
			MInv_->Mult(aux, ezNew);
			ezNew *= opts_.dt;
			ezNew.Add(1.0, Ez_);

			// Update H.
			Kx_->Mult(Ez_, aux);
			MInv_->Mult(aux, hyNew);
			hyNew *= opts_.dt;
			hyNew.Add(1.0, Hy_);

			Ez_ = ezNew;
			Hy_ = hyNew;


			time += opts_.dt;
			cycle++;

			done = (time >= opts_.t_final - 1e-8 * opts_.dt);

			if (done || cycle % opts_.vis_steps == 0) {
				pd_->SetCycle(cycle);
				pd_->SetTime(time);
				pd_->Save();
			}
		}
	}

	void Solver1D::runODESolver()
	{
		FE_Evolution adv(*M_, *Kx_, B_);

		double t = 0.0;
		adv.SetTime(t);

		ODESolver* ode_solver = NULL;
		ode_solver = new RK4Solver;
		ode_solver->Init(adv);

		pd_->SetCycle(0);
		pd_->SetTime(0.0);
		pd_->Save();

		bool done = false;
		for (int cycle = 0; !done;)
		{
			double dt_real = std::min(opts_.dt, opts_.t_final - t);
			ode_solver->Step(Ez_, t, dt_real);
			//ode_solver->Step(Hy_, t, dt_real);
			cycle++;

			done = (t >= opts_.t_final - 1e-8 * opts_.dt);

			if (done || cycle % opts_.vis_steps == 0)
			{

				pd_->SetCycle(cycle);
				pd_->SetTime(t);
				pd_->Save();

			}

		}
	}
}