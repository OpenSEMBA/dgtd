#include "Solver1D.h"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace mfem;

namespace Maxwell {

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

	/** A time-dependent operator for the right-hand side of the ODE. The DG weak
		form of du/dt = -v.grad(u) is M du/dt = K u + b, where M and K are the mass
		and advection matrices, and b describes the flow on the boundary. This can
		be written as a general ODE, du/dt = M^{-1} (K u + b), and this class is
		used to evaluate the right-hand side. */
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

	Solver1D::Solver1D(const Options& opts, const Mesh& mesh)
	{
		checkOptionsAreValid(opts, mesh);

		mesh_ = mfem::Mesh(mesh, true);
		opts_ = opts;

		fec_ = std::make_unique<DG_FECollection>(opts_.order, mesh_.Dimension(), BasisType::GaussLobatto);
		fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());

		ODESolver* ode_solver = NULL;
		ode_solver = new RK4Solver;
		
		MInv_ = buildInverseMassMatrix();
		KxE_ = buildDerivativeAndFluxOperator(X, Electric);
		KxH_ = buildDerivativeAndFluxOperator(X, Magnetic);

		Ez_.SetSpace(fes_.get());
		Ez_.ProjectCoefficient(ConstantCoefficient(0.0));

		Hy_.SetSpace(fes_.get());
		Hy_.ProjectCoefficient(ConstantCoefficient(0.0));

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

	std::unique_ptr<mfem::BilinearForm> Solver1D::buildInverseMassMatrix() const
	{
		auto MInv = std::make_unique<BilinearForm>(fes_.get());
		MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator));
		MInv->Assemble();
		MInv->Finalize();
		return MInv;
	}

	std::unique_ptr<mfem::BilinearForm> Solver1D::buildDerivativeAndFluxOperator(const Direction& d, const FieldType& ft) const
	{
		assert(d == X, "Incorrect argument for direction.");

		auto K = std::make_unique<BilinearForm>(fes_.get());
		
		ConstantCoefficient one(1.0);

		K->AddDomainIntegrator(
			new TransposeIntegrator(
				new DerivativeIntegrator(one, d)));

		std::vector<VectorConstantCoefficient> n = {
			VectorConstantCoefficient(Vector({1.0})),
		};

		double alpha;
		double beta;

		if (ft == Electric)
		{
			alpha = -1.0;
			beta = 0.0;
			K->AddInteriorFaceIntegrator(
				new DGTraceIntegrator(n[d], alpha, beta));
			K->AddBdrFaceIntegrator(
				new DGTraceIntegrator(n[d], alpha, beta));
		}
		else
		{
			alpha = -1.0;
			beta = 0.0;
			K->AddInteriorFaceIntegrator(
				new DGTraceIntegrator(n[d], alpha, beta));
			K->AddBdrFaceIntegrator(
				new DGTraceIntegrator(n[d], alpha, beta));
		}

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
			KxH_->Mult(Hy_, aux);
			MInv_->Mult(aux, ezNew);
			ezNew *= -opts_.dt;
			ezNew.Add(1.0, Ez_);

			// Update H.
			KxE_->Mult(ezNew, aux);
			MInv_->Mult(aux, hyNew);
			hyNew *= -opts_.dt;
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

	//void Solver1D::runODESolver()
	//{
	//	FE_Evolution adv(*M_, *Kx_, B_);

	//	double t = 0.0;
	//	adv.SetTime(t);

	//	ODESolver* ode_solver = NULL;
	//	ode_solver = new RK4Solver;
	//	ode_solver->Init(adv);

	//	pd_->SetCycle(0);
	//	pd_->SetTime(0.0);
	//	pd_->Save();

	//	bool done = false;
	//	for (int cycle = 0; !done;)
	//	{
	//		double dt_real = std::min(opts_.dt, opts_.t_final - t);
	//		ode_solver->Step(Ez_, t, dt_real);
	//		//ode_solver->Step(Hy_, t, dt_real);
	//		cycle++;

	//		done = (t >= opts_.t_final - 1e-8 * opts_.dt);

	//		if (done || cycle % opts_.vis_steps == 0)
	//		{

	//			pd_->SetCycle(cycle);
	//			pd_->SetTime(t);
	//			pd_->Save();

	//		}

	//	}
	//}
}