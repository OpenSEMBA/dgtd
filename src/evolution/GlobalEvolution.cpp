#include "GlobalEvolution.h"

#include <chrono>

namespace maxwell {

GlobalEvolution::GlobalEvolution(
mfem::ParFiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, EvolutionOptions& options) :
mfem::TimeDependentOperator(numberOfFieldComponents* numberOfMaxDimensions* fes.GetNDofs()),
fes_{ fes },
model_{ model },
srcmngr_{ srcmngr },
opts_{ options }
{

	fes_.ExchangeFaceNbrData();

	for (auto d = X; d <= Z; d++){
		eOld_[d].SetSpace(&fes_);
		hOld_[d].SetSpace(&fes_);
	}

	globalOperator_ = std::make_unique<mfem::SparseMatrix>(numberOfFieldComponents * numberOfMaxDimensions * fes_.GetNDofs(), numberOfFieldComponents * numberOfMaxDimensions * (fes_.GetNDofs() + fes_.num_face_nbr_dofs));

#ifdef SHOW_TIMER_INFORMATION
	if (Mpi::WorldRank() == 0){
		std::cout << "------------------------------------------------" << std::endl;
		std::cout << "---------OPERATOR ASSEMBLY INFORMATION----------" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		std::cout << std::endl;
	}
#endif


	Probes probes;
	if (model_.getTotalFieldScatteredFieldToMarker().find(BdrCond::TotalFieldIn) != model_.getTotalFieldScatteredFieldToMarker().end()) {

		srcmngr_.initTFSFPreReqs(model_.getConstMesh(), model_.getTotalFieldScatteredFieldToMarker().at(BdrCond::TotalFieldIn));

		auto globalTFSFfes = srcmngr_.getGlobalTFSFSpace();
		auto tfsfMesh = globalTFSFfes->GetMesh();

		Model tfsfModel = Model(*tfsfMesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(GeomTagToBoundary{}, GeomTagToInteriorBoundary{}));
		
		ProblemDescription tfsfpd(tfsfModel, probes, srcmngr_.sources, opts_);
		DGOperatorFactory<FiniteElementSpace> tfsfops(tfsfpd, *globalTFSFfes);

		TFSFOperator_ = tfsfops.buildTFSFGlobalOperator();

		auto src_sm = static_cast<mfem::SubMesh*>(srcmngr_.getGlobalTFSFSpace()->GetMesh());
		mfem::SubMeshUtils::BuildVdofToVdofMap(*srcmngr_.getGlobalTFSFSpace(), fes_, src_sm->GetFrom(), src_sm->GetParentElementIDMap(), sub_to_parent_ids_);
	}

	ProblemDescription pd(model_, probes, srcmngr_.sources, opts_);
	DGOperatorFactory<mfem::ParFiniteElementSpace> dgops(pd, fes_);

	globalOperator_ = dgops.buildGlobalOperator();

}

const mfem::Vector buildSingleVectorTFSFFunc(const FieldGridFuncs& func)
{
	mfem::Vector res(6 * func[0][0].Size());
	res.UseDevice(true);
	for (auto f : { E, H }) {
		for (auto d : { X, Y, Z }) {
			for (auto v{ 0 }; v < func[f][d].Size(); v++) {
				res[v + (f * 3 + d) * func[f][d].Size()] = func[f][d][v];
			}
		}
	}
	return res;
}

void assertVectorOnDevice(const Vector &v, const std::string &name)
{
    MemoryType mem_type = v.GetMemory().GetMemoryType();
    if (mem_type != MemoryType::DEVICE)
    {
        mfem::out << "Warning: Vector '" << name << "' latest data is NOT on device!\n";
    }
    else
    {
        mfem::out << "Vector '" << name << "' latest data is on DEVICE.\n";
    }
}

void GlobalEvolution::ImplicitSolve(const double dt,
                                             const mfem::Vector& x,
                                             mfem::Vector& k)
{
    // We solve G(k) = k - f(x + dt*k, t+dt) = 0 with Newton–Krylov.
    // Only calls Mult() at t+dt. No need to separate TFSF.

    const double t0 = GetTime();
    SetTime(t0 + dt); // Backward Euler evaluates at t + dt

    const int n = x.Size();

    // Work vectors
    mfem::Vector y(n), fy(n), G(n), delta(n);
    y.UseDevice(true); fy.UseDevice(true); G.UseDevice(true); delta.UseDevice(true);

    // Initial guess: use explicit evaluation at t+dt
    // (MFEM may pass a reused k; you can keep it if you prefer.)
    Mult(x, k); // k0 = f(x, t+dt)

    // Helper to compute G(k) and cache y = x + dt*k and fy = f(y, t+dt)
    auto compute_G = [&](const mfem::Vector& kin) {
        y = x;
        y.Add(dt, kin);   // y = x + dt*k
        Mult(y, fy);      // fy = f(y, t+dt)
        G  = kin;
        G -= fy;          // G = k - f(y)
    };

    // Newton params (tweak if needed)
    const int    newton_max_it = 8;
    const double newton_rtol   = 1e-9;
    const double newton_atol   = 1e-12;

    for (int it = 0; it < newton_max_it; ++it) {
        compute_G(k);

        const double gnorm  = G.Norml2();
        const double target = newton_atol + newton_rtol * std::max(1.0, k.Norml2());
        if (gnorm <= target) { break; }

        // Linear operator: J(v) ≈ v - dt * ( f(y + η v) - f(y) ) / η
        // Finite-diff directional derivative; only calls Mult().
        class JOp : public mfem::Operator {
            GlobalEvolution& ge_;
            const mfem::Vector& y_;
            const mfem::Vector& fy_;
            const double dt_;
            mutable mfem::Vector ypert_, fpert_, diff_;
            const double eps_; // finite-diff scaling
        public:
            JOp(GlobalEvolution& ge, const mfem::Vector& y, const mfem::Vector& fy, double dt)
            : mfem::Operator(y.Size()), ge_(ge), y_(y), fy_(fy), dt_(dt),
              ypert_(y.Size()), fpert_(fy.Size()), diff_(fy.Size()),
              eps_(1e-7) {
                ypert_.UseDevice(true); fpert_.UseDevice(true); diff_.UseDevice(true);
            }
            void Mult(const mfem::Vector& v, mfem::Vector& Jv) const override {
                const double vnorm = v.Norml2();
                if (vnorm == 0.0) {
                    Jv.SetSize(v.Size()); Jv = 0.0; return;
                }
                const double ynorm = y_.Norml2();
                const double eta   = eps_ * (1.0 + ynorm) / vnorm;

                ypert_ = y_;
                ypert_.Add(eta, v);          // y + η v

                ge_.Mult(ypert_, fpert_);    // f(y + η v, t+dt)

                diff_  = fpert_;
                diff_ -= fy_;                // f(y+ηv) - f(y)

                Jv.SetSize(v.Size());
                Jv = v;                      // v
                Jv.Add(-dt_/eta, diff_);     // - dt/η * (f(y+ηv)-f(y))
            }
        } J(*this, y, fy, dt);

        mfem::GMRESSolver gmres;
        gmres.SetOperator(J);
        gmres.SetRelTol(1e-8);
        gmres.SetMaxIter(400);
        gmres.SetPrintLevel(0);
        // If you have a preconditioner, set it here:
        // gmres.SetPreconditioner(P);

        // Solve J * delta = -G
        mfem::Vector rhs = G;
        rhs *= -1.0;
        gmres.Mult(rhs, delta);

        // Update k
        k += delta;
        // (Optional) line search could go here if you see stagnation.
    }

    // Restore original time
    SetTime(t0);
}


void GlobalEvolution::Mult(const mfem::Vector& in, mfem::Vector& out) const
{
    mfem::StopWatch timerTotal, timerExchange, timerAssembleInNew, timerMult, timerTFSF;

	timerTotal.Start();
	
	const auto ndofs = fes_.GetNDofs();
	const auto nbrDofs = fes_.num_face_nbr_dofs;
	const auto blockSize = ndofs + nbrDofs;

    timerExchange.Start();
	#ifdef SEMBA_DGTD_ENABLE_CUDA
	load_in_to_eh_gpu(in, eOld_, hOld_, ndofs);
	#else
	for (auto d = X; d <= Z; d++){
		for (auto v = 0; v < ndofs; v++){
			eOld_[d][v] = in[d * ndofs + v];
			hOld_[d][v] = in[(d+3) * ndofs + v];
		}
	}
	#endif
	for (auto d = X; d <= Z; d++){
		eOld_[d].ExchangeFaceNbrData();
		hOld_[d].ExchangeFaceNbrData();
	}
    timerExchange.Stop();
    
	timerAssembleInNew.Start();
	Vector inNew(6*blockSize);
	#ifdef SEMBA_DGTD_ENABLE_CUDA
	inNew.UseDevice(true);
	load_eh_to_innew_gpu(in, inNew, ndofs, nbrDofs);
	load_nbr_to_innew_gpu(eOld_, hOld_, inNew, ndofs, nbrDofs);
	#else
	for (int d = X; d <= Z; ++d)
	{
		// Copy locals
		MFEM_ASSERT(eOld_[d].Size() == ndofs, "eOld size mismatch");
		MFEM_ASSERT(hOld_[d].Size() == ndofs, "hOld size mismatch");

		inNew.SetVector(eOld_[d],         d      * blockSize);
		inNew.SetVector(hOld_[d],    (3 + d)     * blockSize);

		// Copy neighbor dofs (if any)
		if (nbrDofs)
		{
			mfem::Vector &eNbr = eOld_[d].FaceNbrData();
			mfem::Vector &hNbr = hOld_[d].FaceNbrData();

			MFEM_ASSERT(eNbr.Size() == nbrDofs, "e-nbr size mismatch");
			MFEM_ASSERT(hNbr.Size() == nbrDofs, "h-nbr size mismatch");

			inNew.SetVector(eNbr,      d      * blockSize + ndofs);
			inNew.SetVector(hNbr, (3 + d)     * blockSize + ndofs);
		}
	}
	#endif
	timerAssembleInNew.Stop();

	timerMult.Start();
    globalOperator_->Mult(inNew, out);
	timerMult.Stop();

	timerTFSF.Start();
	bool early_tfsf_deletion = false;
    for (auto& source : srcmngr_.sources) {
        if (dynamic_cast<TotalField*>(source.get()) && srcmngr_.getGlobalTFSFSpace() != nullptr) {
			auto sourcecast = dynamic_cast<TotalField*>(source.get());
			if(opts_.tfsfFinalTime != 0.0 
			&& std::abs(GetTime() - opts_.tfsfFinalTime) <= 1e-5){
				#ifdef SEMBA_DGTD_ENABLE_CUDA
				const auto func = eval_time_var_field_gpu(GetTime(), srcmngr_);
            	mfem::Vector assembledFunc = load_tfsf_into_single_vector_gpu(func);
				#else
				auto func = evalTimeVarFunction(GetTime(), srcmngr_);
            	mfem::Vector assembledFunc = buildSingleVectorTFSFFunc(func);
				#endif
				if(assembledFunc.Norml2() >= 1e-1){
					early_tfsf_deletion = true;
					MFEM_WARNING("TFSF SOURCE DELETED WHEN FIELDS ARE STILL BEING LOADED!");
				}
				source.reset(nullptr);
				break;
			}
			#ifdef SEMBA_DGTD_ENABLE_CUDA
            const auto func = eval_time_var_field_gpu(GetTime(), srcmngr_);
            mfem::Vector assembledFunc = load_tfsf_into_single_vector_gpu(func);
			#else
			auto func = evalTimeVarFunction(GetTime(), srcmngr_);
            mfem::Vector assembledFunc = buildSingleVectorTFSFFunc(func);
			#endif
            mfem::Vector tempTFSF(assembledFunc.Size());
			tempTFSF.UseDevice(true);
            TFSFOperator_->Mult(assembledFunc, tempTFSF);
			#ifdef SEMBA_DGTD_ENABLE_CUDA
			load_tfsf_into_out_vector_gpu(sub_to_parent_ids_, tempTFSF, out, fes_.GetNDofs(), srcmngr_.getGlobalTFSFSpace()->GetNDofs());
			#else
			for (auto f : { E, H }) {
                for (auto d : { X, Y, Z }) {
                    for (int v = 0; v < sub_to_parent_ids_.Size(); v++) {
                        const int outIdx = (f * 3 + d) * ndofs + sub_to_parent_ids_[v];
                        const int tempIdx = (f * 3 + d) * srcmngr_.getGlobalTFSFSpace()->GetNDofs() + v;
                        if (tempTFSF[tempIdx] != 0.0) {
                            out[outIdx] -= tempTFSF[tempIdx];
                        }
                    }
                }
            }
			#endif
		}
	}
	if (early_tfsf_deletion){
		MFEM_WARNING("TFSF SOURCE DELETED WHEN FIELDS ARE STILL BEING LOADED!");
	}
	timerTFSF.Stop();

    timerTotal.Stop();

			
	static double lastPrintTime = -1.0;
	const double printInterval = 0.1;  
	double currentTime = GetTime();

	if (lastPrintTime < 0.0 || currentTime >= lastPrintTime + printInterval) {
		std::cout << "Current time: " << currentTime << std::endl;
		std::cout << "Rank " << Mpi::WorldRank() 
				<< " Mult total: " << timerTotal.RealTime() * 1000 
				<< " ms, exchange: " << timerExchange.RealTime() * 1000 
				<< " ms, assembleIn: " << timerAssembleInNew.RealTime() * 1000 << " ms\n";
		std::cout << "Rank " << Mpi::WorldRank() 
				<< " Mult mult: " << timerMult.RealTime() * 1000 
				<< " ms, tfsf: " << timerTFSF.RealTime() * 1000 << " ms\n";

		lastPrintTime += printInterval;
		if (lastPrintTime < 0.0) lastPrintTime = currentTime;
	}


}


}
