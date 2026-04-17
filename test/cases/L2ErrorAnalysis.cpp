// #include "L2ErrorAnalysis.h"

// #include <iostream>
// #include <fstream>
// #include <regex>
// #include <cmath>
// #include <stdexcept>
// #include <algorithm>
// #include <memory> // For smart pointers

// namespace maxwell {

// const double SQRT2 = std::sqrt(2.0);
// const int m = 5;
// const int n = 5;
// const double w = M_PI * std::sqrt(double(m*m + n*n));

// // Helper for Bessel logic
// inline double bessel_j(int n, double x) {
//     #ifdef _WIN32 
//         return _jn(n, x);
//     #else
//         return jn(n, x);
//     #endif
// }

// // -------------------------------------------------------------------------
// // TM55 EXACT SOLUTIONS (Resonant Box)
// // -------------------------------------------------------------------------

// // Ez = sin(m pi x) * sin(n pi y) * cos(omega t)
// double TM55_Ez_Exact::Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) {
//     mfem::Vector p; 
//     T.Transform(ip, p);
//     return std::sin(double(m)*M_PI*p(0)) * std::sin(double(n)*M_PI*p(1)) * std::cos(w * GetTime());
// }

// // Hx = - (ky / omega) * sin(m pi x) * cos(n pi y) * sin(omega t)
// double TM55_Hx_Exact::Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) {
//     mfem::Vector p; 
//     T.Transform(ip, p);
//     double spatial = std::sin(double(m)*M_PI*p(0)) * std::cos(double(n)*M_PI*p(1));
//     double amplitude = -(double(n) * M_PI) / w;
//     return amplitude * spatial * std::sin(w * GetTime());
// }

// // Hy = + (kx / omega) * cos(m pi x) * sin(n pi y) * sin(omega t)
// double TM55_Hy_Exact::Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) {
//     mfem::Vector p; 
//     T.Transform(ip, p);
//     double spatial = std::cos(double(m)*M_PI*p(0)) * std::sin(double(n)*M_PI*p(1));
//     double amplitude = (double(m) * M_PI) / w;
//     return amplitude * spatial * std::sin(w * GetTime());
// }

// // -------------------------------------------------------------------------
// // BESSEL J6 EXACT SOLUTIONS (2D Ring)
// // -------------------------------------------------------------------------

// class BesselJ6_Ez_Exact : public mfem::Coefficient {
// public:
//     double alpha = 13.589290170541217;
//     double omega = 13.589290170541217;

//     double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override {
//         mfem::Vector p;
//         T.Transform(ip, p);
//         double x = p(0);
//         double y = p(1);
//         double t = GetTime();

//         double r = std::sqrt(x*x + y*y);
//         double theta = std::atan2(y, x);

//         double spatial = bessel_j(6, alpha * r) * std::cos(6.0 * theta);
//         return spatial * std::cos(omega * t);
//     }
// };

// class BesselJ6_Hx_Exact : public mfem::Coefficient {
// public:
//     double alpha = 13.589290170541217;
//     double omega = 13.589290170541217;

//     double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override {
//         mfem::Vector p;
//         T.Transform(ip, p);
//         double x = p(0);
//         double y = p(1);
//         double t = GetTime();

//         double r = std::sqrt(x*x + y*y);
//         double theta = std::atan2(y, x);

//         double sin_theta = y / r;
//         double cos_theta = x / r;
//         double cos_6theta = std::cos(6.0 * theta);
//         double sin_6theta = std::sin(6.0 * theta);

//         double arg = alpha * r;
//         double J6_prime = 0.5 * (bessel_j(5, arg) - bessel_j(7, arg));
//         double dEz_dr = alpha * J6_prime * cos_6theta;

//         double J6 = bessel_j(6, arg);
//         double dEz_dtheta = J6 * (-6.0 * sin_6theta);

//         double dEz_dy = sin_theta * dEz_dr + (cos_theta / r) * dEz_dtheta;
//         double amplitude = -1.0 / omega;
        
//         return amplitude * dEz_dy * std::sin(omega * t);
//     }
// };

// class BesselJ6_Hy_Exact : public mfem::Coefficient {
// public:
//     double alpha = 13.589290170541217;
//     double omega = 13.589290170541217;

//     double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override {
//         mfem::Vector p;
//         T.Transform(ip, p);
//         double x = p(0);
//         double y = p(1);
//         double t = GetTime();

//         double r = std::sqrt(x*x + y*y);
//         double theta = std::atan2(y, x);

//         double sin_theta = y / r;
//         double cos_theta = x / r;
//         double cos_6theta = std::cos(6.0 * theta);
//         double sin_6theta = std::sin(6.0 * theta);

//         double arg = alpha * r;
//         double J6_prime = 0.5 * (bessel_j(5, arg) - bessel_j(7, arg));
//         double dEz_dr = alpha * J6_prime * cos_6theta;

//         double J6 = bessel_j(6, arg);
//         double dEz_dtheta = J6 * (-6.0 * sin_6theta);

//         double dEz_dx = cos_theta * dEz_dr - (sin_theta / r) * dEz_dtheta;
//         double amplitude = 1.0 / omega;

//         return amplitude * dEz_dx * std::sin(omega * t);
//     }
// };

// // -------------------------------------------------------------------------
// // L2 ERROR ANALYSIS IMPLEMENTATION
// // -------------------------------------------------------------------------

// L2ErrorAnalysis::L2ErrorAnalysis(const std::string& case_path) 
//     : case_path_(case_path), 
//       mesh_(1, 1, 1)
// {
//     std::filesystem::path root(case_path_);
//     std::filesystem::path mesh_dir = root / "DomainSnapshotProbes" / "meshes";
//     std::filesystem::path rank0_dir = root / "DomainSnapshotProbes" / "rank_0";
//     std::filesystem::path stats_dir = root / "SimulationStats";

//     std::string mesh_file;
//     for (const auto& entry : std::filesystem::directory_iterator(mesh_dir)) {
//         std::string fname = entry.path().filename().string();
//         if (fname.find("mesh_rank0") != std::string::npos || 
//             fname.find("mesh_rank_0") != std::string::npos) {
//             mesh_file = entry.path().string();
//             break;
//         }
//     }
//     if (mesh_file.empty()) throw std::runtime_error("Could not find rank 0 mesh in " + mesh_dir.string());
    
//     mesh_ = Mesh::LoadFromFile(mesh_file.c_str(), 1, 0, false);

//     t2f_ = buildTimeToFields(rank0_dir.string(), mesh_);

//     if (!std::filesystem::exists(stats_dir)) {
//         std::filesystem::create_directories(stats_dir);
//     }

//     computeAndWriteErrors(stats_dir.string());
// }

// TimeToFields L2ErrorAnalysis::buildTimeToFields(const std::string& rank_root, Mesh& mesh)
// {
//     TimeToFields t2f;
//     std::vector<std::pair<double, std::filesystem::path>> time_dirs;
    
//     for (const auto& entry : std::filesystem::directory_iterator(rank_root)) {
//         if (!entry.is_directory()) continue;
//         if (!checkHasAllFieldFiles(entry.path())) continue;
//         time_dirs.emplace_back(readTimeFile(entry.path()), entry.path());
//     }
    
//     std::sort(time_dirs.begin(), time_dirs.end(),
//               [](const auto& a, const auto& b){ return a.first < b.first; });

//     for (const auto& [t, dir] : time_dirs) {
//         TimeFields tfs;
//         tfs.Ex = loadGridFunction(dir / "Ex.gf", mesh);
//         tfs.Ey = loadGridFunction(dir / "Ey.gf", mesh);
//         tfs.Ez = loadGridFunction(dir / "Ez.gf", mesh);
//         tfs.Hx = loadGridFunction(dir / "Hx.gf", mesh);
//         tfs.Hy = loadGridFunction(dir / "Hy.gf", mesh);
//         tfs.Hz = loadGridFunction(dir / "Hz.gf", mesh);
//         t2f.emplace(t, std::move(tfs));
//     }
//     return t2f;
// }

// void L2ErrorAnalysis::computeAndWriteErrors(const std::string& output_dir)
// {
//     // 1. Define File Paths
//     std::filesystem::path snap_file = std::filesystem::path(output_dir) / "L2ErrorSnapshots.txt";
//     std::filesystem::path total_file = std::filesystem::path(output_dir) / "L2ErrorTotal.txt";

//     std::ofstream out_snap(snap_file);
//     std::ofstream out_total(total_file);
    
//     out_snap << "Time\tL2_Ex\tL2_Ey\tL2_Ez\tL2_Hx\tL2_Hy\tL2_Hz\tL2_Total\n";

//     // --- AUTOMATIC SELECTION LOGIC ---
    
//     // Pointers to the abstract base class 'mfem::Coefficient'
//     mfem::Coefficient *ez_ptr = nullptr;
//     mfem::Coefficient *hx_ptr = nullptr;
//     mfem::Coefficient *hy_ptr = nullptr;

//     // Instance storage
//     TM55_Ez_Exact tm55_ez; 
//     TM55_Hx_Exact tm55_hx; 
//     TM55_Hy_Exact tm55_hy;

//     BesselJ6_Ez_Exact bes_ez;
//     BesselJ6_Hx_Exact bes_hx;
//     BesselJ6_Hy_Exact bes_hy;

//     // Detect mode based on folder name (case_path_)
//     // Usually these are "2D_Bessel..." or "2D_Resonant..."
//     bool is_bessel = (case_path_.find("Bessel") != std::string::npos);

//     if (is_bessel) {
//         std::cout << "Detected Bessel Case. Using BesselJ6 exact solutions." << std::endl;
//         ez_ptr = &bes_ez;
//         hx_ptr = &bes_hx;
//         hy_ptr = &bes_hy;
//     } else {
//         std::cout << "Detected Resonant Box Case (default). Using TM55 exact solutions." << std::endl;
//         ez_ptr = &tm55_ez;
//         hx_ptr = &tm55_hx;
//         hy_ptr = &tm55_hy;
//     }

//     mfem::ConstantCoefficient zero_exact(0.0);

//     double space_time_error_sq = 0.0;
//     double prev_time = t2f_.empty() ? 0.0 : t2f_.begin()->first;

//     std::cout << "Computing L2 Errors for " << t2f_.size() << " snapshots..." << std::endl;

//     for (const auto& [time, fields] : t2f_) {
        
//         // Polymorphism handles the SetTime call for whichever object is pointed to
//         ez_ptr->SetTime(time);
//         hx_ptr->SetTime(time);
//         hy_ptr->SetTime(time);

//         double err_Ex = fields.Ex->ComputeL2Error(zero_exact);
//         double err_Ey = fields.Ey->ComputeL2Error(zero_exact);
//         double err_Ez = fields.Ez->ComputeL2Error(*ez_ptr); // Use pointer

//         double err_Hx = fields.Hx->ComputeL2Error(*hx_ptr); // Use pointer
//         double err_Hy = fields.Hy->ComputeL2Error(*hy_ptr); // Use pointer
//         double err_Hz = fields.Hz->ComputeL2Error(zero_exact);

//         double sum_sq_spatial = err_Ex*err_Ex + err_Ey*err_Ey + err_Ez*err_Ez +
//                                 err_Hx*err_Hx + err_Hy*err_Hy + err_Hz*err_Hz;
        
//         double total_instant = std::sqrt(sum_sq_spatial);

//         // Right Rectangle Rule
//         double dt = time - prev_time;
//         if (dt > 0.0) {
//             space_time_error_sq += sum_sq_spatial * dt;
//         }
//         prev_time = time;

//         out_snap << time << "\t" 
//                  << err_Ex << "\t" << err_Ey << "\t" << err_Ez << "\t" 
//                  << err_Hx << "\t" << err_Hy << "\t" << err_Hz << "\t" 
//                  << total_instant << "\n";
//     }

//     double global_L2_error = std::sqrt(space_time_error_sq);

//     out_total << "Global_Space_Time_L2_Error\n";
//     out_total << global_L2_error << "\n";
// }

// }