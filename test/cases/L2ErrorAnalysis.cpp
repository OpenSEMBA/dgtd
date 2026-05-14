#include "L2ErrorAnalysis.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <stdexcept>

namespace maxwell {

const int m = 5;
const int n = 5;
const double w = M_PI * std::sqrt(double(m * m + n * n));

inline double bessel_j(int order, double x)
{
#ifdef _WIN32
    return _jn(order, x);
#else
    return jn(order, x);
#endif
}

// -------------------------------------------------------------------------
// TM55 EXACT SOLUTIONS (Resonant Box)
// -------------------------------------------------------------------------

double TM55_Ez_Exact::Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip)
{
    mfem::Vector p;
    T.Transform(ip, p);
    return std::sin(double(m) * M_PI * p(0)) * std::sin(double(n) * M_PI * p(1)) * std::cos(w * GetTime());
}

double TM55_Hx_Exact::Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip)
{
    mfem::Vector p;
    T.Transform(ip, p);
    const double spatial = std::sin(double(m) * M_PI * p(0)) * std::cos(double(n) * M_PI * p(1));
    const double amplitude = -(double(n) * M_PI) / w;
    return amplitude * spatial * std::sin(w * GetTime());
}

double TM55_Hy_Exact::Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip)
{
    mfem::Vector p;
    T.Transform(ip, p);
    const double spatial = std::cos(double(m) * M_PI * p(0)) * std::sin(double(n) * M_PI * p(1));
    const double amplitude = (double(m) * M_PI) / w;
    return amplitude * spatial * std::sin(w * GetTime());
}

// -------------------------------------------------------------------------
// BESSEL J6 EXACT SOLUTIONS (2D Ring)
// -------------------------------------------------------------------------

class BesselJ6_Ez_Exact : public mfem::Coefficient {
public:
    double alpha = 13.589290170541217;
    double omega = 13.589290170541217;

    double Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override
    {
        mfem::Vector p;
        T.Transform(ip, p);
        const double x = p(0);
        const double y = p(1);
        const double t = GetTime();

        const double r = std::sqrt(x * x + y * y);
        const double theta = std::atan2(y, x);

        const double spatial = bessel_j(6, alpha * r) * std::cos(6.0 * theta);
        return spatial * std::cos(omega * t);
    }
};

class BesselJ6_Hx_Exact : public mfem::Coefficient {
public:
    double alpha = 13.589290170541217;
    double omega = 13.589290170541217;

    double Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override
    {
        mfem::Vector p;
        T.Transform(ip, p);
        const double x = p(0);
        const double y = p(1);
        const double t = GetTime();

        const double r = std::sqrt(x * x + y * y);
        const double theta = std::atan2(y, x);

        const double sin_theta = y / r;
        const double cos_theta = x / r;
        const double cos_6theta = std::cos(6.0 * theta);
        const double sin_6theta = std::sin(6.0 * theta);

        const double arg = alpha * r;
        const double J6_prime = 0.5 * (bessel_j(5, arg) - bessel_j(7, arg));
        const double dEz_dr = alpha * J6_prime * cos_6theta;

        const double J6 = bessel_j(6, arg);
        const double dEz_dtheta = J6 * (-6.0 * sin_6theta);

        const double dEz_dy = sin_theta * dEz_dr + (cos_theta / r) * dEz_dtheta;
        const double amplitude = -1.0 / omega;

        return amplitude * dEz_dy * std::sin(omega * t);
    }
};

class BesselJ6_Hy_Exact : public mfem::Coefficient {
public:
    double alpha = 13.589290170541217;
    double omega = 13.589290170541217;

    double Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override
    {
        mfem::Vector p;
        T.Transform(ip, p);
        const double x = p(0);
        const double y = p(1);
        const double t = GetTime();

        const double r = std::sqrt(x * x + y * y);
        const double theta = std::atan2(y, x);

        const double sin_theta = y / r;
        const double cos_theta = x / r;
        const double cos_6theta = std::cos(6.0 * theta);
        const double sin_6theta = std::sin(6.0 * theta);

        const double arg = alpha * r;
        const double J6_prime = 0.5 * (bessel_j(5, arg) - bessel_j(7, arg));
        const double dEz_dr = alpha * J6_prime * cos_6theta;

        const double J6 = bessel_j(6, arg);
        const double dEz_dtheta = J6 * (-6.0 * sin_6theta);

        const double dEz_dx = cos_theta * dEz_dr - (sin_theta / r) * dEz_dtheta;
        const double amplitude = 1.0 / omega;

        return amplitude * dEz_dx * std::sin(omega * t);
    }
};

L2ErrorAnalysis::L2ErrorAnalysis(const std::string& case_path)
    : case_path_(case_path),
      mesh_(1, 1, 1)
{
    std::filesystem::path root(case_path_);
    std::filesystem::path mesh_dir = root / "DomainSnapshotProbes" / "meshes";
    std::filesystem::path rank0_dir = root / "DomainSnapshotProbes" / "rank_0";
    std::filesystem::path stats_dir = root / "SimulationStats";

    std::string mesh_file;
    for (const auto& entry : std::filesystem::directory_iterator(mesh_dir)) {
        const std::string fname = entry.path().filename().string();
        if (fname.find("mesh_rank0") != std::string::npos ||
            fname.find("mesh_rank_0") != std::string::npos) {
            mesh_file = entry.path().string();
            break;
        }
    }
    if (mesh_file.empty()) {
        throw std::runtime_error("Could not find rank 0 mesh in " + mesh_dir.string());
    }

    std::ifstream mesh_stream(mesh_file);
    mesh_ = Mesh(mesh_stream, 0, 0, 0);

    t2f_ = buildTimeToFields(rank0_dir.string(), mesh_);

    if (!std::filesystem::exists(stats_dir)) {
        std::filesystem::create_directories(stats_dir);
    }

    computeAndWriteErrors(stats_dir.string());
}

TimeToFields L2ErrorAnalysis::buildTimeToFields(const std::string& rank_root, Mesh& mesh)
{
    return maxwell::buildTimeToFields(rank_root, mesh);
}

void L2ErrorAnalysis::computeAndWriteErrors(const std::string& output_dir)
{
    const std::filesystem::path snap_file = std::filesystem::path(output_dir) / "L2ErrorSnapshots.txt";
    const std::filesystem::path total_file = std::filesystem::path(output_dir) / "L2ErrorTotal.txt";

    std::ofstream out_snap(snap_file);
    std::ofstream out_total(total_file);

    out_snap << "Time\tL2_Ex\tL2_Ey\tL2_Ez\tL2_Hx\tL2_Hy\tL2_Hz\tL2_Total\n";

    mfem::Coefficient* ez_ptr = nullptr;
    mfem::Coefficient* hx_ptr = nullptr;
    mfem::Coefficient* hy_ptr = nullptr;

    TM55_Ez_Exact tm55_ez;
    TM55_Hx_Exact tm55_hx;
    TM55_Hy_Exact tm55_hy;

    BesselJ6_Ez_Exact bes_ez;
    BesselJ6_Hx_Exact bes_hx;
    BesselJ6_Hy_Exact bes_hy;

    const bool is_bessel = (case_path_.find("Bessel") != std::string::npos);

    if (is_bessel) {
        std::cout << "Detected Bessel Case. Using BesselJ6 exact solutions." << std::endl;
        ez_ptr = &bes_ez;
        hx_ptr = &bes_hx;
        hy_ptr = &bes_hy;
    } else {
        std::cout << "Detected Resonant Box Case (default). Using TM55 exact solutions." << std::endl;
        ez_ptr = &tm55_ez;
        hx_ptr = &tm55_hx;
        hy_ptr = &tm55_hy;
    }

    mfem::ConstantCoefficient zero_exact(0.0);

    double space_time_error_sq = 0.0;
    double prev_time = t2f_.empty() ? 0.0 : t2f_.begin()->first;

    std::cout << "Computing L2 Errors for " << t2f_.size() << " snapshots..." << std::endl;

    for (const auto& [time, fields] : t2f_) {
        ez_ptr->SetTime(time);
        hx_ptr->SetTime(time);
        hy_ptr->SetTime(time);

        const double err_Ex = fields.Ex->ComputeL2Error(zero_exact);
        const double err_Ey = fields.Ey->ComputeL2Error(zero_exact);
        const double err_Ez = fields.Ez->ComputeL2Error(*ez_ptr);

        const double err_Hx = fields.Hx->ComputeL2Error(*hx_ptr);
        const double err_Hy = fields.Hy->ComputeL2Error(*hy_ptr);
        const double err_Hz = fields.Hz->ComputeL2Error(zero_exact);

        const double sum_sq_spatial = err_Ex * err_Ex + err_Ey * err_Ey + err_Ez * err_Ez +
                                      err_Hx * err_Hx + err_Hy * err_Hy + err_Hz * err_Hz;

        const double total_instant = std::sqrt(sum_sq_spatial);

        const double dt = time - prev_time;
        if (dt > 0.0) {
            space_time_error_sq += sum_sq_spatial * dt;
        }
        prev_time = time;

        out_snap << time << "\t"
                 << err_Ex << "\t" << err_Ey << "\t" << err_Ez << "\t"
                 << err_Hx << "\t" << err_Hy << "\t" << err_Hz << "\t"
                 << total_instant << "\n";
    }

    const double global_L2_error = std::sqrt(space_time_error_sq);

    out_total << "Global_Space_Time_L2_Error\n";
    out_total << global_L2_error << "\n";
}

}