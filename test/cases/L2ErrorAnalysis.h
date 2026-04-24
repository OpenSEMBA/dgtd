#pragma once

#include "FieldUtils.h"

#include "mfem.hpp"
#include <string>
#include <vector>
#include <map>
#include <filesystem>
#include <memory>

namespace maxwell {

using namespace mfem;

class L2ErrorAnalysis {
public:
    L2ErrorAnalysis(const std::string& case_path);

private:
    void loadData(const std::string& case_path);
    Mesh loadMesh(const std::string& mesh_path);
    TimeToFields buildTimeToFields(const std::string& rank_root, Mesh& mesh);
    
    void computeAndWriteErrors(const std::string& output_dir);

    std::string case_path_;
    Mesh mesh_;
    TimeToFields t2f_;
};

class TM55_Ez_Exact : public mfem::Coefficient {
public:
    double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override;
};

class TM55_Hx_Exact : public mfem::Coefficient {
public:
    double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override;
};

class TM55_Hy_Exact : public mfem::Coefficient {
public:
    double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override;
};

}