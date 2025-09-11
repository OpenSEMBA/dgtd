#pragma once

#include <mfem.hpp>
#include <math.h>
#include <filesystem>
#include <string>
#include "components/Types.h"
#include "math/Function.h"

namespace maxwell {

using namespace mfem;
using Position = mfem::Vector;
using Rank = int;

enum FunctionType 
{
    Gaussian = 0,
    Resonant = 1,
    BesselJ62D = 2,
    BesselJ63D = 3,
    Planewave = 4,
    Dipole = 5
};

class ExcitationCoeffs
{
public:

    ExcitationCoeffs(const std::string& json_path);
    std::array<std::array<double, 3>, 2> FieldCompFactor;

private:

    void initFieldCompFactor();
    void loadInitialPolarizationValues(const FieldType& ft, const Vector& pol);
    void loadExcitedDirectionValues(const FieldType& ft, const Vector& pol);
};

Mesh loadMeshFromFile(const std::string& mesh_path);
GridFunction loadGridFunctionFromFile(const std::string& file_path, Mesh& mesh);

class L2SimDataCalculator
{
public:
    L2SimDataCalculator(const std::string& data_path, const std::string& json_path);

private:

    void loadMeshes(const std::string& data_path);
    void loadFES(const std::string& data_path);
    void loadNodepos(const std::string& data_path);
    void initFunction(const std::string& json_path);

    std::map<Rank, Mesh> meshes_;
    std::map<Rank, std::vector<Position>> nodepos_;
    std::unique_ptr<TimeFunction> function_;
};


}

