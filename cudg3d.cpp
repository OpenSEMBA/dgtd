// OpenSEMBA
// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
//                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
//                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
//                    Daniel Mateos Romero            (damarro@semba.guru)
//
// This file is part of OpenSEMBA.
//
// OpenSEMBA is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//

#include <iostream>
#include <stdlib.h>
#include <sys/time.h>

using namespace std;
using namespace SEMBA;

#include "Solver.h"
#include "solver/Info.h"

int main(int argc, const char *argv[]) {

    Argument::Parser arg;
    arg.args(argc, argv);
    Argument::OptionBase& input = arg.addOption(
            new Argument::Option<std::string>("Input", 'i', "input"));
    arg.addOption(
            (new Argument::Switch("Help", 'h', "help"))->defaultVal(false));
    Argument::Object opts = arg.parseKnownArgs().first;
    if (opts("Input").isNull()) {
        if (opts("Help").getBool()) {
            arg.printHelp();
            exit(EXIT_SUCCESS);
        } else {
            Argument::Error::Required error(input);
            arg.printUsage();
            cerr << arg.getProg() << ": " << error.what() << endl;
            throw error;
        }
    }

    FileSystem::Project inputFile(opts("Input").getString());
    Parser::GiD::Parser parserGiD(inputFile.getFilename());
    Data* smb = parserGiD.read();

    if (!opts("Solver").isNull()) {
        if (opts("Solver").getString() != smb->solver->getName()) {
            delete smb->solver;
            smb->solver = new SEMBA::Solver::Info(opts("Solver").getString());
        }
    }

    Cudg3d::Options solverDGTDOptions;
    solverDGTDOptions.set(opts);
        if (smb->solver->getName() == "cudg3d") {
        Cudg3d::Solver solver(smb);
        if (solverDGTDOptions.isRunSimulation()) {
            solver.run();
        }
    } else {
        throw std::logic_error(std::string("Invalid solver name ") +
                smb->solver->getName());
    }

    exit(EXIT_SUCCESS);
}
