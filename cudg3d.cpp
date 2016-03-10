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

#include "dgtd/Solver.h"

int main(int argc, const char *argv[]) {

    Argument::Parser args;
    args.args(argc, argv);
    Argument::OptionBase& input = args.addOption(
            new Argument::Option<std::string>("Input", 'i', "input"));
    args.addOption(
            (new Argument::Switch("Help", 'h', "help"))->defaultVal(false));
    Argument::Object opts = args.parseKnownArgs().first;
    if (opts("Input").isNull()) {
        if (opts("Help").getBool()) {
            args.printHelp();
            exit(EXIT_SUCCESS);
        } else {
            Argument::Error::Required error(input);
            args.printUsage();
            cerr << args.getProg() << ": " << error.what() << endl;
            throw error;
        }
    }

    FileSystem::Project inputFile(opts("Input").getString());
    Parser::GiD::Parser parserGiD(inputFile.getFilename());
    Data* smb = parserGiD.read();

    if (smb->solver == NULL) {
        if (opts("Solver").isNull()) {
            Argument::Error::Required error(solver);
            args.printUsage();
            std::cerr << args.getProg() << ": " << error.what() << std::endl;
            throw error;
        }
    } else if (!opts("Solver").isNull()) {
        if (opts("Solver").getString() != smb->solver->getName()) {
            delete smb->solver;
            smb->solver = new Solver::Info(opts("Solver").getString());
        }
    }

    Solver::DGTD::Options solverDGTDOptions;
    solverFDTDOptions.set(opts);
    isRunActivated_ = solverFDTDOptions.isRunSimulation();

    if (smb->solver->getName() == "ugrfdtd") {
        solver_ = new Solver::FDTD::Solver(smb, arg);
    } else if (smb->solver->getName() == "cudg3d") {

    } else {
        throw std::logic_error(std::string("Invalid solver name ") +
                smb->solver->getName());
    }

    exit(EXIT_SUCCESS);
}
