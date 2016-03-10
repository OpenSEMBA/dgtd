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

#include "solver/dgtd/SolverDGTD.h"

#define APP_NAME "cudg3d"

int main(int argc, const char *argv[]) {
    Arguments arg(argc, argv);
    arg.printWelcomeMessage(string(APP_NAME), string(APP_VERSION));
    arg.printInfo();

    SmbData* smb;
    ProjectFile inputFile(arg.getFilename());
    ParserGiD parserGiD(inputFile);
    smb = parserGiD.read();
    smb->solverOptions->set(arg);

    {
        SolverDGTD cudg3d(smb);
        cudg3d.run();
    }

    arg.printGoodbyeMessage(string(APP_NAME));
    exit(EXIT_SUCCESS);
}
