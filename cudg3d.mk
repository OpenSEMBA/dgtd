# OpenSEMBA
# Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
#                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
#                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
#                    Daniel Mateos Romero            (damarro@semba.guru)
#
# This file is part of OpenSEMBA.
#
# OpenSEMBA is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with OpenSEMBA. If not, see <http://www.gnu.org/licenses/>.
# =============================================================================
OUT = cudg3d
# =============================================================================
#ifeq ($(FFTW3_SUPPORT),yes)
#	DEFINES +=FFTW3_SUPPORT
#	LIBS += fftw3
#endif
# -------------------- Paths to directories -----------------------------------
SRC_DIRS = $(shell find $(SRC_DIR)/apps/cudg3d/ -type d)

SRCS_CXX := $(shell find $(SRC_DIRS) -maxdepth 1 -type f -name "*.cpp")
OBJS_CXX := $(addprefix $(OBJ_DIR), $(SRCS_CXX:.cpp=.o))

# =============================================================================
LIBS += gidpost opensemba
INCLUDES += $(LIB_DIR)gidpost/include/ \
			$(LIB_DIR)opensemba/include/core/ $(LIB_DIR)opensemba/include/ \
			$(SRC_DIR)/apps/cudg3d/
LIBRARIES += $(LIB_DIR)gidpost/lib/ $(LIB_DIR)opensemba/lib/
# =============================================================================
.PHONY: default check

default: print $(OUT)
	@echo "======================================================="
	@echo "             $(OUT) compilation finished               "
	@echo "======================================================="

$(OBJ_DIR)%.o: %.cpp
	@dirname $@ | xargs mkdir -p
	@echo "Compiling:" $@
	$(CXX) $(CXXFLAGS) $(addprefix -D, $(DEFINES)) $(addprefix -I,$(INCLUDES) ${SOURCE_DIR}) -c -o $@ $<

cudg3d: $(OBJS_CXX) 
	@mkdir -p $(BIN_DIR)
	@echo "Linking:" $@
	${CXX} $^ -o $(BIN_DIR)$(OUT) $(CXXFLAGS) \
	 $(addprefix -D, $(DEFINES)) \
	 $(addprefix -I, ${INCLUDES}) \
	 $(addprefix -L, ${LIBRARIES}) \
	 $(addprefix -l, ${LIBS})

print:
	@echo "======================================================="
	@echo "            ----- Compiling $(OUT) ------              "
	@echo "target:           " $(target)
	@echo "Compiler:         " $(compiler)
	@echo "C++ Compiler:     " `which $(CXX)`
	@echo "C++ Flags:        " $(CXXFLAGS)
	@echo "Defines:          " $(DEFINES)
	@echo "Polynomial order: " $(order)
	@echo "======================================================="
