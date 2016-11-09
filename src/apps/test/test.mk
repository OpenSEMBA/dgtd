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

OUT = test

TEST_CORE_CELL          = yes#

# =============================================================================
SRC_APP_DIR = $(SRC_DIR)apps/test/
# =============================================================================
# --- Core ---
ifeq ($(TEST_CORE_CELL),yes)
	SRC_CORE_CELL_DIRS     := $(shell find $(SRC_DIR)core/cell/ -type d)
	SRC_CORE_CELL_TESTS_DIRS := $(SRC_CORE_CELL_DIRS) \
							    $(shell find $(SRC_APP_DIR)core/cell/ -type d)
endif

SRC_CORE_TESTS_DIRS = $(SRC_CORE_CELL_TESTS_DIRS) 
# ----- Gathers sources ----
SRC_DIRS := $(SRC_APP_DIR) $(SRC_CORE_TESTS_DIRS) 

SRCS_CXX := $(shell find $(SRC_DIRS) -maxdepth 1 -type f -name "*.cpp")
OBJS_CXX := $(addprefix $(OBJ_DIR), $(SRCS_CXX:.cpp=.o))
# =============================================================================
LIBS      += gtest gidpost opensemba 
INCLUDES  += $(SRC_DIR)core/ \
             $(LIB_DIR)gidpost/include/ \
			 $(LIB_DIR)opensemba/include/core/ \
			 $(LIB_DIR)opensemba/include/
LIBRARIES += $(LIB_DIR)gidpost/lib/ $(LIB_DIR)opensemba/lib/

OBJS_LIB  += $(LIB_DIR)opensemba/lib/libopensemba.a \
			 $(LIB_DIR)gidpost/lib/libgidpost.a \
# =============================================================================
.PHONY: default print

default: $(OUT)
	@echo "======================================================="
	@echo "           $(OUT) compilation finished"
	@echo "======================================================="

$(OBJ_DIR)%.o: %.cpp
	@dirname $@ | xargs mkdir -p
	@echo "Compiling:" $@
	$(CXX) $(CXXFLAGS) $(addprefix -D, $(DEFINES)) $(addprefix -I,$(INCLUDES)) -c -o $@ $<

$(BIN_DIR)$(OUT): $(OBJS_CXX) $(OBJS_LIB)
	@mkdir -p $(BIN_DIR)
	@echo "Linking:" $@
	${CXX} $^ \
	-o $@ $(CXXFLAGS) \
	$(addprefix -D, $(DEFINES)) \
	$(addprefix -I, ${INCLUDES}) \
	$(addprefix -L, ${LIBRARIES}) \
	$(addprefix -l, ${LIBS})

$(OUT): print $(BIN_DIR)$(OUT)

print:
	@echo "======================================================="
	@echo "         ----- Compiling $(OUT) ------        "
	@echo "Target:           " $(target)
	@echo "Compiler:         " $(compiler)
	@echo "C++ Compiler:     " `which $(CXX)`
	@echo "C++ Flags:        " $(CXXFLAGS)
	@echo "Defines:          " $(DEFINES)
	@echo "======================================================="

# ------------------------------- END ----------------------------------------
