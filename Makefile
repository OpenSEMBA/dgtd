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

# -- USAGE --------------------------------------------------------------------
# make target   = {debug, release}
#      compiler = {gnu}
# ==================== Default values =========================================
target = release
compiler = gnu
static = yes

CUDG3D_VERSION=\"0.12\"
# =============================================================================
DEFINES +=CUDG3D_VERSION=$(CUDG3D_VERSION)
ifeq ($(demo),yes)
	DEFINES += DEMO
endif
# -------------------- Paths to directories -----------------------------------
BUILD_DIR    = ./build/
OBJ_DIR      = ./obj/
SRC_DIR      = ./src/
EXTERNAL_DIR = ./external/

BIN_DIR = $(BUILD_DIR)bin/
LIB_DIR = $(BUILD_DIR)lib/
#===================== GNU Compiler ===========================================
ifeq ($(compiler),gnu)
	CC = gcc 
	CXX = g++
	CCFLAGS +=
	CXXFLAGS += -static -std=c++0x -fopenmp -pthread
endif # endif choosing the GNU compiler.
# ================= Optimization target =======================================
ifeq ($(target),debug)
	CXXFLAGS +=-O0 -g3 -Wall -Wno-write-strings # -Wconversion #-fprofile-arcs -ftest-coverage
	DEFINES  +=_DEBUG
endif
ifeq ($(target),release)
   	CXXFLAGS += -O2 
endif
# =============================================================================
# -------------------- RULES --------------------------------------------------
default: all
	@echo "======>>>>> Done <<<<<======"

all: cudg3d

cudg3d: opensemba check
	$(MAKE) -f ./src/apps/cudg3d/cudg3d.mk print
	$(MAKE) -f ./src/apps/cudg3d/cudg3d.mk

opensemba: check
	$(MAKE) --directory=$(EXTERNAL_DIR)opensemba/ -f Makefile opensemba
	-mkdir -p $(LIB_DIR)
	-cp -r $(EXTERNAL_DIR)opensemba/$(LIB_DIR)/* $(LIB_DIR)

test: opensemba check
	$(MAKE) -f $(SRC_DIR)apps/test/test.mk print
	$(MAKE) -f $(SRC_DIR)apps/test/test.mk

clean:
	-rm -rf $(OBJ_DIR)
	$(MAKE) --directory=$(EXTERNAL_DIR)opensemba/ -f Makefile clean

clobber: clean
	-rm -rf $(BUILD_DIR)
	$(MAKE) --directory=$(EXTERNAL_DIR)opensemba/ -f Makefile clobber

check:
ifneq ($(target),release) 
ifneq ($(target),debug) 
	@echo "Invalid build target."  
	@echo "Please use 'make target=debug|release'"  
	@exit 1
endif 
endif 
ifneq ($(compiler),gnu)  
	@echo "Invalid build compiler"  
	@echo "Please use 'make compiler= gnu'" 
	@exit 2
endif 

# Exports current variables when other makefiles are called.
export
