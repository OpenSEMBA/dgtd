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

# ============= Makefile for the cudg3d program ===============================
# ==================== Default variables ======================================
target = debug
compiler = gnu
order = 1
OUT = cudg3d
# Default paths.
BINDIR = bin/
OBJDIR = obj/
SRCDIR = src/
# =============================================================================
CXXFLAGS += -fopenmp

INCLUDES +=/usr/local/include
LIBRARIES +=/usr/lib64
LIBS +=

ifeq ($(target),debug)
	DEFINES +=_DEBUG
endif

DEFINES += ORDER_N=$(order) USE_OPENMP
# =============================================================================
# -------------------- Paths to directories -----------------------------------
DIR = $(SRC_CORE_DIR) $(SRC_PARSER_DIR) $(SRC_EXPORTER_DIR) \
 src/apps/cudg3d/ solver/dgtd/ \
 solver/dgtd/core/ solver/dgtd/integrator/ \
 solver/dgtd/DG/ solver/dgtd/DG/dispersives/ solver/dgtd/DG/sources/  
SOURCE_DIR = $(addprefix $(SRCDIR), ${DIR})

SRCS_CXX := $(shell find $(SOURCE_DIR) -maxdepth 1 -type f -name "*.cpp")
SRCS_CXX := $(filter-out $(EXCLUDE), $(SRCS_CXX)) 
OBJS_CXX := $(patsubst $(SRCDIR)%,$(OBJDIR)%,$(SRCS_CXX:.cpp=.o))

SRCS_C := $(shell find $(SOURCE_DIR) -maxdepth 1 -type f -name "*.c")
SRCS_C := $(filter-out $(EXCLUDE), $(SRCS_C)) 
OBJS_C := $(patsubst $(SRCDIR)%,$(OBJDIR)%,$(SRCS_C:.c=.o))

.PHONY: default clean clobber check

default: check cudg3d
	@echo "======================================================="
	@echo "             $(OUT) compilation finished               "
	@echo "======================================================="
	
clean:
	rm -rf *.err *.o *.d $(OBJDIR)

clobber: clean
	rm -rf $(BINDIR)

$(OBJDIR)%.o: $(SRCDIR)%.cpp
	@dirname $@ | xargs mkdir -p
	@echo "Compiling:" $@
	$(CXX) $(CXXFLAGS) $(addprefix -D, $(DEFINES)) $(addprefix -I,$(INCLUDES) ${SOURCE_DIR}) -c -o $@ $<
	
$(OBJDIR)%.o: $(SRCDIR)%.c
	@dirname $@ | xargs mkdir -p
	@echo "Compiling:" $@
	$(CC) $(CCFLAGS) $(addprefix -D, $(DEFINES)) $(addprefix -I,$(INCLUDES) ${SOURCE_DIR}) -c -o $@ $<

cudg3d: $(OBJS_CXX) $(OBJS_C)
	@mkdir -p $(BINDIR)
	@echo "Linking:" $@
	${CXX} $^ -o $(BINDIR)$(OUT) $(CXXFLAGS) \
	 $(addprefix -D, $(DEFINES)) \
	 $(addprefix -I, $(SOURCE_DIR)) $(addprefix -I, ${INCLUDES}) \
	 $(addprefix -L, ${LIBRARIES}) $(addprefix -l, ${LIBS})

check:
	@echo "======================================================="
	@echo "            ----- Compiling $(OUT) ------              "
	@echo "target:           " $(target)
	@echo "Compiler:         " $(compiler)
	@echo "C++ Compiler:     " `which $(CXX)`
	@echo "C++ Flags:        " $(CXXFLAGS)
	@echo "Defines:          " $(DEFINES)
	@echo "Polynomial order: " $(order)
	@echo "======================================================="
	@sleep 1
