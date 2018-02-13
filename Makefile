# -- USAGE --------------------------------------------------------------------
# make target   = {debug, release}
# ==================== Default values =========================================
target = debug
static = yes
# =============================================================================
# -------------------- Paths to directories -----------------------------------
BUILD_DIR    = ./build/
OBJ_DIR      = ./obj/
SRC_DIR      = ./src/
EXTERNAL_DIR = ./external/

BIN_DIR = $(BUILD_DIR)bin/
LIB_DIR = $(BUILD_DIR)lib/
#===================== GNU Compiler ===========================================
CC        = gcc 
CXX       = g++
CCFLAGS  +=
CXXFLAGS += -static -std=c++11 -fopenmp -pthread
# ================= Optimization target =======================================
ifeq ($(target),debug)
	CXXFLAGS +=-O0 -g3 -Wall
	DEFINES  +=_DEBUG
endif
ifeq ($(target),release)
   	CXXFLAGS += -O3 
endif
# =============================================================================
default: all
	@echo "======>>>>> Done <<<<<======"

all: cudg3d test

cudg3d: opensemba check
	$(MAKE) -f $(SRC_DIR)/apps/cudg3d/cudg3d.mk print
	$(MAKE) -f $(SRC_DIR)/apps/cudg3d/cudg3d.mk

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
	@echo "Please use 'make target=release' or 'make target=debug'"  
	@exit 1
endif 
endif 

# Exports current variables when other makefiles are called.
export