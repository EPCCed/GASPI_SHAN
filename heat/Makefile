# Compilers
#CXX=g++
CXX=CC
#MCXX=mcxx
#MPICXX=mpicxx
MPICXX=CC
WRAPPERS=I_MPI_CXX=$(MCXX) MPICH_CXX=$(MCXX) OMPI_CXX=$(MCXX)

# Set the default block size (elements)
BSX?=16
BSY?=$(BSX)
BS=$(BSX)x$(BSY)bs

# Preprocessor flags
CPPFLAGS=-Isrc -DBSX=$(BSX) -DBSY=$(BSY)

# Compiler flags
CXXFLAGS=-O3 -std=c++11
MCXXFLAGS=--ompss-2 $(CXXFLAGS)

# Linker flags
#LDFLAGS=-lrt -lm -libverbs
LDFLAGS=-lrt -lm

GPI2_DIR=$(HOME)/GPI2-1.3.0-2018
#GPI2_DIR=/work/y07/y07/cse/gpi2/1.3.0
#GPI2_DIR=/work/y07/y07/cse/gpi2/1.3.0-shared-notifications
SHAN_DIR=$(HOME)/bsc_workshop/shan

#GPI2_DIR = $(HOME)/GPI-2.openmpi-3.0.0.1S
#MPI_DIR=/sw/laki-SL6x/hlrs/mpi/openmpi/3.0.0-gnu-7.1.0

# Interoperability flags
SHAN_CPPFLAGS = -I$(SHAN_DIR)/include  -I$(GPI2_DIR)/include 
SHAN_LDFLAGS = -L$(SHAN_DIR)/lib64 -lSHAN -L$(GPI2_DIR)/lib64 -lGPI2-dbg

MPI_INTEROP_CPPFLAGS=-I$(MPI_INTEROP_HOME)/include
MPI_INTEROP_LDFLAGS=-L$(MPI_INTEROP_HOME)/lib -L$(MPI_INTEROP_HOME)/lib64 -l:libmpiompss-interop.a

# List of programs
#BIN=	\
#	01.heat_seq.$(BS).bin \
#	02.heat_ompss.$(BS).bin \
#	03.heat_mpi.$(BS).bin \
#	04.heat_mpi_ompss_forkjoin.$(BS).bin \
#	05.heat_mpi_ompss_tasks.$(BS).bin \
#

BIN=	\
	08.heat_mpi_nbuffer.$(BS).bin \
	09.heat_mpi_nbuffer.$(BS).bin \
	010.heat_mpi.$(BS).bin \
	011.heat_gaspi_shan.$(BS).bin

ifdef MPI_INTEROP_HOME
#BIN+=06.heat_mpi_ompss_tasks_interop.$(BS).bin
endif

# Sources
SMP_SRC=src/common/misc.cpp src/smp/main.cpp
MPI_SRC=src/common/misc.cpp src/mpi/main.cpp

all: $(BIN)

01.heat_seq.$(BS).bin: $(SMP_SRC) src/smp/01.solver_seq.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

02.heat_ompss.$(BS).bin: $(SMP_SRC) src/smp/02.solver_ompss.cpp
	$(MCXX) $(CPPFLAGS) $(MCXXFLAGS) -o $@ $^ $(LDFLAGS)

03.heat_mpi.$(BS).bin: $(MPI_SRC) src/mpi/03.solver_mpi.cpp
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS)

04.heat_mpi_ompss_forkjoin.$(BS).bin: $(MPI_SRC) src/mpi/04.solver_mpi_ompss_forkjoin.cpp
	$(WRAPPERS) $(MPICXX) $(CPPFLAGS) $(MCXXFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS)

05.heat_mpi_ompss_tasks.$(BS).bin: $(MPI_SRC) src/mpi/05.06.solver_mpi_ompss_tasks.cpp
	$(WRAPPERS) $(MPICXX) $(CPPFLAGS) $(MCXXFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS)

06.heat_mpi_ompss_tasks_interop.$(BS).bin: $(MPI_SRC) src/mpi/05.06.solver_mpi_ompss_tasks.cpp
	$(WRAPPERS) $(MPICXX) -DINTEROPERABILITY $(CPPFLAGS) $(MCXXFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS) $(MPI_INTEROP_LDFLAGS)

ifdef MPI_INTEROP_HOME
07.heat_mpi_ompss_tasks_interop_async.$(BS).bin: $(MPI_SRC) src/mpi/07.solver_mpi_ompss_tasks_async.cpp
	$(WRAPPERS) $(MPICXX) -DINTEROPERABILITY $(CPPFLAGS) $(MPI_INTEROP_CPPFLAGS) $(MCXXFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS) $(MPI_INTEROP_LDFLAGS)
endif

08.heat_mpi_nbuffer.$(BS).bin: $(MPI_SRC) src/mpi/08.solver_mpi_nbuffer.cpp
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS)

09.heat_mpi_nbuffer.$(BS).bin: $(MPI_SRC) src/mpi/09.solver_mpi_nbuffer.cpp
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS)

010.heat_mpi.$(BS).bin: $(MPI_SRC) src/mpi/010.solver_mpi.cpp
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS)

011.heat_gaspi_shan.$(BS).bin: $(MPI_SRC) src/mpi/011.solver_gaspi_shan.cpp
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) -DUSE_SHAN $(SHAN_CPPFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS) $(SHAN_LDFLAGS)


check: all
	@./scripts/run-tests.sh

clean:
	rm -f *.o *.bin
