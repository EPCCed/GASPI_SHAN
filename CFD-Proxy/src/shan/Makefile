#NETCDF_DIR = $(HOME)/intel64/netcdf-3.6.3
#MPI_DIR= /opt/intel/impi/5.1.2.150
#CC = $(MPI_DIR)/intel64/bin/mpiicc
NETCDF_DIR = /sw/laki-SL6x/hlrs/tools/netcdf/4.3.3-gnu-4.9.1
MPI_DIR = /sw/laki-SL6x/hlrs/mpi/openmpi/3.0.0-gnu-7.1.0
GPI2_DIR = $(HOME)/GPI-2.openmpi-3.0.0.1S
SHAN_DIR=$(HOME)/bsc_workshop/shan
#GPI2_DIR = $(HOME)/GPI2-1.3.0-2018
CC = mpicc
#CC = cc

BIN += hybrid.f6.exe

CFLAGS += -Wall
CFLAGS += -Wextra
CFLAGS += -Wshadow
CFLAGS += -O3 -g 
CFLAGS += -std=c99
#CFLAGS += -DDEBUG 
#CFLAGS += -DUSE_COMM_OVERLAP -DUSE_SHAN
CFLAGS += -DUSE_SHAN -DUSE_COMM_OVERLAP

###############################################################################

INCLUDE_DIR += /usr/include
INCLUDE_DIR += $(SHAN_DIR)/include
INCLUDE_DIR += $(GPI2_DIR)/include
INCLUDE_DIR += $(NETCDF_DIR)/include
INCLUDE_DIR += $(MPI_DIR)/include
#INCLUDE_DIR += $(HOME)/Chaco-2.2/include

LIBRARY_DIR += $(NETCDF_DIR)/lib
LIBRARY_DIR += $(MPI_DIR)/lib
#LIBRARY_DIR += $(HOME)/Chaco-2.2/lib
LIBRARY_DIR += $(SHAN_DIR)/lib64
LIBRARY_DIR += $(GPI2_DIR)/lib64

LDFLAGS = -Wall
LDFLAGS += $(addprefix -L,$(LIBRARY_DIR))
CFLAGS += $(addprefix -I,$(INCLUDE_DIR))

OBJ += comm_data
OBJ += exchange_data_shan
OBJ += solver_data
OBJ += solver
OBJ += error_handling
OBJ += read_netcdf
OBJ += points_of_color
OBJ += rangelist
OBJ += comm_util
OBJ += gradients
OBJ += flux
OBJ += util

#LIB += ibverbs
#LIB += mpich
#LIB += mpi_mt
#LIB += chaco
LIB += netcdf
LIB += m
LIB += SHAN
LIB += GPI2-dbg
LIB += ibverbs
LIB += pthread
LIB += rt

###############################################################################

all: $(BIN)

%.o: %.c
	$(CC) $(CFLAGS) -c $<  -o $@

%.exe: %.o $(addsuffix .o, $(OBJ))
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(addprefix -l, $(LIB))

###############################################################################

.PHONY: clean objclean

objclean:
	rm -f *.o

clean: objclean
	rm -f $(BIN)
