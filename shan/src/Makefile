#GPI2_DIR = $(HOME)/GPI-2.openmpi-3.0.0.1S
#MPI_DIR=/sw/laki-SL6x/hlrs/mpi/openmpi/3.0.0-gnu-7.1.0
#GPI2_DIR = /opt/GPI-2.openmpi-3.0.0
#MPI_DIR = /opt/cray/mpt/7.5.5/gni/mpich-gnu/5.1
#GPI2_DIR = /work/y07/y07/cse/gpi2/1.3.0
GPI2_DIR = $(HOME)/GPI2-1.3.0-2018

CC = cc
FC = ftn
#CC = mpicc

CFLAGS += -Wall
CFLAGS += -Wextra
CFLAGS += -Wshadow
CFLAGS += -O3 -g 
CFLAGS += -std=c99
CFLAGS += -openmp
CFLAGS += -DDEBUG 

#CFLAGS += -DGCC_EXTENSION -DUSE_MPI_SHARED_WIN
CFLAGS += -DGCC_EXTENSION -DUSE_MPI_SHARED_WIN


FFLAGS += -g -O3
# Free-form Fortran source code:
# Next line PGI compiler:
#FFLAGS += -Mfree
# Next line Intel compiler:
FFLAGS += -ffree-form -ffree-line-length-none
# Next line Cray CCE compiler:
# FFLAGS += -f free -m 4



###############################################################################

INCLUDE_DIR += /usr/include
INCLUDE_DIR += $(MPI_DIR)/include 
INCLUDE_DIR += $(GPI2_DIR)/include 
INCLUDE_DIR += ../include
INCLUDE_DIR += ./

LIBRARY_DIR += $(MPI_DIR)/lib
LIBRARY_DIR += $(GPI2_DIR)/lib64

LDFLAGS += $(addprefix -L,$(LIBRARY_DIR))
CFLAGS  += $(addprefix -I,$(INCLUDE_DIR))
CFLAGS  += -D_GNU_SOURCE


OBJ += f_shan
OBJ += shan_util
OBJ += shan_segment
OBJ += shan_comm
OBJ += shan_core
OBJ += shan_type
OBJ += shan_exchange
OBJ += gaspi_util


SHAN = libSHAN.a

#LIB += GPI2-dbg
#LIB += ibverbs
#LIB += mpich
#LIB += mpi_mt
#LIB += mpi


OBJS_C = $(addsuffix .o, $(OBJ))
OBJS_F = F_SHAN.o

###############################################################################

all: $(SHAN) fortran

shan: $(SHAN) 

fortran: $(OBJS_F)
	rm -f ../include/F_SHAN.mod
	cp f_shan.mod ../include

$(SHAN): $(OBJS_C)
	rm -f ../lib64/$(SHAN)
	ar crs ../lib64/$(SHAN) $(OBJS_C)

.F.o:
	$(FC) $(FFLAGS) $(F_INCLUDE) -c $<

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<



###############################################################################

.PHONY: clean objclean

objclean:
	rm -f *.o *.mod 

clean: objclean
