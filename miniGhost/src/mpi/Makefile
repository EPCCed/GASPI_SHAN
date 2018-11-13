#@header
# ************************************************************************
#
#      miniGhost: stencil computations with boundary exchange.
#              Copyright (2012) sandia corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Richard F. Barrett (rfbarre@sandia.gov) or
#                    Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
#@header

# Simple hand-tuned makefile, configured for the gnu compiler and MPICH.
# Modify as necessary for your environment.

#-----------------------------------------------------------------------

# Option: -D_MG_SERIAL
#GPI2_DIR = $(HOME)/GPI-2.openmpi-3.0.0.1S
#MPI_DIR=/sw/laki-SL6x/hlrs/mpi/openmpi/3.0.0-gnu-7.1.0

#MPI_DIR = /opt/cray/mpt/7.5.5/gni/mpich-gnu/5.1
#GPI2_DIR = /work/y07/y07/cse/gpi2/1.3.0
#GPI2_DIR = $(HOME)/GPI2-1.3.0-2018


PROTOCOL = -D_MG_MPI 

# Compilers
#FC=mpif90
#CC=mpicc
#FC=mpif90
FC=ftn
CC=cc
#CC=cc

CFLAGS = $(PROTOCOL)
# C main calling Fortran subroutine:
CFLAGS += -Df2c_ -g

FFLAGS = $(PROTOCOL)

# Variable precision: -D_INT8 and/or -D_REAL8.
FFLAGS += -D_MG_INT4 -D_MG_REAL8
#FFLAGS += -D_MG_INT4 -D_MG_REAL8 

CFLAGS = $(PROTOCOL)
# C main calling Fortran subroutine:
CFLAGS += -Df2c_ -g
#CFLAGS += -Df2c_ -DNUM_SOCKETS=2

# Optimization
# Next line PGI compiler
OPT_F = -O3 -g
#OPT_F = -O3

FFLAGS += $(OPT_F)
# Free-form Fortran source code:
# Next line PGI compiler:
#FFLAGS += -Mfree
# Next line Intel compiler:
FFLAGS += -ffree-form -ffree-line-length-none
# Next line Cray CCE compiler:
# FFLAGS += -f free -m 4

# Array bounds checking: (expensive!)
#FFLAGS += -Mbounds
# Compile to include checkpointing capability.
#FFLAGS += -D_MG_CHECKPT

LD=$(CC)
LDFLAGS=$(CFLAGS)
#LIBS=-lifcore
#LIBS=-lgfortran -lmpi_mpifh
LIBS=-lgfortran

include make_targets

# End makefile
