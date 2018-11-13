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

# Option: -D_SERIAL
PROTOCOL = -D_MPI

# Compilers
#FC=mpif90
FC=ftn
#CC=mpicc
CC=cc $(PROTOCOL)

FFLAGS = $(PROTOCOL)
FFLAGS += -D_TIMER_MPI

# Variable precision: -D_INT8 and/or -D_REAL8.
FFLAGS += -D_MG_INT4 -D_MG_REAL8

# C main calling Fortran subroutine:
CFLAGS = -Df2c_

# Optimization
OPT_F = -O3

FFLAGS += $(OPT_F)
# Free-form Fortran source code:
FFLAGS += -ffree-form
# Array bounds checking: (expensive!)
#FFLAGS += -fbounds-check
# Compile to include checkpointing capability.
#FFLAGS += -D_MG_CHECKPT

LD=$(CC)
LDFLAGS=$(CFLAGS)
LIBS=

include make_targets

# End makefile
