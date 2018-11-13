! ************************************************************************
!
!               miniGhost: stencil computations with boundary exchange.
!                 Copyright (2012) Sandia Corporation
!
! Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
! license for use of this work by or on behalf of the U.S. Government.
!
! This library is free software; you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation; either version 2.1 of the
! License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
! USA
! Questions? Contact Richard F. Barrett (rfbarre@sandia.gov) or
!                    Michael A. Heroux (maherou@sandia.gov)
!
! ************************************************************************
 
MODULE MG_OPTIONS_MOD

   ! Purpose
   ! =======
   ! Declaration of miniGhost options.

   USE MG_CONSTANTS_MOD

   IMPLICIT NONE

   INTEGER :: SCALING

   INTEGER, PARAMETER ::      &
      SCALING_STRONG   = 1,   &
      SCALING_WEAK     = 2           ! Default

   INTEGER :: COMM_METHOD

   INTEGER, PARAMETER ::       &
      COMM_METHOD_SVAF  = 11

   INTEGER :: STENCIL

   INTEGER, PARAMETER ::          &

      STENCIL_NONE         = 20,  &
      STENCIL_2D5PT        = 21,  &  ! Default
      STENCIL_3D7PT        = 23

   ! ---------
   ! Variables
   ! ---------

   INTEGER ::                      &

      NUM_SPIKES,                  &  ! Number of spikes to be inserted, one each NUM_TSTEPS

      NPX, NPY, NPZ,               &  ! Number of processes in direction X, Y, Z

      NX, NY, NZ,                  &  ! Local domain dimensions.

      NUM_VARS,                    &  ! Number of variables defined.
      NUM_TSTEPS,                  &  ! Number of time steps to be executed.
      PERCENT_SUM,                 &  ! Requested percentage of variables to be globally summed.

      DEBUG_GRID,                  &  !
      REPORT_DIFFUSION,            &  !
      REPORT_PERF                     !


END MODULE MG_OPTIONS_MOD
