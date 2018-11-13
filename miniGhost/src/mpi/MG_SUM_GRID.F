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

MODULE MG_SUM_GRID_MOD

   ! Purpose
   ! =======
   ! Coordination of the grid summation on the appropriate variable.

   USE MG_UTILS_MOD
   USE MG_ALLREDUCE_MOD

   IMPLICIT NONE

CONTAINS

   SUBROUTINE MG_SUM_GRID ( IVAR, GSUM, IERR )

      ! ---------------------
      ! Argument Declarations
      ! ---------------------

      INTEGER, INTENT(IN) :: &
         IVAR           ! Variable to be reduced.

      REAL(KIND=MG_REAL), INTENT(OUT) :: & 
         GSUM

      INTEGER, INTENT(OUT) :: &
         IERR           ! Return status

      ! ------------------
      ! Local Declarations
      ! ------------------

      ! ---------------------
      ! Executable Statements
      ! ---------------------

      IERR = 0

      CALL MG_ALLREDUCE_SUM ( IVAR, GRID(:,:,:,IVAR), GSUM, IERR )

      MG_PERF%NUM_SUMGRID = MG_PERF%NUM_SUMGRID + 1

   END SUBROUTINE MG_SUM_GRID

END MODULE MG_SUM_GRID_MOD
