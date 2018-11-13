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

MODULE MG_FLUX_ACCUMULATE_MOD

   USE MG_UTILS_MOD
   USE MG_PROFILING_MOD

   IMPLICIT NONE

   ! Keeps track of heat dissapation out of physical domain.

CONTAINS
   
   SUBROUTINE MG_FLUX_ACCUMULATE ( IVAR, GRID, IERR )
   
      ! ---------------------
      ! Argument Declarations
      ! ---------------------

      INTEGER, INTENT(IN) :: &
         IVAR

      INTEGER, INTENT(OUT) :: &
         IERR           ! Return status

      REAL(KIND=MG_REAL), DIMENSION(0:NX+1, 0:NY+1, 0:NZ+1), INTENT(INOUT) :: &
         GRID

      ! ------------------
      ! Local Declarations
      ! ------------------

      INTEGER ::         &
         I, J, K            ! Counters

      REAL(KIND=MG_REAL) :: &
         DIVISOR

      ! ---------------------
      ! Executable Statements
      ! ---------------------
   
      IERR = 0

      IF ( STENCIL == STENCIL_2D5PT ) THEN 
         DIVISOR = 5.0
      ELSE IF ( STENCIL == STENCIL_3D7PT ) THEN
         DIVISOR = 7.0
      END IF

      IF ( MYPX == 0 )  THEN
         DO K = 1, NZ
            DO J = 1, NY
               FLUX_OUT(IVAR) = FLUX_OUT(IVAR) + ( ( GRID(1,J,K)) / DIVISOR )
            END DO
         END DO
      END IF

      IF ( MYPX == ( NPX - 1 ) )  THEN
         DO K = 1, NZ
            DO J = 1, NY
               FLUX_OUT(IVAR) = FLUX_OUT(IVAR) + ( ( GRID(NX,J,K)) / DIVISOR )
            END DO
         END DO
      END IF

      IF ( MYPY == 0 )  THEN
         DO K = 1, NZ
            DO I = 1, NX
               FLUX_OUT(IVAR) = FLUX_OUT(IVAR) + ( ( GRID(I,1,K)) / DIVISOR )
            END DO
         END DO
      END IF

      IF ( MYPY == ( NPY - 1 ) ) THEN
         DO K = 1, NZ
            DO I = 1, NX
               FLUX_OUT(IVAR) = FLUX_OUT(IVAR) + ( ( GRID(I,NY,K)) / DIVISOR )
            END DO
         END DO
      END IF

      IF ( MYPZ == 0 )  THEN
         DO J = 1, NY
            DO I = 1, NX
               FLUX_OUT(IVAR) = FLUX_OUT(IVAR) + ( ( GRID(I,J,1)) / DIVISOR )
            END DO
         END DO
      END IF

      IF ( MYPZ == ( NPZ - 1 ) ) THEN
         DO J = 1, NY
            DO I = 1, NX
               FLUX_OUT(IVAR) = FLUX_OUT(IVAR) + ( ( GRID(I,J,NZ)) / DIVISOR )
            END DO
         END DO
      END IF

   END SUBROUTINE MG_FLUX_ACCUMULATE
   
END MODULE MG_FLUX_ACCUMULATE_MOD
