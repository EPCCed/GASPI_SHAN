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

MODULE MG_UNPACK_SVAF_MOD

  ! Purpose
  ! =======
  ! Complete incoming non-blocking MPI messages, and unpack data
  ! into appropriate variables.

  USE MG_CONSTANTS_MOD
  USE MG_UTILS_MOD
  USE MG_GET_FACE_MOD
  USE MG_PROFILING_MOD

  IMPLICIT NONE

  ! Note: Using 0-based indexing on arrays.

CONTAINS

  SUBROUTINE MG_UNPACK_SVAF ( GRID, IVAR, IERR )

    ! -------------------------------------
    ! Unpack neighbor recvd data into GRID.
    ! -------------------------------------

    ! ---------------------
    ! Argument Declarations
    ! ---------------------

    INTEGER, INTENT(OUT) :: &
         IERR           ! Return status

    INTEGER, INTENT(IN) :: &
         IVAR   

    REAL(KIND=MG_REAL), DIMENSION(0:NX+1, 0:NY+1, 0:NZ+1), INTENT(INOUT) :: &
         GRID


    ! ------------------
    ! Local Declarations
    ! ------------------

    INTEGER ::        &
         IWHICH,                      &  ! MPI_Wait_any SRC process
         I, L,                        &  ! Counters
         NUM_RECVS_OUTSTANDING,       &  ! Keep track of progress
         NUM_SENDS_OUTSTANDING


#if defined _MG_MPI
    INTEGER ::     &
         ISTAT(MPI_STATUS_SIZE)
#endif

    REAL(KIND=MG_REAL8) ::         &
         TIME_START

    REAL(KIND=MG_REAL)    :: SZ_REAL
    
    ! ---------------------
    ! Executable Statements
    ! ---------------------

    IERR = 0

    IF ( NUMPES == 1 ) &
         RETURN

#if defined _MG_MPI

    TIME_START = MG_TIMER()

    NUM_SENDS_OUTSTANDING = NUM_SENDS
    NUM_RECVS_OUTSTANDING = NUM_RECVS

    ! Complete sends and receives. 
    ! If recv, unpack into user space (domain face).

    DO L = 1, NUM_RECVS + NUM_SENDS


       CALL MPI_WAITANY ( MAX_NUM_SENDS + MAX_NUM_RECVS, MSG_REQS, IWHICH, ISTAT, IERR )
       CALL MG_ASSERT ( IERR, 'MG_UNPACK_SVAF: MPI_WAITANY', ISTAT(MPI_ERROR) )

       IF ( IWHICH > MAX_NUM_RECVS ) THEN

          NUM_SENDS_OUTSTANDING = NUM_SENDS_OUTSTANDING - 1

       ELSE 

          CALL REMOTE_GET_FACE ( IWHICH, GRID, IERR )
          NUM_RECVS_OUTSTANDING = NUM_RECVS_OUTSTANDING - 1

       END IF
    END DO

    CALL MG_ASSERT ( NUM_SENDS_OUTSTANDING, 'MG_UNPACK_SVAF: NUM_SENDS_OUTSTANDING',  &
         NUM_SENDS )
    CALL MG_ASSERT ( NUM_RECVS_OUTSTANDING, 'MG_UNPACK_SVAF: NUM_RECVS_OUTSTANDING',  &
         NUM_RECVS )
    
    NUM_RECVS = 0
    NUM_SENDS = 0

    MG_PERF%TIME_UNPACK_PE = MG_PERF%TIME_UNPACK_PE + MG_TIMER() - &
         TIME_START

#endif _MG_MPI

    RETURN

  END SUBROUTINE MG_UNPACK_SVAF

END MODULE MG_UNPACK_SVAF_MOD
