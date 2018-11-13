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

MODULE MG_PACK_MOD

  ! Purpose
  ! =======
  ! Pack boundary data into send buffer for subsequent transmission,
  ! using the SVAF boundary exchange strategy.
  ! Note that ghost space is sent as well (i.e. rows 0 and NX+1),
  ! necessary for stencils with diagonals and irrelevant for stencils
  ! without diagonals. Could be a separate routine, but no matter.

  USE MG_CONSTANTS_MOD
  USE MG_UTILS_MOD
  USE MG_PROFILING_MOD

  IMPLICIT NONE

CONTAINS

  SUBROUTINE MG_PACK_AND_SEND ( GRID, IVAR, IERR )

    ! ---------------------
    ! Argument Declarations
    ! ---------------------

    INTEGER, INTENT(OUT) :: &
         IERR           ! Return status

    INTEGER, INTENT(IN) :: &
         IVAR           

    REAL(KIND=MG_REAL), DIMENSION(0:NX+1, 0:NY+1, 0:NZ+1), INTENT(IN) :: &
         GRID

    ! ------------------
    ! Local Declarations
    ! ------------------

    INTEGER ::        &
         I, J, K              ! Counters


    INTEGER ::        &
         MSG_REQS_OFFSET
    
    REAL(KIND=MG_REAL8) ::      &
         TIME_START,              &
         TIME_START_DIR
        
    REAL(KIND=MG_REAL)    :: SZ_REAL

    ! ---------------------
    ! Executable Statements
    ! ---------------------

    IERR = 0

    IF ( NUMPES == 1 ) &
         RETURN


    TIME_START = MG_TIMER()
    

    NUM_SENDS = 0
    MSG_REQS_OFFSET = MAX_NUM_RECVS


    COUNT_SEND_BACK  = (NY+2)*(NX+2)
    COUNT_SEND_FRONT = (NY+2)*(NX+2)
    COUNT_SEND_EAST  = (NZ+2)*(NY+2)
    COUNT_SEND_WEST  = (NZ+2)*(NY+2)
    COUNT_SEND_NORTH = (NZ+2)*(NX+2)
    COUNT_SEND_SOUTH = (NZ+2)*(NX+2)
    
    ! Back boundary

    IF ( NEIGHBORS(BACK) /= -1) THEN


       !$OMP PARALLEL DO
       DO J = 0, NY+1
          DO I = 0, NX+1
             SEND_BUFFER_BACK(J*(NX+2) + I + 1) = GRID ( I, J, 1 )
          END DO
       END DO
       !$OMP END PARALLEL DO


       ! Send
       NUM_SENDS = NUM_SENDS + 1
       MSG_REQS_OFFSET = MAX_NUM_RECVS + BACK
       CALL MPI_ISEND ( SEND_BUFFER_BACK, COUNT_SEND_BACK, MG_MPI_REAL, &
            NEIGHBORS(BACK), MSG_TAGS(BACK),                &
            MPI_COMM_MG, MSG_REQS(MSG_REQS_OFFSET), IERR )

    END IF

    ! Front boundary

    IF ( NEIGHBORS(FRONT) /= -1 ) THEN


       !$OMP PARALLEL DO
       DO J = 0, NY+1
          DO I = 0, NX+1
             SEND_BUFFER_FRONT(J*(NX+2) + I + 1) = GRID ( I, J, NZ )
          END DO
       END DO


       !$OMP END PARALLEL DO

       NUM_SENDS = NUM_SENDS + 1
       MSG_REQS_OFFSET = MAX_NUM_RECVS + FRONT
       CALL MPI_ISEND ( SEND_BUFFER_FRONT, COUNT_SEND_FRONT, MG_MPI_REAL, &
            NEIGHBORS(FRONT), MSG_TAGS(FRONT),                &
            MPI_COMM_MG, MSG_REQS(MSG_REQS_OFFSET), IERR )

    END IF

    ! East boundary

    IF ( NEIGHBORS(EAST) /= -1 ) THEN

       !$OMP PARALLEL DO
       DO K = 0, NZ+1
          DO J = 0, NY+1
             SEND_BUFFER_EAST(K*(NY+2) + J + 1) = GRID ( NX, J, K )
          END DO
       END DO


       !$OMP END PARALLEL DO

       ! Send
       NUM_SENDS = NUM_SENDS + 1
       MSG_REQS_OFFSET = MAX_NUM_RECVS + EAST
       CALL MPI_ISEND ( SEND_BUFFER_EAST, COUNT_SEND_EAST, MG_MPI_REAL, &
            NEIGHBORS(EAST), MSG_TAGS(EAST),                &
            MPI_COMM_MG, MSG_REQS(MSG_REQS_OFFSET), IERR )

    END IF

    ! West boundary

    IF ( NEIGHBORS(WEST) /= -1 ) THEN

       !$OMP PARALLEL DO
       DO K = 0, NZ+1
          DO J = 0, NY+1
             SEND_BUFFER_WEST(K*(NY+2) + J + 1) = GRID ( 1, J, K )
          END DO
       END DO


       !$OMP END PARALLEL DO

       ! Send
       NUM_SENDS = NUM_SENDS + 1
       MSG_REQS_OFFSET = MAX_NUM_RECVS + WEST
       CALL MPI_ISEND ( SEND_BUFFER_WEST, COUNT_SEND_WEST, MG_MPI_REAL, &
            NEIGHBORS(WEST), MSG_TAGS(WEST),                &
            MPI_COMM_MG, MSG_REQS(MSG_REQS_OFFSET), IERR )       


    END IF

    ! North boundary

    IF ( NEIGHBORS(NORTH) /= -1 ) THEN

       !$OMP PARALLEL DO
       DO K = 0, NZ+1
          DO I = 0, NX+1
             SEND_BUFFER_NORTH(K*(NX+2) + I + 1) = GRID ( I, NY, K )
          END DO
       END DO


       !$OMP END PARALLEL DO

       ! Send
       NUM_SENDS = NUM_SENDS + 1
       MSG_REQS_OFFSET = MAX_NUM_RECVS + NORTH
       CALL MPI_ISEND ( SEND_BUFFER_NORTH, COUNT_SEND_NORTH, MG_MPI_REAL, &
            NEIGHBORS(NORTH), MSG_TAGS(NORTH),                &
            MPI_COMM_MG, MSG_REQS(MSG_REQS_OFFSET), IERR )

    END IF

    ! South boundary

    IF ( NEIGHBORS(SOUTH) /= -1 ) THEN

       !$OMP PARALLEL DO
       DO K = 0, NZ+1
          DO I = 0, NX+1
             SEND_BUFFER_SOUTH(K*(NX+2) + I + 1) = GRID ( I, 1, K )
          END DO
       END DO


       !$OMP END PARALLEL DO

       NUM_SENDS = NUM_SENDS + 1
       MSG_REQS_OFFSET = MAX_NUM_RECVS + SOUTH
       CALL MPI_ISEND ( SEND_BUFFER_SOUTH, COUNT_SEND_SOUTH, MG_MPI_REAL, &
            NEIGHBORS(SOUTH), MSG_TAGS(SOUTH),                &
            MPI_COMM_MG, MSG_REQS(MSG_REQS_OFFSET), IERR )

    END IF

    MG_PERF%TIME_PACK_PE = MG_PERF%TIME_PACK_PE + MG_TIMER() - TIME_START

    
  END SUBROUTINE MG_PACK_AND_SEND
   
END MODULE MG_PACK_MOD
