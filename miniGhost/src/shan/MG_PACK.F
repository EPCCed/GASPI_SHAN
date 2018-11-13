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
         I, J, K, C_IDX                 ! Counters


    INTEGER ::        &
         MSG_REQS_OFFSET, TYPE_ID
    
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


    TYPE_ID = IVAR - 1
    
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

       ! Pack

       C_IDX = IDX_NEIGHBORS(BACK) - 1
       CALL F_SHAN_COMM_NOTIFY_OR_WRITE(NEIGHBOR_HOOD_ID, SEGMENT_ID, &
            TYPE_ID, C_IDX)
    END IF

    ! Front boundary

    IF ( NEIGHBORS(FRONT) /= -1 ) THEN

       ! Pack
       C_IDX = IDX_NEIGHBORS(FRONT) - 1
       CALL F_SHAN_COMM_NOTIFY_OR_WRITE(NEIGHBOR_HOOD_ID, SEGMENT_ID, &
            TYPE_ID, C_IDX)

    END IF

    ! East boundary

    IF ( NEIGHBORS(EAST) /= -1 ) THEN

       ! Pack
       C_IDX = IDX_NEIGHBORS(EAST) - 1
       CALL F_SHAN_COMM_NOTIFY_OR_WRITE(NEIGHBOR_HOOD_ID, SEGMENT_ID, &
            TYPE_ID, C_IDX)
    END IF

    ! West boundary

    IF ( NEIGHBORS(WEST) /= -1 ) THEN

       ! Pack
       C_IDX = IDX_NEIGHBORS(WEST) - 1
       CALL F_SHAN_COMM_NOTIFY_OR_WRITE(NEIGHBOR_HOOD_ID, SEGMENT_ID, &
            TYPE_ID, C_IDX)

    END IF

    ! North boundary

    IF ( NEIGHBORS(NORTH) /= -1 ) THEN

       ! Pack
       C_IDX = IDX_NEIGHBORS(NORTH) - 1
       CALL F_SHAN_COMM_NOTIFY_OR_WRITE(NEIGHBOR_HOOD_ID, SEGMENT_ID, &
            TYPE_ID, C_IDX)

    END IF

    ! South boundary

    IF ( NEIGHBORS(SOUTH) /= -1 ) THEN

       ! Pack
       C_IDX = IDX_NEIGHBORS(SOUTH) - 1
       CALL F_SHAN_COMM_NOTIFY_OR_WRITE(NEIGHBOR_HOOD_ID, SEGMENT_ID, &
            TYPE_ID, C_IDX)

    END IF

    MG_PERF%TIME_PACK_PE = MG_PERF%TIME_PACK_PE + MG_TIMER() - TIME_START

    
  END SUBROUTINE MG_PACK_AND_SEND
   
END MODULE MG_PACK_MOD
