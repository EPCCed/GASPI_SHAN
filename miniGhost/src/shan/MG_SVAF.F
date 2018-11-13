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

MODULE MG_SVAF_MOD

  ! Purpose
  ! =======
  ! SVAF boundary exchange for 5- and 7-point stencils:
  !    pack, send, receive, unpack.

  USE MG_CONSTANTS_MOD
  USE MG_UTILS_MOD
  USE MG_PACK_MOD
  USE MG_UNPACK_SVAF_MOD

  IMPLICIT NONE

  ! Note: Using 0-based indexing on arrays.

CONTAINS

  SUBROUTINE MG_SVAF ( IVAR, IERR )

    ! -------------------------------------------------------
    ! Pack all variables for single message to each neighbor.
    ! -------------------------------------------------------

    ! ---------------------
    ! Argument Declarations
    ! ---------------------

    INTEGER, INTENT(IN) :: &
         IVAR           ! Return status

    INTEGER, INTENT(OUT) :: &
         IERR           ! Return status

    ! ------------------
    ! Local Declarations
    ! ------------------

    INTEGER :: &
         ICOUNT,                      &  ! Message size.
         ISTAT,                       &  ! MPI return status.
         IWHICH,                      &  ! MPI_Wait_any SRC process
         I, J, K, L, M, N,            &  ! Counters
         NUM_RECVS,                   &
         NUM_RECVS_OUTSTANDING,       &
         MSG_REQ_OFFSET,              &
         OFFSET,                      &
         OFFSET_LAST,                 &
         RECV_BUFFER_SIZE,            &
         SEND_BUFFER_SIZE,            &
         MSG_REQ(12)

    ! ---------------------
    ! Executable Statements
    ! ---------------------

    IERR = 0

    IF (IVAR > 40) THEN
       WRITE(*,*) '** Error ** COMM_MG_SVAF: Unknown variable number ', IVAR, '.'
    ELSE

       ! -------------
       ! Pack and send
       ! -------------
       CALL MG_PACK_AND_SEND ( GRID(:,:,:,IVAR), IVAR, IERR )

       ! ---------------
       ! Wait and unpack
       ! ---------------
       CALL MG_UNPACK_SVAF ( GRID(:,:,:,IVAR), IVAR, IERR )

    ENDIF

  END SUBROUTINE MG_SVAF

END MODULE MG_SVAF_MOD
