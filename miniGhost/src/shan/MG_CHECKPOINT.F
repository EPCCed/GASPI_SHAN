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

MODULE MG_CHECKPOINT_MOD

#if defined _MG_SERIAL

CONTAINS

   SUBROUTINE MG_CHECKPOINT ( IERR )

      ! Current implementation based on MPI IO functionality.

      IMPLICIT NONE

      INTEGER IERR                    ! Return status.

      RETURN

   END SUBROUTINE MG_CHECKPOINT

#elif defined _MG_MPI

   ! Procedures included:
   !
   !   CHECKPOINT

   USE MG_CONSTANTS_MOD
   USE MG_OPTIONS_MOD

   IMPLICIT NONE

   INTEGER(KIND=MG_INT) :: &
         CP_INTERVAL

   CHARACTER*(1024) :: &
         CP_FILE

   INTEGER(KIND=MPI_ADDRESS_KIND) ::  &
      SIZEOFCHAR,       &
      SIZEOFINT,        &
      SIZEOFREAL8

   INTEGER(KIND=MG_INT), PARAMETER :: MAX_GRIDS = 40

   TYPE :: CMDLINE_OPTIONS_TYPE

      INTEGER(KIND=MG_INT) ::  REPORT_DIFFUSION ! boolean
      INTEGER(KIND=MG_INT) ::  SCALING
      INTEGER(KIND=MG_INT) ::  COMM_METHOD
      INTEGER(KIND=MG_INT) ::  STENCIL
      INTEGER(KIND=MG_INT) ::  NUM_SPIKES
      INTEGER(KIND=MG_INT) ::  NPX
      INTEGER(KIND=MG_INT) ::  NPY
      INTEGER(KIND=MG_INT) ::  NPZ
      INTEGER(KIND=MG_INT) ::  NX
      INTEGER(KIND=MG_INT) ::  NY
      INTEGER(KIND=MG_INT) ::  NZ
      INTEGER(KIND=MG_INT) ::  NUM_VARS
      INTEGER(KIND=MG_INT) ::  NUM_TSTEPS
      INTEGER(KIND=MG_INT) ::  PERCENT_SUM
      INTEGER(KIND=MG_INT) ::  CP_INTERVAL
      CHARACTER(1024)      ::  CP_FILE

   END TYPE CMDLINE_OPTIONS_TYPE

   TYPE :: PE_COORDS_TYPE
      INTEGER(KIND=MG_INT) ::  MY_GLOBAL_NX_START
      INTEGER(KIND=MG_INT) ::  MY_GLOBAL_NY_START
      INTEGER(KIND=MG_INT) ::  MY_GLOBAL_NZ_START
      INTEGER(KIND=MG_INT) ::  MY_GLOBAL_NX_END
      INTEGER(KIND=MG_INT) ::  MY_GLOBAL_NY_END
      INTEGER(KIND=MG_INT) ::  MY_GLOBAL_NZ_END
      INTEGER(KIND=MG_INT) ::  NUM_NEIGHS
      INTEGER(KIND=MG_INT) ::  NUM_SUM_GRID
      INTEGER(KIND=MG_INT) ::  MYPE
      INTEGER(KIND=MG_INT) ::  MYPX
      INTEGER(KIND=MG_INT) ::  MYPY
      INTEGER(KIND=MG_INT) ::  MYPZ
   END TYPE PE_COORDS_TYPE

   TYPE :: TSHEADER_TYPE

      INTEGER(KIND=MG_INT)  ::  TSTEP
      INTEGER(KIND=MG_REAL) ::  GSUM_OLD(MAX_GRIDS)

   END TYPE TSHEADER_TYPE

   INTEGER CP_CMDLINE_TYPE
   INTEGER CP_GRIDSTOSUM_TYPE
   INTEGER CP_PECOORDS_TYPE
   INTEGER CP_PECOORDSARRAY_TYPE
   INTEGER CP_GSUMOLD_TYPE
   INTEGER CP_TSHEADER_TYPE
   INTEGER CP_NOGHOST_TYPE
   INTEGER CP_TSGRID_TYPE

   INTEGER(KIND=MPI_ADDRESS_KIND) CP_CMDLINE_EXTENT
   INTEGER(KIND=MPI_ADDRESS_KIND) CP_GRIDSTOSUM_EXTENT
   INTEGER(KIND=MPI_ADDRESS_KIND) CP_PECOORDS_EXTENT
   INTEGER(KIND=MPI_ADDRESS_KIND) CP_PECOORDSARRAY_EXTENT
   INTEGER(KIND=MPI_ADDRESS_KIND) CP_GSUMOLD_EXTENT
   INTEGER(KIND=MPI_ADDRESS_KIND) CP_TSHEADER_EXTENT
   INTEGER(KIND=MPI_ADDRESS_KIND) CP_NOGHOST_EXTENT
   INTEGER(KIND=MPI_ADDRESS_KIND) CP_TSGRID_EXTENT
   INTEGER(KIND=MPI_ADDRESS_KIND) CP_LB

   INTEGER cp_filehandle
CONTAINS
!  =================================================================================

   SUBROUTINE MG_CHECKPOINT ( IERR )

      IMPLICIT NONE

      INTEGER IERR                    ! Return status.

      ! ---------------
      ! Local Variables
      ! ---------------

      INTEGER(KIND=MG_INT) IVAR, cp_num, view_disp

      ! ---------------------
      ! Executable Statements
      ! ---------------------

      IERR = 0

#if defined _MG_CHECKPT

      IF ( CP_INTERVAL == 0 ) THEN
         ! A checkpoint interval of zero means no checkpointing
         RETURN
      END IF

      IF ( MOD(TSTEP,CP_INTERVAL) /= 0 ) THEN
         RETURN
      END IF

      CALL MPI_FILE_OPEN(MPI_COMM_WORLD, CP_FILE, MPI_MODE_CREATE+MPI_MODE_RDWR, MPI_INFO_NULL, cp_filehandle, IERR)

      IF ( TSTEP == CP_INTERVAL ) THEN
         ! This is the first checkpoint.
         ! Initialize the MPI types.
         CALL MG_CHECKPOINT_CREATE_TYPES ( IERR )
         ! Initialize the checkpoint file.
         CALL MG_CHECKPOINT_INIT ( IERR )
      END IF

      ! Before writing the grid vars, write a small header.
      CALL MG_CHECKPOINT_TSTEP_HEADER ( IERR )

      DO IVAR = 1, NUM_VARS
         CALL MG_CHECKPOINT_TSTEP_VAR ( IVAR, IERR )
      END DO

      IF ( TSTEP == CP_INTERVAL ) THEN
         CALL MG_CHECKPOINT_VERIFY ( IERR )
      END IF

      CALL MPI_FILE_CLOSE(cp_filehandle, IERR)

#endif _MG_CHECKPT

      RETURN

   END SUBROUTINE MG_CHECKPOINT

#if defined _MG_CHECKPT

!  ===================================================================================

!  ===================================================================================

   SUBROUTINE MG_CHECKPOINT_CREATE_TYPES ( IERR )

      IMPLICIT NONE

      ! ---------------------
      ! Argument Declarations
      ! ---------------------

      INTEGER IERR                    ! Return status.

      ! ---------------
      ! Local Variables
      ! ---------------

      INTEGER(KIND=MPI_ADDRESS_KIND) ::  &
         disp(3),          &
         view_disp,        &
         local_array_size
      INTEGER ::            &
         blocklen(3),       &
         type(3),           &
         gsizes(3),         &
         psizes(3),         &
         lsizes(3),         &
         global_offsets(3), &
         memsizes(3),       &
         offsets(3),        &
         I,                 &
         status


      ! ---------------------
      ! Executable Statements
      ! ---------------------

      IERR = 0

      CALL MPI_TYPE_GET_EXTENT( MPI_CHARACTER, CP_LB, SIZEOFCHAR,  IERR)
      CALL MPI_TYPE_GET_EXTENT( MPI_INTEGER,   CP_LB, SIZEOFINT,   IERR)
      CALL MPI_TYPE_GET_EXTENT( MPI_REAL8,     CP_LB, SIZEOFREAL8, IERR)

      !
      ! Setup an MPI derived type for writing cmdline options
      !
      blocklen(1) = 16
      blocklen(2) = 1024
      blocklen(3) = 1

      disp(1) = 0
      disp(2) = (SIZEOFINT*16)
      disp(3) = disp(2) + (SIZEOFCHAR*1024)

      type(1) = MPI_INTEGER
      type(2) = MPI_CHARACTER
      type(3) = MPI_UB

      CALL MPI_TYPE_CREATE_STRUCT( 3, blocklen, disp, type, CP_CMDLINE_TYPE, IERR)
      CALL MPI_TYPE_COMMIT( CP_CMDLINE_TYPE, IERR)

      CP_CMDLINE_EXTENT=SIZEOFINT*15
      CALL MPI_TYPE_GET_EXTENT ( CP_CMDLINE_TYPE, CP_LB, CP_CMDLINE_EXTENT, IERR )

      !
      ! Setup an MPI derived type for writing static config options
      !
      CALL MPI_TYPE_CONTIGUOUS( NUM_VARS, MPI_INTEGER, CP_GRIDSTOSUM_TYPE, IERR)
      CALL MPI_TYPE_COMMIT( CP_GRIDSTOSUM_TYPE, IERR)
      CALL MPI_TYPE_GET_EXTENT ( CP_GRIDSTOSUM_TYPE, CP_LB, CP_GRIDSTOSUM_EXTENT, IERR )

      CALL MPI_TYPE_CONTIGUOUS( 12, MPI_INTEGER, CP_PECOORDS_TYPE, IERR)
      CALL MPI_TYPE_COMMIT( CP_PECOORDS_TYPE, IERR)
      CALL MPI_TYPE_GET_EXTENT ( CP_PECOORDS_TYPE, CP_LB, CP_PECOORDS_EXTENT, IERR )

      !
      ! Setup an MPI derived type for writing an array of PE coordinates
      !
      gsizes(1) = 1
      gsizes(2) = NUMPES

      lsizes(1) = 1
      lsizes(2) = 1

      global_offsets(1) = 0
      global_offsets(2) = MYPE

      CALL MPI_TYPE_CREATE_SUBARRAY(2, gsizes, lsizes, global_offsets, MPI_ORDER_FORTRAN, &
                                    CP_PECOORDS_TYPE, CP_PECOORDSARRAY_TYPE, IERR)
      CALL MPI_TYPE_COMMIT(CP_PECOORDSARRAY_TYPE, IERR)
      CALL MPI_TYPE_EXTENT ( CP_PECOORDSARRAY_TYPE, CP_PECOORDSARRAY_EXTENT, IERR )

      !
      ! Setup an MPI derived type for writing a timestep header
      !
      CALL MPI_TYPE_CONTIGUOUS( MAX_GRIDS, MPI_REAL8, CP_GSUMOLD_TYPE, IERR)
      CALL MPI_TYPE_COMMIT( CP_GSUMOLD_TYPE, IERR)
      CALL MPI_TYPE_GET_EXTENT ( CP_GSUMOLD_TYPE, CP_LB, CP_GSUMOLD_EXTENT, IERR )

      blocklen(1) = 1
      blocklen(2) = 1
      blocklen(3) = 1

      disp(1) = 0
      disp(2) = SIZEOFINT
      disp(3) = disp(2)+(SIZEOFREAL8*MAX_GRIDS)

      type(1) = MPI_INTEGER
      type(2) = CP_GSUMOLD_TYPE
      type(3) = MPI_UB

      CALL MPI_TYPE_CREATE_STRUCT( 3, blocklen, disp, type, CP_TSHEADER_TYPE, IERR)
      CALL MPI_TYPE_COMMIT( CP_TSHEADER_TYPE, IERR)

      CP_TSHEADER_EXTENT=SIZEOFINT+(SIZEOFREAL8*MAX_GRIDS)
      CALL MPI_TYPE_GET_EXTENT ( CP_TSHEADER_TYPE, CP_LB, CP_TSHEADER_EXTENT, IERR )

      CALL MPI_TYPE_FREE ( CP_GSUMOLD_TYPE, IERR )

      !
      ! Setup an MPI derived type for writing a grid variable without the ghost cells
      !
      gsizes(1) = NX*NPX
      gsizes(2) = NY*NPY
      gsizes(3) = NZ*NPZ

      psizes(1) = NPX
      psizes(2) = NPY
      psizes(3) = NPZ

      lsizes(1) = NX
      lsizes(2) = NY
      lsizes(3) = NZ

      global_offsets(1) = MYPX*NX
      global_offsets(2) = MYPY*NY
      global_offsets(3) = MYPZ*NZ

      memsizes(1) = lsizes(1)+2
      memsizes(2) = lsizes(2)+2
      memsizes(3) = lsizes(3)+2

      offsets(1) = 1
      offsets(2) = 1
      offsets(3) = 1

      CALL MPI_TYPE_CREATE_SUBARRAY(3, memsizes, lsizes, offsets, MPI_ORDER_FORTRAN, MPI_REAL8, CP_NOGHOST_TYPE, IERR)
      CALL MPI_TYPE_COMMIT(CP_NOGHOST_TYPE, IERR)
      CALL MPI_TYPE_EXTENT ( CP_NOGHOST_TYPE, CP_NOGHOST_EXTENT, IERR )

      CALL MPI_TYPE_CREATE_SUBARRAY(3, gsizes, lsizes, global_offsets, MPI_ORDER_FORTRAN, MPI_REAL8, CP_TSGRID_TYPE, IERR)
      CALL MPI_TYPE_COMMIT(CP_TSGRID_TYPE, IERR)
      CALL MPI_TYPE_EXTENT ( CP_TSGRID_TYPE, CP_TSGRID_EXTENT, IERR )

      RETURN

   END SUBROUTINE MG_CHECKPOINT_CREATE_TYPES

!  ===================================================================================

   SUBROUTINE MG_CHECKPOINT_FREE_TYPES ( IERR )

      IMPLICIT NONE

      ! ---------------------
      ! Argument Declarations
      ! ---------------------

      INTEGER IERR                    ! Return status.

      ! ---------------
      ! Local Variables
      ! ---------------


      ! ---------------------
      ! Executable Statements
      ! ---------------------

      IERR = 0

      CALL MPI_TYPE_FREE ( CP_CMDLINE_TYPE, IERR)
      CALL MPI_TYPE_FREE ( CP_GRIDSTOSUM_TYPE, IERR)
      CALL MPI_TYPE_FREE ( CP_PECOORDS_TYPE, IERR)
      CALL MPI_TYPE_FREE ( CP_PECOORDSARRAY_TYPE, IERR)

      CALL MPI_TYPE_FREE ( CP_TSHEADER_TYPE, IERR )

      CALL MPI_TYPE_FREE ( CP_TSGRID_TYPE, IERR )
      CALL MPI_TYPE_FREE ( CP_NOGHOST_TYPE, IERR )

      RETURN

   END SUBROUTINE MG_CHECKPOINT_FREE_TYPES

!  ===================================================================================

!  ===================================================================================

   SUBROUTINE MG_CHECKPOINT_INIT ( IERR )

      IMPLICIT NONE

      ! ---------------------
      ! Argument Declarations
      ! ---------------------

      INTEGER IERR                    ! Return status.

      ! ---------------
      ! Local Variables
      ! ---------------

      INTEGER(KIND=MPI_ADDRESS_KIND) ::  &
         view_disp,        &
         local_array_size
      INTEGER ::            &
         I

      TYPE ( CMDLINE_OPTIONS_TYPE ) :: CMDLINE_OPTIONS
      INTEGER(KIND=MG_INT)          :: MY_GRIDS_TO_SUM(MAX_GRIDS) ! boolean
      TYPE ( PE_COORDS_TYPE )       :: PE_COORDS

      ! ---------------------
      ! Executable Statements
      ! ---------------------

      IERR = 0

      CMDLINE_OPTIONS%REPORT_DIFFUSION=REPORT_DIFFUSION
      CMDLINE_OPTIONS%SCALING=SCALING
      CMDLINE_OPTIONS%COMM_METHOD=COMM_METHOD
      CMDLINE_OPTIONS%STENCIL=STENCIL
      CMDLINE_OPTIONS%NUM_SPIKES=NUM_SPIKES
      CMDLINE_OPTIONS%NPX=NPX
      CMDLINE_OPTIONS%NPY=NPY
      CMDLINE_OPTIONS%NPZ=NPZ
      CMDLINE_OPTIONS%NX=NX
      CMDLINE_OPTIONS%NY=NY
      CMDLINE_OPTIONS%NZ=NZ
      CMDLINE_OPTIONS%NUM_VARS=NUM_VARS
      CMDLINE_OPTIONS%NUM_TSTEPS=NUM_TSTEPS
      CMDLINE_OPTIONS%PERCENT_SUM=PERCENT_SUM
      CMDLINE_OPTIONS%CP_INTERVAL=CP_INTERVAL
      CMDLINE_OPTIONS%CP_FILE=CP_FILE

      DO I = 1, NUM_VARS
         MY_GRIDS_TO_SUM(I)=GRIDS_TO_SUM(I)
      END DO

      PE_COORDS%MY_GLOBAL_NX_START=MY_GLOBAL_NX_START
      PE_COORDS%MY_GLOBAL_NY_START=MY_GLOBAL_NY_START
      PE_COORDS%MY_GLOBAL_NZ_START=MY_GLOBAL_NZ_START
      PE_COORDS%MY_GLOBAL_NX_END=MY_GLOBAL_NX_END
      PE_COORDS%MY_GLOBAL_NY_END=MY_GLOBAL_NY_END
      PE_COORDS%MY_GLOBAL_NZ_END=MY_GLOBAL_NZ_END
      PE_COORDS%NUM_NEIGHS=NUM_NEIGHS
      PE_COORDS%NUM_SUM_GRID=NUM_SUM_GRID
      PE_COORDS%MYPE=MYPE
      PE_COORDS%MYPX=MYPX
      PE_COORDS%MYPY=MYPY
      PE_COORDS%MYPZ=MYPZ


      view_disp=SIZEOFINT
!      WRITE (*,*) 'view_disp=', view_disp
      call MPI_FILE_SET_VIEW(cp_filehandle, view_disp, CP_CMDLINE_TYPE, &
                           CP_CMDLINE_TYPE, 'native', &
                           MPI_INFO_NULL, IERR)

      IF ( MYPE == ROOT_PE ) THEN
         CALL MPI_FILE_WRITE (cp_filehandle, CMDLINE_OPTIONS, CP_CMDLINE_EXTENT, MPI_BYTE, MPI_STATUS_IGNORE, IERR)
      END IF

      view_disp=SIZEOFINT+        &
                CP_CMDLINE_EXTENT
!      WRITE (*,*) 'view_disp=', view_disp
      call MPI_FILE_SET_VIEW(cp_filehandle, view_disp, CP_GRIDSTOSUM_TYPE, &
                           CP_GRIDSTOSUM_TYPE, 'native', &
                           MPI_INFO_NULL, IERR)

      IF ( MYPE == ROOT_PE ) THEN
         CALL MPI_FILE_WRITE (cp_filehandle, MY_GRIDS_TO_SUM, CP_GRIDSTOSUM_EXTENT, MPI_BYTE, MPI_STATUS_IGNORE, IERR)
      END IF

      view_disp=SIZEOFINT+           &
                CP_CMDLINE_EXTENT+   &
                CP_GRIDSTOSUM_EXTENT
!      WRITE (*,*) 'view_disp=', view_disp
      call MPI_FILE_SET_VIEW(cp_filehandle, view_disp, &
                           CP_PECOORDS_TYPE, &
                           CP_PECOORDSARRAY_TYPE, 'native', &
                           MPI_INFO_NULL, IERR)

      CALL MPI_FILE_WRITE_ALL ( cp_filehandle, PE_COORDS, 1, CP_PECOORDS_TYPE, MPI_STATUS_IGNORE, IERR)

      RETURN

   END SUBROUTINE MG_CHECKPOINT_INIT


!  ===================================================================================

!  ===================================================================================

   SUBROUTINE MG_CHECKPOINT_TSTEP_HEADER ( IERR )

      IMPLICIT NONE

      ! ---------------------
      ! Argument Declarations
      ! ---------------------

      INTEGER IERR                    ! Return status.

      ! ---------------
      ! Local Variables
      ! ---------------

      INTEGER(KIND=MPI_ADDRESS_KIND) ::  &
         cp_num,           & ! the current checkpoint number (0 based)
         cp_size,          & ! the size in bytes of one checkpoint (header+grids)
         view_disp
      INTEGER ::  &
         I,       &
         status

      TYPE ( TSHEADER_TYPE )        :: TSHEADER

      ! ---------------------
      ! Executable Statements
      ! ---------------------

      IERR = 0

!      IF ( MYPE == ROOT_PE ) THEN
!         WRITE(*,*)                                                                  &
!                    'CMDLINE_EXTENT=', CP_CMDLINE_EXTENT,                            &
!                    'GRIDSTOSUM_EXTENT=', CP_GRIDSTOSUM_EXTENT,                      &
!                    'GSUMOLD_EXTENT=', CP_GSUMOLD_EXTENT,                            &
!                    'TSHEADER_EXTENT=', CP_TSHEADER_EXTENT
!      END IF


      TSHEADER%TSTEP=TSTEP
      IF ( REPORT_DIFFUSION ) THEN
         DO I = 1, NUM_VARS
            TSHEADER%GSUM_OLD(I)=GSUM_OLD(I)
         END DO
      END IF

      cp_num =((tstep/cp_interval)-1)
      cp_size=(CP_TSHEADER_EXTENT+(CP_TSGRID_EXTENT*NUM_VARS))

      view_disp=SIZEOFINT+               &
                CP_CMDLINE_EXTENT+       &
                CP_GRIDSTOSUM_EXTENT+    &
                CP_PECOORDSARRAY_EXTENT+ &
                (cp_num*cp_size)
!      WRITE (*,*) 'view_disp=', view_disp
      call MPI_FILE_SET_VIEW(cp_filehandle, view_disp, CP_TSHEADER_TYPE, &
                           CP_TSHEADER_TYPE, 'native', &
                           MPI_INFO_NULL, IERR)

      IF ( MYPE == ROOT_PE ) THEN
         CALL MPI_FILE_WRITE (cp_filehandle, TSHEADER, CP_TSHEADER_EXTENT, MPI_BYTE, MPI_STATUS_IGNORE, IERR)
      END IF

      RETURN

   END SUBROUTINE MG_CHECKPOINT_TSTEP_HEADER

!  ===================================================================================

!  ===================================================================================

   SUBROUTINE MG_CHECKPOINT_TSTEP_VAR ( IVAR, IERR )

      IMPLICIT NONE

      ! ---------------------
      ! Argument Declarations
      ! ---------------------

      INTEGER IVAR
      INTEGER IERR                    ! Return status.

      ! ---------------
      ! Local Variables
      ! ---------------

      INTEGER(KIND=MPI_ADDRESS_KIND) ::  &
         cp_num,           &
         cp_size,          &
         view_disp,        &
         local_array_size
      INTEGER ::            &
         I,                 &
         status

      ! ---------------------
      ! Executable Statements
      ! ---------------------

      cp_num=((tstep/cp_interval)-1)
      cp_size=(CP_TSHEADER_EXTENT+(CP_TSGRID_EXTENT*NUM_VARS))

      view_disp=SIZEOFINT+                  &
                CP_CMDLINE_EXTENT+          &
                CP_GRIDSTOSUM_EXTENT+       &
                CP_PECOORDSARRAY_EXTENT+    &
                (cp_num*cp_size)+           &
                CP_TSHEADER_EXTENT+         &
                (CP_TSGRID_EXTENT*(IVAR-1))

!      WRITE (*,*) 'view_disp=', view_disp
      call MPI_FILE_SET_VIEW(cp_filehandle, view_disp, &
                           MPI_REAL8, &
                           CP_TSGRID_TYPE, 'native', &
                           MPI_INFO_NULL, IERR)

      CALL MG_CHECKPOINT_TSTEP_GRID(GRID(:, :, :, IVAR), IERR)

      RETURN

   END SUBROUTINE MG_CHECKPOINT_TSTEP_VAR

!  ===================================================================================

!  ===================================================================================

   SUBROUTINE MG_CHECKPOINT_TSTEP_GRID ( GRID, IERR )

      ! ---------------------
      ! Argument Declarations
      ! ---------------------

      REAL(KIND=MG_REAL), DIMENSION(0:NX+1, 0:NY+1, 0:NZ+1), INTENT(INOUT) :: GRID
      INTEGER IERR           ! Return status

      ! ---------------------
      ! Executable Statements
      ! ---------------------

      CALL MPI_FILE_WRITE_ALL ( cp_filehandle, GRID, 1, CP_NOGHOST_TYPE, MPI_STATUS_IGNORE, IERR)

      RETURN

   END SUBROUTINE MG_CHECKPOINT_TSTEP_GRID

!  ===================================================================================


!  ===================================================================================

   SUBROUTINE MG_CHECKPOINT_VERIFY ( IERR )

      IMPLICIT NONE

      ! ---------------------
      ! Argument Declarations
      ! ---------------------

      INTEGER IERR                    ! Return status.

      ! ---------------
      ! Local Variables
      ! ---------------

      INTEGER(KIND=MPI_ADDRESS_KIND) ::  &
         cp_num,           &
         cp_size,          &
         view_disp,        &
         local_array_size
      INTEGER ::  &
         I,       &
         J,       &
         K

      TYPE ( CMDLINE_OPTIONS_TYPE ) :: CMDLINE_OPTIONS
      INTEGER(KIND=MG_INT)          :: MY_GRIDS_TO_SUM(MAX_GRIDS) ! boolean
      TYPE ( PE_COORDS_TYPE )       :: PE_COORDS
      TYPE ( TSHEADER_TYPE )        :: TSHEADER

      REAL(KIND=MG_REAL), DIMENSION(0:NX+1, 0:NY+1, 0:NZ+1) :: MYGRID

      ! ---------------------
      ! Executable Statements
      ! ---------------------

      IERR = 0


      cp_num=((tstep/cp_interval)-1)
      cp_size=(CP_TSHEADER_EXTENT+(CP_TSGRID_EXTENT*NUM_VARS))


      view_disp=SIZEOFINT
!      WRITE (*,*) 'view_disp=', view_disp
      call MPI_FILE_SET_VIEW(cp_filehandle, view_disp, CP_CMDLINE_TYPE, &
                           CP_CMDLINE_TYPE, 'native', &
                           MPI_INFO_NULL, IERR)

      IF ( MYPE == ROOT_PE ) THEN
         CALL MPI_FILE_READ (cp_filehandle, CMDLINE_OPTIONS, CP_CMDLINE_EXTENT, MPI_BYTE, MPI_STATUS_IGNORE, IERR)

   !      IF (CMDLINE_OPTIONS%REPORT_DIFFUSION /= REPORT_DIFFUSION ) THEN
   !         WRITE (*,*) 'CMDLINE_OPTIONS%REPORT_DIFFUSION != REPORT_DIFFUSION'
   !      END IF
         IF (CMDLINE_OPTIONS%SCALING         /= SCALING )         THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%SCALING         != SCALING', CMDLINE_OPTIONS%SCALING, '!=', SCALING
         END IF
         IF (CMDLINE_OPTIONS%COMM_METHOD     /= COMM_METHOD )     THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%COMM_METHOD     != COMM_METHOD', CMDLINE_OPTIONS%COMM_METHOD, '!=', COMM_METHOD
         END IF
         IF (CMDLINE_OPTIONS%STENCIL         /= STENCIL )         THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%STENCIL         != STENCIL'
         END IF
         IF (CMDLINE_OPTIONS%NUM_SPIKES      /= NUM_SPIKES )      THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%NUM_SPIKES      != NUM_SPIKES', CMDLINE_OPTIONS%NUM_SPIKES, '!=', NUM_SPIKES
         END IF
         IF (CMDLINE_OPTIONS%NPX             /= NPX )             THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%NPX             != NPX'
         END IF
         IF (CMDLINE_OPTIONS%NPY             /= NPY )             THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%NPY             != NPY'
         END IF
         IF (CMDLINE_OPTIONS%NPZ             /= NPZ )             THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%NPZ             != NPZ'
         END IF
         IF (CMDLINE_OPTIONS%NX              /= NX )              THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%NX              != NX'
         END IF
         IF (CMDLINE_OPTIONS%NY              /= NY )              THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%NY              != NY'
         END IF
         IF (CMDLINE_OPTIONS%NZ              /= NZ )              THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%NZ              != NZ'
         END IF
         IF (CMDLINE_OPTIONS%NUM_VARS        /= NUM_VARS )        THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%NUM_VARS        != NUM_VARS'
         END IF
         IF (CMDLINE_OPTIONS%NUM_TSTEPS      /= NUM_TSTEPS )      THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%NUM_TSTEPS      != NUM_TSTEPS'
         END IF
         IF (CMDLINE_OPTIONS%PERCENT_SUM     /= PERCENT_SUM )     THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%PERCENT_SUM     != PERCENT_SUM'
         END IF
         IF (CMDLINE_OPTIONS%CP_INTERVAL     /= CP_INTERVAL )     THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%CP_INTERVAL     != CP_INTERVAL'
         END IF
         IF (CMDLINE_OPTIONS%CP_FILE         /= CP_FILE )     THEN
            WRITE (*,*) 'CMDLINE_OPTIONS%CP_FILE         != CP_FILE'
         END IF
      END IF


      view_disp=SIZEOFINT+        &
                CP_CMDLINE_EXTENT
!      WRITE (*,*) 'view_disp=', view_disp
      call MPI_FILE_SET_VIEW(cp_filehandle, view_disp, CP_GRIDSTOSUM_TYPE, &
                           CP_GRIDSTOSUM_TYPE, 'native', &
                           MPI_INFO_NULL, IERR)

      IF ( MYPE == ROOT_PE ) THEN
         CALL MPI_FILE_READ (cp_filehandle, MY_GRIDS_TO_SUM, CP_GRIDSTOSUM_EXTENT, MPI_BYTE, MPI_STATUS_IGNORE, IERR)
      END IF

      view_disp=SIZEOFINT+                  &
                CP_CMDLINE_EXTENT+          &
                CP_GRIDSTOSUM_EXTENT
!      WRITE (*,*) 'view_disp=', view_disp
      call MPI_FILE_SET_VIEW(cp_filehandle, view_disp, &
                           CP_PECOORDS_TYPE, &
                           CP_PECOORDSARRAY_TYPE, 'native', &
                           MPI_INFO_NULL, IERR)

      CALL MPI_FILE_READ_ALL ( cp_filehandle, PE_COORDS, 1, CP_PECOORDS_TYPE, MPI_STATUS_IGNORE, IERR)

      IF ( PE_COORDS%MY_GLOBAL_NX_START /= MY_GLOBAL_NX_START ) THEN
         WRITE (*,*) ' PE_COORDS%MY_GLOBAL_NX_START /= MY_GLOBAL_NX_START'
      END IF
      IF ( PE_COORDS%MY_GLOBAL_NY_START /= MY_GLOBAL_NY_START ) THEN
         WRITE (*,*) ' PE_COORDS%MY_GLOBAL_NY_START /= MY_GLOBAL_NY_START'
      END IF
      IF ( PE_COORDS%MY_GLOBAL_NZ_START /= MY_GLOBAL_NZ_START ) THEN
         WRITE (*,*) ' PE_COORDS%MY_GLOBAL_NZ_START /= MY_GLOBAL_NZ_START'
      END IF
      IF ( PE_COORDS%MY_GLOBAL_NX_END /= MY_GLOBAL_NX_END ) THEN
         WRITE (*,*) ' PE_COORDS%MY_GLOBAL_NX_END /= MY_GLOBAL_NX_END'
      END IF
      IF ( PE_COORDS%MY_GLOBAL_NY_END /= MY_GLOBAL_NY_END ) THEN
         WRITE (*,*) ' PE_COORDS%MY_GLOBAL_NY_END /= MY_GLOBAL_NY_END'
      END IF
      IF ( PE_COORDS%MY_GLOBAL_NZ_END /= MY_GLOBAL_NZ_END ) THEN
         WRITE (*,*) ' PE_COORDS%MY_GLOBAL_NZ_END /= MY_GLOBAL_NZ_END'
      END IF
      IF ( PE_COORDS%NUM_NEIGHS /= NUM_NEIGHS ) THEN
         WRITE (*,*) ' PE_COORDS%NUM_NEIGHS /= NUM_NEIGHS'
      END IF
      IF ( PE_COORDS%NUM_SUM_GRID /= NUM_SUM_GRID ) THEN
         WRITE (*,*) ' PE_COORDS%NUM_SUM_GRID /= NUM_SUM_GRID'
      END IF
      IF ( PE_COORDS%MYPE /= MYPE ) THEN
         WRITE (*,*) ' PE_COORDS%MYPE /= MYPE'
      END IF
      IF ( PE_COORDS%MYPX /= MYPX ) THEN
         WRITE (*,*) ' PE_COORDS%MYPX /= MYPX'
      END IF
      IF ( PE_COORDS%MYPY /= MYPY ) THEN
         WRITE (*,*) ' PE_COORDS%MYPY /= MYPY'
      END IF
      IF ( PE_COORDS%MYPZ /= MYPZ ) THEN
         WRITE (*,*) ' PE_COORDS%MYPZ /= MYPZ'
      END IF


      view_disp=SIZEOFINT+                  &
                CP_CMDLINE_EXTENT+          &
                CP_GRIDSTOSUM_EXTENT+       &
                CP_PECOORDSARRAY_EXTENT+    &
                (cp_num*cp_size)
!      WRITE (*,*) 'view_disp=', view_disp
      call MPI_FILE_SET_VIEW(cp_filehandle, view_disp, CP_TSHEADER_TYPE, &
                           CP_TSHEADER_TYPE, 'native', &
                           MPI_INFO_NULL, IERR)

      IF ( MYPE == ROOT_PE ) THEN
         CALL MPI_FILE_READ (cp_filehandle, TSHEADER, CP_TSHEADER_EXTENT, MPI_BYTE, MPI_STATUS_IGNORE, IERR)
         IF ( TSHEADER%TSTEP /= TSTEP ) THEN
            WRITE (*,*) ' TSHEADER%TSTEP != TSTEP', TSHEADER%TSTEP, '!=', TSTEP
         END IF
         IF ( REPORT_DIFFUSION ) THEN
            DO I = 1, NUM_VARS
               IF ( TSHEADER%GSUM_OLD(I) /= GSUM_OLD(I) ) THEN
                     WRITE (*,*) ' TSHEADER%GSUM_OLD(I) != GSUM_OLD(I)'
               END IF
            END DO
         END IF
      END IF



      cp_num =((tstep/cp_interval)-1)
      cp_size=(CP_TSHEADER_EXTENT+(CP_TSGRID_EXTENT*NUM_VARS))

      view_disp=SIZEOFINT+                  &
                CP_CMDLINE_EXTENT+          &
                CP_GRIDSTOSUM_EXTENT+       &
                CP_PECOORDSARRAY_EXTENT+    &
                (cp_num*cp_size)+           &
                CP_TSHEADER_EXTENT
!      WRITE (*,*) 'view_disp=', view_disp
      call MPI_FILE_SET_VIEW(cp_filehandle, view_disp, &
                           MPI_REAL8, &
                           CP_TSGRID_TYPE, 'native', &
                           MPI_INFO_NULL, IERR)

      CALL MPI_FILE_READ_ALL ( cp_filehandle, MYGRID, 1, CP_NOGHOST_TYPE, MPI_STATUS_IGNORE, IERR)

      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               IF ( MYGRID ( I, J, K ) /= GRID ( I, J, K, 1) ) THEN
                  WRITE (*,*) 'MYGRID ( I, J, K ) != GRID ( I, J, K, 1)',&
                       MYGRID (I,J,K), '!=', GRID(I,J,K,1), I,J,K
               END IF
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE MG_CHECKPOINT_VERIFY

#endif _MG_CHECKPT

#endif _MG_MPI

!  ===================================================================================

END MODULE MG_CHECKPOINT_MOD
