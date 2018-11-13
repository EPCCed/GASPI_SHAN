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

MODULE MG_UTILS_MOD

  USE F_SHAN_SEGMENT_MOD
  USE MG_CONSTANTS_MOD
  USE MG_OPTIONS_MOD
  use, intrinsic :: ISO_C_BINDING

  IMPLICIT NONE

CONTAINS

  ! MG_INIT             : Set up MPI environment, problem input.
  ! MG_PRINT_HEADER     : Write configuration to STDIO.
  ! MG_GRID_INIT        : Allocate GRID arrays, call MG_INIT_GRID.
  ! MG_INSERT_SPIKE     : Insert heat sources.
  ! MG_GRID_DEALLOC     : Deallocate GRID arrays.
  ! MG_ASSERT           : Error checking.
  ! MG_INIT_GRID        : Initialize GRID arrays.
  ! MG_TIMER            : Returns time.
  ! MG_COMPUTE_STDDEV   : Compute standard deviation of input profiling data.

  !  =================================================================================

  SUBROUTINE MG_INIT ( IERR )

    IMPLICIT NONE

    INTEGER ::  &
         IERR                    ! Return status.

    ! ---------------
    ! Local Variables
    ! ---------------

    INTEGER ::  &
         I,                    &    ! Counter.
         IX,                   &
         IY,                   &
         IZ,                   &
         J,                    &    ! Counter.
         MYPE_XY,              &    ! tmp var
         NVARS,                &    ! tmp var
         OFFSET_X,             &
         OFFSET_Y,             &
         OFFSET_Z,             &
         REMAINDER,            &
         XLOC,                 &
         YLOC,                 &
         ZLOC

    ! ------------
    ! Local Arrays
    ! ------------

    REAL(KIND=MG_REAL), DIMENSION(:,:), ALLOCATABLE ::   &
         RSPIKE_LOC            ! Temporary for random number generator.

    ! ---------------------
    ! Executable Statements 
    ! ---------------------

    IERR = 0

#if defined _MG_MPI
    ! Parallel processing configuration:

    CALL MPI_COMM_DUP ( MPI_COMM_WORLD, MPI_COMM_MG, IERR )

    CALL MPI_ERRHANDLER_SET ( MPI_COMM_MG, MPI_ERRORS_ARE_FATAL, IERR )

    CALL MPI_COMM_RANK ( MPI_COMM_MG, MYPE, IERR )

    CALL MPI_COMM_SIZE ( MPI_COMM_MG, NUMPES, IERR )

#elif defined _MG_SERIAL
    MYPE = 0    ! All formulas based on these values should work correctly.
    NUMPES = 1
#endif

    ! ---------------------------------
    ! Set position in 3d processor grid
    ! ---------------------------------

    MYPE_XY = MOD ( MYPE, NPX*NPY )
    MYPY = MYPE_XY / NPX
    REMAINDER = MOD ( MYPE_XY, NPX )
    IF ( REMAINDER /= 0 ) THEN
       MYPX = REMAINDER
    ELSE
       MYPX = 0
    END IF
    MYPZ = MYPE / ( NPX*NPY )


    ! --------------
    ! Set neighbors.
    ! --------------
    ALLOCATE ( NEIGHBORS(MAX_NUM_NEIGHBORS), STAT = IERR )      
    NEIGHBORS(1:MAX_NUM_NEIGHBORS) = -1
    CALL MG_ASSERT ( IERR, 'MG_INIT: ALLOCATE ( NEIGHBORS )', MAX_NUM_NEIGHBORS )

    ALLOCATE ( NEIGHBORS_ORIG(MAX_NUM_NEIGHBORS), STAT = IERR )
    NEIGHBORS_ORIG(1:MAX_NUM_NEIGHBORS) = -1
    CALL MG_ASSERT ( IERR, 'MG_INIT: ALLOCATE ( NEIGHBORS_ORIG )', MAX_NUM_NEIGHBORS )

    NUM_NEIGHS = 0
    IF ( MYPY /= 0  ) THEN
       NEIGHBORS(SOUTH) = MYPE - NPX
       NUM_NEIGHS = NUM_NEIGHS + 1
    END IF
    IF ( MYPY /= NPY-1 ) THEN
       NEIGHBORS(NORTH) = MYPE + NPX
       NUM_NEIGHS = NUM_NEIGHS + 1
    END IF
    IF ( MYPX /= 0 ) THEN
       NEIGHBORS(WEST) = MYPE - 1
       NUM_NEIGHS = NUM_NEIGHS + 1
    END IF
    IF ( MYPX /= NPX-1 ) THEN
       NEIGHBORS(EAST) = MYPE + 1
       NUM_NEIGHS = NUM_NEIGHS + 1
    END IF
    IF ( MYPZ /= 0 ) THEN
       NEIGHBORS(BACK) = MYPE - ( NPX*NPY )
       NUM_NEIGHS = NUM_NEIGHS + 1
    END IF
    IF ( MYPZ /= NPZ-1 ) THEN
       NEIGHBORS(FRONT) = MYPE + ( NPX*NPY )
       NUM_NEIGHS = NUM_NEIGHS + 1
    END IF

    NEIGHBORS_ORIG = NEIGHBORS

    ALLOCATE ( SPIKES(NUM_VARS,NUM_SPIKES), STAT = IERR )
    CALL MG_ASSERT ( IERR, 'MG_INIT: ALLOCATE ( SPIKES )', NUM_VARS*NUM_SPIKES )
    SPIKES = -1.0

    ALLOCATE ( SPIKE_LOC(0:3,NUM_SPIKES), STAT = IERR )
    CALL MG_ASSERT ( IERR, 'MG_INIT: ALLOCATE ( SPIKE_LOC )', 4*NUM_SPIKES )
    SPIKE_LOC = -1

    ALLOCATE ( SOURCE_TOTAL(NUM_VARS), STAT = IERR )
    CALL MG_ASSERT ( IERR, 'MG_INIT: ALLOCATE ( SOURCE_TOTAL )', NUM_VARS )
    SOURCE_TOTAL = 0.0

    ! -------------------------------------------------------
    ! ROOT_PE generates and distributes SPIKES and SPIKE_LOC.
    ! Each GRIDi has a unique SPIKE value, but for each set,
    ! all in same location within the GRIDi array.
    ! Global location computed, PEs determines ownership.
    ! -------------------------------------------------------

    ! Determine global indices (excluding ghost)

    GLOBAL_NX = NX * NPX
    GLOBAL_NY = NY * NPY
    GLOBAL_NZ = NZ * NPZ

    MY_GLOBAL_NX_START = NX * MYPX + 1
    MY_GLOBAL_NY_START = NY * MYPY + 1
    MY_GLOBAL_NZ_START = NZ * MYPZ + 1

    MY_GLOBAL_NX_END   = MY_GLOBAL_NX_START + NX - 1
    MY_GLOBAL_NY_END   = MY_GLOBAL_NY_START + NY - 1
    MY_GLOBAL_NZ_END   = MY_GLOBAL_NZ_START + NZ - 1

    IF ( MYPE == ROOT_PE ) THEN

       CALL RANDOM_NUMBER ( SPIKES ( 1:NUM_VARS,1:NUM_SPIKES ) )
       SPIKES = SPIKES * GLOBAL_NX * GLOBAL_NY * GLOBAL_NZ

       ALLOCATE ( RSPIKE_LOC(3,NUM_SPIKES), STAT = IERR )
       CALL RANDOM_NUMBER ( RSPIKE_LOC )

       ! First spike set to problem dimension, inserted center of global grid.
       SPIKES(1,1) = REAL ( GLOBAL_NX * GLOBAL_NY * GLOBAL_NZ )

       SPIKE_LOC(0,1) = -1 
       SPIKE_LOC(1,1) = GLOBAL_NX / 2
       SPIKE_LOC(2,1) = GLOBAL_NY / 2
       SPIKE_LOC(3,1) = GLOBAL_NZ / 2

       ! Set additional spikes randomly about global grid

       DO I = 2, NUM_SPIKES
          SPIKE_LOC(0,I) = -1   ! Owning MPI rank will set to its MYPE.
          SPIKE_LOC(1,I) = INT ( RSPIKE_LOC(1,I) * GLOBAL_NX )
          SPIKE_LOC(2,I) = INT ( RSPIKE_LOC(2,I) * GLOBAL_NY )
          SPIKE_LOC(3,I) = INT ( RSPIKE_LOC(3,I) * GLOBAL_NZ )
       END DO

       DEALLOCATE ( RSPIKE_LOC )

    END IF ! ROOT_PE work.

#if defined _MG_MPI

    ! -----------------
    ! Distribute SPIKES
    ! -----------------

    CALL MPI_BCAST ( SPIKES, NUM_VARS*NUM_SPIKES, MG_MPI_REAL, ROOT_PE,  &
         MPI_COMM_MG, IERR )
    CALL MG_ASSERT ( IERR, 'INIT: MPI_BCAST(SPIKES)', IERR )

    CALL MPI_BCAST ( SPIKE_LOC, 4*NUM_SPIKES, MPI_INTEGER, ROOT_PE,  &
         MPI_COMM_MG, IERR )
    CALL MG_ASSERT ( IERR, 'INIT: MPI_BCAST(SPIKE_LOC)', IERR )

#endif

    ! Owning MPI rank determines converts the global location to its local array indices.

    DO I = 1, NUM_SPIKES

       XLOC = SPIKE_LOC ( 1, I ) ! Global values here.
       YLOC = SPIKE_LOC ( 2, I )
       ZLOC = SPIKE_LOC ( 3, I )

       IF ( ( MY_GLOBAL_NX_START <= XLOC ) .AND. ( XLOC <= MY_GLOBAL_NX_END ) .AND.  &
            ( MY_GLOBAL_NY_START <= YLOC ) .AND. ( YLOC <= MY_GLOBAL_NY_END ) .AND.  &
            ( MY_GLOBAL_NZ_START <= ZLOC ) .AND. ( ZLOC <= MY_GLOBAL_NZ_END ) ) THEN

          SPIKE_LOC ( 0, I ) = MYPE

          SPIKE_LOC ( 1, I ) = SPIKE_LOC ( 1, I ) - MY_GLOBAL_NX_START
          SPIKE_LOC ( 2, I ) = SPIKE_LOC ( 2, I ) - MY_GLOBAL_NY_START 
          SPIKE_LOC ( 3, I ) = SPIKE_LOC ( 3, I ) - MY_GLOBAL_NZ_START 

       ELSE

          ! Not owner, so set to -1 (to make obvious).

          SPIKE_LOC ( 0, I ) = -1

          SPIKE_LOC ( 1, I ) = -1
          SPIKE_LOC ( 2, I ) = -1
          SPIKE_LOC ( 3, I ) = -1

       END IF

    END DO

    RETURN

  END SUBROUTINE MG_INIT

  !  ===================================================================================

  SUBROUTINE MG_PRINT_HEADER ( COMM_METHOD, STENCIL, IERR )

    IMPLICIT NONE

    ! Argument Declarations
    INTEGER, INTENT(IN)  ::       &
         COMM_METHOD, STENCIL

    INTEGER, INTENT(OUT) :: IERR

    !  Purpose
    !  =======
    !  Collate, process, and report performance results.

    ! Local Scalars
    CHARACTER(LEN=30) ::        &
         TEST_DATE
    CHARACTER(LEN=30) ::        &
         TEST_TIME

    INTEGER ::            &
         I,                              &
         ICLOCK_RATE,                    &
         IDUM,                           &
         LEN

    REAL(KIND=MG_REAL4) ::       &
         CLOCK_RES = 0.0

    INTEGER(KIND=MG_INT8), PARAMETER ::  &
         SIZE_OF_DATA = 8


    ! ---------------------
    ! Executable Statements
    ! ---------------------

    IERR = 0

    IF ( MYPE /= ROOT_PE ) &
         RETURN

    CALL DATE_AND_TIME ( TEST_DATE, TEST_TIME )

    WRITE(*,*)
    WRITE(*,100)
    WRITE(*,*) '          Mantevo miniapp MiniGhost experiment'
    WRITE(*,100)

#if defined _MG_MPI
    CLOCK_RES = REAL(MPI_WTICK ( ))

    WRITE(*,*)
    WRITE(*,*) 'Communication strategy: one variable at a time (COMM_METHOD_SVAF)'

#endif

    WRITE(*,*)

    SELECT CASE ( STENCIL )

    CASE ( STENCIL_NONE )

       WRITE(*,*) 'No computation inserted.'

    CASE ( STENCIL_2D5PT )

       WRITE(*,*) 'Computation: 5 pt difference stencil on a 2D grid (STENCIL_2D5PT)'

    CASE ( STENCIL_3D7PT )

       WRITE(*,*) 'Computation: 7 pt difference stencil on a 3D grid (STENCIL_3D7PT)'

    CASE DEFAULT

       WRITE(*,*) '** Warning ** Unknown computation'

    END SELECT

    WRITE(*,*)
    WRITE(*,101) NX * NPX, NY * NPY, NZ * NPZ
    WRITE(*,102) NX, NY, NZ
    WRITE(*,*)
    WRITE(*,103) NUM_VARS
    WRITE(*,*)
    IF ( REPORT_DIFFUSION /= 0 ) THEN
       WRITE(*,104) REPORT_DIFFUSION, ERROR_TOL
    END IF
    WRITE(*,105) NUM_SUM_GRID, PERCENT_SUM
    WRITE(*,*)
    WRITE(*,110) NUM_TSTEPS
    WRITE(*,*)
#if defined _MG_MPI
    WRITE(*,120) NPX, NPY, NPZ
    WRITE(*,*)
    IF ( SCALING == SCALING_STRONG ) THEN   ! Not that it really matters.
       WRITE(*,*) 'MPI version, strong scaling'
    ELSE
       WRITE(*,*) 'MPI version, weak scaling'
    END IF
    WRITE(*,*)
    IF ( NUMPES == 1 ) THEN
       WRITE(*,121) TEST_TIME, TEST_DATE
    ELSE
       WRITE(*,122) NUMPES, TEST_TIME, TEST_DATE
    END IF
#elif defined _MG_SERIAL
    IF ( SCALING == SCALING_STRONG ) THEN   ! Not that it really matters.
       WRITE(*,123) 'Serial version, strong scaling', TEST_TIME, TEST_DATE
    ELSE
       WRITE(*,123) 'Serial version, weak scaling', TEST_TIME, TEST_DATE
    END IF
    WRITE(*,*)
#endif

    ! Format statements

100 FORMAT ( ' =================================================' )

101 FORMAT ( '      Global Grid Dimension: ', I8, ', ', I8, ', ', I8 )
102 FORMAT ( '      Local Grid Dimension : ', I8, ', ', I8, ', ', I8 )

103 FORMAT ( ' Number of variables: ', I2 )
104 FORMAT ( ' Error reported every ', I3, ' time steps. Tolerance is ', 1PE8.2 )
105 FORMAT ( ' Number of variables reduced each time step: ', I2, '; requested  ', I3, '%.')

110 FORMAT ( '      Time steps: ', I6 )

120 FORMAT ( '      Task grid: ', I5, ',', I5, ',', I5 )

121 FORMAT ( ' 1 process executing', // &
         ' Program execution at ', A10, ' on ', A8, '.' )

122 FORMAT ( I4, ' processes executing', // &
         ' Program execution at ', A10, ' on ', A8, '.' )
123 FORMAT ( A32, ', ', A10, ' on ', A8, '.' )

  END SUBROUTINE MG_PRINT_HEADER

  !  ===================================================================================

  SUBROUTINE MG_SET_COMM_TYPE ( IERR )

    IMPLICIT NONE

    INTEGER ::  &
         IERR                    ! Return status.

    INTEGER :: I, J, K, L, M
    TYPE(C_PTR) :: CPTR, RECV_IDX_CPTR, SEND_IDX_CPTR, &
         NELEM_SEND_CPTR, NELEM_RECV_CPTR, &
         SEND_SZ_CPTR, RECV_SZ_CPTR

    REAL(KIND=MG_REAL)    :: SZ_REAL
    INTEGER(KIND=MG_INT4) :: SZ_INT4
    INTEGER(KIND=MG_INT8) :: SZ_INT8
    INTEGER(KIND=MG_INT4) :: NSEND, NRECV, TYPE_ID

    INTEGER(KIND=MG_INT),   DIMENSION(:,:), POINTER :: &
         DUMMY_SEND_SZ, DUMMY_RECV_SZ

    
    ! ---------------
    ! Local Variables
    ! ---------------


    IERR = 0

    DO M = 1, NUM_TYPE

       TYPE_ID = M - 1       
       CALL F_SHAN_TYPE_OFFSET(NEIGHBOR_HOOD_ID, TYPE_ID, &
            NELEM_SEND_CPTR, NELEM_RECV_CPTR, &
            SEND_SZ_CPTR, RECV_SZ_CPTR, &
            SEND_IDX_CPTR, RECV_IDX_CPTR)
       
       CALL C_F_POINTER(NELEM_SEND_CPTR, NELEM_SEND, (/NUM_NEIGHS/))
       CALL C_F_POINTER(NELEM_RECV_CPTR, NELEM_RECV, (/NUM_NEIGHS/))
       CALL C_F_POINTER(SEND_SZ_CPTR, SEND_SZ, (/NUM_NEIGHS/))
       CALL C_F_POINTER(RECV_SZ_CPTR, RECV_SZ, (/NUM_NEIGHS/))
       CALL C_F_POINTER(SEND_IDX_CPTR, OFFSET_SEND, (/TYPE_MAX_NSEND * NUM_NEIGHS/))
       CALL C_F_POINTER(RECV_IDX_CPTR, OFFSET_RECV, (/TYPE_MAX_NRECV * NUM_NEIGHS/))
       
       DO J = 1, NUM_NEIGHS
          I = NEIGHBORS_IDX(J)
          IF (I == BACK .OR. I == FRONT) THEN
             NELEM_SEND(J) = 1
             NELEM_RECV(J) = 1
          END IF
          IF (I == EAST .OR. I == WEST) THEN
             NELEM_SEND(J) = (NZ+2)*(NY+2)
             NELEM_RECV(J) = (NZ+2)*(NY+2)
          END IF
          IF (I == NORTH .OR. I == SOUTH) THEN
             NELEM_SEND(J) = NZ+2
             NELEM_RECV(J) = NZ+2
          ENDIF
       END DO
       
       DO J = 1, NUM_NEIGHS
          I = NEIGHBORS_IDX(J)
          IF (I == BACK .OR. I == FRONT) THEN
             SEND_SZ(J) = (NX+2)*(NY+2)*C_SIZEOF(SZ_REAL)
             RECV_SZ(J) = (NX+2)*(NY+2)*C_SIZEOF(SZ_REAL)
          END IF
          IF (I == WEST .OR. I == EAST) THEN
             SEND_SZ(J) = C_SIZEOF(SZ_REAL)
             RECV_SZ(J) = C_SIZEOF(SZ_REAL)
          END IF
          IF (I == SOUTH .OR. I == NORTH) THEN
             SEND_SZ(J) = (NX+2)*C_SIZEOF(SZ_REAL)
             RECV_SZ(J) = (NX+2)*C_SIZEOF(SZ_REAL)
          END IF

       ENDDO

       DO J = 1, NUM_NEIGHS
          NSEND = (J - 1) * F_SHAN_NELEM_SEND(M) + 1
          NRECV = (J - 1) * F_SHAN_NELEM_RECV(M) + 1
          I = NEIGHBORS_IDX(J)       
          IF (I == BACK) THEN
             OFFSET_SEND(NSEND)  = GET_OFFSET(0,0,1)
             OFFSET_RECV(NRECV)  = GET_OFFSET(0,0,0)
             NSEND = NSEND + 1
             NRECV = NRECV + 1
          ENDIF
          IF (I == FRONT) THEN
             OFFSET_SEND(NSEND)  = GET_OFFSET(0,0,NZ)
             OFFSET_RECV(NRECV)  = GET_OFFSET(0,0,NZ+1)
             NSEND = NSEND + 1
             NRECV = NRECV + 1
          ENDIF
          IF (I == EAST) THEN
             DO K = 0, NZ+1
                DO L = 0, NY+1
                   OFFSET_SEND(NSEND + L)  = GET_OFFSET(NX,L,K)
                   OFFSET_RECV(NRECV + L)  = GET_OFFSET(NX+1,L,K)
                END DO
                NSEND = NSEND + NY+2
                NRECV = NRECV + NY+2                
             END DO
          END IF
          IF (I == WEST) THEN
             DO K = 0, NZ+1
                DO L = 0, NY+1
                   OFFSET_SEND(NSEND + L)  = GET_OFFSET(1,L,K)
                   OFFSET_RECV(NRECV + L)  = GET_OFFSET(0,L,K)
                END DO
                NSEND = NSEND + NY+2
                NRECV = NRECV + NY+2                
             END DO
          END IF
          IF (I == NORTH) THEN
             DO K = 0, NZ+1
                OFFSET_SEND(NSEND + K) = GET_OFFSET(0,NY,K)
                OFFSET_RECV(NRECV + K) = GET_OFFSET(0,NY+1,K)
             END DO
             NSEND = NSEND + NZ+2
             NRECV = NRECV + NZ+2
          END IF
          IF (I == SOUTH) THEN
             DO K = 0, NZ+1
                OFFSET_SEND(NSEND + K) = GET_OFFSET(0,1,K)
                OFFSET_RECV(NRECV + K) = GET_OFFSET(0,0,K)
             END DO
             NSEND = NSEND + NZ+2
             NRECV = NRECV + NZ+2
          END IF
       END DO


    END DO

    
  END SUBROUTINE MG_SET_COMM_TYPE


    SUBROUTINE MG_GRID_INIT ( IERR ) 

    INTEGER, INTENT(OUT) ::  &
         IERR                    ! Return status.

    ! ------------------
    ! Local Declarations
    ! ------------------

    LOGICAL ::               &
         MY_SPIKE

    INTEGER ::  &
         SPIKE_LOC_X,          & !
         SPIKE_LOC_Y,          & !
         SPIKE_LOC_Z

    INTEGER :: I, J
    TYPE(C_PTR) :: CPTR
       INTEGER(KIND=MG_INT8) :: SIZE
    INTEGER, ALLOCATABLE :: SHAPE(:)
    
    ! ---------------
    ! Local Variables
    ! ---------------

    REAL(KIND=MG_REAL)    :: SZ_REAL
    INTEGER(KIND=MG_INT4) :: MAX_NSEND, MAX_NRECV, SEND_SZ, RECV_SZ
    INTEGER(KIND=MG_INT4) :: MAX_NELEM_SEND, MAX_NELEM_RECV


    IERR = 0
    
    ! Determine global indices (excluding ghost)
    IF ( NUM_VARS > 40 ) THEN
       IERR = -1
       CALL MG_ASSERT ( IERR, 'GRID_INIT: TOO MANY VARS', NUM_VARS )
    ENDIF

    IF ( NUM_VARS > 0 ) THEN

       ! ---------------------------------------------
       ! Allocate node-local shared segment
       ! ---------------------------------------------
       SIZE = NUM_VARS*(NX+2)*(NY+2)*(NZ+2)*C_SIZEOF(SZ_REAL)
       CALL F_SHAN_ALLOC_SHARED(SEGMENT_ID, SIZE, CPTR)


       ALLOCATE(SHAPE(4))
       SHAPE = (/ NX+2,NY+2,NZ+2,NUM_VARS /)
       CALL C_F_POINTER(CPTR, GRID, SHAPE)
       GRID(0:, 0:, 0:, 1:) => GRID
       DEALLOCATE ( SHAPE ) 
      
       CALL MG_INIT_GRID ( GRID, IERR )

       ALLOCATE ( MAX_SEND_SZ(MAX_NUM_NEIGHBORS), STAT = IERR )
       CALL MG_ASSERT ( IERR, 'MG_INIT: ALLOCATE ( MAX_SEND_SZ )', MAX_NUM_NEIGHBORS )
       ALLOCATE ( MAX_RECV_SZ(MAX_NUM_NEIGHBORS), STAT = IERR )
       CALL MG_ASSERT ( IERR, 'MG_INIT: ALLOCATE ( MAX_RECV_SZ )', MAX_NUM_NEIGHBORS )

       MAX_SEND_SZ(NORTH) = (NX+2)*(NZ+2)
       MAX_SEND_SZ(SOUTH) = (NX+2)*(NZ+2)
       MAX_SEND_SZ(EAST)  = (NY+2)*(NZ+2)
       MAX_SEND_SZ(WEST)  = (NY+2)*(NZ+2)
       MAX_SEND_SZ(BACK)  = (NX+2)*(NY+2)
       MAX_SEND_SZ(FRONT) = (NX+2)*(NY+2)

       MAX_RECV_SZ(NORTH) = (NX+2)*(NZ+2)
       MAX_RECV_SZ(SOUTH) = (NX+2)*(NZ+2)
       MAX_RECV_SZ(EAST)  = (NY+2)*(NZ+2)
       MAX_RECV_SZ(WEST)  = (NY+2)*(NZ+2)
       MAX_RECV_SZ(BACK)  = (NX+2)*(NY+2)
       MAX_RECV_SZ(FRONT) = (NX+2)*(NY+2)
       
       NUM_NEIGHS = 0
       DO I = 1, MAX_NUM_NEIGHBORS
          IF(NEIGHBORS(I) /= -1) THEN
             NUM_NEIGHS = NUM_NEIGHS + 1
          END IF
       END DO


       ALLOCATE ( IDX_NEIGHBORS(MAX_NUM_NEIGHBORS), STAT = IERR )
       CALL MG_ASSERT ( IERR, 'MG_INIT: ALLOCATE ( IDX_NEIGHBORS )', MAX_NUM_NEIGHBORS )
       ALLOCATE ( NEIGHBORS_IDX(NUM_NEIGHS), STAT = IERR )
       CALL MG_ASSERT ( IERR, 'MG_INIT: ALLOCATE ( NEIGHBORS_IDX )', NUM_NEIGHS )
       ALLOCATE ( F_SHAN_NEIGHBORS(NUM_NEIGHS), STAT = IERR )
       CALL MG_ASSERT ( IERR, 'MG_INIT: ALLOCATE ( F_SHAN_NEIGHBORS )', NUM_NEIGHS )
       

       NUM_TYPE = NUM_VARS
       
       NUM_NEIGHS = 0
       TYPE_SEND_SZ = 0
       TYPE_RECV_SZ = 0
       DO I = 1, MAX_NUM_NEIGHBORS
          IF (NEIGHBORS(I) /= -1) THEN
             NUM_NEIGHS                   = NUM_NEIGHS + 1
             IDX_NEIGHBORS(I)             = NUM_NEIGHS
             NEIGHBORS_IDX(NUM_NEIGHS)    = I
             F_SHAN_NEIGHBORS(NUM_NEIGHS) = NEIGHBORS(I)
             TYPE_SEND_SZ = MAX(TYPE_SEND_SZ, MAX_SEND_SZ(I)*C_SIZEOF(SZ_REAL))
             TYPE_RECV_SZ = MAX(TYPE_RECV_SZ, MAX_RECV_SZ(I)*C_SIZEOF(SZ_REAL))
          END IF
       ENDDO


       MAX_NELEM_SEND = 0
       MAX_NELEM_RECV = 0
       TYPE_MAX_NSEND = 0
       TYPE_MAX_NRECV = 0
       DO J = 1, NUM_NEIGHS
          I = NEIGHBORS_IDX(J)
          IF (I == BACK .OR. I == FRONT) THEN
#ifdef USE_NY_SPLIT
             MAX_NELEM_SEND = NY+2
             MAX_NELEM_RECV = NY+2
#else
             MAX_NELEM_SEND = 1
             MAX_NELEM_RECV = 1
#endif
          END IF
          IF (I == EAST .OR. I == WEST) THEN
             MAX_NELEM_SEND = (NZ+2)*(NY+2)
             MAX_NELEM_RECV = (NZ+2)*(NY+2)
          END IF
          IF (I == NORTH .OR. I == SOUTH) THEN
             MAX_NELEM_SEND = NZ+2
             MAX_NELEM_RECV = NZ+2
          ENDIF
          TYPE_MAX_NSEND = MAX(TYPE_MAX_NSEND, MAX_NELEM_SEND)
          TYPE_MAX_NRECV = MAX(TYPE_MAX_NRECV, MAX_NELEM_RECV)
       END DO

       ALLOCATE ( F_SHAN_SEND_SZ(NUM_TYPE), STAT = IERR )
       CALL MG_ASSERT ( IERR, 'MG_INIT: ALLOCATE ( F_SHAN_SEND_SZ )', NUM_TYPE )       
       ALLOCATE ( F_SHAN_RECV_SZ(NUM_TYPE), STAT = IERR )
       CALL MG_ASSERT ( IERR, 'MG_INIT: ALLOCATE ( F_SHAN_RECV_SZ )', NUM_TYPE )

       ALLOCATE ( F_SHAN_NELEM_SEND(NUM_TYPE), STAT = IERR )
       CALL MG_ASSERT ( IERR, 'MG_INIT: ALLOCATE ( F_SHAN_NELEM_SEND )', NUM_TYPE )       
       ALLOCATE ( F_SHAN_NELEM_RECV(NUM_TYPE), STAT = IERR )
       CALL MG_ASSERT ( IERR, 'MG_INIT: ALLOCATE ( F_SHAN_NELEM_RECV )', NUM_TYPE )

       ! --------------------------------------------------------------------
       ! Identical size and number of elements for all types (same mesh dim)
       ! --------------------------------------------------------------------
       DO J = 1, NUM_TYPE
          F_SHAN_SEND_SZ(J)    = TYPE_SEND_SZ
          F_SHAN_RECV_SZ(J)    = TYPE_RECV_SZ
          F_SHAN_NELEM_SEND(J) = TYPE_MAX_NSEND
          F_SHAN_NELEM_RECV(J) = TYPE_MAX_NRECV
       END DO
          
       ! ------------------------------------------------------------
       ! Sets up the communication neighborhood with max message size 
       ! ------------------------------------------------------------
       CALL F_SHAN_INIT_COMM(NEIGHBOR_HOOD_ID, C_LOC(F_SHAN_NEIGHBORS), &
            NUM_NEIGHS, C_LOC(F_SHAN_SEND_SZ), C_LOC(F_SHAN_RECV_SZ), &
            C_LOC(F_SHAN_NELEM_SEND), C_LOC(F_SHAN_NELEM_RECV), NUM_TYPE)       


       CALL MG_SET_COMM_TYPE(IERR)
       CALL MG_ASSERT ( IERR, 'MG_SET_COMM_TYPE', NUM_TYPE)


    END IF
    
    ALLOCATE ( WORK( 0:NX+1, 0:NY+1, 0:NZ+1 ), STAT = IERR )
    CALL MG_ASSERT ( IERR, 'GRID_INIT: ALLOCATE ( WORK )', (NX+2)*(NY+2)*(NZ+2) )

    RETURN

  END SUBROUTINE MG_GRID_INIT

  !  ===================================================================================

  SUBROUTINE MG_INSERT_SPIKE ( SPIKE_NUM, IERR )

    ! Insert heat source (SPIKES) into arrays. 
    ! This vlaue is added to the total heat applied to the variable.

    INTEGER, INTENT(IN) ::  &
         SPIKE_NUM

    INTEGER, INTENT(OUT) ::  &
         IERR                    ! Return status.

    ! ------------------
    ! Local Declarations
    ! ------------------

    INTEGER ::  &
         IX,      & ! Local index
         IY,      & ! Local index
         IZ,      & ! Local index
         I          ! Counter

    ! ---------------------
    ! Executable Statements
    ! ---------------------

    IERR = 0

    IF ( MYPE == ( SPIKE_LOC ( 0, SPIKE_NUM ) ) ) THEN

       IX = SPIKE_LOC ( 1, SPIKE_NUM )
       IY = SPIKE_LOC ( 2, SPIKE_NUM )
       IZ = SPIKE_LOC ( 3, SPIKE_NUM )

       DO I = 1, NUM_VARS
          GRID( IX, IY, IZ, I ) = SPIKES( I, SPIKE_NUM )
       ENDDO
       IF ( NUM_VARS > 40 ) THEN
          IERR = -1
          CALL MG_ASSERT ( IERR, 'GRID_INIT: TOO MANY VARS', NUM_VARS )
       END IF
    END IF

    IF ( MYPE == ROOT_PE ) THEN
       DO I = 1, NUM_VARS
          SOURCE_TOTAL( I ) = SOURCE_TOTAL( I ) + SPIKES( I, SPIKE_NUM )
       END DO
    END IF

  END SUBROUTINE MG_INSERT_SPIKE

  !  ===================================================================================

  SUBROUTINE MG_GRID_DEALLOC ( IERR )

    USE MG_CONSTANTS_MOD

    IMPLICIT NONE

    ! Variable Declarations

    INTEGER :: &
         IERR                            ! Return status

    ! ---------------------
    ! Executable Statements
    ! ---------------------

    IERR = 0

    ! ---------------------
    ! Free Data Segment
    ! ---------------------
    CALL F_SHAN_FREE_SHARED(SEGMENT_ID)

    ! ---------------------
    ! Free Neighborhood
    ! ---------------------
    CALL F_SHAN_FREE_COMM(NEIGHBOR_HOOD_ID)

  END SUBROUTINE MG_GRID_DEALLOC

  !  ===================================================================================

  SUBROUTINE MG_ASSERT ( IERR, ERROR_MSG, INFO )

    USE MG_CONSTANTS_MOD

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::    &
         ERROR_MSG

    INTEGER, INTENT(IN) ::             &
         IERR,                           & ! Error code.
         INFO

    ! ---------------------
    ! Executable Statements
    ! ---------------------

#if defined _MG_MPI

    CHARACTER*(MPI_MAX_ERROR_STRING) STRING
    INTEGER RESULTLEN, IERROR

    IF ( IERR /= 0 ) THEN
       WRITE(*,80) MYPE, TRIM(ERROR_MSG), IERR, INFO
       call MPI_ERROR_STRING(IERR, STRING, RESULTLEN, IERROR)
       WRITE(*,81) STRING
       CALL MPI_ABORT ( MPI_COMM_MG, -1, IERR )
    END IF

80  FORMAT ( '** Error ** [pe ', I5, '] ', A40, '; CODE = ', I7, &
         '. Additional info:', I4 )
81  FORMAT ( 'MPI error message: ', A80)

#else

    IF ( IERR /= 0 ) THEN
       WRITE(*,90) TRIM(ERROR_MSG), IERR, INFO
       STOP
    END IF

#endif

90  FORMAT ( '** Error ** ', A40, '; CODE = ', I7, '. Additional info:', I4 )


  END SUBROUTINE MG_ASSERT

  !  ===================================================================================

  SUBROUTINE MG_INIT_GRID( GRID, IERR )

    IMPLICIT NONE

    INTEGER, INTENT(OUT) ::  &
         IERR                    ! Return status.

    REAL(KIND=MG_REAL), DIMENSION(0:NX+1, 0:NY+1, 0:NZ+1), INTENT(OUT) :: &
         GRID

    ! ---------------
    ! Local Variables
    ! ---------------

    INTEGER :: I, J, K

    ! ---------------------
    ! Executable Statements
    ! ---------------------

    ! Explicit loops enable first touch affinity for OpenMP.

    IF ( DEBUG_GRID == 1 ) THEN

       !$OMP PARALLEL DO
       DO K = 0, NZ+1
          DO J = 0, NY+1
             DO I = 0, NX+1
                GRID(I, J, K) = 0.0
             END DO
          END DO
       END DO
       !$OMP END PARALLEL DO

    ELSE

       !$OMP PARALLEL DO
       DO K = 0, NZ+1
          DO J = 0, NY+1
             DO I = 0, NX+1
                CALL RANDOM_NUMBER( GRID(I, J, K) )
             END DO
          END DO
       END DO
       !$OMP END PARALLEL DO

    END IF

    RETURN

  END SUBROUTINE MG_INIT_GRID


  
  !  ===================================================================================


  INTEGER(KIND=MG_INT8) FUNCTION GET_OFFSET(I, J, K)

    REAL(KIND=MG_REAL)    :: SZ_REAL
    INTEGER, INTENT(IN) ::  I, J, K
    GET_OFFSET = (I + (NX+2)*J + (NX+2)*(NY+2)*K) * C_SIZEOF(SZ_REAL)    
    RETURN
    
  END FUNCTION GET_OFFSET
  
  !  ===================================================================================

  DOUBLE PRECISION FUNCTION MG_TIMER ()

#if defined _MG_MPI
    include 'mpif.h'

    MG_TIMER = MPI_WTIME ()

#else
    INTEGER COUNT_1, COUNT_RATE, COUNT_MAX

    CALL SYSTEM_CLOCK (COUNT_1, COUNT_RATE, COUNT_MAX)
    MG_TIMER = COUNT_1 * 1.0 / COUNT_RATE

#endif

    RETURN

  END FUNCTION MG_TIMER
  !  ===================================================================================

  DOUBLE PRECISION FUNCTION MG_COMPUTE_STDDEV ( VALUES, MEAN )

    DOUBLE PRECISION :: &
         VALUES(*), MEAN

    INTEGER :: I 
    DOUBLE PRECISION :: TMP

    TMP = 0.0
    DO I = 1, NUMPES
       TMP = TMP + ( VALUES(I) - MEAN )**2
    END DO

    MG_COMPUTE_STDDEV = SQRT ( TMP / REAL( NUMPES ))

    RETURN

  END FUNCTION MG_COMPUTE_STDDEV

  !  ===================================================================================

END MODULE MG_UTILS_MOD
