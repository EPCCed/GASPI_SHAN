      program shmem_mpi

! this is a hello world example based on code examples
! by Jeff Hammond and Jonathan Vincent (big thanks to them)
! in this example master allocates and fills an array with "twos" and then
! each rank prints the sum of all ements, each rank should print 200.0

      use mpi

      use, intrinsic :: iso_c_binding

      implicit none

      integer, parameter :: nsize = 100
      real(8), pointer   :: array(:)
      integer            :: num_procs
      integer            :: ierr
      integer            :: irank
      integer            :: win
      integer            :: comm = MPI_COMM_WORLD
      integer            :: disp_unit
      type(c_ptr)        :: cp1
      type(c_ptr)        :: cp2

      integer(MPI_ADDRESS_KIND) :: win_size
      integer(MPI_ADDRESS_KIND) :: segment_size

      call MPI_Init(ierr)
      call MPI_Comm_size(comm, num_procs, ierr)
      call MPI_Comm_rank(comm, irank, ierr)

      disp_unit = sizeof(1.0d0)

      if (irank == 0) then
         win_size = nsize*disp_unit
      else
         win_size = 0
      endif

      call MPI_Win_allocate_shared(win_size, disp_unit, MPI_INFO_NULL, comm, cp1, win, ierr)
      call MPI_Win_fence(0, win, ierr)

      if (irank == 0) then
         call MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, MPI_MODE_NOCHECK, win, ierr)
         call c_f_pointer(cp1, array, [nsize])
         array = 2.0d0
         call MPI_Win_unlock(0, win, ierr)
         nullify(array)
      endif

      call MPI_Win_fence(0, win, ierr)
      call MPI_Win_shared_query(win, 0, segment_size, disp_unit, cp2, ierr)

      call c_f_pointer(cp2, array, [nsize])

      print *, 'irank=', irank, 'sum=', sum(array)
      nullify(array)

      call MPI_Win_fence(MPI_MODE_NOSUCCEED, win, ierr)
      call MPI_Win_free(win, ierr)
      call MPI_Finalize(ierr)

      end program
      
