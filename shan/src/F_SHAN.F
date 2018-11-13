MODULE F_SHAN

  use, intrinsic :: ISO_C_BINDING

  interface
     subroutine F_SHAN_ALLOC_SHARED(segment_id &
          , dataSz &
          , shm_ptr &
          ) &
          bind(C, name="f_shan_alloc_shared")
       import
       integer(c_int), value  :: segment_id
       integer(c_long), value :: dataSz
       type(c_ptr) :: shm_ptr
     end subroutine F_SHAN_ALLOC_SHARED
  end interface

  interface
     subroutine F_SHAN_FREE_SHARED(segment_id &
          ) &
          bind(C, name="f_shan_free_shared")
       import
       integer(c_int), value  :: segment_id
     end subroutine F_SHAN_FREE_SHARED
  end interface

  interface
     subroutine F_SHAN_INIT_COMM(neighbor_hood_id &
          , neighbors &
          , num_neighbors &
          , maxSendSz &
          , maxRecvSz &
          , max_nelem_send &
          , max_nelem_recv &
          , num_type &
          ) &
          bind(C, name="f_shan_init_comm")
       import
       integer(c_int), value :: neighbor_hood_id
       type(c_ptr), value    :: neighbors
       integer(c_int), value :: num_neighbors
       type(c_ptr), value    :: maxSendSz
       type(c_ptr), value    :: maxRecvSz
       type(c_ptr), value    :: max_nelem_send
       type(c_ptr), value    :: max_nelem_recv
       integer(c_int), value :: num_type
     end subroutine F_SHAN_INIT_COMM
  end interface

  interface
     subroutine F_SHAN_FREE_COMM(neighbor_hood_id) &
          bind(C, name="f_shan_free_comm")
       import
       integer(c_int), value :: neighbor_hood_id
     end subroutine F_SHAN_FREE_COMM
  end interface


  interface
     subroutine F_SHAN_TYPE_OFFSET(neighbor_hood_id &
          , type_id &
          , nelem_send &
          , nelem_recv &
          , send_sz &
          , recv_sz &
          , send_idx &
          , recv_idx &
          ) &
          bind(C, name="f_shan_type_offset")
       import
       integer(c_int), value :: neighbor_hood_id
       integer(c_int), value :: type_id
       type(c_ptr)           :: nelem_send
       type(c_ptr)           :: nelem_recv
       type(c_ptr)           :: send_sz
       type(c_ptr)           :: recv_sz
       type(c_ptr)           :: send_idx
       type(c_ptr)           :: recv_idx
     end subroutine F_SHAN_TYPE_OFFSET
  end interface

  interface
     subroutine F_SHAN_COMM_NOTIFY_OR_WRITE(neighbor_hood_id &
          , segment_id &
          , type_id &
          , idx &
          ) &
          bind(C, name="f_shan_comm_notify_or_write")
       import
       integer(c_int), value :: neighbor_hood_id
       integer(c_int), value :: segment_id
       integer(c_int), value :: type_id
       integer(c_int), value :: idx
     end subroutine F_SHAN_COMM_NOTIFY_OR_WRITE
  end interface


  interface
     subroutine F_SHAN_COMM_WAIT4ALL(neighbor_hood_id &
          , segment_id &
          , type_id &
          ) &
          bind(C, name="f_shan_comm_wait4All")
       import
       integer(c_int), value :: neighbor_hood_id
       integer(c_int), value :: segment_id
       integer(c_int), value :: type_id
     end subroutine F_SHAN_COMM_WAIT4ALL
  end interface


  interface
     subroutine F_SHAN_COMM_WAIT4ALLSEND(neighbor_hood_id &
          , type_id &
          ) &
          bind(C, name="f_shan_comm_wait4AllSend")
       import
       integer(c_int), value :: neighbor_hood_id
       integer(c_int), value :: type_id
     end subroutine F_SHAN_COMM_WAIT4ALLSEND
  end interface


  interface
     subroutine F_SHAN_COMM_WAIT4ALLRECV(neighbor_hood_id &
          , segment_id &
          , type_id &
          ) &
          bind(C, name="f_shan_comm_wait4AllRecv")
       import
       integer(c_int), value :: neighbor_hood_id
       integer(c_int), value :: segment_id
       integer(c_int), value :: type_id
     end subroutine F_SHAN_COMM_WAIT4ALLRECV
  end interface

  
END MODULE F_SHAN
      

