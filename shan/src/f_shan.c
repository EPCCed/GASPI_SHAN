/*
    Copyright (c) T-Systems SfR, C.Simmendinger <christian.simmendinger@t-systems.com>, 2018

    This file is part of SHAN.

    SHAN is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
	    the Free Software Foundation, either version 3 of the License, or
	        (at your option) any later version.

    SHAN is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
	    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	        GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
        along with SHAN.  If not, see <https://www.gnu.org/licenses/>.
	
*/


#include <mpi.h>
#include <fcntl.h>
#include <sys/shm.h>
#include <sys/mman.h>
#include <unistd.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>

#include "SHAN_segment.h"
#include "SHAN_comm.h"
#include "SHAN_type.h"
#include "F_SHAN.h"

#include "assert.h"
#include "shan_util.h"
 
#define MAX_NEIGHBOR_HOOD 32
#define MAX_SEGMENT 32

static shan_neighborhood_t neighbor_hood[MAX_NEIGHBOR_HOOD];
static shan_segment_t data_segment[MAX_SEGMENT];


void f_shan_alloc_shared(const int segment_id
			, const long dataSz
			, void **restrict shm_ptr
			)
{
  int res, iProcLocal;
  MPI_Comm MPI_COMM_SHM;
  MPI_Comm_split_type (MPI_COMM_WORLD
		       , MPI_COMM_TYPE_SHARED
		       , 0
		       , MPI_INFO_NULL
		       , &MPI_COMM_SHM
		       );
  MPI_Comm_rank(MPI_COMM_SHM, &iProcLocal);

  shan_segment_t *dataSegment  = &data_segment[segment_id];
  res = shan_alloc_shared(dataSegment
			  , segment_id
			  , SHAN_DATA
			  , dataSz			  
			  , MPI_COMM_SHM 
			  );
  ASSERT(res == SHAN_SUCCESS);

  res = shan_get_shared_ptr(dataSegment
			    , iProcLocal
			    , shm_ptr
			    );
  ASSERT(res == SHAN_SUCCESS);
}


void f_shan_free_shared(const int segment_id)
{
  int res;
  shan_segment_t *dataSegment  = &data_segment[segment_id];
  res = shan_free_shared(dataSegment);
  ASSERT(res == SHAN_SUCCESS);

  MPI_Comm_free(&(dataSegment ->MPI_COMM_SHM));

}


void f_shan_free_comm(const int neighbor_hood_id)
{
  shan_neighborhood_t *ngbSegment = &neighbor_hood[neighbor_hood_id];
  int res = shan_comm_free_comm(ngbSegment);
  ASSERT(res == SHAN_SUCCESS);

  MPI_Comm_free(&(ngbSegment->MPI_COMM_SHM));
    
}

void f_shan_init_comm(const int neighbor_hood_id
		      , void *neighbors
		      , int num_neighbors
		      , void *maxSendSz
		      , void *maxRecvSz
		      , void *max_nelem_send
		      , void *max_nelem_recv
		      , int num_type
		      )
{
  int res, iProcLocal;
  MPI_Comm MPI_COMM_SHM;
  MPI_Comm_split_type (MPI_COMM_WORLD
		       , MPI_COMM_TYPE_SHARED
		       , 0
		       , MPI_INFO_NULL
		       , &MPI_COMM_SHM
		       );
  MPI_Comm_rank(MPI_COMM_SHM, &iProcLocal);

  shan_neighborhood_t *ngbSegment = &neighbor_hood[neighbor_hood_id];

  res = shan_comm_init_comm(ngbSegment
			    , neighbor_hood_id
			    , (int*) neighbors
			    , num_neighbors
			    , (long*) maxSendSz
			    , (long*) maxRecvSz			    
			    , (int*) max_nelem_send
			    , (int*) max_nelem_recv
			    , num_type
			    , MPI_COMM_SHM
			    , MPI_COMM_WORLD
			    );  
  ASSERT(res == SHAN_SUCCESS);  
}

void f_shan_type_offset(const int neighbor_hood_id
			, const int type_id
			, void **nelem_send
			, void **nelem_recv
			, void **send_sz
			, void **recv_sz
			, void **send_idx
			, void **recv_idx
			)
{
  shan_neighborhood_t *ngbSegment = &neighbor_hood[neighbor_hood_id];
  int res = shan_comm_type_offset(ngbSegment
			       , type_id
			       , (int**) nelem_send
			       , (int**) nelem_recv
			       , (int**) send_sz
			       , (int**) recv_sz
			       , (long**) send_idx
			       , (long**) recv_idx
			       );
  ASSERT(res == SHAN_SUCCESS);  
}


void f_shan_comm_wait4All(const int neighbor_hood_id  
			  , const int segment_id
			  , const int type_id
			 )
{
  shan_neighborhood_t *ngbSegment = &neighbor_hood[neighbor_hood_id];
  shan_segment_t *dataSegment  = &data_segment[segment_id];

  int res = shan_comm_wait4All(ngbSegment
			       , dataSegment
			       , type_id
			       );
 ASSERT(res == SHAN_SUCCESS);
}



void f_shan_comm_wait4AllSend(const int neighbor_hood_id  
			      , const int type_id
    )
{
    shan_neighborhood_t *ngbSegment = &neighbor_hood[neighbor_hood_id];
    
    int res = shan_comm_wait4AllSend(ngbSegment
				     , type_id
	);
    ASSERT(res == SHAN_SUCCESS);
}


void f_shan_comm_wait4AllRecv(const int neighbor_hood_id  
			      , const int segment_id
			      , const int type_id
    )
{
    shan_neighborhood_t *ngbSegment = &neighbor_hood[neighbor_hood_id];
    shan_segment_t *dataSegment  = &data_segment[segment_id];
    
    int res = shan_comm_wait4AllRecv(ngbSegment
				     , dataSegment
				     , type_id
	);
    ASSERT(res == SHAN_SUCCESS);
}


void f_shan_comm_notify_or_write(const int neighbor_hood_id
				 , const int segment_id
				 , const int type_id
				 , int idx
    )
{
  shan_neighborhood_t *ngbSegment = &neighbor_hood[neighbor_hood_id];
  shan_segment_t *dataSegment  = &data_segment[segment_id];

  int res = shan_comm_notify_or_write(ngbSegment
				      , dataSegment
				      , type_id
				      , idx
      );
  ASSERT(res == SHAN_SUCCESS);  
  
}




