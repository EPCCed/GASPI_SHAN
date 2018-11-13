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

#include <GASPI.h>

#include "SHAN_segment.h"
#include "SHAN_comm.h"
#include "SHAN_type.h"

#include "shan_exchange.h"
#include "gaspi_util.h"
#include "shan_util.h"
#include "shan_core.h"
#include "assert.h"


int shan_get_shared_type(type_local_t *type_info
			 , shan_neighborhood_t *neighborhood_id
			 , int local_rank
			 , int num_neighbors
			 , int type_id
			 )
{  
  void *shm_ptr;
  shan_segment_t *shared_segment = &(neighborhood_id->shared_segment);
  shan_get_shared_ptr(shared_segment
		      , local_rank
		      , &shm_ptr);

  shan_comm_get_type(type_info
		     , shm_ptr
		     , num_neighbors
		     , &(neighborhood_id->type_element[type_id])
		     ); 
  
  return SHAN_SUCCESS;
}


int shan_comm_get_type(type_local_t *type_info
		       , void *shm_ptr
		       , int num_neighbors
		       , shan_element_t *type_element
		       )
{
    long typeOffset          = num_neighbors * type_element->elemOffset;
    int max_nelem_send       = type_element->max_nelem_send;
    int max_nelem_recv       = type_element->max_nelem_recv;
    long const maxSz = 
	num_neighbors * MAX_SHARED_NOTIFICATION * sizeof(shan_notification_t)
	+ 4 * num_neighbors * sizeof(int)
	+ num_neighbors * (max_nelem_send + max_nelem_recv) * sizeof(long);
  
    type_info->nid            = (shan_notification_t*) ((char *) shm_ptr + typeOffset);
    typeOffset              += sizeof(shan_notification_t) * num_neighbors * MAX_SHARED_NOTIFICATION;
    type_info->nelem_send     = (int*) ((char*) shm_ptr + typeOffset);
    typeOffset              += num_neighbors * sizeof(int);
    type_info->nelem_recv     = (int*) ((char*) shm_ptr + typeOffset);
    typeOffset              += num_neighbors * sizeof(int);
    type_info->send_sz        = (int*) ((char*) shm_ptr + typeOffset);
    typeOffset              += num_neighbors * sizeof(int);
    type_info->recv_sz        = (int*) ((char*) shm_ptr + typeOffset);
    typeOffset              += num_neighbors * sizeof(int);
    type_info->send_offset    = (long*) ((char*) shm_ptr + typeOffset);
    typeOffset              += max_nelem_send * num_neighbors * sizeof(long);
    type_info->recv_offset    = (long*) ((char*) shm_ptr + typeOffset);
    typeOffset              += max_nelem_recv * num_neighbors * sizeof(long);

    ASSERT(typeOffset == num_neighbors * type_element->elemOffset + maxSz);

    return SHAN_SUCCESS;
  
}

int shan_comm_type_offset(shan_neighborhood_t *neighborhood_id
			  , int type_id
			  , int **nelem_send
			  , int **nelem_recv
			  , int **send_sz
			  , int **recv_sz
			  , long **send_offset
			  , long **recv_offset
			  )
{
    int const iProcLocal = neighborhood_id->iProcLocal;
    int const num_neighbors = neighborhood_id->num_neighbors;
    type_local_t type_info;
    shan_get_shared_type(&type_info
			 , neighborhood_id
			 , iProcLocal
			 , num_neighbors
			 , type_id
	);
            
    *nelem_send    = type_info.nelem_send;
    *nelem_recv    = type_info.nelem_recv;
    *send_sz       = type_info.send_sz;
    *recv_sz       = type_info.recv_sz;
    *send_offset   = type_info.send_offset;
    *recv_offset   = type_info.recv_offset;
        
    return SHAN_SUCCESS;

}


int shan_comm_type_free(shan_segment_t *type_segment)
{
  int res = shan_free_shared(type_segment);
  ASSERT(res == SHAN_SUCCESS);

  return SHAN_SUCCESS;
}



void shan_comm_get_local(shan_neighborhood_t *neighborhood_id
			 , shan_segment_t *data_segment
			 , int const type_id
			 , int const idx
			 )
{
  int const iProcLocal = neighborhood_id->iProcLocal;
  int iProcRemote, i;

  int const rank = neighborhood_id->neighbors[idx];
  ASSERT ((iProcRemote = shan_comm_local_rank(neighborhood_id
					      , rank
					      )) != -1);  
  
  int const RemoteCommIdx = neighborhood_id->RemoteCommIndex[idx];
  int const RemoteNumNeighbors = neighborhood_id->RemoteNumNeighbors[idx];

  type_local_t type_info_src;
  shan_get_shared_type(&type_info_src
		       , neighborhood_id
		       , iProcRemote
		       , RemoteNumNeighbors
		       , type_id
		       );
  
  int  src_nelem_send   = type_info_src.nelem_send[RemoteCommIdx];
  int  src_send_sz      = type_info_src.send_sz[RemoteCommIdx];
  long *src_send_offset = type_info_src.send_offset
      + RemoteCommIdx * neighborhood_id->type_element[type_id].max_nelem_send;		      

  type_local_t type_info_dest;
  shan_get_shared_type(&type_info_dest
		       , neighborhood_id
		       , iProcLocal
		       , neighborhood_id->num_neighbors
		       , type_id
		       );

  int  dest_recv_sz      = type_info_dest.recv_sz[idx];
  long *dest_recv_offset = type_info_dest.recv_offset
    + idx * neighborhood_id->type_element[type_id].max_nelem_recv;    
  ASSERT(src_send_sz    == dest_recv_sz);

  void *send_ptr, *recv_ptr;
  shan_get_shared_ptr(data_segment
		      , iProcRemote
		      , &send_ptr);

  shan_get_shared_ptr(data_segment
		      , iProcLocal
		      , &recv_ptr);

  for (i = 0; i < src_nelem_send; ++i)
    {
      void *restrict src  = (char*) send_ptr + src_send_offset[i];
      void *restrict dest = (char* )recv_ptr + dest_recv_offset[i];
      memcpy(dest, src, src_send_sz);
    }

#ifdef USE_VARIABLE_MESSAGE_LEN
  type_info_dest.nelem_recv[idx] = src_nelem_send; 
  type_info_dest.recv_sz[idx] = src_send_sz; 
#else
  ASSERT(type_info_dest.nelem_recv[idx] == src_nelem_send);
  ASSERT(type_info_dest.recv_sz[idx] == src_send_sz);
#endif  

  /*
   * this triggers notification for 'have read',
   * required in checking 'send' completion.
   */
  shan_increment_local(neighborhood_id
		       , type_id
		       , neighborhood_id->num_neighbors + idx
		       );


  ++(neighborhood_id->type_element[type_id].local_recv_count[idx]);
  
}


void shan_comm_get_remote(shan_neighborhood_t *neighborhood_id
			 , shan_segment_t *data_segment
			 , int const type_id
			 , int const idx
			 )
{
    int i;
    int const iProcLocal    = neighborhood_id->iProcLocal;
    int const num_neighbors = neighborhood_id->num_neighbors;
  
    int const sid 
	= (neighborhood_id->type_element[type_id].local_recv_count[idx]) % 2;  
    shan_remote_t *const remote_segment = &(neighborhood_id->remote_segment);

    type_local_t type_info;
    shan_get_shared_type(&type_info
			 , neighborhood_id
			 , iProcLocal
			 , num_neighbors
			 , type_id
	);

    int const recv_sz = type_info.recv_sz[idx];	  
    int const nelem_recv = type_info.nelem_recv[idx];	  
    long *const recv_offset = type_info.recv_offset
	+ idx * neighborhood_id->type_element[type_id].max_nelem_recv;
    
    long comm_buffer_offset = num_neighbors * neighborhood_id->type_element[type_id].RecvOffset[sid]
        + idx * neighborhood_id->type_element[type_id].maxRecvSz;    
    void *comm_ptr = (char*) remote_segment->shan_ptr + comm_buffer_offset;
    
    void *data_ptr;
    shan_get_shared_ptr(data_segment
			, iProcLocal
			, &data_ptr);
    
    for (i = 0; i < nelem_recv; ++i)
    {
	void *restrict src  = (char*) comm_ptr 
	    + NELEM_COMM_HEADER * sizeof(int) + i * recv_sz;
	void *restrict dest = (char*) data_ptr + recv_offset[i];
	memcpy(dest, src, recv_sz);
    }
    
    ++(neighborhood_id->type_element[type_id].local_recv_count[idx]);
}

