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
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <xmmintrin.h>

#include <GASPI.h>

#include "SHAN_segment.h"
#include "SHAN_comm.h"
#include "SHAN_type.h"

#include "shan_core.h"
#include "shan_exchange.h"
#include "gaspi_util.h"
#include "shan_util.h"
#include "assert.h"


#define GET_NOTIFICATION_ID(sid, num_type, type_id, num_neighbors, idx) \
  ((sid) * ((num_type) * (num_neighbors)) + (type_id) * (num_neighbors) + (idx))  

#ifndef USE_NOCOS
#define GASPI_PROC_LOCAL 0
#endif

void shan_test_shared(shan_neighborhood_t *const neighborhood_id
		      , int const rank_local
		      , int const num_neighbors
		      , int const type_id
		      , int const idx
		      , int *rval
    )
{
    shan_segment_t *const shared_segment = &(neighborhood_id->shared_segment);
    long const typeOffset 
	= num_neighbors * neighborhood_id->type_element[type_id].elemOffset;
  
    void *shm_ptr = NULL; 
    shan_get_shared_ptr(shared_segment
			, rank_local
			, &shm_ptr
	);  
    shan_notify_test_shared((shan_notification_t *) ((char*) shm_ptr + typeOffset)
			    , idx
			    , rval
	);
}

int shan_increment_local(shan_neighborhood_t *const neighborhood_id				  
			 , int const type_id
			 , int const idx
			 ) 
{
    shan_segment_t *const shared_segment = &(neighborhood_id->shared_segment);
    long const typeOffset 
	= neighborhood_id->num_neighbors * neighborhood_id->type_element[type_id].elemOffset;
  
    void *shm_ptr = NULL; 
    shan_get_shared_ptr(shared_segment
			, neighborhood_id->iProcLocal
			, &shm_ptr
	);
    shan_notify_increment_shared((shan_notification_t *)((char*) shm_ptr + typeOffset)
				 , idx
				 , 1
	);        
    return SHAN_SUCCESS;
}


static int shan_alloc_remote(shan_remote_t * const segment   
			     , int const shan_id
			     , const long dataSz
    )
{
    ASSERT(segment != NULL);
    segment->shan_id = shan_id;
    segment->dataSz = dataSz;
    if (dataSz >0)
      {
	  int const page_size = sysconf (_SC_PAGESIZE);
	  void *rem_ptr = NULL; 
	  int res = posix_memalign(&rem_ptr, page_size, dataSz);
	  ASSERT(res == 0);
	  ASSERT(rem_ptr != NULL);

	  segment->shan_ptr = rem_ptr;
      }
    else
      {
	segment->shan_ptr = NULL;
      }
    return SHAN_SUCCESS; 
}


static int shan_free_remote(shan_remote_t * const segment)
{
    ASSERT(segment != NULL);
    check_free(segment->shan_ptr);
  
    return SHAN_SUCCESS; 
}


static void bind_to_segment(shan_neighborhood_t *const neighborhood_id
			    , const gaspi_segment_id_t segment_id
			    )
{
    ASSERT(neighborhood_id != NULL);
    /* 
     * we bind to the entire comm segment at rank 0 
     * includes all (page aligned) send / recv buffers 
     * from all local ranks
     */
  
    shan_remote_t *remote_segment = &(neighborhood_id->remote_segment);
    gaspi_pointer_t seg_ptr = (gaspi_pointer_t *) remote_segment->shan_ptr;
    gaspi_size_t CommSz = (gaspi_size_t) remote_segment->dataSz;
	 
    /* 
     * bind to shared GASPI segment 
     */
    if (CommSz > 0)
    {
	SUCCESS_OR_DIE (gaspi_segment_bind( (gaspi_segment_id_t) segment_id 
					    , seg_ptr
					    , CommSz
					    , GASPI_PROC_LOCAL
			    ));
    }

    MPI_Barrier(neighborhood_id->MPI_COMM_ALL);

    if (CommSz > 0)
    {
	int i;    
	for (i = 0; i < neighborhood_id->num_neighbors; ++i)
	{
	    int const rank = neighborhood_id->neighbors[i];
	    if(shan_comm_local_rank(neighborhood_id
				    , rank
		   ) == -1)
	    {    
		/* 
		 * connect to comm partner and register
		 */
		SUCCESS_OR_DIE( gaspi_connect (rank, GASPI_BLOCK));
		SUCCESS_OR_DIE( gaspi_segment_register((gaspi_segment_id_t) segment_id
						       , rank
						       , GASPI_BLOCK
				    ));
	    }
	}
    }

    MPI_Barrier(neighborhood_id->MPI_COMM_ALL);

}

static void shan_comm_alloc_comm(shan_neighborhood_t *const neighborhood_id)
{
    int i;
    ASSERT(neighborhood_id != NULL);

    gaspi_number_t  segment_max;
    SUCCESS_OR_DIE (gaspi_segment_max (&segment_max));

    gaspi_number_t segment_num;
    SUCCESS_OR_DIE(gaspi_segment_num(&segment_num));  

    gaspi_segment_id_t segment_id = neighborhood_id->neighbor_hood_id;
    ASSERT(segment_id < segment_max);

    /*
     * page alignment
     */


    if (segment_num > 0)
    {
	gaspi_segment_id_t *segment_list =
	    check_malloc(segment_num * sizeof(gaspi_segment_id_t));
	SUCCESS_OR_DIE(gaspi_segment_list (segment_num, segment_list));
	for (i = 0; i < (int) segment_num; ++i)
	{
	    ASSERT(segment_list[i] != segment_id);
	}
	check_free(segment_list);
    }

    int res;
    /*
     * alloc comm space remote
     */
    long sz_remote = neighborhood_id->remoteSz;
    shan_remote_t *remote_segment  = &(neighborhood_id->remote_segment);


    res = shan_alloc_remote(remote_segment
			    , segment_id
			    , sz_remote
	);      
    ASSERT(res == SHAN_SUCCESS);
    /*
     * bind remote comm as GASPI segment
     */
    bind_to_segment(neighborhood_id, segment_id);

    if (sz_remote > 0)
    {
	memset((char*) remote_segment->shan_ptr
	       , 0
	       , sz_remote
	    );
    }

    gaspi_number_t notification_num;
    SUCCESS_OR_DIE(gaspi_notification_num (&notification_num));
    int max_notifications = 2 * neighborhood_id->num_neighbors
	* neighborhood_id->num_type;

    ASSERT(max_notifications < (int) notification_num);
    for (i = 0; i < max_notifications ; ++i)
    {
	gaspi_notification_t nval;
	SUCCESS_OR_DIE(gaspi_notify_reset (remote_segment->shan_id
					   , (gaspi_notification_id_t) i
					   , &nval
			   ));
    }

}

int shan_comm_local_rank(shan_neighborhood_t * const neighborhood_id
		    , int const rank
		    )
{
    int val = -1;
    ASSERT(rank >= 0);
    ASSERT(rank < neighborhood_id->nProcGlobal);
  
    if (neighborhood_id->master == neighborhood_id->remote_master[rank])
    {
	val = rank % neighborhood_id->nProcLocal;
    }    
    return val;
}


#ifndef DEBUG
static void shan_negotiate_meta_data(shan_neighborhood_t * const neighborhood_id)
{
    int i;
    int num_neighbors = neighborhood_id->num_neighbors;
    /*
     * set rank local offsets and notifications
     */
    neighborhood_id->RemoteNumNeighbors = check_malloc(num_neighbors * sizeof(int));
    neighborhood_id->RemoteCommIndex = check_malloc(num_neighbors * sizeof(int));

    for (i = 0; i < num_neighbors; ++i)
    {
	neighborhood_id->RemoteNumNeighbors[i] = 0;
	neighborhood_id->RemoteCommIndex[i] = -1;   
    }
  
    /* 
     * exchange remote comm idx
     */     
    for (i = 0; i < num_neighbors; ++i)
    {
	int const rank = neighborhood_id->neighbors[i];
	local_exchange_comm(neighborhood_id->iProcGlobal
			    , rank
			    , sizeof(int)
			    , &i
			    , &(neighborhood_id->RemoteCommIndex[i])
			    , neighborhood_id->MPI_COMM_ALL
	    );
      
	local_exchange_comm(neighborhood_id->iProcGlobal
			    , rank
			    , sizeof(int)
			    , &num_neighbors
			    , &(neighborhood_id->RemoteNumNeighbors[i])
			    , neighborhood_id->MPI_COMM_ALL
	    );
    }     
}
#else
static void shan_negotiate_meta_data(shan_neighborhood_t * const neighborhood_id)
{
    int i;
    int num_neighbors = neighborhood_id->num_neighbors;
    /* 
     * exchange remote num_neighbors
     */	  
    int *tmp1 = check_malloc (neighborhood_id->nProcGlobal * sizeof(int));  
    MPI_Allgather(&(neighborhood_id->num_neighbors)
		  , 1
		  , MPI_INT
		  , tmp1
		  , 1
		  , MPI_INT
		  , neighborhood_id->MPI_COMM_ALL
	);  

    neighborhood_id->RemoteNumNeighbors = check_malloc(num_neighbors * sizeof(int));
    for (i = 0; i < num_neighbors; ++i)
    {
	int const rank = neighborhood_id->neighbors[i];
	neighborhood_id->RemoteNumNeighbors[i] = tmp1[rank];
    }	  

    /* 
     * exchange remote neighbors
     */	  
    int *tmp2 = check_malloc (neighborhood_id->nProcGlobal * sizeof(int));  
    int ncomm = 0;
    for (i = 0; i < neighborhood_id->nProcGlobal; ++i)
    {
	tmp2[i] = ncomm;
	ncomm += tmp1[i];
    }	  

    int *tmp3 = check_malloc (ncomm * sizeof(int));  
    MPI_Allgatherv(neighborhood_id->neighbors
		   , neighborhood_id->num_neighbors
		   , MPI_INT
		   , tmp3
		   , tmp1
		   , tmp2
		   , MPI_INT
		   , neighborhood_id->MPI_COMM_ALL
	);  


    /* 
     * exchange remote comm idx
     */	  
    neighborhood_id->RemoteCommIndex = check_malloc(num_neighbors * sizeof(int));
    for (i = 0; i < num_neighbors; ++i)
    {
	neighborhood_id->RemoteCommIndex[i] = -1;	  
    }
  
    for (i = 0; i < num_neighbors; ++i)
    {
	int const rank = neighborhood_id->neighbors[i];
	int j, bdir = 0;
	for (j = 0; j < tmp1[rank]; ++j)
	{
	    if (tmp3[tmp2[rank]+j] == neighborhood_id->iProcGlobal)
	    {
		neighborhood_id->RemoteCommIndex[i] = j;	  
		bdir = 1;
		break;
	    }
	}
	ASSERT(bdir == 1);
    }
	    
    check_free(tmp3);
    check_free(tmp2);
    check_free(tmp1);

}
#endif



int shan_comm_free_comm(shan_neighborhood_t *const neighborhood_id)
{
    int i;
    for (i = 0; i < neighborhood_id->num_type; ++i)
    {
	check_free(neighborhood_id->type_element[i].local_recv_count);
	check_free(neighborhood_id->type_element[i].local_send_count);
	check_free(neighborhood_id->type_element[i].local_ack_count);
    }

    check_free(neighborhood_id->type_element);

    check_free(neighborhood_id->neighbors);
    check_free(neighborhood_id->remote_master);
    check_free(neighborhood_id->RemoteNumNeighbors);
    check_free(neighborhood_id->RemoteCommIndex );

    shan_segment_t *shared_segment = &(neighborhood_id->shared_segment);
    shan_free_shared(shared_segment);

    SUCCESS_OR_DIE (gaspi_wait (0, GASPI_BLOCK));
    MPI_Barrier(neighborhood_id->MPI_COMM_ALL);

    SUCCESS_OR_DIE(gaspi_segment_delete(neighborhood_id->neighbor_hood_id));
    shan_remote_t *remote_segment  = &(neighborhood_id->remote_segment);
    shan_free_remote(remote_segment);

    return SHAN_SUCCESS;

}

int shan_comm_init_comm(shan_neighborhood_t *const neighborhood_id
			, int neighbor_hood_id
			, int *neighbors
			, int num_neighbors
			, long *maxSendSz
			, long *maxRecvSz
			, int *max_nelem_send
			, int *max_nelem_recv
			, int num_type 
			, MPI_Comm MPI_COMM_SHM
			, MPI_Comm MPI_COMM_ALL
			)
{
    int i, j, k;
    ASSERT(neighborhood_id != NULL);
    ASSERT(neighbors != NULL);
    ASSERT(num_neighbors > 0);

    ASSERT(maxSendSz != NULL);
    ASSERT(maxRecvSz != NULL);
    ASSERT(max_nelem_send != NULL);
    ASSERT(max_nelem_recv != NULL);
    ASSERT(num_type > 0);

    ASSERT(MPI_COMM_SHM != MPI_COMM_NULL);
    ASSERT(MPI_COMM_ALL != MPI_COMM_NULL);

    neighborhood_id->neighbor_hood_id = neighbor_hood_id;

    neighborhood_id->MPI_COMM_SHM = MPI_COMM_SHM;
    neighborhood_id->MPI_COMM_ALL = MPI_COMM_ALL;
  
    MPI_Comm_rank(MPI_COMM_SHM, &(neighborhood_id->iProcLocal));
    MPI_Comm_size(MPI_COMM_SHM, &(neighborhood_id->nProcLocal));

    MPI_Comm_rank(MPI_COMM_ALL, &(neighborhood_id->iProcGlobal));
    MPI_Comm_size(MPI_COMM_ALL, &(neighborhood_id->nProcGlobal));

    neighborhood_id->master = neighborhood_id->iProcGlobal;
    MPI_Bcast(&(neighborhood_id->master)
	      , 1
	      , MPI_INT
	      , 0
	      , neighborhood_id->MPI_COMM_SHM
	);
  
    neighborhood_id->remote_master
	= check_malloc (neighborhood_id->nProcGlobal * sizeof(int));
  
    MPI_Allgather(&(neighborhood_id->master)
		  , 1
		  , MPI_INT
		  , neighborhood_id->remote_master
		  , 1
		  , MPI_INT
		  , neighborhood_id->MPI_COMM_ALL
	);  
    ASSERT(neighborhood_id->remote_master[neighborhood_id->iProcGlobal]
	   == neighborhood_id->master);
  
    neighborhood_id->num_neighbors = num_neighbors;
    neighborhood_id->neighbors = check_malloc(num_neighbors * sizeof(int));
    neighborhood_id->local_stage_count = check_malloc(num_neighbors * sizeof(int));
  
    for (i = 0; i < num_neighbors; ++i)
    {
	ASSERT(neighbors[i] >= 0);
	neighborhood_id->neighbors[i] = neighbors[i];
	neighborhood_id->local_stage_count[i] = 0;
    }

    neighborhood_id->num_local  = 0;
    for (i = 0; i < neighborhood_id->num_neighbors; ++i)
    {
	int const rank = neighborhood_id->neighbors[i];
	if(shan_comm_local_rank(neighborhood_id
				, rank
	       ) == -1)
	{
	    neighborhood_id->num_local++;
	}
    }  

    neighborhood_id->num_type = num_type;


    /*
     * negotiate remote comm index
     */
    shan_negotiate_meta_data(neighborhood_id);

  
    /* 
     * global max comm sizes 
     */
    MPI_Allreduce( MPI_IN_PLACE
		   , maxSendSz
		   , num_type
		   , MPI_LONG
		   , MPI_MAX
		   , neighborhood_id->MPI_COMM_ALL
	);

    MPI_Allreduce( MPI_IN_PLACE
		   , maxRecvSz
		   , num_type
		   , MPI_LONG
		   , MPI_MAX
		   , neighborhood_id->MPI_COMM_ALL
	);
  
    /* 
     * global max elements
     */
    MPI_Allreduce( MPI_IN_PLACE
		   , max_nelem_send
		   , num_type
		   , MPI_INT
		   , MPI_MAX
		   , neighborhood_id->MPI_COMM_ALL
	);

    MPI_Allreduce( MPI_IN_PLACE
		   , max_nelem_recv
		   , num_type
		   , MPI_INT
		   , MPI_MAX
		   , neighborhood_id->MPI_COMM_ALL
	);

    int const page_size = sysconf (_SC_PAGESIZE);
    neighborhood_id->type_element
	= check_malloc(num_type * sizeof(shan_element_t));    
  

    int const max_header_len = NELEM_COMM_HEADER * sizeof(int);
    long remoteSz = 0;
    for (i = 0; i < num_type; ++i)
    {
	long sz = maxSendSz[i] + max_header_len;
	sz = UP(sz, ALIGNMENT);
	neighborhood_id->type_element[i].maxSendSz    = sz;
	neighborhood_id->type_element[i].max_nelem_send = max_nelem_send[i];
	neighborhood_id->type_element[i].SendOffset[0] = remoteSz;
	neighborhood_id->type_element[i].SendOffset[1] = remoteSz + sz;
	remoteSz += 2 * sz;
    }
    for (i = 0; i < num_type; ++i)
    {
	long sz = maxRecvSz[i] + max_header_len;
	sz = UP(sz, ALIGNMENT);
	neighborhood_id->type_element[i].maxRecvSz    = sz;
	neighborhood_id->type_element[i].max_nelem_recv = max_nelem_recv[i];
	neighborhood_id->type_element[i].RecvOffset[0] = remoteSz;
	neighborhood_id->type_element[i].RecvOffset[1] = remoteSz + sz;
	remoteSz += 2 * sz;
    }	  
    neighborhood_id->remoteSz = UP(num_neighbors * remoteSz, page_size);
  
    long elemOffset = 0;
    for (i = 0; i < num_type; ++i)
    {
	neighborhood_id->type_element[i].elemOffset = elemOffset;
	elemOffset += MAX_SHARED_NOTIFICATION * sizeof(shan_notification_t)
	    + 4 * sizeof(int) + (max_nelem_send[i] + max_nelem_recv[i]) * sizeof(long);
    }
    long const typeOffset = UP(num_neighbors * elemOffset, page_size);

  
    for (i = 0; i < num_type; ++i)
    {
	neighborhood_id->type_element[i].local_recv_count
	    = check_malloc(num_neighbors *sizeof(int));
	neighborhood_id->type_element[i].local_send_count
	    = check_malloc(num_neighbors *sizeof(int));
	neighborhood_id->type_element[i].local_ack_count
	    = check_malloc(num_neighbors *sizeof(int));

      
	for (j = 0; j < num_neighbors; ++j)
	{
	    neighborhood_id->type_element[i].local_recv_count[j]  = 0;
	    neighborhood_id->type_element[i].local_send_count[j]  = 0;
	    neighborhood_id->type_element[i].local_ack_count[j]   = 0;
	}
    }
  

    /*
     * allocate and bind remote comm segment 
     */
    shan_comm_alloc_comm(neighborhood_id);
  
    shan_segment_t *shared_segment = &(neighborhood_id->shared_segment);
    int res = shan_alloc_shared(shared_segment
				, neighborhood_id->neighbor_hood_id
				, SHAN_TYPE
				, typeOffset
				, neighborhood_id->MPI_COMM_SHM
	);


    ASSERT(res == SHAN_SUCCESS);


    for (i = 0; i < num_type; ++i)
    {
	type_local_t type_info;
	shan_get_shared_type(&type_info
			     , neighborhood_id
			     , neighborhood_id->iProcLocal
			     , neighborhood_id->num_neighbors
			     , i
	    );

	for (k = 0; k < num_neighbors * MAX_SHARED_NOTIFICATION; ++k)
	{
	    shan_notify_init_shared(type_info.nid
				    , k
		);
	}
	  
	for (k = 0; k < num_neighbors; ++k)
	{
	    type_info.nelem_send[k] = 0;
	    type_info.nelem_recv[k] = 0;
	    type_info.send_sz[k]    = 0;
	    type_info.recv_sz[k]    = 0;
	}
      
	for (k = 0; k < num_neighbors * max_nelem_send[i]; ++k)
	{
	    type_info.send_offset[k]    = 0;
	}
      
	for (k = 0; k < num_neighbors * max_nelem_recv[i]; ++k)
	{
	    type_info.recv_offset[k]    = 0;
	}

    }
  
    return SHAN_SUCCESS;
}


int shan_comm_waitsome_local(shan_neighborhood_t *const neighborhood_id
			     , int const type_id
			     , int const idx
    )
{
    int iProcRemote, id = -1;  
    int const rank = neighborhood_id->neighbors[idx];
    ASSERT((iProcRemote = shan_comm_local_rank(neighborhood_id
					       , rank
		)) != -1);
    int const RemoteCommIdx = neighborhood_id->RemoteCommIndex[idx];
    int const RemoteNumNeighbors = neighborhood_id->RemoteNumNeighbors[idx];

    int rval = -1;
    shan_test_shared(neighborhood_id
		     , iProcRemote
		     , RemoteNumNeighbors
		     , type_id
		     , RemoteCommIdx
		     , &rval
	);

    int recv_count 
	= neighborhood_id->type_element[type_id].local_recv_count[idx];

    if (rval > recv_count)
    {
	ASSERT(rval == recv_count + 1);
	id = idx;
    }
  
    return (id == -1) ? -1 : SHAN_SUCCESS;
}


int shan_comm_waitsome_remote(shan_neighborhood_t *const neighborhood_id
			      , int const type_id
			      , int const idx
    )
{
    int const iProcLocal    = neighborhood_id->iProcLocal;
    int const num_neighbors = neighborhood_id->num_neighbors;
    int const num_type      = neighborhood_id->num_type;

    int id = -1;      
    int const sid 
	= (neighborhood_id->type_element[type_id].local_recv_count[idx]) % 2;  
    int const nid = GET_NOTIFICATION_ID(sid, num_type, type_id, num_neighbors, idx);
  
    shan_remote_t *const remote_segment = &(neighborhood_id->remote_segment);
    gaspi_notification_id_t tmp_id;
    gaspi_notification_t nval;
    gaspi_return_t ret;
    if (( ret =
	  gaspi_notify_waitsome (remote_segment->shan_id
				 , nid
				 , 1
				 , &tmp_id
				 , GASPI_TEST
	      )
	    ) == GASPI_SUCCESS)
    {
	SUCCESS_OR_DIE(gaspi_notify_reset (remote_segment->shan_id
					   , tmp_id
					   , &nval
			   )); 
	int const remote_rank = nval - 1;

	type_local_t type_info;
	shan_get_shared_type(&type_info
			     , neighborhood_id
			     , iProcLocal
			     , neighborhood_id->num_neighbors
			     , type_id
	    );

	long comm_buffer_offset = num_neighbors * neighborhood_id->type_element[type_id].RecvOffset[sid]
            + idx * neighborhood_id->type_element[type_id].maxRecvSz;

	void *comm_ptr = (char*) remote_segment->shan_ptr + comm_buffer_offset;
	int *const comm_header = (int *) comm_ptr;
	int const nelem_send   = *(comm_header);
	int const send_sz      = *(comm_header + 1);
	int const rval         = *(comm_header + 2);	
	int const rank         = neighborhood_id->neighbors[idx];

	ASSERT(rank == remote_rank);
	ASSERT(rval > neighborhood_id->type_element[type_id].local_recv_count[idx]);
	ASSERT(neighborhood_id->type_element[type_id].local_recv_count[idx] <= rval + 2);

#ifdef USE_VARIABLE_MESSAGE_LEN
	type_info.nelem_recv[idx] = nelem_send;
	type_info.recv_sz[idx]    = send_sz;
#else
	ASSERT(type_info.nelem_recv[idx] == nelem_send);
	ASSERT(type_info.recv_sz[idx] == send_sz);
#endif
	id = idx;

    }
    else
    {
	ASSERT (ret != GASPI_ERROR);
    }

    return (id == -1) ? -1 : SHAN_SUCCESS;
}





int shan_comm_notify_or_write(shan_neighborhood_t *const neighborhood_id
			      , shan_segment_t *const data_segment
			      , int type_id
			      , int idx
    )

{
    int const iProcLocal    = neighborhood_id->iProcLocal;
    int const num_neighbors = neighborhood_id->num_neighbors;
    int const num_type      = neighborhood_id->num_type;
    int iProcRemote;

    ASSERT(idx >= 0);
    ASSERT(idx < num_neighbors);
  
    int const rank = neighborhood_id->neighbors[idx];
    int const sid  
	= (neighborhood_id->type_element[type_id].local_send_count[idx]) % 2;	    
	
    if ((iProcRemote = shan_comm_local_rank(neighborhood_id
					    , rank
	     )) != -1)
    {
	shan_increment_local(neighborhood_id
			     , type_id
			     , idx
	    );      

	++(neighborhood_id->type_element[type_id].local_send_count[idx]);
    }
    else
    {
	int i;      
	type_local_t type_info;
	shan_get_shared_type(&type_info
			     , neighborhood_id
			     , iProcLocal
			     , neighborhood_id->num_neighbors
			     , type_id
	    );
      
	int nelem_send     = type_info.nelem_send[idx];
	int send_sz        = type_info.send_sz[idx];
	long *send_offset  = type_info.send_offset
	    + idx * neighborhood_id->type_element[type_id].max_nelem_send;
      
	void *data_ptr;
	shan_get_shared_ptr(data_segment
			    , iProcLocal
			    , &data_ptr);

	int const RemoteCommIdx = neighborhood_id->RemoteCommIndex[idx];
	int const RemoteNumNeighbors = neighborhood_id->RemoteNumNeighbors[idx];

	const gaspi_offset_t offset_local  = 
	    num_neighbors * neighborhood_id->type_element[type_id].SendOffset[sid]
	    + idx * neighborhood_id->type_element[type_id].maxSendSz;
      
	const gaspi_offset_t offset_remote = 
	    RemoteNumNeighbors * neighborhood_id->type_element[type_id].RecvOffset[sid]
	    + RemoteCommIdx * neighborhood_id->type_element[type_id].maxRecvSz;

	const gaspi_notification_id_t nid
	    = GET_NOTIFICATION_ID(sid, num_type, type_id, RemoteNumNeighbors, RemoteCommIdx);

	gaspi_number_t notification_num;
	SUCCESS_OR_DIE(gaspi_notification_num (&notification_num));
	ASSERT(nid < (int) notification_num);
      
	shan_remote_t *const remote_segment = &(neighborhood_id->remote_segment);
	void *comm_ptr = (char*) remote_segment->shan_ptr + offset_local;
	
	if (nelem_send > 0)
	{
	    for (i = 0; i < nelem_send; ++i)
	    {	      
		void *restrict dest = (char*) comm_ptr 
		    + NELEM_COMM_HEADER * sizeof(int) + i * send_sz;
		void *restrict src  = (char*) data_ptr + send_offset[i];
		memcpy(dest, src, send_sz);
	    }
	}
	  
	int *const comm_header = (int *) ((char*) comm_ptr);
	*(comm_header)      = nelem_send;
	*(comm_header + 1)  = send_sz;
	*(comm_header + 2)  = neighborhood_id->type_element[type_id].local_send_count[idx] + 1;

	long const comm_size = nelem_send * send_sz + NELEM_COMM_HEADER * sizeof(int);		
	write_notify_and_wait ( remote_segment->shan_id
				, offset_local
				, rank
				, offset_remote
				, (gaspi_size_t) comm_size
				, (gaspi_notification_id_t) nid
				, (gaspi_notification_t) neighborhood_id->iProcGlobal + 1
				, 0
	    );

	++(neighborhood_id->type_element[type_id].local_send_count[idx]);
    }

    return SHAN_SUCCESS;

}

