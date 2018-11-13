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
#include <xmmintrin.h>

#include <GASPI.h>

#include "SHAN_segment.h"
#include "SHAN_comm.h"
#include "SHAN_type.h"

#include "shan_exchange.h"
#include "gaspi_util.h"
#include "shan_core.h"
#include "shan_util.h"
#include "assert.h"


void shan_comm_shmemBarrier(shan_neighborhood_t *const neighborhood_id)
{
    MPI_Barrier(neighborhood_id->MPI_COMM_SHM);
}


int shan_comm_test4Send(shan_neighborhood_t *const neighborhood_id
			, int type_id
			, int idx
			) 
{  
    int const rank = neighborhood_id->neighbors[idx];
    volatile int ack_count 
	= neighborhood_id->type_element[type_id].local_ack_count[idx];  

    int iProcRemote, id = -1;
    if ((iProcRemote = shan_comm_local_rank(neighborhood_id
					    , rank
	     )) != -1)
    
    {
	int rval = -1;
	int const RemoteCommIdx = neighborhood_id->RemoteCommIndex[idx];
	int const RemoteNumNeighbors = neighborhood_id->RemoteNumNeighbors[idx];
	shan_test_shared(neighborhood_id
			 , iProcRemote
			 , RemoteNumNeighbors
			 , type_id
			 , RemoteNumNeighbors + RemoteCommIdx
			 , &rval
	    );		  
      
	if (rval > ack_count)
	{
	    ASSERT(ack_count + 1 == rval);
	    ++(neighborhood_id->type_element[type_id].local_ack_count[idx]);
	    id = idx;
	}
    }
    else
    {
	if (neighborhood_id->type_element[type_id].local_recv_count[idx] > 
	    neighborhood_id->type_element[type_id].local_ack_count[idx])
	{	
	    ++(neighborhood_id->type_element[type_id].local_ack_count[idx]);
	    id = idx;
	}
    }
    
    return (id == -1) ? -1 : SHAN_SUCCESS;
}



int shan_comm_wait4Send(shan_neighborhood_t *const neighborhood_id
			, int type_id
			, int idx
			) 
{  
    int res = -1;
    while ((res = shan_comm_test4Send(neighborhood_id
				      , type_id
				      , idx
		)) == -1)
    {	      
	_mm_pause();
    }	      

    return SHAN_SUCCESS;
}



int shan_comm_test4Recv(shan_neighborhood_t *const neighborhood_id
			, shan_segment_t *data_segment
			, int type_id
			, int idx
			) 
{
    int id = -1;
    int const rank = neighborhood_id->neighbors[idx];
    if (shan_comm_local_rank(neighborhood_id
			     , rank
	    ) != -1)
    {
	int res;
	if ((res = shan_comm_waitsome_local(neighborhood_id
					    , type_id
					    , idx
		 )) != -1)
	{
	    shan_comm_get_local(neighborhood_id
				, data_segment
				, type_id
				, idx
		); 
	    id = idx;
	}
    }
    else
    {
	int res;
	if ((res = shan_comm_waitsome_remote(neighborhood_id
					     , type_id
					     , idx
		 )) != -1)
	{
	    shan_comm_get_remote(neighborhood_id
				 , data_segment
				 , type_id
				 , idx
		); 
	    id = idx;
	}
    }
  
    return (id == -1) ? -1 : SHAN_SUCCESS;
}



int shan_comm_wait4Recv(shan_neighborhood_t *const neighborhood_id
			, shan_segment_t *data_segment
			, int type_id
			, int idx
			) 
{  
    int res = -1;
    while ((res = shan_comm_test4Recv(neighborhood_id
				      , data_segment
				      , type_id
				      , idx
		)) == -1)
    {	      
	_mm_pause();
    }	      

    return SHAN_SUCCESS;
}




int shan_comm_wait4AllSend(shan_neighborhood_t *const neighborhood_id
			   , int type_id
    ) 
{
    int i, num_buff = 0;    
    int const num_neighbors  = neighborhood_id->num_neighbors;
    for (i = 0; i < num_neighbors; ++i)
    {
	neighborhood_id->local_stage_count[i] = 0;
    }

    while (num_buff < num_neighbors)
    {
	for (i = 0; i < num_neighbors; ++i)
	{
	    if (!neighborhood_id->local_stage_count[i])
	    {
		int res;
		if ((res = shan_comm_test4Send(neighborhood_id
					       , type_id
					       , i
			 )) != -1)
		{	
		    neighborhood_id->local_stage_count[i] = 1;
		    num_buff++;
		}
	    }
	}
    }
  
    return SHAN_SUCCESS;
}



int shan_comm_wait4AllRecv(shan_neighborhood_t *const neighborhood_id  
			   , shan_segment_t *data_segment
			   , int type_id
    )
{
    int i, num_recv = 0;    
    int const num_neighbors  = neighborhood_id->num_neighbors;
    for (i = 0; i < num_neighbors; ++i)
    {
	neighborhood_id->local_stage_count[i] = 0;
    }

    while (num_recv < num_neighbors)
    {
	for (i = 0; i < num_neighbors; ++i)
	{
	    if (!neighborhood_id->local_stage_count[i])
	    {
		int res;
		if ((res = shan_comm_test4Recv(neighborhood_id
					       , data_segment
					       , type_id
					       , i
			 )) != -1)
		{
		    neighborhood_id->local_stage_count[i] = 1;
		    num_recv++;
		}
	    }
	}
    }
  
    return SHAN_SUCCESS;

}


int shan_comm_wait4All(shan_neighborhood_t *const neighborhood_id  
		       , shan_segment_t *data_segment
		       , int type_id
		       )
{
    int res = shan_comm_wait4AllRecv(neighborhood_id
				     , data_segment 
				     , type_id
	);
    ASSERT(res == SHAN_SUCCESS);


    res = shan_comm_wait4AllSend(neighborhood_id
				 , type_id
	);
    ASSERT(res == SHAN_SUCCESS);  


    return SHAN_SUCCESS;
}

