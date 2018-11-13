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


#ifndef SHAN_TYPE_H
#define SHAN_TYPE_H

#include <stdio.h>
#include <stdlib.h>

#include "GASPI.h"
#include "SHAN_segment.h"


#ifdef __cplusplus
extern "C"
{
#endif

/** \file SHAN_type.h
 *  \brief SHAN_type header. Type conversion in shared memory.
 *   
 */

  
/** Converts shared mem send type in shared mem recv type.
 *  Finalizes receive in shared mem.
 * 
 * @param neighborhood_id - general neighborhood handle
 * @param data_segment    - used data segment
 * @param type_id         - used type id
 * @param idx             - rank index in neighborhood
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
void shan_comm_get_local(shan_neighborhood_t *neighborhood_id
			 , shan_segment_t *data_segment
			 , const int type_id
			 , const int idx
    );

/** Finalizes receive for remote comm
 *
 * @param neighborhood_id - general neighborhood handle
 * @param data_segment    - used data segment
 * @param type_id         - used type id
 * @param idx             - rank index in neighborhood
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
 void shan_comm_get_remote(shan_neighborhood_t *neighborhood_id
			   , shan_segment_t *const data_segment
			   , int const type_id
			   , int const idx
     );

/** Returns type data structure for node local ranks
 *  
 * @param type_info       - type data struct (SHAN_comm.h)
 * @param neighborhood_id - general neighborhood handle
 * @param local_rank      - node local rank 
 * @param num_neighbors   - number of neighbors
 * @param type_id         - used type id
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */  
 int shan_get_shared_type(type_local_t *type_info
			  , shan_neighborhood_t *neighborhood_id
			  , int local_rank
			  , int num_neighbors
			  , int type_id
     );
    
    
/** Gets type data structure for node local ranks
 *  
 * @param neighborhood_id - general neighborhood handle
 * @param type_id         - used type id
 * @param nelem_send      - pointer to number of send elements in shared mem
 * @param nelem_recv      - pointer to number of recv elements in shared mem
 * @param send_sz         - pointer to send size in shared mem
 * @param recv_sz         - pointer to recv size in shared mem
 * @param send_offset     - pointer to offset of send elements in shared mem
 * @param recv_offset     - pointer to offset of recv elements in shared mem
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
 int shan_comm_type_offset(shan_neighborhood_t *neighborhood_id
			   , int type_id
			   , int **nelem_send
			   , int **nelem_recv
			   , int **send_sz
			   , int **recv_sz
			   , long **send_offset
			   , long **recv_offset
     );


/** Getter function for type data
 *  
 * @param type_info       - type data struct (SHAN_comm.h)   
 * @param shm_ptr         - pointer to shared memory 
 * @param num_neighbors   - rank local number of neighbors in neighborhood
 * @param type_element    - type element
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
 int shan_comm_get_type(type_local_t *type_info
			, void *shm_ptr
			, int num_neighbors
			, shan_element_t *type_element
     );


#ifdef __cplusplus
}
#endif

#endif


