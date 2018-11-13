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


#ifndef F_SHAN_H
#define F_SHAN_H

#include <mpi.h>
#include "SHAN_segment.h"


#ifdef __cplusplus
extern "C"
{
#endif

/** \file F_SHAN.h
    \brief Wrapper functions for the SHAN library, mostly targeted at fortran applications.
*/

/** wrapper function for shan_alloc_shared
 *     
 * Note: Memory will be page-aligned.
 *
 * @param segment_id   - segment handle (data)
 * @param dataSz       - segment size in byte
 * @param shm_ptr      - shared mem pointer for allocated memory
 *
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
void f_shan_alloc_shared(const int segment_id
			, const long dataSz
			, void **restrict shm_ptr
    );
    
/** wrapper function for shan_free_shared
 *     
 * @param segment_id    - segment handle (data)
 *
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
void f_shan_free_shared(const int segment_id);


/** wrapper function for shan_free_comm
 *     
 * @param neighbor_hood_id - general neighborhood handle
 *
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
void f_shan_free_comm(const int neighbor_hood_id);


/** wrapper function for shan_init_comm
 *     
 * @param neighbor_hood_id - general neighborhood handle
 * @param neighbors       - comm partners (neighbors)
 * @param num_neighbors   - num comm partners (neighbors)
 * @param maxSendSz       - max send size for every comm type (byte)
 * @param maxRecvSz       - max recv size for every comm type (byte)
 * @param max_nelem_send  - max number of send elements for every comm type
 * @param max_nelem_recv  - max number of recv elements for every comm type
 * @param num_type        - number of types
 *
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
void f_shan_init_comm(const int neighbor_hood_id
		      , void *neighbors
		      , int num_neighbors
		      , void *maxSendSz
		      , void *maxRecvSz
		      , void *max_nelem_send
		      , void *max_nelem_recv
		      , int num_type
    );
    
/** wrapper function for shan_type_offset
 *     
 * @param neighbor_hood_id - general neighborhood handle
 * @param type_id        - used type segment
 * @param nelem_send     - ptr for current number of send elements. (in shared mem, visible node locally)
 * @param nelem_recv     - ptr for current number of recv elements. (in shared mem, visible node locally)
 * @param send_sz        - ptr for current send_size. (in shared mem, visible node locally)
 * @param recv_sz        - ptr for current recv size. (in shared mem, visible node locally)
 * @param send_idx       - ptr for current send offset list. (in shared mem, visible node locally)
 * @param recv_idx       - ptr for current recv offset list. (in shared mem, visible node locally)
 *
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
void f_shan_type_offset(const int neighbor_hood_id
			, const int type_id
			, void **nelem_send
			, void **nelem_recv
			, void **send_sz
			, void **recv_sz
			, void **send_idx
			, void **recv_idx
    );


/** wrapper function for shan_comm_wait4All
 *     
 * @param neighbor_hood_id - general neighborhood handle
 * @param segment_id     - (data) segment handle
 * @param type_id        - used type id
 *
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
void f_shan_comm_wait4All(const int neighbor_hood_id  
			  , const int segment_id
			  , const int type_id
    );


/** wrapper function for shan_comm_wait4AllSend
 *     
 * @param neighbor_hood_id - general neighborhood handle
 * @param type_id        - used type id
 *
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
void f_shan_comm_wait4AllSend(const int neighbor_hood_id  
			      , const int type_id
    );


/** wrapper function for shan_comm_wait4AllRecv
 *     
 * @param neighbor_hood_id - general neighborhood handle
 * @param segment_id     - (data) segment handle
 * @param type_id        - used type id
 *
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
void f_shan_comm_wait4AllRecv(const int neighbor_hood_id  
			      , const int segment_id
			      , const int type_id
    );


/** wrapper function for shan_comm_notify_or_write
 *  
 * @param neighbor_hood_id - general neighborhood handle
 * @param segment_id     - data segment handle
 * @param type_id        - used type id
 * @param idx            - comm index for target rank in neighborhood
 *
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
void f_shan_comm_notify_or_write(const int neighbor_hood_id
				 , const int segment_id
				 , const int type_id
				 , int idx
    );
    















    
#ifdef __cplusplus
}
#endif

#endif
