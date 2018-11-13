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


#ifndef SHAN_COMM_H
#define SHAN_COMM_H


#include <stdio.h>
#include <stdlib.h>

#include "GASPI.h"
#include "SHAN_segment.h"


#ifdef __cplusplus
extern "C"
{
#endif

/** \file SHAN_comm.h
 *  \brief SHAN_comm header for persistant communication in shared memory.
 *   
 */

#define MAX_SHARED_NOTIFICATION 2 //!< 'have written' and 'have read' synchronization


/** Type struct, visible in shared memory
 */
typedef struct
{
    shan_notification_t *nid;    //!< synchronization for types 
    int *nelem_send;             //!< current num send elements per neighbor
    int *nelem_recv;             //!< current num recv elements per neighbor
    int *send_sz;                //!< current send size (in char) per neighbor
    int *recv_sz;                //!< current recv size (in char) per neighbor  
    long *send_offset;           //!< list of send offsets per neighbor
    long *recv_offset;           //!< list of recv offsets per neighbor
} type_local_t;


/** Segment struct, rank_local, holds all segment information.
 */
typedef struct
{
    int shan_id;                //!< shared segment id
    long dataSz;                //!< segment size array
    void *shan_ptr;    //!< local segment pointer    
} shan_remote_t;


/** comm struct, holds all communication information per type.
 */
typedef struct
{
    long maxSendSz;              //!< max send size per type (byte)
    long maxRecvSz;              //!< max recv size per type (byte)
    int  max_nelem_send;         //!< max recv size per type (byte)
    int  max_nelem_recv;         //!< max recv size per type (byte)
    long SendOffset[2];          //!< local offset for send per type (byte)
    long RecvOffset[2];          //!< local offset for recv per type (byte)
    long elemOffset;             //!< element offset in shared mem

    int *local_send_count;      //!< send stage counter array, per type
    int *local_recv_count;      //!< recv stage counter array, per type
    int *local_ack_count;       //!< acknowledge stage counter array, per type
    
} shan_element_t;


/** neighborhood comm struct, holds all communication 
 *  information for the neighborhood.
 */
typedef struct
{
    int neighbor_hood_id;       //!< neighborhood id
    MPI_Comm MPI_COMM_SHM;      //!< shared MPI communicator
    MPI_Comm MPI_COMM_ALL;      //!< global MPI communicator

    int num_neighbors;          //!< num comm partners (neighbors)
    int num_local;              //!< node local number of comm partners
    int *neighbors;             //!< list of neighbors, per rank
    int *RemoteCommIndex;       //!< the remote index corresponding to own rank
    int *RemoteNumNeighbors;    //!< remote number of neighbors for RemoteCommIndex

    int num_type;               //!< num types
    long *typeOffset;           //!< type offsets for all node local ranks

    shan_segment_t shared_segment; //!< shared window for local communication
    shan_element_t *type_element;  //!< local comm data for remote communication.

    long remoteSz;              //!< remote comm size, all types, send + recv (byte)
    shan_remote_t remote_segment;  //!< private segment for remote communication  
    
    int nProcLocal;             //!< num local ranks in shared mem
    int iProcLocal;             //!< local rank id
    int nProcGlobal;            //!< num global ranks
    int iProcGlobal;            //!< global rank id

    int master;                 //!< master of shared segment (local rank 0)
    int *remote_master;         //!< global list of masters

    int *local_stage_count;     //!< generic stage counter for wait4All(Send/Recv)
    
} shan_neighborhood_t;

/** Gets node local rank id.
 *  
 * @param neighborhood_id - handle for neighborhood
 * @param rank            - global rank
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
int shan_comm_local_rank(shan_neighborhood_t * const neighborhood_id
			 , const int rank
    );

/** Increments counter in shared mem
 *  
 * @param neighborhood_id - general neighborhood handle
 * @param type_id        - type index
 * @param idx            - comm index for target rank in neighborhood
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
int shan_increment_local(shan_neighborhood_t *const neighborhood_id				  
			 , int const type_id
			 , int const idx
    );

/** Initialize persistant communication for shared mem and GASPI.
 *  requires bidirectional communication for synchronization in
 *  one-sided communication.
 * 
 *  A zero length messages will work, no message at all will fail.
 * 
 *  - allocates shared and private mem for communication (double buffered).
 *  - figures out local and remote comm partners.
 *  - negotiates remote number f neighbors and comm index
 *
 * @param neighborhood_id - general neighborhood handle
 * @param neighbor_hood_id - neighborhood id
 * @param neighbors     - comm partners (neighbors)
 * @param num_neighbors - num comm partners (neighbors)
 * @param maxSendSz     - max send size for every type.
 * @param maxRecvSz     - max recv size for every type
 * @param max_nelem_send - max number of send elements per type
 * @param max_nelem_recv - max number of recv elements per type
 * @param num_type      - number of types
 * @param MPI_COMM_SHM - MPI shared mem communicator
 * @param MPI_COMM_ALL - embedding of shared communicator (typically MPI_COMM_WORLD) 
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
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
    );


/** Free communication ressources
 *
 * @param neighborhood_id - general neighborhood handle
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
int shan_comm_free_comm(shan_neighborhood_t *const neighborhood_id);


/** Writes data or flags data as readable.
 *  
 *  - aggregates send data into linear buffer or
 *  - flags data as readable 
 *     - number of elements 
 *     - element sizes and 
 *     - element offsets 
 *
 * @param neighborhood_id - general neighborhood handle
 * @param data_segment   - data segment handle
 * @param type_id        - type index
 * @param idx            - comm index for target rank in neighborhood
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
int shan_comm_notify_or_write(shan_neighborhood_t *const neighborhood_id
			      , shan_segment_t *data_segment
			      , int type_id
			      , int idx
    );


/** Waits for entire neighborhood
 *  
 *  - waits for either shared memory notifications
 *    or remote GASPI notifications 
 *  - directly converts send type into recv type in shared memory
 *  - unpacks pipelined remote communictation
 *    into the current receive type.
 *  - waits for all receive requests.
 *  - waits for all send requests 
 *
 * @param neighborhood_id - general neighborhood handle
 * @param data_segment   - data segment handle
 * @param type_id        - type index
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */

 int shan_comm_wait4All(shan_neighborhood_t *const neighborhood_id  
			, shan_segment_t *data_segment
			, int type_id
     );

/** Waits for entire neighborhood
 *  
 *  - waits for either shared memory notifications
 *    or remote GASPI notifications 
 *  - directly converts send type into recv type in shared memory
 *  - unpacks pipelined remote communictation
 *    into the current receive type.
 *  - waits for all receive requests.
 *
 * @param neighborhood_id - general neighborhood handle
 * @param data_segment   - data segment handle
 * @param type_id        - type index
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
int shan_comm_wait4AllRecv(shan_neighborhood_t *const neighborhood_id  
			   , shan_segment_t *data_segment
			   , int type_id
    );



/** waits for specific receive requests.
 *
 * @param neighborhood_id - general neighborhood handle
 * @param data_segment   - data segment handle
 * @param type_id        - type index
 * @param idx            - comm index for target rank in neighborhood
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
 int shan_comm_wait4Recv(shan_neighborhood_t *const neighborhood_id
			 , shan_segment_t *data_segment
			 , int type_id
			 , int idx
     );


/** Tests for specific receive requests.
 *
 * @param neighborhood_id - general neighborhood handle
 * @param data_segment   - data segment handle
 * @param type_id        - type index
 * @param idx            - comm index for target rank in neighborhood
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
int shan_comm_test4Recv(shan_neighborhood_t *const neighborhood_id
			, shan_segment_t *data_segment
			, int type_id
			, int idx
    );

/** Waits for entire neighborhood
 *  
 *  - waits for all send requests
 *
 * @param neighborhood_id - general neighborhood handle
 * @param type_id        - type index
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
 int shan_comm_wait4AllSend(shan_neighborhood_t *const neighborhood_id
			    , int type_id
     );
    
/** Waits for specific send requests
 *  
 * @param neighborhood_id - general neighborhood handle
 * @param type_id        - type index
 * @param idx            - comm index for target rank in neighborhood
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
int shan_comm_wait4Send(shan_neighborhood_t *const neighborhood_id
			, int type_id
			, int idx
    ); 

/** Tests for specific send requests
 *  
 * @param neighborhood_id - general neighborhood handle
 * @param type_id        - type index
 * @param idx            - comm index for target rank in neighborhood
 *
 * @return SHAN_COMM_SUCCESS in case of success, SHAN_COMM_ERROR in case of error.
 */
int shan_comm_test4Send(shan_neighborhood_t *const neighborhood_id
			, int type_id
			, int idx
    ); 


/** Shared mem barrier
 *  
 * @param neighborhood_id - general neighborhood handle
 *
 */
void shan_comm_shmemBarrier(shan_neighborhood_t *const neighborhood_id);

#ifdef __cplusplus
}
#endif

#endif


