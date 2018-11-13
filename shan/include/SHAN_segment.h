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


#ifndef SHAN_segment_H
#define SHAN_segment_H

#include <mpi.h>

#ifdef __cplusplus
extern "C"
{
#endif

/** \file SHAN_segment.h
    \brief SHAN_segment header for notifications in shared memory.

    The SHAN (SHA_red N_otifications) interface is a user-level 
    API which aims at migrating flat MPI (legacy) code towards 
    an asynchronous dataflow model.  SHAN uses the GASPI API 
    and extends ideas from MPI shared windows.

    GASPI is a PGAS communication library which is based on the 
    concept of one-sided, notified communication. The synchronization 
    context here is bundled together with a one-sided message such that 
    a communication target becomes able to test for completion of the 
    received one-sided communication.

    Traditionally the GASPI programming model has been aimed at 
    multithreaded or task-based applications. In GASPI the synchronization
    context is bundled together with a one-sided message such that
    a communication target becomes able to test for completion of the
    received one-sided communication.

    In order to support a migration of legacy applications 
    (with a flat MPI communication model)  towards GASPI, we have extended 
    the concept of shared MPI windows towards a notified communication model 
    in which the processes sharing a common window become able to see all 
    one-sided and notified communication targeted at this window.
    Similarly we have extended the concept of MPI shared windows with
    shared notifications, which are globally visible in shared memory.

    Besides the possibility to entirely avoid node-internal communication 
    and to make use of a much improved overlap of communication and 
    computation the model of notified communication in GASPI shared windows
    will allow legacy SPMD applications a transition towards 
    an asynchronous dataflow model.
*/

enum shan_type  {
     SHAN_DATA   = 0, 
     SHAN_TYPE   = 1,
 };
    
enum shan_return_val  {
    SHAN_ERROR   = -2, 
    SHAN_FAIL    = -1,
    SHAN_SUCCESS =  0
};
    
/** 64 byte aligned notifications struct for shared mem notifications.
 */
typedef struct
{
    volatile int val  __attribute__((aligned(64))); //!< notification value
} shan_notification_t;
    
/** Segment struct, shared, holds all segment information.
 */
typedef struct
{
#ifndef __cplusplus
    void **restrict ptr_array;  //!< shan ptr array
#else
    void ** ptr_array;          //!< shan ptr array
#endif
    int shan_id;                //!< shared segment id
    int shan_type;              //!< shared segment type
    long dataSz;                //!< segment size
    long *localDataSz;          //!< segment size array
    MPI_Comm MPI_COMM_SHM;      //!< MPI shared mem communicator
    
#ifdef USE_MPI_SHARED_WIN
    MPI_Win segment_win;        //!< window handle
#else
    int *fd;                    //!< shmem file descriptor array
    char shan_domain_name[80];  //!< unique shmem name
#endif
    
} shan_segment_t;


/** Local allocation of shared memory of size dataSz
 *     
 * Note: Memory will be page-aligned.
 *
 * @param segment      - segment handle
 * @param shan_id      - (unique) segment id
 * @param shan_type    - type of allocated memory 
 * @param dataSz       - required memory size per rank in byte
 * @param MPI_COMM_SHM - shared mem communicator
 *
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
int shan_alloc_shared(shan_segment_t *const segment
		      , const int shan_id
		      , const int shan_type
		      , const long dataSz
		      , const MPI_Comm MPI_COMM_SHM 
	);

/** Free shared memory.
 *     
 * @param segment      - segment handle
 *
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
int shan_free_shared(shan_segment_t * const segment);


/** Gets shared mem pointer for node local ranks
 *     
 * @param segment      - segment handle
 * @param rank         - node local rank id
 * @param shm_ptr      - required memory size per rank in byte
 *
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
int shan_get_shared_ptr(shan_segment_t * const segment
			, const int rank
			, void **shm_ptr			  
    );  
    
    
/** Resets shared mem notification.
 *
 * @param ptr - pointer to shared notification array
 * @param idx - shared mem notification id
 * @param val - old value of the notification
 *
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
 int shan_notify_reset_shared (shan_notification_t *const ptr
			       , const int idx
			       , int * const val
     );
    
/** Increments shared mem notfication.
 *  Sets write fence such that local result is valid, 
 *  once the incremented value is visible for other local ranks.
 *
 * @param ptr - pointer to shared notification array
 * @param idx - shared mem notification id
 * @param increment - increment value
 * 
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
 int shan_notify_increment_shared (shan_notification_t *const ptr
				   , const int idx
				   , const int increment
     );
    
/** Tests for shared mem notfication.
 *
 * @param ptr - pointer to shared notification array
 * @param idx - shared mem notification id
 * @param val - current notification value
 * 
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
 int shan_notify_test_shared(shan_notification_t *const ptr
			     , const int idx
			     , int * const val
     );
    

/** Tests for shared mem notfication.
 *
 * @param ptr - pointer to shared notification array
 * @param idx - shared mem notification id
 * 
 * @return SHAN_SUCCESS in case of success, SHAN_ERROR in case of error.
 */
 int shan_notify_init_shared(shan_notification_t *const ptr
			     , const int idx
     );
    
#ifdef __cplusplus
}
#endif

#endif
