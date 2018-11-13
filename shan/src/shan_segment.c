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
#include <math.h>
#include <string.h>
#include <GASPI.h>

#include "SHAN_segment.h"
#include "shan_util.h"
#include "assert.h"


#ifndef USE_MPI_SHARED_WIN
static void shan_alloc_shared_root(shan_segment_t *const segment
				   , const MPI_Comm MPI_COMM_SHM 
    )
{
    const int shan_id = segment->shan_id;
    const int shan_type = segment->shan_type;
    int i, j;

    int iProcLocal, nProcLocal;
    MPI_Comm_rank(MPI_COMM_SHM, &iProcLocal);
    MPI_Comm_size(MPI_COMM_SHM, &nProcLocal);

    segment->fd = check_malloc (nProcLocal * sizeof(int));
    for (i = 0; i < nProcLocal; ++i)
    {
	segment->fd[i]     = 0;
    }

    sprintf(segment->shan_domain_name, "/shan_shared_space_%d_%d_%d"
	    , shan_type
	    , shan_id
	    , 0
	);

    void *restrict ptr = NULL;
    int flag;
    for (i = 0; i < nProcLocal; ++i)
    {
	if (i == iProcLocal)
	{
	    segment->fd[0] = shm_open (segment->shan_domain_name
				       , O_CREAT | O_RDWR, 0666);
	    ASSERT(segment->fd[0] != -1);

	    if (iProcLocal == 0)
	    {
		if (ftruncate (segment->fd[0], segment->dataSz))
		{
		    close (segment->fd[0]);
		    exit (EXIT_FAILURE);
		}
	    } 

#if defined(MAP_HUGE_2MB) && defined(MAP_HUGETLB)
	    flag = MAP_SHARED | MAP_HUGETLB | MAP_HUGE_2MB;
	    ptr = (void *restrict) mmap (0, segment->dataSz, PROT_READ | PROT_WRITE
					 , flag, segment->fd[0], 0);	  
	    if( ptr == MAP_FAILED )
	    {
#endif
		flag = MAP_SHARED;
		ptr = (void *restrict) mmap (0, segment->dataSz, PROT_READ | PROT_WRITE
					     , flag, segment->fd[0], 0);		
		if( ptr == MAP_FAILED )
		{
		    close (segment->fd[0]);
		    exit (EXIT_FAILURE);
		}
#if defined(MAP_HUGE_2MB) && defined(MAP_HUGETLB)
	    }
#endif
	    long ptr_offset = 0;
	    for (j = 0; j < nProcLocal; ++j)
	    {
		segment->ptr_array[j] = (char *) ptr + ptr_offset;	    
		ptr_offset += segment->localDataSz[j];
	    }
	}

	__sync_synchronize();
	MPI_Barrier(MPI_COMM_SHM);

    }
}

#endif

int shan_alloc_shared(shan_segment_t *const segment
                      , const int shan_id
                      , const int shan_type
                      , const long dataSz
                      , const MPI_Comm MPI_COMM_SHM 
    )
{
    int i;
    ASSERT(segment != NULL);
    ASSERT(dataSz >= 0);

    int iProcLocal, nProcLocal;
    MPI_Comm_rank(MPI_COMM_SHM, &iProcLocal);
    MPI_Comm_size(MPI_COMM_SHM, &nProcLocal);

    segment->shan_id      = shan_id;
    segment->shan_type    = shan_type;

    segment->ptr_array    = check_malloc (nProcLocal * sizeof(void*));
    segment->localDataSz  = check_malloc (nProcLocal * sizeof(long));
    segment->MPI_COMM_SHM = MPI_COMM_SHM;

    int const page_size = sysconf (_SC_PAGESIZE);
    long sz = UP(dataSz, page_size);
  
    MPI_Allgather(&sz
		  , 1
		  , MPI_LONG
		  , segment->localDataSz
		  , 1
		  , MPI_LONG
		  , MPI_COMM_SHM
	);         

    segment->dataSz = 0;
    for (i = 0; i < nProcLocal; ++i)
    {
	segment->ptr_array[i] = NULL;
	segment->dataSz += segment->localDataSz[i];
    }
  
#ifdef USE_MPI_SHARED_WIN

    void *ptr = NULL;  
    MPI_Aint rsz;
    int dsp;  

    MPI_Info win_info;
    MPI_Info_create(&win_info);
    MPI_Info_set(win_info, "alloc_shared_noncontig", "true");
    MPI_Win_allocate_shared(sz
			    , page_size
			    , win_info
			    , segment->MPI_COMM_SHM
			    , &ptr
			    , &(segment->segment_win)
	);

    for (i = 0; i < nProcLocal; ++i)
    {      
	MPI_Win_shared_query (segment->segment_win 
			      , i
			      , &rsz
			      , &dsp
			      , &ptr
	    );
	segment->ptr_array[i] = ptr;
	ASSERT(rsz == segment->localDataSz[i]);
    }


#else


    shan_alloc_shared_root(segment
			     , MPI_COMM_SHM
	);

#endif

    __sync_synchronize();
    MPI_Barrier(MPI_COMM_SHM);

    ASSERT(segment->ptr_array[iProcLocal] != NULL);
    memset(segment->ptr_array[iProcLocal]
	   , 0
	   , sz
	);

    __sync_synchronize();
    MPI_Barrier(MPI_COMM_SHM);

    return SHAN_SUCCESS; 

}

int shan_get_shared_ptr(shan_segment_t * const segment
			, const int local
			, void **shm_ptr			  
			)
{
    ASSERT(segment->ptr_array != NULL);  
    *shm_ptr = segment->ptr_array[local];
  
    return SHAN_SUCCESS;
}

int shan_free_shared(shan_segment_t * const segment)
{
    MPI_Barrier(segment->MPI_COMM_SHM);

#ifdef USE_MPI_SHARED_WIN

    MPI_Win_free(&segment->segment_win);

#else
    int nProcLocal;
    MPI_Comm_size(segment->MPI_COMM_SHM, &nProcLocal);

    int iProcLocal;
    MPI_Comm_rank(segment->MPI_COMM_SHM, &iProcLocal);


    if (munmap (segment->ptr_array[0], segment->dataSz) || close (segment->fd[0]))
    {
	exit (EXIT_FAILURE);
    }

    MPI_Barrier(segment->MPI_COMM_SHM);  
    if (iProcLocal == 0)
    {
	int res = shm_unlink (segment->shan_domain_name);
	ASSERT(res != -1);
    }

    check_free(segment->fd);

#endif

    check_free(segment->ptr_array);
    check_free(segment->localDataSz);

    MPI_Barrier(segment->MPI_COMM_SHM);  
  
    return SHAN_SUCCESS; 
}

int shan_notify_increment_shared(shan_notification_t *const ptr
				 , const int idx
				 , const int increment
				 )
{
  
  volatile shan_notification_t *nid = ptr + idx;

  __sync_synchronize();
  volatile int res = __sync_add_and_fetch(&(nid->val),increment);

  
  ASSERT(res > 0);
  return SHAN_SUCCESS; 
}

int shan_notify_reset_shared(shan_notification_t *const ptr
			     , const int idx
			     , int * const val)
{
  volatile shan_notification_t *nid = ptr + idx;
  volatile int res = __sync_val_compare_and_swap (&(nid->val), nid->val, 0);

  *val =  res;
  return SHAN_SUCCESS;
}

int shan_notify_init_shared(shan_notification_t *const ptr
			    , const int idx
			    )
{
  volatile shan_notification_t *nid = ptr + idx;
  nid->val = 0;

  return SHAN_SUCCESS;
}



int shan_notify_test_shared(shan_notification_t *const ptr
			    , const int idx
			    , int * const val
			    )
{
  volatile shan_notification_t *nid = ptr + idx;
  volatile int res = nid->val;
  
  *val =  res;
  __sync_synchronize();

  return SHAN_SUCCESS;
}





