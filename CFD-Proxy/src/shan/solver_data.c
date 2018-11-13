/*
 * This file is part of a small exa2ct benchmark kernel
 * The kernel aims at a dataflow implementation for 
 * hybrid solvers which make use of unstructured meshes.
 *
 * Contact point for exa2ct: 
 *                 https://projects.imec.be/exa2ct
 *
 * Contact point for this kernel: 
 *                 christian.simmendinger@t-systems.com
 */


#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

#include "read_netcdf.h"
#include "solver_data.h"
#include "error_handling.h"
#include "rangelist.h"
#include "util.h"

#ifdef USE_SHAN
static void initSharedData(shan_segment_t *dataSegment
			   , long dataSz
			   , void **shm_ptr
    )
{
  int iProcGlobal, nProcGlobal;
  MPI_Comm_rank(MPI_COMM_WORLD, &iProcGlobal);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcGlobal);

  MPI_Comm MPI_COMM_SHM;
  MPI_Comm_split_type (MPI_COMM_WORLD
		       , MPI_COMM_TYPE_SHARED
		       , 0
		       , MPI_INFO_NULL
		       , &MPI_COMM_SHM
		       );
  int iProcLocal;
  MPI_Comm_rank(MPI_COMM_SHM, &iProcLocal);

  int const segment_id = 0;
  int res = shan_alloc_shared(dataSegment
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
#endif

static void init_var(double (*var)[NGRAD], int nallpoints)
{
  int i,j;
  for (i = 0; i< nallpoints; ++i)
    {
      for (j = 0; j < NGRAD; ++j)
	{
	  var[i][j] = 1.0;
	}
    }
}

static void init_grad(double (*grad)[NGRAD][3], int nallpoints)
{
  int i,j,k;
  for (i = 0; i< nallpoints; ++i)
    {
      for (j = 0; j < NGRAD; ++j)
	{
	  for (k = 0; k < 3; ++k)
	    {
	      grad[i][j][k] = 1.0;
	    }
	}
    }
}

static void init_flux(double (*psd_flux)[NFLUX], int nallpoints)
{
  int i,j;
  for (i = 0; i< nallpoints; ++i)
    {
      for (j = 0; j < NFLUX; ++j)
        {         
          psd_flux[i][j] = 1.0;
        }
    }
}

void init_solver_data(solver_data *sd, int NITER)
{
  ASSERT(sd != NULL);
  ASSERT(sd->nallpoints != 0);

  /* initialize var/grad */
  init_var(sd->var, sd->nallpoints);
  init_grad(sd->grad, sd->nallpoints);
  init_flux(sd->psd_flux, sd->nallpoints);

  /* set num iterations */
  sd->niter = NITER;
}


void read_solver_data(int ncid
		      , solver_data *sd
    )
{
  ASSERT(sd != NULL);

  sd->nfaces = 0;
  sd->nallfaces = 0;
  sd->nownpoints = 0;
  sd->nallpoints = 0;
  sd->ncolors = 0;
  sd->fpoint = NULL;
  sd->fnormal = NULL;
  sd->pvolume = NULL;
  sd->var = NULL;
  sd->grad = NULL;
  sd->fcolor = NULL;
  sd->niter = 0;

  /* read val */
  sd->ncolors = get_nc_val(ncid,"ncolors");
  sd->nfaces = get_nc_val(ncid,"nfaces");
  sd->nownpoints = get_nc_val(ncid,"nownpoints");
  sd->nallpoints = get_nc_val(ncid,"nallpoints");

  /* sanity check*/
  ASSERT(sd->ncolors > 0);
  ASSERT(sd->nfaces > 0);
  ASSERT(sd->nownpoints > 0);
  ASSERT(sd->nallpoints > 0);
  ASSERT(NGRAD > 0);

  /* alloc */
  sd->fpoint = check_malloc(sd->nfaces * 2 * sizeof(int));
  sd->fnormal = check_malloc(sd->nfaces * 3 * sizeof(double));
  sd->pvolume = check_malloc(sd->nallpoints * sizeof(double));
  sd->var = check_malloc(sd->nallpoints * NGRAD * sizeof(double));
  sd->psd_flux = check_malloc(sd->nallpoints * NFLUX * sizeof(double));

#ifndef USE_SHAN
  sd->grad = check_malloc(sd->nallpoints * NGRAD * 3 * sizeof(double));
#else
  void *shm_ptr = NULL;
  long dataSz = sd->nallpoints * NGRAD * 3 * sizeof(double);
  initSharedData(&(sd->dataSegment)
		 , dataSz
		 , &shm_ptr
      );
  sd->grad = (double (*)[NGRAD][3]) shm_ptr;
#endif

  /* read data */
  get_nc_int(ncid,"fpoint",&(sd->fpoint[0][0]));
  get_nc_double(ncid,"fnormal",&(sd->fnormal[0][0]));
  get_nc_double(ncid,"pvolume",sd->pvolume);

}


void free_solver_data(solver_data *sd)
{
    check_free(sd->fpoint);
    check_free(sd->fnormal);
    check_free(sd->pvolume);
    check_free(sd->var);
#ifndef USE_SHAN
    check_free(sd->grad);
#else
    int res = shan_free_shared(&(sd->dataSegment));
    ASSERT(res == SHAN_SUCCESS);
#endif

    check_free(sd->psd_flux);
}





