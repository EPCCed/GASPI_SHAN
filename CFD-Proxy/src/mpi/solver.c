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
#include "solver.h"
#include "util.h"

#include "flux.h"
#include "gradients.h"
#include "rangelist.h"
#include "exchange_data_mpi.h"

#define N_MEDIAN 25
#define N_SOLVER 2

void test_solver(comm_data *cd, solver_data *sd)
{
    int k;
    double time, median[N_SOLVER][N_MEDIAN];

    for (k = 0; k < N_MEDIAN; ++k)
    { 
	int i;
	/* comm free */
	time = -now();
	MPI_Barrier(MPI_COMM_WORLD);
	for (i = 0; i < sd->niter; ++i)
	{
	    compute_gradients_gg_comm_free(cd, sd);
	    compute_psd_flux(sd);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	time += now();
	median[0][k] = time;

	/* all remaining tests require ndomains >=1 */
	if (cd->ndomains == 1)
	{
	    continue;
	}

#if 1
	/* MPI bulk sync */
	time = -now();
	MPI_Barrier(MPI_COMM_WORLD);
	for (i = 0; i < sd->niter; ++i)
	{
	    compute_gradients_gg_mpi(cd, sd);
	    compute_psd_flux(sd);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	time += now();
	median[1][k] = time;
#endif
      
	if (cd->iProc == 0)
	{
	    printf(".");
	    fflush(stdout);
	}
    }
  
    if (cd->iProc == 0)
    {

	printf("\n\n*** SETUP\n");
	printf("                                 nProc: %d\n",cd->nProc);
	printf("                                 NITER: %d\n",sd->niter);
	printf("                              N_MEDIAN: %d\n",N_MEDIAN);

	for (k = 0; k < N_SOLVER; ++k)
	{ 
	    sort_median(&median[k][0], &median[k][N_MEDIAN-1]);
	}

	printf("\n*** TIMINGS\n");
	printf("                             comm_free: %10.6f\n",median[0][N_MEDIAN/2]);
	printf("                      exchange_dbl_mpi: %10.6f\n",median[1][N_MEDIAN/2]);

    }
}
