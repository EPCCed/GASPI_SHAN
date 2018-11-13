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

#include "exchange_data_mpi.h"
#include "read_netcdf.h"
#include "comm_data.h"
#include "comm_util.h"
#include "solver_data.h"
#include "error_handling.h"
#include "util.h"
#include "rangelist.h"

#define DATAKEY 4711

static void init_communication_data(int iProc, int nProc, comm_data *cd)
{
    ASSERT(cd != NULL);

    /* init comm data */
    cd->nProc = nProc;
    cd->iProc = iProc;

    cd->ndomains = 0;
    cd->ncommdomains = 0;
    cd->nownpoints = 0;
    cd->naddpoints = 0;

    cd->commpartner = NULL;
    cd->sendcount = NULL;
    cd->recvcount = NULL;
    cd->addpoint_owner = NULL;
    cd->addpoint_id = NULL;
    cd->recvindex = NULL;
    cd->sendindex = NULL;

    cd->nreq = 0;
    cd->req = NULL;
    cd->stat = NULL;
    cd->recvbuf = NULL;
    cd->sendbuf = NULL;

    cd->send_counter = NULL;
    cd->recv_counter = NULL;

}

void read_communication_data(int ncid, comm_data *cd)
{
    ASSERT(cd != NULL);

    /* read val */
    cd->ndomains = get_nc_val(ncid,"ndomains");
    cd->nownpoints = get_nc_val(ncid,"nownpoints");
 
    /* scalar execution */
    if (cd->ndomains == 1)
    {
	return;
    }

    /* read val */
    cd->naddpoints = get_nc_val(ncid,"naddpoints");
    cd->ncommdomains = get_nc_val(ncid,"ncommdomains");

    /* sanity check*/
    ASSERT(cd->ndomains >= 1);
    ASSERT(cd->ndomains == cd->nProc);

    /* sanity check*/
    ASSERT(cd->naddpoints > 0);
    ASSERT(cd->ncommdomains > 0);

    /* alloc */
    cd->commpartner = check_malloc(cd->ncommdomains * sizeof(int));
    cd->sendcount = check_malloc(cd->ndomains * sizeof(int));
    cd->recvcount = check_malloc(cd->ndomains * sizeof(int));
    cd->addpoint_owner = check_malloc(cd->naddpoints * sizeof(int));
    cd->addpoint_id = check_malloc(cd->naddpoints * sizeof(int));

    /* read data */
    get_nc_int(ncid,"commpartner",cd->commpartner);
    get_nc_int(ncid,"sendcount",cd->sendcount);
    get_nc_int(ncid,"recvcount",cd->recvcount);
    get_nc_int(ncid,"addpoint_owner",cd->addpoint_owner);
    get_nc_int(ncid,"addpoint_idx",cd->addpoint_id);

}

static void create_recvsend_index(comm_data *cd)
{
    ASSERT(cd != NULL);

    /* threading model only */
    if (cd->ndomains == 1)
    {
	return;
    }

    int i;
    const int nown   = cd->nownpoints;
    const int nadd   = cd->naddpoints;
    const int nProc  = cd->nProc;
    const int iProc  = cd->iProc;

    ASSERT(cd != NULL);
    ASSERT(cd->ndomains >= 1);
    ASSERT(cd->ncommdomains != 0);
    ASSERT(cd->sendcount != NULL);
    ASSERT(cd->recvcount != NULL);
    ASSERT(cd->addpoint_id != NULL);
    ASSERT(cd->addpoint_owner != NULL);

    /* alloc index tables */
    cd->sendindex = check_malloc(nProc * sizeof(int *));
    cd->recvindex = check_malloc(nProc * sizeof(int *));

    for(i = 0; i < nProc; i++)
    {
	cd->sendindex[i] = NULL;
	cd->recvindex[i] = NULL;
    }

    for(i = 0; i < cd->ncommdomains; i++)
    {
	int k = cd->commpartner[i];
	int j;
	if(cd->sendcount[k] > 0)
	{
	    cd->sendindex[k] = check_malloc(cd->sendcount[k] * sizeof(int));
	    for(j = 0; j < cd->sendcount[k]; j++)
	    {
		cd->sendindex[k][j] = -1;
	    }
	}

	if(cd->recvcount[k] > 0)
	{
	    int count = 0;
	    cd->recvindex[k] = check_malloc(cd->recvcount[k] * sizeof(int));
	    for(j = 0; j < nadd; j++)
	    {
		if(cd->addpoint_owner[j] == k)
		{
		    cd->recvindex[k][count++] = nown + j;
		}
	    }
	}
    }

    int sz = 0;
    for(i = 0; i < cd->ncommdomains; i++)
    {
	int k = cd->commpartner[i];
	sz = MAX(sz, (int)cd->sendcount[k]);
	sz = MAX(sz, (int)cd->recvcount[k]);
    }

    int *ibuf = check_malloc(sz * sizeof(int));
    for(i = 0; i < cd->ncommdomains; i++)
    {      
	int k          = cd->commpartner[i]; 
	int recvcount  = cd->recvcount[k];
	int sendcount  = cd->sendcount[k];
	int *recvindex = cd->recvindex[k];
	int *sendindex = cd->sendindex[k];
	int j;
          
	if(k > iProc) /* first send */
	{
	    for(j = 0; j < recvcount; j++)
	    {	      
		int idx = recvindex[j] - nown;
		ibuf[j] = cd->addpoint_id[idx]; 
	    }

	    MPI_Send(ibuf
		     , recvcount * sizeof(int)
		     , MPI_BYTE
		     , k
		     , DATAKEY
		     , MPI_COMM_WORLD
		);
	    MPI_Recv(ibuf
		     , sendcount * sizeof(int)
		     , MPI_BYTE
		     , k
		     , DATAKEY
		     , MPI_COMM_WORLD
		     , MPI_STATUS_IGNORE
		);

	    for(j = 0; j < sendcount; j++)
	    {
		sendindex[j] = ibuf[j];
	    }
	}
	else  /* first receive */
	{
	    MPI_Recv(ibuf
		     , sendcount * sizeof(int)
		     , MPI_BYTE
		     , k
		     , DATAKEY
		     , MPI_COMM_WORLD
		     , MPI_STATUS_IGNORE
		);
	    for(j = 0; j < sendcount; j++)
	    {
		sendindex[j] = ibuf[j];
	    }
	    for(j = 0; j < recvcount; j++)
	    {
		int idx = recvindex[j] - nown;
		ibuf[j] = cd->addpoint_id[idx];
	    }
	    MPI_Send(ibuf
		     , recvcount * sizeof(int)
		     , MPI_BYTE
		     , k
		     , DATAKEY
		     , MPI_COMM_WORLD
		);
	}
    }

    check_free(ibuf);

}

void init_communication(int argc, char *argv[], comm_data *cd)
{
    ASSERT(cd != NULL);

    /* MPI init */
    int nProc, iProc;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iProc);

    init_communication_data(iProc, nProc, cd);

}


void compute_communication_tables(comm_data *cd)
{
    ASSERT(cd != NULL);

    /* threading model only */
    if (cd->ndomains == 1)
    {
	return;
    }

    ASSERT(cd != NULL);
    ASSERT(cd->naddpoints != 0);
    ASSERT(cd->addpoint_owner != NULL);
    ASSERT(cd->addpoint_id != NULL);
    ASSERT(cd->commpartner != NULL);
    ASSERT(cd->sendcount != NULL);
    ASSERT(cd->recvcount != NULL);

    create_recvsend_index(cd);

#ifdef DEBUG
    int i;
    for(i = 0; i < cd->ncommdomains; i++)
    {
	int k = cd->commpartner[i]; 
	printf(" rank %8d: send %8d to   %8d\n", 
	       cd->iProc, cd->sendcount[k], k);
	printf(" rank %8d: recv %8d from %8d\n", 
	       cd->iProc, cd->recvcount[k], k);
    }
#endif

    /* allocate requests, statuses and buffer */
    const int max_elem_sz = NGRAD * 3;
    init_mpi_requests(cd, max_elem_sz);

}


void free_ressources(comm_data *cd
		     , solver_data *sd
    )
{	
    int i;
    ASSERT(cd != NULL);

    /* solver */
    free_solver_data(sd);

    /* threading model only */
    if (cd->ndomains == 1)
    {
	return;
    }

    /* rangelist */
    free_rangelist_comm(sd);

    /* free */
    check_free(cd->addpoint_owner);
    check_free(cd->addpoint_id);

    for(i = 0; i < cd->ncommdomains; i++)
    {
	int k = cd->commpartner[i];
	if(cd->sendcount[k] > 0)
	{
	    check_free(cd->sendindex[k]);
	}
	if(cd->recvcount[k] > 0)
	{
	    check_free(cd->recvindex[k]);
	}
    }

    check_free(cd->sendindex);
    check_free(cd->recvindex);

    check_free(cd->commpartner);
    check_free(cd->sendcount);
    check_free(cd->recvcount);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

}
