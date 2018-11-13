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
#include "solver_data.h"
#include "comm_data.h"
#include "rangelist.h"
#include "util.h"
#include "error_handling.h"

#define DATAKEY 4712


void init_mpi_requests(comm_data *cd, int dim2)
{
  int i;
  const int max_elem_sz = NGRAD * 3;
  ASSERT(dim2 == max_elem_sz);  
  int szd = sizeof(double);

  ASSERT(cd->nreq == 0);
  ASSERT(cd->ncommdomains > 0);
  
  int szr = 2 * cd->ncommdomains * sizeof(MPI_Request);
  int szs = 2 * cd->ncommdomains * sizeof(MPI_Status);   
  cd->req   = (MPI_Request *)check_malloc(szr);
  cd->stat  = (MPI_Status  *)check_malloc(szs);
  cd->nreq  = cd->ncommdomains;

  /* status flag */
  cd->send_counter = check_malloc(cd->ncommdomains*sizeof(counter_t));
  cd->recv_counter = check_malloc(cd->ncommdomains*sizeof(counter_t));
  for(i = 0; i < cd->ncommdomains; i++)
    {
      cd->send_counter[i].global = 0;
      cd->recv_counter[i].global = 0;
    }

  /* sendbuffer, recvbuffer  */
  cd->sendbuf = check_malloc(cd->ncommdomains * sizeof(double*));
  cd->recvbuf = check_malloc(cd->ncommdomains * sizeof(double*));
  for(i = 0; i < cd->ncommdomains; i++)
    {
      cd->sendbuf[i] = NULL;
      cd->recvbuf[i] = NULL;
    }

  for(i = 0; i < cd->ncommdomains; i++)
    {
      int k = cd->commpartner[i];
      if (cd->sendcount[k] > 0)
	{
	  cd->sendbuf[i] = 
	    check_malloc(dim2 * cd->sendcount[k] * max_elem_sz * szd);
	}
      if (cd->recvcount[k] > 0)
	{
	  cd->recvbuf[i] = 
	    check_malloc(dim2 * cd->recvcount[k] * max_elem_sz * szd);
	}
    }  

}


static void exchange_dbl_copy_in(comm_data *cd
				 , double *sbuf
				 , double *data
				 , int dim2
				 , int i
				 )
{
  int j;
  int *sendcount    = cd->sendcount;
  int **sendindex   = cd->sendindex;
  int k = cd->commpartner[i];
  int count = sendcount[k];

  if(count > 0)
    {
      for(j = 0; j < count; j++)
	{
	  int n1 = dim2 * j;
	  int n2 = dim2 * sendindex[k][j];
	  memcpy(&sbuf[n1], &data[n2], dim2 * sizeof(double));
	}
    }

}


static void exchange_dbl_copy_out(comm_data *cd
				  , double *rbuf
				  , double *data
				  , int dim2
				  , int i
				  )
{
  int j;
  int *recvcount    = cd->recvcount;
  int **recvindex   = cd->recvindex;
  int k = cd->commpartner[i];
  int count = recvcount[k];

  if(count > 0)
    {
      for(j = 0; j < count; j++)
	{
	  int n1 = dim2 * j;
	  int n2 = dim2 * recvindex[k][j];
	  memcpy(&data[n2], &rbuf[n1], dim2 * sizeof(double));
	}
    }

}



static void exchange_dbl_mpi_send(comm_data *cd
				  , double *data
				  , int dim2
				  , int i
				  )
{
    int ncommdomains  = cd->ncommdomains;
    int *commpartner  = cd->commpartner;
    int *sendcount    = cd->sendcount;

    int size, szd = sizeof(double);

    /* send */
    int k = commpartner[i];
    int count = sendcount[k];

    ASSERT(data != NULL);

    if(count > 0)
    {
	double *sbuf = cd->sendbuf[i];
	size = count * szd * dim2;
	MPI_Isend(sbuf
		  , size
		  , MPI_BYTE
		  , k
		  , DATAKEY
		  , MPI_COMM_WORLD
		  , &(cd->req[ncommdomains + i])
	    );
    }
}


void exchange_dbl_mpi_post_recv(comm_data *cd
				, int dim2
				)
{
  int ncommdomains  = cd->ncommdomains;
  int *commpartner  = cd->commpartner;
  int *recvcount    = cd->recvcount;

  int i;
  int size, szd = sizeof(double);

  /* recv */
  for(i = 0; i < ncommdomains; i++)
    {
      int k = commpartner[i];
      int count = recvcount[k];

      if(count > 0)
	{
	  double *rbuf = cd->recvbuf[i];
	  size  = count * szd * dim2;
	  MPI_Irecv(rbuf
		    , size
		    , MPI_BYTE
		    , k
		    , DATAKEY
		    , MPI_COMM_WORLD
		    , &(cd->req[i])
		    );
	}
    }

}

void exchange_dbl_mpi_pack(RangeList *color
			   , comm_data *cd
			   , double *data
			   , int dim2
			   )
{
  int i;
  for(i = 0; i < color->nsendcount; i++)
  {
      /* i1 = commdomain , 0 <= i1 < ncommdomains */
      int i1 = color->sendpartner[i];
      double *sbuf = cd->sendbuf[i1];
      
      /*  k = target rank */
      int k  = cd->commpartner[i1];
      if (color->sendcount[i] > 0 && cd->sendcount[k] > 0)
      {
	  cd->send_counter[i1].global += color->sendcount[i];
	  if(cd->send_counter[i1].global % cd->sendcount[k] == 0)
	  {
	      exchange_dbl_copy_in(cd
				   , sbuf
				   , data
				   , dim2
				   , i1
		  );
#ifdef USE_COMM_OVERLAP
	      exchange_dbl_mpi_send(cd, data, dim2, i1);
#endif
	  }
      }
  }
}


void exchange_dbl_mpi_exchg(comm_data *cd
			    , double *data
			    , int dim2
			    )
{
    int ncommdomains  = cd->ncommdomains;
    int nreq          = cd->nreq;

    int *sendcount    = cd->sendcount;
    int *recvcount    = cd->recvcount;
    int **sendindex   = cd->sendindex;
    int **recvindex   = cd->recvindex;
    int i;

    ASSERT(dim2 > 0);
    ASSERT(ncommdomains != 0);
    ASSERT(nreq > 0);

    ASSERT(sendcount != NULL);
    ASSERT(recvcount != NULL);
    ASSERT(sendindex != NULL);
    ASSERT(recvindex != NULL);

#ifndef USE_COMM_OVERLAP
    for(i = 0; i < ncommdomains; i++)
    {
	int *commpartner  = cd->commpartner;
	int k = commpartner[i];
	int count = sendcount[k];
	if (count > 0)
	{
	    exchange_dbl_mpi_send(cd, data, dim2, i);
	}
    }
#endif

    MPI_Waitall(ncommdomains
		, &(cd->req[0])
		, &(cd->stat[0])
	);
  
    for(i = 0; i < ncommdomains; i++)
    {
	/* copy the data from the recvbuf into out data field */
	double *rbuf = cd->recvbuf[i];
	exchange_dbl_copy_out(cd, rbuf, data, dim2, i);
    }

    MPI_Waitall(ncommdomains
		, &(cd->req[ncommdomains])
		, &(cd->stat[ncommdomains])
	);

}

