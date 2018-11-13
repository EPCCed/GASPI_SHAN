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
#include <string.h>
#include <stdbool.h>

#include "read_netcdf.h"
#include "solver_data.h"
#include "error_handling.h"
#include "util.h"
#include "comm_util.h"
#include "rangelist.h"
#include "points_of_color.h"

#ifdef USE_CHACO
#include <chaco_interface.h>
#endif

#ifdef USE_CHACO
typedef struct
{
  int nelems;
  int *list;
  int *listidx;
} ConnectInfo;
#endif

#ifdef USE_CHACO
static int compute_chaco_graph(ConnectInfo *pnt2pnt
			       , const int nfaces
			       , const int nallpoints
			       , int (*fpoint)[2]
			       )
{
  int i;
  /*----------------------------------------------------------------------------
    | Build up index list for the adjacency graph
    ----------------------------------------------------------------------------*/
  pnt2pnt->listidx = check_malloc((nallpoints + 1)*sizeof(int));
  for(i = 0; i < nallpoints + 1; i++)
    {
      pnt2pnt->listidx[i] = 0;
    }

  for(i = 0; i < nfaces; i++)
    {
      const int p1 = fpoint[i][0];
      const int p2 = fpoint[i][1];

      if(p1 < nallpoints)
	{
	  pnt2pnt->listidx[p1 + 1]++;
	}

      if(p2 < nallpoints)
	{
	  pnt2pnt->listidx[p2 + 1]++;
	}
    }

  /*----------------------------------------------------------------------------
    | Sum up index vector
    ----------------------------------------------------------------------------*/
  for(i = 0; i < nallpoints; i++)
    {
      pnt2pnt->listidx[i + 1] += pnt2pnt->listidx[i];
    }
  /*----------------------------------------------------------------------------
    | Build up the adjacency graph
    ----------------------------------------------------------------------------*/
  pnt2pnt->list = check_malloc((pnt2pnt->listidx[nallpoints])*sizeof(int));
  for(i = 0; i < pnt2pnt->listidx[nallpoints]; i++)
    {
      pnt2pnt->list[i] = -1;
    }

  for(i = 0; i < nfaces; i++)
    {
      const int p1 = fpoint[i][0];
      const int p2 = fpoint[i][1];
      const int pos1 = pnt2pnt->listidx[p1];
      const int pos2 = pnt2pnt->listidx[p2];

      pnt2pnt->list[pos1] = p2 + 1;
      pnt2pnt->listidx[p1]++;

      pnt2pnt->list[pos2] = p1 + 1;
      pnt2pnt->listidx[p2]++;
    }

  for(i = nallpoints; i > 0; i--)
    {
      pnt2pnt->listidx[i] = pnt2pnt->listidx[i - 1];
    }

  pnt2pnt->listidx[0] = 0;
  pnt2pnt->nelems = nallpoints;

  return nallpoints;

}


static void init_pid(int const np
		     , int *pid
		     , solver_data *sd
		     )
{
  int nall = sd->nallpoints;
  int nfaces = sd->nfaces;
  int (*fpoint)[2] = sd->fpoint;
  int *pointweight = NULL;

  short *part;
  float eigent_tol = 0.001;
  int mesh_dims[3] = { 0,0,0};
  int arch   = 0;
  int scheme = 0;
  int coarse = 0;
  int local  = 0;
  int global = 0;
  int ncube  = 0;
  int i;

  ConnectInfo pnt2pnt;
  pnt2pnt.nelems = 0;
  pnt2pnt.listidx = NULL;
  pnt2pnt.list = NULL;

  /* coarsen KL */
  coarse=10;

  /* bi-section etc. */
  scheme=1;

  /* hypercube, flat network .. */
  arch=1;

  /* KL global */
  global=1;

  /* KL local */
  local=1;

  compute_chaco_graph(&pnt2pnt, nfaces, nall, fpoint);

  /* load (0,1,0) */
  pointweight = (int *)check_malloc(nall * sizeof(int));
  for(i = 0; i < nall; i++)
    {
      pointweight[i] = 0;
    }

  for(i = 0; i < nfaces; i++)
    {
      (pointweight[fpoint[i][0]])++;
      (pointweight[fpoint[i][1]])++;
    }

  part = check_malloc(pnt2pnt.nelems * sizeof(short));
  for(i = 0; i < pnt2pnt.nelems; i++)
    {
      part[i] = -1;
    }

  mesh_dims[0] = MIN(sd->nallpoints / np, 65535);

  interface(pnt2pnt.nelems, pnt2pnt.listidx, pnt2pnt.list, pointweight,
	    NULL, NULL, NULL, NULL, NULL, NULL,
	    part, arch, ncube, mesh_dims,
	    NULL, global, local, 0, coarse, scheme, eigent_tol, rand());


  for(i = 0; i < nall; i++)
    {
      ASSERT((part[i] < mesh_dims[0]) && (part[i] >= 0));
      pid[i] = part[i];
    }

  /* re-assign cross edges ending in addpoints (outer halo) */
  int face;
  for (face = 0; face < sd->nfaces; face++)
    {
      int p0 = sd->fpoint[face][0];
      int p1 = sd->fpoint[face][1];
      if (p1 >= sd->nownpoints && pid[p0] != pid[p1])
	{
	  pid[p1] = pid[p0];
	}
      if (p0 >= sd->nownpoints && pid[p1] != pid[p0])
	{
	  pid[p0] = pid[p1];
	}
    }

  sd->ncolors = mesh_dims[0];

}
#endif


static void init_halo_type(int *htype
			   , comm_data *cd
			   , solver_data *sd
			   )
{
  int i, j;
  for(i = 0; i < sd->nallpoints; i++) 
    {
      htype[i] = 0;
    }

  /* set halo type for sendindex data (inner halo) */
  for(i = 0; i < cd->ncommdomains; i++)
    {
      int k = cd->commpartner[i];
      int count = cd->sendcount[k];      
      if(count > 0)
	{
	  for(j = 0; j < count; j++)
	    {
	      int pnt = cd->sendindex[k][j];
	      htype[pnt] = 1;
	    }
	}
    }

  /* addpoints (outer halo) */
  for(i = sd->nownpoints; i < sd->nallpoints; i++) 
    {
      htype[i] = 1;
    }
}


static void local_init_rangelist(RangeList *fcolor)
{  
  // next slice - linked list
  fcolor->succ = NULL;

  // meta data
  fcolor->start = -1;
  fcolor->stop = -1;

  // points of color
  fcolor->nfirst_points_of_color = 0;
  fcolor->first_points_of_color = NULL;
  fcolor->nlast_points_of_color = 0;
  fcolor->last_points_of_color = NULL;

  // comm vars - send
  fcolor->nsendcount = 0;
  fcolor->sendpartner = NULL;
  fcolor->sendcount = NULL;
  fcolor->sendindex = NULL;
  fcolor->sendoffset = NULL;

  // comm vars - recv
  fcolor->nrecvcount = 0;
  fcolor->recvpartner = NULL;
  fcolor->recvcount = NULL;
  fcolor->recvindex = NULL;
  fcolor->recvoffset = NULL;

}


void set_rangelist(comm_data *cd
		   , solver_data *sd
		   , int *htype
		   )
{
  RangeList *tl;
  
  ASSERT(cd != NULL);
  ASSERT(sd != NULL);

  int nfaces          = sd->nfaces;
  int (*fpoint)[2]    = sd->fpoint;
  double(*fnormal)[3] = sd->fnormal;
  int face;

  /* init permutation vector for sorting */
  int *pm    = check_malloc(nfaces * sizeof(int));
  for (face = 0; face < nfaces; face++)
    {
      pm[face]    = face;
    }

  /* sort faces for type and p1/p0 */
  sort_faces(pm, fpoint, htype, nfaces);
  

  /* fix face permutation */
  int    (*fp)[2] = check_malloc(2 * nfaces * sizeof(int));
  double (*fn)[3] = check_malloc(3 * nfaces * sizeof(double));
  for (face = 0; face < nfaces; face++)
    {
      int tf = pm[face];
      memcpy(&(fp[face][0])
	     , &(fpoint[tf][0])
	     , 2 * sizeof(int)
	     );
      memcpy(&(fn[face][0])
	     , &(fnormal[tf][0])
	     , 3 * sizeof(double)
	     );
    }
  check_free(fpoint);
  check_free(fnormal);

  sd->fpoint = fp;
  sd->fnormal = fn;
  
  /* set rangelist */
#define MAX_FACES_IN_COLOR 96

  int count = 0;
  int ncolors = 0;
  for(face = 1 ; face < nfaces; face++)
    {
      if((++count) == MAX_FACES_IN_COLOR)
	{
	  count = 0;
	  ncolors++;
	}
    }
  /* last color */
  ncolors++;

  /* alloc rangelist */
  sd->ncolors = ncolors;
  sd->fcolor  = check_malloc(ncolors * sizeof(RangeList));

  int i;
  for(i = 0; i < sd->ncolors; i++)
    {
      RangeList *rl = &(sd->fcolor[i]); 
      local_init_rangelist(rl);
    }

  /* set color range */
  count = 0;
  int i0 = 0;
  int start = 0;
  for(face = 1 ; face < nfaces; face++)
    {
      if((++count) == MAX_FACES_IN_COLOR)
	{
	  tl = &(sd->fcolor[i0]);
	  tl->start  = start;
	  tl->stop   = face;

	  start = face;
	  count = 0;
	  i0++;
	}
    }
  /* last color */
  tl = &(sd->fcolor[i0]);
  tl->start  = start;
  tl->stop   = nfaces;
  i0++;
    
  /* set rangelist successor */
  for(i = 0; i < sd->ncolors; i++)
    {
      tl = &(sd->fcolor[i]);
      if (i < ncolors - 1 )
	{
	  tl->succ = &(sd->fcolor[i+1]);
	}
      else
	{
	  tl->succ = NULL;
	}
    }
  
    /* first points of color */
  set_first_points_of_color(sd);

  /* last points of color */
  set_last_points_of_color(sd);
  
#if 0
  for (tl = get_color(sd); tl != NULL
	 ; tl = get_next_color(sd,tl)) 
    {
      for (face = tl->start; face < tl->stop; face++)
	{
	  int p0 = sd->fpoint[face][0];
	  int p1 = sd->fpoint[face][1];
	  printf("iProc: %d face: %d p0/p1: %d/%d h0/h1: %d/%d\n",
		 cd->iProc,face,p0,p1,htype[p0],htype[p1]);
	}
    }
#endif
  
}


void init_rangelist(comm_data *cd
		  , solver_data *sd
		  )
{
  /* init halo type */
  int *htype = check_malloc(sd->nallpoints * sizeof(int));
  init_halo_type(htype, cd, sd);

  /* set/sort rangelist */
  set_rangelist(cd, sd, htype);

  /* free metadata */
  check_free(htype);

  /* assign comm for first/last points of color */
  init_rangelist_comm(cd, sd);
  
  if (cd->ndomains == 1)
    {
      return;
    }
}


RangeList* private_get_color_and_exchange(solver_data *sd
					  , RangeList *const prev
					  , pre_fn pre
					  , send_fn send
					  , exch_fn exch
					  , comm_data *cd
					  , double *data
					  , int dim2
					  )
{
    RangeList *color = NULL;

    if(prev != NULL)
    {
	if (send != NULL)
	{
	    send(prev
		 , cd
		 , data
		 , dim2
		);
	}
	color = prev->succ;
    }
    else
    {
	if (pre != NULL)
	{
	    pre(cd
		, dim2
		);
	}      
	color = sd->fcolor;
    }
  
    if (color == NULL)
    {
	if (exch != NULL)
	{
	    exch(cd
		 , data
		 , dim2
		);
	}
    }
    return color;

}

RangeList* private_get_color(solver_data *sd
			     , RangeList *const prev
			     )
{
  RangeList *color = NULL;

  if(prev != NULL)
    {
      color = prev->succ;
    }
  else
    {
      color = sd->fcolor;
    }

  return color;
}

