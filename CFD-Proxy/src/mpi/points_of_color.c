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
#include <omp.h>
#include <stdbool.h>

#include "read_netcdf.h"
#include "solver_data.h"
#include "error_handling.h"
#include "util.h"
#include "rangelist.h"

void set_last_points_of_color(solver_data *sd
			      )
{
  int i, face;
  int    (*fpoint)[2]        = sd->fpoint;
  RangeList *color;

  ASSERT(sd->nallpoints > 0);
  ASSERT(fpoint != NULL);

  int *tmp1   = check_malloc(sd->nallpoints * sizeof(int));
  for(i = 0; i < sd->nallpoints; i++) 
    {
      tmp1[i] = 0;
    }

  /* all finalized/last points of color */
  for (color = get_color(sd); color != NULL
	 ; color = get_next_color(sd,color)) 
    {      
      for(face = color->start ; face < color->stop; face++)
        {
          int p0 = fpoint[face][0];
          int p1 = fpoint[face][1];
	  tmp1[p0]++;
	  tmp1[p1]++;
        }
    }

  int *tmp2   = check_malloc(sd->nallpoints * sizeof(int));
  for(i = 0; i < sd->nallpoints; i++) 
    {
      tmp2[i] = 0;
    }

  sd->first_points_of_color = check_malloc(sd->nallpoints * sizeof(int));
  int *points = sd->first_points_of_color;

  int ncolor = 0;
  int l = 0;
  for (color = get_color(sd); color != NULL
	 ; color = get_next_color(sd,color)) 
    {      
      int npoints = 0;
      for(face = color->start ; face < color->stop; face++)
        {
          int p0 = fpoint[face][0];
          int p1 = fpoint[face][1];
	  if (++tmp2[p0] == tmp1[p0])
	    {
	      points[npoints++] = p0;
	    }           
	  if (++tmp2[p1] == tmp1[p1])
	    {
	      points[npoints++] = p1;
	    }
	}

      l += npoints;
      color->nlast_points_of_color = npoints;
      color->last_points_of_color = points;      
      points += npoints;
      ncolor++;

    }

  ASSERT(l == sd->nallpoints);

  check_free(tmp2);  
  check_free(tmp1);  

}

void set_first_points_of_color(solver_data *sd)
{
  int i, face;
  RangeList *color;
  int    (*fpoint)[2] = sd->fpoint;
  
  /* all first points of color excluding outer halo */
  int *tmp1 = check_malloc(sd->nallpoints * sizeof(int));
  for (i = 0; i < sd->nallpoints; ++i)
    {
      tmp1[i] = -1;
    }
  
  sd->first_points_of_color = check_malloc(sd->nallpoints * sizeof(int));
  int *points = sd->first_points_of_color;

  int ncolor = 0;
  int l = 0;
  for (color = get_color(sd); color != NULL
	 ; color = get_next_color(sd,color)) 
    {      
      int npoints = 0;
      for(face = color->start ; face < color->stop; face++)
        {
          int p0 = fpoint[face][0];
          int p1 = fpoint[face][1];
	  if (tmp1[p0] == -1)
	    {
	      tmp1[p0] = ncolor;
	      points[npoints++] = p0;
	    }
	  if (tmp1[p1] == -1)
	    {
	      tmp1[p1] = ncolor;
	      points[npoints++] = p1;
	    }
	}
      
      l += npoints;
      color->nfirst_points_of_color = npoints;
      color->first_points_of_color = points;      
      points += npoints;
      ncolor++;
    }

  ASSERT(l == sd->nallpoints);

  check_free(tmp1);
}

