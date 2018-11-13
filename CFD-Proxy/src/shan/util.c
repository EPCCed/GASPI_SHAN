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
#include <sys/time.h>
#include "error_handling.h"
#include "solver_data.h"
#include "util.h"

void check_free(void *ptr)
{
  if(ptr != NULL)
  {
    free(ptr);
  }

  return;

}

void *check_malloc(size_t bytes)
{
  void *tmp;
  ASSERT(bytes > 0);
  tmp = malloc(bytes);
  ASSERT(tmp != NULL);

  return tmp;
}



void *check_realloc(void *old, size_t bytes)
{
  void *tmp;
  ASSERT(bytes > 0);
  tmp = realloc(old, bytes);
  ASSERT(tmp != NULL);

  return tmp;
}


static void swap(double *a, double *b)
{
  double tmp = *a;
  *a = *b;
  *b = tmp;
}

void sort_median(double *begin, double *end)
{
  double *ptr;
  double *split;
  if (end - begin <= 1)
    return;
  ptr = begin;
  split = begin + 1;
  while (++ptr != end) {
    if (*ptr < *begin) {
      swap(ptr, split);
      ++split;
    }
  }
  swap(begin, split - 1);
  sort_median(begin, split - 1);
  sort_median(split, end);
}


static void global_msort(int *pm
                         , int *pm_tmp
			 , int *htype
                         , int n
                         , int fp[][2]
                         , int fp_tmp[][2]
                         )
{
  int n1, n2;
  int *pm1, *pm2, *pm0;
  int (*fp1)[2], (*fp2)[2], (*fp0)[2];

  if(n <= 1)
    return;
  
  n1  = n / 2;
  n2  = n - n1;

  pm1 = pm; 
  pm2 = pm + n1; 
  fp1 = fp; 
  fp2 = fp + n1;

  if(n1 > 1)
    global_msort(pm1, pm_tmp, htype, n1, fp1, fp_tmp);

  if(n2 > 1)
    global_msort(pm2, pm_tmp, htype, n2, fp2, fp_tmp);

  pm0 = pm_tmp;
  fp0 = fp_tmp;

  while(n1 > 0 && n2 > 0)
    {
      int diff = 0;

      int fp11 = (*fp1)[1];
      int fp10 = (*fp1)[0];
      int fp21 = (*fp2)[1];
      int fp20 = (*fp2)[0]; 

      int h11  = htype[fp11];
      int h10  = htype[fp10];
      int h21  = htype[fp21];
      int h20  = htype[fp20];
      
      if(diff == 0) 
        diff  = MAX(h21,h20) - MAX(h11,h10);
      
      if(diff == 0) 
        diff  = (fp21 - fp11);
    
      if(diff == 0) 
        diff  = (fp20 - fp10);

      if(diff < 0)
        {
          (*pm0) = (*pm1);
          (*fp0)[0] = (*fp1)[0];
          (*fp0)[1] = (*fp1)[1];
          pm1++; fp1++; n1--;
        }
      else
        {
          (*pm0) = (*pm2);
          (*fp0)[0] = (*fp2)[0];
          (*fp0)[1] = (*fp2)[1];
          pm2++; fp2++; n2--;
        }
      pm0++;
      fp0++;
    }
  
  n2 = n - n2;
  for(n = 0; n < n1; n++)
    {
      pm0[n]    = pm1[n];
      fp0[n][0] = fp1[n][0];
      fp0[n][1] = fp1[n][1];
    }
  for(n = 0; n < n2; n++)
    {
      pm[n]        = pm_tmp[n];
      fp[n][0]     = fp_tmp[n][0];
      fp[n][1]     = fp_tmp[n][1];
    }
}



void sort_faces(int pm[]
		, int fpoint[][2]
                , int htype[]
		, int nfaces
		)
{

  int *pm_tmp      = check_malloc(nfaces * sizeof(int));
  int (*fp)[2]     = check_malloc(2 * nfaces * sizeof(int));
  int (*fp_tmp)[2] = check_malloc(2 * nfaces * sizeof(int));
  int face;

  for(face = 0; face < nfaces; face++) 
    {
      fp[face][0] = fpoint[face][0];
      fp[face][1] = fpoint[face][1];
    }

  global_msort(pm
               , pm_tmp
	       , htype
               , nfaces
               , fp
               , fp_tmp
               );
  
  check_free(pm_tmp);
  check_free(fp);
  check_free(fp_tmp);

}


double now()
{
  struct timeval tp;
  gettimeofday(&tp,NULL);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}




