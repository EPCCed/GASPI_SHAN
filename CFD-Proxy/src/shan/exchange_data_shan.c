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


void init_shan_requests(comm_data *cd)
{
  int i;
  ASSERT(cd->ncommdomains > 0);

  /* status flag */
  cd->send_counter = check_malloc(cd->ncommdomains*sizeof(counter_t));
  cd->recv_counter = check_malloc(cd->ncommdomains*sizeof(counter_t));
  for(i = 0; i < cd->ncommdomains; i++)
    {
      cd->send_counter[i].global = 0;
      cd->recv_counter[i].global = 0;
    }

}



#ifndef USE_COMM_OVERLAP
static void exchange_dbl_shan_send(comm_data *cd
				  , shan_segment_t *const data_segment
				  , int i
				  )
{
    int res, type_id = 0;
    ASSERT((res = shan_comm_notify_or_write(&(cd->neighborhood_id)
					    , data_segment
					    , type_id
					    , i
		)) == SHAN_SUCCESS);
}
#endif



void exchange_dbl_shan_pack(RangeList *color
			   , comm_data *cd
			   , shan_segment_t *const data_segment
			   )
{
  int i;
  for(i = 0; i < color->nsendcount; i++)
  {
      /* i1 = commdomain , 0 <= i1 < ncommdomains */
      int i1 = color->sendpartner[i];
      
      /*  k = target rank */
      int k  = cd->commpartner[i1];
      if (color->sendcount[i] > 0 && cd->sendcount[k] > 0)
      {
	  cd->send_counter[i1].global += color->sendcount[i];
	  if(cd->send_counter[i1].global % cd->sendcount[k] == 0)
	  {
	      int res, type_id = 0;
	      ASSERT((res = shan_comm_notify_or_write(&(cd->neighborhood_id)
						   , data_segment
						   , type_id
						   , i1
			  )) == SHAN_SUCCESS);
	  }
      }
  }
}


void exchange_dbl_shan_exchg(comm_data *cd
			     , shan_segment_t *const data_segment
    )
{
#ifndef USE_COMM_OVERLAP
    int i;
    for(i = 0; i < cd->ncommdomains; i++)
    {
	exchange_dbl_shan_send(cd, data_segment, i);
    }
#endif
    
    int res, type_id = 0;
    ASSERT((res = shan_comm_wait4All(&(cd->neighborhood_id)
				     , data_segment
				     , type_id
		)) == SHAN_SUCCESS);

}


