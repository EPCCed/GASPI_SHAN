#ifndef RANGELIST_H
#define RANGELIST_H

#include <stdbool.h>
#include "comm_data.h"
#include "solver_data.h"
#include "error_handling.h"


void init_rangelist(comm_data *cd
		    , solver_data *sd
		    );

RangeList* private_get_color(solver_data *sd
			     , RangeList *const prev
			     );

RangeList* private_get_color_and_exchange(solver_data *sd
					  , RangeList *const prev
					  , pre_fn pre
					  , send_fn send
					  , exch_fn exch
					  , comm_data *cd
					  , double *data
					  , int dim2
					  );

static inline RangeList* get_color_and_exchange(solver_data *sd
						, pre_fn pre
						, send_fn send
						, exch_fn exch
						, comm_data *cd
						, double *data
						, int dim2
						)
{
  return private_get_color_and_exchange(sd
					, NULL
					, pre
					, send
					, exch
					, cd
					, data
					, dim2
					);
}

static inline RangeList* get_next_color_and_exchange(solver_data *sd
						     , RangeList *const prev
						     , pre_fn pre
						     , send_fn send
						     , exch_fn exch
						     , comm_data *cd
						     , double *data
						     , int dim2
						     )
{
  return private_get_color_and_exchange(sd
					, prev
					, pre
					, send
					, exch
					, cd
					, data
					, dim2
					);
}


static inline RangeList* get_color(solver_data *sd)
{
  return private_get_color(sd, NULL);
}

static inline RangeList* get_next_color(solver_data *sd
					, RangeList *const prev
					)
{
  return private_get_color(sd, prev);
}

#endif
