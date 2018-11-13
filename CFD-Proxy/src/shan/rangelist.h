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
					  , send_fn send
					  , exch_fn exch
					  , comm_data *cd
					  , shan_segment_t *const data_segment
					  );

static inline RangeList* get_color_and_exchange(solver_data *sd
						, send_fn send
						, exch_fn exch
						, comm_data *cd
						, shan_segment_t *const data_segment
						)
{
  return private_get_color_and_exchange(sd
					, NULL
					, send
					, exch
					, cd
					, data_segment
					);
}

static inline RangeList* get_next_color_and_exchange(solver_data *sd
						     , RangeList *const prev
						     , send_fn send
						     , exch_fn exch
						     , comm_data *cd
						     , shan_segment_t *const data_segment
						     )
{
  return private_get_color_and_exchange(sd
					, prev
					, send
					, exch
					, cd
					, data_segment
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
