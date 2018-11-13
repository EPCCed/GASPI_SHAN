#ifndef EXCHANGE_DATA_H
#define EXCHANGE_DATA_H

#include "comm_data.h"

void init_shan_requests(comm_data *cd);

void exchange_dbl_shan_pack(RangeList *color
			   , comm_data *cd
			   , shan_segment_t *const data_segment
    );

void exchange_dbl_shan_exchg(comm_data *cd
			     , shan_segment_t *const data_segment
    );

#endif

