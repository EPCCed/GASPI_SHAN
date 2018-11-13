#ifndef EXCHANGE_DATA_H
#define EXCHANGE_DATA_H

#include "comm_data.h"

void init_mpi_requests(comm_data *cd
		       , int dim2
		       );

void exchange_dbl_mpi_post_recv(comm_data *cd
				, int dim2
				);

void exchange_dbl_mpi_pack(RangeList *color
			   , comm_data *cd
			   , double *data
			   , int dim2
			   );

void exchange_dbl_mpi_exchg(comm_data *cd
			    , double *data
			    , int dim2
			    );

#endif

