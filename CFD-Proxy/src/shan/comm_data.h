#ifndef COMM_DATA_H
#define COMM_DATA_H

#include <mpi.h>

#include "solver_data.h"


typedef void (*send_fn)(RangeList *color
			, comm_data *cd
			, shan_segment_t *const data_segment
			);

typedef void (*exch_fn)(comm_data *cd
			, shan_segment_t *const data_segment
			);


void initSharedNotification(comm_data *cd, solver_data *sd);
void init_communication(int argc, char *argv[], comm_data *cd);
void read_communication_data(int ncid, comm_data *cd);
void compute_communication_tables(comm_data *cd);
void free_ressources(comm_data *cd, solver_data *sd);

#endif
