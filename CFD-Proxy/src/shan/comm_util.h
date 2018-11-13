#ifndef COMM_UTIL_H
#define COMM_UTIL_H

#include <omp.h>
#include <stdbool.h>
#include "comm_data.h"
#include "solver_data.h"
#include "error_handling.h"


void init_rangelist_comm(comm_data *cd
			 , solver_data *sd
    );

void free_rangelist_comm(solver_data *sd
    );

#endif
