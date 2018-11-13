#ifndef GRADIENTS_H
#define GRADIENTS_H

#include "comm_data.h"
#include "solver_data.h"

void compute_gradients_gg_comm_free(comm_data *cd, solver_data *sd);

void compute_gradients_gg_shan(comm_data *cd, solver_data *sd);

#endif
