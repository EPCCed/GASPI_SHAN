#ifndef POINTS_OF_COLOR_H
#define POINTS_OF_COLOR_H

#include <omp.h>
#include <stdbool.h>

#include "comm_data.h"
#include "solver_data.h"
#include "error_handling.h"

void set_last_points_of_color(solver_data *sd
			      );

void set_first_points_of_color(solver_data *sd
			       );

#endif
