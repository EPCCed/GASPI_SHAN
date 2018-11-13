/*
    Copyright (c) T-Systems SfR, C.Simmendinger <christian.simmendinger@t-systems.com>, 2018

    This file is part of SHAN.

    SHAN is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
	    the Free Software Foundation, either version 3 of the License, or
	        (at your option) any later version.

    SHAN is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
	    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	        GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
        along with SHAN.  If not, see <https://www.gnu.org/licenses/>.
	
*/


#ifndef SHAN_CORE_H
#define SHAN_CORE_H


#include <stdio.h>
#include <stdlib.h>

#include "GASPI.h"
#include "SHAN_segment.h"


#define NELEM_COMM_HEADER 3
#define ALIGNMENT 64


void shan_test_shared(shan_neighborhood_t *const neighborhood_id
		      , int const rank_local
		      , int const num_neighbors
		      , int const type_id
		      , int const idx
		      , int *rval
    );

int shan_increment_local(shan_neighborhood_t *const neighborhood_id				  
			 , int const type_id
			 , int const idx
    );

int shan_comm_waitsome_local(shan_neighborhood_t *const neighborhood_id
			     , int const type_id
			     , int const idx
    );

int shan_comm_waitsome_remote(shan_neighborhood_t *const neighborhood_id
			      , int const type_id
			      , int const idx
    );

#endif


