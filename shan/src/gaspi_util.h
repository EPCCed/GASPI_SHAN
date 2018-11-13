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


#ifndef GASPI_UTIL_H
#define GASPI_UTIL_H

#include "GASPI.h"

void 
write_notify_and_wait ( gaspi_segment_id_t segment_id
			, gaspi_offset_t const offset_local
			, gaspi_rank_t const rank
			, gaspi_offset_t const offset_remote
			, gaspi_size_t const size
			, gaspi_notification_id_t const notification_id
			, gaspi_notification_t const notification_value
			, gaspi_queue_id_t const queue
			);

#endif
