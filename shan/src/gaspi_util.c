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


#include <stdio.h>
#include "GASPI.h"
#include "GASPI_Ext.h"
#include "assert.h"

void
write_notify_and_wait ( gaspi_segment_id_t segment_id
			, gaspi_offset_t const offset_local
			, gaspi_rank_t const rank
			, gaspi_offset_t const offset_remote
			, gaspi_size_t const size
			, gaspi_notification_id_t const notification_id
			, gaspi_notification_t const notification_value
			, gaspi_queue_id_t const queue
			)
{
  gaspi_timeout_t const timeout = GASPI_BLOCK;
  gaspi_return_t ret;
  
  /* write, wait if required and re-submit */
  while ((ret = ( gaspi_write_notify( segment_id
				      , offset_local
				      , rank
				      , segment_id
				      , offset_remote
				      , size
				      , notification_id
				      , notification_value
				      , queue
				      , timeout
				      )
		  )) == GASPI_QUEUE_FULL)
    {
      SUCCESS_OR_DIE (gaspi_wait (queue,
				  GASPI_BLOCK));
    }

  ASSERT (ret == GASPI_SUCCESS);
}

