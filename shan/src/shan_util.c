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


#include <mpi.h>
#include <fcntl.h>
#include <sys/shm.h>
#include <sys/mman.h>
#include <unistd.h>
#include <stdio.h>

#include "assert.h"


void check_free(void *ptr)
{
  if(ptr != NULL)
    {
      free(ptr);
    }

  return;
}

void *check_malloc(long bytes)
{
  void *tmp;
  ASSERT(bytes > 0);
  tmp = malloc(bytes);
  ASSERT(tmp != NULL);

  return tmp;
}

