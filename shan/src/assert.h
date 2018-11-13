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


#ifndef ASSERT_H
#define ASSERT_H

#include <stdio.h>
#include <stdlib.h>

#define ASSERT(x...)                                                    \
  if (!(x))                                                             \
  {                                                                     \
    fprintf (stderr, "Error: '%s' [%s:%i]\n", #x, __FILE__, __LINE__);  \
    exit (EXIT_FAILURE);                                                \
  }


#define SUCCESS_OR_DIE(f...)                                            \
  do                                                                    \
  {                                                                     \
    const gaspi_return_t r = f;                                         \
    if (r != GASPI_SUCCESS)                                             \
     {                                                                  \
       printf ("Error: '%s' [%s:%i]: %i\n", #f, __FILE__, __LINE__, r); \
       exit (EXIT_FAILURE);                                             \
     }                                                                  \
  } while (0)


#endif
