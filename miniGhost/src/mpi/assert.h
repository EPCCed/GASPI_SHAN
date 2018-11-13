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
