#include <stdlib.h>
#include <stdio.h>
#include "jkmarray.h"

#define JMAX 5

int
main ()
{
  JKMarray_t *arr;
  int j, k, m;
  unsigned long i;

  arr = JKMarray_ctor (JMAX);
  if (arr == NULL)
    {
      fprintf (stderr, "Error: failed to allocate jkmarray\n");
      exit (EXIT_FAILURE);
    }

  /* This will throw errors if run under valgrind and there's something wrong
     with the memory access. */
  i = 0;
  for (j = 0; j <= JMAX; j++)
    for (k = -j; k <= j; k++)
      for (m = -j; m <= j; m++)
	{
	  int v, w;
	  arr->set (arr, j, k, m, (double) i);
	  v = arr->get (arr, j, k, m);
	  w = arr->data[i];
	  if ((i != (int) v) || (i != (int) w))
	    {
	      fprintf (stderr,
		       "Error with array indexing for (j=%d, k=%d, m=%d), i: %d v: %d w: %d\n",
		       j, k, m, (int) i, (int)v, (int)w);
	      JKMarray_dtor (arr);
	      exit (EXIT_FAILURE);
	    }
	  i++;
	}

  JKMarray_dtor (arr);

  return 0;
}
