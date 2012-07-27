#include <stdlib.h>
#include <stdio.h>
#include "jmarray.h"

#define JMAX 5

int
main ()
{
  JMarray_t *arr;
  int j, m;
  unsigned long i;

  arr = JMarray_ctor (JMAX);
  if (arr == NULL)
    {
      fprintf (stderr, "Error: failed to allocate jmarray\n");
      exit (EXIT_FAILURE);
    }

  /* This will throw errors if run under valgrind and there's something wrong
     with the memory allocation. */
  i = 0;
  for (j = 0; j <= JMAX; j++)
      for (m = -j; m <= j; m++)
	{
	  int v, w;
	  arr->set (arr, j, m, (double) i);
	  v = arr->get (arr, j, m);
	  w = arr->data[i];
	  if ((i != (int) v) || i != (int) w)
	    {
	      fprintf (stderr,
		       "Error with array indexing for (j=%d, m=%d). i: %d v: %d w: %d\n",
		       j, m, (int) i, (int)v, (int)w);
	      JMarray_dtor (arr);
	      exit (EXIT_FAILURE);
	    }
	  i++;
	}

  JMarray_dtor (arr);

  return 0;
}
