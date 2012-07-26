#include <stdlib.h>
#include <stdio.h>
#include "jkmarray.h"

#define JMAX 5

int
main ()
{
  JKMarray arr;
  int ret;
  int j, k, m;
  unsigned long i;

  ret = JKMarray_init (&arr, JMAX);
  if (ret)
    {
      fprintf (stderr, "Error: failed to allocate jkmarray\n");
      exit (ret);
    }

  /* This will throw errors if run under valgrind and there's something wrong
     with the memory allocation. */
  i = 0;
  for (j = 0; j <= JMAX; j++)
    for (k = -j; k <= j; k++)
      for (m = -j; m <= j; m++)
	{
	  int v;
	  arr.set (&arr, j, k, m, (double) i);
	  v = arr.get (&arr, j, k, m);
	  if (i != (int) v)
	    {
	      fprintf (stderr,
		       "Error with array indexing for (j=%d, k=%d, m=%d), i: %d which should be %d\n",
		       j, k, m, v, i);
	      JKMarray_free (&arr);
	      exit (1);
	    }
	  i++;
	}

  JKMarray_free (&arr);

  return ret;
}
