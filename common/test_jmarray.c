#include <stdlib.h>
#include <stdio.h>
#include "jmarray.h"

#define JMAX 5

int
main ()
{
  JMarray arr;
  int ret;
  int j, k, m;
  unsigned long i;

  ret = JMarray_init (&arr, JMAX);
  if (ret)
    {
      fprintf (stderr, "Error: failed to allocate jmarray\n");
      exit (ret);
    }

  /* This will throw errors if run under valgrind and there's something wrong
     with the memory allocation. */
  i = 0;
  for (j = 0; j <= JMAX; j++)
      for (m = -j; m <= j; m++)
	{
	  int v;
	  arr.set (&arr, j, m, (double) i);
	  v = arr.get (&arr, j, m);
	  if (i != (int) v)
	    {
	      fprintf (stderr,
		       "Error with array indexing for (j=%d, m=%d), i: %d which should be %d\n",
		       j, m, v, i);
	      JMarray_free (&arr);
	      exit (1);
	    }
	  i++;
	}

  JMarray_free (&arr);

  return ret;
}
