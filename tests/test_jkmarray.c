/* gcc -O2 test_jkmarray.c  ../src/memory.c -I ../src */

#include <stdlib.h>
#include <stdio.h>
#include "jkmarray.h"
#include "jkmarray_int.h"
#include "jkmarray_int.c"

#define TWO_JMAX_INT 20
#define TWO_JMAX_HALF_INT 21

int
main ()
{
  int two_jmin;

  /* This will throw errors if run under valgrind and there's something wrong
     with the memory access. */

  for (two_jmin = 0; two_jmin <= 1; two_jmin++)
    {
      int two_j, two_k, two_m, two_jmax;
      JKMarray_int_t *arr;
      unsigned long i;

      if (two_jmin == 0)
	two_jmax = TWO_JMAX_INT;
      else
	two_jmax = TWO_JMAX_HALF_INT;

      arr = JKMarray_int_ctor (two_jmax);

      if (arr == NULL)
	{
	  fprintf (stderr, "Error: failed to allocate jkmarray\n");
	  exit (EXIT_FAILURE);
	}

      i = 0;
      for (two_j = two_jmin; two_j <= two_jmax; two_j+= 2)
	for (two_k = -two_j; two_k <= two_j; two_k+= 2)
	  for (two_m = -two_j; two_m <= two_j; two_m+= 2)
	    {
	      int v, w, x;

	      JKMarray_int_set (arr, two_j, two_k, two_m, i);
	      v = JKMarray_int_get (arr, two_j, two_k, two_m);
	      w = arr->data[i];
	      x = JKMarray_idx (two_j, two_k, two_m);

	      if ((i != v) || (i != w) || (i != x))
		{
		  fprintf (stderr,
			   "Error with array indexing for (two_j=%d, two_k=%d, two_m=%d), i: %d v: %d w: %d x: %d\n",
			   two_j, two_k, two_m, i, v, w, x);
		  JKMarray_int_dtor (arr);
		  exit (EXIT_FAILURE);
		}
	      i++;
	    }
      JKMarray_int_dtor (arr);
    }

  return 0;
}

#undef TWO_JMAX_INT
#undef TWO_JMAX_HALF_INT
