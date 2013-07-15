/* gcc -O2 test_jmarray.c  ../src/memory.c -I ../src */

#include <stdlib.h>
#include <stdio.h>
#include "jmarray.h"
#include "jmarray_int.h"
#include "jmarray_int.c"

#define TWO_JMAX_INT 20
#define TWO_JMAX_HALF_INT 21

int
main ()
{
  int two_jmin;

  for (two_jmin = 0; two_jmin <= 1; two_jmin++)
    {
      JMarray_int_t *arr;
      int two_j, two_m, two_jmax;
      unsigned long i;
      
      if (two_jmin == 0)
	two_jmax = TWO_JMAX_INT;
      else
	two_jmax = TWO_JMAX_HALF_INT;

      arr = JMarray_int_ctor (two_jmax);
      if (arr == NULL)
	{
	  fprintf (stderr, "Error: failed to allocate jmarray\n");
	  exit (EXIT_FAILURE);
	}
      
      /* This will throw errors if run under valgrind and there's
	 something wrong with the memory allocation. */
      i = 0;
      for (two_j = two_jmin; two_j <= two_jmax; two_j += 2)
	for (two_m = -two_j; two_m <= two_j; two_m += 2)
	  {
	    int v, w, x;
	    JMarray_int_set (arr, two_j, two_m, i);
	    v = arr->get (arr, two_j, two_m);
	    w = arr->data[i];
	    x = JMarray_idx (two_j, two_m);
	    if ((i != v || i != w || i != x))
	      {
		fprintf (stderr,
			 "Error with array indexing for (two_j=%d, two_m=%d). i: %d v: %d w: %d x: %d\n",
			 two_j, two_m, i, v, w, x);
		JMarray_int_dtor (arr);
		exit (EXIT_FAILURE);
	      }
	    i++;
	  }
      JMarray_int_dtor (arr);
    }

  return 0;
}
