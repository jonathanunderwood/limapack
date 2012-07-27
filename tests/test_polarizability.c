#include <stdlib.h>
#include <stdio.h>

#include "polarizability.h"

#define AXX 1.0
#define AYY 2.0
#define AZZ 3.0

int
main()
{
  polarizability_t *a = polarizability_ctor_from_cart (AXX, AYY, AZZ);

  if (a == NULL)
    {
      fprintf(stderr, "Failed to allocate polarizability.\n");
      exit (EXIT_FAILURE);
    }

  fprintf(stdout, "a_00: %f\n", a->get(a, 0, 0));
  fprintf(stdout, "a_10: %f\n", a->get(a, 1, 0));
  fprintf(stdout, "a_20: %f\n", a->get(a, 2, 0));
  fprintf(stdout, "a_22: %f\n", a->get(a, 2, 2));
  fprintf(stdout, "a_2-2: %f\n", a->get(a, 2, -2));

  polarizability_dtor(a);

  return 0;
}
