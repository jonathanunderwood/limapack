#include <stdlib.h>
#include "memory.h"
#include "jkmarray.h"

#define IDX(J, K, M) (((4 * J * J + 5) * J) / 3 + 2 * J * J + K * (2 * J + 1) + M)
#define DIM(J) (((4 * J * J + 11) * J) / 3 + (4 * J * J) + 1)

JKMarray_t *
JKMarray_ctor ( const int Jmax)
{
  unsigned int dim = DIM (Jmax);
  JKMarray_t *a;

  if (MEMORY_ALLOC(a) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  if (MEMORY_ALLOC_N (a->data, dim) < 0)
    {
      MEMORY_OOMERR;
      MEMORY_FREE(a);
      return NULL;
    }

  a->get = &JKMarray_get;
  a->set = &JKMarray_set;
  a->Jmax = Jmax;
  a->dim = dim;

  return a;
}

void
JKMarray_dtor (JKMarray_t * a)
{
  MEMORY_FREE (a->data);
  MEMORY_FREE (a);
}

double
JKMarray_get (JKMarray_t * a, const int J, const int K, const int M)
{
  return a->data[IDX (J, K, M)];
}

void
JKMarray_set (JKMarray_t * a, const int J, const int K, const int M,
	      const double val)
{
  a->data[IDX (J, K, M)] = val;
}

#undef IDX
#undef DIM
