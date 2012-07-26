#include <stdlib.h>
#include "memory.h"
#include "jmarray.h"

#define IDX(J, M) (J * J + J + M)
#define DIM(J) (J * J + 2 * J + 1)

JMarray_t *
JMarray_ctor (const int Jmax)
{
  unsigned int dim = DIM (Jmax);
  JMarray_t * a;

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

  a->get = &JMarray_get;
  a->set = &JMarray_set;
  a->Jmax = Jmax;
  a->dim = dim;

  return a;
}

void
JMarray_dtor (JMarray_t * a)
{
  MEMORY_FREE (a->data);
  MEMORY_FREE (a);
}

double
JMarray_get (JMarray_t * a, const int J, const int M)
{
  return a->data[IDX (J, M)];
}

void
JMarray_set (JMarray_t * a, const int J, const int M, const double val)
{
  a->data[IDX (J, M)] = val;
}

#undef IDX
#undef DIM
