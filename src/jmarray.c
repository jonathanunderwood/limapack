#include <stdlib.h>
#include "memory.h"
#include "jmarray.h"

#define IDX(J, M) (J * J + J + M)
#define DIM(J) (J * J + 2 * J + 1)

int
JMarray_init (JMarray * a, const int Jmax)
{
  unsigned int dim = DIM (Jmax);
  int ret = MEMORY_ALLOC_N (a->data, dim);

  if (ret == 0)
    {
      a->get = &JMarray_get;
      a->set = &JMarray_set;
      a->Jmax = Jmax;
      a->dim = dim;
      return 0;
    }
  else
    return 1;
}

void
JMarray_free (JMarray * a)
{
  MEMORY_FREE (a->data);
  a->Jmax = 0;
  a->dim = 0;
  a->get = NULL;
  a->set = NULL;
}

double
JMarray_get (JMarray * a, const int J, const int M)
{
  return a->data[IDX (J, M)];
}

void
JMarray_set (JMarray * a, const int J, const int M, const double val)
{
  a->data[IDX (J, M)] = val;
}

#undef IDX
#undef DIM
