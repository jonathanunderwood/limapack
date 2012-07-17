#include <stdlib.h>
#include "memory.h"
#include "jkmarray.h"

#define IDX(J, K, M) (((4 * J * J + 5) * J) / 3 + 2 * J * J + K * (2 * J + 1) + M)
#define DIM(J) (((4 * J * J + 11) * J) / 3 + (4 * J * J) + 1)

int
JKMarray_init (JKMarray * a, const int Jmax)
{
  unsigned int dim = DIM (Jmax);
  int ret = MEMORY_ALLOC_N (a->data, dim);

  if (ret == 0)
    {
      a->get = &JKMarray_get;
      a->set = &JKMarray_set;
      a->Jmax = Jmax;
      a->dim = dim;
      return 0;
    }
  else
    return 1;
}

void
JKMarray_free (JKMarray * a)
{
  MEMORY_FREE (a->data);
  a->Jmax = 0;
  a->dim = 0;
  a->get = NULL;
  a->set = NULL;
}

double
JKMarray_get (JKMarray * a, const int J, const int K, const int M)
{
  return a->data[IDX (J, K, M)];
}

void
JKMarray_set (JKMarray * a, const int J, const int K, const int M,
	      const double val)
{
  a->data[IDX (J, K, M)] = val;
}

#undef IDX
#undef DIM
