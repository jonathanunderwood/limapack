#include <stdlib.h>
#include "memory.h"
#include "jkmarray.h"
#include "jkmarray_double.h"

struct _JKMarray_double
{
  double (*get) (struct _JKMarray_double * self, const int J, const int K,
		 const int M);
  void (*set) (struct _JKMarray_double * self, const int J, const int K,
	       const int M, const double val);
  double *data;
  int Jmax;
  int dim;
};

JKMarray_double_t *
JKMarray_double_ctor (const int Jmax)
{
  unsigned int dim = JKMarray_dim (Jmax);
  JKMarray_double_t *a;

  if (MEMORY_ALLOC (a) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  if (MEMORY_ALLOC_N (a->data, dim) < 0)
    {
      MEMORY_OOMERR;
      MEMORY_FREE (a);
      return NULL;
    }

  a->get = &JKMarray_double_get;
  a->set = &JKMarray_double_set;
  a->Jmax = Jmax;
  a->dim = dim;

  return a;
}

void
JKMarray_double_dtor (JKMarray_double_t * a)
{
  MEMORY_FREE (a->data);
  MEMORY_FREE (a);
}

double
JKMarray_double_get (JKMarray_double_t * a, const int J, const int K,
		     const int M)
{
  return a->data[JKMarray_idx (J, K, M)];
}

void
JKMarray_double_set (JKMarray_double_t * a, const int J, const int K,
		     const int M, const double val)
{
  a->data[JKMarray_idx (J, K, M)] = val;
}
