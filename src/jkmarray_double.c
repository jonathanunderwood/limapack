#include <stdlib.h>
#include "memory.h"
#include "jkmarray.h"
#include "jkmarray_double.h"

struct _JKMarray_double
{
  double (*get) (struct _JKMarray_double * self, const int two_J,
		 const int two_K, const int two_M);
  void (*set) (struct _JKMarray_double * self, const int two_J,
	       const int two_K, const int two_M, const double val);
  double *data;
  int two_Jmax;
  int dim;
};

JKMarray_double_t *
JKMarray_double_ctor (const int two_Jmax)
{
  unsigned int dim;
  JKMarray_double_t *a;

  if (MEMORY_ALLOC (a) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  dim = JKMarray_dim (two_Jmax);

  if (MEMORY_ALLOC_N (a->data, dim) < 0)
    {
      MEMORY_OOMERR;
      MEMORY_FREE (a);
      return NULL;
    }

  a->get = &JKMarray_double_get;
  a->set = &JKMarray_double_set;

  a->two_Jmax = two_Jmax;
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
JKMarray_double_get (JKMarray_double_t * a, const int two_J, const int two_K,
		     const int two_M)
{
  return a->data[JKMarray_idx (two_J, two_K, two_M)];
}

void
JKMarray_double_set (JKMarray_double_t * a, const int two_J, const int two_K,
		     const int two_M, const double val)
{
  a->data[JKMarray_idx (two_J, two_K, two_M)] = val;
}
