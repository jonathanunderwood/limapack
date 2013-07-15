#include <stdlib.h>
#include "memory.h"
#include "jkmarray.h"
#include "jkmarray_int.h"

struct _JKMarray_int
{
  int (*get) (struct _JKMarray_int * self, const int two_J, const int two_K,
	      const int two_M);
  void (*set) (struct _JKMarray_int * self, const int two_J, const int two_K,
	       const int two_M, const int val);
  int *data;
  int two_Jmax;
  int dim;
};

JKMarray_int_t *
JKMarray_int_ctor (const int two_Jmax)
{
  unsigned int dim;
  JKMarray_int_t *a;

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

  a->get = &JKMarray_int_get;
  a->set = &JKMarray_int_set;

  a->two_Jmax = two_Jmax;
  a->dim = dim;

  return a;
}

void
JKMarray_int_dtor (JKMarray_int_t * a)
{
  MEMORY_FREE (a->data);
  MEMORY_FREE (a);
}

int
JKMarray_int_get (JKMarray_int_t * a, const int two_J, const int two_K,
		  const int two_M)
{
  return a->data[JKMarray_idx (two_J, two_K, two_M)];
}

void
JKMarray_int_set (JKMarray_int_t * a, const int two_J, const int two_K,
		  const int two_M, const int val)
{
  a->data[JKMarray_idx (two_J, two_K, two_M)] = val;
}
