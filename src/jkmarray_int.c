#include <stdlib.h>
#include "memory.h"
#include "jkmarray.h"
#include "jkmarray_int.h"

struct _JKMarray_int
{
  int (*get) (struct _JKMarray_int * self, const int J, const int K,
	      const int M);
  void (*set) (struct _JKMarray_int * self, const int J, const int K,
	       const int M, const int val);
  int *data;
  int Jmax;
  int dim;
};

JKMarray_int_t *
JKMarray_int_ctor (const int Jmax)
{
  unsigned int dim = JKMarray_dim (Jmax);
  JKMarray_int_t *a;

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

  a->get = &JKMarray_int_get;
  a->set = &JKMarray_int_set;
  a->Jmax = Jmax;
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
JKMarray_int_get (JKMarray_int_t * a, const int J, const int K, const int M)
{
  return a->data[JKMarray_idx (J, K, M)];
}

void
JKMarray_int_set (JKMarray_int_t * a, const int J, const int K, const int M,
		  const int val)
{
  a->data[JKMarray_idx (J, K, M)] = val;
}
