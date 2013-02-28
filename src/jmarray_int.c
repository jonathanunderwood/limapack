#include <stdlib.h>
#include "memory.h"
#include "jmarray.h"
#include "jmarray_int.h"

struct _JMarray_int
{
  int (*get) (struct _JMarray_int * self, const int J, const int M);
  void (*set) (struct _JMarray_int * self, const int J, const int M,
	       const int val);
  int *data;
  int Jmax;
  int dim;
};

JMarray_int_t *
JMarray_int_ctor (const int Jmax)
{
  unsigned int dim = JMarray_dim (Jmax);
  JMarray_int_t *a;

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

  a->get = &JMarray_int_get;
  a->set = &JMarray_int_set;
  a->Jmax = Jmax;
  a->dim = dim;

  return a;
}

void
JMarray_int_dtor (JMarray_int_t * a)
{
  MEMORY_FREE (a->data);
  MEMORY_FREE (a);
}

int
JMarray_int_get (JMarray_int_t * a, const int J, const int M)
{
  int idx = JMarray_idx (J, M);
  return a->data[idx];
}

void
JMarray_int_set (JMarray_int_t * a, const int J, const int M, const int val)
{
  int idx = JMarray_idx (J, M);
  a->data[idx] = val;
}
