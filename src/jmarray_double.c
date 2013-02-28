#include <stdlib.h>
#include "memory.h"
#include "jmarray.h"
#include "jmarray_double.h"

struct _JMarray_double
{
  double (*get) (struct _JMarray_double * self, const int J, const int M);
  void (*set) (struct _JMarray_double * self, const int J, const int M,
	       const double val);
  double *data;
  int Jmax;
  int dim;
};

JMarray_double_t *
JMarray_double_ctor (const int Jmax)
{
  unsigned int dim = JMarray_dim (Jmax);
  JMarray_double_t *a;

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

  a->get = &JMarray_double_get;
  a->set = &JMarray_double_set;
  a->Jmax = Jmax;
  a->dim = dim;

  return a;
}

void
JMarray_double_dtor (JMarray_double_t * a)
{
  MEMORY_FREE (a->data);
  MEMORY_FREE (a);
}

double
JMarray_double_get (JMarray_double_t * a, const int J, const int M)
{
  int idx = JMarray_idx (J, M);
  return a->data[idx];
}

void
JMarray_double_set (JMarray_double_t * a, const int J, const int M,
		    const double val)
{
  int idx = JMarray_idx (J, M);
  a->data[idx] = val;
}
