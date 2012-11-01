#include <stdlib.h>
#include "memory.h"
#include "jmarray.h"

JMarray_t *
JMarray_@@TYPE@@_ctor (const int Jmax)
{
  unsigned int dim = JMarray_dim (Jmax);
  JMarray_@@TYPE@@_t * a;

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

  a->get = &JMarray_@@TYPE@@_get;
  a->set = &JMarray_@@TYPE@@_set;
  a->Jmax = Jmax;
  a->dim = dim;

  return a;
}

void
JMarray_@@TYPE@@_dtor (JMarray_@@TYPE@@_t * a)
{
  MEMORY_FREE (a->data);
  MEMORY_FREE (a);
}

@@TYPE@@
JMarray_@@TYPE@@_get (JMarray_@@TYPE@@_t * a, const int J, const int M)
{
  int idx = JMarray_idx (J, M);
  return a->data[idx];
}

void
JMarray_@@TYPE@@_set (JMarray_@@TYPE@@_t * a, const int J, const int M, const @@TYPE@@ val)
{
  int idx = JMarray_idx (J, M);
  a->data[idx] = val;
}



