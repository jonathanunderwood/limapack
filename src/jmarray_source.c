#include <stdlib.h>
#include "memory.h"
#include "jmarray.h"
#include "jmarray_@@TYPE@@.h"

struct _JMarray_@@TYPE@@
{
  @@TYPE@@ (*get) (struct _JMarray_@@TYPE@@ * self, const int two_J, const int two_M);
  void (*set) (struct _JMarray_@@TYPE@@ * self, const int two_J, const int two_M,
	       const @@TYPE@@ val);
  @@TYPE@@ *data;
  int two_Jmax;
  int dim;
};

JMarray_@@TYPE@@_t *
JMarray_@@TYPE@@_ctor (const int two_Jmax)
{
  unsigned int dim = JMarray_dim (two_Jmax);
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

  a->two_Jmax = two_Jmax;
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
JMarray_@@TYPE@@_get (JMarray_@@TYPE@@_t * a, const int two_J, const int two_M)
{
  return a->data[JMarray_idx(two_J, two_M)];
}

void
JMarray_@@TYPE@@_set (JMarray_@@TYPE@@_t * a, const int two_J, const int two_M, const @@TYPE@@ val)
{
  a->data[JMarray_idx(two_J, two_M)] = val;
}


