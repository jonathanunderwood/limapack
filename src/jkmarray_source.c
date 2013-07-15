#include <stdlib.h>
#include "memory.h"
#include "jkmarray.h"
#include "jkmarray_@@TYPE@@.h"

struct _JKMarray_@@TYPE@@
{
  @@TYPE@@ (*get) (struct _JKMarray_@@TYPE@@ * self, const int two_J, const int two_K,
		   const int two_M);
  void (*set) (struct _JKMarray_@@TYPE@@ * self, const int two_J, const int two_K,
	       const int two_M, const @@TYPE@@ val);
  @@TYPE@@ *data;
  int two_Jmax;
  int dim;
};

JKMarray_@@TYPE@@_t *
JKMarray_@@TYPE@@_ctor (const int two_Jmax)
{
  unsigned int dim;
  JKMarray_@@TYPE@@_t *a;

  if (MEMORY_ALLOC(a) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  dim = JKMarray_dim (two_Jmax);

  if (MEMORY_ALLOC_N (a->data, dim) < 0)
    {
      MEMORY_OOMERR;
      MEMORY_FREE(a);
      return NULL;
    }

  a->get = &JKMarray_@@TYPE@@_get;
  a->set = &JKMarray_@@TYPE@@_set;

  a->two_Jmax = two_Jmax;
  a->dim = dim;

  return a;
}

void
JKMarray_@@TYPE@@_dtor (JKMarray_@@TYPE@@_t * a)
{
  MEMORY_FREE (a->data);
  MEMORY_FREE (a);
}

@@TYPE@@ 
JKMarray_@@TYPE@@_get (JKMarray_@@TYPE@@_t * a, const int two_J, const int two_K, const int two_M)
{
  return a->data[JKMarray_idx (two_J, two_K, two_M)];
}

void 
JKMarray_@@TYPE@@_set (JKMarray_@@TYPE@@_t * a, const int two_J, const int two_K, const int two_M,
		       const @@TYPE@@ val)
{
  a->data[JKMarray_idx (two_J, two_K, two_M)] = val;
}


