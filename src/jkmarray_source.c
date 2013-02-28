#include <stdlib.h>
#include "memory.h"
#include "jkmarray.h"
#include "jkmarray_@@TYPE@@.h"

struct _JKMarray_@@TYPE@@
{
  @@TYPE@@ (*get) (struct _JKMarray_@@TYPE@@ * self, const int J, const int K,
		 const int M);
  void (*set) (struct _JKMarray_@@TYPE@@ * self, const int J, const int K,
	       const int M, const @@TYPE@@ val);
  @@TYPE@@ *data;
  int Jmax;
  int dim;
};

JKMarray_@@TYPE@@_t *
JKMarray_@@TYPE@@_ctor (const int Jmax)
{
  unsigned int dim = JKMarray_dim (Jmax);
  JKMarray_@@TYPE@@_t *a;

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

  a->get = &JKMarray_@@TYPE@@_get;
  a->set = &JKMarray_@@TYPE@@_set;
  a->Jmax = Jmax;
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
JKMarray_@@TYPE@@_get (JKMarray_@@TYPE@@_t * a, const int J, const int K, const int M)
{
  return a->data[JKMarray_idx (J, K, M)];
}

void
JKMarray_@@TYPE@@_set (JKMarray_@@TYPE@@_t * a, const int J, const int K, const int M,
	      const @@TYPE@@ val)
{
  a->data[JKMarray_idx (J, K, M)] = val;
}

