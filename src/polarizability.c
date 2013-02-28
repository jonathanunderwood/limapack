#include <stdlib.h>
#include <math.h>

#include "polarizability.h"

polarizability_t *
polarizability_ctor()
{
  polarizability_t *alpha = JMarray_double_ctor(2);
  return alpha;
}

void
polarizability_dtor(polarizability_t *alpha)
{
  JMarray_double_dtor(alpha);
}

polarizability_t *
polarizability_from_cart_ctor (const double a_xx, 
			       const double a_yy, 
			       const double a_zz)
{
  polarizability_t *p = polarizability_ctor();
    
  if (p == NULL)
    return p;

  JMarray_double_set(p, 0, 0, -(a_xx + a_yy + a_zz) / sqrt (3.0));
  JMarray_double_set(p, 1, -1, 0.0);
  JMarray_double_set(p, 1, 0, 0.0);
  JMarray_double_set(p, 1, 1, 0.0);
  JMarray_double_set(p, 2, -2, 0.5 * (a_xx - a_yy));
  JMarray_double_set(p, 2, -1, 0.0);
  JMarray_double_set(p, 2, 0, (2.0 * a_zz - a_xx - a_yy) / sqrt (6.0));
  JMarray_double_set(p, 2, 1, 0.0);
  JMarray_double_set(p, 2, 2, 0.5 * (a_xx - a_yy));
  
  return p;
}

double
polarizability_get (polarizability_t *alpha, const int k, const int q)
{
  return JMarray_double_get (alpha, k, q);
}
