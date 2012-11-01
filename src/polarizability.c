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

  p->set(p, 0, 0, -(a_xx + a_yy + a_zz) / sqrt (3.0));
  p->set(p, 1, -1, 0.0);
  p->set(p, 1, 0, 0.0);
  p->set(p, 1, 1, 0.0);
  p->set(p, 2, -2, 0.5 * (a_xx - a_yy));
  p->set(p, 2, -1, 0.0);
  p->set(p, 2, 0, (2.0 * a_zz - a_xx - a_yy) / sqrt (6.0));
  p->set(p, 2, 1, 0.0);
  p->set(p, 2, 2, 0.5 * (a_xx - a_yy));
  
  return p;
}

