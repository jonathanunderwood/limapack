#include <stdlib.h>
#include <math.h>
#include "jmarray.h"

int polarizability_init (JMarray *p, const double a_xx, 
			 const double a_yy, const double a_zz)
{
  int ret = JMarray_init (p, 2);
  if (ret == 0)
    {
      p->set(p, 0, 0, -(a_xx + a_yy + a_zz) / sqrt (3.0));
      p->set(p, 1, -1, 0.0);
      p->set(p, 1, 0, 0.0);
      p->set(p, 1, 1, 0.0);
      p->set(p, 2, -2, 0.5 * (a_xx - a_yy));
      p->set(p, 2, -1, 0.0);
      p->set(p, 2, 0, (2.0 * a_zz - a_xx - a_yy) / sqrt (6.0));
      p->set(p, 2, 1, 0.0);
      p->set(p, 2, 2, 0.5 * (a_xx - a_yy));
    }
  return ret;
}

