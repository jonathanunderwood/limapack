#ifndef __LASER_TYPE1_H__
#define __LASER_TYPE1_H__
#include "laser.h"

typedef struct _laser_type1
{
  laser_t parent;
  /* Repeat these callbacks here (they're already part of the laser_t) so we can
     call them directly without a recast. */
  laser_polzn_vector * (*get_polzn_vector) (struct _laser *self, const double t);
  laser_polzn_tensor * (*get_polzn_tensor) (struct _laser *self, const double t);
  gsl_comlplex * (*get_envelope) (struct _laser *self, const double t);
  /* Params */
  double E0;
  double t0;
  double trise;
  double tfall;
  double envmin;
  double phi, theta, chi;
  gsl_complex ex, ey, ez;
  laser_polzn_vector_t e;
  laser_polzn_tensor_t E;
} laser_type1_t;

#endif /* __LASER_TYPE1_H__ */
