#ifndef __LASER_TYPE1_H__
#define __LASER_TYPE1_H__

#include "laser.h"

typedef struct _laser_type1
{
  laser_t parent;
  /* Repeat these callbacks here (they're already part of the laser_t) so we can
     call them directly without a recast. */
  //  laser_polzn_vector * (*get_polzn_vector) (struct _laser *self, const double t);
  //gsl_comlplex * (*get_envelope) (struct _laser *self, const double t);
  /* Params */
  double freq;
  double E0;
  double t0;
  double trise;
  double tfall;
  double envmin;
  gsl_complex ex, ey, ez;
  double phi, theta, chi;
  laser_polzn_vector_t *e;
} laser_type1_t;

laser_type1_t * laser_type1_ctor();


#endif /* __LASER_TYPE1_H__ */
