#ifndef __LASER_TYPE1_H__
#define __LASER_TYPE1_H__

#include <libconfig.h>
#include "laser.h"

typedef struct _laser_type1
{
  laser_t parent;
  /* Params */
  double freq;
  double E0;
  double t0;
  double trise;
  double tfall;
  double envmin;
  double phi, theta, chi;
  gsl_complex ex, ey, ez;
  laser_polzn_vector_t *e;
} laser_type1_t;

laser_t * laser_type1_cfg_parse_ctor (config_setting_t * element);
//void laser_type1_dtor(laser_t *laser);

#endif /* __LASER_TYPE1_H__ */
