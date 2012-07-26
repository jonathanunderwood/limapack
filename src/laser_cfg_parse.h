#ifndef __LASER_CFG_PARSE_H__
#define __LASER_CFG_PARSE_H__

typedef enum _laser_type {
  NONE = 0, /* Use this for error detection. */
  TYPE1 = 1,
} laser_type_t;

/* Structure to keep track of multiple lasers of different types. */
typedef struct _laser_container 
{
  void **laser;
  int nlasers;
} laser_container_t;

laser_container_t * laser_container_ctor (const int nlasers);
void laser_container_dtor (laser_container_t * self);

laser_container_t *laser_container_cfg_parse_ctor (const config_t * cfg);

#endif /* __LASER_CFG_PARSE_H__ */
