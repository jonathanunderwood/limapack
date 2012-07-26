#include <stdlib.h>
#include <stdio.h>
#include <libconfig.h>

#include "laser.h"
#include "laser_cfg_parse.h"

int
main ()
{
  config_t cfg;
  laser_container_t lasers;
  int ret;

  config_init (&cfg);

  if(! config_read_file(&cfg, "test_laser.cfg"))
  {
    fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
            config_error_line(&cfg), config_error_text(&cfg));
    config_destroy(&cfg);
    exit(EXIT_FAILURE);
  }

  ret = laser_cfg_parse(&cfg, &lasers);

  return ret;
}
