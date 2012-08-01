#include <stdlib.h>
#include <stdio.h>
#include <libconfig.h>

#include "odesys.h"

int
func (double t, const double y[], double f[],
      void *params)
{
  double mu = *(double *)params;
  f[0] = y[1];
  f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
  return ODESYS_SUCCESS;
}

int
main ()
{
  config_t cfg;
  odesys_t *ode;
  double mu = 10;
  double t1 = 0.0, t2 = 100.0;
  double y[2] = { 1.0, 0.0 };
     
  config_init (&cfg);

  if(! config_read_file(&cfg, "test.cfg"))
  {
    fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
            config_error_line(&cfg), config_error_text(&cfg));
    config_destroy(&cfg);
    exit(EXIT_FAILURE);
  }

  ode = odesys_ctor();
  if (ode == NULL)
    {
      fprintf(stderr, "Failed to allocate odesys.\n");
      config_destroy(&cfg);
      exit (EXIT_FAILURE);
    }

  if (odesys_cfg_parse(ode, &cfg))
    {
      fprintf(stderr, "Failed to parse ode configuration.\n");
      config_destroy(&cfg);
      odesys_dtor(ode);
      exit (EXIT_FAILURE);
    }

  if (odesys_init(ode, 2, &mu, func))
    {
      fprintf(stderr, "odesys_init failed.\n");
      exit (EXIT_FAILURE);
    }

  if (odesys_step(ode, t1, t2, y))
    {
      fprintf(stderr, "odesys_step failed.\n");
      config_destroy(&cfg);
      odesys_dtor(ode);
      exit (EXIT_FAILURE);
    }

  /* fprintf(stdout, "y[0]: %f y[1]: %f\n", y[0], y[1]); */
  config_destroy(&cfg);
  odesys_dtor(ode);

  return 0;
}
