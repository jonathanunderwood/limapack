#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <libconfig.h>

#include "odesys.h"
#include "memory.h"
#include "au.h"

int
odesys_cfg_parse(odesys_t *ode, const config_t * cfg)
{
  double hstep;
  config_setting_t *s;
  
  s = config_lookup(cfg, "odesolver");

  if (s == NULL)
    {
      fprintf(stderr, "Failed to find odesolver section in config.\n");
      return -1;
    }

  if (!(config_setting_lookup_float (s, "hstep", &hstep) &&
	config_setting_lookup_float (s, "eps_rel", &(ode->eps_rel)) &&
	config_setting_lookup_float (s, "eps_abs", &(ode->eps_abs)) &&
	config_setting_lookup_float (s, "y_scale", &(ode->y_scale)) &&
	config_setting_lookup_float (s, "dydx_scale", &(ode->dydx_scale))
	))
    {
      fprintf(stderr, "Incomplete/malformed odesolver configuration in file.\n"); 
      return -1;
    }

  /* Convert hstep to atomic units. */
  hstep = PS_TO_AU (hstep);
  ode->hstep = hstep;

  return 0;
}

odesys_t *
odesys_ctor()
{
  odesys_t * ode;

  if (MEMORY_ALLOC(ode) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }
 
  /* Set this to a negative vaule so that odesys_init can check that
     odesys_cfg_parse has been called. */
  ode->hstep = -99.0;

  ode->system.dimension = 0;
  ode->system.function = NULL;
  ode->system.jacobian = NULL;
  ode->system.params = NULL;
  ode->step = NULL;
  ode->evolve = NULL;
  ode->control = NULL;

  return ode;
}

int
odesys_init (odesys_t * ode, const int nvar, void *params,
	     int (*function) (double t, const double y[], double dydt[],
			      void *params)
	     )
{
  /* See odesys_ctor - this checks to see if ode_cfg_parse has been
     called to set up the basic parameters before odesys_init is
     called. */
  if (ode->hstep < 0)
    {
      fprintf(stderr, "odesys_cfg_parse has not been called, hstep < 0.\n");
      return -1;
    }
  
  ode->system.dimension = nvar;
  ode->system.function = function;
  ode->system.jacobian = NULL;
  ode->system.params = params;

  /* Note that we hard-code the ODE step type here. The step type
     chosen is a general purpose type. In the future, should probably
     try others. */
  ode->step = gsl_odeiv_step_alloc
  (gsl_odeiv_step_rkf45, nvar); ode->evolve = gsl_odeiv_evolve_alloc
  (nvar); ode->control = gsl_odeiv_control_standard_new (ode->eps_abs,
  ode->eps_rel, ode->y_scale, ode->dydx_scale);

  /* This allows us to reset hstep if we need to. */
  ode->hinit = ode->hstep;
  
  return 0;
}

void
odesys_dtor (odesys_t * ode)
{
  gsl_odeiv_evolve_free (ode->evolve);
  gsl_odeiv_control_free (ode->control);
  gsl_odeiv_step_free (ode->step);
  MEMORY_FREE(ode);
}

void
odesys_reset (odesys_t * ode)
{
  gsl_odeiv_evolve_reset (ode->evolve);
  gsl_odeiv_step_reset (ode->step);
  ode->hstep = ode->hinit;
}

int
odesys_step (odesys_t * ode, const double t1, const double t2, double *coef)
{
  double t = t1;

  while (t < t2)
    {
      int ode_status =
	gsl_odeiv_evolve_apply (ode->evolve, ode->control, ode->step,
				&(ode->system), &t, t2, &(ode->hstep), coef);
      if (ode_status)
	{
	  fprintf (stderr, "gsl_odeiv_evolve_apply error: %s\n",
		   gsl_strerror (ode_status));
	  return -1;
	}
    }
  return 0;
}


