#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <libconfig.h>

#include "ode.h"
#include "memory.h"

int
ode_cfg_parse(odesys_t *ode, const config_t * cfg)
{

  config_setting_t *s;

  s = config_lookup(cfg, "odesolver");

  if (s == NULL)
    {
      fprintf(stderr, "Failed to find odesolver section in config.\n");
      return -1;
    }

  if (!(config_setting_lookup_float (s, "hstep", &(ode->hstep)) &&
	config_setting_lookup_float (s, "eps_rel", &(ode->eps_rel)) &&
	config_setting_lookup_float (s, "eps_abs", &(ode->eps_abs)) &&
	config_setting_lookup_float (s, "yscale", &(ode->y_scale)) &&
	config_setting_lookup_float (s, "dydx_scale", &(ode->dydx_scale))
	))
    {
      fprintf(stderr, "Incomplete/malformed odesolver configuration  in file.\n"); 
      return -1;
    }

  return 0;
}

odesys_t *
ode_ctor (const int nvar, void *params,
	  int (*function) (double t, const double y[], double dydt[],
			   void *params),
	  int (*jacobian) (double t, const double y[], double *dfdy,
			   double dfdt[], void *params),
	  const gsl_odeiv_step_type * ode_step_type)
{
  odesys_t * ode;

  if (MEMORY_ALLOC(ode) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  ode->system.dimension = nvar;
  ode->system.function = function;
  ode->system.jacobian = jacobian;
  ode->system.params = params;
  ode->step = gsl_odeiv_step_alloc (ode_step_type, nvar);
  ode->evolve = gsl_odeiv_evolve_alloc (nvar);
  ode->control =
    gsl_odeiv_control_standard_new (ode->eps_abs, ode->eps_rel, ode->y_scale,
				    ode->dydx_scale);
  ode->hinit = ode->hstep;

  return ode;
}

void
ode_dtor (odesys_t * ode)
{
  gsl_odeiv_evolve_free (ode->evolve);
  gsl_odeiv_control_free (ode->control);
  gsl_odeiv_step_free (ode->step);
  MEMORY_FREE(ode);
}

void
ode_reset (odesys_t * ode)
{
  gsl_odeiv_evolve_reset (ode->evolve);
  gsl_odeiv_step_reset (ode->step);
  ode->hstep = ode->hinit;
}

int
ode_step (odesys_t * ode, const double t1, const double t2, double *coef)
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

#undef BUFF_LENGTH
