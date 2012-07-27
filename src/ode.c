#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include "slurp.h"
#include "ode.h"

#define BUFF_LENGTH 30

void
ode_pfile_memread (const char *buff, ode_params * ode)
{
  ode->eps_rel = dgrabp (buff, NULL, "eps_rel");
  ode->eps_abs = dgrabp (buff, NULL, "eps_abs");
  ode->y_scale = dgrabp (buff, NULL, "y_scale");
  ode->dydx_scale = dgrabp (buff, NULL, "dydx_scale");
  ode->hstep = dgrabp (buff, NULL, "hstep");
}

void
ode_pfile_read (FILE * fnp, ode_params * ode)
{
  char param_name[BUFF_LENGTH];
  rewind (fnp);

  while (fscanf (fnp, "%s", param_name) != EOF)
    {
      if (!strcmp (param_name, "eps_rel"))
	fscanf (fnp, "%lf", &ode->eps_rel);
      else if (!strcmp (param_name, "eps_abs"))
	fscanf (fnp, "%lf", &ode->eps_abs);
      else if (!strcmp (param_name, "y_scale"))
	fscanf (fnp, "%lf", &ode->y_scale);
      else if (!strcmp (param_name, "dydx_scale"))
	fscanf (fnp, "%lf", &ode->dydx_scale);
      else if (!strcmp (param_name, "hstep"))
	fscanf (fnp, "%lf", &ode->hstep);
      else			/* Ignore all other lines. */
	while ((getc (fnp)) != '\n');
    }
}

void
ode_init (ode_params * ode, const int nvar, void *params,
	  int (*function) (double t, const double y[], double dydt[],
			   void *params),
	  int (*jacobian) (double t, const double y[], double *dfdy,
			   double dfdt[], void *params),
	  const gsl_odeiv_step_type * ode_step_type)
{
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
}

void
ode_free (ode_params * ode)
{
  gsl_odeiv_evolve_free (ode->evolve);
  gsl_odeiv_control_free (ode->control);
  gsl_odeiv_step_free (ode->step);
}

void
ode_reset (ode_params * ode)
{
  gsl_odeiv_evolve_reset (ode->evolve);
  gsl_odeiv_step_reset (ode->step);
  ode->hstep = ode->hinit;
}

void
ode_step (const double t1, const double t2, double *coef, ode_params * ode)
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
	  exit (1);
	}
    }
}

#undef BUFF_LENGTH
