#ifndef __ODE_H__
#define __ODE_H__ 1

#include <gsl/gsl_odeiv.h>

typedef struct ode_params
{
  double hinit, hstep;
  double eps_rel, eps_abs, y_scale, dydx_scale;
  gsl_odeiv_system system;
  gsl_odeiv_step *step;
  gsl_odeiv_control *control;
  gsl_odeiv_evolve *evolve;
  gsl_odeiv_step_type *step_type;
  void *params;
} ode_params;

void ode_pfile_memread (const char *buff, ode_params * ode);
void ode_pfile_read (FILE * fnp, ode_params * ode);
void ode_init (ode_params * ode, const int nvar, void *params,
	       int (*function) (double t, const double y[], double dydt[],
				void *params),
	       int (*jacobian) (double t, const double y[], double *dfdy,
				double dfdt[], void *params),
	       const gsl_odeiv_step_type * ode_step_type);
void ode_free (ode_params * ode);
void ode_reset (ode_params * ode);
void ode_step (const double t1, const double t2, double *coef,
	       ode_params * ode);

#endif /* __ODE_H__ */

/* Some sensible defaults for error calculation:
   double eps_rel = 1.0e-6, eps_abs = 1.0e-6;
   double y_scale = 1.0, dydx_scale = 1.0; */
