#ifndef __ODE_H__
#define __ODE_H__

#include <gsl/gsl_odeiv.h>

typedef struct _odesys
{
  double hinit, hstep;
  double eps_rel, eps_abs, y_scale, dydx_scale;
  gsl_odeiv_system system;
  gsl_odeiv_step *step;
  gsl_odeiv_control *control;
  gsl_odeiv_evolve *evolve;
  gsl_odeiv_step_type *step_type;
  void *params;
} odesys_t;


odesys_t odesys_ctor (const int nvar, void *params,
		      int (*function) (double t, const double y[], 
				       double dydt[], void *params),
		      int (*jacobian) (double t, const double y[], 
				       double *dfdy,double dfdt[], 
				       void *params),
		      const gsl_odeiv_step_type * ode_step_type);
void odesys_dtor (odesys_t * ode);
void odesys_reset (odesys_t * ode);
int odesys_step (odesys_t * ode, const double t1, const double t2, double *coef);
int ode_cfg_parse(odesys_t * ode, const config_t * cfg);

#endif /* __ODE_H__ */

/* Some sensible defaults for error calculation:
   double eps_rel = 1.0e-6, eps_abs = 1.0e-6;
   double y_scale = 1.0, dydx_scale = 1.0; */
