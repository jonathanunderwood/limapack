#ifndef __ODESYS_H__
#define __ODESYS_H__

/* The intention of this module is to hide as much as possible the
   details of the ODE solver - in the future may use something other
   than GSL, but the API should remain the more or less the same. */

#include <gsl/gsl_odeiv.h>

/* The integration function should return ODESYS_SUCCESS. */
#include <gsl/gsl_errno.h>
#define ODESYS_SUCCESS GSL_SUCCESS

typedef struct _odesys
{
  double hstep, hinit;
  double eps_rel, eps_abs, y_scale, dydx_scale;
  gsl_odeiv_system system;
  gsl_odeiv_step *step;
  gsl_odeiv_control *control;
  gsl_odeiv_evolve *evolve;
  gsl_odeiv_step_type *step_type;
  void *params;
} odesys_t;


int odesys_init (odesys_t * ode, const int nvar, void *params,
		 int (*function) (double t, const double y[], 
				  double dydt[], void *params)
		 );


odesys_t * odesys_ctor ();
void odesys_dtor (odesys_t * ode);
void odesys_reset (odesys_t * ode);
int odesys_step (odesys_t * ode, const double t1, const double t2, double *coef);
int odesys_cfg_parse(odesys_t * ode, const config_t * cfg);

#endif /* __ODESYS_H__ */

/* Some sensible defaults for error calculation:
   double eps_rel = 1.0e-6, eps_abs = 1.0e-6;
   double y_scale = 1.0, dydx_scale = 1.0; */
