#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <libconfig.h>

#include "odesys.h"
#include "memory.h"
#include "au.h"

static int
odesys_tdse_function (const double t, const double coefs[], 
		      double derivs[], void *params)
/* Generic wrapper function for the RHS of the TDSE. Performs type
   castings, and dispatches the molecule specific function. */
{
  odeparams_t *p = (odeparams_t *) params;
  laser_collection_t *lasers = p->lasers;
  molecule_t *mol = p->molecule;
  int i;

  /* Ensure all elements of derivative array are zero. This can't be
     done in the main loop of mol->tdse_rhs. */
  for (i = 0; i < 2 * mol->get_ncoef(mol); i++)
    derivs[i] = 0.0;

  /* Dispatch molecule specific function and return. */
  return mol->tdse_rhs(mol, lasers, t, coefs, derivs);
}

int
odesys_cfg_parse(odesys_t *ode, const config_t * cfg)
{
  double hstep, tstart, tend;
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
	config_setting_lookup_float (s, "dydx_scale", &(ode->dydx_scale)) &&
	config_setting_lookup_float (s, "tstart", &tstart) &&
	config_setting_lookup_float (s, "tend", &tend) &&
	config_setting_lookup_int (s, "npoints", &ode->npoints)
	))
    {
      fprintf(stderr, "Incomplete/malformed odesolver configuration in file.\n"); 
      return -1;
    }

  /* Convert to atomic units. */
  ode->hstep = PS_TO_AU (hstep);
  ode->tstart = PS_TO_AU (tstart);
  ode->tend = PS_TO_AU (tend);
  ode->tstep = (ode->tend - ode->tstart) / (ode->npoints - 1);

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

  if (MEMORY_ALLOC(ode->params)  < 0)
    {
      MEMORY_OOMERR;
      MEMORY_FREE(ode);
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

odesys_t *
odesys_cfg_parse_ctor(const config_t * cfg)
{
  odesys_t * ode = odesys_ctor();

  if (ode == NULL)
    return NULL;

  odesys_cfg_parse(ode, cfg);

  return ode;
}

int
odesys_init (odesys_t * ode, molecule_t * molecule, 
	     laser_collection_t *lasers)
{
  int ncoef = molecule->get_ncoef(molecule);
  int nexpval = molecule->get_nexpval(molecule);
  int i;

  /* See odesys_ctor - this checks to see if ode_cfg_parse has been
     called to set up the basic parameters before odesys_init is
     called. */
  if (ode->hstep < 0)
    {
      fprintf(stderr, "odesys_cfg_parse has not been called, hstep < 0.\n");
      return -1;
    }
  
  ode->system.dimension = 2 * ncoef;
  ode->system.function = odesys_tdse_function;
  ode->system.jacobian = NULL;

  /* Wrap up the molecule and lasers structures into a single
     structure so we can pass a single pointer through the ODE
     functions as required by GSL - avoids global variables. */
  ode->params->molecule = molecule;
  ode->params->lasers = lasers;
  ode->system.params = (void *)ode->params;

  /* Note that we hard-code the ODE step type here. The step type
     chosen is a general purpose type. In the future, should probably
     try others. */
  ode->step = gsl_odeiv_step_alloc (gsl_odeiv_step_rkf45, 2 * ncoef); 
  ode->evolve = gsl_odeiv_evolve_alloc (2 * ncoef); 
  ode->control = gsl_odeiv_control_standard_new 
    (ode->eps_abs, ode->eps_rel, ode->y_scale, ode->dydx_scale);

  /* This allows us to reset hstep if we need to. */
  ode->hinit = ode->hstep;

  /* Initialize storge for expectation values. This is done as a
     contiguous memory block intentionally. */
  if (MEMORY_ALLOC_N(ode->expval, ode->npoints) < 0)
    {
      MEMORY_OOMERR;
      //      odesys_dtor();
      return -1;
    }

  if (MEMORY_ALLOC_N((ode->expval)[0], ode->npoints * nexpval) < 0)
    {
      MEMORY_OOMERR;
      MEMORY_FREE(ode->expval);
      return -1;
    }

  for (i = 1; i < ode->npoints; i++)
    (ode->expval)[i] = (ode->expval)[i - 1] + nexpval;
 
  return 0;
}

void
odesys_dtor (odesys_t * ode)
{
  if (ode->expval[0] != NULL)
    MEMORY_FREE(ode->expval[0]);

  if (ode->expval != NULL)
    MEMORY_FREE(ode->expval);

  gsl_odeiv_evolve_free (ode->evolve);
  gsl_odeiv_control_free (ode->control);
  gsl_odeiv_step_free (ode->step);
  MEMORY_FREE(ode->params);
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

int 
odesys_tdse_propagate_simple (odesys_t *odesys)
/* Simple single threaded generic propagator which uses a single
   worker. */
{
  double *coef = NULL;
  double weight;
  molecule_tdse_worker_t *worker = NULL;
  molecule_t *mol = odesys->params->molecule;
  //  laser_collection_t *las = odesys->params->lasers;

  if (MEMORY_ALLOC_N(coef, 2 * mol->get_ncoef(mol)) < 0)
    {
      MEMORY_OOMERR;
      return -1;
    }

  worker = mol->tdse_worker_ctor(mol);
  if (worker == NULL)
    {
      fprintf(stderr, "Failed to allocate worker.\n");
      MEMORY_FREE(coef);
      return -1;
    }

  /* Step through initial states of molecule (eg. the individual
     states in a Boltzmann distribution). For each initial state,
     propagate the TDSE, calculating the expectation values at each
     time-step, and add them to the ensemble averaged expectation
     value, weighted appropriately. */
  while (mol->get_tdse_job(mol, worker, coef, &weight) == 0) 
   { 
     double t1 = odesys->tstart;
     int i;

     odesys_reset(odesys);
     
     for (i = 0; i < odesys->npoints; i++) /* Step through time points. */
       {
	 /* Check that populations in highest levels aren't growing
	    unacceptably. */
	 if (mol->check_populations(mol, coef) != 0)
	   {
	     fprintf(stderr, 
		     "Populations building up unacceptably in TDSE propagation at time=%g ps. Exiting.\n", 
		     AU_TO_PS(t1));
	     MEMORY_FREE(coef);
	     mol->tdse_worker_dtor(mol, worker);
	     return -1;
	   }
	 
	 /* Calculate expectation values at this time. */
	 //some_buffer = mol->calc_exp_values (mol, coef);
	 // add to tally, suitably weighted

	 /* Propagate to the next time point. We use the ODE
	    propagator only if one of the laser fields is
	    non-negligible. Since we're propagating in the interaction
	    picture, if all the laser fields are negligible, the
	    coefficients are unchanged. In this case, we'll need to
	    reset the ODE propagator the next time the lasers are non
	    negligible. */
	 /* if (t1 < odesys->tend) */
	 /*   { */
	 /*     if (laser_collection_all_negligible(las, t1) &&  */
	 /* 	 laser_collection_all_negligible(las, t2)) */
	 /*       need_ode_reset = 1; */
	 /*     else */
	 /*       { */
	 /* 	 if (need_ode_reset) */
	 /* 	   { */
	 /* 	     odesys_reset(odesys); */
	 /* 	     need_ode_reset = 0; */
	 /* 	   } */
		 
	 odesys_step(odesys, t1, t1 + odesys->tstep, coef);
       /* } */
	     
	 t1 += odesys->tstep;
	   /* } */
       }
     mol->set_tdse_job_done(mol, worker);
   }
  
  MEMORY_FREE(coef);
  mol->tdse_worker_dtor(mol, worker);
  return 0;
}
