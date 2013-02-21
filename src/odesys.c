#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <libconfig.h>

#ifdef BUILD_WITH_PTHREADS
#include <pthread.h>
#endif

#include "odesys.h"
#include "molecule.h"
#include "laser.h"
#include "memory.h"
#include "au.h"
#include "slurp.h"

/* Structure to wrap up the molecule and laser parameters to pass to
   the TDSE solver. */
typedef struct _odeparams
{
  molecule_t *molecule;
  laser_collection_t *lasers;
} odeparams_t;

/* This structure is for storing expectation value data for each time
   step. data will be allocated to be an array of pointers to
   molecule_expval_t structures and will have dimension npoints. */
typedef struct _odesys_expval
{
  int npoints;
  molecule_expval_t **data;
#ifdef BUILD_WITH_PTHREADS
  pthread_mutex_t lock;
#endif
} odesys_expval_t;

/* General ODE system container. Contains all parameters, and stuff
   needed for GSL.

   Nb. Some sensible defaults for error calculation: double eps_rel =
   1.0e-6, eps_abs = 1.0e-6; double y_scale = 1.0, dydx_scale =
   1.0. */
struct _odesys
{
  int npoints;
  double tstart, tend, tstep;
  double hstep;//, hinit;
  double eps_rel, eps_abs, y_scale, dydx_scale;
  odeparams_t *params;
  odesys_expval_t *expval;
#ifdef BUILD_WITH_PTHREADS
  int thread_retval;
#endif
};

/* This structure contains all the structures that are needed for the
   GSL ODE propagators. */
typedef struct _odegsl 
{
  gsl_odeiv_system system;
  gsl_odeiv_step *step;
  gsl_odeiv_control *control;
  gsl_odeiv_evolve *evolve;
  gsl_odeiv_step_type *step_type;
  double hstep, hinit;
} odegsl_t;

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

static odegsl_t *
odesys_odegsl_ctor (odesys_t *ode)
{
  odegsl_t *odegsl;
  molecule_t *molecule = ode->params->molecule;
  int ncoef = molecule->get_ncoef(molecule);

  if (MEMORY_ALLOC (odegsl) < 0)
    {
      MEMORY_OOMERR;
      fprintf (stderr, "%s %d: unable to allocate memory for odegsl object.\n",
	       __func__, __LINE__);
      return NULL;
    }

  odegsl->system.dimension = 2 * ncoef;
  odegsl->system.function = odesys_tdse_function;
  odegsl->system.jacobian = NULL;
  odegsl->system.params = (void *) ode->params;

  /* Note that we hard-code the ODE step type here. The step type
     chosen is a general purpose type. In the future, should probably
     try others. */
  odegsl->step = gsl_odeiv_step_alloc (gsl_odeiv_step_rkf45, 2 * ncoef); 
  odegsl->evolve = gsl_odeiv_evolve_alloc (2 * ncoef); 
  odegsl->control = gsl_odeiv_control_standard_new 
    (ode->eps_abs, ode->eps_rel, ode->y_scale, ode->dydx_scale);

  odegsl->hstep = ode->hstep;
  odegsl->hinit = ode->hstep;

  return odegsl;
}

static void
odesys_odegsl_dtor (odegsl_t *odegsl)
{
  MEMORY_FREE (odegsl);
}

static int
odesys_odegsl_step(odegsl_t * odegsl, const double t1, const double t2, 
		   double *coef)
{
  double t = t1;

  while (t < t2)
    {
      int ode_status =
	gsl_odeiv_evolve_apply (odegsl->evolve, odegsl->control, odegsl->step,
				&(odegsl->system), &t, t2, &(odegsl->hstep), coef);
      if (ode_status)
	{
	  fprintf (stderr, "%s %d: gsl_odeiv_evolve_apply error: %s\n",
		   __func__, __LINE__, gsl_strerror (ode_status));
	  return -1;
	}
    }
  return 0;
}

static void
odesys_odegsl_reset (odegsl_t * odegsl)
{
  gsl_odeiv_evolve_reset (odegsl->evolve);
  gsl_odeiv_step_reset (odegsl->step);
  odegsl->hstep = odegsl->hinit;
}

static void
odesys_expval_dtor(const odesys_t *ode, odesys_expval_t *expval)
{
  int i;
  molecule_t * molecule = ode->params->molecule;

  for (i = 0; i < expval->npoints; i++)
    {
      if (expval->data[i] != NULL)
	molecule->expval_dtor(molecule, expval->data[i]);
    }

  if (expval->data != NULL)
    MEMORY_FREE (expval->data);

  MEMORY_FREE (expval);
}

static odesys_expval_t *
odesys_expval_ctor(const odesys_t *ode)
{
  odesys_expval_t * expval;
  int npoints = ode->npoints;
  molecule_t * molecule = ode->params->molecule;
  int i;

  if (MEMORY_ALLOC(expval) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  if (MEMORY_ALLOC_N(expval->data, npoints) < 0)
    {
      MEMORY_OOMERR;
      MEMORY_FREE(expval);
      return NULL;
    }

  for (i = 0; i < npoints; i++)
    {
      expval->data[i] = molecule->expval_ctor (molecule);

      if (expval->data[i] == NULL)
	{
	  fprintf(stderr, 
		  "Failed to allocate storage for expectation values\n");
	  odesys_expval_dtor (ode, expval);
	}
    }
  
#ifdef BUILD_WITH_PTHREADS
  pthread_mutex_init (&(expval->lock), NULL);
#endif

  expval->npoints = npoints;

  return expval;
}

static int
odesys_expval_add_weighted(const odesys_t * ode, odesys_expval_t * a, 
			   odesys_expval_t * b, const double weight)
/* a = a + (weight * b) */
{
  molecule_t * molecule = ode->params->molecule;
  int i;

  if (a->npoints !=b->npoints)
    {
      fprintf (stderr, 
	       "%s %d: dimension missmatch when adding expection values.\n",
	       __func__, __LINE__);
      fprintf (stderr, "a: %d b: %d\n", a->npoints, b->npoints);
      return -1;
    }

  for (i = 0; i < a->npoints; i++)
    molecule->expval_add_weighted 
      (molecule, a->data[i], b->data[i], weight);

  return 0;
}


int
odesys_expval_fwrite(const odesys_t *ode, const char *filename)
/* Write out all expectation values to a file in HDF5 format. This
   assumes the file specified by "filename" doesn't already exist. All
   expectation values are written below a HDF5 group "/expval", with
   each time point being its own group below /expval i.e. /expval/t0,
   expval/t1 etc. The actual time (in ps) of each time point is stored
   as an HDF5 attribute of the /expval/tN group. */
{
  hid_t file_id, root_group_id, dataspace_id;
  int ret = 0;
  const char *root_group = "/expval";
  size_t group_name_size = 5; /* Initial guess. */
  char *group_name = NULL;
  molecule_t *mol = ode->params->molecule;
  int i;

  file_id = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

  /* Failed to open the file, probably because it already exists. */
  if (file_id < 0)
    {
      fprintf(stderr, "%s %d: Failed to open file %s.\n", 
	      __func__, __LINE__, filename);
      return -1;
    }
  
  /* Create the root group /expval. */
  root_group_id = H5Gcreate(file_id, root_group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (root_group_id < 0) /* Failed to create group. */
    {
      fprintf(stderr, "%s %d: Failed to create hdf5 root group in output file: %s.\n", 
	      __func__, __LINE__, filename);
      H5Fclose (file_id);
      return -1;
    }


  /* Create a dataspace that'll be used for the time attribute. */
  dataspace_id = H5Screate(H5S_SCALAR);
  if (dataspace_id < 0)
    {
      fprintf(stderr, "%s %d: Failed to create data space for time attribute.\n",
	      __func__, __LINE__);
      H5Gclose (root_group_id);
      H5Fclose (file_id); 
      return -1;
    }

  /* Allocate a string buffer large enough to hold the group name for
     each time point. In the following loops we take care to grow this buffer if needed. */
  if (MEMORY_ALLOC_N(group_name, group_name_size) < 0) 
    { 
      MEMORY_OOMERR; 
      H5Sclose (dataspace_id);
      H5Gclose (root_group_id);
      H5Fclose (file_id); 
      return -1; 
    }

  /* Create a group for each time step, and write out expvals for that
     time step. */
  for (i = 0; i < ode->npoints; i++)
    {
      hid_t attr_id, group_id;
      int newlength;
      double dt = AU_TO_PS (ode->tstart + i * ode->tstep);

      newlength = snprintf(group_name, group_name_size,
			   "%s%d", "t", i);
      if (newlength > group_name_size)
	{
	  /* The group name for some reason failed to fit in the
	     group_name buffer, so try to repeat with an increased
	     size. If that still fails, give up. */
	  group_name_size = newlength + 1; /* +1 for terminating \0. */
	  if (MEMORY_REALLOC_N(group_name, group_name_size) < 0)
	    {
	      MEMORY_OOMERR;
	      ret = -1;
	      break;
	    }
	  snprintf(group_name, group_name_size, "%s%d", "t", i);
	}

      /* Create a group for this time step. */
      group_id = H5Gcreate(root_group_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (group_id < 0)
	{
	  fprintf(stderr, "%s %d: Failed to create hdf5 group in output file: %s.\n", 
		  __func__, __LINE__, group_name);
	  ret = -1;
	  break;
	}

      /* Store the Delta T value as an attribute of this group. */
      attr_id = H5Acreate (group_id, "time", H5T_NATIVE_DOUBLE, dataspace_id,
			   H5P_DEFAULT, H5P_DEFAULT);
      if (attr_id < 0)
	{
	  fprintf(stderr, "%s %d: Failed to create Delta T attribute.\n",
		  __func__, __LINE__);
	  H5Gclose (group_id);
	  ret= -1;
	  break;
	}

      if (H5Awrite (attr_id, H5T_NATIVE_DOUBLE, &dt) != 0)
	{
	  fprintf(stderr, "%s %d: Failed to write Delta T attribute.\n",
		  __func__, __LINE__);
	  H5Aclose (attr_id);
	  H5Gclose (group_id);
	  ret= -1;
	  break;
	}
      H5Aclose (attr_id);

      /* Write out expval data for this time step. */
      mol->expval_fwrite(mol, ode->expval->data[i], &group_id);

      H5Gclose (group_id);
    }

  H5Sclose (dataspace_id);
  H5Gclose (root_group_id);
  H5Fclose (file_id); 
  
  MEMORY_FREE(group_name);

  return ret;
}

static int
odesys_cfg_parse(odesys_t *ode, const config_t * cfg)
/* Simply parse the config object and return an odesys_t
   object. However, this object is not fully initialized - odesys_init
   would need to be called before use. */
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

static odesys_t *
odesys_ctor()
/* Simply allocate an odesys_t object. Subsequently odesys_cfg_parse and odesys_init
   would need to be called. */
{
  odesys_t * ode;

  if (MEMORY_ALLOC(ode) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  if (MEMORY_ALLOC(ode->params) < 0)
    {
      MEMORY_OOMERR;
      MEMORY_FREE(ode);
      return NULL;
    }

  return ode;
}

static odesys_t *
odesys_parse_from_config_ctor(const config_t * cfg)
{
  odesys_t * ode = odesys_ctor();

  if (ode == NULL)
    return NULL;

  odesys_cfg_parse(ode, cfg);

  return ode;
}


static int
odesys_init (odesys_t * ode, molecule_t * molecule, 
	     laser_collection_t * lasers)
{
  /* Wrap up the molecule and lasers structures into a single
     structure so we can pass a single pointer through the ODE
     functions as required by GSL - avoids global variables. Do this
     before calling odesys_expval_ctor, as that function will call
     molecule->expval_ctor(). */
  ode->params->molecule = molecule;
  ode->params->lasers = lasers;

  ode->expval = odesys_expval_ctor(ode);
  if (ode->expval == NULL)
    {
      fprintf(stderr, 
	      "Failed to allocate memory for expectation value storage.\n");
      return -1;
    }

  return 0;
}

odesys_t *
odesys_parse_cfg_from_buffer_ctor (const char *buffer)
/* Parse the buffer and return an odesys_t object fully initialized
   and ready to use. */
{
  config_t cfg;
  odesys_t *odesys = NULL;
  molecule_t *molecule = NULL; 
  laser_collection_t *lasers = NULL;

  /* Parse configurtion file into a libconfig config_t object. */ 
  config_init (&cfg);

  if (config_read_string (&cfg, buffer) == CONFIG_FALSE)
    {
      switch (config_error_type(&cfg))
	{
	case CONFIG_ERR_PARSE:
	  fprintf(stderr, "Error in configuration at line %d\n",
		  config_error_line(&cfg));
	  fprintf(stderr, "Error was: \"%s\"\n", 
	      config_error_text(&cfg));
	  return NULL; 
	default:
	  fprintf(stderr, "Unknown error in parsing configuration.\n");
	  return NULL;
	}
    }

  /* Parse the config_t object. */
  molecule = molecule_cfg_parse_ctor (&cfg);
  if (molecule == NULL)
    {
      fprintf(stderr, 
	      "Failed to parse molecule information from configuration file.\n");
      config_destroy (&cfg);
      return NULL; 
    }

  lasers = laser_collection_cfg_parse_ctor (&cfg);
  if (lasers == NULL)
    {
      fprintf(stderr, 
	      "Failed to parse laser information from configuration file.\n");
      molecule->dtor(molecule);
      config_destroy (&cfg);
      return NULL;
    }
  
  odesys = odesys_parse_from_config_ctor (&cfg);
  if (odesys == NULL)
    {
      fprintf(stderr, 
	      "Failed to initialize odesys.\n");
      molecule->dtor (molecule);
      laser_collection_dtor (lasers);
      config_destroy (&cfg);
      return NULL;
    }

  config_destroy (&cfg);
  
  odesys_init(odesys, molecule, lasers);

  return odesys;
}

odesys_t *
odesys_parse_cfg_from_file_ctor (const char *file, const int max_cfg_file_size)
/* Parse the file and return an odesys_t object fully initialized
   and ready to use. */
{
  char * buff;
  int buff_size;
  odesys_t *odesys;

  buff_size = slurp_file_to_buffer(file, &buff, max_cfg_file_size);
  
  if (buff_size < 0)
    {
      fprintf(stderr, "Failed to read in config file %s.\n",
	      file);
      return NULL;
    }

  odesys = odesys_parse_cfg_from_buffer_ctor(buff);

  MEMORY_FREE (buff);

  return odesys;
}

void
odesys_dtor (odesys_t * ode)
{
  if (ode->expval != NULL)
    odesys_expval_dtor(ode, ode->expval);

  ode->params->molecule->dtor(ode->params->molecule);
  laser_collection_dtor (ode->params->lasers);
  MEMORY_FREE(ode->params);

  MEMORY_FREE(ode);
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
  odesys_expval_t *buff = NULL;
  odegsl_t *odegsl;
  int retval = 0;

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

  buff = odesys_expval_ctor(odesys);
  if (buff == NULL)
    {
      fprintf(stderr, "%s %d: failed to allocate memory for buff.\n",
	      __func__, __LINE__);
      MEMORY_FREE(coef);
      mol->tdse_worker_dtor(mol, worker);
      return -1;
    }

  odegsl = odesys_odegsl_ctor (odesys);
  if (odegsl == NULL)
    {
      fprintf(stderr, "%s %d: failed to allocate memory for odegsl.\n",
	      __func__, __LINE__);
      odesys_expval_dtor (odesys, buff);
      MEMORY_FREE (coef);
      mol->tdse_worker_dtor (mol, worker);
      return -1;
    }
    
  /* Step through initial states of molecule (eg. the individual
     states in a Boltzmann distribution). For each initial state,
     propagate the TDSE, calculating the expectation values at each
     time-step, and add them to the ensemble averaged expectation
     value, weighted appropriately. */
  while (mol->get_tdse_job(mol, worker) == 0) 
   { 
     double t1 = odesys->tstart;
     int i;

     odesys_odegsl_reset(odegsl);
     mol->get_tdse_job_coef(mol, worker, coef);

     for (i = 0; i < odesys->npoints; i++) /* Step through time points. */
       {
	 /* Check that populations in highest levels aren't growing
	    unacceptably. */
	 if (mol->check_populations(mol, coef) != 0)
	   {
	     fprintf(stderr, 
		     "Populations building up unacceptably in TDSE propagation at t=%g ps. Exiting.\n", 
		     AU_TO_PS(t1));
	     retval = -1;
	     break;
	   }
	 
	 /* Calculate expectation values at this time. */
	 if (mol->expval_calc(mol, coef, t1, buff->data[i]) < 0)
	   {
	     fprintf(stderr, "Error calculating expectation values at t=%g ps. Exiting.\n",
		     AU_TO_PS(t1));
	     retval = -1;
	     break;
	   }

	 /* Propagate to the next time point. */
	 odesys_odegsl_step(odegsl, t1, t1 + odesys->tstep, coef);
	     
	 t1 += odesys->tstep;
       }
     /* Add expectation values for this job to the accumulating
	expectation values, weighted appropriately. */
     weight = mol->get_tdse_job_weight(mol, worker);
     if (odesys_expval_add_weighted(odesys, odesys->expval, buff, weight))
	   {
	     fprintf(stderr, "%s %d: error adding expectation values.\n",
		     __func__, __LINE__);
	     retval = -1;
	     break;
	   }
     
     /* Mark this job as done. */
     mol->set_tdse_job_done(mol, worker);
   }
  
  MEMORY_FREE(coef);
  mol->tdse_worker_dtor(mol, worker);
  odesys_expval_dtor(odesys, buff);
  odesys_odegsl_dtor (odegsl);

  return retval;
}

#ifdef BUILD_WITH_PTHREADS
static odesys_t *
odesys_shallow_copy_ctor (odesys_t * ode)
/* Returns a shallow copy of the input odesys_t object. "Shallow" in
   this context means that all pointers in the output object point to
   the same memory addresses as the input object. In particular, this
   means that the params and expval pointers point to the same memory
   address as the input object. */
{
  odesys_t * copy;

  if (MEMORY_ALLOC (copy) < 0)
    {
      MEMORY_OOMERR;
      fprintf(stderr, "%s %d: Failed to allocate memory for odesys_t copy",
	      __func__, __LINE__);
      return NULL;
    }

  copy->hstep = ode->hstep;
  copy->tstart = ode->tstart;
  copy->tend = ode->tend;
  copy->tstep = ode->tstep;
  copy->npoints = ode->npoints;
  copy->eps_rel = ode->eps_rel;
  copy->eps_abs = ode->eps_abs;
  copy->y_scale = ode->y_scale;
  copy->dydx_scale = ode->dydx_scale;
  copy->params = ode->params;
  copy->expval = ode->expval;

  /* odesys_gsl_init (copy); */

  return copy;
}

static void
odesys_shallow_copy_dtor (odesys_t *copy)
/* Frees memory associated with an odesys_t object created with
   odesys_copy_ctor. This doesn't free the storage associated with the
   params and expval pointers in the odesys_t object. */
{
  MEMORY_FREE (copy);
}

static void *
odesys_tdse_propagate_threaded_child (void *odesystem)
  /* Child thread. TODO: fix up error messages. */
{
  odesys_t *odesys = (odesys_t *) odesystem;
  double *coef = NULL;
  double weight;
  molecule_tdse_worker_t *worker = NULL;
  molecule_t *mol = odesys->params->molecule;
  odesys_expval_t *buff = NULL;
  odegsl_t *odegsl;

  odesys->thread_retval = 0;

  if (MEMORY_ALLOC_N(coef, 2 * mol->get_ncoef(mol)) < 0)
    {
      MEMORY_OOMERR;
      odesys->thread_retval = -1;
      pthread_exit (NULL);
    }

  worker = mol->tdse_worker_ctor(mol);
  if (worker == NULL)
    {
      fprintf(stderr, "Failed to allocate worker.\n");
      MEMORY_FREE(coef);
      odesys->thread_retval = -1;
      pthread_exit (NULL);
    }

  buff = odesys_expval_ctor(odesys);
  if (buff == NULL)
    {
      MEMORY_FREE(coef);
      mol->tdse_worker_dtor(mol, worker);
      odesys->thread_retval = -1;
      pthread_exit (NULL);
    }

  odegsl = odesys_odegsl_ctor (odesys);
  if (odegsl == NULL)
    {
      fprintf(stderr, "%s %d: failed to allocate memory for odegsl.\n",
	      __func__, __LINE__);
      MEMORY_FREE(coef);
      odesys_expval_dtor (odesys, buff);
      mol->tdse_worker_dtor(mol, worker);
      odesys->thread_retval = -1;
      pthread_exit (NULL);
    }

  /* Step through initial states of molecule (eg. the individual
     states in a Boltzmann distribution). For each initial state,
     propagate the TDSE, calculating the expectation values at each
     time-step, and add them to the ensemble averaged expectation
     value, weighted appropriately. */
  do
   { 
     double t1 = odesys->tstart;
     int i, ret;

     /* Get a job to do if any are available. Note that it's important
	here to obtain a mutex so that we don't have a race condition
	in get_tdse_job between looking up a job and assigning a job -
	if two threads access this function at the same time, there's
	a chance they'll both be assigned the same job unless we use a
	mutex. */
     pthread_mutex_lock(&(mol->lock));
     ret = mol->get_tdse_job(mol, worker);
     pthread_mutex_unlock(&(mol->lock));

     if (ret != 0) /* Nothing to do, so exit thread */
	 break;

     odesys_odegsl_reset(odegsl);
     mol->get_tdse_job_coef(mol, worker, coef);

     for (i = 0; i < odesys->npoints; i++) /* Step through time points. */
       {
	 /* Check that populations in highest levels aren't growing
	    unacceptably. */
	 if (mol->check_populations(mol, coef) != 0)
	   {
	     fprintf(stderr, 
		     "Populations building up unacceptably in TDSE propagation at t=%g ps. Exiting.\n", 
		     AU_TO_PS(t1));
	     odesys->thread_retval = -1;
	     break;
	   }
	 
	 /* Calculate expectation values at this time. */
	 if (mol->expval_calc(mol, coef, t1, buff->data[i]) < 0)
	   {
	     fprintf(stderr, "Error calculating expectation values at t=%g ps. Exiting.\n",
		     AU_TO_PS(t1));
	     odesys->thread_retval = -1;
	     break;
	   }

	 /* Propagate to the next time point. */
	 odesys_odegsl_step(odegsl, t1, t1 + odesys->tstep, coef);
	     
	 t1 += odesys->tstep;
       }
     /* Add expectation values for this job to the accumulating
	expectation values, weighted appropriately, and mark this job
	as done. We lock the mol structure here to avoid multiple
	threads writing to the expectation value storeage at the same
	time, though this may not actually be necessary - depends on
	the molecular model. Similarly, locking while updating the
	TDSE job info may not be necessary. */
     weight = mol->get_tdse_job_weight(mol, worker);

     pthread_mutex_lock(&(odesys->expval->lock));
     ret = odesys_expval_add_weighted(odesys, odesys->expval, buff, weight);
     pthread_mutex_unlock(&(odesys->expval->lock));

     if (ret < 0)
       {
	 fprintf(stderr, "%s %d: error adding expectation values.\n",
		 __func__, __LINE__);
	 odesys->thread_retval = -1;
	 break;
       }

     pthread_mutex_lock(&(mol->lock));
     mol->set_tdse_job_done(mol, worker);
     pthread_mutex_unlock(&(mol->lock));
   } while (0);
  
  MEMORY_FREE(coef);
  mol->tdse_worker_dtor(mol, worker);
  odesys_expval_dtor(odesys, buff);
  odesys_odegsl_dtor (odegsl);

  return NULL;
}

int
odesys_tdse_propagate_threaded (odesys_t *odesys, const int nthreads)
/* Simple single threaded generic propagator which uses a single
   worker. */
{
  pthread_t thread[nthreads];
  odesys_t **ode;
  int i, retval = 0;

  /* Each thread process will need a shallow copy of odesys. For
     defintion of shallow copy, see above. */
  if (MEMORY_ALLOC_N(ode, nthreads))
    {
      fprintf(stderr, "%s %d: Failed to allocate memory for ode copy objects.\n",
	      __func__, __LINE__);
      return -1;
    }

  for (i = 0; i < nthreads; i++)
    {
      ode[i] = odesys_shallow_copy_ctor (odesys);
      if (ode[i] == NULL)
	{
	  int j;
	  fprintf(stderr, "%s %d: Failed to create shallow copy of odesys.\n",
		  __func__, __LINE__);

	  for (j = i - 1; j >= 0; j--)
	    odesys_shallow_copy_dtor (ode[i]);

	  MEMORY_FREE (ode);
	  return -1;
	}
    }

  /* Spawn threads. */
  for (i = 0; i < nthreads; i++)
    {
      int ret = pthread_create (&thread[i], NULL, 
				odesys_tdse_propagate_threaded_child,
				(void *) ode[i]);
      if (ret)
	{
	  int j;
	  fprintf(stderr, "%s %d: Failed to create thread.\n", __func__,
		  __LINE__);
	  for (j = i - 1; j >= 0; j--)
	    pthread_cancel(thread[j]);

	  retval = -1;
	}
    }

  /* Wait for threads to finish and check their return status. */
  if (retval == 0)
    {
      for (i = 0; i < nthreads; i++)
	{
	  pthread_join(thread[i], NULL);
	  if (ode[i]->thread_retval < 0)
	    {
	      // TODO: we should get the thread id properly here rather than use i.
	      fprintf (stderr, "%s %d: Thread %d exited uncleanly.\n", __func__, __LINE__, 
		       i);
	      retval = -1;
	    }
	}
    }

  /* Clean up and return. */
  for (i = 0; i < nthreads; i++)
    odesys_shallow_copy_dtor (ode[i]);

  MEMORY_FREE (ode);

  return retval;
}

#endif /* BUILD_WITH_PTHREADS */

#ifdef BUILD_WITH_MPI
#include <mpi.h>

#define NEED_JOB_TAG 0
#define NEW_JOB_TAG 1
#define HAVE_DATA_TAG 2
#define SLAVE_OOM_ERROR_TAG 3
#define SLAVE_CALC_ERROR_TAG 4
#define DIE_TAG 99

#define DO_WORK_STATE 0
#define DIE_STATE 1

int 
odesys_tdse_propagate_mpi_master (odesys_t *odesys)
/* MPI based master process function. This function dispatches tDSE
   jobs to MPI slaves, and takes care of collating the returned
   expectation values. A negative return value indicates a problem,
   and the caller should call MPI_Abort to clean up. The function will
   return 0 otherwise. This function won't call MPI_Abort() or
   MPI_Finalize. */
{
  molecule_t *mol = odesys->params->molecule;
  odesys_expval_t *buff = NULL;
  molecule_tdse_worker_t *worker = NULL;
  double weight;
  int nprocs, nslaves, rank;
  int max_host_length = ODESYS_MAX_HOST_NAME_LENGTH;
  char host[max_host_length];
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(host, &max_host_length);
  if (rank != 0)
    {
      fprintf (stderr, "%s not running as master process.\n", __func__);
      return (-1);
    }

  worker = mol->tdse_worker_ctor(mol);
  if (worker == NULL)
    {
      fprintf(stderr, "Failed to allocate tdse_worker.\n");
      return -1;
    }

  /* Allocate storage for storing expectation values corresponding to
     individual initial states. */
  buff = odesys_expval_ctor(odesys);
  if (buff == NULL)
    {
      mol->tdse_worker_dtor(mol, worker);
      return -1;
    }

  /* Find the number of slave processes running. Our strategy will be
     to loop, dispatching jobs to slave processes as they request
     jobs. Once all jobs have been allocated we'll send a die message
     to any further requests for jobs. Once the number of slave
     processes is zero, we can exit the loop. */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  nslaves = nprocs - 1; /* Since nprocs includes the master process. */

  /* Step through initial states of molecule (eg. the individual
     states in a Boltzmann distribution). For each initial state,
     propagate the TDSE, calculating the expectation values at each
     time-step, and add them to the ensemble averaged expectation
     value, weighted appropriately. */
  do
    {
      MPI_Status stat;
    
     /* Receive an empty message - we'll decide what to do next based
	purely on the message tag. */
     MPI_Recv (NULL, 0, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	       MPI_COMM_WORLD, &stat);

     switch (stat.MPI_TAG)
       {
	 int i;
       case NEED_JOB_TAG:
	 /* Slave is waiting for work so get a new job and send
	    that. If there are no jobs left, tell the slave to
	    shutdown. */
	 if (mol->get_tdse_job(mol, worker) == 0)
	   {
	     MPI_Send(NULL, 0, MPI_INT, stat.MPI_SOURCE, 
		      NEW_JOB_TAG, MPI_COMM_WORLD);
	     mol->tdse_worker_mpi_send(mol, worker, stat.MPI_SOURCE, 
				       NEW_JOB_TAG, MPI_COMM_WORLD);
	     fprintf(stdout, "<%d::%s> Job sent to slave process %d => %s\n",
		     rank, host, stat.MPI_SOURCE, worker->description);
	   }
	 else 
	   {
	     /* No more jobs left to do, so send a terminate message. */
	     MPI_Send(NULL, 0, MPI_INT, stat.MPI_SOURCE, 
		      DIE_TAG, MPI_COMM_WORLD);
	     /* Reduce our tally of number of slaves accordingly. */
	     nslaves--; 
	   }
	 break;
       case HAVE_DATA_TAG:
	 /* Slave is returning expectation value data for a job. So we
	    need to receive (i) the detail of what job the slave has
	    done; (ii) the corresponding expectation values. Once this
	    data is received we need to add it to the accumulating
	    expectation values weighted appropriately. */
	 mol->tdse_worker_mpi_recv(mol, worker, stat.MPI_SOURCE,
				   HAVE_DATA_TAG, MPI_COMM_WORLD);

	 for (i = 0; i < odesys->npoints; i++)
	   mol->expval_mpi_recv(mol, buff->data[i], stat.MPI_SOURCE,
				HAVE_DATA_TAG, MPI_COMM_WORLD);

	 weight = mol->get_tdse_job_weight(mol, worker);
	 
	 if (odesys_expval_add_weighted(odesys, odesys->expval, buff, weight))
	   {
	     fprintf(stderr, "%s %d: error adding expectation values.\n",
		     __func__, __LINE__);
	     mol->tdse_worker_dtor(mol, worker);
	     odesys_expval_dtor(odesys, buff);
	     return -1;
	   }

	 mol->set_tdse_job_done(mol, worker);
	 break;
       case SLAVE_OOM_ERROR_TAG:
	 /* A slave process wasn't able to allocate memory it
	   needed. We'll let this slave process die, but won't kill
	   off the calculation as other slaves may be functioning
	   well. Since all memory allocation in the slave is done
	   before starting jobs, no job has been lost at this point,
	   we simply decrease nslaves by 1. */
	 nslaves--;
	 fprintf(stderr, 
		 "Slave process %d couldn't allocate sufficient memory and quit.\n",
		 stat.MPI_SOURCE);
	 break;
       case SLAVE_CALC_ERROR_TAG:
	 /* A slave process encountered an error in the calculation
	    such that it's not worth continuing the calculation at
	    all. */
	 fprintf(stderr, 
		 "Slave process %d hit a fatal calculation error.\n",
		 stat.MPI_SOURCE);
	 mol->tdse_worker_dtor(mol, worker);
	 odesys_expval_dtor(odesys, buff);
	 return -1;
	 break;
       }

    } while (nslaves > 0);

  odesys_expval_dtor(odesys, buff);

  /* At this point there should be no jobs remaining. However, we'll
     do a sanity check. */
  if (mol->get_tdse_job(mol, worker) == 0)
    {
      fprintf(stderr, "All slaves have exited but there are still jobs left to do.\n");
      mol->tdse_worker_dtor(mol, worker);
      return -1;
    }
      
  mol->tdse_worker_dtor(mol, worker);
  return 0;
}

int 
odesys_tdse_propagate_mpi_slave (odesys_t *odesys)
/* MPI based slave process function. A negative return value
   indicates a problem, and the caller should call MPI_Abort to clean
   up. The function will return 0 otherwise. This function won't call
   MPI_Abort() or MPI_Finalize. */
{
  double *coef = NULL;
  molecule_tdse_worker_t *worker = NULL;
  odegsl_t *odegsl;
  molecule_t *mol = odesys->params->molecule;
  int state;
  int max_host_length = ODESYS_MAX_HOST_NAME_LENGTH;
  char host[max_host_length];
  int rank;

  MPI_Get_processor_name(host, &max_host_length);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (MEMORY_ALLOC_N(coef, 2 * mol->get_ncoef(mol)) < 0)
    {
      MEMORY_OOMERR;
      fprintf(stderr, "Slave process %d on host %s failed to allocate memory for coef.\n",
	      rank, host);
      MPI_Send(NULL, 0, MPI_INT, 0, SLAVE_OOM_ERROR_TAG, MPI_COMM_WORLD);
      return -1;
    }

  worker = mol->tdse_worker_ctor(mol);
  if (worker == NULL)
    {
      fprintf(stderr, "<%d::%s> Failed to allocate memory for worker.\n",
	      rank, host);
      MPI_Send(NULL, 0, MPI_INT, 0, SLAVE_OOM_ERROR_TAG, MPI_COMM_WORLD);
      MEMORY_FREE(coef);
      return -1;
    }

  odegsl = odesys_odegsl_ctor (odesys);
  if (odegsl == NULL)
    {
      fprintf(stderr, "<%d::%s> Failed to allocate memory for odegsl.\n",
	      rank, host);
      MPI_Send(NULL, 0, MPI_INT, 0, SLAVE_OOM_ERROR_TAG, MPI_COMM_WORLD);
      MEMORY_FREE(coef);
      mol->tdse_worker_dtor(mol, worker);
      return -1;
    }

  state = DO_WORK_STATE;

  do
    {
      MPI_Status stat;

      /* Request job. */
      MPI_Send(NULL, 0, MPI_INT, 0, NEED_JOB_TAG, MPI_COMM_WORLD);
      MPI_Recv(NULL, 0, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, 
	       &stat);

      /* Decide what to do based on received message tag. */
      switch (stat.MPI_TAG)
	{
	  int i;

	case (NEW_JOB_TAG):
	  mol->tdse_worker_mpi_recv(mol, worker, stat.MPI_SOURCE, 
				    NEW_JOB_TAG, MPI_COMM_WORLD);
	  fprintf(stdout, "<%d::%s> Job received => %s\n",
		  rank, host, worker->description);

	  mol->get_tdse_job_coef(mol, worker, coef);
	  odesys_odegsl_reset(odegsl);

	  /* Step through time points. */
	  for (i = 0; i < odesys->npoints; i++)
	    {
	      double t1 = odesys->tstart + i * odesys->tstep;
	      double t2 = t1 + odesys->tstep;

	      /* Check that populations in highest levels aren't growing
		 unacceptably. */
	      if (mol->check_populations(mol, coef) != 0)
		{
		  fprintf(stderr, 
			  "<%d::%s> Populations building up unacceptably in TDSE propagation at time=%g ps.\n", 
			  rank, host, AU_TO_PS(t1));
		  MEMORY_FREE(coef);
		  mol->tdse_worker_dtor(mol, worker);
		  odesys_odegsl_dtor (odegsl);
		  MPI_Send(NULL, 0, MPI_INT, 0, SLAVE_CALC_ERROR_TAG, MPI_COMM_WORLD);
		  return -1;
		}
	 
	      /* Calculate expectation values at this time. */
	      if (mol->expval_calc(mol, coef, t1, odesys->expval->data[i]) < 0)
		{
		  fprintf(stderr, "<%d::%s> Error calculating expectation values at t= %g ps.\n",
			  rank, host, AU_TO_PS(t1));
		  MEMORY_FREE(coef);
		  mol->tdse_worker_dtor(mol, worker);
		  odesys_odegsl_dtor (odegsl);
		  MPI_Send(NULL, 0, MPI_INT, 0, SLAVE_CALC_ERROR_TAG, MPI_COMM_WORLD);
		  return -1;
		}

	      /* Propagate to the next time point. */
	      odesys_odegsl_step(odegsl, t1, t2, coef);
	    }
	  
	  /* Send calculated expectation values back to master. */
	  MPI_Send(NULL, 0, MPI_INT, stat.MPI_SOURCE, HAVE_DATA_TAG,
		   MPI_COMM_WORLD);
	  mol->tdse_worker_mpi_send(mol, worker, stat.MPI_SOURCE, 
				    HAVE_DATA_TAG, MPI_COMM_WORLD);
	  for (i = 0; i < odesys->npoints; i++)
	    mol->expval_mpi_send(mol, odesys->expval->data[i], stat.MPI_SOURCE, 
				 HAVE_DATA_TAG, MPI_COMM_WORLD);
	  break;

	case (DIE_TAG):
	  state = DIE_STATE;
	  break;
	}
    } while (state == DO_WORK_STATE);
	
  /* Cleanup. */
  MEMORY_FREE(coef);
  odesys_odegsl_dtor (odegsl);
  mol->tdse_worker_dtor(mol, worker);
      
  return 0;
}

#undef NEED_JOB_TAG 
#undef NEW_JOB_TAG
#undef HAVE_DATA_TAG
#undef SLAVE_OOM_ERROR_TAG
#undef SLAVE_CALC_ERROR_TAG
#undef DIE_TAG

#undef DO_WORK_STATE
#undef DIE_STATE

#endif /* BUILD_WITH_MPI */
