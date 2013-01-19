#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <libconfig.h>

#include "odesys.h"
#include "memory.h"
#include "au.h"


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
  odesys_expval_t *expval;
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
      expval->data[i] = molecule->expval_ctor(molecule);
      if (expval->data[i] == NULL)
	{
	  fprintf(stderr, 
		  "Failed to allocate storage for expectation values\n");
	  odesys_expval_dtor (ode, expval);
	}
    }

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
      fprintf(stderr, 
	    "Dimension missmatch when adding expection values.\n");
      return -1;
    }

  for (i = 0; i < a->npoints; i++)
    molecule->expval_add_weighted(molecule, 
				  a->data[i], b->data[i], weight);

  return 0;
}

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
odesys_expval_fwrite(const odesys_t *ode, const char *filename)
/* Write out all expectation values to a file in HDF5 format. This
   assumes the file specified by "filename" doesn't already exist. If
   it does, this function will return -1, and should be called again
   with a different filename. */
{
  hid_t file_id, group_id;
  int ret = 0;
  const char *root_group = "/expval";
  char *group_name = NULL;
  molecule_t *mol = ode->params->molecule;
  int i;

  file_id = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

  /* Files to open the file, probably because it already exists. */
  if (file_id < 0)
    return -1;

  group_id = H5Gcreate(file_id, root_group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (group_id < 0) /* Failed to create group. */
    {
      fprintf(stderr, "Failed to create hdf5 group in output file: %s.\n", root_group);
      ret = -1;
      goto exit;
    }

  /* Create a group for each time step, and write out expvals for that
     time step. */
  if (MEMORY_ALLOC_N(group_name, 48) < 0)
    {
      MEMORY_OOMERR;
      ret = -1;
      goto exit;
    }

  for (i = 0; i < ode->npoints; i++)
    {
      hid_t gid;
      int newlength;

      newlength = sprintf(group_name, "%s/%s%d", root_group, "t", i);
      if (newlength > 0)
	{
	  /* The group name for some reason failed to fit in the
	     group_name buffer, so try to repeat with an increased
	     size. If that still fails, give up. */
	  if (MEMORY_REALLOC_N(group_name, newlength) < 0)
	    {
	      MEMORY_OOMERR;
	      ret = -1;
	      goto exit;
	    }

	  if (sprintf(group_name, "%s%s%d", "/expval/", "t", i))
	    {
	      fprintf(stderr, "Failed to create group_name.\n");
	      ret = -1;
	      goto exit;
	    }
	}
      /* Create a group for this time step. */
      gid = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (gid < 0)
	{
	  fprintf(stderr, "Failed to create hdf5 group in output file: %s.\n", group_name);
	  ret = -1;
	  goto exit;
	}

      /* Write out expval data for this time step. */
      mol->expval_fwrite(mol, ode->expval->data[i], gid);
    }

 exit:
  if (group_name != NULL)
    MEMORY_FREE(group_name);

  ret += H5Fclose(file_id);

  return ret;
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

  ode->expval = odesys_expval_ctor(ode);
  if (ode->expval == NULL)
    {
      fprintf(stderr, 
	      "Failed to allocate memory for expectation value storage.\n");
      return -1;
    }

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

  return 0;
}

void
odesys_dtor (odesys_t * ode)
{
  if (ode->expval != NULL)
    odesys_expval_dtor(ode, ode->expval);

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
  odesys_expval_t *buff = NULL;

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
      MEMORY_FREE(coef);
      mol->tdse_worker_dtor(mol, worker);
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

     odesys_reset(odesys);
     mol->get_tdse_job_coef(mol, worker, coef);

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
	     odesys_expval_dtor(odesys, buff);
	     return -1;
	   }
	 
	 /* Calculate expectation values at this time. */
	 mol->expval_calc(mol, coef, t1, buff);

	 /* Propagate to the next time point. */
	 odesys_step(odesys, t1, t1 + odesys->tstep, coef);
	     
	 t1 += odesys->tstep;
       }
     /* Add expectation values for this job to the accumulating
	expectation values, weighted appropriately. */
     weight = mol->get_tdse_job_weight(mol, worker);
     odesys_expval_add_weighted(odesys, odesys->expval, buff, weight);
     
     /* Mark this job as done. */
     mol->set_tdse_job_done(mol, worker);
   }
  
  MEMORY_FREE(coef);
  mol->tdse_worker_dtor(mol, worker);
  odesys_expval_dtor(odesys, buff);

  return 0;
}

#ifdef BUILD_WITH_MPI
#include <mpi.h>

#define NEED_JOB_TAG 0
#define NEW_JOB_TAG 1
#define HAVE_DATA_TAG 2
#define DIE_TAG 99

#define DO_WORK_STATE 0
#define DIE_STATE 1

int 
odesys_tdse_propagate_mpi_master (odesys_t *odesys)
/* MPI based master process function. */
{
  molecule_t *mol = odesys->params->molecule;
  //  odesys_expval_t *expval = NULL;
  odesys_expval_t *buff = NULL;
  molecule_tdse_worker_t *worker = NULL;
  double weight;
  int nprocs, nslaves;

  worker = mol->tdse_worker_ctor(mol);
  if (worker == NULL)
    {
      fprintf(stderr, "Failed to allocate worker.\n");
      return -1;
    }

  /* Allocate storage for storing expectation values summed over all
     initial states, weighted appropriately. */
  /* expval = odesys_expval_ctor(odesys); */
  /* if (expval == NULL) */
  /*   { */
  /*     mol->tdse_worker_dtor(mol, worker); */
  /*     return -1; */
  /*   } */

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
	 mol->expval_mpi_recv(mol, buff, stat.MPI_SOURCE,
			      HAVE_DATA_TAG, MPI_COMM_WORLD);

	 weight = mol->get_tdse_job_weight(mol, worker);
	 odesys_expval_add_weighted(odesys, odesys->expval, buff, weight);

	 mol->set_tdse_job_done(mol, worker);
	 break;
       }
    } while (nslaves > 0);

  odesys_expval_dtor(odesys, buff);

  /* At this point there should be no jobs remaining, and no slave
     processes running. However, we'll do a sanity check. */
  if (mol->get_tdse_job(mol, worker) == 0)
    {
      fprintf(stderr, "All slaves have exited but there are still jobs left to do.\n");
      mol->tdse_worker_dtor(mol, worker);
      return -1;
    }
      
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs > 1)
    {
      fprintf(stderr, "Slave processes still running when they shouldn't be!\n");
      return -1;
    }

  mol->tdse_worker_dtor(mol, worker);
  return 0;
}

int 
odesys_tdse_propagate_mpi_slave (odesys_t *odesys)
/* Simple single threaded generic propagator which uses a single
   worker. */
{
  double *coef = NULL;
  molecule_tdse_worker_t *worker = NULL;
  molecule_t *mol = odesys->params->molecule;
  int state;

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

  /* buff = odesys_expval_ctor(odesys); */
  /* if (buff == NULL) */
  /*   { */
  /*     MEMORY_FREE(coef); */
  /*     mol->tdse_worker_dtor(mol, worker); */
  /*     return -1; */
  /*   } */

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
	  mol->get_tdse_job_coef(mol, worker, coef);
	  odesys_reset(odesys);

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
			  "Populations building up unacceptably in TDSE propagation at time=%g ps. Exiting.\n", 
			  AU_TO_PS(t1));
		  MEMORY_FREE(coef);
		  mol->tdse_worker_dtor(mol, worker);
		  //		  odesys_expval_dtor(odesys, buff);

		  // SEND A MESSAGE TO MASTER SAYING CALCULATION SHOULD DIE.
		  //return -1;
		}
	 
	      /* Calculate expectation values at this time. */
	      mol->expval_calc(mol, coef, t1, odesys->expval);

	      /* Propagate to the next time point. */
	      odesys_step(odesys, t1, t2, coef);
	    }
	  
	  /* Send calculated expectation values back to master. */
	  MPI_Send(NULL, 0, MPI_INT, stat.MPI_SOURCE, HAVE_DATA_TAG,
		   MPI_COMM_WORLD);
	  mol->tdse_worker_mpi_send(mol, worker, stat.MPI_SOURCE, 
				    HAVE_DATA_TAG, MPI_COMM_WORLD);
	  mol->expval_mpi_send(mol, odesys->expval, stat.MPI_SOURCE, 
			       HAVE_DATA_TAG, MPI_COMM_WORLD);
	  break;

	case (DIE_TAG):
	  state = DIE_STATE;
	  break;
	}
    } while (state == DO_WORK_STATE);
	
  /* Cleanup. */
  MEMORY_FREE(coef);
  mol->tdse_worker_dtor(mol, worker);
  //  odesys_expval_dtor(odesys, buff);
      
  return 0;
}

#undef NEED_JOB_TAG 
#undef NEW_JOB_TAG
#undef HAVE_DATA_TAG
#undef DIE_TAG
#undef DO_WORK_STATE
#undef DIE_STATE

#endif /* BUILD_WITH_MPI */
