#include <stdlib.h>
#include <stdio.h>
#include <libconfig.h>
#include <string.h>

#include "molecule.h"
#include "molecule_linear.h"
#include "memory.h"
#include "odesys.h"

void
molecule_dispatch_register(molecule_t * molecule,
			   int (*tdse_rhs) (const molecule_t *mol, 
					    const laser_collection_t *lasers, 
					    const double t, const double *coef, double *deriv),
			   int (*get_tdse_job) (molecule_t *self, 
						molecule_tdse_worker_t *worker,
						double *coef, double *weight),
			   void (*set_tdse_job_done) (molecule_t *self, 
						      molecule_tdse_worker_t *worker),
			   int (*check_populations)(const molecule_t *self, 
						    const double * coef),
			   molecule_tdse_worker_t * (*tdse_worker_ctor) (const molecule_t *self),
			   void (*tdse_worker_dtor)(const molecule_t *self, 
						    molecule_tdse_worker_t * worker),
			   int (*get_ncoef) (const molecule_t *self),
			   void (*dtor)(struct _molecule *self)
			   )

{
  /* Register dispatch functions. */
  molecule->tdse_rhs = tdse_rhs;
  molecule->get_tdse_job = get_tdse_job;
  molecule->set_tdse_job_done = set_tdse_job_done;
  molecule->check_populations = check_populations;
  molecule->tdse_worker_ctor = tdse_worker_ctor;
  molecule->tdse_worker_dtor = tdse_worker_dtor;
  molecule->get_ncoef = get_ncoef;
  molecule->dtor = dtor;
}

molecule_t *
molecule_ctor(int (*tdse_rhs) (const molecule_t *mol, const laser_collection_t *lasers, 
			       const double t, const double *coef, double *deriv),
	      int (*get_tdse_job) (molecule_t *self, 
				   molecule_tdse_worker_t *worker,
				   double *coef, double *weight),
	      void (*set_tdse_job_done) (molecule_t *self, 
					 molecule_tdse_worker_t *worker),
	      int (*check_populations)(const molecule_t *self, 
				       const double * coef),
	      molecule_tdse_worker_t * (*tdse_worker_ctor) (const molecule_t *self),
	      void (*tdse_worker_dtor)(const molecule_t *self, 
				       molecule_tdse_worker_t * worker),
	      int (*get_ncoef) (const molecule_t *self),
	      void (*dtor)(molecule_t *self)
	      )
{
  molecule_t *molecule;

  if (MEMORY_ALLOC(molecule) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  /* Register dispatch functions. */
  molecule_dispatch_register(molecule, tdse_rhs, get_tdse_job, 
			     set_tdse_job_done, check_populations,
			     tdse_worker_ctor, tdse_worker_dtor,
			     get_ncoef, dtor);
  return molecule;
}

void
molecule_dtor (molecule_t * molecule)
{
  molecule->dtor(molecule);
}


molecule_t *
molecule_cfg_parse_ctor (const config_t * cfg)
{
  config_setting_t *setting;
  const char *type_string;

  setting = config_lookup(cfg, "molecule");

  if (setting == NULL)
    {
      fprintf(stderr, "Failed to find molecule section in config.\n");
      return NULL;
    }

  if (!config_setting_lookup_string(setting, "type", &type_string))
    {
      fprintf(stderr, "Molecule type not specified in config.\n");
      return NULL;
    }

  if (!strcasecmp(type_string, "linear"))
    {
      molecule_t *mol = linear_molecule_cfg_parse_ctor(setting);
      if (mol == NULL)
	fprintf(stderr, "Error reading linear molecule config.\n");
      return mol;
    }
  else
    {
      fprintf(stderr, "Molecule type %s not recognised in config.\n",
	      type_string);
      return NULL;
    }
}

molecule_tdse_worker_t * 
molecule_tdse_worker_ctor()
{
  molecule_tdse_worker_t * worker;
  
  if (MEMORY_ALLOC(worker) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  worker->state = TW_WAITING;

  return worker;
}

void molecule_tdse_worker_dtor(molecule_tdse_worker_t * worker)
{
  MEMORY_FREE(worker);
}
