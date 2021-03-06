#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <libconfig.h>
#include <hdf5.h>

#ifdef BUILD_WITH_PTHREADS
#include <pthread.h>
#endif

#ifdef BUILD_WITH_MPI
#include <mpi.h>
#endif

#include "molecule.h"
#include "memory.h"

/* Include all molecular model specific header files here only so we
   can properly parse the config file. */
#include "molecule_linear.h"
#include "molecule_symrot.h"
#include "molecule_asymrot.h"

void
molecule_dispatch_register(molecule_t * molecule,
			   int (*tdse_rhs) (const molecule_t *mol, 
					    const laser_collection_t *lasers, 
					    const double t, const double *coef, double *deriv),
			   int (*get_tdse_job) (molecule_t *self, 
						molecule_tdse_worker_t *worker),
			   void (*get_tdse_job_coef) (const molecule_t *self, 
						      const molecule_tdse_worker_t *worker,
						      double *coef),
			   double (*get_tdse_job_weight) (const molecule_t *self, 
							  const molecule_tdse_worker_t *worker),
			   void (*set_tdse_job_done) (molecule_t *self, 
						      molecule_tdse_worker_t *worker),
			   int (*check_populations)(const molecule_t *self, 
						    const double * coef),
			   molecule_tdse_worker_t * (*tdse_worker_ctor) (const molecule_t *self),
			   void (*tdse_worker_dtor)(const molecule_t *self, 
						    molecule_tdse_worker_t * worker),
			   int (*get_ncoef) (const molecule_t *self),
			   molecule_expval_t * (*expval_ctor) (const struct _molecule *self),
			   void (*expval_dtor) (const struct _molecule *self, 
						molecule_expval_t *expval),
			   int (*expval_calc)  (const molecule_t *self, const double *coef,
						 const double t, molecule_expval_t *expval),
			   void (*expval_add_weighted) (const molecule_t *self, 
							molecule_expval_t *a, 
							const molecule_expval_t *b,
							const double weight),
			   int (*expval_fwrite) (const molecule_t *molecule,
						 const molecule_expval_t *expvalue,
						 const hid_t *location),
#ifdef BUILD_WITH_MPI
			   int (*expval_mpi_send) (const molecule_t *self, 
						   const molecule_expval_t *molecule, 
						   int dest, int tag, MPI_Comm comm),
			   int (*expval_mpi_recv) (const molecule_t *self, 
						   molecule_expval_t *molecule, 
						   int dest, int tag, MPI_Comm comm),
			   int (*tdse_worker_mpi_send) (const molecule_t *self, 
							const molecule_tdse_worker_t *worker,
							int dest, int tag, MPI_Comm comm),
			   int (*tdse_worker_mpi_recv) (const molecule_t *self, 
							molecule_tdse_worker_t *worker,
							int dest, int tag, MPI_Comm comm),
#endif
			   void (*dtor)(molecule_t *self)
			   )
{
  molecule->tdse_rhs = tdse_rhs;
  molecule->get_tdse_job = get_tdse_job;
  molecule->set_tdse_job_done = set_tdse_job_done;
  molecule->get_tdse_job_coef = get_tdse_job_coef;
  molecule->get_tdse_job_weight = get_tdse_job_weight;
  molecule->check_populations = check_populations;
  molecule->tdse_worker_ctor = tdse_worker_ctor;
  molecule->tdse_worker_dtor = tdse_worker_dtor;
  molecule->get_ncoef = get_ncoef;
  molecule->expval_ctor = expval_ctor;
  molecule->expval_dtor = expval_dtor;
  molecule->expval_calc = expval_calc;
  molecule->expval_add_weighted = expval_add_weighted;
  molecule->expval_fwrite = expval_fwrite;
#ifdef BUILD_WITH_MPI
  molecule->expval_mpi_send = expval_mpi_send;
  molecule->expval_mpi_recv = expval_mpi_recv;
  molecule->tdse_worker_mpi_send = tdse_worker_mpi_send;
  molecule->tdse_worker_mpi_recv = tdse_worker_mpi_recv;
#endif
  molecule->dtor = dtor;

#ifdef BUILD_WITH_PTHREADS
  pthread_mutex_init(&(molecule->lock), NULL);
#endif

  return;
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
  if (!strcasecmp(type_string, "symrot"))
    {
      molecule_t *mol = symrot_molecule_cfg_parse_ctor(setting);
      if (mol == NULL)
	fprintf(stderr, "Error reading symrot molecule config.\n");
      return mol;
    }
  if (!strcasecmp(type_string, "asymrot"))
    {
      molecule_t *mol = asymrot_molecule_cfg_parse_ctor(setting);
      if (mol == NULL)
	fprintf(stderr, "Error reading asymrot molecule config.\n");
      return mol;
    }
  else
    {
      fprintf(stderr, "Molecule type %s not recognised in config.\n",
	      type_string);
      return NULL;
    }
}

