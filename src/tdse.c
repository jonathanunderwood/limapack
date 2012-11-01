#include <stdlib.h>
#include <stdio.h>

#include "odesys.h"
#include "molecule.h"
#include "tdse.h"

tdse_worker_t * 
tdse_worker_ctor()
{
  tdse_worker_t * worker;
  
  if (MEMORY_ALLOC(worker) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }

  return worker;
      
}

void tdse_worker_dtor(tdse_worker_t * worker)
{
  MEMORY_FREE(worker);
}


tdse_worker_t * 
tdse_worker_ctor()
{
  tdse_worker_t * worker;

  if (MEMORY_ALLOC(worker) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }
  return worker;
}

void 
tdse_worker_dtor(tdse_worker_t * worker)
{
  MEMORY_FREE(worker);
}

tdse_job_queue_t * 
tdse_job_queue_ctor(void (*get_job) (struct _tdse_job_queue, 
				     tdse_worker_t worker),
		    void (*set_job_done) (struct _tdse_job_queue, 
					  tdse_worker_t worker),
		    )
{
  tdse_job_queue_t queue;

  if (MEMORY_ALLOC(queue) < 0)
    {
      MEMORY_OOMERR;
      return NULL;
    }
  /* Register dispatch methods. */
  queue->get_job = get_job;
  queue->set_job_done = set_job_done;

  return queue;
}

void 
tdse_job_queue_dtor (tdse_job_queue_t * queue)
{
  MEMORY_FREE(queue);
}


int 
tdse_propagate_simple (odesys_t *ode, 
		       tdse_job_queue_t *queue, 
		       tdse_worker_t *worker)
/* Simple single threaded generic propagator which uses a single
   worker. */
{
  molecule_t *mol = ode->params->molecule;
  laser_collection_t *las = ode->params->lasers;
  double *coef;
  double weight;

  if (MEMORY_ALLOC_N(coef, 2 * mol->ncoef) < 0)
    {
      MEMORY_OOMERR;
      return -1;
    }

  /* Step through initial states of molecule (eg. the individual
     states in a Boltzmann distribution). For each initial state,
     propagate the TDSE, calculating the expectation values at each
     time-step, and add them to the ensemble averaged expectation
     value, weighted appropriately. */
  while (queue.get_job(queue, worker, &coef, &weight) == 0) 
   { 
     unsigned char need_ode_reset = 0;
     ode_reset(ode);
     t1 = ode->tstart;
     
     for (i = 0; i < ode->npoints; i++) /* Step through time points. */
       {
	 /* Check that populations in highest levels aren't growing
	    unacceptably. */
	 if (mol->check_populations(mol, coef) != 0)
	   {
	     fprintf(stderr, "Populations building up unacceptably in TDSE propagation. Exiting.");
	     MEMORY_FREE(coef);
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
	 if (t1 < ode->tend)
	   {
	     if (las->check_if_negligible(las))
	       need_ode_reset = 1;
	     else
	       {
		 if (need_ode_reset)
		   {
		     ode_reset(ode);
		     need_ode_reset = 0;
		   }
		 
		 odesys_step(ode, t1, t1 + ode->tstep, &coef);
	       }
	     
	     t1 += ode->tstep;
	   }
       }
     queue.set_tdse_job_done(queue, worker);
   }
}
