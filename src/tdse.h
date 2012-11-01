#ifndef __TDSE_H__
#define __TDSE_H__


tdse_worker_t * tdse_worker_ctor();
void tdse_worker_dtor(tdse_worker_t * worker);

/* Infrastructure for keeping track of propagation of initial
   rotational states. The following structure needs to be
   "sub-classed" by the specific molecule type to add extra parameters
   to it which describe the jobs i.e. quantum numbers. */
typedef struct _tdse_job_queue
{
  void (*get_job)  (struct _tdse_job_queue, tdse_worker_t worker);
  void (*set_job_done) (struct _tdse_job_queue, tdse_worker_t worker);
} tdse_job_queue_t;

tdse_job_queue_t * 
tdse_job_queue_ctor(void (*get_job) (struct _tdse_job_queue, 
				     tdse_worker_t worker),
		    void (*set_job_done) (struct _tdse_job_queue, 
					  tdse_worker_t worker));
void tdse_job_queue_dtor (tdse_job_queue_t * queue);


/* This simple propagator function simply propagates each initial
   state serially. */
int tdse_propagate_simple (odesys_t *ode, tdse_job_queue_t *queue, 
			   tdse_worker_t *worker);

#endif
