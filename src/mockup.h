// an instance of this should be in the molecule structure is molecule
// type specific, obviously.  we'll need to modify jmarray and
// jkmarray types to allow int arrays, or even small int.
typedef struct _tdse_queue
{
  jmarray tdse_status; //0 = todo, 1 = started, 2 = done 
} tdse_queue;

// returns suitably initialised coef array and a statistical weight.
// this is a type specific function which needs registering with the molecule_t
// return 1 if no jobs in TODO state
// return 2 if all jobs done
int tdse_get_next(molecule_t *mol, double **coef, double *wt);

// this will be called after each time step to check that population
// isn't growing in the highest energy states - implementation is
// molecule specific.
int tdse_check_populations(molecule_t *mol, double *coef)
