#inclde "laser_polzn.h"

/* Laser definitions employ a crude form of single inheritance implemented using
   structs. The following is the parent class struct for all laser types. This
   struct essentially defines all call-back functions that laser types should
   define and initialize. */

typedef struct _laser {
  laser_polzn_vector * (*get_polzn_vector) (struct _laser *self, const double t);
  laser_polzn_tensor * (*get_polzn_tensor) (struct _laser *self, const double t);
  gsl_comlplex * (*get_envelope) (struct _laser *self, const double t);
} laser_t;

/* Constructor */
void laser_ctor (laser_t *self,
		 laser_polzn_vector * (*get_polzn_vector) (laser_t *self, const double t)
		 laser_polzn_tensor * (*get_polzn_tensor) (laser_t *self, const double t)
		 gsl_comlplex * (*get_envelope) (laser_t *self, const double t)
		 );

/* Destructor */
void laser_xtor (laser_t *self);
