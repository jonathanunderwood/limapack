#ifndef __LASER_H__ 
#define __LASER_H__ 

#include "laser_polzn.h"

/* Laser definitions employ a crude form of single inheritance implemented using
   structs. The following is the parent class struct for all laser types. This
   struct essentially defines all call-back functions that laser types should
   define and initialize. */

typedef struct _laser {
  laser_polzn_vector_t * (*get_polzn_vector) (const struct _laser *self, const double t);
  double (*get_envelope) (const struct _laser *self, const double t);
  double (*get_frequency) (const struct _laser *self, const double t);
  //  laser_type_t type;
} laser_t;

/* Constructor */
laser_t * laser_ctor (laser_polzn_vector_t * (*get_polzn_vector) (const laser_t *self, const double t),
		      double (*get_envelope) (const laser_t *self, const double t),
		      double (*get_frequency) (const laser_t *self, const double t)
		      );

/* Function to register dispatch functions. */
void laser_dispatch_register (laser_t * laser,
			      laser_polzn_vector_t * (*get_polzn_vector) (const laser_t *self, const double t),
			      double (*get_envelope) (const laser_t *self, const double t),
			      double (*get_frequency) (const laser_t *self, const double t)
			      );

/* Destructor */
void laser_dtor (laser_t *self);

/* Macro to convert intensity in W/cm^2 to electric field in V/m in a
   vacuum. See Appendix C of Boyd's Nonlinear Optics. Note extra factor of 2
   compared to Boyd Eq. C3 due to missing half in definition Eq. C1. In other
   words, this is correct if we define E(t) = E0 cos (wt). */
/* #define LASER_ITOE(I) (sqrt (377.0 * 1.0e4 * I * 2.0)) */
#define LASER_I_TO_E(I) (sqrt (7.54e6 * I))


#endif /* __LASER_H__ */
