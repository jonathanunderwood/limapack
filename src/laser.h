#ifndef __LASER_H__ 
#define __LASER_H__ 

#include "laser_polzn.h"

/* Macro to convert intensity in W/cm^2 to electric field in V/m in a
   vacuum. See Appendix C of Boyd's Nonlinear Optics. Note extra factor of 2
   compared to Boyd Eq. C3 due to missing half in definition Eq. C1. In other
   words, this is correct if we define E(t) = E0 cos (wt). */
/* #define LASER_ITOE(I) (sqrt (377.0 * 1.0e4 * I * 2.0)) */
#define LASER_I_TO_E(I) (sqrt (7.54e6 * I))

/* Laser type enumerator - this needs augmenting as new laser types
   are added. */
typedef enum _laser_type {
  LASER_TYPE1 = 1
} laser_type_t;


/* Laser definitions employ a crude form of single inheritance implemented using
   structs. The following is the parent class struct for all laser types. This
   struct essentially defines all call-back functions that laser types should
   define. 

   Note on get_envelope: if this function returns a value < 0, we use
   that to denote "skip this laser" when doing calculations
   eg. evaluating the RHS of the TDSE.
*/

typedef struct _laser {
  laser_polzn_vector_t * (*get_polzn_vector) (const struct _laser *self, const double t);
  double (*get_envelope) (const struct _laser *self, const double t);
  double (*get_frequency) (const struct _laser *self, const double t);
  void (*dtor) (struct _laser *self);
} laser_t;

/* Constructor */
/* laser_t * laser_ctor (laser_polzn_vector_t * (*get_polzn_vector) (const laser_t *self,  */
/* 								  const double t), */
/* 		      double (*get_envelope) (const laser_t *self, const double t), */
/* 		      double (*get_frequency) (const laser_t *self, const double t), */
/* 		      void (*dtor) (laser_t *self) */
/* 		      ); */

/* Function to register dispatch functions. */
void laser_dispatch_register (laser_t * laser,
			      laser_polzn_vector_t * (*get_polzn_vector) (const laser_t *self, 
									  const double t),
			      double (*get_envelope) (const laser_t *self, const double t),
			      double (*get_frequency) (const laser_t *self, const double t),
			      void (*dtor) (laser_t *self)
			      );

/* Destructor */
/* void laser_dtor (laser_t *self); */

/* Structure to keep track of multiple lasers of different types. */
typedef struct _laser_collection
{
  laser_t **laser;
  int nlasers;
} laser_collection_t;

laser_collection_t * laser_collection_ctor (const int nlasers);
void laser_collection_dtor (laser_collection_t * self);

laser_collection_t *laser_collection_cfg_parse_ctor (const config_t * cfg);
int laser_collection_all_negligible (laser_collection_t * self, const double t);
#endif /* __LASER_H__ */
