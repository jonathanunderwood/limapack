/* Polarizability tensor spherical components. */
#ifndef __POLARIZABILITY_H__
#define __POLARIZABILITY_H__

#include "jmarray.h"

typedef JMarray_t polarizability_t;

polarizability_t * polarizability_ctor_from_cart (const double a_xx, 
						  const double a_yy, 
						  const double a_zz);
polarizability_t * polarizability_ctor ();
void polarizability_dtor(polarizability_t *alpha);

#endif /* __POLARIZABILITY_H__ */
