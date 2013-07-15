#ifndef __JKMARRAY_double_H__
#define __JKMARRAY_double_H__

typedef struct _JKMarray_double JKMarray_double_t;

JKMarray_double_t *JKMarray_double_ctor (const int two_Jmax);
void JKMarray_double_dtor (JKMarray_double_t * a);
double JKMarray_double_get (JKMarray_double_t * a, const int two_J,
			    const int two_K, const int two_M);
void JKMarray_double_set (JKMarray_double_t * a, const int two_J,
			  const int two_K, const int two_M, const double val);

#endif /* __JKMARRAY_double_H__ */
