#ifndef __JMARRAY_double_H__
#define __JMARRAY_double_H__

typedef struct _JMarray_double JMarray_double_t;

JMarray_double_t *JMarray_double_ctor (const int Jmax);
void JMarray_double_dtor (JMarray_double_t * a);
double JMarray_double_get (JMarray_double_t * a, const int two_J,
			   const int two_M);
void JMarray_double_set (JMarray_double_t * a, const int two_J,
			 const int two_M, const double val);

#endif /* __JMARRAY_double_H__ */
