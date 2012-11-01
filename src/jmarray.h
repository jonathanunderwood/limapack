#ifndef __JMARRAY_H__
#define __JMARRAY_H__

/* These two utility functions are defined here to allow for use
   outside of the JMarray_foo_t types. */
extern inline int
JMarray_dim (const int Jmax)
{
  return Jmax * Jmax + 2 * Jmax + 1;
}

extern inline int
JMarray_idx (const int J, const int M)
{
  return J * J + J + M;
}

/* Integer variant. */
typedef struct _JMarray_int
{
  int (*get) (struct _JMarray_int * self, const int J, const int M);
  void (*set) (struct _JMarray_int * self, const int J, const int M,
	       const int val);
  double *data;
  int Jmax;
  int dim;
} JMarray_int_t;

JMarray_int_t * JMarray_int_ctor (const int Jmax);
void JMarray_int_dtor (JMarray_int_t * a);
int JMarray_int_get (JMarray_int_t * a, const int J, const int M);
void JMarray_int_set (JMarray_int_t * a, const int J, const int M, const int val);

/* Double variant. */
typedef struct _JMarray_double
{
  double (*get) (struct _JMarray_double * self, const int J, const int M);
  void (*set) (struct _JMarray_double * self, const int J, const int M,
	       const double val);
  double *data;
  int Jmax;
  int dim;
} JMarray_double_t;

JMarray_double_t * JMarray_double_ctor (const int Jmax);
void JMarray_double_dtor (JMarray_double_t * a);
double JMarray_double_get (JMarray_double_t * a, const int J, const int M);
void JMarray_double_set (JMarray_double_t * a, const int J, const int M, const double val);

#endif /* __JMARRAY_H__ */
