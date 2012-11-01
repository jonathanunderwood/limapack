#ifndef __JMARRAY_@@TYPE@@_H__
#define __JMARRAY_@@TYPE@@_H__

typedef struct _JMarray_@@TYPE@@
{
  double (*get) (struct _JMarray_@@TYPE@@ * self, const int J, const int M);
  void (*set) (struct _JMarray_@@TYPE@@ * self, const int J, const int M,
	       const double val);
  double *data;
  int Jmax;
  int dim;
} JMarray_@@TYPE@@_t;

JMarray_@@TYPE@@_t * JMarray_@@TYPE@@_ctor (const int Jmax);
void JMarray_@@TYPE@@_dtor (JMarray_@@TYPE@@_t * a);
@@TYPE@@ JMarray_@@TYPE@@_get (JMarray_@@TYPE@@_t * a, const int J, const int M);
void JMarray_@@TYPE@@_set (JMarray_@@TYPE@@_t * a, const int J, const int M, const @@TYPE@@ val);

#endif /* __JMARRAY_@@TYPE@@_H__ */
