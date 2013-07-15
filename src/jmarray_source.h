#ifndef __JMARRAY_@@TYPE@@_H__
#define __JMARRAY_@@TYPE@@_H__

typedef struct _JMarray_@@TYPE@@ JMarray_@@TYPE@@_t;

JMarray_@@TYPE@@_t * JMarray_@@TYPE@@_ctor (const int Jmax);
void JMarray_@@TYPE@@_dtor (JMarray_@@TYPE@@_t * a);
@@TYPE@@ JMarray_@@TYPE@@_get (JMarray_@@TYPE@@_t * a, const int two_J, const int two_M);
void JMarray_@@TYPE@@_set (JMarray_@@TYPE@@_t * a, const int two_J, const int two_M, const @@TYPE@@ val);

#endif /* __JMARRAY_@@TYPE@@_H__ */
