#ifndef __JKMARRAY_@@TYPE@@_H__
#define __JKMARRAY_@@TYPE@@_H__

typedef struct _JKMarray_@@TYPE@@ JKMarray_@@TYPE@@_t;

JKMarray_@@TYPE@@_t * JKMarray_@@TYPE@@_ctor (const int Jmax);
void JKMarray_@@TYPE@@_dtor (JKMarray_@@TYPE@@_t * a);
@@TYPE@@ JKMarray_@@TYPE@@_get (JKMarray_@@TYPE@@_t * a, const int J, const int K, const int M);
void JKMarray_@@TYPE@@_set (JKMarray_@@TYPE@@_t * a, const int J, const int K, const int M,
		            const @@TYPE@@ val);

#endif /* __JKMARRAY_@@TYPE@@_H__ */
