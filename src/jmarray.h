#ifndef __JMARRAY_H__
#define __JMARRAY_H__

/* These two utility functions are defined here to allow for use
   outside of the JMarray_foo_t types. The first maps a (J,M) pair to
   a 1D array index. The second returns the size needed for such a 1D
   array for a given maximum value of J. */
extern inline int
JMarray_idx (const int J, const int M)
{
  return J * J + J + M;
}

extern inline int
JMarray_dim (const int Jmax)
{
  return Jmax * Jmax + 2 * Jmax + 1;
}

#endif /* __JMARRAY_H__ */
