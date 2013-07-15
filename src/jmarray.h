#ifndef __JMARRAY_H__
#define __JMARRAY_H__

/* These two utility functions are defined here to allow for use
   outside of the JMarray_foo_t types. The first maps a (J,M) pair to
   a 1D array index. The second returns the size needed for such a 1D
   array for a given maximum value of J. */

extern inline int
JMarray_idx (const int two_J, const int two_M)
{
  return ((two_J * two_J) / 2 + two_J + two_M) / 2;
}


extern inline int
JMarray_dim (const int two_Jmax)
{
  return JMarray_idx (two_Jmax, two_Jmax) + 1;
}

#endif /* __JMARRAY_H__ */
