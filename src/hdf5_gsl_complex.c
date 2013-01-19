#include <hdf5.h>
#include <gsl/gsl_complex.h>

#include "hdf5_gsl_complex.h"

hid_t
hdf5_gsl_complex_type_create()
{
  hid_t complex_id = H5Tcreate (H5T_COMPOUND, sizeof (gsl_complex));
  int ret;

  if (complex_id >= 0)
    {
      ret = H5Tinsert (complex_id, "real", HOFFSET(gsl_complex, dat[0]), H5T_NATIVE_DOUBLE);
      ret += H5Tinsert (complex_id, "imaginary", HOFFSET(gsl_complex, dat[1]),
			H5T_NATIVE_DOUBLE);
      if (ret >= 0)
	return complex_id;
      else 
	  H5Tclose (complex_id);
    }

  /* If we get here, something went wrong. */
  return -1;
}
