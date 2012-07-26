#include <stdlib.h>
#include <libconfig.h>
#include <gsl/gsl_complex.h>

#include "config_gsl_complex.h"

int
config_setting_lookup_gsl_complex(const config_setting_t * setting, const char * name, 
				  gsl_complex * value)
{
  config_setting_t * c = config_setting_get_member(setting, name);
  
  if (c == NULL || (config_setting_length(c) != 2) )
    return CONFIG_FALSE;
  else
    {
      double re = config_setting_get_float_elem(c, 0);
      double im = config_setting_get_float_elem(c, 1);
      GSL_SET_COMPLEX(value, re, im);
      return CONFIG_TRUE;
    }

}
