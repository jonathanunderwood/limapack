#include "config_gsl_complex.h"

int
config_setting_lookup_gsl_complex(const config_setting_t * setting, const char * name, 
				  gsl_complex * value)
{
  c = config_setting_get_member(element, name);
  
  if (c == NULL || (config_setting_length(c) != 2) )
    return CONFIG_FALSE;
  else
    {
      double re = config_setting_get_float_elem(c, 0);
      double im = config_setting_get_float_elem(c, 1);
      gsl_set_complex(value, re, im);
      return CONFIG_TRUE;
    }

}
