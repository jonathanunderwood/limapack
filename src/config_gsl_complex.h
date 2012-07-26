#ifndef __CONFIG_GSL_COMPLEX_H__
#define __CONFIG_GSL_COMPLEX_H__

#include <gsl/gsl_complex.h>
#include <libconfig.h>

int config_setting_lookup_gsl_complex(const config_setting_t * setting, 
				      const char * name, 
				      gsl_complex * value);

#endif /* __CONFIG_GSL_COMPLEX_H__ */
