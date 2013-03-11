#ifndef __MOLECULE_SYMROT_H__
#define __MOLECULE_SYMROT_H__

#include "molecule.h"

typedef struct _symrot_molecule symrot_molecule_t;

symrot_molecule_t * symrot_molecule_ctor();
void symrot_molecule_dtor(molecule_t * mol);
molecule_t * symrot_molecule_cfg_parse_ctor(const config_setting_t *cfg);

#endif /* __MOLECULE_SYMROT_H__ */
