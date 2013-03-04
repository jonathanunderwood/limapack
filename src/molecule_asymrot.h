#ifndef __MOLECULE_ASYMROT_H__
#define __MOLECULE_ASYMROT_H__

#include "molecule.h"

typedef struct _asymrot_molecule asymrot_molecule_t;

asymrot_molecule_t * asymrot_molecule_ctor();
void asymrot_molecule_dtor(molecule_t * mol);
molecule_t * asymrot_molecule_cfg_parse_ctor(const config_setting_t *cfg);

#endif /* __MOLECULE_ASYMROT_H__ */
