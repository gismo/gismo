/** @file gsDofMapper_.cpp.in

    @brief instantiation of gsDofMapper

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Bressan
**/

#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsDofMapper.h>
#include <gsCore/gsDofMapper.hpp>


namespace gismo {

    TEMPLATE_INST void gsDofMapper::init(
         const gsMultiBasis<real_t> & bases, index_t nComp);

    TEMPLATE_INST void gsDofMapper::init(
            std::vector<const gsMultiBasis<real_t> *> const & bases);

      TEMPLATE_INST void gsDofMapper::init(
        const gsMultiBasis<real_t>         &basis,
        const gsBoundaryConditions<real_t> &bc, int unk
        );

    TEMPLATE_INST void gsDofMapper::initSingle(
        const gsBasis<real_t> & bases, index_t nComp);
}


