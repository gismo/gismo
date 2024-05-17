/** @file gsWeightMapper_.cpp

    @brief instantiation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsCore/gsTemplateTools.h>

#include <gsMSplines/gsWeightMapper.h>
#include <gsMSplines/gsWeightMapper.hpp>

#include <gsMSplines/gsWeightMapperUtils.h>
#include <gsMSplines/gsWeightMapperUtils.hpp>

namespace gismo
{

    CLASS_TEMPLATE_INST gsWeightMapper<real_t> ;

    //Utils

    TEMPLATE_INST
    index_t reorderMapperTarget (gsWeightMapper<real_t> &mapper, const std::vector<index_t>& permutation, gsPermutationMatrix* permMatrix);
    TEMPLATE_INST
    gsWeightMapper<real_t>* combineMappers (const std::vector<gsWeightMapper<real_t>*> &mappers, std::vector<index_t> &shifts, bool needShifting);
    TEMPLATE_INST
    void combineMappers (const std::vector<gsWeightMapper<real_t>*> &mappers, bool needShifting, std::vector<index_t> &shifts, gsWeightMapper<real_t>& result);


} // end namespace gismo
