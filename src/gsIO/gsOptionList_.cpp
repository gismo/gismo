/** @file gsOptionList_.cpp

    @brief instantiation of gsOptionList

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
**/

#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsOptionList.h>
#include <gsIO/gsOptionList.hpp>

namespace gismo {

TEMPLATE_INST
real_t gsOptionList::getReal<real_t>(const std::string & label) const;

TEMPLATE_INST
std::vector<real_t> gsOptionList::getMultiReal<real_t>(const std::string & gn) const;

TEMPLATE_INST
real_t gsOptionList::askReal<real_t>(const std::string & label, const real_t & value) const;

TEMPLATE_INST
void gsOptionList::setReal<real_t>(const std::string & label, const real_t & value);

TEMPLATE_INST
void gsOptionList::addReal<real_t>(const std::string & label,
                                   const std::string & desc,
                                   const real_t & value);

}


