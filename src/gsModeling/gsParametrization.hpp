/** @file gsParametrization.hpp

    @brief Provides implementation gsParametrization class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl

*/

#include <gsIO/gsOptionList.h>

namespace gismo
{


template<class T>
gsParametrization<T>::gsParametrization(gismo::gsOptionList optionList) : m_mesh(MeshInfo(optionList.getString("filenameIn")))
{
    std::string filenameIn = optionList.getString("filenameIn");
    std::string filenameOut = optionList.getString("filenameOut");
    std::string boundaryMethod = optionList.getString("boundaryMethod");
    std::string paraMethod = optionList.getString("parametrizationMethod");
    std::vector<int> cornersInput = optionList.getMultiInt("corners");
    double range = optionList.getReal("range");
    index_t number = optionList.getInt("number");
}

} // namespace gismo
