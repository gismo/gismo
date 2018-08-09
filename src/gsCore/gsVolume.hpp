/** @file gsVolume.hpp

    @brief Provides implementation of Volume common operations.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsUtils/gsPointGrid.h>
#include <gsUtils/gsMesh/gsMesh.h>

namespace gismo
{

template<class T> 
gsGeometryEvaluator<T> *
gsVolume<T>::evaluator(unsigned flags) const
{
    switch ( this->coDim() )
    {
    case 0:
        return new gsGenericGeometryEvaluator<T,3,0 >(*this, flags);
    case 1:
        return new gsGenericGeometryEvaluator<T,3,1 >(*this, flags);
    case -2:
        return new gsGenericGeometryEvaluator<T,3,-2>(*this, flags);
    default:
        GISMO_ERROR("Codimension problem.( codim="<<this->coDim() );
    }
}

}; // namespace gismo

