/** @file gsBoundingBox.h

    @brief Provides gsBoundingBox class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D.-M. Nguyen
*/

#ifndef _BOUNDINGBOX_H_
#define _BOUNDINGBOX_H_

#include <gsCore/gsLinearAlgebra.h>

namespace gismo {


template <class T>
class gsBoundingBox {
public:
    gsBoundingBox() { };
    gsBoundingBox(const gsVector3d<T>& pMin, const gsVector3d<T>& pMax) {
        this->pMin = pMin;
        this->pMax = pMax;
    }
/// data members
    gsVector3d<T> pMin, pMax;
    
    /// get the maximum length of the edges of the bounding box
    T getMaxSize() const {return (pMax-pMin).maxCoeff();};
};

} // namespace gismo

#endif

