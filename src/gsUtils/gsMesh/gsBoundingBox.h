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

}; // namespace gismo

#endif

