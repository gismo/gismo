// Tensor Nurbs geometry with with weights constaint to product of 1D Nurbs-component weights

#pragma once

#include <gsNurbs/gsQuasiTensorNurbsBasis.h>

namespace gismo
{

template<class T, unsigned d> class gsQuasiTensorNurbs;


template<class T, unsigned d>
class gsQuasiTensorNurbs : 
public gsGenericGeometry< gsQuasiTensorNurbsBasis<d,T> > 
{

// TO DO later

};





}; // namespace gismo

