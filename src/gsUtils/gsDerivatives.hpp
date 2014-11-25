
#pragma once

#include <gsUtils/gsCombinatorics.h>

namespace gismo {

template<int Ddim, int Tdim, int Order>
class computedDerivativeSizeT;

// INTERFACE FUNCTIONS FOR COMPUTING NUMBER OF PARTIAL DERIVATIVES
// BOTH COMPILE AND EXECUTION TIME
// computedDerivativeSize (int Ddim, int Tdim, int Order)   -> execution time
// computedDerivativeSize<int Ddim, int Tdim, int Order> () -> compile time

int computedDerivativeSize(int Ddim, int Tdim, int Order);

template <int Ddim, int Tdim, int Order>
int computedDerivativeSize ()
{
    return computedDerivativeSizeT<Ddim,Tdim,Order>::size;
}

int computedDerivativeSize(int Ddim, int Tdim, int Order)
{
    return Tdim*binomial(Ddim+Order-1,Order);
}

// INTERFACE FUNCTIONS TO CONVERT BETWEEN SINGLE
// AND MULTI INDEX NOTATION
// derivatives of multivariate functions are stored in a vector
// following some rules:
// pure derivatives first
// mixed derivatives in lexicographic order afterwards

// this functions convert a sequence of non decreasing indexes
// to the position of the corresponding derivative in the vector
// (0,0)-> second derivative with respect to the first coordinate
// (1,2)->second mixed derivative with respect of the second and third coordinate
// (0,0,1)->third order derivative 2 times with respect to the first coordinate and
//          one time with respect to the second coordinate
// IF THE SEQUENCE IS NOT NON DECREASING GARBAGE OR SEGMENTATION ERRORS WILL HAPPEN!

// second order derivatives
// modified lexicographic order
template <int Ddim, int Tdim>
inline int derivativeToIndex (int a,int b)
{
    return Tdim*( a==b ? a : Ddim+(2*Ddim-a-1)*a/2+b-1-a );
}


// third order derivatives
// lexicographic order
template <int Ddim, int Tdim>
inline int derivativeToIndex (int a,int b, int c)
{
    return Tdim*(a*(2+a*a-3*a*(1+Ddim)+3*Ddim*(2+Ddim))/6 +(2*Ddim-2*a-(b-a)+1)*(b-a)/2+(c-b));
}


// TEMPLATE IMPLEMENTATION
template<int Ddim, int Tdim, int Order>
class computedDerivativeSizeT
{
public:
   enum { size= Tdim*binomialT<Ddim+Order-1,Order>::value};
};

}; // gismoo namespace
