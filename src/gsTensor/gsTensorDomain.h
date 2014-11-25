//
// Purpose: iterates over elements, samples domain, 
// produces a quadrature rule on the domain
//
#pragma once

#include <ostream>

#include <gsCore/gsDomain.h>

namespace gismo
{

  /** 
      Class 
  */

  
template<class T, unsigned d>
class gsTensorDomain : public gsDomain<T>
{

public:

  /// Default empty constructor
  gsTensorDomain() { };

  ~gsTensorDomain() { }; //destructor


public:



// Data members
private:


}; // TensorDomain gsClass


//////////////////////////////////////////////////
//////////////////////////////////////////////////




}; // namespace gismo
