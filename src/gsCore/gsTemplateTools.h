
#pragma once

#include <gsCore/gsExport.h>

namespace gismo
{

template<class T>
class gsKnotVector;

template <bool flag, class IsTrue, class IsFalse>
struct choose;

template <class IsTrue, class IsFalse>
struct choose<true, IsTrue, IsFalse> {
   typedef IsTrue type;
};

template <class IsTrue, class IsFalse>
struct choose<false, IsTrue, IsFalse> {
   typedef IsFalse type;
};

} // end namespace gismo


/*
  Macros for instantization
*/

#define STRUCT_TEMPLATE_INST template struct GISMO_EXPORT
#define CLASS_TEMPLATE_INST  template class  GISMO_EXPORT
#define TEMPLATE_INST        template        GISMO_EXPORT
