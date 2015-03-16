/** @file gsTemplateTools.h

    @brief Utilities related to template programming.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo
{

/*
  Compile time type switching.
  Example: The type

  gismo::choose<k==1,double, gsVector<double> >::type

  evaluates (at compile time) to "double" if k is equal to 1, 
  and to "gsVector<double>" otherwise
*/

template <bool flag, class IsTrue, class IsFalse>
struct choose;

template <class IsTrue, class IsFalse>
struct choose<true, IsTrue, IsFalse> { typedef IsTrue type; };

template <class IsTrue, class IsFalse>
struct choose<false, IsTrue, IsFalse> {
   typedef IsFalse type;
};

} // end namespace gismo


/*
  Macros for instantiating templated structs, classes or functions
*/

#if defined(_WIN32)
#define STRUCT_TEMPLATE_INST template struct __declspec(dllexport)
#define CLASS_TEMPLATE_INST  template class  __declspec(dllexport)
#define TEMPLATE_INST        template        __declspec(dllexport)
#else
#define STRUCT_TEMPLATE_INST template struct
#define CLASS_TEMPLATE_INST  template class 
#define TEMPLATE_INST        template       
#endif
