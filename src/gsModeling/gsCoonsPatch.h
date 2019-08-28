/** @file gsCoonsPatch.h

    @brief Provides Coons's patch construction from boundary data.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Haberleitner, A. Mantzaflaris
*/

#pragma once

#include<gsModeling/gsPatchGenerator.h>

namespace gismo
{


/**
   \brief Computes a Coons' patch parametrization given a set of
   boundary geometries.  Parametrization is not guaranteed to be
   non-singular. Works for surface, volumes, or any dimension.

   \tparam T Coefficient type

   \ingroup Modeling
*/
template <typename T>
class gsCoonsPatch : public gsPatchGenerator<T>
{
public:
    typedef gsPatchGenerator<T> Base;
public:

    /**
       \brief Constructs a Coon's patch object by a collection of
       tensor-product patches defining the boundaries of a patch.  
       
       The number of boundaries is expected to be 2*d, where d is the
       space dimension. They are expected to form a closed "shell"

       \param boundary a set of boundary curves or patches
    */
    gsCoonsPatch(const gsMultiPatch<T> & boundary) : Base(boundary)
    { }

public:

    // Look at gsPatchGenerator
    const gsGeometry<T> & compute();

private:

    template<short_t d> void compute_impl();

protected:

    using Base::m_boundary;
   
    using Base::m_result;
    
}; // gsCoonsPatch


}// namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCoonsPatch.hpp)
#endif
