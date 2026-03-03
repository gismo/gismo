/** @file gsSpringPatch.h

    @brief Provides spring patch construction from boundary data.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include<gsModeling/gsPatchGenerator.h>

namespace gismo
{

/**
   \brief Computes a parametrization based on the spring patch
   technique, given a set of boundary geometries.
*/
template <typename T>
class gsSpringPatch : public gsPatchGenerator<T>
{
public:
    typedef gsPatchGenerator<T> Base;

public:

    /**
       \brief Constructs a spring patch object by a collection of
       tensor-product patches defining the boundaries of a patch.  
       
       The number of boundaries is expected to be 2*d, where d is the
       space dimension. They are expected to form a closed "shell"

       \param boundary a set of boundary curves or patches
    */
    gsSpringPatch(const gsMultiPatch<T> & boundary): Base(boundary)
    { }

public:

    /// \brief Main routine that performs the computation
    const gsGeometry<T> & compute();

private:

    template<unsigned d> void compute_impl();

protected:

    using Base::m_boundary;
   
    using Base::m_result;

}; // gsSpringPatch


}// namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsSpringPatch.hpp)
#endif
