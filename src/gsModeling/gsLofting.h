/** @file gsLofting.h

    @brief Provides declaration of data fitting algorithms by least
    squares approximation.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
# include <gsNurbs/gsKnotVector.h>
# include <gsNurbs/gsBSpline.h>
# include <gsCore/gsMultiPatch.h>

namespace gismo
{

/**
  @brief 
   Class for performing a lofting gsGeometry.
    
   \ingroup Modeling
**/
template<class T>
class gsLofting
{
    
public:// typedefs

    typedef std::vector<gsBSpline<T>> curveContainer;

    
    

public:// constructors
    /// default constructor
    gsLofting()
    {
        m_basis = NULL;
        m_result= NULL ;
    }

    /// destructor
    ~gsLofting();

    // explicit gsLofting(curveContainer container);
    explicit gsLofting(gsMultiPatch<T> container);

    explicit gsLofting(gsMultiPatch<T> container, index_t deg_v);

    explicit gsLofting(gsMultiPatch<T> container, gsKnotVector<T> kv);

    /// constructor
    // gsLofting(std::vector<gsBSpline<T>> const & container);


public:

    /// Compues the isoparameters in v-direction.
    void isovparameters(T alpha, gsMatrix<T> & v_parameters);

    /// Computes the least squares fit for a gsBasis
    void make_compatible();

    /// compute the lofting surface.
    void compute();

    /// gives back the computed approximation
    gsGeometry<T> * result() const { return m_result; }

    // /// gives back the computed approximation for multipatch geometry
    // const gsMappedSpline<2,T> & mresult() const { return m_mresult; }

    // /// Returns the basis of the approximation
    // const gsBasis<T> & getBasis() const {return *static_cast<const gsBasis<T>*>(m_basis);}

    

protected:

    //gsOptionList

    /// the container of B-spline curves to be lofted
    // std::vector<gsBSpline<T>> m_container;
    gsMultiPatch<T> m_container;

    gsKnotVector<T> m_ku;
    gsKnotVector<T> m_kv;

    index_t m_vdeg;

    index_t m_N;

    index_t m_m;

    index_t m_DIM;

    
    

    /// Pointer keeping the basis
    // gsFunctionSet<T> * m_basis;

    gsTensorBSplineBasis <2, T> * m_basis;

    /// Pointer keeping the resulting geometry
    gsGeometry<T> * m_result;

    /// Maximum point-wise error
    T m_max_error;

    /// Minimum point-wise error
    T m_min_error;

    // All point-wise errors
    std::vector<T> m_pointErrors;

    mutable T m_last_lambda;


private:
   

}; // class gsLofting


// #ifdef GISMO_WITH_PYBIND11

//   /**
//    * @brief Initializes the Python wrapper for the class: gsKnotVector
//    */
//   void pybind11_init_gsFitting(pybind11::module &m);

// #endif // GISMO_WITH_PYBIND11


}// namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsLofting.hpp)
#endif
