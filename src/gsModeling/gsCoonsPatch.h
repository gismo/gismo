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

    template<unsigned d> void compute_impl();

protected:

    using Base::m_boundary;
   
    using Base::m_result;
    
}; // gsCoonsPatch



template <typename T> const 
gsGeometry<T> & gsCoonsPatch<T>::compute()
{
    const int dim = m_boundary.parDim();

    delete m_result;// delete any previous result
    m_result = NULL;

    switch ( dim ) // dispatch to implementation
    {
    case 1:
        compute_impl<2>();
        break;
    case 2:
        compute_impl<3>();
        break;
    case 3:
        compute_impl<4>();
        break;
    default:
        GISMO_ERROR("Dimension "<< dim << "is invalid.");
        break;
    }

    return *m_result;
}

template <typename T> template <unsigned d>
void gsCoonsPatch<T>::compute_impl()
{
    gsTensorBSplineBasis<d,T> resultBasis; // Basis for the Coon's patch
    gsMatrix<T> coefs;                     // CPs of the Coon's patch

    // Resolve boundary configuration and set boundary coefficients
    this->preparePatch(resultBasis, coefs);

    //-------- Compute interior control points    
    
    // Note: assuming that the param. domain is [0,1]^d
    gsMatrix<T> gr[d]; // component-wise Greville points

    gsVector<unsigned,d> vend,  // Multi-index of the furthest vertex
        stride, // The strides of the tensor-basis
        tmp;
    
    for (unsigned k = 0; k < d; ++k)
    {
        gr    [k] = resultBasis.component(k).anchors();
        vend  [k] = resultBasis.size(k) - 1;
        stride[k] = resultBasis.stride(k);
    }
    
    // Check whether there are any interior points to fill in
    if ( (vend.array() < 2).any() )
    {
        gsWarn<<"There where no interior control points.\n";
        m_result = resultBasis.makeGeometry( give(coefs) ).release();
        return;
    }
    
    // Multi-index of current interior CP (in iteration)
    gsGridIterator<unsigned,CUBE,d> grid(gsVector<unsigned,d>::Ones(), vend);
    // Multi-index of current cube element stencil (in iteration) ( [-1,1]^d --> [0,2]^d )
    static const gsVector<unsigned,d> twos = gsVector<unsigned,d>::Constant(2);
    gsGridIterator<unsigned,CUBE,d> cf(twos, false);

    for(; grid; ++grid) // loop over all interior coefficients
    {
        // grab current coefficient
        typename gsMatrix<T>::RowXpr curCoef = coefs.row(stride.dot(*grid));
        
        cf.reset();
        ++cf;// skip (0..0), ie. the coordinates of CP "*grid"
        
        for(; cf; ++cf) // loop over all cube elements (stencil around *grid)            
        {   
            // cf:  ***
            //      *+*
            //      ***
            
            // Initialize weight sign
            T w = -1.0;
            
            for(unsigned k=0; k!=d; k++)
            {
                // Compute weight
                if (0 != cf->at(k) ) // coordinate not equal to cf->at(k) ?
                {
                    // fetch Greville point
                    const T x = gr[k](0, grid->at(k));
                    // update weight
                    w *= - ( 2 == cf->at(k) ? x : 1.0-x );
                }
                
                // Compute index of contributing coefficient
                // 0 : CP index coordinate of current CP grid->at(k)
                // 1 : fixed to corner coordinate 0
                // 2 : fixed to corner coordinate vend[k]
                // tmp: *--*--*
                //      |     |
                //      *  +  *
                //      |     |
                //      *--*--*
                tmp[k] = ( cf->at(k)==0 ? grid->at(k)  :
                           cf->at(k)==1 ? 0 : vend[k] );
            }

            // Add contribution of current cube element "*cf" to curCoef
            curCoef += w * coefs.row( stride.dot(tmp) ) ;

        }//cf
    }//grid

    // Coons' patch is ready
    m_result = resultBasis.makeGeometry( give(coefs) ).release();
}


}// namespace gismo
