/** @file gsTensorBSplineBasis.hpp

    @brief Provides declaration of tensor-product B-spline basis class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsTemplateTools.h>
#include <gsTensor/gsTensorTools.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsBoehm.h>

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>


namespace gismo
{


template<short_t d, class T>
void gsTensorBSplineBasis<d,T>::
active_cwise(const gsMatrix<T> & u, 
             gsVector<index_t,d>& low,
             gsVector<index_t,d>& upp ) const
{
    for (index_t j = 0; j < u.cols(); ++j)
    {
        for (short_t i = 0; i < d; ++i)
        {
            low[i] = component(i).firstActive( u(i,j) );
            upp[i] = low[i] + component(i).degree();
        }
    }
}

template<short_t d, class T>
void gsTensorBSplineBasis<d,T>::
refine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, 
                    const std::vector< std::vector<T> >& refineKnots)
{
    GISMO_ASSERT( refineKnots.size() == d, "refineKnots vector has wrong size" );
    gsSparseMatrix<T,RowMajor> B[d];
    
    // refine component bases and obtain their transfer matrices
    for (unsigned i = 0; i < d; ++i)
    {
        this->component(i).refine_withTransfer( B[i], refineKnots[i] );
    }
    
    tensorCombineTransferMatrices<d, T>( B, transfer );
}

template<short_t d, class T>
void gsTensorBSplineBasis<d,T>::
refine_withCoefs(gsMatrix<T> & coefs,const std::vector< std::vector<T> >& refineKnots)
{
    GISMO_ASSERT( refineKnots.size() == d, "refineKnots vector has wrong size" );
    gsVector<unsigned> strides(d);
    for (unsigned j = 0; j < d; ++j) //calculate the strides for each direction
    {
        strides[j]=this->stride(j);
    }
    for (unsigned i = 0; i < d; ++i)
    {
        if(refineKnots[i].size()>0)
        {
            gsTensorBoehmRefine(this->component(i).knots(), coefs, i, strides,
                                refineKnots[i].begin(), refineKnots[i].end(), true);
            
            for (index_t j = i+1; j<strides.rows(); ++j)
                strides[j]=this->stride(j); //new stride for this direction
        }
    }
}


template<short_t d, class T>
void gsTensorBSplineBasis<d,T>::refine(gsMatrix<T> const & boxes, int)
{
    // Note: refExt parameter is ignored
    GISMO_ASSERT( boxes.rows() == this->dim() , 
                  "Number of rows of refinement boxes must equal dimension of parameter space.");
    GISMO_ASSERT( boxes.cols() % 2 == 0, 
                  "Refinement boxes must have even number of columns.");
    
    const T tol = 0.000000001;
    
    // for each coordinate direction of the parameter domain:
    for( int di = 0; di < this->dim(); di++)
    {
        
        // for simplicity, get the corresponding knot vector.
        KnotVectorType kold_di = Self_t::component(di).knots();

        // vector of flags for refining knotspans
        std::vector<bool> flagInsertKt( kold_di.size(), false);

        // This is a very crude and unelegant test, but
        // conveniently straightforward to implement (and maybe to):
        // For each knot span, check if its midpoint is
        // contained in any of the refinement boxes.
        // If yes, set the corresponding flag to 1.
        for( size_t i=1; i < kold_di.size(); i++ ) // loop over spans
            if( kold_di[i]-kold_di[i-1] > tol)  // check for empty spans
            {
                const T midpt = (kold_di[i] + kold_di[i-1])/2; // midpoint of knot span
                for( index_t j=0; j < boxes.cols(); j+=2 )     // loop over all boxes
                {
                    // if the box contains the midpoint, mark it
                    flagInsertKt[i] = boxes(di,j) < midpt && midpt < boxes(di,j+1);
                }
            }

        // now, with the flags set, loop over all knots spans and
        // insert midpoint in each knot span which is marked for refinement.
        for( size_t i=1; i < kold_di.size(); i++ )
            if( flagInsertKt[i] )
            {
                T midpt = (kold_di[i] + kold_di[i-1])/2;
                Self_t::component(di).insertKnot( midpt );
            }

    } // for( int di )


} // refine()


/*
 * NB:
 * This implementation assumes that the component bases implement the firstActive() / numActive() protocol.
 * In particular, their active basis functions must always be continuous intervals.
 * This is the case for all current component bases, so we only keep this version for now.
 * Above, commented out, is the generic version which is quite a bit slower.
 */
template<short_t d, class T>
void gsTensorBSplineBasis<d,T>::
active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const
{
    GISMO_ASSERT( u.rows() == static_cast<index_t>(d), "Invalid point dimension: "<<u.rows()<<", expected "<< d);

    unsigned firstAct[d];
    gsVector<unsigned, d> v, size;

    // count active functions in each tensor direction
    unsigned numAct = 1;
    for (unsigned i = 0; i < d; ++i)
    {
        size[i] = component(i).numActive();
        numAct *= size[i];
    }

    result.resize( numAct, u.cols() );  
  
    // Fill with active bases indices
    for (index_t j = 0; j < u.cols(); ++j)
    {
        // get the active basis indices for the component bases at u(:,j)
        for (short_t i = 0; i < d; ++i)
        {
            firstAct[i] = component(i).firstActive( u(i,j) );
        }

        // iterate over all tensor product active functions
        unsigned r = 0;
        v.setZero();
        do
        {
            index_t gidx = firstAct[d-1] + v(d-1);    //compute global index in the tensor product
            for ( short_t i=d-2; i>=0; --i )
                gidx = gidx * this->size(i) + firstAct[i] + v(i);

            result(r, j) = gidx;
            ++r ;
        } while (nextLexicographic(v, size));
    }
}


namespace internal
{

/// @brief Get a TensorBSplineBasis from XML data
template<short_t d, class T>
class gsXml< gsTensorBSplineBasis<d,T> >
{
private:
    gsXml() { }
    typedef gsTensorBSplineBasis<d,T> Object;
public:
    GSXML_COMMON_FUNCTIONS(Object);
    static std::string tag () { return "Basis"; }
    static std::string type () { return "TensorBSplineBasis"+to_string(d); }

    static Object * get (gsXmlNode * node)
    {
        if (d>1)
        {
            GISMO_ASSERT( ( !strcmp( node->name(),"Basis") )
                          &&  ( !strcmp(node->first_attribute("type")->value(),
                                        internal::gsXml<Object>::type().c_str() ) ),
                          "Something is wrong with the XML data: There should be a node with a "
                          <<internal::gsXml<Object>::type().c_str()<<" Basis.");
        }
        else
        {
            GISMO_ASSERT( ( !strcmp( node->name(),"Basis") )
                          &&  ( !strcmp(node->first_attribute("type")->value(),
                                        internal::gsXml<Object>::type().c_str())
                                ||  !strcmp(node->first_attribute("type")->value(), "BSplineBasis") ),
                          "Something is wrong with the XML data: There should be a node with a "
                          <<internal::gsXml<Object>::type().c_str()<<" Basis.");
        }

        return getTensorBasisFromXml<Object >( node );
    }
    
    static gsXmlNode * put (const Object & obj, 
                            gsXmlTree & data )
    {
        return putTensorBasisToXml<Object >(obj,data);
    }
};

} // internal

} // gismo

