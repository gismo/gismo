/** @file gsPatchGenerator.h

    @brief Provides an interface for patch generators.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include<gsCore/gsMultiPatch.h>

namespace gismo
{


/**
   \brief Abstract class that accepts a set of input boundaries and
   computes a new geometry

   \tparam T Coefficient type

   \ingroup Modeling
*/
template <typename T>
class gsPatchGenerator
{
public:

    gsPatchGenerator() : m_result(NULL)
    { }

    /**
       \brief Constructs a patch generator object by a collection of
       geometries defining the boundaries of a patch.  
       
       \param boundary a set of boundary curves or patches
    */
    gsPatchGenerator(const gsMultiPatch<T> & boundary)
    : m_boundary(boundary), m_result(NULL)
    {
        GISMO_ASSERT( static_cast<int>(m_boundary.nPatches()) == 2*(m_boundary.parDim()+1), 
                      "Expecting "<<2*(m_boundary.parDim()+1)
                      <<" boundaries, received "<< m_boundary.nPatches()<<".");
    }

    virtual ~gsPatchGenerator()
    {
        delete m_result;
    }

public:

    /// \brief Main routine that performs the computation (to be
    /// implemented in derived classes)
    virtual const gsGeometry<T> & compute() = 0;

    /// \brief Main routine that receives input and performs the computation
    const gsGeometry<T> & compute(const gsMultiPatch<T> & boundary)
    {
        m_boundary = boundary;
        return compute();
    }
    
    /// \brief Returns the resulting patch. Assumes that compute() has
    /// been called before.
    const gsGeometry<T> & result() const
    { 
        GISMO_ASSERT(m_result != NULL, "Result not computed.");
        return *m_result; 
    }

    /// \brief Returns the input boundaries
    const gsMultiPatch<T> & input() 
    {
        return m_boundary;
    }

protected:

    /// \brief Resolves the configuration of the input boundaries and
    /// creates a patch filled with the boundary coefficients
    template<short_t d> void preparePatch(gsTensorBSplineBasis<d,T> & resultBasis,
                                           gsMatrix<T> & coefs);

protected:

    /// Input boundaries
    gsMultiPatch<T> m_boundary;
    
    /// Resulting patch
    gsGeometry<T> * m_result;
    

}; // gsPatchGenerator



template <typename T> template <short_t d>
void gsPatchGenerator<T>::preparePatch(gsTensorBSplineBasis<d,T> & resultBasis, gsMatrix<T> & coefs)
{
    GISMO_ASSERT(m_boundary.nPatches()  == 2*d, 
                 "Expecting "<<2*d<<" boundaries");

    typedef typename gsBSplineTraits<static_cast<short_t>(d-1),T>::Geometry Boundary_t;

    //-------- 1. Find the pairs of facing boundaries

    // A. Compute the topology of the input boundaries based on corners
    m_boundary.computeTopology(1e-3);
    GISMO_ASSERT(m_boundary.nBoundary()==0,"The input boundary is not closed.");
    //gsDebugVar( m_boundary.detail() );

    // Make sure that the shell is water-tight (remove small gaps)
    m_boundary.closeGaps(1e-3);

    // B. Permute boundaries so that we get a sequence of facing
    //    geometries: 0-1, 2-3, 4-5, ..

    std::vector<Boundary_t*> input;// permutation of the input boundaries
    input.reserve(2*d);
    std::vector<short_t> perm;// permutation to be computed
    perm.reserve(2*d);

    for (unsigned k = 0; k!=2*d; ++k) //for all boundaries
    {
        for (unsigned l = k+1; l!=2*d; ++l)
            if ( (k!=l) && (m_boundary.findInterface(k,l)==NULL) )
            {
                input.push_back( dynamic_cast<Boundary_t*>(&m_boundary.patch(k)) );
                perm.push_back(k);
                GISMO_ASSERT(input.back()!=NULL,
                             "Could not convert to the expected geometry type.");
                input.push_back( dynamic_cast<Boundary_t*>(&m_boundary.patch(l)) );
                perm.push_back(l);
                GISMO_ASSERT(input.back()!=NULL,
                             "Could not convert to the expected geometry type.");
                break;
            }
    }

    //gsDebugVar(gsAsMatrix<int>(perm));//debug
    GISMO_ASSERT(input.size()==2*d,"Pairs not correctly identified.");    

    //-------- 2. Resolve the configuration: Force orientation
    //    according to the positive orthant. Origin (0,..,0) is
    //    pinched to the first vertex of boundary 0, and all other
    //    boundaries are oriented accordingly.

    // A. Origin corner

    gsMatrix<T> v = input[0]->coef(0);
    for (unsigned k = 2; k<2*d; k+=2) //for all pairs but the first
    {
        // Find joint with direction k/2
        if ( ! input[k]->isPatchCorner(v) )
        {
            std::swap(input[k], input[k+1]);
            std::swap(perm [k], perm [k+1]);
            GISMO_ASSERT(input[k]->isPatchCorner(v),"No joint found");
        }
        // Set the origin of the parametrization to v
        input[k]->setOriginCorner(v);
    }

    // Apply the same permutation to m_boundary
    m_boundary.permute(perm);
    //gsDebugVar(gsAsMatrix<int>(perm));

    // B. Furthest corner
    int count = 0;// becomes true if a common corner is found for all odd boundaries
    for (boxCorner c = boxCorner::getFirst(d-1); c<boxCorner::getEnd(d-1); ++c)
    {
        v = input[1]->coefAtCorner(c);
        count = 1;
        for (unsigned k = 3; k<2*d; k+=2)
            if ( input[k]->isPatchCorner(v) )
                ++count;
        
        if (count == d)
            break;
    }

    GISMO_ASSERT(count==d, "Furthest corner not found");
    // Set furthest corner of all parametrizations to v
    for (unsigned k = 1; k<2*d; k+=2)
        input[k]->setFurthestCorner(v); 
    
    // C. Fix orientation
    m_boundary.computeTopology();
    //gsDebugVar( m_boundary.detail() );
    GISMO_ASSERT(m_boundary.nBoundary()==0, "Something went wrong with boundary identification.");

    if ( d>2)
    {
        std::fill(perm.begin(),perm.end(),0);
        for (unsigned k = 0; k!=2*(d-1); ++k) //for all pairs
        {
            if ( typename gsMultiPatch<T>::InterfacePtr bi = m_boundary.findInterface(k,k+2) )
            {
                if (bi->first() .side() - k-1 != 0 )//side!=k+1
                    perm[k]   = 1;
                if (bi->second().side() - k-1 != 0 )//side!=k+1
                    perm[k+2] = 1;
            }
            else
            {
                gsWarn<<"Topology is not the expected one.";
            }
        }

        for (unsigned k = 0; k!=2*d; ++k)
            if ( perm[k] )
                input[k]->swapDirections(0,1);

        m_boundary.computeTopology();
        // gsDebugVar(gsAsMatrix<int>(perm));
        // gsDebugVar( m_boundary.detail() );

        GISMO_ASSERT(m_boundary.nBoundary()==0,
                     "Something went wrong with boundary identification.");
    }

    //-------- 3. At this point the configuration of the boundaries is
    //  correct.  Decide on a common knot vector and degree for each
    //  coordinate.

    std::vector<gsBSplineBasis<T>*> cBases(d);
    // Issues: 1. find the correct knot-vectors per direction
    //         2. do the input patches meet in the right way at facets ? flips needed ?

    // A. Check initial compatibility of knot-vectors

    for (unsigned k = 0; k<2*d; k+=2) //for all pairs
    {
        for (unsigned l = 0; l!=d-1; ++l)
        {
            if ( input[k]->basis().degree(l) !=  input[k+1]->basis().degree(l) )
                // to do
                gsWarn<< "*** Degrees: "<< input[k]->basis().degree(l) <<", "
                      << input[k+1]->basis().degree(l) << "differ.\n";
        
            if ( input[k]->knots(l) !=  input[k+1]->knots(l) )
            {
                // to do
                gsWarn<< "*** Knots differ.\n";
                //1. Normalize both and check again
                //2. Merge
            }
        }
    }

    // B. Try to reorder the input so that we get a compatible ordering for in dD.

    for (unsigned k = 0; k!=d-1; ++k) // For now use the input knot-vectors
    {
        cBases[k] = input[2*d-2]->basis().component( k ).clone().release();//
    }
    // Last one
    cBases[d-1] = input.front()->basis().component(d-2).clone().release();

/* For debugging
    for (unsigned k = 0; k!=d; ++k)
    {
        gsDebug<<"*** Boundaries "<<k<<":\n";
        gsDebugVar( input[2*k  ]->basis() );
        gsDebugVar( input[2*k+1]->basis() );
    }
//*/

    // Create the final tensor basis
    resultBasis = gsTensorBSplineBasis<d,T>(cBases); //Note: constructor consumes pointers

    //gsDebug<<"*** Choice "<<resultBasis<<"\n";

    //-------- 4. Fill in the boundary of the patch

    gsMatrix<index_t> bdr; // indices of the boundary control points
    coefs.setZero(resultBasis.size(), m_boundary.geoDim());
    
    // Fill in boundary coefficients
    for ( short_t i = 0; i<d; i++ )
    {
        for ( int s = 0; s<2; s++ )
        {
            // Get the input coefficients for this boundary
            const gsMatrix<T> & bCoefs = input[2*i+s]->coefs();
            
            //Get the boundary indices ("1+" is because boxSides start from 1)
            bdr = resultBasis.boundary(static_cast<boxSide>(1+2*i+s));

            GISMO_ASSERT( bdr.rows() == bCoefs.rows(),
                          "Wrong dim "<<s <<": "<<bdr.rows()<<" != "<<bCoefs.rows() 
                          <<" at patch "<< 2*i+s <<", side "<<static_cast<boxSide>(1+2*i+s)
                          <<".\n Component:\n"<< input[2*i+s]->basis()
                          << "\n ResultBasis:\n"<< resultBasis//.component(2*i+s)
                          );
            
            // Fill in resulting patch CPs
            for ( index_t k = 0; k< bCoefs.rows(); k++ )
                coefs.row( bdr(k,0) ) = bCoefs.row(k);
        }
    }

}


}// namespace gismo
