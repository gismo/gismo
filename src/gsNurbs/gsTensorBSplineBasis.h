/** @file gsTensorBSplineBasis.h

    @brief Provides declaration of TensorBSplineBasis abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsTensor/gsTensorBasis.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsTensor/gsTensorTools.h>

namespace gismo
{

// forward declaration
template<unsigned d, class T, class KnotVectorType>
class gsTensorBSplineBasis;


/// Traits for TensorBSplineBasis
template<unsigned d, class T, class KnotVectorType>
struct gsTraits< gsTensorBSplineBasis<d,T,KnotVectorType>, d>
{
    typedef typename 
    gsTraits<gsBSplineBasis<T,KnotVectorType>,d>::RationalBasisType
    RationalBasisType;
  
    typedef typename 
    gsTraits<gsBSplineBasis<T,KnotVectorType>,d>::RationalGeometryType
    RationalGeometryType;  

    typedef typename 
    gsTraits<gsBSplineBasis<T,KnotVectorType>,d-1>::RationalBasisType
        RationalBoundaryType;
};

/** 
    @brief A tensor product B-spline basis.

    \param T coefficient type
    \param d dimension of the parameter domain

    \ingroup basis
    \ingroup Nurbs
*/

  
template<unsigned d, class T, class KnotVectorType>
class gsTensorBSplineBasis : public gsTensorBasis< d, gsBSplineBasis<T,KnotVectorType> >  
{

public: 
    /// Base type
    typedef gsTensorBasis< d, gsBSplineBasis<T,KnotVectorType> > Base;

    /// Coordinate basis type
    typedef gsBSplineBasis<T,KnotVectorType> Basis_t;

    /// Coefficient type
    typedef T Scalar_t;

    /// Associated geometry type
    typedef gsTensorBSpline<d, T, KnotVectorType> GeometryType;

    /// Associated Boundary basis type
    //typedef typename Base::BoundaryBasisType BoundaryBasisType;

    typedef typename Base::iterator        iterator;
    typedef typename Base::const_iterator  const_iterator;

    /// Shared pointer for gsTensorBasis
    typedef memory::shared_ptr< gsTensorBSplineBasis > Ptr;

public:
    /// Constructors for gsTensorBSplineBasis
    gsTensorBSplineBasis( const KnotVectorType& KV1, const KnotVectorType& KV2 )
    : Base( new Basis_t(KV1), new Basis_t(KV2) )
    { m_isPeriodic = -1; }

    // TO DO: temporary, REMOVE
    gsTensorBSplineBasis( const gsTensorBasis<d,gsBSplineBasis<T,KnotVectorType> > & bb)
    : Base( bb )
    { }

    gsTensorBSplineBasis( const KnotVectorType& KV1, 
                          const KnotVectorType& KV2, 
                          const KnotVectorType& KV3 )
    : Base( new Basis_t(KV1), new Basis_t(KV2), new Basis_t(KV3) )
    { m_isPeriodic = -1; }

    // TO DO: more constructors



    // Constructors forwarded from the base class
    gsTensorBSplineBasis() : Base() { }

    // TODO: Can we store the error message somewhere centrally to be sure that it is the same for all the constructors?

    explicit gsTensorBSplineBasis( Basis_t* x) : Base(x) 
    {
        GISMO_ENSURE(d==1,"Invalid constructor." );
        setIsPeriodic();
    }

    gsTensorBSplineBasis( Basis_t* x,  Basis_t*  y) : Base(x,y) 
    { 
        GISMO_ENSURE(d==2,"Invalid constructor." );
        setIsPeriodic();
    }
    
    gsTensorBSplineBasis( Basis_t* x,  Basis_t* y, Basis_t* z ) : Base(x,y,z) 
    { 
        GISMO_ENSURE(d==3,"Invalid constructor." );
        setIsPeriodic();
        
    }
    
    gsTensorBSplineBasis( Basis_t* x,  Basis_t* y, Basis_t* z, Basis_t* w ) : Base(x,y,z,w) 
    { 
        GISMO_ENSURE(d==4,"Invalid constructor." );
        setIsPeriodic();
    }
    
    
    gsTensorBSplineBasis( std::vector<Basis_t* > const & bb ) : Base(bb)
    {
        GISMO_ENSURE( d == bb.size(), "Wrong d in the constructor of gsTensorBSplineBasis." );
        setIsPeriodic();
    }
    
    gsTensorBSplineBasis( const gsTensorBSplineBasis & o) : Base(o)
    {
        m_isPeriodic = o.periodicDirection();
    }
    
    gsTensorBSplineBasis * clone() const
    { return new gsTensorBSplineBasis(*this); }
    
public:

    KnotVectorType & knots (int i)
    { return this->component(i).knots(); }

    const KnotVectorType & knots (int i) const 
    { return this->component(i).knots(); }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "TensorBSplineBasis: dim=" << this->dim()<< ", size="<< this->size() <<".";
        if( m_isPeriodic != -1 )
            os << "Periodic in " << m_isPeriodic << "-th direction.\n";
        for ( unsigned i = 0; i!=d; ++i )
            os << "\n  Direction "<< i <<": "<< this->component(i).knots() <<" ";
        os << "\n";
        return os;
    }

    /**
     * \brief Perform k-refinement coordinate-wise
     *
     * \param[in] i number of k-refinement steps to perform
     */
    void k_refine(int const & i = 1) 
    { 
        for (unsigned j = 0; j < d; ++j)
            this->m_bases[j]->k_refine(i);
    }

    /**
     * \brief
     * Takes a vector of coordinate wise knot values and inserts these values to the basis.
     * Also constructs and returns the transfer matrix that transfers coefficients to the new basis.
     *
     * \param u     refineKnots Coordinate-wise knot values to be inserted
     * \param[out]  transfer Transfer matrix
     */
    void refine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, const std::vector< std::vector<T> >& refineKnots)
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

    /**
     * \brief
     * Takes a vector of coordinate wise knot values and inserts these values to the basis.
     * Also takes the old coefficients and changes them to reflect the new coefficients.
     *
     * \param u     refineKnots Coordinate-wise knot values to be inserted
     * \param[out]  new coefficients
     \todo rename to insertKnots_withCoefs
     */
    void refine_withCoefs(gsMatrix<T> & coefs,const std::vector< std::vector<T> >& refineKnots)
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

    /**
     * \brief
     * Takes a vector of coordinate wise knot values and inserts these values to the basis.
     *
     * \param u     refineKnots Coordinate-wise knot values to be inserted
     */
    void insertKnots(const std::vector< std::vector<T> >& refineKnots)
    {
        GISMO_ASSERT( refineKnots.size() == d, "refineKnots vector has wrong size" );
        for (unsigned j = 0; j < d; ++j) // refine basis in each direction
            this->knots(j).insert(refineKnots[j]);
    }

    /** \brief Refinement of the tensor basis on the area defined by \em boxes.
     *
     * Applies "local" refinement within the tensor-product structure of
     * gsTensorBSplineBasis. The areas for refinement are specified in \em boxes.\n
     * \em boxes is a gsMatrix of size <em>d</em> x <em>(2*N)</em>, where\n
     * \em d is the dimension of the parameter domain and\n
     * \em N is the number of refinement boxes.\n
     * \n
     * Every two successive columns in \em boxes correspond to the coordinates
     * of the lower and upper corners of one refinement box, respectively (see example below).
     * Note that a new knot will be inserted in every knot span contained in this area.
     * If some of the given boxes overlap, the refinement will only be done once.
     *
     * <b>Example</b>, let
     * \verbatim
     d = 2
     knotvector1 = knotvector2 = [ 0 0 0  0.25  0.5  0.75  1 1 1 ]

     boxes = [ 0.25  0.75  0     0.5 ]
     [ 0     0.25  0.75  1   ]
     \endverbatim
     * The areas
     * <em>[ 0.25, 0.75 ] x [ 0, 0.25 ]</em> and
     * <em>[ 0, 0.5 ] x [ 0.75, 1 ]</em> will be refined. The knots \em 0.125, \em 0.375, and \em 0.625 will
     * be inserted in \em knotvector1, and the knots \em 0.125 and \em 0.875 in \em knotvector2.
     *
     * \param[in] boxes gsMatrix of size <em>d</em> x <em>(2*N)</em>;
     * specifies areas for refinement.\n
     * See above for details and format.
     * \param NOTE This function directly modifies the basis (by inserting
     * knots in the underlying univariate B-spline bases).
     *
     * \ingroup Nurbs
     */
    void refine( gsMatrix<T> const & boxes )
    {
        GISMO_ASSERT( boxes.rows() == this->dim() , "Number of rows of refinement boxes must equal dimension of parameter space.");
        GISMO_ASSERT( boxes.cols() % 2 == 0, "Refinement boxes must have even number of columns.");

        const T tol = 0.000000001;

        // for each coordinate direction of the parameter domain:
        for( int di = 0; di < this->dim(); di++)
        {

            // for simplicity, get the corresponding knot vector.
            KnotVectorType kold_di = this->component(di).knots();

            // vector of flags for refining knotspans
            gsVector<int> flagInsertKt( kold_di.size() );
            flagInsertKt.setZero();

            // This is a very crude and unelegant test, but
            // conveniently straightforward to implement (and maybe to):
            // For each knot span, check if its midpoint is
            // contained in any of the refinement boxes.
            // If yes, set the corresponding flag to 1.
            for( int i=1; i < kold_di.size(); i++ ) // loop over spans
                if( kold_di[i]-kold_di[i-1] > tol)  // check for empty spans
                {
                    T midpt = (kold_di[i] + kold_di[i-1])/2; // midpoint of knot span
                    for( int j=0; j < boxes.cols(); j+=2 ) // loop over all boxes
                    {
                        if( boxes(di,j) < midpt && midpt < boxes(di,j+1) )
                            flagInsertKt[i] = 1; // if the box contains the midpoint, mark it
                    }
                }

            // now, with the flags set, loop over all knots spans and
            // insert midpoint in each knot span which is marked for refinement.
            for( int i=1; i < kold_di.size(); i++ )
                if( flagInsertKt[i] == 1)
                {
                    T midpt = (kold_di[i] + kold_di[i-1])/2;
                    this->component(di).insertKnot( midpt );
                }

        } // for( int di )


    } // refine()

    GISMO_MAKE_GEOMETRY_NEW

    /// Reduces spline continuity (in all directions) at interior knots by \a i
    void reduceContinuity(int const & i = 1) 
    { 
        for (unsigned j = 0; j < d; ++j)
            this->component(j).reduceContinuity(i);
    }

    /// Returns knot indices of the beginning and end of the support of the i-th
    /// basis function.
    void elementSupport_into(const unsigned& i,
                             gsMatrix<unsigned, d, 2>& result) const
    {
        gsMatrix<unsigned> tmp_vec;
        gsVector<unsigned, d> ti = this->tensorIndex(i);

        for (unsigned dim = 0; dim < d; ++dim)
        {
            this->component(dim).knots().supportIndex_into(ti[dim], tmp_vec);
            result.row(dim) = tmp_vec.row(0);
        }
    }


    /// Returns knot indices of the beginning and end of the support of the i-th
    /// basis function.
    gsMatrix<unsigned, d, 2> elementSupport(const unsigned & i) const
    {
        gsMatrix<unsigned, d, 2> result(d, 2);
        elementSupport_into(i, result);
        return result;
    }

    /// Tells, whether there is a coordinate direction in which the basis is periodic.
    inline bool isPeriodic() const { return m_isPeriodic != -1; }

    /// Gives the value of m_isPeriodic.
    inline int periodicDirection() const { return m_isPeriodic; }

    /// Converts \param dir -th basis to periodic.
    inline void setPeriodic( const int dir )
    {
        this->component(dir).setPeriodic();
        if( this->component(dir).isPeriodic() ) // Only when succeeded when converting to periodic.
            m_isPeriodic = dir;
    }

    /// Sets the coefficients so that the resulting TensorBSpline is periodic in direction dir.
    gsMatrix<T> perCoefs( const gsMatrix<T>& originalCoefs, int dir ) const
    {
        // Identify which coefficients to copy and where to copy them.
        std::vector<index_t> sourceSliceIndices;
        std::vector<index_t> targetSliceIndices;
        int numPeriodic = this->component(dir).numCrossingFunctions();
        for( int i = 0; i < numPeriodic; i++ )
        {
            gsMatrix<unsigned> currentSourceSlice = *(this->slice(dir,i));
            gsMatrix<unsigned> currentTargetSlice = *(this->slice(dir, this->size(dir)  + i ));

            for( index_t j = 0; j < currentSourceSlice.size(); j++ )
            {
                sourceSliceIndices.push_back( static_cast<index_t>( currentSourceSlice(j) ) );
                targetSliceIndices.push_back( static_cast<index_t>( currentTargetSlice(j) ) );
            }
        }

        // Copy the chosen coefficients.
        gsMatrix<T> result = originalCoefs;
        for( std::size_t i = 0; i < sourceSliceIndices.size(); i++ )
        {
            //std::cout << "source: " << sourceSliceIndices[i]  << "\n";
            //std::cout << "target: " << targetSliceIndices[i]  << "\n";
            result.row( targetSliceIndices[ i ] ) = originalCoefs.row( sourceSliceIndices[ i ] );
        }

        return result;
    }

    /// Evaluate an element of the space given by coefs at points u
    void eval_into_new(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const
    {

    }


private:

    /// Repeated code from the constructors is held here.
    /// Sets m_isPeriodic to either -1 (if none of the underlying bases is periodic) or the index of the one basis that is periodic.
    void setIsPeriodic()
    {
        m_isPeriodic = -1;
        for( int i = 0; i < this->dim(); i++ )
        {
            if( this->component(i).isPeriodic() )
            {
                if( m_isPeriodic == -1 )
                    m_isPeriodic = i;
                else
                    gsWarn << "Cannot handle a basis that is periodic in more than one direction.\n";
            }
        }
    }

    /////////////
    // Members //
    /////////////

    // Coordinate direction, where the basis is periodic (equal to -1 if there is no such direction).
    int m_isPeriodic;

};


} // namespace gismo




//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorBSplineBasis.hpp)
#endif
