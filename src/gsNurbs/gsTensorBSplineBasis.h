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
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsTensorBSpline.h>

namespace gismo
{

/** 
    @brief A tensor product B-spline basis.

    \tparam T coefficient type
    \tparam d dimension of the parameter domain
    \tparam KnotVectorType type of the knot-vector

    \ingroup basis
    \ingroup Nurbs
*/
  
template<short_t d, class T>
class gsTensorBSplineBasis : public gsTensorBasis<d,T>
{
public: 
    typedef gsKnotVector<T> KnotVectorType;

    /// Base type
    typedef gsTensorBasis<d,T> Base;
    
    /// Family type
    typedef gsBSplineBasis<T>  Family_t;

    typedef gsTensorBSplineBasis Self_t;

    /// Coordinate basis type
    typedef gsBSplineBasis<T> CoordinateBasis;
    typedef CoordinateBasis                  Basis_t;

    /// Coefficient type
    typedef T Scalar_t;

    /// Associated Boundary basis type
    typedef typename gsBSplineTraits<d,T>::Geometry GeometryType;

    /// Associated Boundary basis type
    typedef typename gsBSplineTraits<static_cast<short_t>(d-1),T>::Basis BoundaryBasisType;

    typedef typename Base::iterator        iterator;
    typedef typename Base::const_iterator  const_iterator;

    /// Smart pointer for gsTensorBSplineBasis
    typedef memory::shared_ptr< Self_t > Ptr;

    /// Smart pointer for gsTensorBSplineBasis
    typedef memory::unique_ptr< Self_t > uPtr;

public:

    /// \brief Default constructor
    gsTensorBSplineBasis() : Base()
    {
        for(short_t i = 0; i!=d; ++i)
            this->m_bases[i] = new Basis_t();
    }

    void swap(gsTensorBSplineBasis & other)
    {
        this->Base::swap(static_cast<Base&>(other));
        std::swap(m_isPeriodic, other.m_isPeriodic);
    }
    
#if !EIGEN_HAS_RVALUE_REFERENCES
    gsTensorBSplineBasis & operator=(gsTensorBSplineBasis other)
    { gsTensorBSplineBasis::swap(other); return *this;}
#endif
    
    /**
       \brief Constructs a 2D tensor product B-spline basis. Assumes
       that the tamplate parameter \a d is equal to 2.
       
       \param KV1 knot-vector with respect to the first parameter dimension
       \param KV2 knot-vector with respect to the second parameter dimension
     */
    template<typename U>
    gsTensorBSplineBasis( KnotVectorType KV1, gsKnotVector<U> KV2,
                          typename util::enable_if<d==2,U>::type * = NULL )
    : Base( new Basis_t(give(KV1)), new Basis_t(give(KV2)) )
    { m_isPeriodic = -1; }

    /**
       \brief Constructs a 3D tensor product B-spline basis. Assumes
       that the tamplate parameter \a d is equal to 3.
       
       \param KV1 knot-vector with respect to the first dimension
       \param KV2 knot-vector with respect to the second dimension
       \param KV3 knot-vector with respect to the third dimension
     */
    gsTensorBSplineBasis( KnotVectorType KV1, 
                          KnotVectorType KV2, 
                          KnotVectorType KV3 )
    : Base( new Basis_t(give(KV1)), new Basis_t(give(KV2)), new Basis_t(give(KV3)) )
    { m_isPeriodic = -1; }

    gsTensorBSplineBasis( KnotVectorType KV1, 
                          KnotVectorType KV2, 
                          KnotVectorType KV3,
                          KnotVectorType KV4)
    : Base( new Basis_t(give(KV1)), new Basis_t(give(KV2)),
            new Basis_t(give(KV3)), new Basis_t(give(KV4)) )
    { m_isPeriodic = -1; }

    explicit gsTensorBSplineBasis(std::vector<KnotVectorType> KV)
    { 
        GISMO_ENSURE(d == KV.size(), "Invalid number of knot-vectors given." );
        for(short_t i = 0; i!=d; ++i)
            this->m_bases[i] = new Basis_t( give(KV[i]) );
        m_isPeriodic = -1; 
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
    
    gsTensorBSplineBasis(std::vector< gsBasis<T>*> & bb ) : Base(bb.data())
    {
        GISMO_ASSERT( checkVectorPtrCast<Basis_t>(bb), "Invalid vector of basis pointers.");
        GISMO_ENSURE( d == bb.size(), "Wrong d in the constructor of gsTensorBSplineBasis." );
        bb.clear();
        setIsPeriodic();
    }
    
    explicit gsTensorBSplineBasis(std::vector< Basis_t*> & bb ) 
    : Base( castVectorPtr<gsBasis<T> >(bb).data() )
    {
        GISMO_ENSURE( d == bb.size(), "Wrong d in the constructor of gsTensorBSplineBasis." );
        bb.clear();
        setIsPeriodic();
    }

#ifdef __DOXYGEN__
    /// \brief Returns the boundary basis for side s.
    typename gsBSplineTraits<static_cast<short_t>(d-1),T>::Basis::uPtr boundaryBasis(boxSide const & s);
#endif
    GISMO_UPTR_FUNCTION_DEF(BoundaryBasisType, boundaryBasis, boxSide const &)
    {
        std::vector<gsBasis<T>*> rr;
        this->getComponentsForSide(n1,rr);
        return BoundaryBasisType::New(rr);
    }

    GISMO_CLONE_FUNCTION(gsTensorBSplineBasis)
    
    static Self_t * New(std::vector<gsBasis<T>*> & bb )
    { return new Self_t(bb); }

    static Self_t * New(std::vector<Basis_t*> & bb )
    { return new Self_t(bb); }

    static uPtr make(std::vector<gsBasis<T>*> & bb )
    { return uPtr( new Self_t(bb) ); }

    static uPtr make(std::vector<Basis_t*> & bb )
    { return uPtr( new Self_t(bb) ); }
    
public:

    KnotVectorType & knots (int i)
    { return Self_t::component(i).knots(); }

    const KnotVectorType & knots (int i) const 
    { return Self_t::component(i).knots(); }

    // knot \a k of direction \a i
    T knot(int i, int k) const 
    { return Self_t::component(i).knots()[k]; }


    const Basis_t & component(short_t dir) const
    {
        return static_cast<const Basis_t &>(Base::component(dir));
    }

    Basis_t & component(short_t dir)
    {
        return static_cast<Basis_t &>(Base::component(dir));
    }

    // Look at gsBasis class for a description
    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const;

    /// Returns a box with the coordinate-wise active functions
    /// \param u evaluation points
    /// \param low lower left corner of the box
    /// \param upp upper right corner of the box
    void active_cwise(const gsMatrix<T> & u, gsVector<index_t,d>& low,
                      gsVector<index_t,d>& upp ) const;

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "TensorBSplineBasis: dim=" << this->dim()<< ", size="<< this->size() <<".";
        if( m_isPeriodic != -1 )
            os << "Periodic in " << m_isPeriodic << "-th direction.\n";
        for ( short_t i = 0; i!=d; ++i )
            os << "\n  Direction "<< i <<": "<< Self_t::component(i).knots() <<" ";
        os << "\n";
        return os;
    }

    //
    // param other parent/reference mesh determining the smoothness at the inner knots.
    // param i number of k-refinement steps to perform

    /// \brief Perform k-refinement coordinate-wise, in all directions.
    ///
    /// \copydetails gsBSplineBasis::refine_k
    void k_refine(Self_t & other, int const & i = 1)
    { 
        for (short_t j = 0; j < d; ++j)
            Self_t::component(j).refine_k(other.component(j), i);
    }

    /// \brief p-refinement (essentially degree elevation in all
    /// directions)
    void refine_p(int const & i = 1)
    {
        for (short_t j = 0; j < d; ++j)
            Self_t::component(j).refine_p(i);
    }
    
    /// \brief Uniform h-refinement (placing \a i new knots inside
    /// each knot-span, for all directions
    void refine_h(int const & i = 1)
    {
        for (short_t j = 0; j < d; ++j)
            Self_t::component(j).refine_h(i);
    }

    /**
     * \brief Takes a vector of coordinate wise knot values and
     * inserts these values to the basis.  
     *
     * Also constructs and returns
     * the transfer matrix that transfers coefficients to the new
     * basis.
     *
     * \param[out]  transfer Transfer matrix
     * \param[in] refineKnots Coordinate-wise knot values to be inserted
     */
    void refine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, const std::vector< std::vector<T> >& refineKnots);


    /**
     * \brief Takes a vector of coordinate wise knot values and
     * inserts these values to the basis.
     *
     * Also takes the old coefficients and changes them to reflect the new coefficients.
     *
     * \param[out]  coefs new coefficients
     * \param refineKnots Coordinate-wise knot values to be inserted
     *
     * \todo rename to insertKnots_withCoefs
     */
    void refine_withCoefs(gsMatrix<T> & coefs,const std::vector< std::vector<T> >& refineKnots);

    /// Inserts the knot \em knot with multiplicity \em mult in the knot
    /// vector of direction \a dir.
    void insertKnot(T knot, index_t dir, int mult=1)
    { this->knots(dir).insert( knot, mult); }

    /**
     * \brief Takes a vector of coordinate wise knot values and
     * inserts these values to the basis.
     *
     * \param refineKnots Coordinate-wise knot values to be inserted
     */
    void insertKnots(const std::vector< std::vector<T> >& refineKnots)
    {
        GISMO_ASSERT( refineKnots.size() == d, "refineKnots vector has wrong size" );
        for (short_t j = 0; j < d; ++j) // refine basis in each direction
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
     * \param[in] refExt is ignored
     *
     * \remarks NOTE This function directly modifies the basis (by inserting
     * knots in the underlying univariate B-spline bases).
     *
     * \ingroup Nurbs
     *
     *    \note: the \a refExt parameter is ignored in this implementation
     */
    void refine( gsMatrix<T> const & boxes, int refExt = 0);

    GISMO_MAKE_GEOMETRY_NEW

    /// \brief Reduces spline continuity (in all directions) at
    /// interior knots by \a i
    void reduceContinuity(int const & i = 1) 
    { 
        for (short_t j = 0; j < d; ++j)
            Self_t::component(j).reduceContinuity(i);
    }

    /// \brief Returns span (element) indices of the beginning and end
    /// of the support of the i-th basis function.
    template <int _Rows>
    void elementSupport_into(const index_t i, gsMatrix<index_t, _Rows, 2> & result) const
    {
        //result.resize(d,2);
        gsMatrix<index_t> tmp_vec;
        const gsVector<index_t, d> ti = this->tensorIndex(i);

        for (short_t dim = 0; dim < d; ++dim)
        {
            const gsKnotVector<T> & kv = Self_t::component(dim).knots();
            kv.supportIndex_into(ti[dim], tmp_vec);
            result.row(dim).noalias() =
                tmp_vec.cwiseMax(0).cwiseMin(kv.numElements());
        }
    }

    /// \brief Returns span (element) indices of the beginning and end
    /// of the support of the i-th basis function.
    gsMatrix<index_t, d, 2> elementSupport(const index_t & i) const
    {
        gsMatrix<index_t, d, 2> result(d, 2);
        elementSupport_into(i, result);
        return result;
    }

    /// \brief Returns the indices of active basis functions in the
    /// given input element box
    template <int _Rows>
    void elementActive_into(const gsMatrix<index_t,_Rows,2> & box,
                             gsMatrix<index_t> & result) const
    {
        GISMO_ASSERT( box.rows() == static_cast<index_t>(d), "Invalid input box");
        gsMatrix<index_t,d,2> tmp;
        
        gsVector<index_t,d> str;
        this->stride_cwise(str);

        for (short_t dm = 0; dm != d; ++dm)
        {
            tmp(dm,0) = Self_t::component(dm).knots().lastKnotIndex (box(dm,0)) - this->degree(dm);
            tmp(dm,1) = Self_t::component(dm).knots().firstKnotIndex(box(dm,1)) - 1;
        }

        const gsVector<index_t,d> sz = tmp.col(1)- tmp.col(0) + gsVector<index_t,d>::Ones();

        result = gsVector<index_t>::LinSpaced(sz[0], tmp(0,0), tmp(0,1));
        for (short_t dm = 1; dm != d; ++dm)
            result = result.replicate(1,sz[dm]) + 
                gsVector<index_t>::Constant(result.rows(), str[dm] )
                * gsVector<index_t>::LinSpaced(sz[dm], tmp(dm,0), tmp(dm,1))
                .transpose();
    }

    template <int _Rows>
    gsMatrix<T> elementDom(const gsMatrix<index_t,_Rows,2> & box) const
    {
        GISMO_ASSERT( box.rows() == static_cast<index_t>(d), "Invalid input box");
        gsMatrix<T> rvo;
        rvo.resize(d,2);
        for (short_t dm = 0; dm != d; ++dm)
        {
            rvo(dm,0) = Self_t::component(dm).knots().uValue(box(dm,0));
            rvo(dm,1) = Self_t::component(dm).knots().uValue(box(dm,1));
        }
        return rvo;
    }
    
    /// Tells, whether there is a coordinate direction in which the basis is periodic.
    inline bool isPeriodic() const { return m_isPeriodic != -1; }

    /// Gives the value of m_isPeriodic.
    inline int periodicDirection() const { return m_isPeriodic; }

    /// Converts \param dir -th basis to periodic.
    inline void setPeriodic( const int dir )
    {
        Self_t::component(dir).setPeriodic();
        if( Self_t::component(dir).isPeriodic() ) // Only when succeeded when converting to periodic.
            m_isPeriodic = dir;
    }

    /// Sets the coefficients so that the resulting TensorBSpline is periodic in direction dir.
    gsMatrix<T> perCoefs( const gsMatrix<T>& originalCoefs, short_t dir ) const
    {
        // Identify which coefficients to copy and where to copy them.
        std::vector<index_t> sourceSliceIndices;
        std::vector<index_t> targetSliceIndices;
        index_t numPeriodic = Self_t::component(dir).numCrossingFunctions();

        const index_t sz = this->size(dir) - numPeriodic;
        for( index_t i = 0; i < numPeriodic; i++ )
        {
            gsMatrix<index_t> currentSourceSlice = this->coefSlice(dir,i);
            gsMatrix<index_t> currentTargetSlice = this->coefSlice(dir, sz  + i );

            for( index_t j = 0; j != currentSourceSlice.size(); j++ )
            {
                sourceSliceIndices.push_back( static_cast<index_t>( currentSourceSlice(j) ) );
                targetSliceIndices.push_back( static_cast<index_t>( currentTargetSlice(j) ) );
            }
        }

        // Copy the chosen coefficients.
        gsMatrix<T> result = originalCoefs;
        for( size_t i = 0; i != sourceSliceIndices.size(); i++ )
        {
            //gsDebug << "source: " << sourceSliceIndices[i]  << "\n";
            //gsDebug << "target: " << targetSliceIndices[i]  << "\n";
            result.row( targetSliceIndices[ i ] ) = originalCoefs.row( sourceSliceIndices[ i ] );
        }

        return result;
    }

private:

    /// Repeated code from the constructors is held here.
    /// Sets m_isPeriodic to either -1 (if none of the underlying bases is periodic) or the index of the one basis that is periodic.
    void setIsPeriodic()
    {
        m_isPeriodic = -1;
        for( short_t i = 0; i < this->dim(); i++ )
        {
            if( Self_t::component(i).isPeriodic() )
            {
                if( m_isPeriodic == -1 )
                    m_isPeriodic = i;
                else
                    gsWarn << "Cannot handle a basis that is periodic in more than one direction.\n";
            }
        }
    }

protected:

    /// Coordinate direction, where the basis is periodic (when equal
    /// to -1 if there is no such direction).
    short_t m_isPeriodic;

};


} // namespace gismo


// *****************************************************************
#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorBSplineBasis.hpp)
/*
#else
#ifdef gsTensorBSplineBasis_EXPORT
#include GISMO_HPP_HEADER(gsTensorBSplineBasis.hpp)
#undef  EXTERN_CLASS_TEMPLATE
#define EXTERN_CLASS_TEMPLATE CLASS_TEMPLATE_INST
#endif
namespace gismo
{
EXTERN_CLASS_TEMPLATE gsTensorBSplineBasis<2,real_t>;
EXTERN_CLASS_TEMPLATE gsTensorBSplineBasis<3,real_t>;
EXTERN_CLASS_TEMPLATE gsTensorBSplineBasis<4,real_t>;
}
//*/
#endif
// *****************************************************************
