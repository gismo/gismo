/** @file gsBSplineBasis.h

    @brief Represents a B-spline basis with one parameter

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D. Mokris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsConstantBasis.h>

#include <gsTensor/gsTensorDomainIterator.h>
#include <gsTensor/gsTensorDomainBoundaryIterator.h>

#include <gsNurbs/gsBSplineAlgorithms.h>

namespace gismo
{

/// \brief Traits for BSplineBasis in more dimensions
template<unsigned d, class T, class KnotVectorType>
struct gsTraits<gsBSplineBasis<T, KnotVectorType>,d>
{
    typedef gsBSplineBasis<T, KnotVectorType> Family_t;
    typedef gsBSpline<T,KnotVectorType>       GeometryType;
    typedef gsConstantBasis<T>                BoundaryBasisType;

    typedef gsTensorBSplineBasis<d,T,KnotVectorType>  TensorBasisType;
    typedef gsTensorBSpline<d, T, KnotVectorType>     TensorGeometryType;
    typedef typename 
    choose<d==1, gsConstantBasis<T>, gsTensorBSplineBasis<d-1,T,KnotVectorType>
           >::type TensorBoundaryType;
    typedef typename 
    choose<d==1, gsConstantFunction<T>, gsTensorBSpline<d-1,T,KnotVectorType>
           >::type TBoundaryGeometryType;

    typedef typename 
    choose<d==1, gsNurbsBasis<T,KnotVectorType>, gsTensorNurbsBasis<d, T, KnotVectorType>
           >::type RationalBasisType;
    typedef typename 
    choose<d==1, gsNurbs<T,KnotVectorType>, gsTensorNurbs<d, T, KnotVectorType>
           >::type RationalGeometryType;
    typedef typename 
    choose<d==1, gsConstantBasis<T>       , gsTensorNurbsBasis<d-1, T, KnotVectorType>
           >::type RationalBoundaryType;
};

/** \brief
    A univariate B-spline basis.

    \tparam T coefficient type
    \tparam KnotVectorType the type of knot vector to use

    \ingroup basis
*/

template<class T, class KnotVectorType >
class gsBSplineBasis : public gsBasis<T>
{
public:
    typedef gsBasis<T> Base;

    typedef gsBSplineBasis<T,KnotVectorType> Self_t;

    /// This object identifies the family of B-spline bases
    typedef Self_t Family_t;

    /// Coefficient type
    typedef T Scalar_t;

    /// Associated geometry type
    typedef typename gsTraits<Family_t,1>::GeometryType GeometryType;

    /// Associated Boundary basis type
    typedef typename gsTraits<Family_t,1>::BoundaryBasisType BoundaryBasisType;

    /// Dimension of the parameter domain
    static const int Dim = 1;

    /// Shared pointer for gsBSplineBasis
    typedef memory::shared_ptr< gsBSplineBasis > Ptr;

    static Ptr makeShared ( const KnotVectorType & KV )
    { return Ptr( new gsBSplineBasis(KV) ); }
  
public:

    /// Default empty constructor
    explicit gsBSplineBasis( bool periodic = false ) :
    Base(), m_p(0), m_periodic(0)
    {
        m_knots.initClamped(0);
        if( periodic )
        {
            gsWarn << "Converting your basis to periodic.\n";
            _convertToPeriodic();
        }

        if( ! check()  )
            gsWarn << "Warning: Insconsistent "<< *this<< "\n";
    }


    /// Construct BSpline basis of a knot vector
    gsBSplineBasis( const KnotVectorType & KV, bool periodic = false) :
    m_p(KV.degree()), m_knots(KV), m_periodic(0)
    { 
        if( periodic )
        {
            gsWarn << "Converting your basis to periodic.\n";
            _convertToPeriodic();
        }

        if( ! check()  )
            gsWarn << "Warning: Insconsistent "<< *this<< "\n";

    }

    /// Construct a BSpline basis
    /// \param u0 starting parameter
    /// \param u1 end parameter parameter
    /// \param interior number of interior knots
    /// \param degree degree of the spline space
    /// \param mult_interior multiplicity at the interior knots
    gsBSplineBasis(T u0, T u1, unsigned interior, 
                   int degree, unsigned mult_interior=1,
                   bool periodic = false ) :
    m_p(degree), 
    m_knots(u0, u1, interior, m_p+1, mult_interior, m_p), 
    m_periodic(0)
    {
        if( periodic )
        {
            _convertToPeriodic();
            gsWarn << "Converting your basis to periodic.\n";
        }

        if( ! check()  )
            gsWarn << "Warning: Insconsistent "<< *this<< "\n";
    }

    
    /// Copy Constructor
    gsBSplineBasis( const gsBSplineBasis & o)
    : m_p(o.m_p), m_knots( o.knots() ), m_periodic(o.m_periodic)
    {

    }
    
    // gsBSplineBasis( const Base & o)
    // { 
    //     const gsBSplineBasis * a;
    //     if ( ( a = dynamic_cast<const gsBSplineBasis *>( &o )) )
    //     {
    //         m_p        = a->degree() ; 
    //         m_knots    = KnotVectorType( a->knots() );
    //         m_periodic = a->m_periodic;
    //     }
    //     else
    //         GISMO_ERROR("Cannot convert "<<o<<" to gsBSplineBasis\n");
    // }

public:

    void swap(gsBSplineBasis& other)
    {
        std::swap(m_p, other.m_p);
        std::swap(m_periodic, other.m_periodic);
        m_knots.swap(other.m_knots);
    }

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

    // Look at gsBasis class for a description
    int dim() const { return Dim; }

    // Look at gsBasis class for a description
    int size() const { return m_knots.size() - m_p - 1 - m_periodic; }
    
    // Look at gsBasis class for a description
    int numElements() const { return m_knots.numKnotSpans(); }

    // Look at gsBasis class for a description
    int elementIndex(const gsVector<T> & u ) const;

    // Same as gsBasis::elementIndex but argument is a value instead of a vector
    int elementIndex(T u ) const;

    // Look at gsBasis class for a description
    gsBSplineBasis & component(unsigned i) const;

    /// Returns the anchors (greville points) of the basis
    void anchors_into(gsMatrix<T> & result) const 
    { 
        m_knots.greville_into(result); 
    }

    /// Returns the anchors (greville points) of the basis
    void anchor_into(unsigned i, gsMatrix<T> & result) const 
    { 
        result.resize(1,1);
        result(0,0) = m_knots.greville(i);
    }

    // Look at gsBasis class for a description
    void connectivity(const gsMatrix<T> & nodes, 
                      gsMesh<T> & mesh) const;

    // Look at gsBasis class for a description
    void active_into(const gsMatrix<T> & u, gsMatrix<unsigned>& result) const;

    // Look at gsBasis class for a description
    gsMatrix<unsigned> * allBoundary( ) const ;

    // Look at gsBasis class for a description
    gsMatrix<unsigned> * boundaryOffset(boxSide const & s,unsigned offset) const;

    // Look at gsBasis class for a description
    gsConstantBasis<T> * boundaryBasis(boxSide const & s ) const;

    // Look at gsBasis class for a description
    gsMatrix<T> support() const ;

    // Look at gsBasis class for a description
    gsMatrix<T> support( const unsigned & i ) const ;

    /// Only meaningfull for periodic basis: For basis members that have
    /// a twin, this function returns the other twin index, otherwise it
    /// returns the same index as the argument
    unsigned twin(unsigned i) const ;

    // Look at gsBasis class for a description
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    virtual void evalSingle_into(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    virtual void evalFunc_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void derivSingle_into(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void deriv_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void deriv2Single_into(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void deriv2_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    gsMatrix<T> * laplacian(const gsMatrix<T> & u ) const ;

    // Look at gsBasis class for a description
    gsBSplineBasis * clone() const;

    // Look at gsBasis class for a description
    gsBasis<T> * tensorize(const gsBasis<T> & other) const;
    
    // Look at gsBasis class for a description
    virtual gsGeometry<T> * makeGeometry( const gsMatrix<T> & coefs ) const;

    // Look at gsBasis class for a description
    virtual gsGeometry<T> * makeGeometry( gsMovable< gsMatrix<T> > coefs ) const;

    /// Check the BSplineBasis for consistency
    bool check() const
    { 
        if ( m_periodic > 0 )
        {
            // Periodicity check wrt knot values
            return ( 
                m_knots.degree()      == m_p &&
                (int)m_knots.size()   >  2*m_p+1
                );           
        }
        else
        {
            return ( 
                m_knots.degree()      == m_p &&
                (int)m_knots.size()   >  2*m_p+1
                ); 
        }
    }
  
    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const;

    /// Return a string with detailed information on the basis.
    std::string detail() const;

    // Look at gsBasis class for a description
    virtual void evalDerSingle_into(unsigned i, const gsMatrix<T> & u, 
                                    int n, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    virtual void evalAllDers_into(const gsMatrix<T> & u, int n, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    virtual void evalAllDersSingle_into(unsigned i, const gsMatrix<T> & u, 
                                        int n, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    int degree(int i) const 
    { 
        GISMO_ASSERT(i==0,"Asked for degree(i) in 1D basis.");
        return m_p; 
    }

    int degree() const {return m_p;}

    // Look at gsBasis class for a description
    int maxDegree()   const { return m_p; }

    // Look at gsBasis class for a description
    int minDegree()   const { return m_p; }

    // Look at gsBasis class for a description
    int totalDegree() const { return m_p; }
 
    /// Returns the order of the B-spline  basis
    inline unsigned order() const { return m_p+1; }

    /// True iff the point \a pp is in the domain of the basis
    inline bool inDomain(T const & pp) const 
    { return ( (pp >= m_knots[m_p]) &&  (pp <= *(m_knots.end()-m_p-1) ) ); }

    /// Returns the starting value of the domain of the basis
    T domainStart() const { return m_knots[m_p]; }

    /// Returns the ending value of the domain of the basis
    T domainEnd() const { return m_knots[m_knots.size()-m_p-1]; }

    /// Returns length of the ``active" part of the knot vector.
    T _activeLength() const { return domainEnd() - domainStart(); }

    /// Returns the index of the first active (ie. non-zero) basis function at point u
    /// Takes into account non-clamped knots.
    inline unsigned firstActive(T u) const { 
        return ( inDomain(u) ? m_knots.findspan(u)-m_p : 0 );
    }

    // Number of active functions at any point of the domain
    inline unsigned numActive() const { return m_p + 1; }

    /// Returns the index of the first active (ie. non-zero) basis
    /// function at all columns (points) of u
    inline gsMatrix<unsigned,1> * firstActive(const gsMatrix<T,1> & u) const { 
        gsMatrix<unsigned,1> * fa = new gsMatrix<unsigned,1>(1, u.cols() );
        for ( index_t i = 0; i < u.cols(); i++ )
            fa->at(0,i) = ( inDomain(u(0,i)) ? m_knots.findspan(u(0,i))-m_p : 0 );
        return fa;
    }

    // Look at gsBasis class for a description
    gsDomain<T> * domain() const { return const_cast<KnotVectorType *>(&m_knots); }

    /// Returns the knot vector of the basis
    const KnotVectorType & knots () const { return m_knots;}
    KnotVectorType & knots ()       { return m_knots;}

    T knot(index_t const & i) const { return m_knots[i];}

    /// Inserts the knot \em knot in the underlying knot vector.
    void insertKnot(T knot, int mult=1)
    { m_knots.insert( knot, mult); }

    // Look at gsBasis class for a description
    void refineElements(std::vector<unsigned> const & elements)
    { m_knots.refineSpans(elements); }

    // Look at gsBasis class for a description
    void uniformRefine(int numKnots = 1, int mul=1)
    { m_knots.uniformRefine(numKnots,mul); }

    // Look at gsBasis class for a description
    void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots = 1, int mul=1);

    // Look at gsBasis class for a description
    void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots = 1, int mul=1);

    /// Refine the basis by inserting the given knots and perform knot
    /// refinement for the given coefficient matrix.
    void refine_withCoefs(gsMatrix<T>& coefs, const std::vector<T>& knots);

    /// Refine the basis by inserting the given knots and produce a sparse matrix which maps coarse coefficient vectors to refined ones.
    void refine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, const std::vector<T>& knots);

    /// \brief Increases the degree without adjusting the smoothness at inner
    /// knots, except from the knot values in \a knots (constrained
    /// knots of initial geometry)
    /// 
    /// This type of refinement is known as k-refinement.  Note that
    /// this type of refinement is ment to be performed after one (or
    /// more) with one h-refinement steps the parent mesh (\a other),
    /// otherwise this function is equivalent to p-refinement of \a other.
    /// Note: not tested yet
    void refine_k(const Self_t & other, int const & i = 1)
    { 
        GISMO_ASSERT( m_p >= other.m_p, "Degree of other knot-vector should be lower.");
        //for (typename std::vector<T>::iterator it = 
        //         m_knots.begin(); it != m_knots.end(); ++it)
        //    GISMO_ASSERT( has(*it), "Knot "<< *it<<" is not in the knot vector.");

        // grab unique knots
        std::vector<T> knots = other.m_knots.unique();
        // Increase the degree without adjusting any knot
        m_p += i; 
        m_knots.set_degree(m_p);
        // Adjust (reduce) smoothness to satisfy initial constraint knots
        m_knots.insert(knots,i);
    }

    /// \brief p-refinement (essentially degree elevation)
    void refine_p(int const & i = 1)
    { degreeElevate(i); }

    /// \brief Uniform h-refinement (placing \a i new knots inside each knot-span
    void refine_h(int const & i = 1)
    { uniformRefine(i); }
  
    /// Elevate the degree of the basis and preserve the smoothness
    void degreeElevate(int const & i = 1, int const dir = -1)
    { 
        GISMO_ASSERT( dir == -1 || dir == 0, "Invalid direction");
        m_p+=i; m_knots.degreeElevate(i); 
    }

    /// Elevate the degree of the basis and preserve the knot multiplicity
    void degreeIncrese(int const & i = 1, int const dir = -1)
    {
        GISMO_ASSERT( dir == -1 || dir == 0, "Invalid direction");
        m_p+=i; m_knots.degreeIncrease(i);
    }

    void setDegree(int const & i);

    /// Reduce the degree of the basis and preserve the smoothness
    void degreeReduce (int const & i = 1) 
    { 
        GISMO_ASSERT( i<=m_p, "Cannot reduce degree to negative");
        m_p-=i; m_knots.degreeReduce(i);
        //m_periodic =
    }

    /// Reduces spline continuity at interior knots by \a i
    void reduceContinuity(int const & i = 1) 
    { 
        GISMO_ASSERT( i>=0 && ( m_knots.size()>2*(m_p+1) || i<=m_p ), 
                      "Cannot achieve continuity less than C^{-1} at interior knots.");
        // TODO check: max interior mult + i <= m_p+1
        m_knots.increaseMultiplicity(i);
    }

    /// Tells, whether the basis is periodic.
    bool isPeriodic() const
    {
        return (m_periodic > 0);
    }

    // Look at gsBasis class for a description
    unsigned functionAtCorner(boxCorner const & c) const;

    /// Returns number of functions crossing the boundary of the knot vector.
    int numCrossingFunctions () const
    {
        return m_periodic;
    }

    /// Checks, if both endknots have multiplicity m_p + 1.
    bool isClamped() const
    {
        if( m_knots[0] != m_knots[m_p])
            return false;

        else if( m_knots[m_knots.size() - m_p -1] != m_knots[m_knots.size()-1])
            return false;

        else
            return true;
    }

    /// If flag is true, tries to convert the basis to periodic (succeeds only if the knot vector is suitable).
    void setPeriodic(bool flag = true)
    {
        if ( flag )
            _convertToPeriodic();
        else
            m_periodic = 0;
    }

    /// Returns the multiplicity of the first ``significant" knot
    /// (i.e., the m_p+1st). If it is different from the multiplicity
    /// of the corresponding knot at the end, returns zero.
    int borderKnotMult() const;

    typename gsBasis<T>::domainIter makeDomainIterator() const
    {
        return typename gsBasis<T>::domainIter(new gsTensorDomainIterator<T,1>(*this));
    }

    typename gsBasis<T>::domainIter makeDomainIterator(const boxSide & s) const
    {
        return ( s == boundary::none ? 
                 typename gsBasis<T>::domainIter(new gsTensorDomainIterator<T,1>(*this)) :
                 typename gsBasis<T>::domainIter(new gsTensorDomainBoundaryIterator<T,1>(*this, s))
                );
    }

    /// Moves the knot vectors to enforce periodicity.
    void enforceOuterKnotsPeriodic();

    // Look at gsBasis class for a description
    void reverse() { m_knots.reverse(); }


    /// Returns the size of the basis ignoring the bureaucratic way of
    /// turning the basis into periodic.
    int trueSize() const
    { return size() + m_periodic; }

private:

    /// Tries to convert the basis into periodic
    void _convertToPeriodic();

    /// Adjusts endknots so that the knot vector can be made periodic.
    void _stretchEndKnots();

public:

    /// \brief Helper function for evaluation with periodic basis.
    ///
    /// \param coefs coefficients (control points, one per row) before
    /// converting the basis into periodic.  
    ///
    ///\return copy of coefs with the first m_periodic rows copied to
    /// the last m_periodic rows.
    gsMatrix<T> perCoefs( const gsMatrix<T>& coefs ) const
    {
        gsMatrix<T> per_coefs = coefs;
        per_coefs.bottomRows( m_periodic ) = coefs.topRows( m_periodic );
        return per_coefs;
    }

    /// \brief Helper function for transforming periodic coefficients
    /// to full coefficients
    void expandCoefs(gsMatrix<T> & coefs) const
    {
        const index_t sz = coefs.rows();
        coefs.conservativeResize(sz+m_periodic, Eigen::NoChange);
        coefs.bottomRows( m_periodic ) = coefs.topRows( m_periodic );
    }

    /// \brief Helper function for transforming full coefficients to
    /// periodic coefficients
    void trimCoefs(gsMatrix<T> & coefs) const
    {
        const index_t sz = coefs.rows();
        coefs.conservativeResize(sz-m_periodic, Eigen::NoChange);
    }
 
// Data members
private:

    /// Degree
    int m_p;

    /// Knot vector
    KnotVectorType m_knots;
    
    /// Denotes whether the basis is periodic, ( 0 -- non-periodic, >0 -- number of ``crossing" functions)
    int m_periodic;

    /*/// Multiplicity of the p+1st knot from the beginning and from the end.
      int m_bordKnotMulti;*/

}; // class gsBSplineBasis


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBSplineBasis.hpp)
#endif
