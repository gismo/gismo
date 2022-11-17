/** @file gsFupBasis.h

    @brief Provides declaration of FupBasis class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

*/


#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsConstantBasis.h>

#include <gsTensor/gsTensorBasis.h>
#include <gsTensor/gsTensorDomainIterator.h>
#include <gsTensor/gsTensorDomainBoundaryIterator.h>

#include <gsNurbs/gsKnotVector.h>

// Start  Fortran declarations
extern "C" double __fup_0_16_d_MOD_fupn(int * deg, double * anchor, double * u,
                                        double * ch_knot_length, int * deriv_order);

extern "C" void __fup_0_16_d_MOD_racun();
// End  Fortran declarations

namespace gismo
{

/// @brief Traits for FupBasis in more dimensions
template<short_t d, class T>
struct gsFupTraits
{
    typedef gsKnotVector<T> KnotVectorType;

    //typedef gsTensorFupBasis<d,T> Basis;
    //typedef gsTensorRatBasis<d,T> RatBasis;
    //typedef gsTensorFup<d,T>      Geometry;
    //typedef gsTensorNurbs<d,T>    RatGeometry;
};

/** \brief
    A univariate B-spline basis.

    \tparam T coefficient type
    \tparam KnotVectorType the type of knot vector to use

    \ingroup basis
    \ingroup Nurbs
*/
template<class T>
class gsFupBasis : public gsBasis<T>
{
    mutable index_t m_p;
    gsKnotVector<T> m_knots;
    mutable T m_knot_length;

public:
    typedef gsBasis<T> Base;

    typedef gsFupBasis<T> Self_t;

    /// @brief Coefficient type
    typedef T Scalar_t;

    /// @brief Associated geometry type
    //typedef typename gsFupTraits<T>::Geometry GeometryType;

    /// @brief Associated Boundary basis type
    //typedef typename gsFupTraits<T>::Basis BoundaryBasisType;

    /// @brief Dimension of the parameter domain
    static const short_t Dim = 1;

    /// @brief Smart pointer for gsFupBasis
    typedef memory::shared_ptr< Self_t > Ptr;

    /// @brief Smart pointer for gsFupBasis
    typedef memory::unique_ptr< Self_t > uPtr;

public:
    GISMO_CLONE_FUNCTION(gsFupBasis)

    gsFupBasis()
    {
        m_p = 0;
        m_knots.initClamped(0);
    }

    gsFupBasis(const T u0, const T u1,
               const unsigned interior,
               const int degree)
    {
        m_p = degree;
        //m_knots.initUniform(u0, u1, interior, m_p+2, 1, m_p+1);
        m_knots.initUniform(u0, u1, interior, 1, 1, m_p+1);
        m_knot_length = (T)(u1-u0)/(m_knots.numElements());

        //Precompute values at rational points
        __fup_0_16_d_MOD_racun();
    }

private:

    void fupn_eval(const gsMatrix<T> & u, int deriv_order,
                   gsMatrix<T>& result) const
    {
        result.resize(m_p+2, u.cols() );
        index_t act;
        for ( index_t k = 0; k!=u.cols(); ++k)
        {
            act = firstActive(u(0,k));
            for ( index_t i = 0; i!=m_p+2; ++i)
                result(i,k) = fupn_eval_single(u(0,k), deriv_order, act++);
        }
    }

    T fupn_eval_single(T u, int deriv_order, index_t i) const
    {
        T anc = m_knots.greville(i);
        return __fup_0_16_d_MOD_fupn(&m_p, &anc, &u,
                                     &m_knot_length, &deriv_order);
                //gsInfo<<std::fixed<<std::setw(6) << "fupn("<<m_p<<","<<anc<<","<<u(k)<<"," <<m_knot_length<<","<<deriv_order<<") = " << result(i,k) <<"\n";
    }

public:
         
    memory::unique_ptr<gsGeometry<T> > makeGeometry( gsMatrix<T> coefs ) const
    { return nullptr; }

    // Look at gsBasis class for a description
    short_t domainDim() const { return Dim; }

    // Look at gsBasis class for a description
    index_t size() const { return m_knots.size() - m_p - 2; }

    /// @brief Returns the anchors (greville points) of the basis
    void anchors_into(gsMatrix<T> & result) const
    {
        m_knots.greville_into(result);
    }

    /// @brief Returns the anchors (greville points) of the basis
    void anchor_into(index_t i, gsMatrix<T> & result) const
    {
        result.resize(1,1);
        result(0,0) = m_knots.greville(i);
    }

    // Look at gsBasis class for a description
    void connectivity(const gsMatrix<T> & nodes,
                      gsMesh<T> & mesh) const { }


    inline index_t firstActive(T u) const
    {
        return ( inDomain(u) ? (m_knots.iFind(u)-m_knots.begin()) - m_p - 1 : 0 );
    }

    // Look at gsBasis class for a description
    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const
    {
        result.resize(m_p+2, u.cols());
        for (index_t j = 0; j < u.cols(); ++j)
        {
            index_t first = firstActive(u(0,j));
            for (int i = 0; i != m_p+2; ++i)
                result(i,j) = first++;
        }
    }

    // Look at gsBasis class for a description
    gsMatrix<T> support() const
    {
        gsMatrix<T> res(1,2);
        res << domainStart() , domainEnd() ;
        return res;
    }

    // Look at gsBasis class for a description
    gsMatrix<T> support(const index_t & i ) const
    {
        GISMO_ASSERT( static_cast<size_t>(i) < m_knots.size()-m_p-1,
                      "Invalid index of basis function." );
        gsMatrix<T> res(1,2);
        res << ( i > m_p ? m_knots[i] : m_knots[m_p] ),
            ( static_cast<size_t>(i) < (m_knots.size()-2*m_p-2) ? m_knots[i+m_p+1] :
              m_knots[m_knots.size()-m_p-1] );
        return res ;
    }

    // Look at gsBasis class for a description
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
    {
        fupn_eval(u,0,result);
    }

    // Look at gsBasis class for a description
    void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const
    {
        result.resize(1, u.cols() );
        for ( index_t k = 0; k!=u.cols(); ++k)
            result(0,k) = fupn_eval_single(u(0,k), 0, i);
    }

    // Look at gsBasis class for a description
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
    {
        fupn_eval(u,1,result);
    }

    // Look at gsBasis class for a description
    void derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const {
        result.resize(1, u.cols() );
        for ( index_t k = 0; k!=u.cols(); ++k)
            result(0,k) = fupn_eval_single(u(0,k), 1, i);
    }

    // Look at gsBasis class for a description
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const { }

    // Look at gsBasis class for a description
    void deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const { }

    // Look at gsBasis class for a description
    gsMatrix<T> laplacian(const gsMatrix<T> & u ) const { return gsMatrix<T>(); }

    /// @brief Check the FupBasis for consistency
    bool check() const
    {
        return true;
    }

    /// @brief Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "FupBasis: deg=" <<degree()
           << ", size="<< size() << ", knot vector:\n";
        os << m_knots;
        os << "\nCharacteristic length: "<< m_knot_length <<"\n";
        return os;
    }

    // Look at gsBasis class for a description
    virtual void evalDerSingle_into(index_t i, const gsMatrix<T> & u,
                                    int n, gsMatrix<T>& result) const{ }

    // Look at gsBasis class for a description
    virtual void evalAllDers_into(const gsMatrix<T> & u, int n,
                                  std::vector<gsMatrix<T> >& result) const
    {
        result.resize(n+1);
        for ( index_t j = 0; j<=n; ++j)
            fupn_eval(u,j,result[j]);
    }

    // Look at gsBasis class for a description
    virtual void evalAllDersSingle_into(index_t i, const gsMatrix<T> & u,
                                        int n, gsMatrix<T>& result) const
    {
        result.resize(n+1, u.cols() );
        for ( index_t j = 0; j<=n; ++j)
            for ( index_t k = 0; k!=u.cols(); ++k)
                result(j,k) = fupn_eval_single(u(0,k), j, i);
    }

    // Look at gsBasis class for a description
    short_t degree(short_t i) const
    {
        return m_p;
    }

    short_t degree() const {return m_p;}

    // Look at gsBasis class for a description
    short_t maxDegree()   const { return m_p; }

    // Look at gsBasis class for a description
    short_t minDegree()   const { return m_p; }

    // Look at gsBasis class for a description
    short_t totalDegree() const { return m_p; }

    /// @brief Returns the order of the B-spline  basis
    inline unsigned order() const { return m_p+1; }

    /// @brief True iff the point \a pp is in the domain of the basis
    inline bool inDomain(T const & pp) const
    { return true; }

    /// @brief Returns the starting value of the domain of the basis
    T domainStart() const { return 0; }

    /// @brief Returns the ending value of the domain of the basis
    T domainEnd() const { return 1; }

    // Number of active functions at any point of the domain
    inline index_t numActive() const { return m_p + 1; }

    // Look at gsBasis class for a description
    void uniformRefine(int numKnots = 1, int mul=1)
    {
        GISMO_ASSERT(1==mul, "multiple knot ?");
        m_knots.uniformRefine(numKnots,mul);
    }

    /// @brief Elevate the degree of the basis and preserve the smoothness
    void degreeElevate(short_t const & i = 1, short_t const dir = -1)
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT( dir == -1 || dir == 0, "Invalid direction");

    }

    // Look at gsBasis for documentation
    void degreeReduce (short_t const & i = 1, short_t const dir = -1)
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT( dir == -1 || dir == 0, "Invalid direction");
        GISMO_ASSERT( i<=m_p, "Cannot reduce degree to negative");

    }

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

}; // class gsFupBasis<1>


} // namespace gismo


//#ifndef GISMO_BUILD_LIB
//#include GISMO_HPP_HEADER(gsFupBasis.hpp)
//#endif
