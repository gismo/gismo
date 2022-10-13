/** @file gsFupBasis.h

    @brief Provides declaration of FupBasis class

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D. Mokris
*/


#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsConstantBasis.h>

#include <gsTensor/gsTensorBasis.h>
#include <gsTensor/gsTensorDomainIterator.h>
#include <gsTensor/gsTensorDomainBoundaryIterator.h>

#include <gsNurbs/gsKnotVector.h>

namespace gismo
{

/// @brief Traits for FupBasis in more dimensions
template<short_t d, class T>
struct gsFupTraits
{
    typedef gsKnotVector<T> KnotVectorType;

    //typedef gsTensorFupBasis<d,T> Basis;
    //typedef gsTensorNurbsBasis<d,T>   RatBasis;
    //typedef gsTensorFup<d,T>      Geometry;
    //typedef gsTensorNurbs<d,T>        RatGeometry;
};
template<class T>
struct gsFupTraits<1,T>
{
    typedef gsKnotVector<T> KnotVectorType;
    typedef gsFupBasis<T>         Basis;
    typedef gsNurbsBasis<T>           RatBasis;
    typedef gsFup<T>              Geometry;
    typedef gsNurbs<T>                RatGeometry;
};
template<class T>
struct gsFupTraits<0,T>
{
    typedef gsKnotVector<T> KnotVectorType;
    typedef gsConstantBasis<T>                       Basis;
    typedef gsConstantBasis<T>                       RatBasis;
    typedef gsConstantFunction<T>                    Geometry;
    typedef gsConstantFunction<T>                    RatGeometry;
};

/** \brief
    A univariate B-spline basis.

    \tparam T coefficient type
    \tparam KnotVectorType the type of knot vector to use

    \ingroup basis
    \ingroup Nurbs
*/
template<class T>
class gsFupBasis<T> : public gsBasis<T>
{
    /// @brief Degree
    short_t m_p;

    /// @brief Knot vector
    KnotVectorType m_knots;
    
    /// @brief Denotes whether the basis is periodic, ( 0 -- non-periodic, >0 -- number of ``crossing" functions)
    int m_periodic;

    /*/// @brief Multiplicity of the p+1st knot from the beginning and from the end.
      int m_bordKnotMulti;*/

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
    private: virtual gsFupBasis * clone_impl() const = 0;
    public: uPtr clone() const { return uPtr(dynamic_cast<Self_t*>(clone_impl())); }

public:

    // Look at gsBasis class for a description
    short_t domainDim() const { return Dim; }

    // Look at gsBasis class for a description
    index_t size() const { return 10; }

    // Look at gsBasis class for a description
    size_t numElements() const { return m_knots.numElements(); }

    /// @brief Returns the anchors (greville points) of the basis
    void anchors_into(gsMatrix<T> & result) const 
    { 
        //CV points
    }

    /// @brief Returns the anchors (greville points) of the basis
    void anchor_into(index_t i, gsMatrix<T> & result) const
    { 
        result.resize(1,1);
        result(0,0) = m_knots.greville(i);
    }

    // Look at gsBasis class for a description
    void connectivity(const gsMatrix<T> & nodes, 
                      gsMesh<T> & mesh) const;

    // Look at gsBasis class for a description
    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const;

    // Look at gsBasis class for a description
    gsMatrix<T> support() const ;

    // Look at gsBasis class for a description
    gsMatrix<T> support(const index_t & i ) const ;

    // Look at gsBasis class for a description
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    virtual void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void deriv_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void deriv2_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    gsMatrix<T> laplacian(const gsMatrix<T> & u ) const ;
    
    /// @brief Check the FupBasis for consistency
    bool check() const
    { 
        return true;
    }
  
    /// @brief Prints the object as a string.
    std::ostream &print(std::ostream &os) const;

    // Look at gsBasis class for a description
    virtual void evalDerSingle_into(index_t i, const gsMatrix<T> & u,
                                    int n, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    virtual void evalAllDers_into(const gsMatrix<T> & u, int n,
                                  std::vector<gsMatrix<T> >& result) const;

    // Look at gsBasis class for a description
    virtual void evalAllDersSingle_into(index_t i, const gsMatrix<T> & u,
                                        int n, gsMatrix<T>& result) const;

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
    { return ( (pp >= *(m_knots.begin()+m_p)) &&  (pp <= *(m_knots.end()-m_p-1) ) ); }

    /// @brief Returns the starting value of the domain of the basis
    T domainStart() const { return *(m_knots.begin()+m_p); }

    /// @brief Returns the ending value of the domain of the basis
    T domainEnd() const { return *(m_knots.end()-m_p-1); }

    // Number of active functions at any point of the domain
    inline index_t numActive() const { return m_p + 1; }

    // Look at gsBasis class for a description
    void uniformRefine(int numKnots = 1, int mul=1)
    { m_knots.uniformRefine(numKnots,mul); }

    /// @brief Elevate the degree of the basis and preserve the smoothness
    void degreeElevate(short_t const & i = 1, short_t const dir = -1)
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT( dir == -1 || dir == 0, "Invalid direction");
        m_p+=i; m_knots.degreeElevate(i); 
    }

    // Look at gsBasis for documentation
    void degreeReduce (short_t const & i = 1, short_t const dir = -1)
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT( dir == -1 || dir == 0, "Invalid direction");
        GISMO_ASSERT( i<=m_p, "Cannot reduce degree to negative");
        m_p-=i; m_knots.degreeReduce(i);
        //m_periodic =
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
