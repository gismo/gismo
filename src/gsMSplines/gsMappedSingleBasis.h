/** @file gsMappedSingleBasis.h

    @brief Implementation of a piece of the gsMappedBasis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, P. Weinmueller
*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsUtils/gsMesh/gsMesh.h>
#include <gsMSplines/gsMappedBasis.h>

namespace gismo
{

template<short_t d,class T> class gsMappedBasis;

/**
   @brief Class gsMappedSingleBasis represents an indivisual .....of a

   Note that it does not own the underlying basis: if you delete
   the object the basis is not deleted. The lifetime of the
   underlying basis should be at least the lifetime of gsMappedSingleBasis.

*/
template<short_t d,class T>
class gsMappedSingleBasis : public gsBasis<T>
{
private:
    typedef memory::unique_ptr< gsDomainIterator<T> > domainIter;
    typedef gsBasis<T> Base;

public:
    /// Shared pointer for gsMappedSingleBasis
    typedef memory::shared_ptr< gsMappedSingleBasis > Ptr;

    /// Unique pointer for gsMappedSingleBasis
    typedef memory::unique_ptr< gsMappedSingleBasis > uPtr;

private:
    /// Default empty constructor
    gsMappedSingleBasis() : m_basis(nullptr) { }

public:

    /// Dimension of the parameter domain
    static const int Dim = d;

    /// Construct a basis function by a pointer to a basis and an index i
    gsMappedSingleBasis(gsMappedBasis<d,T> * basis, unsigned const & i = 0)
    : m_basis(basis),m_index(i)
    {
        GISMO_ASSERT( i<unsigned(m_basis->nPatches()),"Invalid basis function index" );
    }

    gsMappedSingleBasis( const gsMappedSingleBasis& other ) : Base( other )
    {
        m_basis = other.m_basis;
        m_index = other.m_index;
    }

    static uPtr make(   const gsMappedSingleBasis& other)
    { return uPtr( new gsMappedSingleBasis( other ) ); }

    ~gsMappedSingleBasis() { } //destructor

public:

    short_t domainDim() const
    {
        return d;
    }

    void connectivity(const gsMatrix<T> & nodes,
                      gsMesh<T>   & mesh) const
    {
        GISMO_UNUSED(nodes); GISMO_UNUSED(mesh);
        GISMO_NO_IMPLEMENTATION;
        //m_basis->connectivity(nodes,mesh);
    }

    // Look at gsBasis class for a description
    size_t numElements() const { return m_basis->getBase(m_index).numElements(); }

    /*
      void refine(gsMatrix<T> const & boxes)
      {  m_basis->refine(m_index,boxes);  }

      void refineElements(std::vector<unsigned> const & boxes)
      {  m_basis->refineElements(m_index,boxes);  }
    */

    /// Returns the indices of active (non zero) basis functions at points (columns of) u, as a list of indices, in result
    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const
    {
        m_basis->active_into(m_index,u,result);
    }

    /// Returns the number of active (nonzero) basis functions at points \a u in \a result.
    void numActive_into(const gsMatrix<T> & u, gsVector<unsigned>& result) const
    {
        // Assuming all patches have the same degree
        m_basis->numActive_into(m_index,u,result);
    }

    /// Returns a bounding box for the basis' domain
    gsMatrix<T> support() const
    {
        return m_basis->getBase(m_index).support();
    }

    /// Returns a bounding box for the basis' domain on the domain of *this
    gsMatrix<T> support(const index_t & i) const
    {
        typename gsMappedBasis<d,T>::IndexContainer sourceIndices;
        m_basis->getMapper().targetToSource(i,sourceIndices);
        // Get the support on the whole patch
        gsMatrix<T> supp = gsMatrix<T>::Zero(d,2);
        gsMatrix<T> localSupp;
        for (typename gsMappedBasis<d,T>::IndexContainer::iterator i = sourceIndices.begin(); i!=sourceIndices.end(); i++)
        {
            // Only consider local basis functions on the same patch
            if (m_basis->getPatch(*i)!=m_index) break;
            // Get the support of the basis function
            localSupp = m_basis->getBase(m_index).support(m_basis->getPatchIndex(*i));
            // If no support is available, we assign it
            if (supp.rows()==0 && supp.cols()==0)
            {
                supp = localSupp;
                continue;
            }
            // If a support is available, we increase it if needd
            for (index_t dim=0; dim!=d; dim++)
            {
                if (localSupp(dim,0) < supp(dim,0))
                    supp(dim,0) = localSupp(dim,0);
                if (localSupp(dim,1) > supp(dim,1))
                    supp(dim,1) = localSupp(dim,1);
            }
        }
        return supp;
        // return m_basis->getBase(m_index).support();
    }
    /// Returns the boundary basis on side s
    gsBasis<T>* boundaryBasis_impl(boxSide const & s) const
    {
        return m_basis->getBase(m_index).boundaryBasis(s).release(); // Wrong, Should return 1-D mappedSingleBasis
    }

    /// Evaluates the non-zero basis functions at value u.
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
    {
        // m_basis->evalGlobal_into(m_index,u,result);
        m_basis->eval_into(m_index,u,result);
    }

    /// Evaluates i-th basis functions at value u.
    void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const
    {
        m_basis->evalSingle_into(m_index,i,u,result);
    }

    /// Evaluates the (partial) derivatives of non-zero basis functions at (the columns of) u.
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
    {
        m_basis->deriv_into(m_index,u,result);
    }

    /// Evaluates the (partial)derivatives of the i-th basis function at (the columns of) u.
    void derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const
    {
        //GISMO_UNUSED(i); GISMO_UNUSED(u); GISMO_UNUSED(result);
        //GISMO_NO_IMPLEMENTATION;
        m_basis->derivSingle_into(m_index,i,u,result);
    }

    /// Evaluates the (partial) derivatives of the nonzero basis functions at points \a u into \a result.
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
    {
        m_basis->deriv2_into(m_index,u,result);
    }

    /// @brief Evaluate the (partial) derivatives of the \a i-th basis function
    /// at points \a u into \a result.
    void deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const
    {
        GISMO_UNUSED(i); GISMO_UNUSED(u); GISMO_UNUSED(result);
        GISMO_NO_IMPLEMENTATION;
    }

    /// @brief Evaluate the nonzero basis functions and their derivatives up
    /// to order \a n at points \a u into \a result.
    void evalAllDers_into(const gsMatrix<T> & u, int n, std::vector<gsMatrix<T> >& result) const
    {
        m_basis->evalAllDers_into(m_index,u,n,result);
    }

    /// @brief Evaluate the basis function \a i and its derivatives up
    /// to order \a n at points \a u into \a result.
    void evalAllDersSingle_into(index_t i, const gsMatrix<T> & u, int n, gsMatrix<T>& result) const
    {
        GISMO_UNUSED(i); GISMO_UNUSED(u); GISMO_UNUSED(n); GISMO_UNUSED(result);
        GISMO_NO_IMPLEMENTATION;
    }

    /// @brief Evaluate the (partial) derivative(s) of order \a n the \a i-th basis function
    /// at points \a u into \a result.
    void evalDerSingle_into(index_t i, const gsMatrix<T> & u, int n, gsMatrix<T>& result) const
    {
        GISMO_UNUSED(i); GISMO_UNUSED(u); GISMO_UNUSED(n); GISMO_UNUSED(result);
        GISMO_NO_IMPLEMENTATION;
    }

    GISMO_CLONE_FUNCTION(gsMappedSingleBasis)

    memory::unique_ptr<gsGeometry<T> > makeGeometry( gsMatrix<T> coefs ) const
    {
        GISMO_UNUSED(coefs);
        GISMO_NO_IMPLEMENTATION;
    }

    std::ostream &print(std::ostream &os) const
    {
        GISMO_UNUSED(os);
        os << "Mapped basis function "<< m_index << " / "<< m_basis->size()-1 <<"\n";
        return os;
    }

    /// Prints the object as a string with extended details.
    std::string detail() const
    {
        // By default just uses print(..)
        std::ostringstream os;
        print(os);
        return os.str();
    }

    //////////////////////////////////////////////////
    // Virtual member that may be implemented or not by the derived class
    //////////////////////////////////////////////////

    index_t size() const
    {
        return m_basis->size(m_index);
    }

    /// Returns the polynomial degree.
    short_t maxDegree() const
    {
        // TODO Not always working: make it more general
        return degree();
    }

    /// Returns the polynomial degree.
    short_t minDegree() const
    {
        // TODO Not always working: make it more general
        return degree();
    }

    /// Returns the polynomial degree.
    short_t degree() const
    {
        // TODO Not always working: make it more general
        return m_basis->maxDegree();                                   // must fix this (just took max_degree)
    }

    /// Returns the polynomial degree.
    short_t degree(short_t i) const
    {
        return m_basis->degree(m_index,i);
    }

    /// Check the LagrangeBasis for consistency
    bool check() const
    {
        GISMO_NO_IMPLEMENTATION;
    }

    /// The gsBasisFun points to the i-th basis function of m_basis
    /// after calling this setter.
    void setBasis( unsigned const & i ) const
    {
        GISMO_ASSERT( i<unsigned(m_basis->nPatches()),"Invalid basis index" );
        m_index = i;
    }


    /// Point to the first basis function of the basis
    void first() const
    {
        m_index = 0;
    }

    /// Points to the first basis function of the basis
    /// or return false if this is the last basis function
    bool next() const
    {
        return ( ++m_index < (index_t)m_basis->nPatches() );
    }

    /// Return the 1-d basis of the underlying tensor product basis for the \a i-th parameter component.
    const gsBasis<T>& component(short_t i) const
    {
        // TODO Not always working: make it more general
        return m_basis->getBase(m_index).component(i);
    }

    gsBasis<T>& component(short_t i)
    {
        // TODO Not always working: make it more general
        // return gsMappedSingleBasisComponent<d-1,T> (this, i);
        return m_basis->getBase(m_index).component(i);
    }

    typename gsBasis<T>::domainIter makeDomainIterator() const
    {
        // TODO Not always working: make it more general
        return m_basis->getBase(m_index).makeDomainIterator();
    }

    typename gsBasis<T>::domainIter makeDomainIterator(const boxSide & s) const
    {
        // TODO Not always working: make it more general
        return m_basis->getBase(m_index).makeDomainIterator(s);
    }


    gsMatrix<index_t> boundaryOffset(boxSide const & s, index_t offset) const
    {
        std::vector<index_t> temp, rtemp;
        m_basis->addLocalIndicesOfPatchSide(patchSide(m_index,s),offset,temp);
        m_basis->getMapper().sourceToTarget(temp,rtemp);

        // Better way for offset one: compute (anchors()) the normal derivatives at the boundary and return the indices
        if (offset == 1) // Small fix
        {
            GISMO_ASSERT(offset==1, "The indices of boundaryOffset(s,1) "
                                    "will be substract from boundaryOffset(s,0)");

            std::vector<index_t> diff, temp2, rtemp2;

            m_basis->addLocalIndicesOfPatchSide(patchSide(m_index,s),0,temp2);
            m_basis->getMapper().sourceToTarget(temp2,rtemp2);
            // Subtract the indizes of Offset = 0
            std::set_difference(rtemp.begin(), rtemp.end(), rtemp2.begin(), rtemp2.end(),
                        std::inserter(diff, diff.begin()));
            rtemp = diff;
        }

        return makeMatrix<index_t>(rtemp.begin(),rtemp.size(),1 );
    }

    index_t functionAtCorner(boxCorner const & c) const
    {
        index_t cindex = m_basis->getBase(m_index).functionAtCorner(c);
        cindex = m_basis->_getLocalIndex(m_index,cindex);
        GISMO_ENSURE(m_basis->getMapper().sourceIsId(cindex),"Corner function has no identity map, i.e. there are more than 1 functions associated to the corner?");
        std::vector<index_t> indices;
        m_basis->getMapper().sourceToTarget(cindex,indices);
        GISMO_ASSERT(indices.size()==1,"Size of the indices returned for the corner basis function should be 1 but is "<<indices.size()<<". Otherwise, there are more than 1 functions associated to the corner");
        return indices.front();
    }

    
// Data members
private:
    gsMappedBasis<d,T> * m_basis;
    mutable index_t m_index;

}; // class gsMappedSingleBasis


} // namespace gismo
