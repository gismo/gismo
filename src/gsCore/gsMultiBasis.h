/** @file gsMultiBasis.h

    @brief Provides declaration of MultiBasis class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsGeometry.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsDofMapper.h>
#include <gsCore/gsBoxTopology.h>
#include <gsPde/gsBoundaryConditions.h>

namespace gismo
{

/** @brief
    Holds a set of patch-wise bases and their
    topology information.

    \tparam T coefficient type
*/
template<class T>
class gsMultiBasis
{

public:
    typedef memory::shared_ptr<gsMultiBasis> Ptr;

    typedef std::vector<gsBasis<T> *> BasisContainer;

public:

    /// Type definitions
    typedef typename BasisContainer::size_type size_t;
    typedef typename BasisContainer::iterator iterator;
    typedef typename BasisContainer::const_iterator const_iterator;

    typedef          gsBasis<T> & reference;
    typedef          const gsBasis<T> & const_reference;

public:

    /// Default empty constructor
    gsMultiBasis() { }

    /// Create a multi-basis instance grom a gsMultiPatch
    explicit gsMultiBasis( const gsMultiPatch<T> & mpatch );

    /// Create from a vector of bases and topology
    gsMultiBasis( const BasisContainer& bases, const gsBoxTopology & topology)
    : m_topology( topology )
    {
       m_bases.resize(bases.size());
       cloneAll( bases.begin(), bases.end(), m_bases.begin() );
    }
    
    /// Create a single-basis instance
    explicit gsMultiBasis( const gsBasis<T> & geo );


    /// Create from bases and boundary/interface information
    gsMultiBasis( const BasisContainer& bases,
                  const std::vector<patchSide>& boundary,
                  const std::vector<boundaryInterface>& interfaces )
    : m_topology( bases[0]->dim(), bases.size(), boundary, interfaces )      
    { 
       m_bases.resize(bases.size());
       cloneAll( bases.begin(), bases.end(), m_bases.begin() );
    }
    
    /// Destructor
    ~gsMultiBasis();

    /// Copy constructor (makes deep copy)
    gsMultiBasis( const gsMultiBasis& other );

    /// Assignment operator (uses copy-and-swap idiom)
    gsMultiBasis& operator= ( gsMultiBasis other )
    {
        this->swap( other );
        return *this;
    }

    /// Clone function. Used to make a copy of the object
    gsMultiBasis* clone() const {
        return new gsMultiBasis( *this );
    }

public:

    /// Get a const-iterator to the patches
    /// \return an iterator to the beginning of the patches
    const_iterator begin() const 
    {
        return m_bases.begin();
    }

    /// Get a const iterator to the end of the patches
    /// \return an iterator to the end of the patches
    const_iterator end() const 
    {
        return m_bases.end();
    }

    /// Get an iterator to the beginning of the  patches
    /// \return an iterator to the beginning of the  patches
    iterator begin() {
        return m_bases.begin();
    }

    /// Get an iterator to the end of the  patches
    /// \return an iterator to the end of the  patches
    iterator end() {
        return m_bases.end();
    }

    const gsBoxTopology & topology() const { return m_topology; }

public:

    /// Assess i-th parametric basis
    const_reference at(size_t i) const
	{return *m_bases.at(i);}

    //reference at(size_t i)
    //{return *m_bases.at(i);}

	const_reference operator[](size_t i) const
    {return *m_bases[i];}
    
	reference operator[](size_t i)			
    {return *m_bases[i];}

    // reference front()
    // {return *m_bases.front();}

	const_reference front() const
    {return *m_bases.front();}
    
    // reference back()						
    // {return *m_bases.back();}

    const_reference back() const
    {return *m_bases.back();}

    void pop_back()
    {m_bases.pop_back();}

public:

    /// Swap with another gsMultiBasis.
    void swap(gsMultiBasis& other)
    {
        m_topology.swap( m_topology );

        m_bases.swap( other.m_bases );
    }

    /// Prints the object as a string.
    std::ostream& print( std::ostream& os ) const;

    /// Dimension of the parameter domain (must match for all bases).
    int dim() const
    { return m_bases[0]->dim();}

    /// @brief Returns the polynomial degree of basis \a i in component \a j,
    /// if the basis is of polynomial or piecewise polynomial type.
    int degree(int i = 0, int comp = 0) const
    { return m_bases[i]->degree(comp);}

    /// Maximum degree with respect to variable \a k.
    int maxDegree(int k) const;

    /// Minimum degree with respect to variable \a k.
    int minDegree(int k) const;

    /// The number of basis functions in basis \a i.
    int size(int i) const
    { return m_bases[i]->size();}

    /// The total number of basis functions in all bases
    size_t totalSize() const
    {
        size_t sum = 0;
        for (size_t k = 0; k < m_bases.size(); ++k)
            sum += m_bases[k]->size();
        return sum;
    }

    /// The total number of elements in all patches
    size_t totalElements() const
    {
        size_t sum = 0;
        for (size_t k = 0; k < m_bases.size(); ++k)
            sum += m_bases[k]->numElements();
        return sum;
    }

    /// Number of bases
    size_t nBases() const          { return m_bases.size(); }

    /// Returns a vector of paches // to do : replace by copies
    BasisContainer const& bases() const { return m_bases; }

    /// Makes a deep copy of all bases and puts them in a vector
    std::vector<gsBasis<T> *> basesCopy() const;

    /// Return the \a i-th basis block.
    const gsBasis<T> & basis( std::size_t i ) const
    {
        GISMO_ASSERT( i < m_bases.size(), "Invalid patch index requested from gsMultiBasis" );
        return *m_bases[i];
    }

    /// Add a basis
    void addBasis( gsBasis<T>* g );

    /// Search for the given basis and return its index.
    int findBasisIndex( gsBasis<T>* g ) const;
    
    /// @brief Add an interface joint betweeen side \a s1 of geometry
    /// \a g1 side \a s2 of geometry \a g2.
    ///
    /// \todo add orientation information
    void addInterface( gsBasis<T>* g1, boundary::side s1,
                       gsBasis<T>* g2, boundary::side s2 );

    /// Add side s of patch g to the outer boundary of the domain
    void addPatchBoundary( gsBasis<T>* g, boundary::side s ) 
    {
        const int p =findBasisIndex( g );
        m_topology.addBoundary( patchSide( p, s ) );
    }
    
    /// Refine every basis uniformly by inserting \a numKnots new knots on each knot span
    void uniformRefine(int numKnots = 1)
    {
        for (size_t k = 0; k < m_bases.size(); ++k)
        {
            m_bases[k]->uniformRefine(numKnots);
        }
    }

    //add to the domain structure the boxes defined in "boxes"
    void refine(int k, gsMatrix<T> const & boxes)
    {
        m_bases[k]->refine(boxes);
    }


    /// Elevate the degree of every basis by the given amount.
    void degreeElevate(int const& i = 1)
    {
        for (size_t k = 0; k < m_bases.size(); ++k)
        {
            m_bases[k]->degreeElevate(i);
        }
    }

    /// Elevate the degree of every basis by the given amount.
    void degreeElevateComponent(unsigned dir, int const& i = 1)
    {
        for (size_t k = 0; k < m_bases.size(); ++k)
        {
            m_bases[k]->degreeElevateComponent(dir,i);
        }
    }

    /// Reduce the degree of the basis by the given amount.
    void degreeReduce(int const& i = 1)
    {
        for (size_t k = 0; k < m_bases.size(); ++k)
        {
            m_bases[k]->degreeReduce(i);
        }
    }

    /// Set the degree of the basis.
    void setDegree(int const& i)
    {
        for (size_t k = 0; k < m_bases.size(); ++k)
        {
            m_bases[k]->setDegree(i);
        }
    }

   BasisContainer & patchBases() 
   {
       return m_bases;
   }

    void setTopology(const gsBoxTopology & tpl)
	{
        m_topology = tpl;
    }


    void getMapper(bool conforming, 
                   const gsBoundaryConditions<T> & bc, 
                   gsDofMapper & mapper) const;


    void getMapper(bool conforming, gsDofMapper & mapper) const;

    // to remove
    gsDofMapper * makeMapper(bool conforming) const
    {
        gsDofMapper * mapper = new gsDofMapper;
        getMapper(conforming, *mapper);
        return mapper;
    }

    // to remove
    gsDofMapper * makeMapper(bool conforming, 
                             const gsBoundaryConditions<T> & bc) const
    {
        gsDofMapper * mapper = new gsDofMapper;
        getMapper(conforming, bc, *mapper);
        return mapper;
    }

    gsDofMapper * makeIdMapper()  const
	{
        gsDofMapper * mapper = new gsDofMapper();

        int nDofs = 0;
        for (size_t k = 0; k < m_bases.size(); ++k)
            nDofs += m_bases[k]->size();

        mapper->setIdentity(m_bases.size(), nDofs );       
        return mapper;
    }


    

    // Data members
private:

    BasisContainer m_bases;

    gsBoxTopology m_topology;

}; // class gsMultiBasis


//////////////////////////////////////////////////
//////////////////////////////////////////////////

/// Print (as string) a multibasis structure
template<class T>
std::ostream& operator<<( std::ostream& os, const gsMultiBasis<T>& b )
{
    return b.print( os );
}


} // namespace gismo


#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsMultiBasis.hpp)
#endif
