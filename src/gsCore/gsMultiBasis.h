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
    
    \ingroup Core
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

    /// Create a multi-basis instance from a gsMultiPatch
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

    /// Clear (delete) all patches
    void clear()
    {
        m_topology.clearAll();
        m_bases   .clear();
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
    int degree(size_t i = 0, int comp = 0) const
    { 
        GISMO_ASSERT( i < m_bases.size(),
        "Invalid patch index "<<i<<" requested from gsMultiBasis" );
        return m_bases[i]->degree(comp);
    }

    /// @brief Maximum degree with respect to variable \a k.
    int maxDegree(int k) const;

    /// @brief Minimum degree with respect to variable \a k.
    int minDegree(int k) const;

    /// @brief Maximum degree with respect to all variables
    int maxCwiseDegree() const;

    /// @brief Minimum degree with respect to all variables
    int minCwiseDegree() const;

    /// @brief The number of basis functions in basis \a i.
    int size(size_t i) const
    {
        GISMO_ASSERT( i < m_bases.size(),
        "Invalid patch index "<<i<<" requested from gsMultiBasis" ); 
        return m_bases[i]->size();
    }

    /// @brief The total number of basis functions in all bases
    size_t totalSize() const
    {
        size_t sum = 0;
        for (size_t k = 0; k < m_bases.size(); ++k)
            sum += m_bases[k]->size();
        return sum;
    }

    /// @brief The total number of elements in all patches
    size_t totalElements() const
    {
        size_t sum = 0;
        for (size_t k = 0; k < m_bases.size(); ++k)
            sum += m_bases[k]->numElements();
        return sum;
    }

    /// @brief Number of patch-wise bases
    size_t nBases() const          { return m_bases.size(); }

    /// Return the \a i-th basis block.
    const gsBasis<T> & basis( std::size_t i ) const
    {
        GISMO_ASSERT( i < m_bases.size(), 
        "Invalid patch index"<<i<<" requested from gsMultiBasis" );
        return *m_bases[i];
    }

    /// @brief Add a basis (ownership of the pointer is also acquired)
    void addBasis( gsBasis<T>* g );

    /// @brief Search for the given basis and return its index.
    int findBasisIndex( gsBasis<T>* g ) const;
    
    /// @brief Add an interface joint between side \a s1 of geometry
    /// \a g1 side \a s2 of geometry \a g2.
    ///
    /// \todo add orientation information
    void addInterface( gsBasis<T>* g1, boxSide s1,
                       gsBasis<T>* g2, boxSide s2 );

    /// @brief Add side s of patch g to the outer boundary of the domain
    void addPatchBoundary( gsBasis<T>* g, boxSide s ) 
    {
        const int p =findBasisIndex( g );
        m_topology.addBoundary( patchSide( p, s ) );
    }
    
    /// @brief Refine every basis uniformly by inserting \a numKnots
    /// new knots on each knot span
    void uniformRefine(int numKnots = 1, int mul=1)
    {
        for (size_t k = 0; k < m_bases.size(); ++k)
        {
            m_bases[k]->uniformRefine(numKnots,mul);
        }
    }

    // @brief Refine the boxes defined by "boxes"
    void refine(int k, gsMatrix<T> const & boxes)
    {
        m_bases[k]->refine(boxes);
    }

    /// @brief Refine the are defined by \em boxes
    /// on patch \em k.
    ///
    /// See gsHTensorBasis::refineElements() for further documentation.
    void refineElements(int k, std::vector<unsigned> const & boxes)
    {
        m_bases[k]->refineElements(boxes);
    }

    /// @brief Refine the are defined by \em boxes
    /// on patch \em k with extension \em refExt.
    ///
    /// See gsHTensorBasis::refineWithExtension() for further documentation.
    void refine(size_t k, gsMatrix<T> const & boxes, int refExt)
    {
        GISMO_ASSERT( k < m_bases.size(),
                      "Invalid patch index "<<k<<" requested from gsMultiBasis" );
        m_bases[k]->refine( boxes, refExt);
    }

    /** @brief Checks if the interfaces \em bivec are fully matching, and if not, repairs them, i.e., makes them fully matching.
    *
    * \remarks Designed for gsHTensorBasis and derived bases.
    * Assumes that the meshes on all levels of the gsHTensorBasis
    * are fully matching.
    *
    * Calls repairInterface() for each boundaryInterface in \em bivec.
    */
    void repairInterfaces( const std::vector< boundaryInterface > & bivec )
    {
        size_t kmax = 2*bivec.size();
        size_t k = 0;
        bool somethingChanged = false;
        do
        {
            somethingChanged = false;
            for( size_t i = 0; i < bivec.size(); ++i )
                somethingChanged = ( somethingChanged || repairInterface( bivec[i] ) );
            k++; // just to be sure this can't go wrong
        }
        while( somethingChanged && k <= kmax );
    }

    /** @brief Checks if the interface is fully matching, and if not, repairs it.
    *
    * \remarks Designed for gsHTensorBasis and derived bases.
    * Assumes that the meshes on all levels of the gsHTensorBasis
    * are fully matching.
    *
    * \returns true, if something was repaired, i.e., if the mesh on the interface was changed.
    */
    bool repairInterface( const boundaryInterface & bi );

    /// @brief Elevate the degree of every basis by the given amount.
    void degreeElevate(int const& i = 1, int const dir = -1)
    {
        for (size_t k = 0; k < m_bases.size(); ++k)
            m_bases[k]->degreeElevate(i,dir);
    }

    /// Reduce the degree of the basis by the given amount.
    void degreeReduce(int const& i = 1)
    {
        for (size_t k = 0; k < m_bases.size(); ++k)
            m_bases[k]->degreeReduce(i);
    }

    /// Set the degree of the basis.
    void setDegree(int const& i)
    {
        for (size_t k = 0; k < m_bases.size(); ++k)
        {
            m_bases[k]->setDegree(i);
        }
    }

    /// Reduce the continuity by i
    void reduceContinuity(int const i)
    {
        for (size_t k = 0; k < m_bases.size(); ++k)
        {
            m_bases[k]->reduceContinuity(i);
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
                   int unk,
                   gsDofMapper & mapper, 
                   bool finalize = true) const;

    void getMapper(bool conforming, 
                   const gsBoundaryConditions<T> & bc,
                   gsDofMapper & mapper,
                   bool finalize = true) const
    { getMapper(conforming, bc, 0, mapper, finalize); }


    void getMapper(bool conforming, gsDofMapper & mapper, bool finalize = true) const;

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
        getMapper(conforming, bc, 0, *mapper);// using unknown 0
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


//private: // to do

    /**
     * @brief Matches the degrees of freedom along an interface.
     *
     * The boundaryInterface specifying the interface
     * is passed as argument.\n
     * The degrees of freedom (DOFs) along the interface from
     * both patches are matched in the sense that they are
     * then treated as one global DOF by the mapper.
     *
     * @todo Check if this description is accurate.
     *
     * @remarks This function assumes
     * tensor-product-structure with matching mesh along
     * the interface! If the gsMultiBasis contains
     * bases of class gsHTensorBasis (or derived), it
     * calls the function matchInterfaceHTensor().
     *
     * @param bi specifies the interface to be matched
     * @param mapper the gsDofMapper which should know that
     * these interface-DOFs are matched.
     */
    void matchInterface(const boundaryInterface & bi,
                        gsDofMapper & mapper) const;

    /**
     * @brief Matches the degrees of freedom along an interface.
     *
     * Same as matchInterface(), but for bases of
     * class gsHTensorBasis or derived classes.
     *
     * \remarks Assumes that the meshes on all levels of
     * the gsHTensorBasis are fully matching at the interface.
     *
     * @param bi specifies the interface to be matched
     * @param mapper the gsDofMapper which should know that
     * these interface-DOFs are matched.
     */
    void matchInterfaceHTensor(const boundaryInterface & bi,
                               gsDofMapper & mapper) const;

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


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMultiBasis.hpp)
#endif
