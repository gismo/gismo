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
#include <gsAssembler/gsAssemblerOptions.h>


namespace gismo
{

/** @brief
    Holds a set of patch-wise bases and their
    topology information.

    \tparam T coefficient type
    
    \ingroup Core
*/
template<class T>
class gsMultiBasis : public gsFunctionSet<T>
{
public:
    typedef memory::shared_ptr<gsMultiBasis> Ptr;
    typedef memory::unique_ptr<gsMultiBasis> uPtr;

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
    gsMultiBasis(BasisContainer& bases, const gsBoxTopology & topology)
        : m_topology( topology )
    {
        m_bases.swap(bases);// consumes the pointers
    }
    
    /// Create a single-basis instance
    gsMultiBasis( const gsBasis<T> & geo );

    /// Create from bases and boundary/interface information
    gsMultiBasis( BasisContainer & bases,
                  const std::vector<patchSide>& boundary,
                  const std::vector<boundaryInterface>& interfaces )
        : m_topology( bases[0]->dim(), bases.size(), boundary, interfaces )
    {
        m_bases.swap(bases);// consumes the pointers
    }
    
    /// Destructor
    ~gsMultiBasis();

    /// Copy constructor (makes deep copy)
    gsMultiBasis( const gsMultiBasis& other );
    
#if EIGEN_HAS_RVALUE_REFERENCES
    /// Move constructor
    gsMultiBasis(gsMultiBasis&& other)
    {
        m_bases.swap(other.m_bases);
        m_topology.swap(other.m_topology);
    }
#endif

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

    int domainDim () const {return m_bases.front()->domainDim();}
    
    int targetDim () const {return m_bases.front()->targetDim();}

    /// Swap with another gsMultiBasis.
    void swap(gsMultiBasis& other)
    {
        m_topology.swap( other.m_topology );

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

    index_t size() const { GISMO_ERROR("call gsMultiBasis::nBases() instead.");}

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
    const gsBasis<T> & basis(const  std::size_t i ) const
    {
        GISMO_ASSERT( i < m_bases.size(),
                      "Invalid patch index"<<i<<" requested from gsMultiBasis" );
        return *m_bases[i];
    }

    const gsBasis<T> & piece(const index_t i) const 
    {
        GISMO_ASSERT( static_cast<size_t>(i) < m_bases.size(),
                      "Invalid patch index"<<i<<" requested from gsMultiBasis" );
        return *m_bases[i]; 
    }

    /// Return the \a i-th basis block.
    gsBasis<T> & basis(const std::size_t i )
    {
        GISMO_ASSERT( i < m_bases.size(),
                      "Invalid patch index"<<i<<" requested from gsMultiBasis" );
        return *m_bases[i];
    }

    /// @brief Add a basis (ownership of the pointer is also acquired)
    void addBasis( gsBasis<T> * g );

    /// @brief Add a basis (ownership of the pointer is also acquired)
    void addBasis(typename gsBasis<T>::uPtr g);    

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

    /// @brief Refine the component \a comp of every basis uniformly
    /// by inserting \a numKnots new knots on each knot span
    void uniformRefineComponent(int comp, int numKnots = 1, int mul=1)
    {
        for (size_t k = 0; k < m_bases.size(); ++k)
        {
            m_bases[k]->component(comp).uniformRefine(numKnots,mul);
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
        bool sthChanged = false;
        bool change = false;
        do
        {
            sthChanged = false;
            for( size_t i = 0; i < bivec.size(); i++ )
            {
                change = repairInterface( bivec[i] );
                sthChanged = sthChanged || change;
            }
            k++; // just to be sure this cannot go on infinitely
        }
        while( sthChanged && k <= kmax );
    }

    /** @brief Checks if the 2D-interface is fully matching, and if not, repairs it.
    *
    * Same as repairInterface(), but only for 2D and a bit more efficient.
    *
    * \returns true, if something was repaired, i.e., if the mesh on the interface was changed.
    */
    bool repairInterface2d( const boundaryInterface & bi );

    /** @brief Checks if the interface is fully matching, and if not, repairs it.
    *
    * \remarks Designed for gsHTensorBasis and derived bases.
    * Assumes that the respective meshes on all levels of the
    * gsHTensorBasis are fully matching.
    *
    * \returns true, if something was repaired, i.e., if the mesh on the interface was changed.
    */
    bool repairInterface( const boundaryInterface & bi );

    /** @brief Finds the elements that need to be refined in order to repair an interface.
     *
     * This function compares the hierarchical meshes on both patches associated with
     * the boundaryInterface \em bi. The elements that need to be refined on <em>bi.first()</em>
     * and <em>bi.second()</em> are stored in \em refEltsFirst and \em refEltsSecond,
     * respectively.
     *
     * Subsequent calls of the functions\n
     * m_bases[ bi.first().patch ]->refineElements( refEltsFirst )\n
     * m_bases[ bi.second().patch ]->refineElements( refEltsSecond )\n
     * will repair the interface in the sense that the resulting meshes are fully matching
     * (this is done in repairInterface()).
     *
     * \param[in] bi bondaryInterface to be checked.
     * \param[out] refEltsFirst Contains elements (and levels) specifying needed refinement on patch bi.first().
     * \param[out] refEltsSecond Contains elements (and levels) specifying needed refinement on patch bi.second().
     *
     * \returns <em>True</em>, if anything needs to be refined
     * (i.e., if refEltsFirst.size() > 0 or refEltsSecond.size() > 0).\n
     * <em>False</em> if the patches are already fully matching at interface \em bi.
     *
     * Is called by repairInterface(), templated over dimension.
     */
    template<int d>
    bool repairInterfaceFindElements( const boundaryInterface & bi,
                                      std::vector<unsigned> & refEltsFirst,
                                      std::vector<unsigned> & refEltsSecond );

    /// @brief Elevate the degree of every basis by the given amount. (keeping the smoothness)
    void degreeElevate(int const& i = 1, int const dir = -1)
    {
        for (size_t k = 0; k < m_bases.size(); ++k)
            m_bases[k]->degreeElevate(i,dir);
    }

    /// @brief Increase the degree of every basis by the given amount. (keeping the multiplicity)
    void degreeIncrease(int const& i = 1, int const dir = -1)
    {
        for (size_t k = 0; k < m_bases.size(); ++k)
            m_bases[k]->degreeIncrease(i,dir);
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
            m_bases[k]->setDegree(i);
    }

    /// Reduce the continuity by i
    void reduceContinuity(int const i)
    {
        for (size_t k = 0; k < m_bases.size(); ++k)
            m_bases[k]->reduceContinuity(i);
    }

    BasisContainer & patchBases()
    {
        return m_bases;
    }

    const BasisContainer & patchBases() const
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

    void getMapper(iFace::strategy is,
                   const gsBoundaryConditions<T> & bc,
                   gsDofMapper & mapper,
                   int unk,
                   bool finalize = true) const
    { getMapper(is==iFace::glue, bc, unk, mapper, finalize); }

    void getMapper(dirichlet::strategy ds,
                   iFace::strategy is,
                   const gsBoundaryConditions<T> & bc,
                   gsDofMapper & mapper,
                   int unk,
                   bool finalize = true) const
    {
        if ( ds == dirichlet::elimination )
            getMapper(is==iFace::glue, bc, unk, mapper, finalize); 
        else
            getMapper(is==iFace::glue,        mapper, finalize); 
    }
    
    gsDofMapper getMapper(dirichlet::strategy ds,
                          iFace::strategy is,
                          const gsBoundaryConditions<T> & bc,
                          int unk,
                          bool finalize = true) const
    {
        gsDofMapper mapper;
        if ( ds == dirichlet::elimination )
            getMapper(is==iFace::glue, bc, unk, mapper, finalize); 
        else
            getMapper(is==iFace::glue,        mapper, finalize); 
        return mapper;
    }


    // to remove
    void getMapper(bool conforming, gsDofMapper & mapper, bool finalize = true) const;

    void getMapper(iFace::strategy is, gsDofMapper & mapper, bool finalize = true) const
    { getMapper(is==iFace::glue, mapper, finalize); }


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

//    // OUTDATED since implementation of matchWith
//    /**
//     * @brief Matches the degrees of freedom along an interface.
//     *
//     * Same as matchInterface(), but for bases of
//     * class gsHTensorBasis or derived classes.
//     *
//     * \remarks Assumes that the meshes on all levels of
//     * the gsHTensorBasis are fully matching at the interface.
//     *
//     * @param bi specifies the interface to be matched
//     * @param mapper the gsDofMapper which should know that
//     * these interface-DOFs are matched.
//     */
//    template<unsigned d>
//    void matchInterfaceHTensor(const boundaryInterface & bi,
//                               gsDofMapper & mapper) const;

    // Data members
private:

    BasisContainer m_bases;

    gsBoxTopology m_topology;

}; // class gsMultiBasis


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
