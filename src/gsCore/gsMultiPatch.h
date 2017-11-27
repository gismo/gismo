/** @file gsMultiPatch.h

    @brief Provides declaration of the MultiPatch class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsBoxTopology.h>
#include <gsCore/gsGeometry.h>

#include <gsCore/gsMultiBasis.h>

namespace gismo
{

/** @brief Container class for a set of geometry patches and their
    topology, that is, the interface connections and outer boundary
    faces.

    \tparam T coefficient type
    
    \ingroup Core
*/
template<class T>
class gsMultiPatch : public gsBoxTopology, public gsFunctionSet<T>
{

public:
    typedef gsBoxTopology Base;
    /// Shared pointer for gsMultiPatch
    typedef memory::shared_ptr< gsMultiPatch > Ptr;
    /// Unique pointer for gsMultiPatch
    typedef memory::unique_ptr< gsMultiPatch > uPtr;

    typedef std::vector<gsGeometry<T> *> PatchContainer;

    typedef typename PatchContainer::size_type size_t;
    typedef typename PatchContainer::iterator iterator;
    typedef typename PatchContainer::const_iterator const_iterator;

public:

    /// Default empty constructor
    gsMultiPatch() { }

    /// Copy constructor (makes deep copy)
    gsMultiPatch( const gsMultiPatch& other );

#if ! EIGEN_HAS_RVALUE_REFERENCES

    /// Assignment operator (uses copy-and-swap idiom)
    gsMultiPatch& operator= ( gsMultiPatch other )
    {
        this->swap(other);
        return *this;
    }

#else
    /// Move constructor
    gsMultiPatch( gsMultiPatch&& other )
        : gsBoxTopology( give(other) ), m_patches( give(other.m_patches) )
    {}

    /// Assignment operator 
    gsMultiPatch& operator= ( const gsMultiPatch& other );

    /// Move assignment operator 
    gsMultiPatch& operator= ( gsMultiPatch&& other );

#endif

    /// Create from a vector of patches
    explicit gsMultiPatch( PatchContainer & patches );

    /// Create a single-patch instance
    gsMultiPatch( const gsGeometry<T> & geo );

    /// Create from patches and boundary/interface information
    gsMultiPatch( PatchContainer & patches,
                  const std::vector<patchSide>& boundary,
                  const std::vector<boundaryInterface>& interfaces );
    
    /// Destructor
    ~gsMultiPatch();

    GISMO_CLONE_FUNCTION(gsMultiPatch)

public:

    /// Get a const-iterator to the patches
    /// \return an iterator to the beginning of the patches
    const_iterator begin() const 
    { return m_patches.begin(); }

    /// Get a const iterator to the end of the patches
    /// \return an iterator to the end of the patches
    const_iterator end() const 
    { return m_patches.end(); }

    /// Get an iterator to the beginning of the  patches
    /// \return an iterator to the beginning of the  patches
    iterator begin() 
    { return m_patches.begin(); }

    /// Get an iterator to the end of the  patches
    /// \return an iterator to the end of the  patches
    iterator end()  
    { return m_patches.end(); }

public:
    
    const gsGeometry<T> & piece(const index_t i) const
    { return *m_patches[i]; }

    virtual index_t size() const
    { return m_patches.size(); }

public:
    /**
     * @brief construct the affine map that places bi.first() next to bi.second() and
     *        identifies the two matching sides.
     *        The optional scaling specifies the scaling of first() in the direction
     *        orthogonal to the interface.
     *        By default it preserves the size of bi.first() in that direction
     * @param bi
     * @param scaling
     * @return
     */
    gsAffineFunction<T> getMapForInterface(const boundaryInterface &bi, T scaling=0) const;

    /// \brief Swap with another gsMultiPatch
    void swap(gsMultiPatch& other)
    {
        gsBoxTopology::swap( other );
        m_patches.swap( other.m_patches );
    }

    /// \brief Prints the object as a string
    std::ostream& print( std::ostream& os ) const;

    /// \brief Prints the object as a string with extended details
    std::string detail() const;

    /// \brief Dimension of the parameter domain (must match for all patches).
    int parDim() const 
    {
        //GISMO_ASSERT( m_patches.size() > 0 , "Empty multipatch object.");
        return m_dim;
    }
    int domainDim () const {return parDim();}
    
    /// \brief Dimension of the geometry (must match for all patches).
    int geoDim() const;
    int targetDim () const {return geoDim();}
    
    /// \brief Co-dimension of the geometry (must match for all patches).
    int coDim() const;

    /// \brief Returns true if the multipatch object is a closed
    /// manifold (ie. it has no boundaries)
    bool isClosed() { return this->nBoundary() == 0; }

    /// \brief Returns the range of parameter
    gsMatrix<T> parameterRange(int i = 0) const;

    /// \brief Number of patches
    std::size_t nPatches() const          { return m_patches.size(); }

    /// \brief Returns a vector of patches // to do : replace by copies
    PatchContainer const& patches() const { return m_patches; }

    const gsBoxTopology & topology() const { return *this; }
    
    /// \brief Return the \a i-th patch.
    const gsGeometry<T> & operator []( size_t i ) const { return *m_patches[i]; }

    /// \brief Makes a deep copy of all bases and puts them in a vector
    std::vector<gsBasis<T> *> basesCopy() const;

    /// Return the \a i-th patch.
    gsGeometry<T>& patch( std::size_t i ) const
    {
        GISMO_ASSERT( i < m_patches.size(), "Invalid patch index requested from gsMultiPatch" );
        return *m_patches[i];
    }

    ///\brief Permutes the patches according to \a perm
    void permute(const std::vector<int> & perm);

    ///\brief Return the basis of the \a i-th patch.
    gsBasis<T> & basis( std::size_t i ) const;

    ///\brief Add a patch (pointer is consumed)
    void addPatch( gsGeometry<T> * g );

    ///\brief Add a patch from a gsGeometry<T>::uPtr
    void addPatch(typename gsGeometry<T>::uPtr g);

    /// Add a patch by copying argument
    void addPatch(const gsGeometry<T> & g);

    /// \brief Search for the given geometry and return its patch index.
    int findPatchIndex( gsGeometry<T>* g ) const;

    /// @brief Add an interface joint between side \a s1 of geometry
    /// \a g1 side \a s2 of geometry \a g2.
    ///
    /// \todo add orientation information
    GISMO_DEPRECATED void addInterface( gsGeometry<T>* g1, boxSide s1,
            gsGeometry<T>* g2, boxSide s2 );

    using Base::addInterface; // unhide base function
    
    /// Add side s of patch g to the outer boundary of the domain
    void addPatchBoundary( gsGeometry<T>* g, boxSide s ) {
        int p = findPatchIndex( g );
        gsBoxTopology::addBoundary( patchSide( p, s ) );
    }

    /// \brief Refine uniformly all patches by inserting \a numKnots
    /// in each knot-span with multipliplicity \a mul
    void uniformRefine(int numKnots = 1, int mul = 1);

    /// \brief Elevate the degree of all patches by \a elevationSteps.
    void degreeElevate(int elevationSteps = 1);

    void embed(const index_t N)
    {
        for ( typename PatchContainer::const_iterator it = m_patches.begin();
              it != m_patches.end(); ++it )
        {
            ( *it )->embed(N);
        }
    }
    
    /// \brief Attempt to compute interfaces and boundaries
    /// automatically.
    /// \param tol The tolerance to test for matching points
    /// \param cornersOnly When set to true an interface is accepted
    /// if the patch corners match, even if the parameterization does
    /// not agree
    bool computeTopology( T tol = 1e-4, bool cornersOnly = false);

    /// \brief Attempt to close gaps between the interfaces. Assumes
    /// that the topology is computed, ie. computeTopology() has been
    /// called.
    void closeGaps( T tol = 1e-4 );

    /// Clear (delete) all patches
    void clear()
    {
        Base::clearAll();
        freeAll(m_patches);
        m_patches.clear();
    }

    /// \brief Returns a bounding box for the multipatch domain. The
    /// output \a result is a matrix with two columns, corresponding
    /// to two points: the lower and upper corner of the bounding box.
    void boundingBox(gsMatrix<T> & result) const;

    /// \brief Splits each patch uniformly in each direction into two new patches,
    /// giving a total number of 2^d new patches. This method allocated new
    /// space for each new geometry, the original one stays unchanged.
    gsMultiPatch<T> uniformSplit() const;


    /** @brief Checks if all patch-interfaces are fully matching, and if not, repairs them, i.e., makes them fully matching.
    *
    * \remarks Designed for gsHTensorBasis and derived bases.
    * Assumes that the meshes on all levels of the gsHTensorBasis
    * are fully matching.
    */
    void repairInterfaces();

protected:

    void setIds();

    // Data members
private:

    PatchContainer m_patches;

private:
    // implementation functions

    // match the vertices in ci1 starting from start to the end with the vertices
    // in ci2 that are still non matched
    // cc1 and cc2 are the physical coordinates of the vertices
    // ci1 and ci2 are the indices corner indices
    // start is the index in ci1 of the vertex to match
    // reference is the image of the vertex ci1(0) and it is used to compute orientation
    // tol is the allowed distance between two vertexes in physical domain
    // matched keeps track of the already matched vertices
    // dirMap and dirO are the output orientation of the match
    // return true if all the vertices starting from start are matched
    static bool matchVerticesOnSide (
           const gsMatrix<T> &cc1, const std::vector<boxCorner> &ci1, index_t start,
           const gsMatrix<T> &cc2, const std::vector<boxCorner> &ci2, 
           const gsVector<bool> &matched,
           gsVector<index_t> &dirMap, gsVector<bool>    &dirO,
           T tol, index_t reference=0);

}; // class gsMultiPatch


/// Print (as string) a multipatch structure
template<class T>
std::ostream& operator<<( std::ostream& os, const gsMultiPatch<T>& b )
{
    return b.print( os );
}


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMultiPatch.hpp)
#endif
