/** @file gsBoxTopology.h

    @brief Provides declaration of the BoxTopology class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): F. Buchegger, A. Mantzaflaris
*/


#pragma once

#include <gsCore/gsExport.h>

#include <gsCore/gsBoundary.h>

namespace gismo
{

/** @brief
    Defines a topological arrangement of a collection of "boxes"
    (e.g., parameter domains that map to physical-domain patches).

    The information on outer boundaries is stored as a list of patchSide structs,
    each one defining the corresponding patch side to lie on the boundary.

    The topological arrangement is stored as a list of
    boundaryInterface structs, each one defining an interface between
    two patch sides.
    
    \ingroup Core
*/
class GISMO_EXPORT gsBoxTopology
{

public:
    /// Shared pointer for gsBoxTopology
    typedef memory::shared_ptr< gsBoxTopology > Ptr;
    typedef memory::unique_ptr< gsBoxTopology > uPtr;

    typedef std::vector< patchSide > bContainer;
    typedef bContainer::iterator biterator;
    typedef bContainer::const_iterator const_biterator;

    typedef std::vector<boundaryInterface> ifContainer;
    typedef ifContainer::iterator iiterator;
    typedef ifContainer::const_iterator const_iiterator;

    typedef const boundaryInterface * InterfacePtr;
public:

    /// Default constructor
    gsBoxTopology(short_t d = -1, index_t n = 0) : m_dim(d), nboxes(n) { }

    gsBoxTopology( short_t d, index_t boxes,
            const bContainer & boundary,
            const ifContainer & interfaces )
        : m_dim(d), nboxes(boxes), m_boundary(boundary), m_interfaces(interfaces) { }

    // Default copy constructor does the same as the following:
    //gsBoxTopology(const gsBoxTopology & other) : dim(other.dim), nboxes(other.nboxes), 
    //    m_boundary(other.m_boundary), m_interfaces(other.m_interfaces)
    // { }
    
    /// Clone function. Used to make a copy of the object
    gsBoxTopology * clone() const
    {
        return new gsBoxTopology(*this);
    }

public:

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const;

    /// Print (as string) a boxTopology object
    friend std::ostream& operator<<( std::ostream& os, const gsBoxTopology& b )
    {
        return b.print( os );
    }
    
    /// Number of boxes
    index_t nBoxes() const       { return nboxes; }

    /// Dimension of the boxes
    short_t dim  () const       { return m_dim; }

    /// Set the dimension of the boxes
    void setDim  (short_t i)
    { 
        GISMO_ASSERT(m_dim==-1 || i==m_dim, "Changing box dimension.");
        m_dim = i; 
    }

    /// Number of interfaces
    size_t nInterfaces() const { return m_interfaces.size(); }

    /// Number of boundaries
    size_t nBoundary() const   { return m_boundary.size(); }

/*
 * Additional members for Multipatch geometries
 */

	/// Get a const-iterator to the interfaces
	/// \return an iterator to the beginning of the interfaces
	const_iiterator iBegin() const
	{ return m_interfaces.begin(); }

	/// Get a const iterator to the end of the interfaces
	/// \return an iterator to the end of the interfaces
	const_iiterator iEnd() const
	{ return m_interfaces.end(); }

	/// Get an iterator to the beginning of the  interfaces
	/// \return an iterator to the beginning of the  interfaces
	iiterator iBegin()
	{ return m_interfaces.begin(); }

	/// Get an iterator to the end of the  interfaces
	/// \return an iterator to the end of the  interfaces
	iiterator iEnd()
	{ return m_interfaces.end(); }

	/// Get a const-iterator to the beginning of the boundaries
	/// \return an iterator to the beginning of the boundaries
	const_biterator bBegin() const
	{ return m_boundary.begin(); }

	/// Get a const-iterator to the end of the boundaries
	/// \return an iterator to the end of the boundaries
	const_biterator bEnd() const
	{ return m_boundary.end(); }

	/// Get an iterator to the beginning of the boundaries
	/// \return an iterator to the beginning of the boundaries
	biterator bBegin()
	{ return m_boundary.begin(); }

	/// Get an iterator to the end of the boundaries
	/// \return an iterator to the end of the knotvector
	biterator bEnd()
	{ return m_boundary.end(); }

    /// Clear all boundary and interface data.
    void clearTopology()
    {
        m_boundary  .clear();
        m_interfaces.clear();
    }

    /// Clear all boxes, boundary and interface data.
    void clearAll()
    {
        clearTopology();
        m_dim  = -1;
        nboxes =  0;
    }

    /// Swap with another gsBoxTopology.
    void swap(gsBoxTopology& other)
    {
        std::swap( m_dim, other.m_dim );
        std::swap( nboxes, other.nboxes );
        m_boundary.swap( other.m_boundary );
        m_interfaces.swap( other.m_interfaces );
    }

    /// Add an interface between side \a s1 of box \a p1 and side \a s2 of box \a p2.
    void addInterface(index_t p1, boxSide s1,
                      index_t p2, boxSide s2)
    {
        addInterface(boundaryInterface(patchSide(p1, s1), patchSide(p2, s2), m_dim));
    }

    /// Add an interface described by \a bi.
    void addInterface( const boundaryInterface& bi )
    {
        m_interfaces.push_back( bi );
    }

    /// Add \a i new boxes.
    void addBox(index_t i = 1)
    {
        nboxes += i;
    }

    /// Set side \a s of box \a p to a boundary.
    void addBoundary(index_t p, boxSide s)
    {
        addBoundary(patchSide(p, s));
    }

    /// Set patch side \a ps to a boundary.
    void addBoundary(const patchSide& ps)
    {
        m_boundary.push_back( ps );
    }

    /// Make all patch sides which are not yet declared as interface or boundary to a boundary.
    void addAutoBoundaries();

    /// Is the given patch side \a ps set to a boundary?
    bool isBoundary(const patchSide& ps) const
    {
        return std::find(m_boundary.begin(), m_boundary.end(), ps) != m_boundary.end();
    }

    /// Returns true if side \a s on patch \a p is a boundary.
    bool isBoundary( index_t p, boxSide s )
    {
        return isBoundary( patchSide(p,s) );
    }

    /// Is the given patch side \a ps set to an interface?
    bool isInterface(const patchSide& ps) const;

    /// Return the vector of boundaries.
    const bContainer & boundaries() const { return m_boundary;}
    bContainer & boundaries() { return m_boundary;}

    /// Return the vector of interfaces.
    const ifContainer & interfaces() const { return m_interfaces; }
    ifContainer & interfaces() { return m_interfaces; }

    /// Check that boundaries and interfaces are consistent.
    void checkConsistency() const;

    /// Access i-th boundary interface
    const boundaryInterface & bInterface(int i) const {return m_interfaces[i];}

    /// set \a result to the associated patchSide of \a ps, returns
    /// false if it is a boundary patchSide
    bool getNeighbour(const patchSide& ps ,patchSide& result, int & ii) const;

    /// set \a result to the associated patchSide of \a ps, returns
    /// false if it is a boundary patchSide
    bool getNeighbour(const patchSide& ps ,patchSide& result) const;

    /// Returns a pointer to the interface between boxes \a b1 and \a
    /// b2, if one exists, otherwise it returns a null pointer
    InterfacePtr findInterface(const index_t b1, const index_t b2) const;

    /// set \a result to the associated interface of \a ps, returns
    /// false if it is a boundary patchSide
    bool getInterface(const patchSide& ps,boundaryInterface & result) const
    {
        for ( unsigned i = 0; i < m_interfaces.size(); ++i )
            if ( m_interfaces[i].first() == ps || m_interfaces[i].second() == ps )
            {
                result = m_interfaces[i];
                return true;
            }
        return false;
    }

    /// takes a patchCorner \a start and gives back all other patchCorners,
    /// that represent the same point in the vector \a cornerList
    bool getCornerList(const patchCorner& start, std::vector<patchCorner> & cornerList) const;

    /// returns the maximal valence of a vertex of this topology.
    int getMaxValence() const;

    /// @brief Returns all components representing the topology
    ///
    /// Each entry of the outer vector represents one component (patch-interior, face,
    /// edge, corner, etc.). Since the components refering to one interface can be
    /// addressed as belonging to different patches, each component itself is represented
    /// by an inner vector which contains all \a patchComponent objects that refer
    /// to the particular component.
    ///
    /// @param combineCorners If this is set, all corners are treated as one component
    std::vector< std::vector<patchComponent> > allComponents(bool combineCorners = false) const;

    /// gives back all the extraordinary vertices (3 faces or more than 4) of the topology
    /// each EV is represented by a vector of patchCorners, which represent the same vertex
    /// all the vectors are put in the vector \a cornerLists. It will only find vertices on
    /// the inside.
    /// CAREFUL: works only for 2D
    void getEVs(std::vector<std::vector<patchCorner> > & cornerLists) const;

    /// gives back all the ordinary vertices (4 faces) of the topology
    /// each OV is represented by a vector of patchCorners, which represent the same vertex
    /// all the vectors are put in the vector \a cornerLists It will only find vertices on
    /// the inside.
    /// CAREFUL: works only for 2D
    void getOVs(std::vector<std::vector<patchCorner> > & cornerLists) const;

protected:
    // Data members

    /// Dimension of the boxes held
    short_t m_dim;

    /// Number of boxes held
    index_t nboxes;

    /// List of boundaries of the boxes
    bContainer m_boundary;

    /// List of intefaces between boxes
    ifContainer m_interfaces ;

}; // class gsBoxTopology


} // namespace gismo

