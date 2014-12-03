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

    The information on outer boundaries is stored as a list of patch_side structs,
    each one defining the corresponding patch side to lie on the boundary.

    The topological arrangement is stored as a list of
    boundaryInterface structs, each one defining an interface between
    two patch sides.
*/
class GISMO_EXPORT gsBoxTopology
{

public:
    /// Shared pointer for gsBoxTopology
    typedef memory::shared_ptr< gsBoxTopology > Ptr;
    //typedef memory::unique_ptr< gsBoxTopology > LocalPtr;

    typedef std::vector< patch_side >::iterator biterator;
    typedef std::vector< patch_side >::const_iterator const_biterator;

    typedef std::vector< gismo::boundaryInterface >::iterator iiterator;
    typedef std::vector< gismo::boundaryInterface >::const_iterator const_iiterator;

public:

    /// Default constructor
    gsBoxTopology(int d = -1, int n = 0) : m_dim(d), nboxes(n) { }

    gsBoxTopology( int d, unsigned boxes,
            const std::vector< patch_side > & boundary,
            const std::vector< boundaryInterface > & interfaces )
        : m_dim(d), nboxes(boxes), m_boundary(boundary), m_interfaces(interfaces) { }

    // Default copy constructor does the same as the following:
    //gsBoxTopology(const gsBoxTopology & other) : dim(other.dim), nboxes(other.nboxes), 
    //    m_boundary(other.m_boundary), m_interfaces(other.m_interfaces)
    // { }

    ~gsBoxTopology()
    { }

    /// Clone function. Used to make a copy of the object
    gsBoxTopology * clone() const
    {
        return new gsBoxTopology(*this);
    }

public:

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        if ( this->size() > 0 )
        {
            os << "gsBoxTopology (" << this->size() << ").\n";
        }
        else
            os << "gsBoxTopology ( empty! ).\n";

        if ( m_boundary.size() )
            os << "Boundaries:\n";
        for( std::vector< patch_side >::const_iterator bit =
                m_boundary.begin(); bit != m_boundary.end(); ++bit )
        {
            os <<"("<< *bit <<") ";
        }

        if ( m_interfaces.size() )
            os << "\nInterfaces:\n";
        for( std::vector< boundaryInterface >::const_iterator iit =
                m_interfaces.begin(); iit != m_interfaces.end(); ++iit )
        {
            os <<"["<< *iit <<"] ";
        }
        return os;
    }

    /// Print (as string) a boxTopology object
    friend std::ostream& operator<<( std::ostream& os, const gsBoxTopology& b )
    {
        return b.print( os );
    }
    
    /// Number of boxes
    int size () const       { return nboxes; }

    /// Dimension of the boxes
    int dim  () const       { return m_dim; }

    /// Set the dimension of the boxes
    void setDim  (int i)
    { 
        GISMO_ASSERT(m_dim==-1 || i==m_dim, "Changing box dimension.");
        m_dim = i; 
    }

    /// Number of interfaces
    int nInterfaces() const { return m_interfaces.size(); }

    /// Number of boundaries
    int nBoundary() const   { return m_boundary.size(); }

//////////////////////////////////////////////////
// Additional members for Multipatch geometries
//////////////////////////////////////////////////

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
    void clear()
    {
        m_boundary.clear();
        m_interfaces.clear();
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
    void addInterface( int p1, boundary::side s1,
                       int p2, boundary::side s2)
    {
        addInterface(  boundaryInterface( patch_side(p1,s1),patch_side(p2,s2), m_dim ));
    }

    /// Add an interface between side \a s1 of box \a p1 and side \a s2 of box \a p2 with the given orientation.
    void addInterface( int p1, boundary::side s1,
                       int p2, boundary::side s2, const gsVector<bool>& orient)
    {
        addInterface(  boundaryInterface(patch_side(p1,s1),patch_side(p2,s2), orient ));
    }

    /// Add an interface described by \a bi.
    void addInterface( const boundaryInterface& bi )
    {
        m_interfaces.push_back( bi );
    }

    /// Add \a i new boxes.
    void addBox( int i = 1 )
    {
        nboxes +=i;
    }

    /// Set side \a s of box \a p to a boundary.
    void addBoundary( int p, boundary::side s )
    {
        addBoundary( patch_side(p,s) );
    }

    /// Set patch side \a ps to a boundary.
    void addBoundary(const patch_side& ps)
    {
        m_boundary.push_back( ps );
    }

    /// Make all patch sides which are not yet declared as interface or boundary to a boundary.
    void addAutoBoundaries();

    /// Is the given patch side \a ps set to a boundary?
    bool isBoundary(const patch_side& ps) const
    {
        return std::find(m_boundary.begin(), m_boundary.end(), ps) != m_boundary.end();
    }

    /// Is the given patch side \a ps set to an interface?
    bool isInterface(const patch_side& ps) const;

    /// Return the vector of boundaries.
    std::vector< patch_side >        boundaries() const { return m_boundary;   }

    /// Return the vector of interfaces.
    std::vector< boundaryInterface > interfaces() const { return m_interfaces; }

    /// Check that boundaries and interfaces are consistent.
    void checkConsistency() const;

    /// Iteration: set \a result to the first patch side of the first box.
    void firstPatchSide(patch_side& result);

    /// Iteration: increment \a result to the next patch side, iterating over all sides of all boxes.
    bool nextPatchSide(patch_side& result);

    /// Access i-th boundary interface
    const boundaryInterface & bInterface(int i) const {return m_interfaces[i];}

    /// set \a result to the associated patchside of \a ps, returns false if it is a boundary patch_side
    bool getNeighbour(const patch_side& ps ,patch_side& result, boundaryInterface* iface=NULL) const;

    /// set \a result to the associated interface of \a ps, returns false if it is a boundary patch_side
    bool getInterface(const patch_side& ps,boundaryInterface & result) const
    {
        for ( unsigned i = 0; i < m_interfaces.size(); ++i )
            if ( m_interfaces[i].ps1 == ps || m_interfaces[i].ps2 == ps )
            {
                result = m_interfaces[i];
                return true;
            }
        return false;
    }

    /// set \a result to the orientation vector of the interface, returns false if it is a boundary patch_side
    bool getOrientationOfInterface(const patch_side& ps, gsVector<bool> & result) const
    {
        for ( unsigned i = 0; i < m_interfaces.size(); ++i )
            if ( m_interfaces[i].ps1 == ps || m_interfaces[i].ps2 == ps )
            {
                result = m_interfaces[i].orient();
                return true;
            }
        return false;
    }

    /// takes a patch_corner \a start and gives back all other patch_corners,
    /// that represent the same point in the vector \a cornerList
    /// CAREFUL: works for 2D only
    bool getCornerList(const patch_corner& start,std::vector<patch_corner> & cornerList) const;

    /// gives back all the extraordinary vertices (3 faces or more than 4) of the topology
    /// each EV is represented by a vector of patch_corners, which represent the same vertex
    /// all the vectors are put in the vector \a cornerLists
    /// CAREFUL: works for 2D only
    void getEVs(std::vector<std::vector<patch_corner> > & cornerLists) const;

    /// gives back all the ordinary vertices (4 faces) of the topology
    /// each OV is represented by a vector of patch_corners, which represent the same vertex
    /// all the vectors are put in the vector \a cornerLists
    /// CAREFUL: works for 2D only
    void getOVs(std::vector<std::vector<patch_corner> > & cornerLists) const;

protected:
    // Data members

    int m_dim;
    int nboxes;
    std::vector< patch_side > m_boundary;
    std::vector< boundaryInterface > m_interfaces ;

}; // class gsBoxTopology


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsBoxTopology.cpp)
#endif
