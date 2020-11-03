/** @file gsPlanarDomain.h

    @brief Provides declaration of gsPlanarDomain class. The outer boundary
    (m_loops[0]) is a loop of curves, listed in anticlockwise order. Inner
    boundaries (m_loops[i], i > 0) are loops of curves, listed in clockwise
    order.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Falini, A. Mantzaflaris, D. Mayer, D.-M. Nguyen, M. Pauley
*/

#pragma once

#include <iostream>
#include <deque>

#include <gsModeling/gsCurveLoop.h> 
#include <gsModeling/gsTemplate.h>

// #include <gsUtils/gsSortedVector.h>


namespace gismo
{
template <class T>
class gsBemSolution;

/**
      @brief Class representing a Planar domain with an outer boundary and a number of holes.

      The outer boundary is oriented CCW and the holes are CW.
      
      \ingroup Modeling
  */

// TODO Add gsWriteParaview using PolygonalData....
template<class T>
class gsPlanarDomain
{
public:

    typedef memory::shared_ptr<gsPlanarDomain> Ptr;
    typedef memory::unique_ptr<gsPlanarDomain> uPtr;
    
    /// Default empty constructor
    gsPlanarDomain() { }

    /// Construct planar domain by giving the outer boundary
    gsPlanarDomain( gsCurveLoop<T> * boundary)
    {
        if ( boundary->is_ccw() )
        {
            m_loops.push_back(boundary);
        }
        else
        {
            boundary->reverse();
            m_loops.push_back(boundary);
        }

        updateBoundingBox();
    }

    /// Construct planar domain by a list of loops
    gsPlanarDomain( std::vector< gsCurveLoop<T> *> const & loops);

    /// Construct planar domain by an outer boundary given by a curve
    gsPlanarDomain( gsCurve<T> * boundary)
    {
        m_loops.push_back( new gsCurveLoop<T>(boundary) );
        updateBoundingBox();
    }

    ~gsPlanarDomain()
    {
        freeAll( m_loops );
    }

    /// Copy constructor
    gsPlanarDomain( const gsPlanarDomain & other)
        : m_loops( other.m_loops.size() )
    {
        m_bbox = other.m_bbox;
        cloneAll( other.m_loops.begin(), other.m_loops.end(),
                  this->m_loops.begin() );
    }


    /// Assignment operator
    gsPlanarDomain& operator=( const gsPlanarDomain& other)
    {
        freeAll( m_loops );
        m_loops.resize( other.m_loops.size() );
        m_bbox = other.m_bbox;
        cloneAll( other.m_loops.begin(), other.m_loops.end(),
                  this->m_loops.begin() );
        return *this;
    }

    /// Clone function. Used to make a copy of the (derived) geometry
    //GISMO_CLONE_FUNCTION(gsPlanarDomain)
    uPtr clone() const
    {
        return uPtr(new gsPlanarDomain(*this));
    }

public:

    void check()
    {
        if ( !m_loops[0]->is_ccw() )
            gsWarn<< "Wrong orientation in outer loop of planar domain.";
        for(size_t i=1; i< m_loops.size(); i++)
        {
            if( m_loops[i]->is_ccw())
                gsWarn<< "Wrong orientation in loop["<< i <<"] of planar domain.";
        }
    }

    void insertHole( gsCurveLoop<T> * hole )
    {
        if ( hole->is_ccw() )
        {
            hole->reverse();
        }
        m_loops.push_back( hole );
    }

    int numLoops() const { return m_loops.size();    }
    int numHoles() const { return m_loops.size() -1; }

    gsCurveLoop<T> & outer()             { return loop(0); }
    const gsCurveLoop<T> & outer() const { return loop(0); }

    gsCurveLoop<T> & loop(unsigned loopNumber)
    {
        GISMO_ASSERT( loopNumber<m_loops.size(), "Loop does not exist" );
        return *m_loops[loopNumber];
    }
    const gsCurveLoop<T> & loop(unsigned loopNumber) const
    {
        GISMO_ASSERT( loopNumber<m_loops.size(), "Loop does not exist" );
        return *m_loops[loopNumber];
    }

    gsCurve<T> & curve(unsigned loopNumber, unsigned curveNumber)
    {
        GISMO_ASSERT( loopNumber<m_loops.size(), "Loop does not exist" );
        return m_loops[loopNumber]->curve(curveNumber);
    }
    const gsCurve<T> & curve(unsigned loopNumber, unsigned curveNumber) const
    {
        GISMO_ASSERT( loopNumber<m_loops.size(), "Loop does not exist" );
        return m_loops[loopNumber]->curve(curveNumber);
    }

    bool contains( gsVector<T> const &, T)
    {
        GISMO_NO_IMPLEMENTATION
    }

    gsMatrix<T> boundingBox() const
    {
        gsMatrix<T> res(2,2);

        res = m_loops[0]->getBoundingBox();
        return res;
    }

    void translate(gsVector<T> const & v)
    {
        for ( typename std::vector< gsCurveLoop<T> *>::iterator it =
              m_loops.begin();  it != m_loops.end(); ++it)
            (*it)->translate(v);
    }

    //gsMatrix<T> averageValue( std::vector<gsFunction<T>*> const &f, std::vector<T> const & breaks);

    /// @name inDomain
    /// given a matrix of points \param u, returns true if they are inside the planar domain
    ///\param direction sets to 0 states we are performing our checking running parallel to the x-axis
    bool inDomain( gsMatrix<T> const & u, int direction = 0);

    ///@name onBoundary
    ///given a matrix of points \param u, returns true if they lie on the boundary of the planar domain
    bool onBoundary(gsMatrix<T> const & u);

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const;

    friend std::ostream& operator<<(std::ostream& os, const gsPlanarDomain& pd) { return pd.print(os); }

    /// linearly discriti
    void sampleLoop_into( int loopID, int npoints, int numEndPoints, gsMatrix<T> & u );

    gsMatrix<T> sampleLoop( int loopID, int npoints=50, int numEndPoints=2)
    {
        gsMatrix<T> u;
        sampleLoop_into( loopID, npoints, numEndPoints, u );
        return u;
    }

    void sampleCurve_into( int loopID, int curveID, int npoints, gsMatrix<T> & u );

    gsMatrix<T> sampleCurve(int loopID, int curveID, int npoints = 50)
    {
        gsMatrix<T> u;
        sampleCurve_into( loopID, curveID, npoints, u );
        return u;
    }

    /// Return a triangulation of the planar domain
    memory::unique_ptr<gsMesh<T> > toMesh(int npoints = 50) const;     // FOR NOW ONLY ONE LOOP

    /// split this planar domain in two, returning the new planar domain created
    /// as a result.
    uPtr split(int startIndex, int endIndex,
                          gsCurve<T> * newCurveThisFace, gsCurve<T> * newCurveNewFace)
    {
        typename gsCurveLoop<T>::uPtr newCurveLoop = this->m_loops[0]->split(startIndex, endIndex, newCurveThisFace, newCurveNewFace);
        updateBoundingBox();

        return uPtr(new gsPlanarDomain<T>(newCurveLoop.release()));
    }

    /// Update the bounding box. Needs to be called after any operation that
    /// modifies the outer loop.
    // TODO: make this function private or protected. In order to do this we
    // would need to replace the functions "outer", "loop" and "curve" with
    // new public functions for modifying the loops that call this function
    // when they are done.
    void updateBoundingBox()
    {
        assert(!m_loops.empty()); // outer loop does not exist
        m_bbox = m_loops[0]->getBoundingBox();
    }

    /// split the \a curveId^th curve in the \a loopId^th loop of the planar domain into two curves
    /// \param loopId specifies the loop
    /// \param curveId specifies the curve in the loop
    /// \param lengthRatio   ratio of the lengths of the first new curve and of the original curve
    gsMatrix<T> splitCurve(size_t loopId, size_t curveId, T lengthRatio=.5)
    {
        return m_loops[loopId]->splitCurve(curveId,lengthRatio);
    }



    // Data members
private:

    // m_loops[0] is CCW, all others CW holes
    std::vector< gsCurveLoop<T> *> m_loops;

    // The lower left and upper right corner of a square bounding the
    // domain
    gsMatrix<T,2,2> m_bbox;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

}; // class gsPlanarDomain

} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPlanarDomain.hpp)
#endif
