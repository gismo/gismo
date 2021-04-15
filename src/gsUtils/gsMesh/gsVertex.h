/** @file gsVertex.h

    @brief Provides gsVertex class for a vertex of a gsMesh

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Mayer
*/

#pragma once

#include <gsUtils/gsMesh/gsMeshElement.h>
#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{

/**
 * @brief gsVertex class that represents a 3D vertex for a gsMesh.
 */
template <class T>
class gsVertex  : public gsMeshElement<T>, public gsVector3d<T>
{
public:
    /// @brief Shared pointer for gsVertex
    typedef memory::shared_ptr< gsVertex > Ptr;

    /// @brief Unique pointer for gsVertex
    typedef memory::unique_ptr< gsVertex > uPtr;

    typedef typename gsVector3d<T>::Scalar       scalar_t;
    typedef gsMeshElement<T>                     MeshElement;
    typedef typename MeshElement::gsFaceHandle   gsFaceHandle;
    typedef typename MeshElement::gsVertexHandle gsVertexHandle;

public:
    /// @brief Constructor
    gsVertex() : MeshElement(), gsVector3d<T>() { }

    template<typename OtherDerived>
    gsVertex(const Eigen::MatrixBase<OtherDerived>& other) :
    MeshElement(), gsVector3d<T>(other) { }

    /// @brief Constructor, take 3 scalars.
    /// \param x, y, z Coordinates of position in 3D space.
    gsVertex(scalar_t x, scalar_t y, scalar_t z = 0) :
        MeshElement(), gsVector3d<T>(x,y,z),sharp(0)
    { }

    /// @brief Constructor, takes a gsVector3d
    /// \param u the gsVector3d
    gsVertex( gsVector3d<T> const & u) :
        MeshElement(), gsVector3d<T>(u),sharp(0)
    { }

    /// @brief Constructor, takes a gsVector.
    /// \param u gsVector of dimension 1, 2 or 3. Fills with zero.
    gsVertex( gsVector<T> const & u) :
        MeshElement(), gsVector3d<T>(),sharp(0)
    {
        // vertex is always a 3-cooordinate vector.
        const index_t r = u.rows();
        GISMO_ASSERT(r<=3, "Invalid input in gsVertex constructor");

        if ( r < 3 )
            this->tail(3-r).setZero();

        this->head(r) = u;
    }

    virtual ~gsVertex() { };

    gsVertex & operator=(const gsVertex & other)
    {
        this->gsVector3d<T>::operator=(other);
        return *this;
    }
    
    //GISMO_CLONE_FUNCTION(gsVertex)
    /// @brief Clone Function (deep copy)
    uPtr clone() const { return uPtr(new gsVertex(*this)); }

    /// @brief Moves a gsVertex relatively.
    /// \param dx, dy, dz values added to x, y and z
    void move(scalar_t dx, scalar_t dy, scalar_t dz)
    {
        this->x() += dx;
        this->y() += dy;
        this->z() += dz;
    }

    /// @brief Adds a gsFaceHandle \a f to the list of faces adjaent to this vertex.
    /// \param f gsFaceHandle
    inline void addFace(gsFaceHandle const& f)
    { faces.push_back( f ); }

    /// \return x
    inline T   x () const { return (*this)(0); }
    /// \return y
    inline T   y () const { return (*this)(1); }
    /// \return z
    inline T   z () const { return (*this)(2); }
    /// \return &x
    inline T & x () { return (*this)(0); }
    /// \return &y
    inline T & y () { return (*this)(1); }
    /// \return &z
    inline T & z () { return (*this)(2); }

    std::ostream &print(std::ostream &os) const
    {
        os << "Vertex( " << this->x() << " " << this->y() << " " << this->z() << " )\n";
        return os;
    }

public:

    std::vector<gsVertexHandle> nVertices;
    //gsVector3d<T> coords;
    bool sharp;

    /// @brief List of faces adjacent to this vertex
    std::vector<gsFaceHandle> faces;
    int numEdges;

    T data;
};

  template<class T>
  static bool Xless( typename gsVertex<T>::gsVertexHandle const & a,
		     typename gsVertex<T>::gsVertexHandle const & b )
  {return a->x()< b->x() || ( a->x()==b->x() && a->y()<b->y() )
      || ( a->x()==b->x() && a->y()==b->y() && a->z()<b->z()); }

  template<class T>
  static bool Yless( typename gsVertex<T>::gsVertexHandle const & a,
		     typename gsVertex<T>::gsVertexHandle const & b )
  {return a->y()< b->y() || ( a->y()==b->y() && a->x()<b->x() ); }
  template<class T>
  static bool Zless( typename gsVertex<T>::gsVertexHandle const & a,
		     typename gsVertex<T>::gsVertexHandle const & b )
  {return a->z()< b->z(); }

// Lexicographic comparison for vertex handles

template<class T>
struct lexCompareVHandle
{
    bool operator() (typename gsVertex<T>::gsVertexHandle const & lhs,
		     typename gsVertex<T>::gsVertexHandle const & rhs) const
      {return lhs->x()< rhs->x() || ( lhs->x()==rhs->x() && lhs->y()<rhs->y() )
      || ( lhs->x()==rhs->x() && lhs->y()==rhs->y() && lhs->z()<rhs->z()); }
};

template<class T>
T length(gsVertex<T> const & vert)
{
    //return (sqrt(vert.x()*vert.x()+vert.y()*vert.y()+vert.z()*vert.z()));
    return vert.norm();
}

template<class T>
bool operator < (typename gsVertex<T>::gsVertexHandle const & lhs,
		typename gsVertex<T>::gsVertexHandle const & rhs)
{
    return !(lhs->x() < rhs->x() || (lhs->x() == rhs->x() && lhs->y() < rhs->y())
        || (lhs->x() == rhs->x() && lhs->y() == rhs->y() && lhs->z() < rhs->z()));
}
template<class T>
bool operator > (typename gsVertex<T>::gsVertexHandle const & lhs,
        typename gsVertex<T>::gsVertexHandle const & rhs)
{
    return !(lhs->x() > rhs->x() || (lhs->x() == rhs->x() && lhs->y() > rhs->y())
        || (lhs->x() == rhs->x() && lhs->y() == rhs->y() && lhs->z() > rhs->z()));
}
//template<class T>
//bool operator == (typename gsVertex<T> const & lhs,
//        typename gsVertex<T> const & rhs)
//{
//  return (lhs->x()== rhs->x()&& lhs->y()==rhs->y()&& lhs->z()==rhs->z())
//    ;}

/**
 * Compare LHS == RHS
 * @tparam T
 * @param lhs
 * @param rhs
 * @return
 */
template<class T>
bool operator == (gsVertex<T> const & lhs,gsVertex<T> const & rhs)
{
    return (lhs.x()==rhs.x())&&
           (lhs.y()==rhs.y())&&
           (lhs.z()==rhs.z());
//    return lhs.Eigen::template Matrix<T,3,1>::operator==(rhs); /slower
}

template<class T>
bool operator < (gsVertex<T> const & lhs,gsVertex<T> const & rhs)
{
    return (lhs.x()> rhs.x() || ( lhs.x()==rhs.x() && lhs.y()>rhs.y() )
            || ( lhs.x()==rhs.x() && lhs.y()==rhs.y() && lhs.z()>rhs.z()));
}

/**
 * Compares LHS.x < RHS.X
 * Compares first x, if equal y, if equal z.
 * @tparam T
 * @param lhs
 * @param rhs
 * @return true if lhs.x < rhs.x; lhs.x==rhs.x and lhs.y < rhs.y; lhs.x==rhs.x, lhs.y==rhs.y and lhs.z < rhs.z;
 * false if lhs.x >= rhs.x and lhs.y >= rhs.y and lhs.z >= rhs.y
 */
template<class T>
bool operator > (gsVertex<T> const & lhs,gsVertex<T> const & rhs)
{
    return (lhs.x() < rhs.x() ||
        (lhs.x() == rhs.x() && lhs.y() < rhs.y()) ||
        (lhs.x() == rhs.x() && lhs.y() == rhs.y() && lhs.z() < rhs.z()));
}

/**
 * Compares to gsVertex
 * @param lhs LHS gsVertex
 * @param rhs RHS gsVertex
 * @return true if LHS != RHS, false if LHS == RHS
 */
template<class T>
bool operator != (gsVertex<T> const & lhs, gsVertex<T> const & rhs)
{
    //return lhs.Eigen::template Matrix<T,3,1>::operator!=(rhs);
    return !(lhs.x()== rhs.x()&& lhs.y()==rhs.y()&& lhs.z()==rhs.z());
}

} // namespace gismo
