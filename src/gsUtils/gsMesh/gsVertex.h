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

template <class T> 
class gsVertex  : public gsMeshElement<T>, public gsVector3d<T>
{
public:
    /// Shared pointer for gsVertex
    typedef memory::shared_ptr< gsVertex > Ptr;

    /// Unique pointer for gsVertex
    typedef memory::unique_ptr< gsVertex > uPtr;

    typedef typename gsVector3d<T>::Scalar       scalar_t;
    typedef gsMeshElement<T>                     MeshElement;
    typedef typename MeshElement::gsFaceHandle   gsFaceHandle;
    typedef typename MeshElement::gsVertexHandle gsVertexHandle;

public:
    gsVertex() : MeshElement(), gsVector3d<T>() { }

    gsVertex(scalar_t x, scalar_t y, scalar_t z = 0) : 
        MeshElement(), gsVector3d<T>(x,y,z),sharp(0)
    { }

    gsVertex( gsVector3d<T> const & u) : 
        MeshElement(), gsVector3d<T>(u),sharp(0)
    { }

    gsVertex( gsVector<T> const & u) : 
        MeshElement(), gsVector3d<T>(u),sharp(0)
    { 
        const index_t r = u.rows();

        if ( r < 3 )
        {
            this->bottomRows(3-r).setZero();
        }
    }

    virtual ~gsVertex() { };

    // clone function
    //GISMO_CLONE_FUNCTION(gsVertex)
    uPtr clone() const { return uPtr(new gsVertex(*this)); }

    void move(scalar_t dx, scalar_t dy, scalar_t dz) 
    {
        this->x() += dx;
        this->y() += dy;
        this->z() += dz;
    }
    
    inline void addFace(gsFaceHandle const& f)
    { faces.push_back( f ); }


    inline T   x () const { return this->operator()(0); }
    inline T   y () const { return this->operator()(1); }
    inline T   z () const { return this->operator()(2); }
    inline T & x () { return this->operator()(0); }
    inline T & y () { return this->operator()(1); }
    inline T & z () { return this->operator()(2); }

    std::ostream &print(std::ostream &os) const
    {
        os << "Vertex( " << this->x() << " " << this->y() << " " << this->z() << " )\n";
        return os;
    }

public:

    std::vector<gsVertexHandle> nVertices;
    //gsVector3d<T> coords;
    bool sharp;
    
    //List of faces adjacent to this vertex
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
template<class T>
bool operator == (gsVertex<T> const & lhs,gsVertex<T> const & rhs)
{
    return (lhs.x()==rhs.x())&&
           (lhs.y()==rhs.y())&&
           (lhs.z()==rhs.z());
//    return lhs.Eigen::template Matrix<T,3,1>::operator==(rhs); /slower
}
//void operator = (gsVertex<T> & lhs,gsVertex<T> const & rhs)
//{
//    lhs.coords=rhs.coords;
//    lhs.nVertices=rhs.nVertices;
//    lhs.numEdges=rhs.numEdges;
//    lhs.sharp=rhs.sharp;
//    lhs.faces=rhs.faces;
//    return 0;
//}

template<class T>
bool operator < (gsVertex<T> const & lhs,gsVertex<T> const & rhs)
{
    return (lhs.x()> rhs.x() || ( lhs.x()==rhs.x() && lhs.y()>rhs.y() )
            || ( lhs.x()==rhs.x() && lhs.y()==rhs.y() && lhs.z()>rhs.z()));
}
template<class T>
bool operator > (gsVertex<T> const & lhs,gsVertex<T> const & rhs)
{
    return (lhs.x()< rhs.x() || ( lhs.x()==rhs.x() && lhs.y()<rhs.y() )
            || ( lhs.x()==rhs.x() && lhs.y()==rhs.y() && lhs.z()<rhs.z()));
}

template<class T>
T operator *(gsVertex<T> const & lhs,gsVertex<T> const & rhs)
{
//    return (lhs.x()*rhs.x()+
//            lhs.y()*rhs.y()+
//            lhs.z()*rhs.z());
    return lhs.Eigen::template Matrix<T,3,1>::operator*(rhs);
}

template<class T>
gsVertex<T> operator -(gsVertex<T> const & lhs,gsVertex<T> const & rhs)
{
    //return gsVertex<T>(lhs.x()-rhs.x(),lhs.y()-rhs.y(),lhs.z()-rhs.z());
    return (gsVector3d<T>)lhs.Eigen::template Matrix<T,3,1>::operator-(rhs);
}

template<class T>
gsVertex<T> operator +(gsVertex<T> const & lhs,gsVertex<T> const & rhs)
{
    //return gsVertex<T>(lhs.x()+rhs.x(),lhs.y()+rhs.y(),lhs.z()+rhs.z());
        return (gsVector3d<T>)lhs.Eigen::template Matrix<T,3,1>::operator+(rhs);
}

template<class T>
bool operator != (gsVertex<T> const & lhs, gsVertex<T> const & rhs)
{
    //return lhs.Eigen::template Matrix<T,3,1>::operator!=(rhs);
    return !(lhs.x()== rhs.x()&& lhs.y()==rhs.y()&& lhs.z()==rhs.z());
}

} // namespace gismo

