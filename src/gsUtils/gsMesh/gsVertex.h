#pragma once

#include <gsUtils/gsMesh/gsMeshElement.h>
#include <gsCore/gsLinearAlgebra.h>
#include <set>
namespace gismo {

template <class T> 
class gsVertex  : public gsMeshElement<T> 
{
public:
    typedef typename gsVector3d<T>::Scalar       scalar_t;
    typedef gsMeshElement<T>                     MeshElement;
    typedef typename MeshElement::gsFaceHandle   gsFaceHandle;
    typedef typename MeshElement::gsVertexHandle gsVertexHandle;

public:
    gsVertex(scalar_t x, scalar_t y, scalar_t z = 0) : 
        MeshElement(), coords(x,y,z),sharp(0) 
    { }

    gsVertex( gsVector3d<T> const & u) : 
        MeshElement(),sharp(0) 
    { }

    gsVertex( gsVector<T> const & u) : 
        MeshElement(),sharp(0) 
    { 
        const index_t r = u.rows();

        if ( r == 3 )
            coords = u;
        else if ( r < 3 )
        {
            coords.topRows(r) = u;
            coords.bottomRows(3-r).setZero();
        }
        else
            coords = u.topRows(3);
    }

    virtual ~gsVertex() { };

    void move(scalar_t dx, scalar_t dy, scalar_t dz) 
    {
        coords.x() += dx;
        coords.y() += dy;
        coords.z() += dz;
    }
    
    inline void addFace(gsFaceHandle const& f)
    { faces.push_back( f ); }


    inline T   x () const { return coords(0); }
    inline T   y () const { return coords(1); }
    inline T   z () const { return coords(2); }
    inline T & x () { return coords(0); }
    inline T & y () { return coords(1); }
    inline T & z () { return coords(2); }

    std::ostream &print(std::ostream &os) const
    {
        os<<"Vertex( "<<coords.x()<<" "<<coords.y()<<" "<<coords.z()<<" )\n";
        return os;
    }

public:

    std::vector<gsVertexHandle> nVertices;
    gsVector3d<T> coords;
    bool sharp;
    
    //List of faces adjacent to this vertex
    std::vector<gsFaceHandle> faces;
    int numEdges;
    
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
    return (sqrt(vert.x()*vert.x()+vert.y()*vert.y()+vert.z()*vert.z()));
}

template<class T>
bool operator < (typename gsVertex<T>::gsVertexHandle const & lhs,
		typename gsVertex<T>::gsVertexHandle const & rhs)
{
  return !(lhs->x()< rhs->x() || ( lhs->x()==rhs->x() && lhs->y()<rhs->y() )
	   || ( lhs->x()==rhs->x() && lhs->y()==rhs->y() && lhs->z()<rhs->z()) 
    );}
template<class T>
bool operator > (typename gsVertex<T>::gsVertexHandle const & lhs,
        typename gsVertex<T>::gsVertexHandle const & rhs)
{
  return !(lhs->x()> rhs->x() || ( lhs->x()==rhs->x() && lhs->y()>rhs->y() )
       || ( lhs->x()==rhs->x() && lhs->y()==rhs->y() && lhs->z()>rhs->z())
    );}
//template<class T>
//bool operator == (typename gsVertex<T> const & lhs,
//        typename gsVertex<T> const & rhs)
//{
//  return (lhs->x()== rhs->x()&& lhs->y()==rhs->y()&& lhs->z()==rhs->z())
//    ;}
template<class T>
bool operator == (gsVertex<T> const & lhs,gsVertex<T> const & rhs)
{
    return (((lhs.x())==rhs.x())&&
            ((lhs.y())==rhs.y())&&
            ((lhs.z())==rhs.z()));
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
            || ( lhs.x()==rhs.x() && lhs.y()==rhs.y() && lhs.z()>rhs.z()))
         ;
}
template<class T>
bool operator > (gsVertex<T> const & lhs,gsVertex<T> const & rhs)
{
    return (lhs.x()< rhs.x() || ( lhs.x()==rhs.x() && lhs.y()<rhs.y() )
            || ( lhs.x()==rhs.x() && lhs.y()==rhs.y() && lhs.z()<rhs.z()))
         ;
}

template<class T>
T operator *(gsVertex<T> const & lhs,gsVertex<T> const & rhs)
{
    return (lhs.x()*rhs.x()+
            lhs.y()*rhs.y()+
            lhs.z()*rhs.z());
}

template<class T>
gsVertex<T> operator -(gsVertex<T> const & lhs,gsVertex<T> const & rhs)
{
    return gsVertex<T>(lhs.x()-rhs.x(),lhs.y()-rhs.y(),lhs.z()-rhs.z());
}
template<class T>
gsVertex<T> operator +(gsVertex<T> const & lhs,gsVertex<T> const & rhs)
{
    return gsVertex<T>(lhs.x()+rhs.x(),lhs.y()+rhs.y(),lhs.z()+rhs.z());
}
template<class T>
bool operator != (gsVertex<T> const & lhs, gsVertex<T> const & rhs)
{
  return !(lhs.x()== rhs.x()&& lhs.y()==rhs.y()&& lhs.z()==rhs.z())
    ;}

};// namespace gismo


