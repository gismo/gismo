/** @file gsMeshTopology.h

    @brief Combinatorial core shared by the surface and the volumetric meshes
    of gsMesh2: an oriented 2-map with handles, connectivity and properties.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

//#define Eigen gsEigen
//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(gismo::Point)
//#undef Eigen

#include <gsMesh2/gsProperty.h>


namespace gismo
{

/** \brief The combinatorial core of the gsMesh2 mesh classes: an oriented
    2-map with property arrays.

    gsMeshTopology carries everything that is common to the surface and to the
    volumetric mesh and that depends neither on the coordinate type nor on the
    dimension in which the mesh is read: the Vertex / Halfedge / Edge / Face
    handles and their connectivity, the property containers, the iterators and
    circulators, incremental construction through add_face(), deletion and
    garbage collection.

    Its four entity kinds are the cells of an oriented 2-map, and the index
    arithmetic they rely on -- \c edge(h)==h>>1 and \c opposite_halfedge(h)==h^1,
    i.e. exactly two halfedges per edge -- is what fixes their meaning:

    - in gsSurfMeshTopology they are read as the cells of a surface;
    - in gsVolMeshTopology the very same arithmetic makes them the \em cell-local
      cells of a volumetric mesh (a corner, an edge-use, a dart and a half-face),
      because a geometric edge of a volume mesh carries two halfedges \em per
      incident face and therefore cannot be a gsMeshTopology::Edge.

    Only operations that are meaningful under both readings live here; the Euler
    operators of a surface (collapse, flip, split, triangulate) belong to
    gsSurfMeshTopology.

    Being non-templated, all of this compiles exactly once into the library
    instead of once per scalar type.

    \sa gsSurfMeshTopology, gsVolMeshTopology
*/
class GISMO_EXPORT gsMeshTopology
{
public:
    using Self = gsMeshTopology;

public:

    /// Base class for all topology types (internally it is basically an index)
    /// \sa Vertex, Halfedge, Edge, Face
    class Base_handle
    {
    public:

        /// constructor
        explicit Base_handle(int _idx=-1) : idx_(_idx) {}

        /// Get the underlying index of this handle
        int idx() const { return idx_; }

        /// reset handle to be invalid (index=-1)
        void reset() { idx_=-1; }

        /// return whether the handle is valid, i.e., the index is not equal to -1.
        bool is_valid() const { return idx_ != -1; }

        /// are two handles equal?
        bool operator==(const Base_handle& _rhs) const {
            return idx_ == _rhs.idx_;
        }

        /// are two handles different?
        bool operator!=(const Base_handle& _rhs) const {
            return idx_ != _rhs.idx_;
        }

        /// compare operator useful for sorting handles
        bool operator<(const Base_handle& _rhs) const {
            return idx_ < _rhs.idx_;
        }

    private:
        friend class Vertex_iterator;
        friend class Halfedge_iterator;
        friend class Edge_iterator;
        friend class Face_iterator;
        friend class gsMeshTopology;
        int idx_;
    };


    /// this type represents a vertex (internally it is basically an index)
    ///  \sa Halfedge, Edge, Face
    struct GISMO_EXPORT Vertex : public Base_handle
    {
        /// default constructor (with invalid index)
        explicit Vertex(int _idx=-1) : Base_handle(_idx) {}
    };


    /// this type represents a halfedge (internally it is basically an index)
    /// \sa Vertex, Edge, Face
    struct GISMO_EXPORT Halfedge : public Base_handle
    {
        /// default constructor (with invalid index)
        explicit Halfedge(int _idx=-1) : Base_handle(_idx) {}
    };


    /// this type represents an edge (internally it is basically an index)
    /// \sa Vertex, Halfedge, Face
    struct Edge : public Base_handle
    {
        /// default constructor (with invalid index)
        explicit Edge(int _idx=-1) : Base_handle(_idx) {}
    };


    /// this type represents a face (internally it is basically an index)
    /// \sa Vertex, Halfedge, Edge
    struct Face : public Base_handle
    {
        /// default constructor (with invalid index)
        explicit Face(int _idx=-1) : Base_handle(_idx) {}
    };

public:

    /// This type stores the vertex connectivity
    /// \sa Halfedge_connectivity, Face_connectivity
    struct Vertex_connectivity
    {
        /// an outgoing halfedge per vertex (it will be a bounday halfedge for boundary vertices)
        Halfedge  halfedge_;
    };


    /// This type stores the halfedge connectivity
    /// \sa Vertex_connectivity, Face_connectivity
    struct Halfedge_connectivity
    {
        /// face incident to halfedge
        Face      face_;
        /// vertex the halfedge points to
        Vertex    vertex_;
        /// next halfedge within a face (or along a boundary)
        Halfedge  next_halfedge_;
        /// previous halfedge within a face (or along a boundary)
        Halfedge  prev_halfedge_;
    };


    /// This type stores the face connectivity
    /// \sa Vertex_connectivity, Halfedge_connectivity
    struct Face_connectivity
    {
        /// a halfedge that is part of the face
        Halfedge  halfedge_;

    };




public:

    /// Vertex property of type T
    /// \sa Halfedge_property, Edge_property, Face_property
    template <class T> class Vertex_property : public gsProperty<T>
    {
    public:

        /// default constructor
        explicit Vertex_property() {}
        explicit Vertex_property(gsProperty<T> p) : gsProperty<T>(p) {}

        /// access the data stored for vertex \c v
        typename gsProperty<T>::reference operator[](Vertex v)
        {
            return gsProperty<T>::operator[](v.idx());
        }

        /// access the data stored for vertex \c v
        typename gsProperty<T>::const_reference operator[](Vertex v) const
        {
            return gsProperty<T>::operator[](v.idx());
        }
    };


    /// Halfedge property of type T
    /// \sa Vertex_property, Edge_property, Face_property
    template <class T> class Halfedge_property : public gsProperty<T>
    {
    public:

        /// default constructor
        explicit Halfedge_property() {}
        explicit Halfedge_property(gsProperty<T> p) : gsProperty<T>(p) {}

        /// access the data stored for halfedge \c h
        typename gsProperty<T>::reference operator[](Halfedge h)
        {
            return gsProperty<T>::operator[](h.idx());
        }

        /// access the data stored for halfedge \c h
        typename gsProperty<T>::const_reference operator[](Halfedge h) const
        {
            return gsProperty<T>::operator[](h.idx());
        }
    };


    /// Edge property of type T
    /// \sa Vertex_property, Halfedge_property, Face_property
    template <class T> class Edge_property : public gsProperty<T>
    {
    public:

        /// default constructor
        explicit Edge_property() {}
        explicit Edge_property(gsProperty<T> p) : gsProperty<T>(p) {}

        /// access the data stored for edge \c e
        typename gsProperty<T>::reference operator[](Edge e)
        {
            return gsProperty<T>::operator[](e.idx());
        }

        /// access the data stored for edge \c e
        typename gsProperty<T>::const_reference operator[](Edge e) const
        {
            return gsProperty<T>::operator[](e.idx());
        }
    };


    /// Face property of type T
    /// \sa Vertex_property, Halfedge_property, Edge_property
    template <class T> class Face_property : public gsProperty<T>
    {
    public:

        /// default constructor
        explicit Face_property() {}
        explicit Face_property(gsProperty<T> p) : gsProperty<T>(p) {}

        /// access the data stored for face \c f
        typename gsProperty<T>::reference operator[](Face f)
        {
            return gsProperty<T>::operator[](f.idx());
        }

        /// access the data stored for face \c f
        typename gsProperty<T>::const_reference operator[](Face f) const
        {
            return gsProperty<T>::operator[](f.idx());
        }
    };


    /// Mesh property of type T
    /// \sa Vertex_property, Halfedge_property, Edge_property
    template <class T> class Mesh_property : public gsProperty<T>
    {
    public:

        /// default constructor
        explicit Mesh_property() {}
        explicit Mesh_property(gsProperty<T> p) : gsProperty<T>(p) {}

        /// access the data stored for the mesh
        typename gsProperty<T>::reference operator[](size_t idx)
        {
            return gsProperty<T>::operator[](idx);
        }

        /// access the data stored for the mesh
        typename gsProperty<T>::const_reference operator[](size_t idx) const
        {
            return gsProperty<T>::operator[](idx);
        }
    };



public:

    /// this class iterates linearly over all vertices
    /// \sa vertices_begin(), vertices_end()
    /// \sa Halfedge_iterator, Edge_iterator, Face_iterator
    class Vertex_iterator
    {
    public:

        /// Default constructor
        Vertex_iterator(Vertex v=Vertex(), const Self* m=NULL) : hnd_(v), mesh_(m)
        {
            if (mesh_ && mesh_->garbage()) while (mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
        }

        /// get the vertex the iterator refers to
        Vertex operator*()  const { return  hnd_; }

        const Vertex * operator->()  const { return  &hnd_; }

        // (!)
        size_t operator-(const Vertex_iterator& rhs) const
        {
            //assumes no deleted vertices..
            return (hnd_.idx_ - rhs.hnd_.idx_);
        }

        // (!)
        Vertex_iterator& operator+=(const size_t i)
        {
            //assumes no deleted vertices..
            hnd_.idx_+= i;
            return *this;
        }

        /// are two iterators equal?
        bool operator==(const Vertex_iterator& rhs) const
        {
            return (hnd_==rhs.hnd_);
        }

        /// how do two iterators compare?
        bool operator<(const Vertex_iterator& rhs) const
        {
            return (hnd_<rhs.hnd_);
        }

        /// are two iterators different?
        bool operator!=(const Vertex_iterator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment iterator
        Vertex_iterator& operator++()
        {
            ++hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
            return *this;
        }

        /// pre-decrement iterator
        Vertex_iterator& operator--()
        {
            --hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) --hnd_.idx_;
            return *this;
        }

    private:
        Vertex  hnd_;
        const Self* mesh_;
    };


    /// this class iterates linearly over all halfedges
    /// \sa halfedges_begin(), halfedges_end()
    /// \sa Vertex_iterator, Edge_iterator, Face_iterator
    class Halfedge_iterator
    {
    public:

        /// Default constructor
        Halfedge_iterator(Halfedge h=Halfedge(), const Self* m=NULL) : hnd_(h), mesh_(m)
        {
            if (mesh_ && mesh_->garbage()) while (mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
        }

        /// get the halfedge the iterator refers to
        Halfedge operator*()  const { return  hnd_; }

        const Halfedge * operator->()  const { return  &hnd_; }

        /// are two iterators equal?
        bool operator==(const Halfedge_iterator& rhs) const
        {
            return (hnd_==rhs.hnd_);
        }

        /// are two iterators different?
        bool operator!=(const Halfedge_iterator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment iterator
        Halfedge_iterator& operator++()
        {
            ++hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
            return *this;
        }

        /// pre-decrement iterator
        Halfedge_iterator& operator--()
        {
            --hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) --hnd_.idx_;
            return *this;
        }

    private:
        Halfedge  hnd_;
        const Self* mesh_;
    };


    /// this class iterates linearly over all edges
    /// \sa edges_begin(), edges_end()
    /// \sa Vertex_iterator, Halfedge_iterator, Face_iterator
    class Edge_iterator
    {
    public:

        /// Default constructor
        Edge_iterator(Edge e=Edge(), const Self* m=NULL) : hnd_(e), mesh_(m)
        {
            if (mesh_ && mesh_->garbage()) while (mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
        }

        /// get the edge the iterator refers to
        Edge operator*()  const { return  hnd_; }

        const Edge * operator->()  const { return  &hnd_; }

        /// are two iterators equal?
        bool operator==(const Edge_iterator& rhs) const
        {
            return (hnd_==rhs.hnd_);
        }

        /// are two iterators different?
        bool operator!=(const Edge_iterator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment iterator
        Edge_iterator& operator++()
        {
            ++hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
            return *this;
        }

        /// pre-decrement iterator
        Edge_iterator& operator--()
        {
            --hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) --hnd_.idx_;
            return *this;
        }

    private:
        Edge  hnd_;
        const Self* mesh_;
    };


    /// this class iterates linearly over all faces
    /// \sa faces_begin(), faces_end()
    /// \sa Vertex_iterator, Halfedge_iterator, Edge_iterator
    class Face_iterator
    {
    public:

        /// Default constructor
        Face_iterator(Face f=Face(), const Self* m=NULL) : hnd_(f), mesh_(m)
        {
            if (mesh_ && mesh_->garbage()) while (mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
        }

        /// get the face the iterator refers to
        Face operator*()  const { return  hnd_; }

        const Face * operator->()  const { return  &hnd_; }

        // (!)
        size_t operator-(const Face_iterator& rhs) const
        {
            //assumes no deleted faces..
            return (hnd_.idx_ - rhs.hnd_.idx_);
        }

        // (!)
        Face_iterator& operator+=(const size_t i)
        {
            //assumes no deleted faces..
            hnd_.idx_+= i;
            return *this;
        }

        /// are two iterators equal?
        bool operator==(const Face_iterator& rhs) const
        {
            return (hnd_==rhs.hnd_);
        }

        /// how do two iterators compare?
        bool operator<(const Face_iterator& rhs) const
        {
            return (hnd_<rhs.hnd_);
        }

        /// are two iterators different?
        bool operator!=(const Face_iterator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment iterator
        Face_iterator& operator++()
        {
            ++hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
            return *this;
        }

        /// pre-decrement iterator
        Face_iterator& operator--()
        {
            --hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) --hnd_.idx_;
            return *this;
        }

    private:
        Face  hnd_;
        const Self* mesh_;
    };



public:

    /// this helper class is a container for iterating through all
    /// vertices using C++11 range-based for-loops.
    /// \sa vertices()
    class Vertex_container
    {
    public:
        Vertex_container(Vertex_iterator _begin, Vertex_iterator _end) : begin_(_begin), end_(_end) {}
        Vertex_iterator begin() const { return begin_; }
        Vertex_iterator end()   const { return end_;   }
    private:
        Vertex_iterator begin_, end_;
    };



    /// this helper class is a container for iterating through all
    /// halfedge using C++11 range-based for-loops.
    /// \sa halfedges()
    class Halfedge_container
    {
    public:
        Halfedge_container(Halfedge_iterator _begin, Halfedge_iterator _end) : begin_(_begin), end_(_end) {}
        Halfedge_iterator begin() const { return begin_; }
        Halfedge_iterator end()   const { return end_;   }
    private:
        Halfedge_iterator begin_, end_;
    };



    /// this helper class is a container for iterating through all
    /// edges using C++11 range-based for-loops.
    /// \sa edges()
    class Edge_container
    {
    public:
        Edge_container(Edge_iterator _begin, Edge_iterator _end) : begin_(_begin), end_(_end) {}
        Edge_iterator begin() const { return begin_; }
        Edge_iterator end()   const { return end_;   }
    private:
        Edge_iterator begin_, end_;
    };



    /// this helper class is a container for iterating through all
    /// faces using C++11 range-based for-loops.
    /// \sa faces()
    class Face_container
    {
    public:
        Face_container(Face_iterator _begin, Face_iterator _end) : begin_(_begin), end_(_end) {}
        Face_iterator begin() const { return begin_; }
        Face_iterator end()   const { return end_;   }
    private:
        Face_iterator begin_, end_;
    };





public:

    /// this class circulates through all one-ring neighbors of a vertex.
    /// it also acts as a container-concept for C++11 range-based for loops.
    /// \sa Halfedge_around_vertex_circulator, Face_around_vertex_circulator, vertices(Vertex)
    class Vertex_around_vertex_circulator
    {
    public:

        /// default constructor
        Vertex_around_vertex_circulator(const Self* m=NULL, Vertex v=Vertex())
        : mesh_(m), active_(true)
        {
            if (mesh_)
                halfedge_ = mesh_->halfedge(v);
        }

        /// are two circulators equal?
        bool operator==(const Vertex_around_vertex_circulator& rhs) const
        {
            assert(mesh_);
            return (active_ && (mesh_==rhs.mesh_) && (halfedge_==rhs.halfedge_));
        }

        /// are two circulators different?
        bool operator!=(const Vertex_around_vertex_circulator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment (rotate counter-clockwise)
        Vertex_around_vertex_circulator& operator++()
        {
            assert(mesh_);
            halfedge_ = mesh_->ccw_rotated_halfedge(halfedge_);
            active_ = true;
            return *this;
        }

        /// pre-decrement (rotate clockwise)
        Vertex_around_vertex_circulator& operator--()
        {
            assert(mesh_);
            halfedge_ = mesh_->cw_rotated_halfedge(halfedge_);
            return *this;
        }

        /// get the vertex the circulator refers to
        Vertex operator*()  const
        {
            assert(mesh_);
            return mesh_->to_vertex(halfedge_);
        }

        /// cast to bool: true if vertex is not isolated
        operator bool() const { return halfedge_.is_valid(); }

        /// return current halfedge
        Halfedge halfedge() const { return halfedge_; }

        // helper for C++11 range-based for-loops
        Vertex_around_vertex_circulator& begin()
        {
            active_=!halfedge_.is_valid();
            return *this;
        }

        // helper for C++11 range-based for-loops
        Vertex_around_vertex_circulator& end()   { active_=true;  return *this; }

    private:
        const Self*  mesh_;
        Halfedge         halfedge_;
        // helper for C++11 range-based for-loops
        bool active_;
    };


    /// this class circulates through all outgoing halfedges of a vertex.
    /// it also acts as a container-concept for C++11 range-based for loops.
    /// \sa Vertex_around_vertex_circulator, Face_around_vertex_circulator, halfedges(Vertex)
    class Halfedge_around_vertex_circulator
    {
    public:

        /// default constructor
        Halfedge_around_vertex_circulator(const Self* m=NULL, Vertex v=Vertex())
        : mesh_(m), active_(true)
        {
            if (mesh_) halfedge_ = mesh_->halfedge(v);
        }

        /// are two circulators equal?
        bool operator==(const Halfedge_around_vertex_circulator& rhs) const
        {
            assert(mesh_);
            return (active_ && (mesh_==rhs.mesh_) && (halfedge_==rhs.halfedge_));
        }

        /// are two circulators different?
        bool operator!=(const Halfedge_around_vertex_circulator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment (rotate couter-clockwise)
        Halfedge_around_vertex_circulator& operator++()
        {
            assert(mesh_);
            halfedge_ = mesh_->ccw_rotated_halfedge(halfedge_);
            active_ = true;
            return *this;
        }

        /// pre-decrement (rotate clockwise)
        Halfedge_around_vertex_circulator& operator--()
        {
            assert(mesh_);
            halfedge_ = mesh_->cw_rotated_halfedge(halfedge_);
            return *this;
        }

        /// get the halfedge the circulator refers to
        Halfedge operator*() const { return halfedge_; }

        const Halfedge * operator->()  const { return  &halfedge_; }

        /// cast to bool: true if vertex is not isolated
        operator bool() const { return halfedge_.is_valid(); }

        // helper for C++11 range-based for-loops
        Halfedge_around_vertex_circulator& begin() { active_=!halfedge_.is_valid(); return *this; }
        // helper for C++11 range-based for-loops
        Halfedge_around_vertex_circulator& end()   { active_=true;  return *this; }

    private:
        const Self*  mesh_;
        Halfedge         halfedge_;
        // helper for C++11 range-based for-loops
        bool active_;
    };


    /// this class circulates through all incident faces of a vertex.
    /// it also acts as a container-concept for C++11 range-based for loops.
    /// \sa Vertex_around_vertex_circulator, Halfedge_around_vertex_circulator, faces(Vertex)
    class Face_around_vertex_circulator
    {
    public:

        /// construct with mesh and vertex (vertex should not be isolated!)
        Face_around_vertex_circulator(const Self* m=NULL, Vertex v=Vertex())
        : mesh_(m), active_(true)
        {
            if (mesh_)
            {
                halfedge_ = mesh_->halfedge(v);
                if (halfedge_.is_valid() && mesh_->is_boundary(halfedge_))
                    operator++();
            }
        }

        /// are two circulators equal?
        bool operator==(const Face_around_vertex_circulator& rhs) const
        {
            assert(mesh_);
            return (active_ && (mesh_==rhs.mesh_) && (halfedge_==rhs.halfedge_));
        }

        /// are two circulators different?
        bool operator!=(const Face_around_vertex_circulator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment (rotates counter-clockwise)
        Face_around_vertex_circulator& operator++()
        {
            assert(mesh_ && halfedge_.is_valid());
            do {
                if (halfedge_ == mesh_->ccw_rotated_halfedge(halfedge_))
                    break;
                halfedge_ = mesh_->ccw_rotated_halfedge(halfedge_);
            } while (mesh_->is_boundary(halfedge_));
            active_ = true;
            return *this;
        }

        /// pre-decrement (rotate clockwise)
        Face_around_vertex_circulator& operator--()
        {
            assert(mesh_ && halfedge_.is_valid());
            do
                halfedge_ = mesh_->cw_rotated_halfedge(halfedge_);
            while (mesh_->is_boundary(halfedge_));
            return *this;
        }

        /// get the face the circulator refers to
        Face operator*() const
        {
            assert(mesh_ && halfedge_.is_valid());
            return mesh_->face(halfedge_);
        }

        /// cast to bool: true if vertex is not isolated
        operator bool() const { return halfedge_.is_valid(); }

        // helper for C++11 range-based for-loops
        Face_around_vertex_circulator& begin()
        {
            active_=!(halfedge_.is_valid() );
            return *this;
        }
        // helper for C++11 range-based for-loops
        Face_around_vertex_circulator& end()   { active_=true;  return *this; }

    private:
        const Self*  mesh_;
        Halfedge         halfedge_;
        // helper for C++11 range-based for-loops
        bool active_;
    };


    /// this class circulates through the vertices of a face.
    /// it also acts as a container-concept for C++11 range-based for loops.
    /// \sa Halfedge_around_face_circulator, vertices(Face)
    class Vertex_around_face_circulator
    {
    public:

        /// default constructor
        Vertex_around_face_circulator(const Self* m=NULL, Face f=Face())
        : mesh_(m), active_(true)
        {
            if (mesh_) halfedge_ = mesh_->halfedge(f);
        }

        /// are two circulators equal?
        bool operator==(const Vertex_around_face_circulator& rhs) const
        {
            assert(mesh_);
            return (active_ && (mesh_==rhs.mesh_) && (halfedge_==rhs.halfedge_));
        }

        /// are two circulators different?
        bool operator!=(const Vertex_around_face_circulator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment (rotates counter-clockwise)
        Vertex_around_face_circulator& operator++()
        {
            assert(mesh_ && halfedge_.is_valid());
            halfedge_ = mesh_->next_halfedge(halfedge_);
            active_ = true;
            return *this;
        }

        /// pre-decrement (rotates clockwise)
        Vertex_around_face_circulator& operator--()
        {
            assert(mesh_ && halfedge_.is_valid());
            halfedge_ = mesh_->prev_halfedge(halfedge_);
            return *this;
        }

        /// get the vertex the circulator refers to
        Vertex operator*() const
        {
            assert(mesh_ && halfedge_.is_valid());
            return mesh_->to_vertex(halfedge_);
        }

        /// get the halfedge the circulator refers to
        Halfedge he() const
        {
            assert(mesh_ && halfedge_.is_valid());
            return halfedge_;
        }

        // helper for C++11 range-based for-loops
        Vertex_around_face_circulator& begin() { active_=false; return *this; }
        // helper for C++11 range-based for-loops
        Vertex_around_face_circulator& end()   { active_=true;  return *this; }

    private:
        const Self*  mesh_;
        Halfedge         halfedge_;
        // helper for C++11 range-based for-loops
        bool active_;
    };


    /// this class circulates through all halfedges of a face.
    /// it also acts as a container-concept for C++11 range-based for loops.
    /// \sa Vertex_around_face_circulator, halfedges(Face)
    class Halfedge_around_face_circulator
    {
    public:

        /// default constructur
        Halfedge_around_face_circulator(const Self* m=NULL, Face f=Face())
        : mesh_(m), active_(true)
        {
            if (mesh_) halfedge_ = mesh_->halfedge(f);
        }

        /// are two circulators equal?
        bool operator==(const Halfedge_around_face_circulator& rhs) const
        {
            assert(mesh_);
            return (active_ && (mesh_==rhs.mesh_) && (halfedge_==rhs.halfedge_));
        }

        /// are two circulators different?
        bool operator!=(const Halfedge_around_face_circulator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment (rotates counter-clockwise)
        Halfedge_around_face_circulator& operator++()
        {
            assert(mesh_ && halfedge_.is_valid());
            halfedge_ = mesh_->next_halfedge(halfedge_);
            active_ = true;
            return *this;
        }

        /// pre-decrement (rotates clockwise)
        Halfedge_around_face_circulator& operator--()
        {
            assert(mesh_ && halfedge_.is_valid());
            halfedge_ = mesh_->prev_halfedge(halfedge_);
            return *this;
        }

        /// get the halfedge the circulator refers to
        Halfedge operator*() const { return halfedge_; }

        const Halfedge * operator->()  const { return  &halfedge_; }

        // helper for C++11 range-based for-loops
        Halfedge_around_face_circulator& begin() { active_=false; return *this; }
        // helper for C++11 range-based for-loops
        Halfedge_around_face_circulator& end()   { active_=true;  return *this; }

    private:
        const Self*  mesh_;
        Halfedge         halfedge_;
        // helper for C++11 range-based for-loops
        bool active_;
    };

public:

    /// \name Construct, destruct, assignment
    //@{

    /// default constructor
    gsMeshTopology();

    // destructor (is virtual, since we inherit from Geometry_representation)
    virtual ~gsMeshTopology();

    /// copy constructor: copies \c rhs to \c *this. performs a deep copy of all properties.
    gsMeshTopology(const gsMeshTopology& rhs) { operator=(rhs); }

    /// assign \c rhs to \c *this. performs a deep copy of all properties.
    gsMeshTopology& operator=(const gsMeshTopology& rhs);

    /// move \c rhs to \c *this
    gsMeshTopology& operator=(gsMeshTopology&& rhs) noexcept;

    /// assign \c rhs to \c *this. does not copy custom properties.
    gsMeshTopology& assign(const gsMeshTopology& rhs);

    //@}

public:

    /// add \a nverts new vertex at once, with position \a 0
    void add_batch_vertices(size_t nverts);

    /// add a new face with vertex list \c vertices
    /// \sa add_triangle, add_quad
    Face add_face(const std::vector<Vertex>& vertices);

    /// add a new triangle connecting vertices \c v1, \c v2, \c v3
    /// \sa add_face, add_quad
    Face add_triangle(Vertex v1, Vertex v2, Vertex v3);

    /// add a new quad connecting vertices \c v1, \c v2, \c v3, \c v4
    /// \sa add_triangle, add_face
    Face add_quad(Vertex v1, Vertex v2, Vertex v3, Vertex v4);

    // Insert a new quad by four consecutive border edges
    Face add_quad(Edge e0, Edge e1, Edge e2, Edge e3);

    // Insert a new edge between two vertices, even if one exists
    Edge add_edge(Vertex start, Vertex end)
    {
        if (end.idx()>start.idx())
            std::swap(start,end);

        Halfedge h = new_edge(start,end);
        set_halfedge(start, h);
        set_halfedge(end, opposite_halfedge(h));
        set_next_halfedge(h, opposite_halfedge(h));
        set_next_halfedge(opposite_halfedge(h), h);
        return edge(h);
    }

    /// Insert new edge between two vertices, if one does not exist
    /// Complexity proportional to the number of halfedges
    Edge find_or_add_edge(Vertex start, Vertex end)
    {
        if (end.idx()>start.idx())
            std::swap(start,end);

        //brut search
        for( Halfedge a : halfedges() )
            if ( from_vertex(a) == start && to_vertex(a)== end)
                return edge(a);
        return add_edge(start,end);

        /*
        Halfedge h = find_halfedge(start,end);
        if (!h.is_valid())
        {   //try the opposite as well
            h = find_halfedge(end,start);
            if ( h.is_valid() )
                h = opposite_halfedge(h);
        }
        return h.is_valid() ? edge(h) : add_edge(start,end);
        */
    }

public:

    /// returns number of (deleted and valid) vertices in the mesh
    unsigned int vertices_size() const { return (unsigned int) vprops_.size(); }
    /// returns number of (deleted and valid)halfedge in the mesh
    unsigned int halfedges_size() const { return (unsigned int) hprops_.size(); }
    /// returns number of (deleted and valid)edges in the mesh
    unsigned int edges_size() const { return (unsigned int) eprops_.size(); }
    /// returns number of (deleted and valid)faces in the mesh
    unsigned int faces_size() const { return (unsigned int) fprops_.size(); }

    /// returns number of vertices in the mesh
    unsigned int n_vertices() const { return vertices_size() - deleted_vertices_; }
    /// returns number of halfedge in the mesh
    unsigned int n_halfedges() const { return halfedges_size() - 2*deleted_edges_; }
    /// returns number of edges in the mesh
    unsigned int n_edges() const { return edges_size() - deleted_edges_; }
    /// returns number of faces in the mesh
    unsigned int n_faces() const { return faces_size() - deleted_faces_; }

    /// returns true iff the mesh is empty, i.e., has no vertices
    unsigned int empty() const { return n_vertices() == 0; }

    /// clear mesh: remove all vertices, edges, faces
    void clear();

    /// remove unused memory from vectors
    void free_memory();

    /// reserve memory (mainly used in file readers)
    void reserve(unsigned int nvertices,
                 unsigned int nedges,
                 unsigned int nfaces );


    /// remove deleted vertices/edges/faces
    virtual void garbage_collection();

    /// returns whether vertex \c v is deleted
    /// \sa garbage_collection()
    bool is_deleted(Vertex v) const {return vdeleted_[v];}

    /// returns whether halfedge \c h is deleted
    /// \sa garbage_collection()
    bool is_deleted(Halfedge h) const {return edeleted_[edge(h)];}

    /// returns whether edge \c e is deleted
    /// \sa garbage_collection()
    bool is_deleted(Edge e) const { return edeleted_[e];}

    /// returns whether face \c f is deleted
    /// \sa garbage_collection()
    bool is_deleted(Face f) const {return fdeleted_[f];}

public:

    /// returns an outgoing halfedge of vertex \c v.
    /// if \c v is a boundary vertex this will be a boundary halfedge.
    Halfedge halfedge(Vertex v) const { return vconn_[v].halfedge_; }

    /// returns the \c i'th halfedge of edge \c e. \c i has to be 0 or 1.
    Halfedge halfedge(Edge e, unsigned int i) const
    {
        assert(i<=1);
        return Halfedge((e.idx() << 1) + i);
    }

    /// returns a halfedge of face \c f
    Halfedge halfedge(Face f) const { return fconn_[f].halfedge_; }

    /// return the edge that contains halfedge \c h as one of its two halfedges.
    inline Edge edge(Halfedge h) const { return Edge(h.idx() >> 1); }

    /// returns the \c i'th vertex of edge \c e. \c i has to be 0 or 1.
    Vertex vertex(Edge e, unsigned int i) const
    {
        assert(i<=1);
        return to_vertex(halfedge(e, i));
    }

    /// returns the vertex the halfedge \c h emanates from
    Vertex from_vertex(Halfedge h) const { return to_vertex(opposite_halfedge(h)); }

    /// returns the vertex the halfedge \c h points to
    Vertex to_vertex(Halfedge h) const { return hconn_[h].vertex_; }

    /// returns the face incident to halfedge \c h
    Face face(Halfedge h) const
    {
        return hconn_[h].face_;
    }

    /// returns the face incident to the \c i'th halfedge of edge \c e. \c i has to be 0 or 1.
    Face face(Edge e, unsigned int i) const
    {
        assert(i<=1);
        return face(halfedge(e, i));
    }

public:

    /// set the outgoing halfedge of vertex \c v to \c h
    void set_halfedge(Vertex v, Halfedge h) { vconn_[v].halfedge_ = h; }

    /// sets the vertex the halfedge \c h points to to \c v
    void set_vertex(Halfedge h, Vertex v) { hconn_[h].vertex_ = v; }

    /// sets the incident face to halfedge \c h to \c f
    void set_face(Halfedge h, Face f) { hconn_[h].face_ = f; }

    /// sets the next halfedge of \c h within the face to \c nh
    void set_next_halfedge(Halfedge h, Halfedge nh)
    {
        hconn_[h].next_halfedge_ = nh;
        hconn_[nh].prev_halfedge_ = h;
    }

    /// sets the halfedge of face \c f to \c h
    void set_halfedge(Face f, Halfedge h) { fconn_[f].halfedge_ = h; }

public:

    /// returns the next halfedge within the incident face
    Halfedge next_halfedge(Halfedge h) const { return hconn_[h].next_halfedge_; }

    /// returns the previous halfedge within the incident face
    Halfedge prev_halfedge(Halfedge h) const  { return hconn_[h].prev_halfedge_; }

    /// returns the opposite halfedge of \c h
    Halfedge opposite_halfedge(Halfedge h) const
    { return Halfedge((h.idx() & 1) ? h.idx()-1 : h.idx()+1); }

    /// returns the halfedge that is rotated counter-clockwise around the
    /// start vertex of \c h. it is the opposite halfedge of the previous halfedge of \c h.
    inline Halfedge ccw_rotated_halfedge(Halfedge h) const
    { return opposite_halfedge(prev_halfedge(h)); }

    /// returns the halfedge that is rotated clockwise around the
    /// start vertex of \c h. it is the next halfedge of the opposite halfedge of \c h.
    inline Halfedge cw_rotated_halfedge(Halfedge h) const
    { return next_halfedge(opposite_halfedge(h)); }

    inline Halfedge forward_halfedge(Self::Halfedge he) const
    { return next_halfedge(opposite_halfedge(next_halfedge(he))); }

    inline Halfedge backward_halfedge(Self::Halfedge he) const
    { return prev_halfedge(opposite_halfedge(prev_halfedge(he))); }

    inline Halfedge forward_halfedge(gsMeshTopology::Halfedge he, int r) const
    {
        for (int k = 0; k < r; k++)
            he = forward_halfedge(he);
        return he;
    }

    inline Halfedge backward_halfedge(gsMeshTopology::Halfedge he, int r) const
    {
        for (int k = 0; k < r; k++)
            he = backward_halfedge(he);
        return he;
    }

    inline Halfedge next_halfedge(gsMeshTopology::Halfedge he, int r) const
    {
        for (int k = 0; k < r; k++)
            he = next_halfedge(he);
        return he;
    }

    inline Halfedge prev_halfedge(gsMeshTopology::Halfedge he, int r) const
    {
        for (int k = 0; k < r; k++)
            he = prev_halfedge(he);
        return he;
    }


public:

    /// return whether vertex \c v is valid, i.e. the index it stores is within the array bounds.
    bool is_valid(Vertex v) const
    { return (0 <= v.idx()) && (v.idx() < (int)vertices_size()); }
    /// return whether halfedge \c h is valid, i.e. the index it stores is within the array bounds.
    bool is_valid(Halfedge h) const { return (0 <= h.idx()) && (h.idx() < (int)halfedges_size());}
    /// return whether edge \c e is valid, i.e. the index it stores is within the array bounds.
    bool is_valid(Edge e) const { return (0 <= e.idx()) && (e.idx() < (int)edges_size()); }
    /// return whether face \c f is valid, i.e. the index it stores is within the array bounds.
    bool is_valid(Face f) const { return (0 <= f.idx()) && (f.idx() < (int)faces_size()); }

    /// returns whether \c v is a manifold vertex (not incident to several patches)
    bool is_manifold(Vertex v) const
    {
        // The vertex is non-manifold if more than one gap exists, i.e.
        // more than one outgoing boundary halfedge.
        int n(0);
        Halfedge_around_vertex_circulator hit = halfedges(v), hend=hit;
        if (hit) do
        {
            if (is_boundary(*hit))
                ++n;
        }
        while (++hit!=hend);
        return n<2;
    }

    /// returns whether \c v is a boundary vertex
    bool is_boundary(Vertex v) const
    {
        Halfedge h(halfedge(v));
        return (!(h.is_valid() && face(h).is_valid()));
    }

    /// returns whether \c e is a boundary edge, i.e., if one of its
    /// halfedges is a boundary halfedge.
    bool is_boundary(Edge e) const
    { return (is_boundary(halfedge(e, 0)) || is_boundary(halfedge(e, 1))); }

    /// returns whether \c f is a boundary face, i.e., it one of its edges is a boundary edge.
    bool is_boundary(Face f) const
    {
        Halfedge h  = halfedge(f);
        Halfedge hh = h;
        do
        {
            if (is_boundary(opposite_halfedge(h)))
                return true;
            h = next_halfedge(h);
        }
        while (h != hh);
        return false;
    }

    /// returns whether h is a boundary halfege, i.e., if its face does not exist.
    bool is_boundary(Halfedge h) const { return !face(h).is_valid(); }

    /// returns whether h touches the boundary, i.e., its opposite is a boundary halfedge.
    bool touches_boundary(Halfedge h) const { return is_boundary(opposite_halfedge(h)); }

    /// returns whether \c v is isolated, i.e., not incident to any face
    bool is_isolated(Vertex v) const { return !halfedge(v).is_valid(); }

    inline bool is_near_boundary(Vertex v) const
    {
        for (auto u : vertices(v))
            if ( is_boundary(u) )
                return true;
        return false;
    }

public:

    /** add a vertex property of type \c T with name \c name and default value \c t.
     fails if a property named \c name exists already, since the name has to be unique.
     in this case it returns an invalid property */
    template <class T> Vertex_property<T> add_vertex_property(const std::string& name, T t=T())
    { return Vertex_property<T>(vprops_.add<T>(name, give(t))); }
    /** add a halfedge property of type \c T with name \c name and default value \c t.
     fails if a property named \c name exists already, since the name has to be unique.
     in this case it returns an invalid property */
    template <class T> Halfedge_property<T> add_halfedge_property(const std::string& name, T t=T())
    { return Halfedge_property<T>(hprops_.add<T>(name, give(t))); }
    /** add a edge property of type \c T with name \c name and default value \c t.
     fails if a property named \c name exists already, since the name has to be unique.
     in this case it returns an invalid property */
    template <class T> Edge_property<T> add_edge_property(const std::string& name, T t=T())
    {
        return Edge_property<T>(eprops_.add<T>(name, t));
    }
    /** add a face property of type \c T with name \c name and default value \c t.
     fails if a property named \c name exists already, since the name has to be unique.
     in this case it returns an invalid property */
    template <class T> Face_property<T> add_face_property(const std::string& name, const T t=T())
    { return Face_property<T>(fprops_.add<T>(name, t)); }
    /** add a mesh property of type \c T with name \c name and default value \c t.
     fails if a property named \c name exists already, since the name has to be unique.
     in this case it returns an invalid property */
    template <class T> Mesh_property<T> add_mesh_property(const std::string& name, const T t=T())
    { return Mesh_property<T>(mprops_.add<T>(name, t)); }

    /** get the vertex property named \c name of type \c T. returns an invalid
     Vertex_property if the property does not exist or if the type does not match. */
    template <class T> Vertex_property<T> get_vertex_property(const std::string& name) const
    { return Vertex_property<T>(vprops_.get<T>(name)); }
    /** get the halfedge property named \c name of type \c T. returns an invalid
     Vertex_property if the property does not exist or if the type does not match. */
    template <class T> Halfedge_property<T> get_halfedge_property(const std::string& name) const
    { return Halfedge_property<T>(hprops_.get<T>(name)); }
    /** get the edge property named \c name of type \c T. returns an invalid
     Vertex_property if the property does not exist or if the type does not match. */
    template <class T> Edge_property<T> get_edge_property(const std::string& name) const
    { return Edge_property<T>(eprops_.get<T>(name)); }
    /** get the face property named \c name of type \c T. returns an invalid
     Vertex_property if the property does not exist or if the type does not match. */
    template <class T> Face_property<T> get_face_property(const std::string& name) const
    { return Face_property<T>(fprops_.get<T>(name)); }
    /** get the mesh property named \c name of type \c T. returns an invalid
     Mesh_property if the property does not exist or if the type does not match. */
    template <class T> Mesh_property<T> get_mesh_property(const std::string& name) const
    { return Mesh_property<T>(mprops_.get<T>(name)); }

    /** if a vertex property of type \c T with name \c name exists, it is returned.
     otherwise this property is added (with default value \c t) */
    template <class T> Vertex_property<T> vertex_property(const std::string& name, const T t=T())
    { return Vertex_property<T>(vprops_.get_or_add<T>(name, give(t))); }
    /** if a halfedge property of type \c T with name \c name exists, it is returned.
     otherwise this property is added (with default value \c t) */
    template <class T> Halfedge_property<T> halfedge_property(const std::string& name, T t=T())
    { return Halfedge_property<T>(hprops_.get_or_add<T>(name, give(t))); }
    /** if an edge property of type \c T with name \c name exists, it is returned.
     otherwise this property is added (with default value \c t) */
    template <class T> Edge_property<T> edge_property(const std::string& name, const T t=T())
    { return Edge_property<T>(eprops_.get_or_add<T>(name, t)); }
    /** if a face property of type \c T with name \c name exists, it is returned.
     otherwise this property is added (with default value \c t) */
    template <class T> Face_property<T> face_property(const std::string& name, const T t=T())
    { return Face_property<T>(fprops_.get_or_add<T>(name, t)); }

    /** if a mesh property of type \c T with name \c name exists, it is returned.
     otherwise this property is added (with default value \c t) */
    template <class T> Mesh_property<T> mesh_property(const std::string& name, const T t=T())
    { return Mesh_property<T>(mprops_.get_or_add<T>(name, t)); }

    /// rename the vertex property \c p
    template <class T> void rename_vertex_property(Vertex_property<T>& p,
                                                   std::string newname)
    { vprops_.rename(p, give(newname)); }

    /// swaps (the values of) two vertex properties of the same type
    void swap_vertex_property(const std::string & name1,
                              const std::string & name2)
    { vprops_.swap(name1,name2); }

    /// remove the vertex property \c p
    template <class T> void remove_vertex_property(Vertex_property<T>& p)
    { vprops_.remove(p); }
    /// remove the halfedge property \c p
    template <class T> void remove_halfedge_property(Halfedge_property<T>& p)
    { hprops_.remove(p); }
    /// remove the edge property \c p
    template <class T> void remove_edge_property(Edge_property<T>& p)
    { eprops_.remove(p); }
    /// remove the face property \c p
    template <class T> void remove_face_property(Face_property<T>& p)
    { fprops_.remove(p); }
    /// remove the mesh property \c p
    template <class T> void remove_mesh_property(Mesh_property<T>& p)
    { mprops_.remove(p); }

    /** get the type_info \c T of vertex property named \c. returns an typeid(void)
     if the property does not exist or if the type does not match. */
    const std::type_info& get_vertex_property_type(const std::string& name) const
    { return vprops_.get_type(name); }
    /** get the type_info \c T of halfedge property named \c. returns an typeid(void)
     if the property does not exist or if the type does not match. */
    const std::type_info& get_halfedge_property_type(const std::string& name) const
    { return hprops_.get_type(name); }
    /** get the type_info \c T of edge property named \c. returns an typeid(void)
     if the property does not exist or if the type does not match. */
    const std::type_info& get_edge_property_type(const std::string& name) const
    { return eprops_.get_type(name); }
    /** get the type_info \c T of face property named \c. returns an typeid(void)
     if the property does not exist or if the type does not match. */
    const std::type_info& get_face_property_type(const std::string& name) const
    { return fprops_.get_type(name); }
    /** get the type_info \c T of face property named \c. returns an typeid(void)
     if the property does not exist or if the type does not match. */
    const std::type_info& get_mesh_property_type(const std::string& name) const
    { return mprops_.get_type(name); }

    /// returns the names of all vertex properties
    std::vector<std::string> vertex_properties() const
    { return vprops_.properties(); }
    /// returns the names of all halfedge properties
    std::vector<std::string> halfedge_properties() const
    { return hprops_.properties(); }
    /// returns the names of all edge properties
    std::vector<std::string> edge_properties() const
    { return eprops_.properties(); }
    /// returns the names of all face properties
    std::vector<std::string> face_properties() const
    { return fprops_.properties(); }
    /// returns the names of all mesh properties
    std::vector<std::string> mesh_properties() const
    { return mprops_.properties(); }

    /// prints the names of all properties
    void property_stats() const;

public:
    /// returns start iterator for vertices
    Vertex_iterator vertices_begin() const
    { return Vertex_iterator(Vertex(0), this); }

    /// returns end iterator for vertices
    Vertex_iterator vertices_end() const
    { return Vertex_iterator(Vertex(vertices_size()), this); }

    /// returns vertex container for C++11 range-based for-loops
    Vertex_container vertices() const
    { return Vertex_container(vertices_begin(), vertices_end()); }

    /// returns start iterator for halfedges
    Halfedge_iterator halfedges_begin() const
    { return Halfedge_iterator(Halfedge(0), this); }

    /// returns end iterator for halfedges
    Halfedge_iterator halfedges_end() const
    { return Halfedge_iterator(Halfedge(halfedges_size()), this); }

    /// returns halfedge container for C++11 range-based for-loops
    Halfedge_container halfedges() const
    { return Halfedge_container(halfedges_begin(), halfedges_end()); }

    /// returns start iterator for edges
    Edge_iterator edges_begin() const
    { return Edge_iterator(Edge(0), this); }

    /// returns end iterator for edges
    Edge_iterator edges_end() const
    { return Edge_iterator(Edge(edges_size()), this); }

    /// returns edge container for C++11 range-based for-loops
    Edge_container edges() const
    { return Edge_container(edges_begin(), edges_end()); }

    /// returns start iterator for faces
    Face_iterator faces_begin() const
    { return Face_iterator(Face(0), this); }

    /// returns end iterator for faces
    Face_iterator faces_end() const
    { return Face_iterator(Face(faces_size()), this); }

    /// returns face container for C++11 range-based for-loops
    Face_container faces() const
    { return Face_container(faces_begin(), faces_end()); }

    /// returns circulator for vertices around vertex \c v
    Vertex_around_vertex_circulator vertices(Vertex v) const
    { return Vertex_around_vertex_circulator(this, v); }

    /// returns circulator for outgoing halfedges around vertex \c v
    Halfedge_around_vertex_circulator halfedges(Vertex v) const
    { return Halfedge_around_vertex_circulator(this, v); }

    /// returns circulator for faces around vertex \c v
    Face_around_vertex_circulator faces(Vertex v) const
    { return Face_around_vertex_circulator(this, v); }

    /// returns circulator for vertices of face \c f
    Vertex_around_face_circulator vertices(Face f) const
    { return Vertex_around_face_circulator(this, f); }

    /// returns circulator for halfedges of face \c f
    Halfedge_around_face_circulator halfedges(Face f) const
    { return Halfedge_around_face_circulator(this, f); }


public:

    /** returns the valence (number of incident edges or neighboring vertices)
     of vertex \c v. */
    unsigned int valence(Vertex v) const;

    /// returns the valence of face \c f (its number of vertices)
    unsigned int valence(Face f) const;

    /// find the halfedge from start to end
    Halfedge find_halfedge(Vertex start, Vertex end) const;

    /// find the edge (a,b)
    Edge find_edge(Vertex a, Vertex b) const;

    /// deletes the vertex \c v from the mesh
    void delete_vertex(Vertex v);

    /// deletes the edge \c e from the mesh
    void delete_edge(Edge e);

    /// deletes the face \c f from the mesh
    void delete_face(Face f);

protected:

    Face new_quad(Self::Halfedge * h, int sz);

    /// allocate a new vertex, resize vertex properties accordingly.
    Vertex new_vertex()
    {
        vprops_.push_back();
        return Vertex(vertices_size()-1);
    }

    /// allocate a new edge, resize edge and halfedge properties accordingly.
    Halfedge new_edge(Vertex start, Vertex end)
    {
        assert(start != end);

        eprops_.push_back();
        hprops_.push_back();
        hprops_.push_back();

        Halfedge h0(halfedges_size()-2);
        Halfedge h1(halfedges_size()-1);

        set_vertex(h0, end);   // h0: start->end
        set_vertex(h1, start); // h1: end-->start
        return h0;
    }

    /// allocate a new face, resize face properties accordingly.
    Face new_face()
    {
        fprops_.push_back();
        return Face(faces_size()-1);
    }

    friend std::ostream& operator<<(std::ostream& os, const gsMeshTopology & sm)
    {
        os<<"gsMeshTopology with "<<sm.n_vertices()<<" vertices, "<<sm.n_edges()<<
            " edges and "<<sm.n_faces()<<" faces.\n";
        return os;
    }

protected: //------------------------------------------------ helper functions

    /** make sure that the outgoing halfedge of vertex v is a boundary halfedge
     if v is a boundary vertex. */
    void adjust_outgoing_halfedge(Vertex v);

    /** Called by garbage_collection() once the arrays have been compacted and
        while the index maps are still alive, so that a derived class can remap
        its own handle-valued properties.  The default implementation does
        nothing.

        \param vmap old vertex index -> new vertex handle
        \param hmap old halfedge index -> new halfedge handle
        \param fmap old face index -> new face handle
        \note the edge map is not passed because it is derivable:
        <code>emap[e] = Edge(hmap[halfedge(e,0)].idx()>>1)</code>
    */
    virtual void remap_handles(const Vertex_property<Vertex>     & vmap,
                               const Halfedge_property<Halfedge> & hmap,
                               const Face_property<Face>         & fmap)
    { GISMO_UNUSED(vmap); GISMO_UNUSED(hmap); GISMO_UNUSED(fmap); }

    /// are there deleted vertices, edges or faces?
    bool garbage() const { return garbage_; }

protected:

    gsProperty_container vprops_;
    gsProperty_container hprops_;
    gsProperty_container eprops_;
    gsProperty_container fprops_;
    gsProperty_container mprops_;

    Vertex_property<Vertex_connectivity>      vconn_;
    Halfedge_property<Halfedge_connectivity>  hconn_;
    Face_property<Face_connectivity>          fconn_;

    Vertex_property<bool>  vdeleted_;
    Edge_property<bool>    edeleted_;
    Face_property<bool>    fdeleted_;

    unsigned int deleted_vertices_;
    unsigned int deleted_edges_;
    unsigned int deleted_faces_;
    bool garbage_;

    // helper data for add_face()
    typedef std::pair<Halfedge, Halfedge>  NextCacheEntry;
    typedef std::vector<NextCacheEntry>    NextCache;
    std::vector<Vertex>      add_face_vertices_;
    std::vector<Halfedge>    add_face_halfedges_;
    std::vector<bool>        add_face_is_new_;
    std::vector<bool>        add_face_needs_adjust_;
    NextCache                add_face_next_cache_;
};

/// print a vertex handle
inline std::ostream& operator<<(std::ostream& os, gsMeshTopology::Vertex v)
{ return (os << 'v' << v.idx()); }

/// print a halfedge handle
inline std::ostream& operator<<(std::ostream& os, gsMeshTopology::Halfedge h)
{ return (os << 'h' << h.idx()); }

/// print an edge handle
inline std::ostream& operator<<(std::ostream& os, gsMeshTopology::Edge e)
{ return (os << 'e' << e.idx()); }

/// print a face handle
inline std::ostream& operator<<(std::ostream& os, gsMeshTopology::Face f)
{ return (os << 'f' << f.idx()); }

} // namespace gismo
