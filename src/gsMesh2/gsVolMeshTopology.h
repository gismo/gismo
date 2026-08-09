/** @file gsVolMeshTopology.h

    @brief Half-face mesh topology: the scalar-independent (non-templated) base
    of gsVolMesh, a combinatorial 3-map for polyhedral volume meshes.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsMesh2/gsMeshTopology.h>

#include <string>
#include <vector>

namespace gismo
{

/** \brief A half-face data structure for polyhedral volume meshes: topology only.

    gsVolMeshTopology is a combinatorial 3-map \f$M=(D,\beta_1,\beta_2,\beta_3)\f$
    in the sense of Feng, Wang, Weng and Tong, <em>Compact combinatorial maps: a
    volume mesh data structure</em>, Graphical Models 75 (2013) 149--156, and of
    the CMap3 class of the CGoGN library.  It is built as a sibling of
    gsSurfMeshTopology on the shared 2-map core gsMeshTopology.

    <b>The inherited 2-map is the cell-boundary complex.</b> gsMeshTopology here
    stores the disjoint union of the oriented boundary surfaces of all cells.
    Each cell boundary is a closed oriented genus-0 surface and therefore a legal
    component of the base class, so its darts, its face and vertex circulators,
    its property machinery, add_face() and garbage collection all remain valid
    \em inside a cell.  The paper puts it the same way: "each 3-cell is treated
    locally as a 2-manifold cell complex, which can be represented by a local
    half-edge structure, i.e. a 2D combinatorial map".

    <b>Two tiers of cells.</b> A geometric edge of a volume mesh is surrounded by
    \f$k\f$ faces and so carries \f$2k\f$ darts, which the two-darts-per-edge
    arithmetic of the base (<code>edge(h)==h>>1</code>) cannot express.  The
    inherited cells are therefore the \em cell-local ones and the geometric cells
    are separate:

    | tier       | name       | is                                             | CGoGN     |
    |------------|------------|------------------------------------------------|-----------|
    | geometric  | Vertex     | a point of the mesh                             | Vertex    |
    | geometric  | Edge       | an edge shared by all incident cells            | Edge      |
    | geometric  | Face       | a face shared by one or two cells               | Face      |
    | volumetric | Cell       | a polyhedral cell                               | Volume    |
    | cell-local | Corner     | one use of a Vertex by one Cell                 | Vertex2   |
    | cell-local | Edge_use   | one use of an Edge by one Cell                  | Edge2     |
    | cell-local | Halfedge   | a \em dart: one oriented use of an edge in one half-face | HalfEdge |
    | cell-local | Halfface   | an oriented face, as seen from one Cell         | Face2     |

    Every cell-local entity belongs to exactly one cell; that is what keeps the
    inherited 2-map valid.  The geometric names deliberately hide the inherited
    ones, which stay reachable through the aliases above.

    <b>The three permutations.</b> On the darts,
    - \f$\beta_1\f$ = next_halfedge(), the next dart in the half-face cycle;
    - \f$\beta_2\f$ = opposite_halfedge(), the dart on the adjacent face of the
      \em same cell;
    - \f$\beta_3\f$ = mate(), the dart on the same geometric edge in the opposite
      half-face, i.e. in the neighbouring cell.

    \f$\beta_1\f$ and \f$\beta_2\f$ are inherited; \f$\beta_3\f$ is the only new
    topological relation and everything volumetric follows from it in constant
    time -- see mate(), radial_next(), opposite_halfface().  \f$\beta_3\f$ is
    partial: a dart with no mate is a boundary dart, exactly as a halfedge with
    no face is a boundary halfedge in the surface mesh.

    \sa gsMeshTopology, gsSurfMeshTopology, gsVolMesh
*/
class GISMO_EXPORT gsVolMeshTopology : public gsMeshTopology
{
public:

    /// Dimension-agnostic topology base, holding the cell-boundary complex
    typedef gsMeshTopology Base;

    using Self = gsVolMeshTopology;

public: //------------------------------------------ cell-local (inherited) tier

    /// One use of a Vertex by one Cell: the corner of that cell at that vertex.
    /// There are as many corners at a vertex as there are incident cells.
    /// \sa Vertex
    typedef Base::Vertex   Corner;

    /// One use of an Edge by one Cell: the two darts of that edge on the two
    /// faces of the cell that meet along it.
    /// \sa Edge
    typedef Base::Edge     Edge_use;

    /// A \em dart: one oriented use of an edge inside one half-face.
    typedef Base::Halfedge Halfedge;

    /// A \em dart (combinatorial-map spelling of Halfedge).
    typedef Base::Halfedge Dart;

    /// An oriented face as seen from one Cell.  A geometric Face carries two of
    /// them in the interior of the mesh and one on its boundary.
    /// \sa Face
    typedef Base::Face     Halfface;

    /// Property of type T attached to the corners
    template <class T> using Corner_property   = Base::Vertex_property<T>;
    /// Property of type T attached to the edge-uses
    template <class T> using Edge_use_property = Base::Edge_property<T>;
    /// Property of type T attached to the half-faces
    template <class T> using Halfface_property = Base::Face_property<T>;

    // Halfedge_property and Mesh_property are inherited unchanged

public: //------------------------------------------------- geometric tier

    /** Base class of the geometric and volumetric handles.

        Unlike gsMeshTopology::Base_handle the comparison operators are \em typed.
        With two tiers of 0-, 1- and 2-cells in play, comparing handles of
        different kinds is always a bug, and here it fails to compile rather than
        silently comparing indices.
    */
    template <class D>
    class Vol_handle
    {
    public:
        /// constructor
        explicit Vol_handle(int _idx=-1) : idx_(_idx) {}

        /// Get the underlying index of this handle
        int idx() const { return idx_; }

        /// reset handle to be invalid (index=-1)
        void reset() { idx_=-1; }

        /// return whether the handle is valid, i.e. the index is not -1
        bool is_valid() const { return idx_ != -1; }

        /// are two handles equal?
        bool operator==(const D& rhs) const { return idx_ == rhs.idx_; }
        /// are two handles different?
        bool operator!=(const D& rhs) const { return idx_ != rhs.idx_; }
        /// compare operator, useful for sorting handles
        bool operator< (const D& rhs) const { return idx_ <  rhs.idx_; }

    protected:
        int idx_;
    };

    /// A geometric vertex, shared by all incident cells.  This is the entity
    /// that carries the point coordinates in gsVolMesh.
    /// \sa Corner
    struct GISMO_EXPORT Vertex : public Vol_handle<Vertex>
    { explicit Vertex(int _idx=-1) : Vol_handle<Vertex>(_idx) {} };

    /// A geometric edge, shared by all incident cells.
    /// \sa Edge_use
    struct GISMO_EXPORT Edge : public Vol_handle<Edge>
    { explicit Edge(int _idx=-1) : Vol_handle<Edge>(_idx) {} };

    /// A geometric face, carrying one (boundary) or two (interior) half-faces.
    /// \sa Halfface
    struct GISMO_EXPORT Face : public Vol_handle<Face>
    { explicit Face(int _idx=-1) : Vol_handle<Face>(_idx) {} };

    /// A polyhedral cell (a 3-cell).
    struct GISMO_EXPORT Cell : public Vol_handle<Cell>
    { explicit Cell(int _idx=-1) : Vol_handle<Cell>(_idx) {} };

public: //-------------------------------------------- geometric connectivity

    /// Connectivity of a geometric vertex
    struct Vertex_connectivity
    {
        /// one corner at this vertex; the others follow the corner ring
        Corner corner_;
    };

    /// Connectivity of a geometric edge
    struct Edge_connectivity
    {
        /// one dart on this edge; the others follow the radial cycle
        Halfedge halfedge_;
    };

    /// Connectivity of a geometric face
    struct Face_connectivity
    {
        /// one of the (at most two) half-faces of this face
        Halfface halfface_;
    };

    /// Connectivity of a cell
    struct Cell_connectivity
    {
        /// one half-face of the cell; the others follow the half-face ring
        Halfface halfface_;
        /// one corner of the cell; the others follow the cell corner ring
        Corner   corner_;
        /// one edge-use of the cell; the others follow the cell edge-use ring
        Edge_use edge_use_;
    };

public: //--------------------------------------------- geometric properties

    /// Property of type \a T attached to the geometric entities of kind \a H
    template <class H, class T>
    class Entity_property : public gsProperty<T>
    {
    public:
        /// default constructor
        explicit Entity_property() {}
        explicit Entity_property(gsProperty<T> p) : gsProperty<T>(p) {}

        /// access the data stored for entity \c h
        typename gsProperty<T>::reference operator[](H h)
        { return gsProperty<T>::operator[](h.idx()); }

        /// access the data stored for entity \c h
        typename gsProperty<T>::const_reference operator[](H h) const
        { return gsProperty<T>::operator[](h.idx()); }
    };

    /// Vertex property of type T.  Hides gsMeshTopology::Vertex_property, which
    /// stays available as Corner_property.
    template <class T> using Vertex_property = Entity_property<Vertex, T>;
    /// Edge property of type T.  Hides gsMeshTopology::Edge_property, which
    /// stays available as Edge_use_property.
    template <class T> using Edge_property   = Entity_property<Edge  , T>;
    /// Face property of type T.  Hides gsMeshTopology::Face_property, which
    /// stays available as Halfface_property.
    template <class T> using Face_property   = Entity_property<Face  , T>;
    /// Cell property of type T
    template <class T> using Cell_property   = Entity_property<Cell  , T>;

public: //------------------------------------------------------- iterators

    /// Linear iterator over all geometric entities of kind \a H, skipping the
    /// deleted ones.
    /// \sa vertices(), edges(), faces(), cells()
    template <class H>
    class Entity_iterator
    {
    public:
        /// Default constructor
        Entity_iterator(H h=H(), const gsVolMeshTopology* m=NULL) : hnd_(h), mesh_(m)
        {
            if (mesh_ && mesh_->garbage())
                while (mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) step(1);
        }

        /// get the entity the iterator refers to
        H operator*() const { return hnd_; }

        /// are two iterators equal?
        bool operator==(const Entity_iterator& rhs) const { return hnd_==rhs.hnd_; }
        /// are two iterators different?
        bool operator!=(const Entity_iterator& rhs) const { return hnd_!=rhs.hnd_; }
        /// is this iterator before \c rhs ?
        bool operator< (const Entity_iterator& rhs) const { return hnd_< rhs.hnd_; }

        /// pre-increment
        Entity_iterator& operator++() { advance( 1); return *this; }
        /// pre-decrement
        Entity_iterator& operator--() { advance(-1); return *this; }

    private:
        void step(int d) { hnd_ = H(hnd_.idx()+d); }

        void advance(int d)
        {
            step(d);
            GISMO_ASSERT(mesh_, "gsVolMeshTopology: iterator without a mesh");
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_))
                step(d);
        }

        H                        hnd_;
        const gsVolMeshTopology* mesh_;
    };

    typedef Entity_iterator<Vertex> Vertex_iterator;
    typedef Entity_iterator<Edge>   Edge_iterator;
    typedef Entity_iterator<Face>   Face_iterator;
    typedef Entity_iterator<Cell>   Cell_iterator;

    /// Helper class acting as a container for C++11 range-based for-loops
    template <class I>
    class Entity_container
    {
    public:
        Entity_container(I _begin, I _end) : begin_(_begin), end_(_end) {}
        I begin() const { return begin_; }
        I end()   const { return end_;   }
    private:
        I begin_, end_;
    };

    typedef Entity_container<Vertex_iterator> Vertex_container;
    typedef Entity_container<Edge_iterator>   Edge_container;
    typedef Entity_container<Face_iterator>   Face_container;
    typedef Entity_container<Cell_iterator>   Cell_container;

public: //--------------------------------------- circulators: around a Cell

    /// Circulates through the half-faces of a cell, following its half-face ring.
    /// It also acts as a container-concept for C++11 range-based for loops.
    /// \sa halffaces(Cell)
    class Halfface_around_cell_circulator
    {
    public:
        Halfface_around_cell_circulator(const gsVolMeshTopology* m=NULL, Cell c=Cell())
        : mesh_(m), active_(true)
        { if (mesh_ && c.is_valid()) halfface_ = mesh_->halfface(c); }

        bool operator==(const Halfface_around_cell_circulator& rhs) const
        { return (active_ && mesh_==rhs.mesh_ && halfface_==rhs.halfface_); }
        bool operator!=(const Halfface_around_cell_circulator& rhs) const
        { return !operator==(rhs); }

        /// pre-increment
        Halfface_around_cell_circulator& operator++()
        {
            GISMO_ASSERT(mesh_ && halfface_.is_valid(), "invalid circulator");
            halfface_ = mesh_->next_halfface(halfface_);
            active_ = true;
            return *this;
        }

        /// get the half-face the circulator refers to
        Halfface operator*() const { return halfface_; }
        const Halfface* operator->() const { return &halfface_; }

        /// cast to bool: true if the cell has any half-face
        operator bool() const { return halfface_.is_valid(); }

        // helpers for C++11 range-based for-loops
        Halfface_around_cell_circulator& begin() { active_=!halfface_.is_valid(); return *this; }
        Halfface_around_cell_circulator& end()   { active_=true;  return *this; }

    private:
        const gsVolMeshTopology* mesh_;
        Halfface                 halfface_;
        bool                     active_;
    };

    /// Circulates through the corners of a cell, following its corner ring.
    /// Every vertex of the cell is visited exactly once.
    /// \sa corners(Cell), vertices(Cell)
    class Corner_around_cell_circulator
    {
    public:
        Corner_around_cell_circulator(const gsVolMeshTopology* m=NULL, Cell c=Cell())
        : mesh_(m), active_(true)
        { if (mesh_ && c.is_valid()) corner_ = mesh_->corner(c); }

        bool operator==(const Corner_around_cell_circulator& rhs) const
        { return (active_ && mesh_==rhs.mesh_ && corner_==rhs.corner_); }
        bool operator!=(const Corner_around_cell_circulator& rhs) const
        { return !operator==(rhs); }

        /// pre-increment
        Corner_around_cell_circulator& operator++()
        {
            GISMO_ASSERT(mesh_ && corner_.is_valid(), "invalid circulator");
            corner_ = mesh_->next_corner_in_cell(corner_);
            active_ = true;
            return *this;
        }

        /// get the corner the circulator refers to
        Corner operator*() const { return corner_; }
        const Corner* operator->() const { return &corner_; }

        operator bool() const { return corner_.is_valid(); }

        Corner_around_cell_circulator& begin() { active_=!corner_.is_valid(); return *this; }
        Corner_around_cell_circulator& end()   { active_=true;  return *this; }

    private:
        const gsVolMeshTopology* mesh_;
        Corner                   corner_;
        bool                     active_;
    };

    /// Circulates through the edge-uses of a cell, following its edge-use ring.
    /// Every edge of the cell is visited exactly once.
    /// \sa edge_uses(Cell), edges(Cell)
    class Edge_use_around_cell_circulator
    {
    public:
        Edge_use_around_cell_circulator(const gsVolMeshTopology* m=NULL, Cell c=Cell())
        : mesh_(m), active_(true)
        { if (mesh_ && c.is_valid()) edge_use_ = mesh_->edge_use(c); }

        bool operator==(const Edge_use_around_cell_circulator& rhs) const
        { return (active_ && mesh_==rhs.mesh_ && edge_use_==rhs.edge_use_); }
        bool operator!=(const Edge_use_around_cell_circulator& rhs) const
        { return !operator==(rhs); }

        /// pre-increment
        Edge_use_around_cell_circulator& operator++()
        {
            GISMO_ASSERT(mesh_ && edge_use_.is_valid(), "invalid circulator");
            edge_use_ = mesh_->next_edge_use_in_cell(edge_use_);
            active_ = true;
            return *this;
        }

        /// get the edge-use the circulator refers to
        Edge_use operator*() const { return edge_use_; }
        const Edge_use* operator->() const { return &edge_use_; }

        operator bool() const { return edge_use_.is_valid(); }

        Edge_use_around_cell_circulator& begin() { active_=!edge_use_.is_valid(); return *this; }
        Edge_use_around_cell_circulator& end()   { active_=true;  return *this; }

    private:
        const gsVolMeshTopology* mesh_;
        Edge_use                 edge_use_;
        bool                     active_;
    };

    /** Circulates through the cells adjacent to a cell across its half-faces.

        Boundary half-faces are skipped, so the circulator is empty for an
        isolated cell.  A neighbour that shares several faces with the cell is
        reported once per shared face.
        \sa cells(Cell)
    */
    class Cell_around_cell_circulator
    {
    public:
        Cell_around_cell_circulator(const gsVolMeshTopology* m=NULL, Cell c=Cell())
        : mesh_(m), active_(true)
        {
            if (mesh_ && c.is_valid())
            {
                halfface_ = mesh_->halfface(c);
                if (halfface_.is_valid() && mesh_->is_boundary(halfface_))
                {
                    const Halfface start = halfface_;
                    do halfface_ = mesh_->next_halfface(halfface_);
                    while (halfface_ != start && mesh_->is_boundary(halfface_));
                    if (mesh_->is_boundary(halfface_)) halfface_.reset(); // all boundary
                }
            }
        }

        bool operator==(const Cell_around_cell_circulator& rhs) const
        { return (active_ && mesh_==rhs.mesh_ && halfface_==rhs.halfface_); }
        bool operator!=(const Cell_around_cell_circulator& rhs) const
        { return !operator==(rhs); }

        /// pre-increment
        Cell_around_cell_circulator& operator++()
        {
            GISMO_ASSERT(mesh_ && halfface_.is_valid(), "invalid circulator");
            do halfface_ = mesh_->next_halfface(halfface_);
            while (mesh_->is_boundary(halfface_));
            active_ = true;
            return *this;
        }

        /// get the neighbouring cell the circulator refers to
        Cell operator*() const
        {
            GISMO_ASSERT(mesh_ && halfface_.is_valid(), "invalid circulator");
            return mesh_->opposite_cell(halfface_);
        }

        /// the half-face of the current cell that is being crossed
        Halfface halfface() const { return halfface_; }

        operator bool() const { return halfface_.is_valid(); }

        Cell_around_cell_circulator& begin() { active_=!halfface_.is_valid(); return *this; }
        Cell_around_cell_circulator& end()   { active_=true;  return *this; }

    private:
        const gsVolMeshTopology* mesh_;
        Halfface                 halfface_;
        bool                     active_;
    };

public: //----------------------------------- circulators: around a Halfface

    /// Circulates through the darts of a half-face, i.e. along its oriented
    /// boundary.  This is the face circulator of the inherited 2-map.
    /// \sa halfedges(Halfface)
    typedef Base::Halfedge_around_face_circulator Halfedge_around_halfface_circulator;

    /// Circulates through the corners of a half-face.
    /// \sa corners(Halfface)
    typedef Base::Vertex_around_face_circulator   Corner_around_halfface_circulator;

    /// Circulates through the geometric vertices of a half-face.
    /// \sa vertices(Halfface)
    class Vertex_around_halfface_circulator
    {
    public:
        Vertex_around_halfface_circulator(const gsVolMeshTopology* m=NULL,
                                          Halfface hf=Halfface())
        : mesh_(m), active_(true)
        { if (mesh_ && hf.is_valid()) halfedge_ = mesh_->halfedge(hf); }

        bool operator==(const Vertex_around_halfface_circulator& rhs) const
        { return (active_ && mesh_==rhs.mesh_ && halfedge_==rhs.halfedge_); }
        bool operator!=(const Vertex_around_halfface_circulator& rhs) const
        { return !operator==(rhs); }

        /// pre-increment (advances along the half-face cycle)
        Vertex_around_halfface_circulator& operator++()
        {
            GISMO_ASSERT(mesh_ && halfedge_.is_valid(), "invalid circulator");
            halfedge_ = mesh_->next_halfedge(halfedge_);
            active_ = true;
            return *this;
        }

        /// pre-decrement
        Vertex_around_halfface_circulator& operator--()
        {
            GISMO_ASSERT(mesh_ && halfedge_.is_valid(), "invalid circulator");
            halfedge_ = mesh_->prev_halfedge(halfedge_);
            return *this;
        }

        /// get the vertex the circulator refers to
        Vertex operator*() const
        {
            GISMO_ASSERT(mesh_ && halfedge_.is_valid(), "invalid circulator");
            return mesh_->to_vertex(halfedge_);
        }

        /// the dart the circulator currently sits on
        Halfedge halfedge() const { return halfedge_; }

        operator bool() const { return halfedge_.is_valid(); }

        Vertex_around_halfface_circulator& begin() { active_=false; return *this; }
        Vertex_around_halfface_circulator& end()   { active_=true;  return *this; }

    private:
        const gsVolMeshTopology* mesh_;
        Halfedge                 halfedge_;
        bool                     active_;
    };

public: //--------------------------------------- circulators: around an Edge

    /** Circulates through the darts on a geometric edge: the \em radial cycle.

        \f$\beta_2\f$ and \f$\beta_3\f$ each reverse the orientation of a dart,
        so their composition radial_next() = \f$\beta_2\circ\beta_3\f$ preserves
        it and steps into the next cell around the edge.  The orbit therefore
        reports exactly one dart per incident \em cell, all oriented the same
        way, and never repeats a cell.

        For an interior edge the orbit is a cycle of length valence(Edge).  For a
        boundary edge it is an open chain: the constructor rewinds with
        radial_prev() to the free end, and the circulator stops once the current
        dart has no mate.
        \sa halfedges(Edge), halffaces(Edge), cells(Edge)
    */
    class Halfedge_around_edge_circulator
    {
    public:
        Halfedge_around_edge_circulator(const gsVolMeshTopology* m=NULL, Edge e=Edge())
        : mesh_(m), active_(true), done_(false)
        {
            if (mesh_ && e.is_valid())
            {
                halfedge_ = mesh_->halfedge(e);
                // rewind to the free end of an open (boundary) radial chain
                const Halfedge start = halfedge_;
                Halfedge p = mesh_->radial_prev(halfedge_);
                while (p.is_valid() && p != start)
                {
                    halfedge_ = p;
                    p = mesh_->radial_prev(halfedge_);
                }
            }
        }

        bool operator==(const Halfedge_around_edge_circulator& rhs) const
        { return (active_ && mesh_==rhs.mesh_ && (done_ || halfedge_==rhs.halfedge_)); }
        bool operator!=(const Halfedge_around_edge_circulator& rhs) const
        { return !operator==(rhs); }

        /// pre-increment (rotates into the next cell around the edge)
        Halfedge_around_edge_circulator& operator++()
        {
            GISMO_ASSERT(mesh_ && halfedge_.is_valid(), "invalid circulator");
            const Halfedge n = mesh_->radial_next(halfedge_);
            if (n.is_valid()) halfedge_ = n; else done_ = true;
            active_ = true;
            return *this;
        }

        /// pre-decrement (rotates into the previous cell around the edge)
        Halfedge_around_edge_circulator& operator--()
        {
            GISMO_ASSERT(mesh_ && halfedge_.is_valid(), "invalid circulator");
            const Halfedge p = mesh_->radial_prev(halfedge_);
            if (p.is_valid()) halfedge_ = p;
            return *this;
        }

        /// get the dart the circulator refers to
        Halfedge operator*() const { return halfedge_; }
        const Halfedge* operator->() const { return &halfedge_; }

        operator bool() const { return halfedge_.is_valid() && !done_; }

        Halfedge_around_edge_circulator& begin() { active_=!halfedge_.is_valid(); return *this; }
        Halfedge_around_edge_circulator& end()   { active_=true;  return *this; }

    private:
        const gsVolMeshTopology* mesh_;
        Halfedge                 halfedge_;
        bool                     active_, done_;
    };

    /** Circulates through the half-faces incident to a geometric edge.

        It follows the radial cycle and therefore reports one half-face per
        incident cell: the face by which that cell is left.  The remaining
        half-faces at the edge are their opposites and, on a boundary edge, the
        two boundary half-faces.
        \sa halffaces(Edge)
    */
    class Halfface_around_edge_circulator
    {
    public:
        Halfface_around_edge_circulator(const gsVolMeshTopology* m=NULL, Edge e=Edge())
        : mesh_(m), circ_(m,e) {}

        bool operator==(const Halfface_around_edge_circulator& rhs) const { return circ_==rhs.circ_; }
        bool operator!=(const Halfface_around_edge_circulator& rhs) const { return circ_!=rhs.circ_; }

        Halfface_around_edge_circulator& operator++() { ++circ_; return *this; }
        Halfface_around_edge_circulator& operator--() { --circ_; return *this; }

        /// get the half-face the circulator refers to
        Halfface operator*() const
        {
            GISMO_ASSERT(mesh_, "invalid circulator");
            return mesh_->halfface(*circ_);
        }

        operator bool() const { return static_cast<bool>(circ_); }

        Halfface_around_edge_circulator& begin() { circ_.begin(); return *this; }
        Halfface_around_edge_circulator& end()   { circ_.end();   return *this; }

    private:
        const gsVolMeshTopology*        mesh_;
        Halfedge_around_edge_circulator circ_;
    };

    /// Circulates through the cells incident to a geometric edge.
    /// \sa cells(Edge), valence(Edge)
    class Cell_around_edge_circulator
    {
    public:
        Cell_around_edge_circulator(const gsVolMeshTopology* m=NULL, Edge e=Edge())
        : mesh_(m), circ_(m,e) {}

        bool operator==(const Cell_around_edge_circulator& rhs) const { return circ_==rhs.circ_; }
        bool operator!=(const Cell_around_edge_circulator& rhs) const { return circ_!=rhs.circ_; }

        Cell_around_edge_circulator& operator++() { ++circ_; return *this; }
        Cell_around_edge_circulator& operator--() { --circ_; return *this; }

        /// get the cell the circulator refers to
        Cell operator*() const
        {
            GISMO_ASSERT(mesh_, "invalid circulator");
            return mesh_->cell(mesh_->halfface(*circ_));
        }

        operator bool() const { return static_cast<bool>(circ_); }

        Cell_around_edge_circulator& begin() { circ_.begin(); return *this; }
        Cell_around_edge_circulator& end()   { circ_.end();   return *this; }

    private:
        const gsVolMeshTopology*        mesh_;
        Halfedge_around_edge_circulator circ_;
    };

public: //-------------------------------------- circulators: around a Vertex

    /** Circulates through the corners at a geometric vertex, i.e. through its
        uses by the incident cells.

        Each incident cell contributes exactly one corner, so no de-duplication
        is involved and cells(Vertex) is a faithful traversal of the cell star.
        \sa corners(Vertex), cells(Vertex), valence(Vertex)
    */
    class Corner_around_vertex_circulator
    {
    public:
        Corner_around_vertex_circulator(const gsVolMeshTopology* m=NULL, Vertex v=Vertex())
        : mesh_(m), active_(true)
        { if (mesh_ && v.is_valid()) corner_ = mesh_->corner(v); }

        bool operator==(const Corner_around_vertex_circulator& rhs) const
        { return (active_ && mesh_==rhs.mesh_ && corner_==rhs.corner_); }
        bool operator!=(const Corner_around_vertex_circulator& rhs) const
        { return !operator==(rhs); }

        /// pre-increment
        Corner_around_vertex_circulator& operator++()
        {
            GISMO_ASSERT(mesh_ && corner_.is_valid(), "invalid circulator");
            corner_ = mesh_->next_corner_at_vertex(corner_);
            active_ = true;
            return *this;
        }

        /// get the corner the circulator refers to
        Corner operator*() const { return corner_; }
        const Corner* operator->() const { return &corner_; }

        operator bool() const { return corner_.is_valid(); }

        Corner_around_vertex_circulator& begin() { active_=!corner_.is_valid(); return *this; }
        Corner_around_vertex_circulator& end()   { active_=true;  return *this; }

    private:
        const gsVolMeshTopology* mesh_;
        Corner                   corner_;
        bool                     active_;
    };

    /// Circulates through the cells incident to a geometric vertex.
    /// \sa cells(Vertex)
    class Cell_around_vertex_circulator
    {
    public:
        Cell_around_vertex_circulator(const gsVolMeshTopology* m=NULL, Vertex v=Vertex())
        : mesh_(m), circ_(m,v) {}

        bool operator==(const Cell_around_vertex_circulator& rhs) const { return circ_==rhs.circ_; }
        bool operator!=(const Cell_around_vertex_circulator& rhs) const { return circ_!=rhs.circ_; }

        Cell_around_vertex_circulator& operator++() { ++circ_; return *this; }

        /// get the cell the circulator refers to
        Cell operator*() const
        {
            GISMO_ASSERT(mesh_, "invalid circulator");
            return mesh_->cell(*circ_);
        }

        operator bool() const { return static_cast<bool>(circ_); }

        Cell_around_vertex_circulator& begin() { circ_.begin(); return *this; }
        Cell_around_vertex_circulator& end()   { circ_.end();   return *this; }

    private:
        const gsVolMeshTopology*        mesh_;
        Corner_around_vertex_circulator circ_;
    };

    /** Circulates through all darts emanating from a geometric vertex.

        For every corner at the vertex it walks the outgoing darts of that
        corner inside its cell, then steps to the next corner.  A geometric edge
        at the vertex is therefore reported once per incident half-face; use
        edges(Vertex) for the de-duplicated list of edges.
        \sa halfedges(Vertex), edges(Vertex)
    */
    class Halfedge_around_vertex_circulator
    {
    public:
        Halfedge_around_vertex_circulator(const gsVolMeshTopology* m=NULL,
                                          Vertex v=Vertex());

        bool operator==(const Halfedge_around_vertex_circulator& rhs) const
        { return (active_ && mesh_==rhs.mesh_ && halfedge_==rhs.halfedge_); }
        bool operator!=(const Halfedge_around_vertex_circulator& rhs) const
        { return !operator==(rhs); }

        /// pre-increment: advance inside the star of the current corner, and
        /// step to the next corner once that star is exhausted
        Halfedge_around_vertex_circulator& operator++();

        /// get the dart the circulator refers to
        Halfedge operator*() const { return halfedge_; }
        const Halfedge* operator->() const { return &halfedge_; }

        /// the corner the current dart emanates from
        Corner corner() const { return corner_; }

        operator bool() const { return halfedge_.is_valid(); }

        Halfedge_around_vertex_circulator& begin() { active_=!halfedge_.is_valid(); return *this; }
        Halfedge_around_vertex_circulator& end()   { active_=true;  return *this; }

    private:
        const gsVolMeshTopology* mesh_;
        Vertex                   vertex_;
        Corner                   corner_, corner_start_;
        Halfedge                 halfedge_, corner_first_;
        bool                     active_;
    };

public:

    /// \name Construct, destruct, assignment
    //@{

    /// default constructor
    gsVolMeshTopology();

    virtual ~gsVolMeshTopology();

    /// copy constructor: copies \c rhs to \c *this. performs a deep copy of all properties.
    gsVolMeshTopology(const gsVolMeshTopology& rhs) : Base() { operator=(rhs); }

    /// assign \c rhs to \c *this. performs a deep copy of all properties.
    gsVolMeshTopology& operator=(const gsVolMeshTopology& rhs);

    /// move \c rhs to \c *this
    gsVolMeshTopology& operator=(gsVolMeshTopology&& rhs) noexcept;

    /// assign \c rhs to \c *this. does not copy custom properties.
    gsVolMeshTopology& assign(const gsVolMeshTopology& rhs);

    //@}

public:

    /// \name Adding vertices and cells
    //@{

    /// add a new (isolated) vertex
    Vertex add_vertex();

    /// add \a nverts new vertices at once
    void add_batch_vertices(size_t nverts);

    /** @brief add a new cell bounded by the given faces.

        Each entry of \a faces is the vertex loop of one boundary face of the
        cell, oriented so that its normal points \em out of the cell.  The
        cell's corners, half-faces, edge-uses and darts are created; each face
        and each edge is matched against the existing geometry and, whenever a
        free opposite half-face is found, the two half-faces are sewn.

        \returns the new cell, or an invalid handle if the input is not a closed
        orientable polyhedron, or if gluing would give a face three incident
        cells.
        \sa add_tet, add_hex, add_prism, add_pyramid
    */
    Cell add_cell(const std::vector< std::vector<Vertex> >& faces);

    /// add a tetrahedron; seen from \a v3 the triangle \a v0 \a v1 \a v2 is
    /// counter-clockwise (VTK_TETRA ordering)
    Cell add_tet(Vertex v0, Vertex v1, Vertex v2, Vertex v3);

    /// add a hexahedron in VTK_HEXAHEDRON / CGNS ordering: \a v0..\a v3 is the
    /// bottom quad seen from outside, \a v4..\a v7 the matching top quad
    Cell add_hex(Vertex v0, Vertex v1, Vertex v2, Vertex v3,
                 Vertex v4, Vertex v5, Vertex v6, Vertex v7);

    /// add a triangular prism in VTK_WEDGE ordering: \a v0..\a v2 is the bottom
    /// triangle seen from outside, \a v3..\a v5 the matching top triangle
    Cell add_prism(Vertex v0, Vertex v1, Vertex v2,
                   Vertex v3, Vertex v4, Vertex v5);

    /// add a pyramid in VTK_PYRAMID ordering: \a v0..\a v3 is the quad base
    /// seen from outside and \a v4 the apex
    Cell add_pyramid(Vertex v0, Vertex v1, Vertex v2, Vertex v3, Vertex v4);

    /** glue the two half-faces \a hf0 and \a hf1.

        Both must be free (on the boundary) and must carry the same vertex loop
        in opposite orientations.  Their geometric faces are merged and
        \f$\beta_3\f$ is set on every dart of the two loops.
        \sa unsew
    */
    void sew(Halfface hf0, Halfface hf1);

    /// undo the gluing of \a hf and its opposite, splitting their geometric
    /// face into two boundary faces
    void unsew(Halfface hf);

    //@}

public:

    /// \name Sizes
    //@{

    /// number of (deleted and valid) vertices in the mesh
    unsigned int vertices_size()  const { return (unsigned int) vertprops_.size(); }
    /// number of (deleted and valid) edges in the mesh
    unsigned int edges_size()     const { return (unsigned int) edgeprops_.size(); }
    /// number of (deleted and valid) faces in the mesh
    unsigned int faces_size()     const { return (unsigned int) faceprops_.size(); }
    /// number of (deleted and valid) cells in the mesh
    unsigned int cells_size()     const { return (unsigned int) cellprops_.size(); }
    /// number of (deleted and valid) corners in the mesh
    unsigned int corners_size()   const { return Base::vertices_size(); }
    /// number of (deleted and valid) edge-uses in the mesh
    unsigned int edge_uses_size() const { return Base::edges_size(); }
    /// number of (deleted and valid) half-faces in the mesh
    unsigned int halffaces_size() const { return Base::faces_size(); }

    /// number of vertices in the mesh
    unsigned int n_vertices()  const { return vertices_size() - deleted_verts_; }
    /// number of edges in the mesh
    unsigned int n_edges()     const { return edges_size()    - deleted_geoedges_; }
    /// number of faces in the mesh
    unsigned int n_faces()     const { return faces_size()    - deleted_geofaces_; }
    /// number of cells in the mesh
    unsigned int n_cells()     const { return cells_size()    - deleted_cells_; }
    /// number of corners in the mesh
    unsigned int n_corners()   const { return Base::n_vertices(); }
    /// number of edge-uses in the mesh
    unsigned int n_edge_uses() const { return Base::n_edges(); }
    /// number of half-faces in the mesh
    unsigned int n_halffaces() const { return Base::n_faces(); }

    /// returns true iff the mesh is empty, i.e. has no vertices
    bool empty() const { return n_vertices() == 0; }

    /// clear mesh: remove all vertices, edges, faces and cells
    void clear();

    /// remove unused memory from the arrays
    void free_memory();

    /// reserve memory (mainly used in file readers)
    void reserve(unsigned int nvertices, unsigned int nedges,
                 unsigned int nfaces,    unsigned int ncells);

    //@}

public:

    /// \name Validity and deletion flags
    //@{

    // keep the cell-local overloads of the base next to the geometric ones
    using Base::is_valid;
    using Base::is_deleted;

    bool is_valid(Vertex v) const { return (0<=v.idx()) && (v.idx()<(int)vertices_size()); }
    bool is_valid(Edge   e) const { return (0<=e.idx()) && (e.idx()<(int)edges_size());    }
    bool is_valid(Face   f) const { return (0<=f.idx()) && (f.idx()<(int)faces_size());    }
    bool is_valid(Cell   c) const { return (0<=c.idx()) && (c.idx()<(int)cells_size());    }

    bool is_deleted(Vertex v) const { return vertdeleted_[v]; }
    bool is_deleted(Edge   e) const { return edgedeleted_[e]; }
    bool is_deleted(Face   f) const { return facedeleted_[f]; }
    bool is_deleted(Cell   c) const { return celldeleted_[c]; }

    /// are there deleted entities in either tier?
    bool garbage() const { return Base::garbage() || volgarbage_; }

    /// remove all deleted entities and compact both tiers
    virtual void garbage_collection();

    //@}

public:

    /// \name Cross-tier accessors
    //@{

    // the base overloads halfedge(Corner), halfedge(Edge_use,i) and
    // halfedge(Halfface); keep them next to halfedge(Edge)
    using Base::halfedge;

    /// the geometric vertex that corner \c cn is a use of
    Vertex   vertex(Corner cn)      const { return cvertex_[cn]; }
    /// one corner at vertex \c v; the others follow the corner ring
    Corner   corner(Vertex v)       const { return vertconn_[v].corner_; }
    /// the cell that corner \c cn belongs to
    Cell     cell(Corner cn)        const { return cell(Base::face(Base::halfedge(cn))); }

    /// the geometric edge that dart \c h runs along
    Edge     edge(Halfedge h)       const { return uedge_[Base::edge(h)]; }
    /// the geometric edge that edge-use \c u is a use of
    Edge     edge(Edge_use u)       const { return uedge_[u]; }
    /// the edge-use that dart \c h belongs to
    Edge_use edge_use(Halfedge h)   const { return Base::edge(h); }
    /// one dart on the geometric edge \c e; the others follow the radial cycle
    Halfedge halfedge(Edge e)       const { return edgeconn_[e].halfedge_; }

    /// the half-face that dart \c h belongs to
    Halfface halfface(Halfedge h)   const { return Base::face(h); }
    /// the geometric face that half-face \c hf is a use of
    Face     face(Halfface hf)      const { return hfface_[hf]; }
    /** the \c i'th half-face of face \c f; \c i has to be 0 or 1.

        Half-face 0 always exists.  For a boundary face half-face 1 is an
        invalid handle.
    */
    Halfface halfface(Face f, unsigned int i=0) const
    {
        GISMO_ASSERT(i<=1, "A face has at most two half-faces");
        return 0==i ? faceconn_[f].halfface_
                    : opposite_halfface(faceconn_[f].halfface_);
    }

    /// the cell incident to half-face \c hf
    Cell     cell(Halfface hf)      const { return hfcell_[hf]; }
    /// the \c i'th cell incident to face \c f; \c i has to be 0 or 1
    Cell     cell(Face f, unsigned int i=0) const
    {
        const Halfface hf = halfface(f,i);
        return hf.is_valid() ? cell(hf) : Cell();
    }

    /// one half-face of cell \c c; the others follow the half-face ring
    Halfface halfface(Cell c)       const { return cellconn_[c].halfface_; }
    /// one corner of cell \c c; the others follow the corner ring
    Corner   corner(Cell c)         const { return cellconn_[c].corner_; }
    /// one edge-use of cell \c c; the others follow the edge-use ring
    Edge_use edge_use(Cell c)       const { return cellconn_[c].edge_use_; }

    /// the corner that dart \c h points to
    Corner   to_corner(Halfedge h)   const { return Base::to_vertex(h); }
    /// the corner that dart \c h emanates from
    Corner   from_corner(Halfedge h) const { return Base::from_vertex(h); }
    /// the geometric vertex that dart \c h points to
    Vertex   to_vertex(Halfedge h)   const { return cvertex_[Base::to_vertex(h)]; }
    /// the geometric vertex that dart \c h emanates from
    Vertex   from_vertex(Halfedge h) const { return cvertex_[Base::from_vertex(h)]; }

    /// the \c i'th vertex of geometric edge \c e; \c i has to be 0 or 1
    Vertex   vertex(Edge e, unsigned int i) const
    {
        GISMO_ASSERT(i<=1, "An edge has two vertices");
        const Halfedge h = halfedge(e);
        return 0==i ? from_vertex(h) : to_vertex(h);
    }

    /// the \c i'th corner of edge-use \c u; \c i has to be 0 or 1
    Corner   corner(Edge_use u, unsigned int i) const { return Base::vertex(u,i); }

    //@}

public:

    /// \name Navigation
    //@{

    /** the dart on the same geometric edge in the opposite half-face, i.e. in
        the neighbouring cell: the map \f$\beta_3\f$.

        Invalid for a boundary dart.  \f$\beta_3\f$ is an involution and reverses
        the half-face cycle, i.e. \f$(\beta_1\circ\beta_3)^2 = \mathrm{id}\f$.
    */
    Halfedge mate(Halfedge h) const { return hmate_[h]; }

    /// rotate around the geometric edge of \a h into the next cell, preserving
    /// the orientation of the dart: \f$\beta_2\circ\beta_3\f$
    Halfedge radial_next(Halfedge h) const
    {
        const Halfedge m = hmate_[h];
        return m.is_valid() ? Base::opposite_halfedge(m) : Halfedge();
    }

    /// rotate around the geometric edge of \a h into the previous cell
    Halfedge radial_prev(Halfedge h) const
    { return hmate_[Base::opposite_halfedge(h)]; }

    /// the opposite half-face of \a hf; invalid iff \a hf is on the boundary
    Halfface opposite_halfface(Halfface hf) const
    {
        const Halfedge m = hmate_[Base::halfedge(hf)];
        return m.is_valid() ? Base::face(m) : Halfface();
    }

    /// the cell on the other side of \a hf; invalid iff \a hf is on the boundary
    Cell opposite_cell(Halfface hf) const
    {
        const Halfface o = opposite_halfface(hf);
        return o.is_valid() ? cell(o) : Cell();
    }

    /// next half-face in the half-face ring of the cell
    Halfface next_halfface(Halfface hf) const { return hfnext_[hf]; }
    /// next corner in the corner ring of the cell
    Corner   next_corner_in_cell(Corner cn) const { return ccnext_[cn]; }
    /// next corner in the corner ring of the geometric vertex
    Corner   next_corner_at_vertex(Corner cn) const { return cvnext_[cn]; }
    /// next edge-use in the edge-use ring of the cell
    Edge_use next_edge_use_in_cell(Edge_use u) const { return ucnext_[u]; }

    /// find the dart of half-face \a hf running from \a a to \a b, if any
    Halfedge find_halfedge(Halfface hf, Vertex a, Vertex b) const;

    /// find the geometric edge between \a a and \a b, if any
    Edge find_edge(Vertex a, Vertex b) const;

    /// find the geometric face carrying exactly the vertices \a verts, if any
    Face find_face(const std::vector<Vertex>& verts) const;

    /// find a free (boundary) half-face whose vertex loop is \a loop reversed,
    /// if any; this is the half-face that add_cell() would sew \a loop to
    Halfface find_free_halfface(const std::vector<Vertex>& loop) const;

    //@}

public:

    /// \name Predicates
    //@{

    /// a half-face is on the boundary iff it has no opposite
    bool is_boundary(Halfface hf) const
    { return !hmate_[Base::halfedge(hf)].is_valid(); }

    /// a dart is on the boundary iff it has no mate
    bool is_boundary(Halfedge h)  const { return !hmate_[h].is_valid(); }

    /// a geometric face is on the boundary iff it carries a single half-face
    bool is_boundary(Face f)      const { return is_boundary(halfface(f,0)); }

    /// a geometric edge is on the boundary iff one of its faces is
    bool is_boundary(Edge e)      const;

    /// a geometric vertex is on the boundary iff one of its faces is
    bool is_boundary(Vertex v)    const;

    /// a cell is on the boundary iff one of its half-faces is
    bool is_boundary(Cell c)      const;

    /// returns whether \c v is isolated, i.e. not used by any cell
    bool is_isolated(Vertex v) const { return !corner(v).is_valid(); }

    /// returns whether the cells around \c e form a single fan or a single
    /// cycle, i.e. whether the radial orbit reaches all of them
    bool is_manifold(Edge e) const;

    /// returns whether the cell star of \c v is connected through the faces at
    /// \c v, i.e. whether \c v is not a pinch point
    bool is_manifold(Vertex v) const;

    bool is_tet    (Cell c) const;
    bool is_hex    (Cell c) const;
    bool is_prism  (Cell c) const;
    bool is_pyramid(Cell c) const;

    /// returns whether every cell is a tetrahedron
    bool is_tet_mesh() const;
    /// returns whether every cell is a hexahedron
    bool is_hex_mesh() const;

    /** @brief check every topological invariant of the structure.

        Verifies that \f$\beta_3\f$ is a fixed-point-free involution on the
        non-boundary darts, that it reverses the half-face cycles and the dart
        orientations, that the two tiers agree on the edge, the face and the
        vertices of a dart and of its mate, that a geometric face carries at most
        two half-faces, that the four rings are consistent with the orbits they
        stand for, and that every cell boundary is a closed surface.

        \param msg if non-null, receives a description of the first violation
        \returns true if the structure is sound
    */
    bool is_valid_topology(std::string* msg = NULL) const;

    //@}

public:

    /// \name Counts and valences
    //@{

    // the base provides valence(Corner) and valence(Halfface)
    using Base::valence;

    /// the number of cells incident to vertex \c v
    unsigned int valence(Vertex v) const;
    /// the number of cells incident to edge \c e
    unsigned int valence(Edge e) const;
    /// the number of half-faces of face \c f: 1 on the boundary, 2 inside
    unsigned int valence(Face f) const { return is_boundary(f) ? 1u : 2u; }
    /// the number of half-faces of cell \c c
    unsigned int valence(Cell c) const;

    /// the number of vertices of cell \c c
    unsigned int n_vertices(Cell c) const;
    /// the number of edges of cell \c c
    unsigned int n_edges(Cell c) const;
    /// the number of faces of cell \c c
    unsigned int n_faces(Cell c) const { return valence(c); }

    /// the number of vertices of face \c f
    unsigned int n_vertices(Face f) const { return Base::valence(halfface(f,0)); }

    /** Mesh statistics.

        Prints the number of vertices, edges, faces, cells and darts, the
        boundary counts, the cell-type histogram and the extremal edge and
        vertex valences.
    */
    void mesh_statistics() const;

    //@}

public:

    /// \name Iterators and containers
    //@{

    Vertex_iterator vertices_begin() const { return Vertex_iterator(Vertex(0), this); }
    Vertex_iterator vertices_end()   const { return Vertex_iterator(Vertex(vertices_size()), this); }
    Edge_iterator   edges_begin()    const { return Edge_iterator(Edge(0), this); }
    Edge_iterator   edges_end()      const { return Edge_iterator(Edge(edges_size()), this); }
    Face_iterator   faces_begin()    const { return Face_iterator(Face(0), this); }
    Face_iterator   faces_end()      const { return Face_iterator(Face(faces_size()), this); }
    Cell_iterator   cells_begin()    const { return Cell_iterator(Cell(0), this); }
    Cell_iterator   cells_end()      const { return Cell_iterator(Cell(cells_size()), this); }

    /// returns vertex container for C++11 range-based for-loops
    Vertex_container vertices() const { return Vertex_container(vertices_begin(), vertices_end()); }
    /// returns edge container for C++11 range-based for-loops
    Edge_container   edges()    const { return Edge_container(edges_begin(), edges_end()); }
    /// returns face container for C++11 range-based for-loops
    Face_container   faces()    const { return Face_container(faces_begin(), faces_end()); }
    /// returns cell container for C++11 range-based for-loops
    Cell_container   cells()    const { return Cell_container(cells_begin(), cells_end()); }

    /// returns corner container for C++11 range-based for-loops
    Base::Vertex_container corners()   const { return Base::vertices(); }
    /// returns edge-use container for C++11 range-based for-loops
    Base::Edge_container   edge_uses() const { return Base::edges();    }
    /// returns half-face container for C++11 range-based for-loops
    Base::Face_container   halffaces() const { return Base::faces();    }

    // halfedges() and halfedges(Corner) / halfedges(Halfface) come from the base
    using Base::halfedges;

    //@}

public:

    /// \name Circulators
    //@{

    /// returns circulator for the half-faces of cell \c c
    Halfface_around_cell_circulator halffaces(Cell c) const
    { return Halfface_around_cell_circulator(this, c); }

    /// returns circulator for the corners of cell \c c
    Corner_around_cell_circulator corners(Cell c) const
    { return Corner_around_cell_circulator(this, c); }

    /// returns circulator for the edge-uses of cell \c c
    Edge_use_around_cell_circulator edge_uses(Cell c) const
    { return Edge_use_around_cell_circulator(this, c); }

    /// returns circulator for the cells adjacent to cell \c c
    Cell_around_cell_circulator cells(Cell c) const
    { return Cell_around_cell_circulator(this, c); }

    /// returns circulator for the corners of half-face \c hf
    Corner_around_halfface_circulator corners(Halfface hf) const
    { return Corner_around_halfface_circulator(this, hf); }

    /// returns circulator for the vertices of half-face \c hf
    Vertex_around_halfface_circulator vertices(Halfface hf) const
    { return Vertex_around_halfface_circulator(this, hf); }

    /// returns circulator for the darts on geometric edge \c e (radial cycle)
    Halfedge_around_edge_circulator halfedges(Edge e) const
    { return Halfedge_around_edge_circulator(this, e); }

    /// returns circulator for the half-faces around geometric edge \c e
    Halfface_around_edge_circulator halffaces(Edge e) const
    { return Halfface_around_edge_circulator(this, e); }

    /// returns circulator for the cells around geometric edge \c e
    Cell_around_edge_circulator cells(Edge e) const
    { return Cell_around_edge_circulator(this, e); }

    /// returns circulator for the corners at vertex \c v
    Corner_around_vertex_circulator corners(Vertex v) const
    { return Corner_around_vertex_circulator(this, v); }

    /// returns circulator for the cells at vertex \c v
    Cell_around_vertex_circulator cells(Vertex v) const
    { return Cell_around_vertex_circulator(this, v); }

    /// returns circulator for the darts emanating from vertex \c v
    Halfedge_around_vertex_circulator halfedges(Vertex v) const
    { return Halfedge_around_vertex_circulator(this, v); }

    //@}

public:

    /// \name Collected adjacencies
    ///
    /// These orbits are not cycles: an entity would be visited once per
    /// incident cell or per incident half-face.  They are therefore returned as
    /// de-duplicated vectors rather than as circulators, so that the cost is
    /// visible at the call site.
    //@{

    /// the vertices of cell \c c, each exactly once
    std::vector<Vertex> vertices(Cell c) const;
    /// the edges of cell \c c, each exactly once
    std::vector<Edge>   edges(Cell c)    const;
    /// the faces of cell \c c, each exactly once
    std::vector<Face>   faces(Cell c)    const;
    /// the edges incident to vertex \c v, each exactly once
    std::vector<Edge>   edges(Vertex v)  const;
    /// the faces incident to vertex \c v, each exactly once
    std::vector<Face>   faces(Vertex v)  const;
    /// the vertices of face \c f, in the order of its first half-face
    std::vector<Vertex> vertices(Face f) const;
    /// the two vertices of edge \c e
    std::vector<Vertex> vertices(Edge e) const;

    //@}

public:

    /// \name Deletion
    //@{

    /// deletes cell \c c together with the cell-local entities it owns; the
    /// geometric entities left without any use are deleted as well
    void delete_cell(Cell c);

    /// deletes the (isolated) vertex \c v
    void delete_vertex(Vertex v);

    //@}

public:

    /// \name Property handling
    //@{

    template <class T> Vertex_property<T> add_vertex_property(const std::string& name, const T t=T())
    { return Vertex_property<T>(vertprops_.add<T>(name, t)); }
    template <class T> Edge_property<T> add_edge_property(const std::string& name, const T t=T())
    { return Edge_property<T>(edgeprops_.add<T>(name, t)); }
    template <class T> Face_property<T> add_face_property(const std::string& name, const T t=T())
    { return Face_property<T>(faceprops_.add<T>(name, t)); }
    template <class T> Cell_property<T> add_cell_property(const std::string& name, const T t=T())
    { return Cell_property<T>(cellprops_.add<T>(name, t)); }

    template <class T> Corner_property<T> add_corner_property(const std::string& name, const T t=T())
    { return Base::add_vertex_property<T>(name, t); }
    template <class T> Edge_use_property<T> add_edge_use_property(const std::string& name, const T t=T())
    { return Base::add_edge_property<T>(name, t); }
    template <class T> Halfface_property<T> add_halfface_property(const std::string& name, const T t=T())
    { return Base::add_face_property<T>(name, t); }

    template <class T> Vertex_property<T> get_vertex_property(const std::string& name) const
    { return Vertex_property<T>(vertprops_.get<T>(name)); }
    template <class T> Edge_property<T> get_edge_property(const std::string& name) const
    { return Edge_property<T>(edgeprops_.get<T>(name)); }
    template <class T> Face_property<T> get_face_property(const std::string& name) const
    { return Face_property<T>(faceprops_.get<T>(name)); }
    template <class T> Cell_property<T> get_cell_property(const std::string& name) const
    { return Cell_property<T>(cellprops_.get<T>(name)); }

    template <class T> Corner_property<T> get_corner_property(const std::string& name) const
    { return Base::get_vertex_property<T>(name); }
    template <class T> Edge_use_property<T> get_edge_use_property(const std::string& name) const
    { return Base::get_edge_property<T>(name); }
    template <class T> Halfface_property<T> get_halfface_property(const std::string& name) const
    { return Base::get_face_property<T>(name); }

    template <class T> Vertex_property<T> vertex_property(const std::string& name, const T t=T())
    { return Vertex_property<T>(vertprops_.get_or_add<T>(name, t)); }
    template <class T> Edge_property<T> edge_property(const std::string& name, const T t=T())
    { return Edge_property<T>(edgeprops_.get_or_add<T>(name, t)); }
    template <class T> Face_property<T> face_property(const std::string& name, const T t=T())
    { return Face_property<T>(faceprops_.get_or_add<T>(name, t)); }
    template <class T> Cell_property<T> cell_property(const std::string& name, const T t=T())
    { return Cell_property<T>(cellprops_.get_or_add<T>(name, t)); }

    template <class T> void remove_vertex_property(Vertex_property<T>& p) { vertprops_.remove(p); }
    template <class T> void remove_edge_property  (Edge_property<T>&   p) { edgeprops_.remove(p); }
    template <class T> void remove_face_property  (Face_property<T>&   p) { faceprops_.remove(p); }
    template <class T> void remove_cell_property  (Cell_property<T>&   p) { cellprops_.remove(p); }

    /// prints the names of all properties of both tiers
    void property_stats() const;

    //@}

protected: //--------------------------------------------------- construction

    /// allocate a new geometric vertex
    Vertex new_vertex()
    {
        vertprops_.push_back();
        return Vertex(vertices_size()-1);
    }

    /// allocate a new geometric edge
    Edge new_edge()
    {
        edgeprops_.push_back();
        return Edge(edges_size()-1);
    }

    /// allocate a new geometric face
    Face new_face()
    {
        faceprops_.push_back();
        return Face(faces_size()-1);
    }

    /// allocate a new cell
    Cell new_cell()
    {
        cellprops_.push_back();
        return Cell(cells_size()-1);
    }

    /// allocate a corner as a use of \a v by cell \a c, and link it into the
    /// corner ring of \a v and the corner ring of \a c
    Corner new_corner(Vertex v, Cell c);

    /// register the (already created) edge-use \a u of cell \a c with the
    /// geometric edge \a e, and link it into the edge-use ring of \a c
    void register_edge_use(Edge_use u, Edge e, Cell c);

    /// register the (already created) half-face \a hf of cell \a c with the
    /// geometric face \a f, and link it into the half-face ring of \a c
    void register_halfface(Halfface hf, Face f, Cell c);

    /** set \f$\beta_3\f$ between the two half-faces \a hf0 and \a hf1.

        The two must carry the same vertex loop in opposite orientations, agree
        on the geometric edge of every dart, and both be free.  Only the darts
        are touched; merging the geometric faces is left to the caller, because
        add_cell() reuses the existing face instead of merging two.
    */
    void sew_darts(Halfface hf0, Halfface hf1);

    /// remap the handle-valued properties after the base has compacted its
    /// arrays; see gsMeshTopology::remap_handles
    virtual void remap_handles(const Base::Vertex_property<Base::Vertex>     & vmap,
                               const Base::Halfedge_property<Base::Halfedge> & hmap,
                               const Base::Face_property<Base::Face>         & fmap);

    /// bind the standard properties of both tiers to the (freshly copied)
    /// property containers
    void bind_properties();

    /// the face templates of the four standard cell types, as vertex-index
    /// loops oriented outwards; \a type is the number of vertices of the cell
    static const std::vector< std::vector<int> > & cell_type(int nverts);

protected:

    // ---- geometric tier ---------------------------------------------------
    gsProperty_container vertprops_;
    gsProperty_container edgeprops_;
    gsProperty_container faceprops_;
    gsProperty_container cellprops_;

    Vertex_property<Vertex_connectivity> vertconn_;
    Edge_property<Edge_connectivity>     edgeconn_;
    Face_property<Face_connectivity>     faceconn_;
    Cell_property<Cell_connectivity>     cellconn_;

    Vertex_property<bool> vertdeleted_;
    Edge_property<bool>   edgedeleted_;
    Face_property<bool>   facedeleted_;
    Cell_property<bool>   celldeleted_;

    // ---- cell-local tier: the extra connectivity of the 3-map -------------

    /// "h:mate" -- \f$\beta_3\f$, the only new topological relation
    Base::Halfedge_property<Halfedge> hmate_;

    /// "v:vertex"    -- corner -> its geometric vertex (the Cv2V table of the paper)
    Corner_property<Vertex>   cvertex_;
    /// "v:vertexnext" -- next corner in the corner ring of the geometric vertex
    Corner_property<Corner>   cvnext_;
    /// "v:cellnext"   -- next corner in the corner ring of the cell
    Corner_property<Corner>   ccnext_;

    /// "e:edge"     -- edge-use -> its geometric edge
    Edge_use_property<Edge>     uedge_;
    /// "e:cellnext" -- next edge-use in the edge-use ring of the cell
    Edge_use_property<Edge_use> ucnext_;

    /// "f:cell"     -- half-face -> its cell
    Halfface_property<Cell>     hfcell_;
    /// "f:face"     -- half-face -> its geometric face
    Halfface_property<Face>     hfface_;
    /// "f:cellnext" -- next half-face in the half-face ring of the cell
    Halfface_property<Halfface> hfnext_;

    unsigned int deleted_verts_;
    unsigned int deleted_geoedges_;
    unsigned int deleted_geofaces_;
    unsigned int deleted_cells_;
    bool         volgarbage_;

    /// scratch space for add_cell(), so that it does not allocate per cell
    std::vector<Vertex>   add_cell_verts_;
    std::vector<Halfface> add_cell_halffaces_;
};

/// print a geometric vertex handle
inline std::ostream& operator<<(std::ostream& os, gsVolMeshTopology::Vertex v)
{ return (os << 'V' << v.idx()); }

/// print a geometric edge handle
inline std::ostream& operator<<(std::ostream& os, gsVolMeshTopology::Edge e)
{ return (os << 'E' << e.idx()); }

/// print a geometric face handle
inline std::ostream& operator<<(std::ostream& os, gsVolMeshTopology::Face f)
{ return (os << 'F' << f.idx()); }

/// print a cell handle
inline std::ostream& operator<<(std::ostream& os, gsVolMeshTopology::Cell c)
{ return (os << 'C' << c.idx()); }

} // namespace gismo
