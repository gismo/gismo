/** @file gsVolMeshTopology.cpp

    @brief Implementation of gsVolMeshTopology, the non-templated half-face
    topology base of gsVolMesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#include <gsMesh2/gsVolMeshTopology.h>

#include <gsCore/gsDebug.h>

#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <utility>

namespace gismo {

// ===========================================================================
//  construction, assignment
// ===========================================================================

gsVolMeshTopology::
gsVolMeshTopology() : Base()
{
    // allocate the standard properties of both tiers; the same list is used in
    // bind_properties(), which is called after the containers are copied
    vertconn_    = add_vertex_property<Vertex_connectivity>("V:connectivity");
    edgeconn_    = add_edge_property  <Edge_connectivity>  ("E:connectivity");
    faceconn_    = add_face_property  <Face_connectivity>  ("F:connectivity");
    cellconn_    = add_cell_property  <Cell_connectivity>  ("C:connectivity");

    vertdeleted_ = add_vertex_property<bool>("V:deleted", false);
    edgedeleted_ = add_edge_property  <bool>("E:deleted", false);
    facedeleted_ = add_face_property  <bool>("F:deleted", false);
    celldeleted_ = add_cell_property  <bool>("C:deleted", false);

    hmate_       = Base::add_halfedge_property<Halfedge>("h:mate");

    cvertex_     = Base::add_vertex_property<Vertex>("v:vertex");
    cvnext_      = Base::add_vertex_property<Corner>("v:vertexnext");
    ccnext_      = Base::add_vertex_property<Corner>("v:cellnext");

    uedge_       = Base::add_edge_property<Edge>    ("e:edge");
    ucnext_      = Base::add_edge_property<Edge_use>("e:cellnext");

    hfcell_      = Base::add_face_property<Cell>    ("f:cell");
    hfface_      = Base::add_face_property<Face>    ("f:face");
    hfnext_      = Base::add_face_property<Halfface>("f:cellnext");

    deleted_verts_ = deleted_geoedges_ = deleted_geofaces_ = deleted_cells_ = 0;
    volgarbage_ = false;
}

gsVolMeshTopology::
~gsVolMeshTopology()
{
}

void
gsVolMeshTopology::
bind_properties()
{
    vertconn_    = vertex_property<Vertex_connectivity>("V:connectivity");
    edgeconn_    = edge_property  <Edge_connectivity>  ("E:connectivity");
    faceconn_    = face_property  <Face_connectivity>  ("F:connectivity");
    cellconn_    = cell_property  <Cell_connectivity>  ("C:connectivity");

    vertdeleted_ = vertex_property<bool>("V:deleted");
    edgedeleted_ = edge_property  <bool>("E:deleted");
    facedeleted_ = face_property  <bool>("F:deleted");
    celldeleted_ = cell_property  <bool>("C:deleted");

    hmate_       = Base::halfedge_property<Halfedge>("h:mate");

    cvertex_     = Base::vertex_property<Vertex>("v:vertex");
    cvnext_      = Base::vertex_property<Corner>("v:vertexnext");
    ccnext_      = Base::vertex_property<Corner>("v:cellnext");

    uedge_       = Base::edge_property<Edge>    ("e:edge");
    ucnext_      = Base::edge_property<Edge_use>("e:cellnext");

    hfcell_      = Base::face_property<Cell>    ("f:cell");
    hfface_      = Base::face_property<Face>    ("f:face");
    hfnext_      = Base::face_property<Halfface>("f:cellnext");
}

gsVolMeshTopology&
gsVolMeshTopology::
operator=(const gsVolMeshTopology& rhs)
{
    if (this != &rhs)
    {
        Base::operator=(rhs);

        vertprops_ = rhs.vertprops_;
        edgeprops_ = rhs.edgeprops_;
        faceprops_ = rhs.faceprops_;
        cellprops_ = rhs.cellprops_;

        // property handles hold pointers into the containers, so they have to
        // be looked up again in the copies
        bind_properties();

        deleted_verts_    = rhs.deleted_verts_;
        deleted_geoedges_ = rhs.deleted_geoedges_;
        deleted_geofaces_ = rhs.deleted_geofaces_;
        deleted_cells_    = rhs.deleted_cells_;
        volgarbage_       = rhs.volgarbage_;
    }
    return *this;
}

gsVolMeshTopology&
gsVolMeshTopology::
operator=(gsVolMeshTopology&& rhs) noexcept
{
    if (this != &rhs)
    {
        Base::operator=(std::move(rhs));

        vertprops_ = std::move(rhs.vertprops_);
        edgeprops_ = std::move(rhs.edgeprops_);
        faceprops_ = std::move(rhs.faceprops_);
        cellprops_ = std::move(rhs.cellprops_);

        bind_properties();

        deleted_verts_    = rhs.deleted_verts_;
        deleted_geoedges_ = rhs.deleted_geoedges_;
        deleted_geofaces_ = rhs.deleted_geofaces_;
        deleted_cells_    = rhs.deleted_cells_;
        volgarbage_       = rhs.volgarbage_;
    }
    return *this;
}

gsVolMeshTopology&
gsVolMeshTopology::
assign(const gsVolMeshTopology& rhs)
{
    if (this != &rhs)
    {
        Base::assign(rhs);

        vertprops_.clear();
        edgeprops_.clear();
        faceprops_.clear();
        cellprops_.clear();

        vertconn_    = add_vertex_property<Vertex_connectivity>("V:connectivity");
        edgeconn_    = add_edge_property  <Edge_connectivity>  ("E:connectivity");
        faceconn_    = add_face_property  <Face_connectivity>  ("F:connectivity");
        cellconn_    = add_cell_property  <Cell_connectivity>  ("C:connectivity");
        vertdeleted_ = add_vertex_property<bool>("V:deleted", false);
        edgedeleted_ = add_edge_property  <bool>("E:deleted", false);
        facedeleted_ = add_face_property  <bool>("F:deleted", false);
        celldeleted_ = add_cell_property  <bool>("C:deleted", false);

        vertconn_.array()    = rhs.vertconn_.array();
        edgeconn_.array()    = rhs.edgeconn_.array();
        faceconn_.array()    = rhs.faceconn_.array();
        cellconn_.array()    = rhs.cellconn_.array();
        vertdeleted_.array() = rhs.vertdeleted_.array();
        edgedeleted_.array() = rhs.edgedeleted_.array();
        facedeleted_.array() = rhs.facedeleted_.array();
        celldeleted_.array() = rhs.celldeleted_.array();

        vertprops_.resize(rhs.vertices_size());
        edgeprops_.resize(rhs.edges_size());
        faceprops_.resize(rhs.faces_size());
        cellprops_.resize(rhs.cells_size());

        // Base::assign() only carries the standard properties of the 2-map, so
        // the 3-map relations have to be re-created and copied here
        hmate_   = Base::add_halfedge_property<Halfedge>("h:mate");
        cvertex_ = Base::add_vertex_property<Vertex>("v:vertex");
        cvnext_  = Base::add_vertex_property<Corner>("v:vertexnext");
        ccnext_  = Base::add_vertex_property<Corner>("v:cellnext");
        uedge_   = Base::add_edge_property<Edge>    ("e:edge");
        ucnext_  = Base::add_edge_property<Edge_use>("e:cellnext");
        hfcell_  = Base::add_face_property<Cell>    ("f:cell");
        hfface_  = Base::add_face_property<Face>    ("f:face");
        hfnext_  = Base::add_face_property<Halfface>("f:cellnext");

        hmate_.array()   = rhs.hmate_.array();
        cvertex_.array() = rhs.cvertex_.array();
        cvnext_.array()  = rhs.cvnext_.array();
        ccnext_.array()  = rhs.ccnext_.array();
        uedge_.array()   = rhs.uedge_.array();
        ucnext_.array()  = rhs.ucnext_.array();
        hfcell_.array()  = rhs.hfcell_.array();
        hfface_.array()  = rhs.hfface_.array();
        hfnext_.array()  = rhs.hfnext_.array();

        deleted_verts_    = rhs.deleted_verts_;
        deleted_geoedges_ = rhs.deleted_geoedges_;
        deleted_geofaces_ = rhs.deleted_geofaces_;
        deleted_cells_    = rhs.deleted_cells_;
        volgarbage_       = rhs.volgarbage_;
    }
    return *this;
}

void
gsVolMeshTopology::
clear()
{
    Base::clear();

    vertprops_.resize(0);
    edgeprops_.resize(0);
    faceprops_.resize(0);
    cellprops_.resize(0);

    vertprops_.free_memory();
    edgeprops_.free_memory();
    faceprops_.free_memory();
    cellprops_.free_memory();

    deleted_verts_ = deleted_geoedges_ = deleted_geofaces_ = deleted_cells_ = 0;
    volgarbage_ = false;
}

void
gsVolMeshTopology::
free_memory()
{
    Base::free_memory();
    vertprops_.free_memory();
    edgeprops_.free_memory();
    faceprops_.free_memory();
    cellprops_.free_memory();
}

void
gsVolMeshTopology::
reserve(unsigned int nvertices, unsigned int nedges,
        unsigned int nfaces,    unsigned int ncells)
{
    // one corner per (cell,vertex) incidence and one dart per (half-face,edge)
    // incidence; the factors are the averages of a hexahedral mesh and only
    // steer the initial allocation
    Base::reserve(8*ncells, 12*ncells, 6*ncells);

    vertprops_.reserve(nvertices);
    edgeprops_.reserve(nedges);
    faceprops_.reserve(nfaces);
    cellprops_.reserve(ncells);
}

void
gsVolMeshTopology::
property_stats() const
{
    gsInfo << "cell-local (2-map) properties:\n";
    Base::property_stats();

    std::vector<std::string> props;

    gsInfo << "vertex properties:\n";
    props = vertprops_.properties();
    for (size_t i=0; i<props.size(); ++i) gsInfo << "\t" << props[i] << "\n";

    gsInfo << "edge properties:\n";
    props = edgeprops_.properties();
    for (size_t i=0; i<props.size(); ++i) gsInfo << "\t" << props[i] << "\n";

    gsInfo << "face properties:\n";
    props = faceprops_.properties();
    for (size_t i=0; i<props.size(); ++i) gsInfo << "\t" << props[i] << "\n";

    gsInfo << "cell properties:\n";
    props = cellprops_.properties();
    for (size_t i=0; i<props.size(); ++i) gsInfo << "\t" << props[i] << "\n";
}

// ===========================================================================
//  low level construction
// ===========================================================================

gsVolMeshTopology::Vertex
gsVolMeshTopology::
add_vertex()
{
    const Vertex v = new_vertex();
    vertconn_[v].corner_ = Corner();
    vertdeleted_[v] = false;
    return v;
}

void
gsVolMeshTopology::
add_batch_vertices(size_t nverts)
{
    for (size_t i=0; i!=nverts; ++i)
        add_vertex();
}

gsVolMeshTopology::Corner
gsVolMeshTopology::
new_corner(Vertex v, Cell c)
{
    const Corner cn = Base::new_vertex();

    Base::set_halfedge(cn, Halfedge());
    cvertex_[cn] = v;

    // link into the corner ring of the vertex
    const Corner vfirst = vertconn_[v].corner_;
    if (vfirst.is_valid())
    {
        cvnext_[cn] = cvnext_[vfirst];
        cvnext_[vfirst] = cn;
    }
    else
    {
        cvnext_[cn] = cn;
        vertconn_[v].corner_ = cn;
    }

    // link into the corner ring of the cell
    const Corner cfirst = cellconn_[c].corner_;
    if (cfirst.is_valid())
    {
        ccnext_[cn] = ccnext_[cfirst];
        ccnext_[cfirst] = cn;
    }
    else
    {
        ccnext_[cn] = cn;
        cellconn_[c].corner_ = cn;
    }

    return cn;
}

void
gsVolMeshTopology::
register_edge_use(Edge_use u, Edge e, Cell c)
{
    uedge_[u] = e;

    const Edge_use ufirst = cellconn_[c].edge_use_;
    if (ufirst.is_valid())
    {
        ucnext_[u] = ucnext_[ufirst];
        ucnext_[ufirst] = u;
    }
    else
    {
        ucnext_[u] = u;
        cellconn_[c].edge_use_ = u;
    }
}

void
gsVolMeshTopology::
register_halfface(Halfface hf, Face f, Cell c)
{
    hfcell_[hf] = c;
    hfface_[hf] = f;

    const Halfface hffirst = cellconn_[c].halfface_;
    if (hffirst.is_valid())
    {
        hfnext_[hf] = hfnext_[hffirst];
        hfnext_[hffirst] = hf;
    }
    else
    {
        hfnext_[hf] = hf;
        cellconn_[c].halfface_ = hf;
    }
}

// ===========================================================================
//  add_cell
// ===========================================================================

namespace {

// key of an unordered pair of indices
inline std::pair<int,int> upair(int a, int b)
{ return a<b ? std::make_pair(a,b) : std::make_pair(b,a); }

} // anonymous namespace

gsVolMeshTopology::Cell
gsVolMeshTopology::
add_cell(const std::vector< std::vector<Vertex> >& faces)
{
    // --- 1. validate the input -------------------------------------------
    //
    // The boundary of a cell has to be a closed, orientable, manifold surface.
    // That is exactly the statement that every directed vertex pair occurs in
    // the loops exactly once and that its reverse occurs exactly once.  Testing
    // it up front means the construction below cannot fail half way and leave
    // the mesh in a broken state.
    if (faces.size() < 4)
    {
        gsWarn << "gsVolMeshTopology::add_cell: a cell needs at least 4 faces, got "
               << faces.size() << ".\n";
        return Cell();
    }

    std::map<std::pair<int,int>, int> darts; // directed pair -> multiplicity
    add_cell_verts_.clear();
    for (size_t k=0; k!=faces.size(); ++k)
    {
        const std::vector<Vertex> & loop = faces[k];
        if (loop.size() < 3)
        {
            gsWarn << "gsVolMeshTopology::add_cell: face " << k
                   << " has fewer than 3 vertices.\n";
            return Cell();
        }
        for (size_t i=0; i!=loop.size(); ++i)
        {
            const Vertex a = loop[i], b = loop[(i+1)%loop.size()];
            if (!is_valid(a) || is_deleted(a))
            {
                gsWarn << "gsVolMeshTopology::add_cell: invalid vertex "<<a<<".\n";
                return Cell();
            }
            if (a == b)
            {
                gsWarn << "gsVolMeshTopology::add_cell: face " << k
                       << " repeats vertex " << a << ".\n";
                return Cell();
            }
            ++darts[std::make_pair(a.idx(), b.idx())];
            add_cell_verts_.push_back(a);
        }
    }

    for (std::map<std::pair<int,int>,int>::const_iterator it = darts.begin();
         it != darts.end(); ++it)
    {
        if (1 != it->second)
        {
            gsWarn << "gsVolMeshTopology::add_cell: the directed edge "
                   << Vertex(it->first.first) << "->" << Vertex(it->first.second)
                   << " occurs " << it->second << " times; the cell boundary is "
                   << "not a manifold surface.\n";
            return Cell();
        }
        const std::pair<int,int> rev(it->first.second, it->first.first);
        if (darts.end() == darts.find(rev))
        {
            gsWarn << "gsVolMeshTopology::add_cell: the edge "
                   << Vertex(it->first.first) << "-" << Vertex(it->first.second)
                   << " is used by a single face; the cell boundary is not closed "
                   << "or not consistently oriented.\n";
            return Cell();
        }
    }

    // distinct vertices of the cell
    std::sort(add_cell_verts_.begin(), add_cell_verts_.end());
    add_cell_verts_.erase(std::unique(add_cell_verts_.begin(), add_cell_verts_.end()),
                          add_cell_verts_.end());

    // Euler's formula for a genus-0 cell; darts.size() is twice the edge count
    const int nV = (int)add_cell_verts_.size();
    const int nE = (int)darts.size()/2;
    const int nF = (int)faces.size();
    if (2 != nV-nE+nF)
    {
        gsWarn << "gsVolMeshTopology::add_cell: the cell boundary has Euler "
               << "characteristic " << nV-nE+nF << " instead of 2.\n";
        return Cell();
    }

    // --- 2. the cell and its corners --------------------------------------
    const Cell c = new_cell();
    celldeleted_[c] = false;
    cellconn_[c].halfface_ = Halfface();
    cellconn_[c].corner_   = Corner();
    cellconn_[c].edge_use_ = Edge_use();

    std::map<int,Corner> corner_of; // vertex index -> corner of this cell
    for (size_t i=0; i!=add_cell_verts_.size(); ++i)
        corner_of[add_cell_verts_[i].idx()] = new_corner(add_cell_verts_[i], c);

    // --- 3. the half-faces, their darts and their edge-uses ---------------
    //
    // The corners are private to this cell and the input has been checked to be
    // a closed manifold surface, so every cell-local edge is used by exactly two
    // of the loops and the wiring below is unambiguous.  This is why add_cell
    // does not go through gsMeshTopology::add_face(): there is no boundary to
    // re-link and no complex vertex to guard against.
    std::map<std::pair<int,int>, Edge_use> use_of;
    add_cell_halffaces_.clear();
    std::vector<Halfedge> loop_darts;

    for (size_t k=0; k!=faces.size(); ++k)
    {
        const std::vector<Vertex> & loop = faces[k];
        const size_t n = loop.size();

        loop_darts.assign(n, Halfedge());
        for (size_t i=0; i!=n; ++i)
        {
            const Corner ca = corner_of[loop[i].idx()];
            const Corner cb = corner_of[loop[(i+1)%n].idx()];
            const std::pair<int,int> key = upair(ca.idx(), cb.idx());

            std::map<std::pair<int,int>, Edge_use>::iterator it = use_of.find(key);
            if (use_of.end() == it)
            {
                const Halfedge h = Base::new_edge(ca, cb); // h: ca -> cb
                use_of[key] = Base::edge(h);
                ucnext_[Base::edge(h)] = Edge_use();
                uedge_[Base::edge(h)]  = Edge();
                loop_darts[i] = h;
            }
            else
            {
                const Halfedge h0 = Base::halfedge(it->second, 0);
                loop_darts[i] = (Base::from_vertex(h0) == ca)
                              ? h0 : Base::opposite_halfedge(h0);
            }
            hmate_[loop_darts[i]] = Halfedge();
        }

        const Halfface hf = Base::new_face();
        hfnext_[hf] = Halfface();
        Base::set_halfedge(hf, loop_darts[0]);
        for (size_t i=0; i!=n; ++i)
        {
            Base::set_face(loop_darts[i], hf);
            Base::set_next_halfedge(loop_darts[i], loop_darts[(i+1)%n]);
            Base::set_halfedge(Base::from_vertex(loop_darts[i]), loop_darts[i]);
        }
        add_cell_halffaces_.push_back(hf);
    }

    // --- 4. attach the cell-local edges to the geometric edges ------------
    for (std::map<std::pair<int,int>, Edge_use>::const_iterator it = use_of.begin();
         it != use_of.end(); ++it)
    {
        const Edge_use u = it->second;
        const Halfedge h = Base::halfedge(u, 0);
        const Vertex a = cvertex_[Base::from_vertex(h)];
        const Vertex b = cvertex_[Base::to_vertex(h)];

        Edge e = find_edge(a, b);
        if (!e.is_valid())
        {
            e = new_edge();
            edgedeleted_[e] = false;
            edgeconn_[e].halfedge_ = h;
        }
        register_edge_use(u, e, c);
    }

    // --- 5. attach the half-faces to the geometric faces, gluing where the
    //        neighbouring cell is already there ---------------------------
    for (size_t k=0; k!=add_cell_halffaces_.size(); ++k)
        register_halfface(add_cell_halffaces_[k], Face(), c);

    for (size_t k=0; k!=add_cell_halffaces_.size(); ++k)
    {
        const Halfface hf = add_cell_halffaces_[k];
        const Halfface other = find_free_halfface(faces[k]);
        if (other.is_valid())
        {
            hfface_[hf] = hfface_[other];
            sew_darts(hf, other);
        }
        else
        {
            const Face f = new_face();
            facedeleted_[f] = false;
            faceconn_[f].halfface_ = hf;
            hfface_[hf] = f;
        }
    }

    return c;
}

// ---------------------------------------------------------------------------
//  the standard cell types, after Section 3.2.1 of Feng et al. (2013): a cell
//  type is described by its faces as local vertex loops oriented outwards.  The
//  local vertex numbering is that of VTK / CGNS.
// ---------------------------------------------------------------------------
const std::vector< std::vector<int> > &
gsVolMeshTopology::
cell_type(int nverts)
{
    struct Tables
    {
        std::vector< std::vector<int> > tet, pyr, pri, hex, none;

        static void push(std::vector< std::vector<int> > & to,
                         const int * f, int n)
        { to.push_back(std::vector<int>(f, f+n)); }

        Tables()
        {
            static const int t0[]={0,2,1}, t1[]={0,1,3}, t2[]={1,2,3}, t3[]={2,0,3};
            push(tet,t0,3); push(tet,t1,3); push(tet,t2,3); push(tet,t3,3);

            static const int y0[]={0,3,2,1}, y1[]={0,1,4}, y2[]={1,2,4},
                             y3[]={2,3,4},   y4[]={3,0,4};
            push(pyr,y0,4); push(pyr,y1,3); push(pyr,y2,3); push(pyr,y3,3);
            push(pyr,y4,3);

            static const int p0[]={0,2,1}, p1[]={3,4,5}, p2[]={0,1,4,3},
                             p3[]={1,2,5,4}, p4[]={2,0,3,5};
            push(pri,p0,3); push(pri,p1,3); push(pri,p2,4); push(pri,p3,4);
            push(pri,p4,4);

            static const int h0[]={0,3,2,1}, h1[]={4,5,6,7}, h2[]={0,1,5,4},
                             h3[]={1,2,6,5}, h4[]={2,3,7,6}, h5[]={3,0,4,7};
            push(hex,h0,4); push(hex,h1,4); push(hex,h2,4); push(hex,h3,4);
            push(hex,h4,4); push(hex,h5,4);
        }
    };

    static const Tables tables;

    switch (nverts)
    {
    case 4: return tables.tet;
    case 5: return tables.pyr;
    case 6: return tables.pri;
    case 8: return tables.hex;
    default: return tables.none;
    }
}

namespace {

inline std::vector< std::vector<gsVolMeshTopology::Vertex> >
expand_type(const std::vector< std::vector<int> > & type,
            const std::vector<gsVolMeshTopology::Vertex> & v)
{
    std::vector< std::vector<gsVolMeshTopology::Vertex> > faces(type.size());
    for (size_t i=0; i!=type.size(); ++i)
    {
        faces[i].reserve(type[i].size());
        for (size_t j=0; j!=type[i].size(); ++j)
            faces[i].push_back(v[type[i][j]]);
    }
    return faces;
}

} // anonymous namespace

gsVolMeshTopology::Cell
gsVolMeshTopology::
add_tet(Vertex v0, Vertex v1, Vertex v2, Vertex v3)
{
    std::vector<Vertex> v;
    v.push_back(v0); v.push_back(v1); v.push_back(v2); v.push_back(v3);
    return add_cell(expand_type(cell_type(4), v));
}

gsVolMeshTopology::Cell
gsVolMeshTopology::
add_pyramid(Vertex v0, Vertex v1, Vertex v2, Vertex v3, Vertex v4)
{
    std::vector<Vertex> v;
    v.push_back(v0); v.push_back(v1); v.push_back(v2); v.push_back(v3);
    v.push_back(v4);
    return add_cell(expand_type(cell_type(5), v));
}

gsVolMeshTopology::Cell
gsVolMeshTopology::
add_prism(Vertex v0, Vertex v1, Vertex v2, Vertex v3, Vertex v4, Vertex v5)
{
    std::vector<Vertex> v;
    v.push_back(v0); v.push_back(v1); v.push_back(v2);
    v.push_back(v3); v.push_back(v4); v.push_back(v5);
    return add_cell(expand_type(cell_type(6), v));
}

gsVolMeshTopology::Cell
gsVolMeshTopology::
add_hex(Vertex v0, Vertex v1, Vertex v2, Vertex v3,
        Vertex v4, Vertex v5, Vertex v6, Vertex v7)
{
    std::vector<Vertex> v;
    v.push_back(v0); v.push_back(v1); v.push_back(v2); v.push_back(v3);
    v.push_back(v4); v.push_back(v5); v.push_back(v6); v.push_back(v7);
    return add_cell(expand_type(cell_type(8), v));
}

// ===========================================================================
//  sewing
// ===========================================================================

void
gsVolMeshTopology::
sew_darts(Halfface hf0, Halfface hf1)
{
    const Halfedge h0 = Base::halfedge(hf0);
    const Vertex a = from_vertex(h0), b = to_vertex(h0);

    // hf1 carries the same loop reversed, so it holds the dart b->a
    const Halfedge h1 = find_halfedge(hf1, b, a);
    GISMO_ENSURE(h1.is_valid(),
                 "gsVolMeshTopology::sew: the two half-faces do not carry the "
                 "same vertex loop in opposite orientations.");

    Halfedge d0 = h0, d1 = h1;
    do
    {
        GISMO_ENSURE(from_vertex(d0) == to_vertex(d1) &&
                     to_vertex(d0)   == from_vertex(d1),
                     "gsVolMeshTopology::sew: the two vertex loops do not match.");
        GISMO_ENSURE(edge(d0) == edge(d1),
                     "gsVolMeshTopology::sew: the two half-faces disagree on the "
                     "geometric edge "<<from_vertex(d0)<<"-"<<to_vertex(d0)<<".");
        GISMO_ENSURE(!hmate_[d0].is_valid() && !hmate_[d1].is_valid(),
                     "gsVolMeshTopology::sew: half-face is already glued.");

        hmate_[d0] = d1;
        hmate_[d1] = d0;

        // beta3 reverses the face cycle: walking forwards on hf0 corresponds to
        // walking backwards on hf1, which is (beta1 o beta3)^2 = id
        d0 = Base::next_halfedge(d0);
        d1 = Base::prev_halfedge(d1);
    } while (d0 != h0);

    GISMO_ENSURE(d1 == h1,
                 "gsVolMeshTopology::sew: the two half-faces have different "
                 "numbers of vertices.");
}

void
gsVolMeshTopology::
sew(Halfface hf0, Halfface hf1)
{
    GISMO_ENSURE(is_boundary(hf0) && is_boundary(hf1),
                 "gsVolMeshTopology::sew: both half-faces have to be free.");
    GISMO_ENSURE(hf0 != hf1,
                 "gsVolMeshTopology::sew: cannot glue a half-face to itself.");

    const Face f0 = hfface_[hf0], f1 = hfface_[hf1];
    sew_darts(hf0, hf1);

    // merge the two geometric faces, keeping f0
    if (f0 != f1)
    {
        hfface_[hf1] = f0;
        faceconn_[f0].halfface_ = hf0;
        facedeleted_[f1] = true;
        ++deleted_geofaces_;
        volgarbage_ = true;
    }
}

void
gsVolMeshTopology::
unsew(Halfface hf)
{
    if (is_boundary(hf)) return;

    const Halfface other = opposite_halfface(hf);
    const Face     f     = hfface_[hf];

    Halfedge d = Base::halfedge(hf);
    do
    {
        const Halfedge m = hmate_[d];
        hmate_[d] = Halfedge();
        if (m.is_valid()) hmate_[m] = Halfedge();
        d = Base::next_halfedge(d);
    } while (d != Base::halfedge(hf));

    // the geometric face splits into two boundary faces
    faceconn_[f].halfface_ = hf;
    const Face nf = new_face();
    facedeleted_[nf] = false;
    faceconn_[nf].halfface_ = other;
    hfface_[other] = nf;
}

// ===========================================================================
//  searching
// ===========================================================================

gsVolMeshTopology::Halfedge
gsVolMeshTopology::
find_halfedge(Halfface hf, Vertex a, Vertex b) const
{
    const Halfedge start = Base::halfedge(hf);
    Halfedge h = start;
    do
    {
        if (from_vertex(h) == a && to_vertex(h) == b) return h;
        h = Base::next_halfedge(h);
    } while (h != start);
    return Halfedge();
}

gsVolMeshTopology::Edge
gsVolMeshTopology::
find_edge(Vertex a, Vertex b) const
{
    if (!a.is_valid() || !b.is_valid()) return Edge();

    for (auto cn : corners(a))
    {
        if (!Base::halfedge(cn).is_valid()) continue;
        for (auto h : Base::halfedges(cn))
        {
            if (cvertex_[Base::to_vertex(h)] != b) continue;
            const Edge e = uedge_[Base::edge(h)];
            // an edge-use created moments ago by add_cell() has no geometric
            // edge yet; it is not a match
            if (e.is_valid()) return e;
        }
    }
    return Edge();
}

gsVolMeshTopology::Halfface
gsVolMeshTopology::
find_free_halfface(const std::vector<Vertex>& loop) const
{
    if (loop.size() < 3) return Halfface();

    const Vertex a = loop[0], b = loop[1];

    for (auto cn : corners(a))
    {
        if (!Base::halfedge(cn).is_valid()) continue;
        for (auto hfit : Base::faces(cn))
        {
            const Halfface hf = hfit;
            if (!is_boundary(hf)) continue;
            if (Base::valence(hf) != (unsigned int)loop.size()) continue;

            // the candidate carries the loop reversed
            const Halfedge h = find_halfedge(hf, b, a);
            if (!h.is_valid()) continue;

            bool match = true;
            Halfedge d = h;
            for (size_t i=0; i!=loop.size(); ++i)
            {
                if (to_vertex(d) != loop[(loop.size()-i)%loop.size()])
                { match = false; break; }
                d = Base::next_halfedge(d);
            }
            if (match) return hf;
        }
    }
    return Halfface();
}

gsVolMeshTopology::Face
gsVolMeshTopology::
find_face(const std::vector<Vertex>& verts) const
{
    if (verts.empty()) return Face();

    std::vector<Vertex> want(verts);
    std::sort(want.begin(), want.end());

    for (auto cn : corners(verts[0]))
    {
        if (!Base::halfedge(cn).is_valid()) continue;
        for (auto hfit : Base::faces(cn))
        {
            const Halfface hf = hfit;
            if (Base::valence(hf) != (unsigned int)verts.size()) continue;

            std::vector<Vertex> have;
            for (auto v : vertices(hf)) have.push_back(v);
            std::sort(have.begin(), have.end());
            if (have == want) return hfface_[hf];
        }
    }
    return Face();
}

// ===========================================================================
//  predicates
// ===========================================================================

bool
gsVolMeshTopology::
is_boundary(Edge e) const
{
    const Halfedge start = halfedge(e);
    if (!start.is_valid()) return true;

    Halfedge h = start;
    do
    {
        const Halfedge n = radial_next(h);
        if (!n.is_valid()) return true;
        h = n;
    } while (h != start);

    return false;
}

bool
gsVolMeshTopology::
is_boundary(Vertex v) const
{
    for (auto cn : corners(v))
    {
        if (!Base::halfedge(cn).is_valid()) continue;
        for (auto hf : Base::faces(cn))
            if (is_boundary(hf)) return true;
    }
    return false;
}

bool
gsVolMeshTopology::
is_boundary(Cell c) const
{
    for (auto hf : halffaces(c))
        if (is_boundary(hf)) return true;
    return false;
}

bool
gsVolMeshTopology::
is_manifold(Edge e) const
{
    // the radial orbit reaches one cell per step; a non-manifold edge carries
    // several such fans, so compare with the total number of its edge-uses
    unsigned int reached = valence(e);

    unsigned int total = 0;
    const Vertex a = vertex(e,0), b = vertex(e,1);
    for (auto cn : corners(a))
    {
        if (!Base::halfedge(cn).is_valid()) continue;
        for (auto h : Base::halfedges(cn))
            if (cvertex_[Base::to_vertex(h)] == b && uedge_[Base::edge(h)] == e)
                ++total;
    }

    return reached == total;
}

bool
gsVolMeshTopology::
is_manifold(Vertex v) const
{
    // breadth-first search over the corners at v, stepping from a corner to the
    // corners of the cells that share a face with it at v
    std::vector<Corner> all;
    for (auto cn : corners(v)) all.push_back(cn);
    if (all.size() < 2) return true;

    std::set<int> seen;
    std::vector<Corner> stack;
    stack.push_back(all[0]);
    seen.insert(all[0].idx());

    while (!stack.empty())
    {
        const Corner cn = stack.back(); stack.pop_back();
        if (!Base::halfedge(cn).is_valid()) continue;

        for (auto h : Base::halfedges(cn))
        {
            const Halfedge m = hmate_[h];
            if (!m.is_valid()) continue;
            // beta3 keeps the dart on the same geometric vertex
            const Corner nb = Base::to_vertex(m);
            if (seen.insert(nb.idx()).second) stack.push_back(nb);

            const Corner nb2 = Base::from_vertex(m);
            if (cvertex_[nb2] == v && seen.insert(nb2.idx()).second)
                stack.push_back(nb2);
        }
    }

    return seen.size() == all.size();
}

bool gsVolMeshTopology::is_tet(Cell c) const
{ return 4 == valence(c) && 4 == n_vertices(c); }

bool gsVolMeshTopology::is_hex(Cell c) const
{ return 6 == valence(c) && 8 == n_vertices(c); }

bool gsVolMeshTopology::is_prism(Cell c) const
{ return 5 == valence(c) && 6 == n_vertices(c); }

bool gsVolMeshTopology::is_pyramid(Cell c) const
{ return 5 == valence(c) && 5 == n_vertices(c); }

bool
gsVolMeshTopology::
is_tet_mesh() const
{
    for (auto c : cells())
        if (!is_tet(c)) return false;
    return true;
}

bool
gsVolMeshTopology::
is_hex_mesh() const
{
    for (auto c : cells())
        if (!is_hex(c)) return false;
    return true;
}

// ===========================================================================
//  valences and counts
// ===========================================================================

unsigned int
gsVolMeshTopology::
valence(Vertex v) const
{
    unsigned int n = 0;
    for (auto cn : corners(v)) { ++n; GISMO_UNUSED(cn); }
    return n;
}

unsigned int
gsVolMeshTopology::
valence(Edge e) const
{
    unsigned int n = 0;
    for (auto h : halfedges(e)) { ++n; GISMO_UNUSED(h); }
    return n;
}

unsigned int
gsVolMeshTopology::
valence(Cell c) const
{
    unsigned int n = 0;
    for (auto hf : halffaces(c)) { ++n; GISMO_UNUSED(hf); }
    return n;
}

unsigned int
gsVolMeshTopology::
n_vertices(Cell c) const
{
    unsigned int n = 0;
    for (auto cn : corners(c)) { ++n; GISMO_UNUSED(cn); }
    return n;
}

unsigned int
gsVolMeshTopology::
n_edges(Cell c) const
{
    unsigned int n = 0;
    for (auto u : edge_uses(c)) { ++n; GISMO_UNUSED(u); }
    return n;
}

// ===========================================================================
//  canonical cell ordering
// ===========================================================================

int
gsVolMeshTopology::
cell_vtk_order(Cell c, std::vector<Vertex>& out) const
{
    out.clear();

    // classify the cell by its face structure; counting vertices and faces is
    // not enough, since e.g. a 6-face 8-vertex cell need not be a hexahedron
    unsigned int ntri = 0, nquad = 0, nother = 0;
    Halfface tri, quad;
    for (auto hf : halffaces(c))
    {
        const unsigned int k = Base::valence(hf);
        if      (3 == k) { ++ntri;  if (!tri.is_valid())  tri  = hf; }
        else if (4 == k) { ++nquad; if (!quad.is_valid()) quad = hf; }
        else             { ++nother; }
    }

    const unsigned int nv = n_vertices(c), nf = valence(c);

    const bool is_tet_ = (4==nv && 4==nf && 4==ntri && 0==nquad && 0==nother);
    const bool is_hex_ = (8==nv && 6==nf && 0==ntri && 6==nquad && 0==nother);
    const bool is_pri_ = (6==nv && 5==nf && 2==ntri && 3==nquad && 0==nother);
    const bool is_pyr_ = (5==nv && 5==nf && 4==ntri && 1==nquad && 0==nother);

    if (!is_tet_ && !is_hex_ && !is_pri_ && !is_pyr_)
    {
        for (auto cn : corners(c)) out.push_back(cvertex_[cn]);
        return 42;                                          // VTK_POLYHEDRON
    }

    // the vertices of a half-face, in its outward cyclic order
    std::vector<Vertex> loop;
    const Halfface seed = (is_pyr_ || is_hex_) ? quad : tri;
    for (auto v : vertices(seed)) loop.push_back(v);

    // vertex -> its corner in this cell, so that the neighbours of a vertex
    // inside the cell can be walked
    std::map<int,Corner> cor;
    for (auto cn : corners(c)) cor[cvertex_[cn].idx()] = cn;

    // the neighbour of \a v inside the cell that is not one of \a skip; for a
    // corner of a hexahedron, a prism cap or a pyramid base that is the one
    // edge leaving the seed face
    struct Local
    {
        static Vertex opposite(const gsVolMeshTopology& m,
                               const std::map<int,Corner>& cor, Vertex v,
                               const std::vector<Vertex>& skip)
        {
            const std::map<int,Corner>::const_iterator it = cor.find(v.idx());
            if (cor.end() == it) return Vertex();
            for (auto h : m.Base::halfedges(it->second))
            {
                const Vertex w = m.to_vertex(h);
                if (skip.end() == std::find(skip.begin(), skip.end(), w))
                    return w;
            }
            return Vertex();
        }
    };

    if (is_tet_)
    {
        // the outward loop of face 0 of the tetrahedron template is (v0,v2,v1)
        const Vertex v0 = loop[0], v2 = loop[1], v1 = loop[2];
        Vertex v3;
        for (auto cn : corners(c))
        {
            const Vertex w = cvertex_[cn];
            if (w!=v0 && w!=v1 && w!=v2) { v3 = w; break; }
        }
        out.push_back(v0); out.push_back(v1); out.push_back(v2); out.push_back(v3);
        return 10;                                               // VTK_TETRA
    }

    if (is_pyr_)
    {
        // the outward loop of the base of the pyramid template is (v0,v3,v2,v1)
        const Vertex v0 = loop[0], v3 = loop[1], v2 = loop[2], v1 = loop[3];
        Vertex v4;
        for (auto cn : corners(c))
        {
            const Vertex w = cvertex_[cn];
            if (w!=v0 && w!=v1 && w!=v2 && w!=v3) { v4 = w; break; }
        }
        out.push_back(v0); out.push_back(v1); out.push_back(v2);
        out.push_back(v3); out.push_back(v4);
        return 14;                                               // VTK_PYRAMID
    }

    if (is_pri_)
    {
        // the outward loop of the bottom triangle of the prism template is (v0,v2,v1)
        const Vertex v0 = loop[0], v2 = loop[1], v1 = loop[2];
        std::vector<Vertex> bottom;
        bottom.push_back(v0); bottom.push_back(v1); bottom.push_back(v2);
        out = bottom;
        out.push_back(Local::opposite(*this, cor, v0, bottom));
        out.push_back(Local::opposite(*this, cor, v1, bottom));
        out.push_back(Local::opposite(*this, cor, v2, bottom));
        return 13;                                               // VTK_WEDGE
    }

    // the outward loop of the bottom of the hexahedron template is (v0,v3,v2,v1)
    const Vertex v0 = loop[0], v3 = loop[1], v2 = loop[2], v1 = loop[3];
    std::vector<Vertex> bottom;
    bottom.push_back(v0); bottom.push_back(v1); bottom.push_back(v2); bottom.push_back(v3);
    out = bottom;
    out.push_back(Local::opposite(*this, cor, v0, bottom));
    out.push_back(Local::opposite(*this, cor, v1, bottom));
    out.push_back(Local::opposite(*this, cor, v2, bottom));
    out.push_back(Local::opposite(*this, cor, v3, bottom));
    return 12;                                                   // VTK_HEXAHEDRON
}

// ===========================================================================
//  collected adjacencies
// ===========================================================================

std::vector<gsVolMeshTopology::Vertex>
gsVolMeshTopology::
vertices(Cell c) const
{
    std::vector<Vertex> out;
    for (auto cn : corners(c)) out.push_back(cvertex_[cn]);
    return out;
}

std::vector<gsVolMeshTopology::Edge>
gsVolMeshTopology::
edges(Cell c) const
{
    std::vector<Edge> out;
    for (auto u : edge_uses(c)) out.push_back(uedge_[u]);
    return out;
}

std::vector<gsVolMeshTopology::Face>
gsVolMeshTopology::
faces(Cell c) const
{
    std::vector<Face> out;
    for (auto hf : halffaces(c)) out.push_back(hfface_[hf]);
    return out;
}

std::vector<gsVolMeshTopology::Edge>
gsVolMeshTopology::
edges(Vertex v) const
{
    std::vector<Edge> out;
    for (auto h : halfedges(v)) out.push_back(edge(h));
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

std::vector<gsVolMeshTopology::Face>
gsVolMeshTopology::
faces(Vertex v) const
{
    std::vector<Face> out;
    for (auto cn : corners(v))
    {
        if (!Base::halfedge(cn).is_valid()) continue;
        for (auto hf : Base::faces(cn)) out.push_back(hfface_[hf]);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

std::vector<gsVolMeshTopology::Vertex>
gsVolMeshTopology::
vertices(Face f) const
{
    std::vector<Vertex> out;
    for (auto v : vertices(halfface(f,0))) out.push_back(v);
    return out;
}

std::vector<gsVolMeshTopology::Vertex>
gsVolMeshTopology::
vertices(Edge e) const
{
    std::vector<Vertex> out;
    out.push_back(vertex(e,0));
    out.push_back(vertex(e,1));
    return out;
}

// ===========================================================================
//  the two-level vertex circulator
// ===========================================================================

gsVolMeshTopology::Halfedge_around_vertex_circulator::
Halfedge_around_vertex_circulator(const gsVolMeshTopology* m, Vertex v)
: mesh_(m), vertex_(v), active_(true)
{
    if (!mesh_ || !v.is_valid()) return;

    // settle on the first corner at v that actually carries a dart
    const Corner first = mesh_->corner(v);
    Corner cn = first;
    while (cn.is_valid())
    {
        const Halfedge h = mesh_->Base::halfedge(cn);
        if (h.is_valid())
        {
            corner_ = corner_start_ = cn;
            halfedge_ = corner_first_ = h;
            return;
        }
        cn = mesh_->next_corner_at_vertex(cn);
        if (cn == first) break;
    }
}

gsVolMeshTopology::Halfedge_around_vertex_circulator&
gsVolMeshTopology::Halfedge_around_vertex_circulator::
operator++()
{
    GISMO_ASSERT(mesh_ && halfedge_.is_valid(), "invalid circulator");
    active_ = true;

    // next outgoing dart within the star of the current corner
    halfedge_ = mesh_->Base::ccw_rotated_halfedge(halfedge_);
    if (halfedge_ != corner_first_) return *this;

    // that star is exhausted: step to the next corner carrying a dart
    Corner cn = mesh_->next_corner_at_vertex(corner_);
    while (cn != corner_start_)
    {
        const Halfedge h = mesh_->Base::halfedge(cn);
        if (h.is_valid())
        {
            corner_   = cn;
            halfedge_ = corner_first_ = h;
            return *this;
        }
        cn = mesh_->next_corner_at_vertex(cn);
    }

    // back where we started: the orbit is complete
    corner_   = corner_start_;
    halfedge_ = corner_first_ = mesh_->Base::halfedge(corner_start_);
    return *this;
}

// ===========================================================================
//  validation
// ===========================================================================

bool
gsVolMeshTopology::
is_valid_topology(std::string* msg) const
{
    std::ostringstream os;

#   define GSVOL_FAIL(cond, what)                                       \
    if (!(cond)) { os << what; if (msg) *msg = os.str(); return false; }

    // --- beta3 -----------------------------------------------------------
    for (auto h : Base::halfedges())
    {
        if (Base::is_deleted(h)) continue;

        const Halfface hf = Base::face(h);
        GSVOL_FAIL(hf.is_valid(),
                   "dart "<<h<<" has no half-face: the cell shell is not closed");

        const Halfedge m = hmate_[h];
        if (!m.is_valid()) continue;                    // boundary dart

        GSVOL_FAIL(m != h, "dart "<<h<<" is its own mate");
        GSVOL_FAIL(hmate_[m] == h,
                   "beta3 is not an involution at dart "<<h);
        GSVOL_FAIL(Base::face(m) != hf,
                   "dart "<<h<<" is glued inside its own half-face");
        GSVOL_FAIL(hmate_[Base::next_halfedge(h)] == Base::prev_halfedge(m),
                   "beta3 does not reverse the face cycle at dart "<<h
                   <<": (beta1 o beta3)^2 != id");
        GSVOL_FAIL(to_vertex(m) == from_vertex(h) && from_vertex(m) == to_vertex(h),
                   "dart "<<h<<" and its mate disagree on the vertices "
                   <<from_vertex(h)<<"-"<<to_vertex(h));
        GSVOL_FAIL(edge(m) == edge(h),
                   "dart "<<h<<" and its mate disagree on the geometric edge");
        GSVOL_FAIL(hfface_[Base::face(m)] == hfface_[hf],
                   "half-faces "<<hf.idx()<<" and "<<Base::face(m).idx()
                   <<" are glued but belong to different geometric faces");
    }

    // --- cells and their rings -------------------------------------------
    for (auto c : cells())
    {
        unsigned int nhf = 0, ncn = 0, nu = 0;

        for (auto hf : halffaces(c))
        {
            ++nhf;
            GSVOL_FAIL(hfcell_[hf] == c,
                       "half-face "<<hf.idx()<<" is in the ring of cell "<<c
                       <<" but reports cell "<<hfcell_[hf]);
            GSVOL_FAIL(hfface_[hf].is_valid(),
                       "half-face "<<hf.idx()<<" has no geometric face");
            GSVOL_FAIL(nhf <= halffaces_size(), "the half-face ring of cell "<<c<<" does not close");
        }
        GSVOL_FAIL(nhf >= 4, "cell "<<c<<" has only "<<nhf<<" faces");

        for (auto cn : corners(c))
        {
            ++ncn;
            GSVOL_FAIL(cell(cn) == c,
                       "corner "<<cn.idx()<<" is in the ring of cell "<<c
                       <<" but reports cell "<<cell(cn));
            GSVOL_FAIL(ncn <= corners_size(), "the corner ring of cell "<<c<<" does not close");
        }

        for (auto u : edge_uses(c))
        {
            ++nu;
            GSVOL_FAIL(uedge_[u].is_valid(),
                       "edge-use "<<u.idx()<<" of cell "<<c<<" has no geometric edge");
            GSVOL_FAIL(nu <= edge_uses_size(), "the edge-use ring of cell "<<c<<" does not close");
        }

        // Euler characteristic of the cell boundary
        GSVOL_FAIL(2 == (int)ncn - (int)nu + (int)nhf,
                   "the boundary of cell "<<c<<" has Euler characteristic "
                   <<(int)ncn-(int)nu+(int)nhf<<" instead of 2");
    }

    // --- geometric faces --------------------------------------------------
    for (auto f : faces())
    {
        const Halfface hf0 = halfface(f,0);
        GSVOL_FAIL(hf0.is_valid(), "face "<<f<<" has no half-face");
        GSVOL_FAIL(hfface_[hf0] == f,
                   "face "<<f<<" points at half-face "<<hf0.idx()
                   <<", which belongs to face "<<hfface_[hf0]);

        const Halfface hf1 = halfface(f,1);
        if (hf1.is_valid())
            GSVOL_FAIL(hfface_[hf1] == f,
                       "face "<<f<<" and its second half-face disagree");
    }

    // --- geometric vertices and their corner rings ------------------------
    for (auto v : vertices())
    {
        unsigned int n = 0;
        for (auto cn : corners(v))
        {
            ++n;
            GSVOL_FAIL(cvertex_[cn] == v,
                       "corner "<<cn.idx()<<" is in the ring of vertex "<<v
                       <<" but reports vertex "<<cvertex_[cn]);
            GSVOL_FAIL(n <= corners_size(), "the corner ring of vertex "<<v<<" does not close");
        }
    }

    // --- geometric edges and their radial cycles --------------------------
    for (auto e : edges())
    {
        const Halfedge h = halfedge(e);
        GSVOL_FAIL(h.is_valid(), "edge "<<e<<" has no dart");
        GSVOL_FAIL(uedge_[Base::edge(h)] == e,
                   "edge "<<e<<" points at a dart whose edge-use belongs to "
                   <<uedge_[Base::edge(h)]);

        unsigned int n = 0;
        for (auto d : halfedges(e))
        {
            ++n;
            GSVOL_FAIL(uedge_[Base::edge(d)] == e,
                       "the radial cycle of edge "<<e<<" leaves the edge");
            GSVOL_FAIL(n <= halfedges_size(), "the radial cycle of edge "<<e<<" does not close");
        }
        GSVOL_FAIL(n > 0, "edge "<<e<<" has an empty radial cycle");
    }

#   undef GSVOL_FAIL

    if (msg) msg->clear();
    return true;
}

// ===========================================================================
//  statistics
// ===========================================================================

void
gsVolMeshTopology::
mesh_statistics() const
{
    unsigned int ntet=0, nhex=0, npri=0, npyr=0, nother=0, nbcell=0;
    for (auto c : cells())
    {
        if      (is_tet(c))     ++ntet;
        else if (is_hex(c))     ++nhex;
        else if (is_prism(c))   ++npri;
        else if (is_pyramid(c)) ++npyr;
        else                    ++nother;
        if (is_boundary(c))     ++nbcell;
    }

    unsigned int nbface = 0;
    for (auto f : faces())
        if (is_boundary(f)) ++nbface;

    unsigned int emin = 0, emax = 0;
    bool first = true;
    for (auto e : edges())
    {
        const unsigned int k = valence(e);
        if (first) { emin = emax = k; first = false; }
        else { emin = std::min(emin,k); emax = std::max(emax,k); }
    }

    unsigned int vmin = 0, vmax = 0;
    first = true;
    for (auto v : vertices())
    {
        if (is_isolated(v)) continue;
        const unsigned int k = valence(v);
        if (first) { vmin = vmax = k; first = false; }
        else { vmin = std::min(vmin,k); vmax = std::max(vmax,k); }
    }

    gsInfo << "gsVolMeshTopology statistics\n"
           << "  vertices  : " << n_vertices()  << "\n"
           << "  edges     : " << n_edges()     << "\n"
           << "  faces     : " << n_faces()     << " (" << nbface << " on the boundary)\n"
           << "  cells     : " << n_cells()     << " (" << nbcell << " touching the boundary)\n"
           << "  corners   : " << n_corners()   << "\n"
           << "  edge-uses : " << n_edge_uses() << "\n"
           << "  half-faces: " << n_halffaces() << "\n"
           << "  darts     : " << Base::n_halfedges() << "\n"
           << "  cell types: " << ntet << " tet, " << nhex << " hex, "
           << npri << " prism, " << npyr << " pyramid, " << nother << " other\n"
           << "  edge valence   (cells per edge)  : [" << emin << ", " << emax << "]\n"
           << "  vertex valence (cells per vertex): [" << vmin << ", " << vmax << "]\n";
}

// ===========================================================================
//  deletion and garbage collection
// ===========================================================================

void
gsVolMeshTopology::
delete_cell(Cell c)
{
    if (celldeleted_[c]) return;

    std::vector<Edge_use> us;
    for (auto u : edge_uses(c)) us.push_back(u);

    // 1. a surviving geometric edge must not keep its representative dart
    //    inside the cell that is about to disappear.  This has to happen while
    //    beta3 is still intact, i.e. before the unsewing below.
    for (size_t i=0; i!=us.size(); ++i)
    {
        const Edge e = uedge_[us[i]];
        if (!e.is_valid() || edgedeleted_[e]) continue;
        if (hfcell_[Base::face(edgeconn_[e].halfedge_)] != c) continue;

        Halfedge keep;
        for (auto d : halfedges(e))
            if (hfcell_[Base::face(d)] != c) { keep = d; break; }
        if (keep.is_valid()) edgeconn_[e].halfedge_ = keep;
    }

    // 2. release the neighbours
    std::vector<Halfface> hfs;
    for (auto hf : halffaces(c)) hfs.push_back(hf);
    for (size_t i=0; i!=hfs.size(); ++i)
        if (!is_boundary(hfs[i])) unsew(hfs[i]);

    // 3. drop the cell-local entities, all of which belong to this cell alone
    for (size_t i=0; i!=hfs.size(); ++i)
    {
        Base::fdeleted_[hfs[i]] = true;
        ++Base::deleted_faces_;
        facedeleted_[hfface_[hfs[i]]] = true;
        ++deleted_geofaces_;
    }

    for (size_t i=0; i!=us.size(); ++i)
    {
        Base::edeleted_[us[i]] = true;
        ++Base::deleted_edges_;
    }

    std::vector<Corner> cns;
    for (auto cn : corners(c)) cns.push_back(cn);
    for (size_t i=0; i!=cns.size(); ++i)
    {
        // unlink from the corner ring of the geometric vertex
        const Vertex v = cvertex_[cns[i]];
        Corner p = cvnext_[cns[i]];
        while (cvnext_[p] != cns[i]) p = cvnext_[p];
        if (p == cns[i]) vertconn_[v].corner_ = Corner();     // last corner
        else
        {
            cvnext_[p] = cvnext_[cns[i]];
            if (vertconn_[v].corner_ == cns[i]) vertconn_[v].corner_ = p;
        }

        Base::vdeleted_[cns[i]] = true;
        ++Base::deleted_vertices_;
    }

    // 4. geometric edges left without any use
    for (size_t i=0; i!=us.size(); ++i)
    {
        const Edge e = uedge_[us[i]];
        if (!e.is_valid() || edgedeleted_[e]) continue;
        const Vertex a = vertex(e,0), b = vertex(e,1);
        if (!find_edge(a,b).is_valid())
        {
            edgedeleted_[e] = true;
            ++deleted_geoedges_;
        }
    }

    celldeleted_[c] = true;
    ++deleted_cells_;
    volgarbage_ = true;
    Base::garbage_ = true;
}

void
gsVolMeshTopology::
delete_vertex(Vertex v)
{
    if (vertdeleted_[v]) return;
    GISMO_ENSURE(is_isolated(v),
                 "gsVolMeshTopology::delete_vertex: vertex "<<v<<" is still used "
                 "by "<<valence(v)<<" cells; delete those first.");
    vertdeleted_[v] = true;
    ++deleted_verts_;
    volgarbage_ = true;
}

void
gsVolMeshTopology::
remap_handles(const Base::Vertex_property<Base::Vertex>     & vmap,
              const Base::Halfedge_property<Base::Halfedge> & hmap,
              const Base::Face_property<Base::Face>         & fmap)
{
    // The base has compacted its arrays but not yet resized them; the surviving
    // entities occupy the first n_* slots and the maps translate old indices
    // into new ones.
    const unsigned int nCn = Base::n_vertices();
    const unsigned int nH  = Base::n_halfedges();
    const unsigned int nU  = Base::n_edges();
    const unsigned int nHf = Base::n_faces();

    for (unsigned int i=0; i!=nH; ++i)
    {
        const Halfedge h(i);
        if (hmate_[h].is_valid()) hmate_[h] = hmap[hmate_[h]];
    }

    for (unsigned int i=0; i!=nCn; ++i)
    {
        const Corner cn(i);
        if (cvnext_[cn].is_valid()) cvnext_[cn] = vmap[cvnext_[cn]];
        if (ccnext_[cn].is_valid()) ccnext_[cn] = vmap[ccnext_[cn]];
    }

    for (unsigned int i=0; i!=nU; ++i)
    {
        const Edge_use u(i);
        if (ucnext_[u].is_valid())
            ucnext_[u] = Edge_use(hmap[Base::halfedge(ucnext_[u],0)].idx()>>1);
    }

    for (unsigned int i=0; i!=nHf; ++i)
    {
        const Halfface hf(i);
        if (hfnext_[hf].is_valid()) hfnext_[hf] = fmap[hfnext_[hf]];
    }

    // the geometric tier points back into the cell-local tier
    for (auto v : vertices())
        if (vertconn_[v].corner_.is_valid())
            vertconn_[v].corner_ = vmap[vertconn_[v].corner_];

    for (auto e : edges())
        if (edgeconn_[e].halfedge_.is_valid())
            edgeconn_[e].halfedge_ = hmap[edgeconn_[e].halfedge_];

    for (auto f : faces())
        if (faceconn_[f].halfface_.is_valid())
            faceconn_[f].halfface_ = fmap[faceconn_[f].halfface_];

    for (auto c : cells())
    {
        Cell_connectivity & cc = cellconn_[c];
        if (cc.halfface_.is_valid()) cc.halfface_ = fmap[cc.halfface_];
        if (cc.corner_.is_valid())   cc.corner_   = vmap[cc.corner_];
        if (cc.edge_use_.is_valid())
            cc.edge_use_ = Edge_use(hmap[Base::halfedge(cc.edge_use_,0)].idx()>>1);
    }
}

void
gsVolMeshTopology::
garbage_collection()
{
    if (!garbage()) return;

    // --- 1. compact the geometric tier -----------------------------------
    //
    // The scheme is the one of gsMeshTopology::garbage_collection(): rows are
    // brought to the front by disjoint transpositions, so the resulting
    // permutation is an involution and one map serves both directions.
    Vertex_property<Vertex> Vmap = add_vertex_property<Vertex>("V:garbage-collection");
    Edge_property<Edge>     Emap = add_edge_property<Edge>    ("E:garbage-collection");
    Face_property<Face>     Fmap = add_face_property<Face>    ("F:garbage-collection");
    Cell_property<Cell>     Cmap = add_cell_property<Cell>    ("C:garbage-collection");

    int nV = (int)vertices_size(), nE = (int)edges_size();
    int nF = (int)faces_size(),    nC = (int)cells_size();

    for (int i=0; i<nV; ++i) Vmap[Vertex(i)] = Vertex(i);
    for (int i=0; i<nE; ++i) Emap[Edge(i)]   = Edge(i);
    for (int i=0; i<nF; ++i) Fmap[Face(i)]   = Face(i);
    for (int i=0; i<nC; ++i) Cmap[Cell(i)]   = Cell(i);

    if (nV > 0)
    {
        int i0=0, i1=nV-1;
        while (true)
        {
            while (!vertdeleted_[Vertex(i0)] && i0 < i1) ++i0;
            while ( vertdeleted_[Vertex(i1)] && i0 < i1) --i1;
            if (i0 >= i1) break;
            vertprops_.swap(i0, i1);
        }
        nV = vertdeleted_[Vertex(i0)] ? i0 : i0+1;
    }
    if (nE > 0)
    {
        int i0=0, i1=nE-1;
        while (true)
        {
            while (!edgedeleted_[Edge(i0)] && i0 < i1) ++i0;
            while ( edgedeleted_[Edge(i1)] && i0 < i1) --i1;
            if (i0 >= i1) break;
            edgeprops_.swap(i0, i1);
        }
        nE = edgedeleted_[Edge(i0)] ? i0 : i0+1;
    }
    if (nF > 0)
    {
        int i0=0, i1=nF-1;
        while (true)
        {
            while (!facedeleted_[Face(i0)] && i0 < i1) ++i0;
            while ( facedeleted_[Face(i1)] && i0 < i1) --i1;
            if (i0 >= i1) break;
            faceprops_.swap(i0, i1);
        }
        nF = facedeleted_[Face(i0)] ? i0 : i0+1;
    }
    if (nC > 0)
    {
        int i0=0, i1=nC-1;
        while (true)
        {
            while (!celldeleted_[Cell(i0)] && i0 < i1) ++i0;
            while ( celldeleted_[Cell(i1)] && i0 < i1) --i1;
            if (i0 >= i1) break;
            cellprops_.swap(i0, i1);
        }
        nC = celldeleted_[Cell(i0)] ? i0 : i0+1;
    }

    // --- 2. the cell-local tier points at the geometric one ---------------
    for (unsigned int i=0; i!=Base::vertices_size(); ++i)
    {
        const Corner cn(i);
        if (!Base::is_deleted(cn) && cvertex_[cn].is_valid())
            cvertex_[cn] = Vmap[cvertex_[cn]];
    }
    for (unsigned int i=0; i!=Base::edges_size(); ++i)
    {
        const Edge_use u(i);
        if (!Base::is_deleted(u) && uedge_[u].is_valid())
            uedge_[u] = Emap[uedge_[u]];
    }
    for (unsigned int i=0; i!=Base::faces_size(); ++i)
    {
        const Halfface hf(i);
        if (Base::is_deleted(hf)) continue;
        if (hfface_[hf].is_valid()) hfface_[hf] = Fmap[hfface_[hf]];
        if (hfcell_[hf].is_valid()) hfcell_[hf] = Cmap[hfcell_[hf]];
    }

    remove_vertex_property(Vmap);
    remove_edge_property(Emap);
    remove_face_property(Fmap);
    remove_cell_property(Cmap);

    vertprops_.resize(nV); vertprops_.free_memory();
    edgeprops_.resize(nE); edgeprops_.free_memory();
    faceprops_.resize(nF); faceprops_.free_memory();
    cellprops_.resize(nC); cellprops_.free_memory();

    deleted_verts_ = deleted_geoedges_ = deleted_geofaces_ = deleted_cells_ = 0;
    volgarbage_ = false;

    // --- 3. compact the cell-local tier; remap_handles() fixes the rest ----
    Base::garbage_collection();
}

} // namespace gismo
