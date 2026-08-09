/** @file gsMeshTopology.cpp

    @brief Implementation of gsMeshTopology, the combinatorial core shared by
    gsSurfMeshTopology and gsVolMeshTopology.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#include <gsMesh2/gsMeshTopology.h>

#include <gsCore/gsDebug.h>

#include <algorithm>
#include <map>

namespace gismo {

gsMeshTopology::
gsMeshTopology()
{
    // allocate standard properties
    // same list is used in operator=() and assign()
    vconn_    = add_vertex_property<Vertex_connectivity>("v:connectivity");
    hconn_    = add_halfedge_property<Halfedge_connectivity>("h:connectivity");
    fconn_    = add_face_property<Face_connectivity>("f:connectivity");
    vdeleted_ = add_vertex_property<bool>("v:deleted", false);
    edeleted_ = add_edge_property<bool>("e:deleted", false);
    fdeleted_ = add_face_property<bool>("f:deleted", false);

    mprops_.push_back();

    deleted_vertices_ = deleted_edges_ = deleted_faces_ = 0;
    garbage_ = false;
}

gsMeshTopology::
~gsMeshTopology()
{
}

gsMeshTopology&
gsMeshTopology::
operator=(const gsMeshTopology& rhs)
{
    if (this != &rhs)
    {
        // deep copy of property containers
        vprops_ = rhs.vprops_;
        hprops_ = rhs.hprops_;
        eprops_ = rhs.eprops_;
        fprops_ = rhs.fprops_;
        mprops_ = rhs.mprops_;

        // property handles contain pointers, have to be reassigned
        vconn_    = vertex_property<Vertex_connectivity>("v:connectivity");
        hconn_    = halfedge_property<Halfedge_connectivity>("h:connectivity");
        fconn_    = face_property<Face_connectivity>("f:connectivity");
        vdeleted_ = vertex_property<bool>("v:deleted");
        edeleted_ = edge_property<bool>("e:deleted");
        fdeleted_ = face_property<bool>("f:deleted");

        // how many elements are deleted?
        deleted_vertices_ = rhs.deleted_vertices_;
        deleted_edges_    = rhs.deleted_edges_;
        deleted_faces_    = rhs.deleted_faces_;
        garbage_          = rhs.garbage_;
    }

    return *this;
}

gsMeshTopology&
gsMeshTopology::
assign(const gsMeshTopology& rhs)
{
    if (this != &rhs)
    {
        // clear properties
        vprops_.clear();
        hprops_.clear();
        eprops_.clear();
        fprops_.clear();
        mprops_.clear();

        // allocate standard properties
        vconn_    = add_vertex_property<Vertex_connectivity>("v:connectivity");
        hconn_    = add_halfedge_property<Halfedge_connectivity>("h:connectivity");
        fconn_    = add_face_property<Face_connectivity>("f:connectivity");
        vdeleted_ = add_vertex_property<bool>("v:deleted", false);
        edeleted_ = add_edge_property<bool>("e:deleted", false);
        fdeleted_ = add_face_property<bool>("f:deleted", false);

        // copy properties from other mesh
        vconn_.array()     = rhs.vconn_.array();
        hconn_.array()     = rhs.hconn_.array();
        fconn_.array()     = rhs.fconn_.array();
        vdeleted_.array()  = rhs.vdeleted_.array();
        edeleted_.array()  = rhs.edeleted_.array();
        fdeleted_.array()  = rhs.fdeleted_.array();

        // resize (needed by property containers)
        vprops_.resize(rhs.vertices_size());
        hprops_.resize(rhs.halfedges_size());
        eprops_.resize(rhs.edges_size());
        fprops_.resize(rhs.faces_size());
        mprops_.resize(1);

        // how many elements are deleted?
        deleted_vertices_ = rhs.deleted_vertices_;
        deleted_edges_    = rhs.deleted_edges_;
        deleted_faces_    = rhs.deleted_faces_;
        garbage_          = rhs.garbage_;
    }

    return *this;
}

gsMeshTopology&
gsMeshTopology::
operator=(gsMeshTopology&& rhs) noexcept
{
    if (this != &rhs)
    {
        // deep copy of property containers
        vprops_ = rhs.vprops_;
        hprops_ = rhs.hprops_;
        eprops_ = rhs.eprops_;
        fprops_ = rhs.fprops_;
        mprops_ = rhs.mprops_;

        // property handles contain pointers, have to be reassigned
        vconn_    = vertex_property<Vertex_connectivity>("v:connectivity");
        hconn_    = halfedge_property<Halfedge_connectivity>("h:connectivity");
        fconn_    = face_property<Face_connectivity>("f:connectivity");
        vdeleted_ = vertex_property<bool>("v:deleted");
        edeleted_ = edge_property<bool>("e:deleted");
        fdeleted_ = face_property<bool>("f:deleted");

        // how many elements are deleted?
        deleted_vertices_ = rhs.deleted_vertices_;
        deleted_edges_    = rhs.deleted_edges_;
        deleted_faces_    = rhs.deleted_faces_;
        garbage_          = rhs.garbage_;

        rhs.clear();
    }
    return *this;
}


void
gsMeshTopology::
clear()
{
    vprops_.resize(0);
    hprops_.resize(0);
    eprops_.resize(0);
    fprops_.resize(0);
    mprops_.resize(0);

    free_memory();

    deleted_vertices_ = deleted_edges_ = deleted_faces_ = 0;
    garbage_ = false;
}

void
gsMeshTopology::
free_memory()
{
    vprops_.free_memory();
    hprops_.free_memory();
    eprops_.free_memory();
    fprops_.free_memory();
    mprops_.free_memory();
}

void
gsMeshTopology::
reserve(unsigned int nvertices,
        unsigned int nedges,
        unsigned int nfaces )
{
    vprops_.reserve(nvertices);
    hprops_.reserve(2*nedges);
    eprops_.reserve(nedges);
    fprops_.reserve(nfaces);
    mprops_.reserve(1);
}

void
gsMeshTopology::
property_stats() const
{
    std::vector<std::string> props;

    std::cout << "vertex properties:\n";
    props = vertex_properties();
    for (unsigned int i=0; i<props.size(); ++i)
        std::cout << "\t" << props[i] << std::endl;

    std::cout << "halfedge properties:\n";
    props = halfedge_properties();
    for (unsigned int i=0; i<props.size(); ++i)
        std::cout << "\t" << props[i] << std::endl;

    std::cout << "edge properties:\n";
    props = edge_properties();
    for (unsigned int i=0; i<props.size(); ++i)
        std::cout << "\t" << props[i] << std::endl;

    std::cout << "face properties:\n";
    props = face_properties();
    for (unsigned int i=0; i<props.size(); ++i)
        std::cout << "\t" << props[i] << std::endl;
}

gsMeshTopology::Halfedge
gsMeshTopology::
find_halfedge(Vertex start, Vertex end) const
{
    assert(is_valid(start) && is_valid(end));

    Halfedge h  = halfedge(start);
    const Halfedge hh = h;

    if (h.is_valid())
    {
        do
        {
            if (to_vertex(h) == end)
                return h;
            h = cw_rotated_halfedge(h);
        }
        while (h != hh);
    }

    return Halfedge();
}

gsMeshTopology::Edge
gsMeshTopology::
find_edge(Vertex a, Vertex b) const
{
    Halfedge h = find_halfedge(a,b);
    return h.is_valid() ? edge(h) : Edge();
}

void
gsMeshTopology::
adjust_outgoing_halfedge(Vertex v)
{
    Halfedge h  = halfedge(v);
    const Halfedge hh = h;

    if (h.is_valid() && (h != cw_rotated_halfedge(h)) )
    {
        do
        {
            if (is_boundary(h))
            {
                set_halfedge(v, h);
                return;
            }
            h = cw_rotated_halfedge(h);
        }
        while (h != hh);
    }
}

gsMeshTopology::Face
gsMeshTopology::
add_triangle(Vertex v0, Vertex v1, Vertex v2)
{
    add_face_vertices_.resize(3);
    add_face_vertices_[0] = v0;
    add_face_vertices_[1] = v1;
    add_face_vertices_[2] = v2;
    return add_face(add_face_vertices_);
}

gsMeshTopology::Face
gsMeshTopology::
add_quad(Vertex v0, Vertex v1, Vertex v2, Vertex v3)
{
    add_face_vertices_.resize(4);
    add_face_vertices_[0] = v0;
    add_face_vertices_[1] = v1;
    add_face_vertices_[2] = v2;
    add_face_vertices_[3] = v3;
    return add_face(add_face_vertices_);
}


gsMeshTopology::Face
gsMeshTopology::
add_face(const std::vector<Vertex>& vertices)
{
    const unsigned int n(vertices.size());
    assert (n > 2);

    Vertex        v;
    unsigned int  i, ii, id;
    Halfedge      inner_next, inner_prev, outer_next, outer_prev, boundary_next, boundary_prev, patch_start, patch_end;


    // use global arrays to avoid new/delete of local arrays!!!
    std::vector<Halfedge>&  halfedges    = add_face_halfedges_;
    std::vector<bool>&      is_new       = add_face_is_new_;
    std::vector<bool>&      needs_adjust = add_face_needs_adjust_;
    NextCache&              next_cache   = add_face_next_cache_;
    halfedges.clear();
    halfedges.resize(n);
    is_new.clear();
    is_new.resize(n);
    needs_adjust.clear();
    needs_adjust.resize(n, false);
    next_cache.clear();
    next_cache.reserve(3*n);

    // test for topological errors
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
        if ( !is_boundary(vertices[i]) )
        {
            std::cerr << "gsMeshTopology::add_face: complex vertex v"<<vertices[i].idx()<<".\n";
            return Face();
        }

        halfedges[i] = find_halfedge(vertices[i], vertices[ii]);
        is_new[i]    = !halfedges[i].is_valid();

        if (!is_new[i] && !is_boundary(halfedges[i]))
        {
            std::cerr << "gsMeshTopology::add_face: complex edge at v"<<from_vertex(halfedges[i]).idx()<<" -- v"<<to_vertex(halfedges[i]).idx()<<".\n";
            return Face();
        }
    }

    // re-link patches if necessary
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
        if (!is_new[i] && !is_new[ii])
        {
            inner_prev = halfedges[i];
            inner_next = halfedges[ii];

            if (next_halfedge(inner_prev) != inner_next)
            {
                // here comes the ugly part... we have to relink a whole patch

                // search a free gap
                // free gap will be between boundary_prev and boundary_next
                outer_prev = opposite_halfedge(inner_next);
                outer_next = opposite_halfedge(inner_prev);
                boundary_prev = outer_prev;
                do
                {
                    boundary_prev = opposite_halfedge(next_halfedge(boundary_prev));
                }
                while (!is_boundary(boundary_prev) || boundary_prev==inner_prev);
                boundary_next = next_halfedge(boundary_prev);
                assert(is_boundary(boundary_prev));
                assert(is_boundary(boundary_next));


                // ok ?
                if (boundary_next == inner_next)
                {
                    std::cerr << "gsMeshTopology::add_face: patch re-linking failed\n";
                    return Face();
                }

                // other halfedges' handles
                patch_start = next_halfedge(inner_prev);
                patch_end   = prev_halfedge(inner_next);

                // relink
                next_cache.push_back(NextCacheEntry(boundary_prev, patch_start));
                next_cache.push_back(NextCacheEntry(patch_end, boundary_next));
                next_cache.push_back(NextCacheEntry(inner_prev, inner_next));
            }
        }
    }



    // create missing edges
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
        if (is_new[i])
        {
            halfedges[i] = new_edge(vertices[i], vertices[ii]);
        }
    }



    // create the face
    Face f(new_face());
    set_halfedge(f, halfedges[n-1]);



    // setup halfedges
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
        v          = vertices[ii];
        inner_prev = halfedges[i];
        inner_next = halfedges[ii];

        id = 0;
        if (is_new[i])  id |= 1;
        if (is_new[ii]) id |= 2;

        if (id)
        {
            outer_prev = opposite_halfedge(inner_next);
            outer_next = opposite_halfedge(inner_prev);

            // set outer links
            switch (id)
            {
            case 1: // prev is new, next is old
                boundary_prev = prev_halfedge(inner_next);
                next_cache.push_back(NextCacheEntry(boundary_prev, outer_next));
                set_halfedge(v, outer_next);
                break;

            case 2: // next is new, prev is old
                boundary_next = next_halfedge(inner_prev);
                next_cache.push_back(NextCacheEntry(outer_prev, boundary_next));
                set_halfedge(v, boundary_next);
                break;

            case 3: // both are new
                if (!halfedge(v).is_valid())
                {
                    set_halfedge(v, outer_next);
                    next_cache.push_back(NextCacheEntry(outer_prev, outer_next));
                }
                else
                {
                    boundary_next = halfedge(v);
                    boundary_prev = prev_halfedge(boundary_next);
                    next_cache.push_back(NextCacheEntry(boundary_prev, outer_next));
                    next_cache.push_back(NextCacheEntry(outer_prev, boundary_next));
                }
                break;
            }

            // set inner link
            next_cache.push_back(NextCacheEntry(inner_prev, inner_next));
        }
        else needs_adjust[ii] = (halfedge(v) == inner_next);


        // set face handle
        set_face(halfedges[i], f);
    }



    // process next halfedge cache
    typename NextCache::const_iterator ncIt(next_cache.begin()), ncEnd(next_cache.end());
    for (; ncIt != ncEnd; ++ncIt)
    {
        set_next_halfedge(ncIt->first, ncIt->second);
    }



    // adjust vertices' halfedge handle
    for (i=0; i<n; ++i)
    {
        if (needs_adjust[i])
        {
            adjust_outgoing_halfedge(vertices[i]);
        }
    }


    return f;
}

unsigned int
gsMeshTopology::
valence(Vertex v) const
{
    unsigned int count(0);
    Halfedge h  = halfedge(v);
    if (h.is_valid())
    {
        const Halfedge hh = h;
        do
        {
            ++count;
            h = cw_rotated_halfedge(h);
        }
        while (h != hh);
    }
    return count;
}

unsigned int
gsMeshTopology::
valence(Face f) const
{
    unsigned int count(0);
    Halfedge h  = halfedge(f);
    const Halfedge hh = h;
    do
    {
        ++count;
        h = next_halfedge(h);
    }
    while (h != hh);

    return count;
}

void
gsMeshTopology::
delete_vertex(Vertex v)
{
    if (vdeleted_[v])  return;

    // collect incident faces
    std::vector<Face> incident_faces;
    incident_faces.reserve(6);

    Face_around_vertex_circulator fc, fc_end;
    fc = fc_end = faces(v);

    if (fc)
        do
        {
            incident_faces.push_back(*fc);
        } while (++fc != fc_end);

    // delete incident faces
    typename std::vector<Face>::iterator fit(incident_faces.begin()),
        fend(incident_faces.end());

    for (; fit != fend; ++fit)
        delete_face(*fit);

    // mark v as deleted if not yet done by delete_face()
    if (!vdeleted_[v])
    {
        vdeleted_[v] = true;
        deleted_vertices_++;
        garbage_ = true;
    }
}

void
gsMeshTopology::
delete_edge(Edge e)
{
    if (edeleted_[e])  return;

    Face f0 = face(halfedge(e, 0));
    Face f1 = face(halfedge(e, 1));

    if (f0.is_valid()) delete_face(f0);
    if (f1.is_valid()) delete_face(f1);
}

void
gsMeshTopology::
delete_face(Face f)
{
    if (fdeleted_[f])  return;

    // mark face deleted
    if (!fdeleted_[f])
    {
        fdeleted_[f] = true;
        deleted_faces_++;
    }

    // boundary edges of face f to be deleted
    std::vector<Edge> deleted_edges;
    deleted_edges.reserve(3);


    // vertices of face f for updating their outgoing halfedge
    std::vector<Vertex> vertices;
    vertices.reserve(3);


    // for all halfedges of face f do:
    //   1) invalidate face handle.
    //   2) collect all boundary halfedges, set them deleted
    //   3) store vertex handles
    Halfedge_around_face_circulator hc, hc_end;
    hc = hc_end = halfedges(f);

    do
    {
        set_face(*hc, Face());

        if (is_boundary(opposite_halfedge(*hc)))
            deleted_edges.push_back(edge(*hc));

        vertices.push_back(to_vertex(*hc));

    } while (++hc != hc_end);


    // delete all collected (half)edges
    // delete isolated vertices
    if (!deleted_edges.empty())
    {
        typename std::vector<Edge>::iterator del_it(deleted_edges.begin()),
            del_end(deleted_edges.end());

        Halfedge h0, h1, next0, next1, prev0, prev1;
        Vertex   v0, v1;

        for (; del_it!=del_end; ++del_it)
        {
            h0    = halfedge(*del_it, 0);
            v0    = to_vertex(h0);
            next0 = next_halfedge(h0);
            prev0 = prev_halfedge(h0);

            h1    = halfedge(*del_it, 1);
            v1    = to_vertex(h1);
            next1 = next_halfedge(h1);
            prev1 = prev_halfedge(h1);

            // adjust next and prev handles
            set_next_halfedge(prev0, next1);
            set_next_halfedge(prev1, next0);

            // mark edge deleted
            if (!edeleted_[*del_it])
            {
                edeleted_[*del_it] = true;
                deleted_edges_++;
            }

            // update v0
            if (halfedge(v0) == h1)
            {
                if (next0 == h1)
                {
                    if (!vdeleted_[v0])
                    {
                        vdeleted_[v0] = true;
                        deleted_vertices_++;
                    }
                }
                else set_halfedge(v0, next0);
            }

            // update v1
            if (halfedge(v1) == h0)
            {
                if (next1 == h0)
                {
                    if (!vdeleted_[v1])
                    {
                        vdeleted_[v1] = true;
                        deleted_vertices_++;
                    }
                }
                else  set_halfedge(v1, next1);
            }
        }
    }


    // update outgoing halfedge handles of remaining vertices
    typename std::vector<Vertex>::iterator v_it(vertices.begin()),
        v_end(vertices.end());
    for (; v_it!=v_end; ++v_it)
        adjust_outgoing_halfedge(*v_it);

    garbage_ = true;
}

void
gsMeshTopology::
garbage_collection()
{
    int  i, i0, i1,
        nV(vertices_size()),
        nE(edges_size()),
        nH(halfedges_size()),
        nF(faces_size());

    Vertex    v;
    Halfedge  h;
    Face      f;


    // setup handle mapping
    Vertex_property<Vertex>      vmap = add_vertex_property<Vertex>("v:garbage-collection");
    Halfedge_property<Halfedge>  hmap = add_halfedge_property<Halfedge>("h:garbage-collection");
    Face_property<Face>          fmap = add_face_property<Face>("f:garbage-collection");
    for (i=0; i<nV; ++i)
        vmap[Vertex(i)] = Vertex(i);
    for (i=0; i<nH; ++i)
        hmap[Halfedge(i)] = Halfedge(i);
    for (i=0; i<nF; ++i)
        fmap[Face(i)] = Face(i);



    // remove deleted vertices
    if (nV > 0)
    {
        i0=0;  i1=nV-1;

        while (1)
        {
            // find first deleted and last un-deleted
            while (!vdeleted_[Vertex(i0)] && i0 < i1)  ++i0;
            while ( vdeleted_[Vertex(i1)] && i0 < i1)  --i1;
            if (i0 >= i1) break;

            // swap
            vprops_.swap(i0, i1);
        };

        // remember new size
        nV = vdeleted_[Vertex(i0)] ? i0 : i0+1;
    }


    // remove deleted edges
    if (nE > 0)
    {
        i0=0;  i1=nE-1;

        while (1)
        {
            // find first deleted and last un-deleted
            while (!edeleted_[Edge(i0)] && i0 < i1)  ++i0;
            while ( edeleted_[Edge(i1)] && i0 < i1)  --i1;
            if (i0 >= i1) break;

            // swap
            eprops_.swap(i0, i1);
            hprops_.swap(2*i0,   2*i1);
            hprops_.swap(2*i0+1, 2*i1+1);
        };

        // remember new size
        nE = edeleted_[Edge(i0)] ? i0 : i0+1;
        nH = 2*nE;
    }


    // remove deleted faces
    if (nF > 0)
    {
        i0=0;  i1=nF-1;

        while (1)
        {
            // find 1st deleted and last un-deleted
            while (!fdeleted_[Face(i0)] && i0 < i1)  ++i0;
            while ( fdeleted_[Face(i1)] && i0 < i1)  --i1;
            if (i0 >= i1) break;

            // swap
            fprops_.swap(i0, i1);
        };

        // remember new size
        nF = fdeleted_[Face(i0)] ? i0 : i0+1;
    }


    // update vertex connectivity
    for (i=0; i<nV; ++i)
    {
        v = Vertex(i);
        if (!is_isolated(v))
            set_halfedge(v, hmap[halfedge(v)]);
    }


    // update halfedge connectivity
    for (i=0; i<nH; ++i)
    {
        h = Halfedge(i);
        set_vertex(h, vmap[to_vertex(h)]);
        set_next_halfedge(h, hmap[next_halfedge(h)]);
        if (!is_boundary(h))
            set_face(h, fmap[face(h)]);
    }


    // update handles of faces
    for (i=0; i<nF; ++i)
    {
        f = Face(i);
        set_halfedge(f, hmap[halfedge(f)]);
    }


    // let derived classes remap their own handle-valued properties
    remap_handles(vmap, hmap, fmap);

    // remove handle maps
    remove_vertex_property(vmap);
    remove_halfedge_property(hmap);
    remove_face_property(fmap);


    // finally resize arrays
    vprops_.resize(nV); vprops_.free_memory();
    hprops_.resize(nH); hprops_.free_memory();
    eprops_.resize(nE); eprops_.free_memory();
    fprops_.resize(nF); fprops_.free_memory();

    deleted_vertices_ = deleted_edges_ = deleted_faces_ = 0;
    garbage_ = false;
}

// e(v1,v0): h0(v1->v0) and h1(v0->v1)
// v0 = vertex(e,0) ==   to_vertex(h0)  == from_vertex(h1)
// v1 = vertex(e,1) == from_vertex (h0) ==  to_vertex(h1)
// Insert a new quad by four consecutive border edges
gsMeshTopology::Face gsMeshTopology::add_quad(Edge e0, Edge e1,
                                      Edge e2, Edge e3)
{
    Halfedge h[4];
    for (int i : {0,1})
        for (int j : {0,1})
        {
            if ( vertex(e0,i) == vertex(e1,j) )
            {
                h[0] = halfedge(e0, i);
                h[1] = halfedge(e1,!j);
                break;
            }
        }

    GISMO_ASSERT( is_valid(h[0]) && is_valid(h[1]), "Error in gsSurfMesh::add_quad");
    h[2] = ( to_vertex(h[1]) == vertex(e2,0) ? halfedge(e2,1) : halfedge(e2,0) );
    h[3] = ( to_vertex(h[2]) == vertex(e3,0) ? halfedge(e3,1) : halfedge(e3,0) );
    GISMO_ASSERT( to_vertex(h[1]) == from_vertex(h[2]), "Error in gsSurfMesh::add_quad");
    GISMO_ASSERT( to_vertex(h[2]) == from_vertex(h[3]), "Error in gsSurfMesh::add_quad");
    GISMO_ASSERT( to_vertex(h[3]) == from_vertex(h[0]), "Error in gsSurfMesh::add_quad");

    for(int i = 0; i!=4; ++i)
    {
        if (is_valid(face(h[i])))
        {
            std::cerr << "gsSurfMesh::add_face: discovered a complex edge f"
                      << face(h[i]).idx()<<": "
                      << from_vertex(h[i]).idx()<<"~"<<to_vertex(h[i]).idx()<<".\n";

            // Try to add new twin edges for the extra faces
            auto start = from_vertex(h[i]);
            auto  end  = to_vertex(h[i]);
            for( Halfedge hh : halfedges() )
                if ( from_vertex(hh) == start && to_vertex(hh)== end
                     && hh !=h[i] && !is_valid(face(hh)) )
                {
                    h[i] = hh;
                    return add_quad(edge(h[0]),edge(h[1]),edge(h[2]),edge(h[3]));
                }

            Edge e = add_edge( from_vertex(h[i]), to_vertex(h[i]) );
            h[i] = halfedge(e, 0);
            return add_quad(edge(h[0]),edge(h[1]),edge(h[2]),edge(h[3]));
        }

        if ( !is_boundary( from_vertex(h[i]) ) )
        {
            std::cerr << "gsSurfMesh::add_face: discovered a complex vertex v"
                      <<from_vertex(h[i]).idx()<<".\n";
            return Face();
        }
    }

    return new_quad(h,4);
}

gsMeshTopology::Face gsMeshTopology::new_quad(Halfedge * h, int sz)
{
    GISMO_UNUSED(sz);
    Face f = new_face();
    set_halfedge(f, h[0]);
    for (index_t i = 0; i!=4; ++i)
    {
        set_face(h[i], f);
        set_next_halfedge(h[i], h[(i+1)%4]);
    }

    Halfedge ho, hl, hr; // relink boundary
    for (index_t i = 0; i!=4; ++i)
    {
        ho = opposite_halfedge(h[i]);
        hr = ho;
        do //CW
        {
            if (is_boundary(hr))
            {
                hl = h[i];
                do //CCW
                {
                    hl = prev_halfedge(hl);
                    hl = opposite_halfedge(hl);
                    if (is_boundary(hl))
                    {
                        set_next_halfedge(hr,hl);
                        break;
                    }
                }
                while (hl!=h[i]);                
            }
            hr = next_halfedge(hr);
            hr = opposite_halfedge(hr);
        }
        while ( hr!=ho );
    }

    for (index_t i = 0; i!=4; ++i)
    {
        set_halfedge(from_vertex(h[i]), h[i]);
        adjust_outgoing_halfedge(from_vertex(h[i]));
    }

    return f;
}
  

void gsMeshTopology::add_batch_vertices(size_t nverts)
{
    for(size_t i = 0; i!=nverts; ++i)
        new_vertex();
}

} // namespace gismo
