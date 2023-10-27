/** @file gsSurfMesh.cpp

    @brief Half edge mesh structure

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#include <gsCore/gsTemplateTools.h>

#include <gsMesh2/gsSurfMesh.h>
#include <gsMesh2/IO.h>

#include <gsCore/gsMultiPatch.h>
#include <gsNurbs/gsTensorBSpline.h>

//== NAMESPACE ================================================================


namespace gismo {

//== IMPLEMENTATION ===========================================================


gsSurfMesh::
gsSurfMesh()
{
    // allocate standard properties
    // same list is used in operator=() and assign()
    vconn_    = add_vertex_property<Vertex_connectivity>("v:connectivity");
    hconn_    = add_halfedge_property<Halfedge_connectivity>("h:connectivity");
    fconn_    = add_face_property<Face_connectivity>("f:connectivity");
    vpoint_   = add_vertex_property<Point>("v:point",Point(0,0,0));
    vdeleted_ = add_vertex_property<bool>("v:deleted", false);
    edeleted_ = add_edge_property<bool>("e:deleted", false);
    fdeleted_ = add_face_property<bool>("f:deleted", false);

    mprops_.push_back();

    deleted_vertices_ = deleted_edges_ = deleted_faces_ = 0;
    garbage_ = false;
}


gsSurfMesh::
gsSurfMesh(const gsMatrix<Scalar> & pts)
{
    for( auto & col : pts.colwise() )
        this->add_vertex(col);
}

//-----------------------------------------------------------------------------


gsSurfMesh::
~gsSurfMesh()
{
}


//-----------------------------------------------------------------------------


gsSurfMesh&
gsSurfMesh::
operator=(const gsSurfMesh& rhs)
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
        vpoint_   = get_vertex_property<Point>("v:point");

        // normals might be there, therefore use get_property
        vnormal_  = get_vertex_property<Point>("v:normal");
        fnormal_  = get_face_property<Point>("f:normal");

        // how many elements are deleted?
        deleted_vertices_ = rhs.deleted_vertices_;
        deleted_edges_    = rhs.deleted_edges_;
        deleted_faces_    = rhs.deleted_faces_;
        garbage_          = rhs.garbage_;
    }

    return *this;
}


//-----------------------------------------------------------------------------


gsSurfMesh&
gsSurfMesh::
assign(const gsSurfMesh& rhs)
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
        vpoint_   = add_vertex_property<Point>("v:point",Point(0,0,0));
        vdeleted_ = add_vertex_property<bool>("v:deleted", false);
        edeleted_ = add_edge_property<bool>("e:deleted", false);
        fdeleted_ = add_face_property<bool>("f:deleted", false);

        // normals might be there, therefore use get_property
        vnormal_  = get_vertex_property<Point>("v:normal");
        fnormal_  = get_face_property<Point>("f:normal");

        // copy properties from other mesh
        vconn_.array()     = rhs.vconn_.array();
        hconn_.array()     = rhs.hconn_.array();
        fconn_.array()     = rhs.fconn_.array();
        vpoint_.array()    = rhs.vpoint_.array();
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


//-----------------------------------------------------------------------------


bool
gsSurfMesh::
read(const std::string& filename)
{
    return read_mesh(*this, filename);
}


//-----------------------------------------------------------------------------


bool
gsSurfMesh::
write(const std::string& filename) const
{
    return write_mesh(*this, filename);
}


//-----------------------------------------------------------------------------


void
gsSurfMesh::
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


//-----------------------------------------------------------------------------


void
gsSurfMesh::
free_memory()
{
    vprops_.free_memory();
    hprops_.free_memory();
    eprops_.free_memory();
    fprops_.free_memory();
    mprops_.free_memory();
}


//-----------------------------------------------------------------------------


void
gsSurfMesh::
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


//-----------------------------------------------------------------------------


void
gsSurfMesh::
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


//-----------------------------------------------------------------------------


gsSurfMesh::Vertex
gsSurfMesh::
add_vertex(const Point& p)
{
    Vertex v = new_vertex();
    vpoint_[v] = p;
    return v;
}


//-----------------------------------------------------------------------------


gsSurfMesh::Halfedge
gsSurfMesh::
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


//-----------------------------------------------------------------------------


gsSurfMesh::Edge
gsSurfMesh::
find_edge(Vertex a, Vertex b) const
{
    Halfedge h = find_halfedge(a,b);
    return h.is_valid() ? edge(h) : Edge();
}


//-----------------------------------------------------------------------------


void
gsSurfMesh::
adjust_outgoing_halfedge(Vertex v)
{
    Halfedge h  = halfedge(v);
    const Halfedge hh = h;

    if (h.is_valid())
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


//-----------------------------------------------------------------------------


gsSurfMesh::Face
gsSurfMesh::
add_triangle(Vertex v0, Vertex v1, Vertex v2)
{
    add_face_vertices_.resize(3);
    add_face_vertices_[0] = v0;
    add_face_vertices_[1] = v1;
    add_face_vertices_[2] = v2;
    return add_face(add_face_vertices_);
}


//-----------------------------------------------------------------------------


gsSurfMesh::Face
gsSurfMesh::
add_quad(Vertex v0, Vertex v1, Vertex v2, Vertex v3)
{
    add_face_vertices_.resize(4);
    add_face_vertices_[0] = v0;
    add_face_vertices_[1] = v1;
    add_face_vertices_[2] = v2;
    add_face_vertices_[3] = v3;
    return add_face(add_face_vertices_);
}


//-----------------------------------------------------------------------------


gsSurfMesh::Face
gsSurfMesh::
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
            std::cerr << "gsSurfMesh::add_face: complex vertex\n";
            return Face();
        }

        halfedges[i] = find_halfedge(vertices[i], vertices[ii]);
        is_new[i]    = !halfedges[i].is_valid();

        if (!is_new[i] && !is_boundary(halfedges[i]))
        {
            std::cerr << "gsSurfMesh::add_face: complex edge\n";
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
                    std::cerr << "gsSurfMeshT::add_face: patch re-linking failed\n";
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
    NextCache::const_iterator ncIt(next_cache.begin()), ncEnd(next_cache.end());
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


//-----------------------------------------------------------------------------


unsigned int
gsSurfMesh::
valence(Vertex v) const
{
    unsigned int count(0);

    Vertex_around_vertex_circulator vvit = vertices(v);
    Vertex_around_vertex_circulator vvend = vvit;
    if (vvit) do
              {
                  ++count;
              } while (++vvit != vvend);

    return count;
}


//-----------------------------------------------------------------------------


unsigned int
gsSurfMesh::
valence(Face f) const
{
    unsigned int count(0);

    Vertex_around_face_circulator fvit = vertices(f);
    Vertex_around_face_circulator fvend = fvit;
    do {
        ++count;
    } while (++fvit != fvend);

    return count;
}


unsigned int
gsSurfMesh::
face_valence_sum() const
{
    unsigned int count = 0;
#   pragma omp parallel for reduction(+:count)
    for (auto fit = faces_begin(); fit < faces_end(); ++fit)
        count += valence(*fit);
    return count;
}

//-----------------------------------------------------------------------------


bool
gsSurfMesh::
is_triangle_mesh() const
{
    Face_iterator fit=faces_begin(), fend=faces_end();
    for (; fit!=fend; ++fit)
        if (valence(*fit) != 3)
            return false;

    return true;
}


//-----------------------------------------------------------------------------


bool
gsSurfMesh::
is_quad_mesh() const
{
    Face_iterator fit=faces_begin(), fend=faces_end();
    for (; fit!=fend; ++fit)
        if (valence(*fit) != 4)
            return false;

    return true;
}


//-----------------------------------------------------------------------------


void
gsSurfMesh::
triangulate()
{
    /* The iterators will stay valid, even though new faces are added,
       because they are now implemented index-based instead of
       pointer-based.
    */
    Face_iterator fit=faces_begin(), fend=faces_end();
    for (; fit!=fend; ++fit)
        triangulate(*fit);
}


//-----------------------------------------------------------------------------


void
gsSurfMesh::
triangulate(Face f)
{
    /*
      Split an arbitrary face into triangles by connecting
      each vertex of fh after its second to vh.

      - fh will remain valid (it will become one of the
      triangles)
      - the halfedge handles of the new triangles will
      point to the old halfedges
    */

    Halfedge base_h  = halfedge(f);
    Vertex   start_v = from_vertex(base_h);
    Halfedge next_h  = next_halfedge(base_h);

    while (to_vertex(next_halfedge(next_h)) != start_v)
    {
        Halfedge next_next_h(next_halfedge(next_h));

        Face new_f = new_face();
        set_halfedge(new_f, base_h);

        Halfedge new_h = new_edge(to_vertex(next_h), start_v);

        set_next_halfedge(base_h, next_h);
        set_next_halfedge(next_h, new_h);
        set_next_halfedge(new_h,  base_h);

        set_face(base_h, new_f);
        set_face(next_h, new_f);
        set_face(new_h,  new_f);

        base_h = opposite_halfedge(new_h);
        next_h = next_next_h;
    }
    set_halfedge(f, base_h);  //the last face takes the handle _fh

    set_next_halfedge(base_h, next_h);
    set_next_halfedge(next_halfedge(next_h), base_h);

    set_face(base_h, f);
}


//-----------------------------------------------------------------------------


void
gsSurfMesh::
update_face_normals()
{
    if (!fnormal_)
        fnormal_ = face_property<Point>("f:normal",Point(0,0,0));

    Face_iterator fit, fend=faces_end();

    for (fit=faces_begin(); fit!=fend; ++fit)
        fnormal_[*fit] = compute_face_normal(*fit);
}


//-----------------------------------------------------------------------------


Normal
gsSurfMesh::
compute_face_normal(Face f) const
{
    Halfedge h = halfedge(f);
    Halfedge hend = h;

    Point p0 = vpoint_[to_vertex(h)];
    h = next_halfedge(h);
    Point p1 = vpoint_[to_vertex(h)];
    h = next_halfedge(h);
    Point p2 = vpoint_[to_vertex(h)];

    if (next_halfedge(h) == hend) // face is a triangle
    {
        p2-=p1; p0-=p1;
        return p2.cross(p1).normalized();
    }

    else // face is a general polygon
    {
        Normal n(0,0,0);

        hend = h;
        do
        {
            n += (p2-p1).cross(p0-p1);
            h  = next_halfedge(h);
            p0 = p1;
            p1 = p2;
            p2 = vpoint_[to_vertex(h)];
        }
        while (h != hend);

        return n.normalized();
    }
}


//-----------------------------------------------------------------------------


void
gsSurfMesh::
update_vertex_normals()
{
    if (!vnormal_)
        vnormal_ = vertex_property<Point>("v:normal",Point(0,0,0));

    Vertex_iterator vit, vend=vertices_end();

    for (vit=vertices_begin(); vit!=vend; ++vit)
        vnormal_[*vit] = compute_vertex_normal(*vit);
}


//-----------------------------------------------------------------------------


Normal
gsSurfMesh::
compute_vertex_normal(Vertex v) const
{
    Point     nn(0,0,0);
    Halfedge  h = halfedge(v);

    if (h.is_valid())
    {
        const Halfedge hend = h;
        const Point p0 = vpoint_[v];

        Point   n, p1, p2;
        Scalar  cosine, angle, denom;

        do
        {
            if (!is_boundary(h))
            {
                p1 = vpoint_[to_vertex(h)];
                p1 -= p0;

                p2 = vpoint_[from_vertex(prev_halfedge(h))];
                p2 -= p0;

                // check whether we can robustly compute angle
                denom = sqrt(p1.squaredNorm()*p2.squaredNorm());
                if (denom > std::numeric_limits<Scalar>::min())
                {
                    cosine = p1.dot(p2) / denom;
                    if      (cosine < -1.0) cosine = -1.0;
                    else if (cosine >  1.0) cosine =  1.0;
                    angle = acos(cosine);

                    n   = p1.cross(p2);

                    // check whether normal is != 0
                    denom = n.norm();
                    if (denom > std::numeric_limits<Scalar>::min())
                    {
                        n  *= angle/denom;
                        nn += n;
                    }
                }
            }

            h  = cw_rotated_halfedge(h);
        }
        while (h != hend);

        nn.normalize();
    }

    return nn;
}


//-----------------------------------------------------------------------------


gsSurfMesh::Scalar
gsSurfMesh::
edge_length(Edge e) const
{
    return (vpoint_[vertex(e,0)] - vpoint_[vertex(e,1)]).norm();
}


//-----------------------------------------------------------------------------


void
gsSurfMesh::
split(Face f, Vertex v)
{
    /*
      Split an arbitrary face into triangles by connecting each vertex of fh to vh.
      - fh will remain valid (it will become one of the triangles)
      - the halfedge handles of the new triangles will point to the old halfeges
    */

    Halfedge hend = halfedge(f);
    Halfedge h    = next_halfedge(hend);

    Halfedge hold = new_edge(to_vertex(hend), v);

    set_next_halfedge(hend, hold);
    set_face(hold, f);

    hold = opposite_halfedge(hold);

    while (h != hend)
    {
        Halfedge hnext = next_halfedge(h);

        Face fnew = new_face();
        set_halfedge(fnew, h);

        Halfedge hnew = new_edge(to_vertex(h), v);

        set_next_halfedge(hnew, hold);
        set_next_halfedge(hold, h);
        set_next_halfedge(h,    hnew);

        set_face(hnew, fnew);
        set_face(hold, fnew);
        set_face(h,    fnew);

        hold = opposite_halfedge(hnew);

        h = hnext;
    }

    set_next_halfedge(hold, hend);
    set_next_halfedge(next_halfedge(hend), hold);

    set_face(hold, f);

    set_halfedge(v, hold);
}


void
gsSurfMesh::
quad_split(Face f, Vertex v, Halfedge s)
{
    /*
      Split an arbitrary face into quads by connecting each vertex of fh to vh.
      - fh will remain valid (it will become one of the quads)
    */

    //assert: number of vertices is even (4,6,..)
    //assert: vertex s is on face
    //assert: vertex is isolated

    //std::cout<< "s: "<< from_vertex(s) <<"->"<<to_vertex(s) <<std::endl;
    set_halfedge(f,s);
    Halfedge hnext = next_halfedge(next_halfedge(s));
    Halfedge h  = s;
    //std::cout<< "h: "<< from_vertex(h) <<"->"<<to_vertex(h) <<std::endl;
        
    //first step
    Halfedge e0 = new_edge(v,from_vertex(h));
    set_halfedge(v, e0);
    set_face(e0, f);
    set_next_halfedge(e0, h);//sets next(e0) and also prev(h)

    // ---------------
    Halfedge e = new_edge(from_vertex(hnext),v);
    set_face(e, f);

    set_next_halfedge(prev_halfedge(hnext),e);//sets next(h) and also prev(last)
    set_next_halfedge(e,e0);//sets next(h) and also prev(last)

    //std::cout <<"face: ";for (auto fv : vertices(f)) std::cout <<" "<< fv.idx(); std::cout <<"\n";

    //std::cout<< "h: "<< from_vertex(h) <<"->"<<to_vertex(h) <<std::endl;
    e0 = opposite_halfedge(e);//v->hnext
    h = hnext;
    hnext  = next_halfedge(next_halfedge(hnext));

    while (h != s) // face containing h2
    {
        //std::cout<< "e0: "<< from_vertex(e0) <<"->"<<to_vertex(e0) <<std::endl;
        //std::cout<< "h: "<< from_vertex(h) <<"->"<<to_vertex(h) <<std::endl;
        f = new_face();
        e = ( hnext!=s ? new_edge(from_vertex(hnext),v) : 
              opposite_halfedge(halfedge(v)) );
        //std::cout<< "e: "<< from_vertex(e) <<"->"<<to_vertex(e) <<std::endl;
        set_halfedge(f, e0);
        set_face(e0, f); // 1

        set_next_halfedge(e0,h);
        set_face(h, f); // 2

        set_face(next_halfedge(h), f); // 3
        set_next_halfedge(next_halfedge(h),e);

        set_face(e, f); //4
        set_next_halfedge(e,e0);

        //std::cout <<"face: ";for (auto fv : vertices(f)) std::cout <<" "<< fv.idx(); std::cout <<"\n";

        e0 = opposite_halfedge(e);
        h = hnext;
        hnext  = next_halfedge(next_halfedge(hnext));
    }
}

void gsSurfMesh::quad_split()
{
    gsSurfMesh::Vertex v;
    gsSurfMesh::Halfedge he;

    // reserve vertices, edges, faces
    reserve(n_vertices() + n_edges() + n_faces(),
            2 * n_edges(), 4 * n_faces());

    auto points = get_vertex_property<Point>("v:point");

    index_t env = n_vertices(); // edge vertices start here

    // loop over all edges, add edge points
    Point tmp;
    for (auto eit : edges())
    {
        he = halfedge(eit, 0);
        tmp = (points[from_vertex(he)] + points[to_vertex(he)]) / 2;
        v = add_vertex(tmp);
        insert_vertex(he, v);
    }

    index_t fnv = n_vertices(); // face vertices start here

    // loop over all faces, add face points
    for (auto fit : faces())
    {
        auto fv = vertices(fit);
        tmp.setZero();
        for (auto vc = fv.begin(); vc != fv.end(); ++vc, ++vc)
            tmp += points[*vc];
        tmp /= 4;
        add_vertex(tmp);  // vertex gets shifted face id
    }

    int i = 0;
    for (auto fit : faces())
    {
        v = gsSurfMesh::Vertex(fnv + (i++));//face vertex id ?
        //Start from an original vertex
        auto fv = vertices(fit).begin();
        if ((*fv).idx() >= env) ++fv; //todo: add -> operator
        //assert ( (*fv).idx() < nv )
        quad_split(fit, v, fv.he());
    }

}

//-----------------------------------------------------------------------------


gsSurfMesh::Halfedge
gsSurfMesh::
split(Edge e, Vertex v)
{
    Halfedge h0 = halfedge(e, 0);
    Halfedge o0 = halfedge(e, 1);

    Vertex   v2 = to_vertex(o0);

    Halfedge e1 = new_edge(v, v2);
    Halfedge t1 = opposite_halfedge(e1);

    Face     f0 = face(h0);
    Face     f3 = face(o0);

    set_halfedge(v, h0);
    set_vertex(o0, v);

    if (!is_boundary(h0))
    {
        Halfedge h1 = next_halfedge(h0);
        Halfedge h2 = next_halfedge(h1);

        Vertex   v1 = to_vertex(h1);

        Halfedge e0 = new_edge(v, v1);
        Halfedge t0 = opposite_halfedge(e0);

        Face f1 = new_face();
        set_halfedge(f0, h0);
        set_halfedge(f1, h2);

        set_face(h1, f0);
        set_face(t0, f0);
        set_face(h0, f0);

        set_face(h2, f1);
        set_face(t1, f1);
        set_face(e0, f1);

        set_next_halfedge(h0, h1);
        set_next_halfedge(h1, t0);
        set_next_halfedge(t0, h0);

        set_next_halfedge(e0, h2);
        set_next_halfedge(h2, t1);
        set_next_halfedge(t1, e0);
    }
    else
    {
        set_next_halfedge(prev_halfedge(h0), t1);
        set_next_halfedge(t1, h0);
        // halfedge handle of _vh already is h0
    }


    if (!is_boundary(o0))
    {
        Halfedge o1 = next_halfedge(o0);
        Halfedge o2 = next_halfedge(o1);

        Vertex v3 = to_vertex(o1);

        Halfedge e2 = new_edge(v, v3);
        Halfedge t2 = opposite_halfedge(e2);

        Face f2 = new_face();
        set_halfedge(f2, o1);
        set_halfedge(f3, o0);

        set_face(o1, f2);
        set_face(t2, f2);
        set_face(e1, f2);

        set_face(o2, f3);
        set_face(o0, f3);
        set_face(e2, f3);

        set_next_halfedge(e1, o1);
        set_next_halfedge(o1, t2);
        set_next_halfedge(t2, e1);

        set_next_halfedge(o0, e2);
        set_next_halfedge(e2, o2);
        set_next_halfedge(o2, o0);
    }
    else
    {
        set_next_halfedge(e1, next_halfedge(o0));
        set_next_halfedge(o0, e1);
        set_halfedge(v, e1);
    }

    if (halfedge(v2) == h0)
        set_halfedge(v2, t1);

    return t1;
}


//-----------------------------------------------------------------------------


gsSurfMesh::Halfedge
gsSurfMesh::
insert_vertex(Halfedge h0, Vertex v)
{
    // before:
    //
    // v0      h0       v2
    //  o--------------->o
    //   <---------------
    //         o0
    //
    // after:
    //
    // v0  h0   v   h1   v2
    //  o------>o------->o
    //   <------ <-------
    //     o0       o1

    Halfedge h2 = next_halfedge(h0);
    Halfedge o0 = opposite_halfedge(h0);
    Halfedge o2 = prev_halfedge(o0);
    Vertex   v2 = to_vertex(h0);
    Face     fh = face(h0);
    Face     fo = face(o0);

    Halfedge h1 = new_edge(v, v2);
    Halfedge o1 = opposite_halfedge(h1);

    // adjust halfedge connectivity
    set_next_halfedge(h1, h2);
    set_next_halfedge(h0, h1);
    set_vertex(h0, v);
    set_vertex(h1, v2);
    set_face(h1, fh);

    set_next_halfedge(o1, o0);
    set_next_halfedge(o2, o1);
    set_vertex(o1, v);
    set_face(o1, fo);

    // adjust vertex connectivity
    set_halfedge(v2, o1);
    adjust_outgoing_halfedge(v2);
    set_halfedge(v, h1);
    adjust_outgoing_halfedge(v);

    // adjust face connectivity
    if (fh.is_valid()) set_halfedge(fh, h0);
    if (fo.is_valid()) set_halfedge(fo, o1);

    return o1;
}


//-----------------------------------------------------------------------------


gsSurfMesh::Halfedge
gsSurfMesh::
insert_edge(Halfedge h0, Halfedge h1)
{
    assert(face(h0) == face(h1));
    assert(face(h0).is_valid());

    Vertex   v0 = to_vertex(h0);
    Vertex   v1 = to_vertex(h1);

    Halfedge h2 = next_halfedge(h0);
    Halfedge h3 = next_halfedge(h1);

    Halfedge h4 = new_edge(v0, v1);
    Halfedge h5 = opposite_halfedge(h4);

    Face     f0 = face(h0);
    Face     f1 = new_face();

    set_halfedge(f0, h0);
    set_halfedge(f1, h1);

    set_next_halfedge(h0, h4);
    set_next_halfedge(h4, h3);
    set_face(h4, f0);

    set_next_halfedge(h1, h5);
    set_next_halfedge(h5, h2);
    Halfedge h = h2;
    do
    {
        set_face(h, f1);
        h = next_halfedge(h);
    }
    while (h != h2);

    return h4;
}


//-----------------------------------------------------------------------------


bool
gsSurfMesh::
is_flip_ok(Edge e) const
{
    // boundary edges cannot be flipped
    if (is_boundary(e)) return false;

    // check if the flipped edge is already present in the mesh

    Halfedge h0 = halfedge(e, 0);
    Halfedge h1 = halfedge(e, 1);

    Vertex v0 = to_vertex(next_halfedge(h0));
    Vertex v1 = to_vertex(next_halfedge(h1));

    if (v0 == v1)   // this is generally a bad sign !!!
        return false;

    if (find_halfedge(v0, v1).is_valid())
        return false;

    return true;
}


//-----------------------------------------------------------------------------


void
gsSurfMesh::
flip(Edge e)
{
    // CAUTION : Flipping a halfedge may result in
    // a non-manifold mesh, hence check for yourself
    // whether this operation is allowed or not!

    //let's make it sure it is actually checked
    assert(is_flip_ok(e));

    Halfedge a0 = halfedge(e, 0);
    Halfedge b0 = halfedge(e, 1);

    Halfedge a1 = next_halfedge(a0);
    Halfedge a2 = next_halfedge(a1);

    Halfedge b1 = next_halfedge(b0);
    Halfedge b2 = next_halfedge(b1);

    Vertex   va0 = to_vertex(a0);
    Vertex   va1 = to_vertex(a1);

    Vertex   vb0 = to_vertex(b0);
    Vertex   vb1 = to_vertex(b1);

    Face     fa  = face(a0);
    Face     fb  = face(b0);

    set_vertex(a0, va1);
    set_vertex(b0, vb1);

    set_next_halfedge(a0, a2);
    set_next_halfedge(a2, b1);
    set_next_halfedge(b1, a0);

    set_next_halfedge(b0, b2);
    set_next_halfedge(b2, a1);
    set_next_halfedge(a1, b0);

    set_face(a1, fb);
    set_face(b1, fa);

    set_halfedge(fa, a0);
    set_halfedge(fb, b0);

    if (halfedge(va0) == b0)
        set_halfedge(va0, a1);
    if (halfedge(vb0) == a0)
        set_halfedge(vb0, b1);
}


//-----------------------------------------------------------------------------


bool
gsSurfMesh::
is_collapse_ok(Halfedge v0v1)
{
    Halfedge  v1v0(opposite_halfedge(v0v1));
    Vertex    v0(to_vertex(v1v0));
    Vertex    v1(to_vertex(v0v1));
    Vertex    vv, vl, vr;
    Halfedge  h1, h2;


    // the edges v1-vl and vl-v0 must not be both boundary edges
    if (!is_boundary(v0v1))
    {
        vl = to_vertex(next_halfedge(v0v1));
        h1 = next_halfedge(v0v1);
        h2 = next_halfedge(h1);
        if (is_boundary(opposite_halfedge(h1)) && is_boundary(opposite_halfedge(h2)))
            return false;
    }


    // the edges v0-vr and vr-v1 must not be both boundary edges
    if (!is_boundary(v1v0))
    {
        vr = to_vertex(next_halfedge(v1v0));
        h1 = next_halfedge(v1v0);
        h2 = next_halfedge(h1);
        if (is_boundary(opposite_halfedge(h1)) && is_boundary(opposite_halfedge(h2)))
            return false;
    }


    // if vl and vr are equal or both invalid -> fail
    if (vl == vr) return false;


    // edge between two boundary vertices should be a boundary edge
    if ( is_boundary(v0) && is_boundary(v1) &&
         !is_boundary(v0v1) && !is_boundary(v1v0))
        return false;


    // test intersection of the one-rings of v0 and v1
    Vertex_around_vertex_circulator vv_it, vv_end;
    vv_it = vv_end = vertices(v0);
    do
    {
        vv = *vv_it;
        if (vv != v1 && vv != vl && vv != vr)
            if (find_halfedge(vv, v1).is_valid())
                return false;
    }
    while (++vv_it != vv_end);


    // passed all tests
    return true;
}


//-----------------------------------------------------------------------------


void
gsSurfMesh::
collapse(Halfedge h)
{
    Halfedge h0 = h;
    Halfedge h1 = prev_halfedge(h0);
    Halfedge o0 = opposite_halfedge(h0);
    Halfedge o1 = next_halfedge(o0);

    // remove edge
    remove_edge(h0);

    // remove loops
    if (next_halfedge(next_halfedge(h1)) == h1)
        remove_loop(h1);
    if (next_halfedge(next_halfedge(o1)) == o1)
        remove_loop(o1);
}


//-----------------------------------------------------------------------------


void
gsSurfMesh::
remove_edge(Halfedge h)
{
    Halfedge  hn = next_halfedge(h);
    Halfedge  hp = prev_halfedge(h);

    Halfedge  o  = opposite_halfedge(h);
    Halfedge  on = next_halfedge(o);
    Halfedge  op = prev_halfedge(o);

    Face      fh = face(h);
    Face      fo = face(o);

    Vertex    vh = to_vertex(h);
    Vertex    vo = to_vertex(o);



    // halfedge -> vertex
    Halfedge_around_vertex_circulator vh_it, vh_end;
    vh_it = vh_end = halfedges(vo);
    do
    {
        set_vertex(opposite_halfedge(*vh_it), vh);
    }
    while (++vh_it != vh_end);


    // halfedge -> halfedge
    set_next_halfedge(hp, hn);
    set_next_halfedge(op, on);


    // face -> halfedge
    if (fh.is_valid())  set_halfedge(fh, hn);
    if (fo.is_valid())  set_halfedge(fo, on);


    // vertex -> halfedge
    if (halfedge(vh) == o)  set_halfedge(vh, hn);
    adjust_outgoing_halfedge(vh);
    set_halfedge(vo, Halfedge());


    // delete stuff
    if (!vdeleted_) vdeleted_ = vertex_property<bool>("v:deleted", false);
    if (!edeleted_) edeleted_ = edge_property<bool>("e:deleted", false);
    vdeleted_[vo]      = true; ++deleted_vertices_;
    edeleted_[edge(h)] = true; ++deleted_edges_;
    garbage_ = true;
}


//-----------------------------------------------------------------------------


void
gsSurfMesh::
remove_loop(Halfedge h)
{
    Halfedge  h0 = h;
    Halfedge  h1 = next_halfedge(h0);

    Halfedge  o0 = opposite_halfedge(h0);
    Halfedge  o1 = opposite_halfedge(h1);

    Vertex    v0 = to_vertex(h0);
    Vertex    v1 = to_vertex(h1);

    Face      fh = face(h0);
    Face      fo = face(o0);



    // is it a loop ?
    assert ((next_halfedge(h1) == h0) && (h1 != o0));


    // halfedge -> halfedge
    set_next_halfedge(h1, next_halfedge(o0));
    set_next_halfedge(prev_halfedge(o0), h1);


    // halfedge -> face
    set_face(h1, fo);


    // vertex -> halfedge
    set_halfedge(v0, h1);  adjust_outgoing_halfedge(v0);
    set_halfedge(v1, o1);  adjust_outgoing_halfedge(v1);


    // face -> halfedge
    if (fo.is_valid() && halfedge(fo) == o0)
        set_halfedge(fo, h1);


    // delete stuff
    if (!edeleted_) edeleted_ = edge_property<bool>("e:deleted", false);
    if (!fdeleted_) fdeleted_ = face_property<bool>("f:deleted", false);
    if (fh.is_valid()) { fdeleted_[fh] = true; ++deleted_faces_; }
    edeleted_[edge(h0)] = true; ++deleted_edges_;
    garbage_ = true;
}


//-----------------------------------------------------------------------------


void
gsSurfMesh::
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
    std::vector<Face>::iterator fit(incident_faces.begin()),
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


//-----------------------------------------------------------------------------


void
gsSurfMesh::
delete_edge(Edge e)
{
    if (edeleted_[e])  return;

    Face f0 = face(halfedge(e, 0));
    Face f1 = face(halfedge(e, 1));

    if (f0.is_valid()) delete_face(f0);
    if (f1.is_valid()) delete_face(f1);
}


//-----------------------------------------------------------------------------

void
gsSurfMesh::
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
        std::vector<Edge>::iterator del_it(deleted_edges.begin()),
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
    std::vector<Vertex>::iterator v_it(vertices.begin()),
        v_end(vertices.end());
    for (; v_it!=v_end; ++v_it)
        adjust_outgoing_halfedge(*v_it);

    garbage_ = true;
}


//-----------------------------------------------------------------------------


void
gsSurfMesh::
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


void gsSurfMesh::cc_subdivide()
{
    gsSurfMesh::Vertex v;
    gsSurfMesh::Halfedge he;

    // reserve vertices, edges, faces
    reserve( n_vertices()+n_edges()+n_faces(),
             2*n_edges(), 4*n_faces() );

    auto points = get_vertex_property<Point>("v:point");

    index_t env = n_vertices(); // edge vertices start here

    // loop over all edges, add edge points
    Point tmp;
    for (auto eit : edges())
    {
        he  = halfedge(eit,0);
        tmp = (points[from_vertex(he)]+points[to_vertex(he)]) / 2;
        v   = add_vertex(tmp);
        insert_vertex(he,v);
    }

    index_t fnv = n_vertices(); // face vertices start here

    // loop over all faces, add face points
    for (auto fit : faces())
    {
        auto fv = vertices(fit);
        tmp.setZero();
        for (auto vc = fv.begin(); vc!=fv.end(); ++vc, ++vc)
            tmp += points[*vc];
        tmp /= 4;
        add_vertex(tmp);  // vertex gets shifted face id
    }

    int i = 0;
    for (auto fit : faces())
    {
        v = gsSurfMesh::Vertex(fnv+(i++));//face vertex id ?
        //Start from an original vertex
        auto fv = vertices(fit).begin();
        if ( (*fv).idx() >= env) ++fv; //todo: add -> operator
        //assert ( (*fv).idx() < nv )
        quad_split(fit,v,fv.he());
    }

#   pragma omp parallel for default(none) shared(points,env,fnv) private(v,tmp)
    for (i = env; i<fnv;++i)
    {
        v = gsSurfMesh::Vertex(i); //edge points
        if (is_boundary(v))
        {
            //gsWarn<< "Boundary vertex "<< v.idx() <<"\n";
            continue;
        }
        auto vit = vertices(v);
        auto vcp = vit;
        tmp.setZero();
        if (vit) do
                 {
                     tmp += points[*vit];
                 } while (++vit != vcp);
        tmp /= 4 ; // =valence(v);
        points[v] = tmp;
    }

#   pragma omp parallel for default(none) shared(env,points) private(v)
    for (i = 0; i<env;++i)
    {
        v = gsSurfMesh::Vertex(i); // original vertices
        auto n = valence(v);
        auto & pt = points[v];//original vertex positions are computed using new edge/face points only
        if (is_boundary(v))
        {
            //gsWarn<< "Boundary vertex "<< v.idx() <<"\n";
            if (2>n)
            {
                auto vv = vertices(v);
                pt *= 2;
                pt += points[*vv]; // first boundary neighbor
                pt += points[*(--vv.end())]; // last boundary neighbor
                pt /= 4;               
            }
            continue;
        }

        auto vit = halfedges(v);
        auto vcp = vit;
        //formula: pt = ( (n*(n-3))*points[v] + 4*E - F ) / (n*n);
        pt *= n*(n-3);
        if (vit)
            do
            { //pt += 4*E-F
                pt += 4*points[ to_vertex(*vit) ]
                    - points[ to_vertex(next_halfedge(*vit)) ];
            } while (++vit != vcp);
        pt /= n*n;
    }
}

gsSurfMesh::Vertex_property<Point>
gsSurfMesh::cc_limit_points(std::string label)
{
    auto points = get_vertex_property<Point>("v:point");
    auto limits = add_vertex_property<Point>(
        (label == "v:point" ? "v:limit_points_2022" : label),Point(0,0,0));
    real_t n;
#   pragma omp parallel for default(none) shared(std::cout,points,limits) private(n)
    for (auto vit = vertices_begin(); vit < vertices_end(); ++vit)
    {
        n = valence(*vit);
        if (is_boundary(*vit))
        {
            gsWarn<< "Boundary vertex is ignored.\n";
            
            if (2>n)
            {
                
            }
            continue;
        }

        auto & pt = limits[*vit];
        pt = n*n*points[*vit];
        for ( auto he : halfedges(*vit) )
        {
            if (is_boundary(he))
            {
                gsWarn<< "Boundary halfedge is ignored.\n";                
            }

            pt += 4 * points[ to_vertex(he) ] +
                points[ to_vertex(next_halfedge(he)) ];
        }
        pt /= (n*(n+5));
    }

    if (label == "v:point") //vertices are replaced by their limit positions
    {
        rename_vertex_property(points,"v:point_original");
        rename_vertex_property(limits,"v:point");
    }
    return limits;
}


gsSurfMesh::Vertex_property<Point>
gsSurfMesh::cc_limit_normals(std::string label, bool normalize)
{
    auto points = get_vertex_property<Point>("v:point");
    //todo: check if label exists
    auto limits = add_vertex_property<Point>(label,Point(0,0,0));
    Point t1, t2;
    real_t c1, c2, cc1, cc2;
    index_t i;
    gsSurfMesh::Halfedge h2;
#   pragma omp parallel for default(none) shared(limits,points,normalize) private(h2,t1,t2,c1,c2,cc1,cc2,i)
    for (auto vit = vertices_begin(); vit < vertices_end(); ++vit)
    {
        const real_t n = valence(*vit);
        const real_t cospin = math::cos(EIGEN_PI/n);
        cc2 = 1 / ( n * math::sqrt(4+cospin*cospin) );
        cc1 = 1/n + cospin*cc2;
        t1.setZero();
        t2.setZero();
        i = 0;
        for ( auto he : halfedges(*vit) )
        {
            h2 = ccw_rotated_halfedge(he);
            c1  = math::cos( 2*i   *EIGEN_PI/n)*cc1;
            c2  = math::cos((2*i+1)*EIGEN_PI/n)*cc2;
            t1 += c1 * points[ to_vertex(he ) ]
                + c2 * points[ to_vertex(next_halfedge(he)) ];
            t2 += c1 * points[ to_vertex(h2 ) ]
                + c2 * points[ to_vertex(next_halfedge(h2)) ];
            ++i;
        }
        if (normalize)
            limits[*vit] = t1.cross(t2).normalized();
        else
            limits[*vit] = t1.cross(t2);
    }
    return limits;
}

gsSurfMesh::Vertex_property<Point>
gsSurfMesh::cc_limit_tangent_vec(std::string label, bool normalize)
{
    gsSurfMesh::Vertex v;
    gsSurfMesh::Halfedge h2;

    auto points = get_vertex_property<Point>("v:point");
    //todo: check if label exists
    auto limits = add_vertex_property<Point>(label,Point(0,0,0));
    Point t1, t2;
    real_t c1, c2, cc1, cc2;
    index_t i;
#   pragma omp parallel for default(none) shared(limits,points,normalize) private(v,h2,t1,t2,c1,c2,cc1,cc2,i)
    for (auto vit = vertices_begin(); vit < vertices_end(); ++vit)
    {
        const real_t n = valence(*vit);
        const real_t cospin = math::cos(EIGEN_PI/n);
        cc2 = 1 / ( n * math::sqrt(4+cospin*cospin) );
        cc1 = 1/n + cospin*cc2;
        t1.setZero();
        t2.setZero();
        i = 0;
        for ( auto he : halfedges(*vit) )
        {
            h2 = ccw_rotated_halfedge(he);
            c1  = math::cos( 2*i   *EIGEN_PI/n)*cc1;
            c2  = math::cos((2*i+1)*EIGEN_PI/n)*cc2;
            t1 += c1 * points[ to_vertex(he ) ]
                + c2 * points[ to_vertex(next_halfedge(he)) ];
            ++i;
        }
        if (normalize)
            limits[*vit] = t1.normalized();
        else
            limits[*vit] = t1;
    }
    return limits;
}


namespace {
// Flat index of a tensor index (\a i,\a a j) of a gir of size \a sz
// per direction, where (i,j) is given in coordinates rotated by \a s
// quadrants (ie. s=0 are the original coordinates)
inline index_t face_pt_idx(index_t i, index_t j, index_t s, index_t sz)
{
    switch (s)
    {
    case 0:
        return sz * j + i;
    case 1:
        return sz * i + (sz-1-j);
    case 2:
        return sz * (sz-1-j) + (sz-1-i);
    case 3:
        return sz * (sz-1-i) + j;
    default:
        GISMO_ERROR("idx error");
    }
}
}

gsMultiPatch<real_t> gsSurfMesh::cc_acc3(bool comp_topology) const
{
    auto points = get_vertex_property<Point>("v:point");
    gsMultiPatch<real_t> mp;
    gsKnotVector<> kv(0,1,0,4);//cubic degree
    gsTensorBSplineBasis<2> bb(kv,kv);
    gsMatrix<real_t> coefs;
    gsSurfMesh::Halfedge h2;
    gsSurfMesh::Vertex v;
    real_t n;
#   pragma omp parallel for default(none) shared(std::cout,points,bb) private(n,v,h2,coefs) shared(mp)
    for (auto fit = faces_begin(); fit < faces_end(); ++fit)
    {
        //gsInfo << "face id: "<< fit->idx() <<"\n"; 
        coefs.resize(16,3);//thread privates must be initialized for each thread
        index_t s = 0;

        for ( auto he : halfedges(*fit) )
        {
            v = from_vertex(he);
            n = valence(v);
            auto c00 = coefs.row(face_pt_idx(0,0,s,4)).transpose();
            auto c11 = coefs.row(face_pt_idx(1,1,s,4)).transpose();
            auto c10 = coefs.row(face_pt_idx(1,0,s,4)).transpose();
            auto c01 = coefs.row(face_pt_idx(0,1,s,4)).transpose();
            h2 = opposite_halfedge(prev_halfedge(he));
            if (is_boundary(v))
            {
                c10 = ( 2*points[v] + points[to_vertex(he)] ) / 3;
                c01 = ( 2*points[v] + points[to_vertex(h2)] ) / 3;

                if (n>2)
                {
                    h2 = he; //find boundary edge
                    while ( !is_boundary(opposite_halfedge(h2)) )
                        h2 = cw_rotated_halfedge(h2);
                    c00 = (4*points[v] + points[to_vertex(h2)] +
                           points[to_vertex(next_halfedge((opposite_halfedge(h2))))] )/6;
                }
                else//single face at corner
                {
                    c00 = points[v];
                    n=3; // for c11
                }

                c11 = 2*(n-1) * points[v] +
                    2 * points[to_vertex(he)] +
                    points[to_vertex(next_halfedge(he))] +
                    2 * points[to_vertex(ccw_rotated_halfedge(he))];
                c11 /= 2*n+3;
            }
            else // interior vertex
            {
                c00 = n*n*points[v];
                for (auto h : halfedges(v))
                    c00 += 4 * points[ to_vertex(h) ] +
                        points[ to_vertex(next_halfedge(h)) ] ;
                c00 /= n*(n+5);
                c11 = n * points[v] +
                    2 * points[to_vertex(he)] +
                    points[to_vertex(next_halfedge(he))] +
                    2 * points[to_vertex(ccw_rotated_halfedge(he))];
                c11 /= n+5;
                c10 = n * points[v]
                    + 2 * points[to_vertex(he)]+
                    points[to_vertex(ccw_rotated_halfedge(he))]+
                    points[to_vertex(cw_rotated_halfedge(he)) ]+
                    0.5 * points[to_vertex(next_halfedge(he))] +
                    0.5 * points[to_vertex(next_halfedge(cw_rotated_halfedge(he))) ];
                c10 /= n+5;
                c01 = n * points[v]
                    + 2 * points[to_vertex(h2)]+
                    points[to_vertex(ccw_rotated_halfedge(h2))]+
                    points[to_vertex(cw_rotated_halfedge(h2)) ]+
                    0.5 * points[to_vertex(next_halfedge(h2))] +
                    0.5 * points[to_vertex(next_halfedge(cw_rotated_halfedge(h2))) ];
                c01 /= n+5;
            }

            /*
            h2 = opposite_halfedge(prev_halfedge(he));
            n = valence(v);
            if ( is_boundary( opposite_halfedge(he) ) )
            {
                c10 = ( 2*points[v] + points[to_vertex(he)] ) / 3;
            }
            else //interior edge
            {
                c10 = n * points[v] + 2 * points[to_vertex(he)] +
                    points[to_vertex(ccw_rotated_halfedge(he))]+
                    points[to_vertex(cw_rotated_halfedge(he)) ]+
                    0.5 * points[to_vertex(next_halfedge(he))] +
                    0.5 * points[to_vertex(next_halfedge(cw_rotated_halfedge(he))) ];
                c10 /= n+5;
            }

            h2 = opposite_halfedge(prev_halfedge(he));
            //if ( is_boundary(h2) )
            if (is_boundary(v))
            {
                c01 = ( 2 * points[v] + points[to_vertex(h2)] ) / 3;
            }
            else //interior edge
            {
                c01 = n * points[v] + 2 * points[to_vertex(h2)] +
                    points[to_vertex(ccw_rotated_halfedge(h2))]+
                    points[to_vertex(cw_rotated_halfedge(h2)) ]+
                    0.5 * points[to_vertex(next_halfedge(h2))] +
                    0.5 * points[to_vertex(next_halfedge(cw_rotated_halfedge(h2))) ];
                c01 /= n+5;
            }
            */

            ++s;//next halfedge
        }

#       pragma omp critical (mp_addPatch)
        mp.addPatch( bb.makeGeometry(give(coefs)) );
    }
    if (comp_topology)
        mp.computeTopology();
    return mp;
}


gsMultiPatch<real_t> gsSurfMesh::linear_patches() const
{
    gsMultiPatch<> mp;
    gsSurfMesh HEmesh(*this);

    // Counts for each vertex the number of passes
    auto vpassed = HEmesh.vertex_property<gsSurfMesh::Scalar>("v:passed", 0 );
    // Index of the curve loop (negative) or of the patch (positive)
    auto hindex = HEmesh.halfedge_property<int>("h:index", 0 );
    // Patch index of each face (seems unused)
    //auto findex = HEmesh.face_property<int>("f:index", -1 );

    // Create a stack of EVs and boundary EVs
    std::list<gsSurfMesh::Vertex> EVs;
    for (auto v : HEmesh.vertices())
    {
        if (HEmesh.valence(v) != 4 && !HEmesh.is_boundary(v))//interior
            EVs.push_back(v);
        else if (HEmesh.valence(v) > 3 && HEmesh.is_boundary(v)) //boundary
            EVs.push_back(v);
        // else if (HEmesh.valence(v) == 2 && HEmesh.is_boundary(v)) //boundary corner
        //     EVs.push_back(v);
    }

    // For all EVs, find the curve-loop over all adjacent edges
    index_t curveloopIdx = -1;
    for (auto EV = EVs.begin(); EV!=EVs.end(); EV++)
    {
        // std::list<gsSurfMesh::Vertex>::iterator EV = EVs.begin();
        for ( auto he : HEmesh.halfedges(*EV) ) // iterate over all HE that come from a EV
        {
            auto h = he;
            // Mark the vertex from which we departed as passed
            auto v = HEmesh.from_vertex(h);
            auto vold = HEmesh.from_vertex(h);
            vpassed[v] += 1;

            // stopping conditions:
            // - boundary is hit
            // - HE is already assigned to another loop
            // - hit another V that has been crossed,
            while (true)
            {
                v = HEmesh.to_vertex(h);
                vold = HEmesh.from_vertex(h);
                // If h is already assigned to a curve loop, we stop
                if (hindex[h]!=0)
                    break;

                hindex[h] =
                    hindex[HEmesh.opposite_halfedge(h)] = curveloopIdx;
                if (!HEmesh.is_boundary(v))
                {
                    if (HEmesh.valence(v)==4) // interior ordinary vertex
                    {
                        h = HEmesh.next_halfedge(h);
                        h = HEmesh.opposite_halfedge(h);
                        h = HEmesh.next_halfedge(h);
                    }
                    else
                        break; // EV, thus stop this curveloop
                }
                else // is boundary
                {
                    if (HEmesh.valence(v) > 3)
                        break; // EV, thus stop this curveloop
                }
                // Mark the to-vertex as passed, if it is not an EV
                vpassed[v] += 1;

                // Check if the to-vertex is a boundary vertex. If yes, we stop
                if (HEmesh.is_boundary(v))
                    break;
            }
            curveloopIdx--;
        }
    }

    // gsWriteParaview(HEmesh,"HEmesh",{"v:passed"});
    // Collect intersections and EVs (points that have been passed more than once)
    EVs.clear();
    EVs.resize(0);
    for (auto v : HEmesh.vertices())
        if (vpassed[v]>1) // Intersection or EV
            EVs.push_back(v);
    // Probably not needed:
    // else if (vpassed[v]==1 && HEmesh.is_boundary(v)) // Point that splits the boundary
    //     EVs.push_back(v);

    // Patch index which will be assigned to half-edges, starting from 1 now, because 0 is reserved for non-assigned half edges
    index_t patchIdx = 1;
    // Stores the number of element in both directions for each patch
    std::vector<std::pair<index_t,index_t>> dirSizes;
    // From each starting point, we take each half edge again and we assign patch indices to the incident faces
    for (auto EV = EVs.begin(); EV!=EVs.end(); EV++)
    {
        for ( auto he : HEmesh.halfedges(*EV) ) // iterate over all HE that come from a EV
        {
            if (!HEmesh.is_valid(he))
                continue;
            // Check if the half-edge is already assigned to a patch
            if (hindex[he] > 0)
                continue;

            auto h = he;
            if (HEmesh.is_boundary(h))
                continue;
            // Store other patch direction sizes for a check
            std::pair<index_t,index_t> dirSize(1,1);
            bool dir = 0;
            // Stopping conditions:
            // 1. Startpoint reached
            // 2. The half-edge h has a positive index, meaning it is assigned to a patch.
            while (true)
            {
                if (hindex[h] > 0)
                {
                    gsWarn<<"Half-edge has already been passed\n";
                    break;
                }
                // ---- Assign the patch index to the half-edge
                hindex[h] = patchIdx;
                // ---- Count on the direction
                if (dir)
                    dirSize.second++;
                else
                    dirSize.first++;

                // If h goes to the original EV, we stop
                if (HEmesh.to_vertex(h)==*EV)
                    break;

                // Move on (goes around the corner on the same face)
                h = HEmesh.next_halfedge(h);

                // Check if the next half-edge has a negative index OR is a boundary OR if the opposite has a face with patchIdx.
                // If so, it is a patch boundary and we ended up in a corner
                // ---- Go around the corner
                if (hindex[h]<0
                    || ( ( HEmesh.is_boundary(h) || HEmesh.is_boundary(HEmesh.opposite_halfedge(h)) ) && HEmesh.is_boundary(HEmesh.to_vertex(h)) )
                    )
                {
                    // if dir was 1, pushback dirSize and empty temporary one
                    if (dir==1)
                    {
                        dirSizes.push_back(dirSize);
                        dirSize.first = 1;
                        dirSize.second= 1;
                    }
                    dir = !dir;
                    continue;
                }
                // ---- Stopping condition 2
                else if (hindex[h] > 0)
                {
                    gsWarn<<"Something went wrong. Half-edge is already assigned\n";
                    break;
                }
                // check if the next half-edge has an index of 0. If so, it is an interior half-edge
                // ---- Go straight
                else if (hindex[h]==0)
                {
                    h = HEmesh.opposite_halfedge(h);
                    h = HEmesh.next_halfedge(h);
                    continue;
                }
                else
                    gsWarn<<"Something went wrong?\n";
            }

            // Make a linear tensor basis of size dir.first,dir.second
            gsKnotVector<real_t> KV0(0,dirSize.first -1,dirSize.first -2,2,1,1);
            gsKnotVector<real_t> KV1(0,dirSize.second-1,dirSize.second-2,2,1,1);
            gsTensorBSplineBasis<2,real_t> tbasis(KV0,KV1);
            gsMatrix<> coefs(dirSize.first*dirSize.second,3);
            coefs.setZero();

            // Start from the EV again. We start in dirsize
            h = he;
            if (HEmesh.is_boundary(h))
                continue;
            if (!HEmesh.is_valid(he))
                continue;

            auto h0 = h;
            if (HEmesh.is_boundary(h0))
                gsWarn<<"Edge is already a boundary...\n";
            bool flip = false;
            gsMatrix<> tmpcoefs(dirSize.first,3);
            for (index_t i=0; i!=dirSize.second; i++)
            {
                for (index_t j=0; j!=dirSize.first-2; j++) //
                {
                    tmpcoefs.row(j) = HEmesh.position(HEmesh.from_vertex(h)).transpose();

                    h = HEmesh.next_halfedge(h);
                    if (!HEmesh.is_boundary(h))
                    {
                        h = HEmesh.opposite_halfedge(h);
                        h = HEmesh.next_halfedge(h);
                    }
                }
                tmpcoefs.row(dirSize.first-2) = HEmesh.position(HEmesh.from_vertex(h)).transpose();
                tmpcoefs.row(dirSize.first-1) = HEmesh.position(HEmesh.to_vertex(h)).transpose();

                if (flip)
                    h = HEmesh.opposite_halfedge(h);

                h = HEmesh.next_halfedge(h);
                h = HEmesh.next_halfedge(h);

                if (flip)
                {
                    h = HEmesh.opposite_halfedge(h);
                    tmpcoefs = tmpcoefs.colwise().reverse().eval();
                }
                coefs.block(i*dirSize.first,0,dirSize.first,3) = tmpcoefs;

                flip = !flip;
            }
            mp.addPatch(tbasis.makeGeometry(give(coefs)));
            patchIdx++;
        }
    }

    /*
      for (size_t k = 0; k!=dirSizes.size(); k++)
      gsDebug<<"Patch "<<k<<": Dir 0: "<<dirSizes[k].first<<"; Dir 1: "<<dirSizes[k].second<<"\n";
*/
    return mp;
}

namespace internal {


void gsXml<gsSurfMesh>::get_into(gsXmlNode * node, gsSurfMesh & result)
{
    assert( ( !strcmp( node->name(),"SurfMesh") )
            &&  ( !strcmp(node->first_attribute("type")->value(),"off") ) );

    result = gsSurfMesh();

    /*
      if ( !strcmp(node->first_attribute("type")->value(),"off") )
      {
      read_off_ascii(result,node->value());
      return;
      }
    */

    // !strcmp(node->first_attribute("type")->value(),"poly")
    // !strcmp(node->first_attribute("type")->value(),"stl")
    //!strcmp(node->first_attribute("type")->value(),"obj")
    //!strcmp(node->first_attribute("type")->value(),"vtk")


    std::istringstream str;
    str.str( node->value() );

    unsigned nv  = atoi ( node->first_attribute("vertices")->value() ) ;
    unsigned nf  = atoi ( node->first_attribute("faces")->value() ) ;
    unsigned ne  = atoi ( node->first_attribute("edges")->value() ) ;
    result.reserve(nv, std::max(3*nv, ne), nf);
    real_t x(0), y(0), z(0); // T?
    for (unsigned i=0; i<nv; ++i)
    {
        gsGetReal(str, x);
        gsGetReal(str, y);
        gsGetReal(str, z);
        result.add_vertex(Point(x,y,z));
    }

    unsigned k, c = 0;
    std::vector<gsSurfMesh::Vertex> face;
    for (unsigned i=0; i<nf; ++i)
    {
        gsGetInt(str, c);
        face.resize(c);
        for (unsigned j=0; j<c; ++j)
        {
            gsGetInt(str, k);
            face[j] = gsSurfMesh::Vertex(k);
        }
        result.add_face(face);
    }
}

gsXmlNode *
gsXml<gsSurfMesh>::put (const gsSurfMesh & obj, gsXmlTree & data)
{

    return nullptr;
};

}//namespace internal

//=============================================================================
} // namespace gismo
//=============================================================================
