/** @file gsSurfMesh.hpp

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
#include <gsIO/gsXml.h>

#include <algorithm>
#include <fstream>
#include <map>

namespace gismo {

template <class Scalar>
gsSurfMesh<Scalar>::
gsSurfMesh() : Base()
{
    // the topology properties are allocated by the base class constructor;
    // here we only add the geometry
    vpoint_   = add_vertex_property<Point>("v:point",Point(0,0,0));
}

template <class Scalar>
gsSurfMesh<Scalar>::
gsSurfMesh(const gsMatrix<Scalar> & pts) : gsSurfMesh()
{
    for( auto & col : pts.colwise() )
        this->add_vertex(col);
}

template <class Scalar>
gsSurfMesh<Scalar>::
~gsSurfMesh()
{
}

template <class Scalar>
gsSurfMesh<Scalar>&
gsSurfMesh<Scalar>::
operator=(const gsSurfMesh<Scalar>& rhs)
{
    if (this != &rhs)
    {
        // deep copy of the topology (property containers, connectivity, counters)
        Base::operator=(rhs);

        // property handles contain pointers, have to be reassigned
        vpoint_   = get_vertex_property<Point>("v:point");

        // normals might be there, therefore use get_property
        vnormal_  = get_vertex_property<Point>("v:normal");
        fnormal_  = get_face_property<Point>("f:normal");
    }

    return *this;
}

template <class Scalar>
gsSurfMesh<Scalar>&
gsSurfMesh<Scalar>::
assign(const gsSurfMesh<Scalar>& rhs)
{
    if (this != &rhs)
    {
        // clears the property containers and re-creates the topology properties
        Base::assign(rhs);

        // re-create and copy the geometry
        vpoint_   = add_vertex_property<Point>("v:point",Point(0,0,0));
        vpoint_.array()    = rhs.vpoint_.array();

        // normals might be there, therefore use get_property
        vnormal_  = get_vertex_property<Point>("v:normal");
        fnormal_  = get_face_property<Point>("f:normal");
    }

    return *this;
}

template <class Scalar>
gsSurfMesh<Scalar>&
gsSurfMesh<Scalar>::
operator=(gsSurfMesh<Scalar>&& rhs) noexcept
{
    if (this != &rhs)
    {
        Base::operator=(std::move(rhs));

        // property handles contain pointers, have to be reassigned
        vpoint_  = get_vertex_property<Point>("v:point");

        // normals might be there, therefore use get_property
        vnormal_ = get_vertex_property<Point>("v:normal");
        fnormal_ = get_face_property<Point>("f:normal");
    }
    return *this;
}


template <class Scalar>
bool
gsSurfMesh<Scalar>::
read(const std::string& filename)
{
    // extract file extension
    std::string::size_type dot(filename.rfind("."));
    if (dot == std::string::npos) return false;
    std::string ext = filename.substr(dot+1, filename.length()-dot-1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    // extension determines reader
    if (ext == "off")
    {
        return read_off<Scalar>(*this, filename);
    }
    else if (ext == "obj")
    {
        std::ifstream in(filename.c_str());
        if (!in) return false;
        return read_obj<Scalar>(*this, in);
    }
    else if (ext == "stl")
    {
        return read_stl<Scalar>(*this, filename);
    }

    // we didn't find a reader module
    return false;
}

template <class Scalar>
bool
gsSurfMesh<Scalar>::
write(const std::string& filename) const
{
    return write_mesh(*this, filename);
}

template <class Scalar>
typename gsSurfMesh<Scalar>::Vertex
gsSurfMesh<Scalar>::
add_vertex(const Point& p)
{
    Vertex v = new_vertex();
    vpoint_[v] = p;
    return v;
}

template <class Scalar>
void
gsSurfMesh<Scalar>::
update_face_normals()
{
    if (!fnormal_)
        fnormal_ = face_property<Point>("f:normal",Point(0,0,0));

    Face_iterator fit, fend=faces_end();

    for (fit=faces_begin(); fit!=fend; ++fit)
        fnormal_[*fit] = compute_face_normal(*fit);
}

template <class Scalar>
typename gsSurfMesh<Scalar>::Normal
gsSurfMesh<Scalar>::compute_face_normal(Face f) const
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
        return p2.cross(p0).normalized();
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

template <class Scalar>
void
gsSurfMesh<Scalar>::
update_vertex_normals()
{
    if (!vnormal_)
        vnormal_ = vertex_property<Point>("v:normal",Point(0,0,0));

    Vertex_iterator vit, vend=vertices_end();

    for (vit=vertices_begin(); vit!=vend; ++vit)
        vnormal_[*vit] = compute_vertex_normal(*vit);
}


template <class Scalar>
typename gsSurfMesh<Scalar>::Normal
gsSurfMesh<Scalar>::compute_vertex_normal(Vertex v) const
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

template <class Scalar>
Scalar
gsSurfMesh<Scalar>::
edge_length(Edge e) const
{
    return (vpoint_[vertex(e,0)] - vpoint_[vertex(e,1)]).norm();
}

template <class Scalar>
void gsSurfMesh<Scalar>::quad_split()
{
    typename gsSurfMesh<Scalar>::Vertex v;
    typename gsSurfMesh<Scalar>::Halfedge he;

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
        v = typename gsSurfMesh<Scalar>::Vertex(fnv + (i++));//face vertex id ?
        //Start from an original vertex
        auto fv = vertices(fit).begin();
        if ((*fv).idx() >= env) ++fv; //todo: add -> operator
        //assert ( (*fv).idx() < nv )
        quad_split(fit, v, fv.he());
    }

}

template <class Scalar>
inline memory::unique_ptr<gsGeometry<Scalar> > gsSurfMesh<Scalar>::asPatch(typename gsSurfMesh<Scalar>::Halfedge h, int deg) const
{
    static gsKnotVector<Scalar> kv;
    kv.initUniform(0, 1, 0, 1, 1, deg);
    static gsTensorBSplineBasis<2> bs;
    bs = { kv,kv };
    static gsMatrix<Scalar> coefs;
    coefs.resize( (deg+1)*(deg+1), 3);

    Halfedge hh;
    index_t c = 0;

    // Finding the first point and take its halfedge
    for (index_t i = 0; i<(deg/2); ++i)
    {
        h = backward_halfedge(h);
        h = next_halfedge(next_halfedge(opposite_halfedge(h)));
    }

    // Going through the control points for each row
    for (int j = 0; j<deg; ++j)
    {
        coefs.row(c++) = vpoint_[from_vertex(h)];
        coefs.row(c++) = vpoint_[  to_vertex(h)];
        hh = h;
        for (int i = 1; i<deg; ++i)
        {
            hh = forward_halfedge(hh);
            coefs.row(c++) = vpoint_[to_vertex(hh)];
        }
        h = ccw_rotated_halfedge((prev_halfedge(h)));
    }
    //last row
    coefs.row(c++) = vpoint_[from_vertex(h)];
    coefs.row(c++) = vpoint_[  to_vertex(h)];
    hh = h;
    for (int i = 1; i<deg; ++i)
    {
        hh = opposite_halfedge( prev_halfedge( opposite_halfedge(
             prev_halfedge( opposite_halfedge(hh) ))));
        coefs.row(c++) = vpoint_[to_vertex(hh)];
    }
    return bs.makeGeometry(coefs);
}

template <class Scalar>
gsMultiPatch<Scalar> gsSurfMesh<Scalar>::asSpline(int deg) const
{
    gsMultiPatch<Scalar> res;
    int n;
    Halfedge he, hh;
    Vertex ve;
    bool evface;
    if ( 0 == deg%2 )
        for (auto v : vertices())
        {
            if (is_boundary(v)) continue; // handling boundaries
            n = valence(v);

            he = halfedge(v);
            
            evface = false;
            hh = he;
            // exluding ev faces
            do
            {
                if (valence(face(hh)) != 4)
                {
                    evface = true;
                    break;
                }
                hh = ccw_rotated_halfedge(hh);
            } while (hh != he);

            if (evface)
                continue;

            if (4==n)  // exluding ev points
                res.addPatch( asPatch(he, deg) );
        }
    else // 1 == deg%2
        for (auto f : faces())
        {
            if (valence(f) != 4) continue;

            evface = false;
            if (is_boundary(f)) continue; // handling boundaries

            for ( auto v : vertices(f) )
            {
                n = valence(v);
                                
                if ( n!=4 ) break;

                he = halfedge(v);
                hh = he;

                // exluding ev faces
                do
                {
                    if (valence(face(hh)) != 4)
                    {
                        evface = true;
                        break;
                    }
                    hh = ccw_rotated_halfedge(hh);
                } while (hh != he);

                if (evface) break;
            }

            if (evface) continue;

            if (4==n)
                res.addPatch( asPatch(halfedge(f), deg) );
        }

    return res;
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

template <class Scalar>
gsMultiPatch<Scalar> gsSurfMesh<Scalar>::linear_patches() const
{
    gsMultiPatch<Scalar> mp;
    gsSurfMesh<Scalar> HEmesh(*this);

    // Counts for each vertex the number of passes
    auto vpassed = HEmesh.vertex_property<Scalar>("v:passed", 0 );
    // Index of the curve loop (negative) or of the patch (positive)
    auto hindex = HEmesh.halfedge_property<int>("h:index", 0 );
    // Patch index of each face (seems unused)
    //auto findex = HEmesh.face_property<int>("f:index", -1 );

    // Create a stack of EVs and boundary EVs
    std::list<typename gsSurfMesh<Scalar>::Vertex> EVs;
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
        // std::list<typename gsSurfMesh<Scalar>::Vertex>::iterator EV = EVs.begin();
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
            gsKnotVector<Scalar> KV0(0,dirSize.first -1,dirSize.first -2,2,1,1);
            gsKnotVector<Scalar> KV1(0,dirSize.second-1,dirSize.second-2,2,1,1);
            gsTensorBSplineBasis<2,Scalar> tbasis(KV0,KV1);
            gsMatrix<Scalar> coefs(dirSize.first*dirSize.second,3);
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
            gsMatrix<Scalar> tmpcoefs(dirSize.first,3);
            for (index_t i=0; i!=dirSize.second; i++)
            {
                for (index_t j=0; j!=dirSize.first-2; j++) //
                {
                    tmpcoefs.row(j) = HEmesh.position(HEmesh.from_vertex(h)).transpose().template cast<Scalar>();

                    h = HEmesh.next_halfedge(h);
                    if (!HEmesh.is_boundary(h))
                    {
                        h = HEmesh.opposite_halfedge(h);
                        h = HEmesh.next_halfedge(h);
                    }
                }
                tmpcoefs.row(dirSize.first-2) = HEmesh.position(HEmesh.from_vertex(h)).transpose().template cast<Scalar>();
                tmpcoefs.row(dirSize.first-1) = HEmesh.position(HEmesh.to_vertex(h)).transpose().template cast<Scalar>();

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


template <class Scalar>
typename gsSurfMesh<Scalar>::Point
gsSurfMesh<Scalar>::face_barycenter(Face f)
{
    unsigned int f_val{ 0 };
    Point coords, tmp;

    f_val = valence(f);
    coords.setZero();
    for (auto vit : vertices(f)) {
        coords += position(vit);
    }
    coords = coords / static_cast<Scalar>(f_val);
    tmp[0] = coords(0);
    tmp[1] = coords(1);
    tmp[2] = coords(2);

    return tmp;

}

    
template <class Scalar>
void gsSurfMesh<Scalar>::mergeDoubleVertices()
{
    GISMO_NO_IMPLEMENTATION
}

template <class Scalar>
gsVector<typename gsSurfMesh<Scalar>::Vertex>
gsSurfMesh<Scalar>::add_mesh(gsSurfMesh<Scalar>& subMesh)
{

    gsVector<Vertex> idmap(subMesh.n_vertices()); // local vertex of subMesh mapping with global to new mesh

    // Adding vertices to current mesh
    for (auto vit : subMesh.vertices())
        idmap[vit.idx()] = this->add_vertex(subMesh.position(vit));

    std::vector<Vertex> vv;
    // Adding faces to current mesh
    for (auto fit : subMesh.faces())
    {
        vv.clear();
        for (auto vit : subMesh.vertices(fit))
           vv.push_back(idmap[vit.idx()]);

        this->add_face(vv);
    }

    return idmap;
}

template <class Scalar>
void gsSurfMesh<Scalar>::quad_split(index_t w)
{

    if (w==0) // original faces (dummy)
    {
        return;
    }
    else if (w == 1) // uniform split at half of each edge
    {
        quad_split();
        return;
    }
    else // general cases for w >=2
    {

        GISMO_ASSERT(w < 3, "NOT TESTED for w>=3!");

        Vertex v, vs, ve;
        Halfedge he, hh, hb;
        // reserve vertices, edges, faces
        reserve(n_vertices() + n_edges() + n_faces(),
            2 * n_edges(), 4 * n_faces());


        gsSurfMesh<Scalar> nm;

        // loop over all edges, add edge points
        Point tmp, dx, tmpA, tmpB, tmpC;

        Vertex_property<bool> master_verts;

        for (auto vit : vertices())
            nm.add_vertex(position(vit));

        if (!vprops_.has("v:master"))
            master_verts = nm.add_vertex_property<bool>("v:master", false);
        else
            master_verts = nm.get_vertex_property<bool>("v:master");


        if (!nm.vprops_.has("v:neighval"))
            nm.add_vertex_property<int>("v:neighval", 4);

        for (auto eit : edges())
        {

            he = halfedge(eit, 0);
            dx = (position(to_vertex(he)) - position(from_vertex(he))) / (Scalar)(w + 1);
            tmp = position(from_vertex(he));
            master_verts[from_vertex(he)] = true;

            hh = he;
            for (index_t i = 0; i < w; i++)
            {
                tmp += dx;
                nm.add_vertex(tmp);
                hh = prev_halfedge(opposite_halfedge(insert_vertex(hh, add_vertex(tmp))));
                hh = next_halfedge(hh);
            }


        }

        auto points = get_vertex_property<Point>("v:point");
        std::vector<Vertex> intverts;
        index_t n = 0, count = 0;
        for (auto fit : faces())
        {


            // Find one master vertex in the face
            for (auto hit : halfedges(fit))
                if (master_verts[from_vertex(hit)])
                {
                    he = hit;
                    break;
                }

            // Count the number of initial edges
            hh = he;
            count = 0;
            do
            {
                count++;
                hh = next_halfedge(hh);
            } while (he != hh);

            // Create corner subfaces

            n = count / (w + 1); // # of master edges (initial ones)
            count = 0;
            hh = he;
            hb = prev_halfedge(prev_halfedge(he));
            intverts.clear();
            do
            {
                count++;
                tmpA = points[from_vertex(hh)];
                tmpB = points[to_vertex(hh)];
                tmpC = points[from_vertex(prev_halfedge(hh))];
                tmp = (tmpB - tmpA) + (tmpC - tmpA) + tmpA;
                v = nm.add_vertex(tmp);
                intverts.push_back(v);
                nm.add_quad(Vertex(to_vertex(hh)), Vertex(from_vertex(hh)), Vertex(from_vertex(prev_halfedge(hh))), v);
                hh = next_halfedge(hh);
                if (count < n)
                {
                    hh = next_halfedge(hh);
                    hh = next_halfedge(hh);

                }
            } while (hb != hh);

            // Create interior face(s)
            hh = he;
            count = 0;
            do
            {

                hh = next_halfedge(hh);
                vs = intverts[count];
                if (count == n - 1)
                    ve = intverts[0];
                else
                    ve = intverts[count + 1];
                nm.add_quad(Vertex(from_vertex(hh)), vs, ve, Vertex(to_vertex(hh)));
                hh = next_halfedge(hh);
                hh = next_halfedge(hh);
                count++;
            } while (count != n);



            nm.add_quad(intverts[0], intverts[3], intverts[2], intverts[1]);
        }

        *this = std::move(nm);

    }
}

template <class Scalar>
void gsSurfMesh<Scalar>::dual_mesh_inplace()
{
    // Dual-mesh instance
    Self dm = dual_mesh();

    *this = std::move(dm);

}

template <class Scalar>
typename gsSurfMesh<Scalar>::Self gsSurfMesh<Scalar>::dual_mesh()
{
    // Dual-mesh instance
    gsSurfMesh<Scalar> dm;

    //Instances of the original mesh
    typename gsSurfMesh<Scalar>::Vertex v;
    typename gsSurfMesh<Scalar>::Face f;
    typename gsSurfMesh<Scalar>::Face fop;
    typename gsSurfMesh<Scalar>::Edge e;


    // Calculate the dual vertices (from barycenter of original faces)

    std::map<Face, Vertex> FVMap;

    // For each face take the barycenter
    for (auto fit : faces()) {
        v = dm.add_vertex(face_barycenter(fit));
        FVMap[fit] = v;
    }

    // For the connected vertices in the original mesh create the dual faces
    std::vector<Vertex> df;
    for (auto vit : vertices()) {
        if (is_boundary(vit)) { continue; }
        df.clear();
        for (auto fit : faces(vit)) {
            df.push_back(FVMap[fit]);
        }
        dm.add_face(df);


    }



    return dm;

}

template <class Scalar>
real_t gsSurfMesh<Scalar>::
angle(Halfedge h1, Halfedge h2)
{
    real_t result = 0.0;
    Point v1 = position(to_vertex(h1)) - position(from_vertex(h1));
    Point v2 = position(to_vertex(h2)) - position(from_vertex(h2));
    result = static_cast<real_t>(math::acos(v1.dot(v2)/(v1.norm()*v2.norm())));

    return result;
}

template <class Scalar>
void gsSurfMesh<Scalar>::
display_halfedge()
{
    Point tmp;
    auto hpp = add_vertex_property<Point>("v:halfedge", tmp.setZero());
    Halfedge he;
    for (auto fit : faces())
    {
        he = halfedge(fit);
        hpp[from_vertex(he)] = (position(to_vertex(he)) -
                position(from_vertex(he))).normalized();
    }
}

template <class Scalar>
typename gsSurfMesh<Scalar>::Self gsSurfMesh<Scalar>::flip_orientation()
{

    // Build flipped face connectivity
    std::vector<std::vector<Vertex>> F;
    F.reserve(n_faces());

    std::vector<Vertex> face;

    for (auto fit : faces())
    {
        face.clear();
        for (auto fv : vertices(fit))
            face.push_back(fv);

        // Flip orientation (reverse vertices in the face)
        std::reverse(face.begin() + 1, face.end());

        F.push_back(face);
    }

    // Reconstruct the mesh
    Self out;
    for (auto vit : vertices())
    {
        out.add_vertex(position(vit));
    }
    for (auto vec : F)
    {
        out.add_face(vec);
    }

    //////// Recompute normals
    out.update_face_normals();
    // out.update_vertex_normals();

    return out;
}

template <class Scalar>
void gsSurfMesh<Scalar>::polyhedral_modification_boundary()
{
    // Current implementation only for regular boundary (vertex valence = 3) and
       // conrners (vertex valence = 2).
       // TODO: General case for EF in boundary by using Chebysev points (see A.Nashri 1987).

    std::map<Vertex, Point> bvmap; // New positions for boundary vertices
    Vertex bv;
    auto pts = points();
    // Compute the new positions for boundary vertices.
    for (auto hit : halfedges())
    {
        if (touches_boundary(hit))
        {
            bv = from_vertex(hit);

            if (valence(bv) == 3) // Regular boundary case
            {
                bvmap[bv] = 2 * pts[bv] - pts[from_vertex(prev_halfedge(hit))];
            }
            else if (valence(bv) == 2) // Corner boundary case
            {
                bvmap[bv] = 4 * pts[bv] - 2 * pts[from_vertex(prev_halfedge(hit))]
                    - 2 * pts[to_vertex(hit)] + pts[to_vertex(next_halfedge(hit))];
            }
            else // irregular case
            {
                gsWarn << "Irregular boundary stop process\n";
                return;
            }

        }
    }

    // Modify mesh boundary
    for (auto vit : vertices())
        if (is_boundary(vit))
            position(vit) = bvmap[vit];
}


namespace internal
{

template <class Scalar>
void gsXml< gsSurfMesh<Scalar> >::get_into(gsXmlNode * node, gsSurfMesh<Scalar> & result)
{
    GISMO_ASSERT( !strcmp( node->name(),"SurfMesh") || !strcmp( node->name(),"Mesh"),
                  "gsXml<gsSurfMesh>: expected a <Mesh> or <SurfMesh> node, got <"
                  << node->name() << ">." );

    result = gsSurfMesh<Scalar>();

    // Readers that hand the raw file through verbatim tag the node with a
    // "format" attribute and store the file contents as the node value.
    gsXmlAttribute * fmt = node->first_attribute("format");
    if ( NULL != fmt )
    {
        if ( !strcmp(fmt->value(),"off") )
        {
            read_off_ascii(result, node->value());
            return;
        }
        if ( !strcmp(fmt->value(),"obj") )
        {
            std::istringstream is( node->value() );
            read_obj(result, is);
            return;
        }
        GISMO_ERROR("gsXml<gsSurfMesh>: unsupported mesh format \""
                    << fmt->value() << "\".");
    }

    // Otherwise the node carries the counts as attributes and the coordinates
    // and face indices as its value.
    gsXmlAttribute * av = node->first_attribute("vertices");
    gsXmlAttribute * af = node->first_attribute("faces");
    gsXmlAttribute * ae = node->first_attribute("edges");
    GISMO_ASSERT( av && af && ae, "gsXml<gsSurfMesh>: <Mesh> node has neither a "
                  "\"format\" attribute nor \"vertices\"/\"faces\"/\"edges\" counts." );

    std::istringstream str;
    str.str( node->value() );

    unsigned nv  = atoi ( av->value() ) ;
    unsigned nf  = atoi ( af->value() ) ;
    unsigned ne  = atoi ( ae->value() ) ;
    result.reserve(nv, std::max(3*nv, ne), nf);
    real_t x(0), y(0), z(0); // reading as real_t, cast to Scalar below if needed
    for (unsigned i=0; i<nv; ++i)
    {
        gsGetReal(str, x);
        gsGetReal(str, y);
        gsGetReal(str, z);
        result.add_vertex(typename gsSurfMesh<Scalar>::Point(static_cast<Scalar>(x),static_cast<Scalar>(y),static_cast<Scalar>(z)));
    }

    /* //Alternative for reading quads only (with complex topolog)
   unsigned k, c = 0;
    std::vector<typename gsSurfMesh<Scalar>::Vertex> face(4);
    std::vector<typename gsSurfMesh<Scalar>::Edge> e(4);
    for (unsigned i=0; i<nf; ++i)
    {
        gsGetInt(str, c);
        GISMO_ASSERT(4==c, "quads?");
        for (unsigned j=0; j<c; ++j)
        {
            gsGetInt(str, k);
            face[j] = typename gsSurfMesh<Scalar>::Vertex(k);
        }
        for (unsigned j=0; j<c; ++j)
            e[j] = result.find_or_add_edge(face[j],face[(j+1)%c]);
        result.add_quad(e[0], e[1], e[2], e[3]);
    }
    //*/

    unsigned k, c = 0;
    std::vector<typename gsSurfMesh<Scalar>::Vertex> face;
    for (unsigned i=0; i<nf; ++i)
    {
        gsGetInt(str, c);
        face.resize(c);
        for (unsigned j=0; j<c; ++j)
        {
            gsGetInt(str, k);
            face[j] = typename gsSurfMesh<Scalar>::Vertex(k);
        }
        result.add_face(face);
    }

    if (0!=ne)
    {
        typename gsSurfMesh<Scalar>::template Halfedge_property<bool> sharp = result.template add_halfedge_property<bool>("h:sharp");
        face.resize(2);
        typename gsSurfMesh<Scalar>::Halfedge he;
        for(unsigned i = 0; i!=ne; ++i)
        {
            gsGetInt(str, k);
            face[0] = typename gsSurfMesh<Scalar>::Vertex(k);
            gsGetInt(str, k);
            face[1] = typename gsSurfMesh<Scalar>::Vertex(k);
            he = result.find_halfedge(face[0], face[1]);
            sharp[he] = true;
            he = result.opposite_halfedge(he);
            sharp[he] = true;
        }
    }
}

template <class Scalar>
gsXmlNode *
gsXml< gsSurfMesh<Scalar> >::put (const gsSurfMesh<Scalar> & obj, gsXmlTree & data)
{
    GISMO_UNUSED(obj);
    GISMO_UNUSED(data);
    return nullptr;
};

}//namespace internal

} // namespace gismo
