
#include <gsMesh2/IO.h>

#include <cstdio>


//== NAMESPACES ===============================================================


namespace gismo {


//== IMPLEMENTATION ===========================================================


// helper function
template <typename T> size_t read(FILE* in, T& t)
{
    return fread((char*)&t, 1, sizeof(t), in);
}


// helper function
template <typename T> size_t write(FILE* out, T& t)
{
    return fwrite((char*)&t, 1, sizeof(t), out);
}


//-----------------------------------------------------------------------------


bool read_poly(gsSurfMesh& mesh, const std::string& filename)
{
    // open file (in binary mode)
    FILE* in = fopen(filename.c_str(), "rb");
    if (!in) return false;


    // clear mesh
    mesh.clear();


    // how many elements?
    unsigned int nv, ne, nh, nf;
    read(in, nv);
    read(in, ne);
    read(in, nf);
    nh = 2*ne;


    // resize containers
    mesh.vprops_.resize(nv);
    mesh.hprops_.resize(nh);
    mesh.eprops_.resize(ne);
    mesh.fprops_.resize(nf);


    // get properties
    gsSurfMesh::Vertex_property<gsSurfMesh::Vertex_connectivity>      vconn = mesh.vertex_property<gsSurfMesh::Vertex_connectivity>("v:connectivity");
    gsSurfMesh::Halfedge_property<gsSurfMesh::Halfedge_connectivity>  hconn = mesh.halfedge_property<gsSurfMesh::Halfedge_connectivity>("h:connectivity");
    gsSurfMesh::Face_property<gsSurfMesh::Face_connectivity>          fconn = mesh.face_property<gsSurfMesh::Face_connectivity>("f:connectivity");
    gsSurfMesh::Vertex_property<Point>                                  point = mesh.vertex_property<Point>("v:point",Point(0,0,0));

    // read properties from file
    size_t result;
    result = fread((char*)vconn.data(), sizeof(gsSurfMesh::Vertex_connectivity),   nv, in);
    GISMO_ENSURE(result==sizeof(gsSurfMesh::Vertex_connectivity),"Vertex connectivity reading error");
    result = fread((char*)hconn.data(), sizeof(gsSurfMesh::Halfedge_connectivity), nh, in);
    GISMO_ENSURE(result==sizeof(gsSurfMesh::Halfedge_connectivity),"Vertex connectivity reading error");
    result = fread((char*)fconn.data(), sizeof(gsSurfMesh::Face_connectivity),     nf, in);
    GISMO_ENSURE(result==sizeof(gsSurfMesh::Face_connectivity),"Vertex connectivity reading error");
    result = fread((char*)point.data(), sizeof(Point),                               nv, in);
    GISMO_ENSURE(result==sizeof(Point),"Vertex connectivity reading error");

    fclose(in);
    return true;
}


//-----------------------------------------------------------------------------


bool write_poly(const gsSurfMesh& mesh, const std::string& filename)
{
    // check for colors
    auto color = mesh.get_vertex_property<Color>("v:color");
    bool has_colors = color;


    // open file (in binary mode)
    FILE* out = fopen(filename.c_str(), "wb");
    if (!out) return false;


    // how many elements?
    unsigned int nv, ne, nh, nf;
    nv = mesh.n_vertices();
    ne = mesh.n_edges();
    nh = mesh.n_halfedges();
    nf = mesh.n_faces();

    write(out, nv);
    write(out, ne);
    write(out, nf);
    write(out, has_colors);
    nh = 2*ne;


    // get properties
    gsSurfMesh::Vertex_property<gsSurfMesh::Vertex_connectivity>      vconn = mesh.get_vertex_property<gsSurfMesh::Vertex_connectivity>("v:connectivity");
    gsSurfMesh::Halfedge_property<gsSurfMesh::Halfedge_connectivity>  hconn = mesh.get_halfedge_property<gsSurfMesh::Halfedge_connectivity>("h:connectivity");
    gsSurfMesh::Face_property<gsSurfMesh::Face_connectivity>          fconn = mesh.get_face_property<gsSurfMesh::Face_connectivity>("f:connectivity");
    gsSurfMesh::Vertex_property<Point>                                  point = mesh.get_vertex_property<Point>("v:point");


    // write properties to file
    fwrite((char*)vconn.data(), sizeof(gsSurfMesh::Vertex_connectivity),   nv, out);
    fwrite((char*)hconn.data(), sizeof(gsSurfMesh::Halfedge_connectivity), nh, out);
    fwrite((char*)fconn.data(), sizeof(gsSurfMesh::Face_connectivity),     nf, out);
    fwrite((char*)point.data(), sizeof(Point),                               nv, out);

    if (has_colors) fwrite((char*)color.data(), sizeof(Color), nv, out);

    fclose(out);

    return true;
}


//=============================================================================
} // namespace gismo
//=============================================================================
