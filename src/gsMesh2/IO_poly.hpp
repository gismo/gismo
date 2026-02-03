
#include <gsMesh2/IO.h>

#include <cstdio>
#include <vector>
#include <string>


//== NAMESPACES ===============================================================


namespace gismo {


//== IMPLEMENTATION ===========================================================


// helper function
static inline size_t io_poly_read_file_binary(FILE* in, void* t, size_t size)
{
    return fread((char*)t, 1, size, in);
}


// helper function
static inline size_t io_poly_write_file_binary(FILE* out, const void* t, size_t size)
{
    return fwrite((char*)t, 1, size, out);
}

// compatibility wrappers for prior 'read'/'write' template names
namespace {
template <typename T>
static inline size_t io_poly_read(FILE* in, T& t)
{
    return fread((char*)&t, 1, sizeof(t), in);
}

template <typename T>
static inline size_t io_poly_write(FILE* out, const T& t)
{
    return fwrite((char*)&t, 1, sizeof(t), out);
}
} // anonymous namespace


//-----------------------------------------------------------------------------


template<class Scalar>
bool read_poly(gsSurfMesh<Scalar>& mesh, const std::string& filename)
{
    // open file (in binary mode)
    FILE* in = fopen(filename.c_str(), "rb");
    if (!in) return false;


    // clear mesh
    mesh.clear();


    // how many elements?
    unsigned int nv, ne, nh, nf;
    io_poly_read(in, nv);
    io_poly_read(in, ne);
    io_poly_read(in, nf);
    nh = 2*ne;


    // resize containers
    mesh.vprops_.resize(nv);
    mesh.hprops_.resize(nh);
    mesh.eprops_.resize(ne);
    mesh.fprops_.resize(nf);


    // get properties
    typename gsSurfMesh<Scalar>::Vertex_property<typename gsSurfMesh<Scalar>::Vertex_connectivity>      vconn = mesh.template vertex_property<typename gsSurfMesh<Scalar>::Vertex_connectivity>("v:connectivity");
    typename gsSurfMesh<Scalar>::Halfedge_property<typename gsSurfMesh<Scalar>::Halfedge_connectivity>  hconn = mesh.template halfedge_property<typename gsSurfMesh<Scalar>::Halfedge_connectivity>("h:connectivity");
    typename gsSurfMesh<Scalar>::Face_property<typename gsSurfMesh<Scalar>::Face_connectivity>          fconn = mesh.template face_property<typename gsSurfMesh<Scalar>::Face_connectivity>("f:connectivity");
    typename gsSurfMesh<Scalar>::Vertex_property<typename gsSurfMesh<Scalar>::Point>                                  point = mesh.template vertex_property<typename gsSurfMesh<Scalar>::Point>("v:point",typename gsSurfMesh<Scalar>::Point(0,0,0));

    // read properties from file
    size_t result;
    result = fread((char*)vconn.data(), sizeof(typename gsSurfMesh<Scalar>::Vertex_connectivity),   nv, in);
    GISMO_ENSURE(result==sizeof(typename gsSurfMesh<Scalar>::Vertex_connectivity),"Vertex connectivity reading error");
    result = fread((char*)hconn.data(), sizeof(typename gsSurfMesh<Scalar>::Halfedge_connectivity), nh, in);
    GISMO_ENSURE(result==sizeof(typename gsSurfMesh<Scalar>::Halfedge_connectivity),"Vertex connectivity reading error");
    result = fread((char*)fconn.data(), sizeof(typename gsSurfMesh<Scalar>::Face_connectivity),     nf, in);
    GISMO_ENSURE(result==sizeof(typename gsSurfMesh<Scalar>::Face_connectivity),"Vertex connectivity reading error");
    result = fread((char*)point.data(), sizeof(typename gsSurfMesh<Scalar>::Point),                               nv, in);
    GISMO_ENSURE(result==sizeof(typename gsSurfMesh<Scalar>::Point),"Vertex connectivity reading error");

    fclose(in);
    return true;
}


//-----------------------------------------------------------------------------

template<class Scalar>
bool write_poly(const gsSurfMesh<Scalar>& mesh, const std::string& filename)
{
    // check for colors
    auto color = mesh.template get_vertex_property<typename gsSurfMesh<Scalar>::Point>("v:color");
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

    io_poly_write(out, nv);
    io_poly_write(out, ne);
    io_poly_write(out, nf);
    io_poly_write(out, has_colors);
    nh = 2*ne;


    // get properties
    typename gsSurfMesh<Scalar>::Vertex_property<typename gsSurfMesh<Scalar>::Vertex_connectivity>      vconn = mesh.template get_vertex_property<typename gsSurfMesh<Scalar>::Vertex_connectivity>("v:connectivity");
    typename gsSurfMesh<Scalar>::Halfedge_property<typename gsSurfMesh<Scalar>::Halfedge_connectivity>  hconn = mesh.template get_halfedge_property<typename gsSurfMesh<Scalar>::Halfedge_connectivity>("h:connectivity");
    typename gsSurfMesh<Scalar>::Face_property<typename gsSurfMesh<Scalar>::Face_connectivity>          fconn = mesh.template get_face_property<typename gsSurfMesh<Scalar>::Face_connectivity>("f:connectivity");
    typename gsSurfMesh<Scalar>::Vertex_property<typename gsSurfMesh<Scalar>::Point>                                  point = mesh.template get_vertex_property<typename gsSurfMesh<Scalar>::Point>("v:point");


    // write properties to file
    fwrite((char*)vconn.data(), sizeof(typename gsSurfMesh<Scalar>::Vertex_connectivity),   nv, out);
    fwrite((char*)hconn.data(), sizeof(typename gsSurfMesh<Scalar>::Halfedge_connectivity), nh, out);
    fwrite((char*)fconn.data(), sizeof(typename gsSurfMesh<Scalar>::Face_connectivity),     nf, out);
    fwrite((char*)point.data(), sizeof(typename gsSurfMesh<Scalar>::Point),                               nv, out);

    if (has_colors) fwrite((char*)color.data(), sizeof(typename gsSurfMesh<Scalar>::Point), nv, out);

    fclose(out);

    return true;
}


//=============================================================================
} // namespace gismo
//=============================================================================
