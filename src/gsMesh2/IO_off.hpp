
#include <gsMesh2/IO.h>

#include <cstdio>


//== NAMESPACE ================================================================


namespace gismo {


//== IMPLEMENTATION ===========================================================

namespace {
// helper function
template <typename T> int read(FILE* in, T& t)
{
    int err = 0;
    err = fread(&t, 1, sizeof(t), in);
    return err;
}
} // anonymous namespace


// avoid copying in stringstreams !!
struct membuf : std::streambuf
{
    membuf(const char* begin)
    {
        char* p(const_cast<char*>(begin));
        // note: the buffer length is the string length; computing it as
        // sizeof(begin)/... would measure the *pointer*, giving 0 and thus an
        // immediately-empty stream.
        this->setg(p, p, p + std::strlen(begin));
    }
};

struct imemstream: virtual membuf, std::istream
{
    imemstream(char const* base)
    : membuf(base), std::istream(static_cast<std::streambuf*>(this)) { }
};


template <class Scalar>
bool read_off_ascii(gsSurfMesh<Scalar>& mesh,
                    char * node)
{
    typedef typename gsSurfMesh<Scalar>::Point Point;
//    std::fstream fs;
//    fs.getline(in, sizeof buffer );

    char line[200];
    //int                  nc;
    unsigned int         i, j, idx;
    unsigned int         nV, nF, nE;
    Point                 p, n, c;
    gsVector<real_t,2>   t;
    typename gsSurfMesh<Scalar>::Vertex v;

    imemstream is(node);

/*
    // properties
    typename gsSurfMesh<Scalar>::Vertex_property<Normal>              normals;
    typename gsSurfMesh<Scalar>::Vertex_property<Point>  texcoords;
    typename gsSurfMesh<Scalar>::Vertex_property<Point>               colors;
    if (has_normals)   normals   = mesh.template vertex_property<Normal>("v:normal",Point(0,0,0));
    if (has_texcoords) texcoords = mesh.template vertex_property<Point>("v:texcoord",Point(0,0,0));
    if (has_colors)    colors    = mesh.template vertex_property<Point>("v:color",Point(0,0,0));
*/

    // header: [ST][C][N][4][n]OFF
    if ( !is.getline(line,200) ) return false;
    if ( NULL == strstr(line, "OFF") ) return false;

    // #Vertice, #Faces, #Edges
    is >> nV >> nF >> nE;
    if ( !is ) return false;

    mesh.clear();
    mesh.reserve(nV, std::max(3*nV, nE), nF);

    // read vertices: pos [normal] [color] [texcoord]
    for (i=0; i<nV && !is.eof(); ++i)
    {
        // read line
        //lp = is.getline(line, 200);
        //lp = line;

        // position
        //items = sscanf(lp, "%f %f %f%n", (float*)&p[0], (float*)&p[1], (float*)&p[2], &nc);
        //assert(items==3);

        is >> p[0] >> p[1] >> p[2];
        v = mesh.add_vertex(p.template cast<Scalar>());
        //lp += nc;
/*
        // normal
        if (has_normals)
        {
            if (sscanf(lp, "%f %f %f%n", (float*)&n[0], (float*)&n[1], (float*)&n[2], &nc) == 3)
            {
                normals[v] = n;
            }
            lp += nc;
        }

        // color
        if (has_colors)
        {
            if (sscanf(lp, "%f %f %f%n", (float*)&c[0], (float*)&c[1], (float*)&c[2], &nc) == 3)
            {
                if (c[0]>1.0f || c[1]>1.0f || c[2]>1.0f) c *= (1.0/255.0);
                colors[v] = c;
            }
            lp += nc;
        }

        // tex coord
        if (has_texcoords)
        {
            items = sscanf(lp, "%f %f%n", (float*)&t[0], (float*)&t[1], &nc);
            assert(items == 2);
            texcoords[v][0] = t[0];
            texcoords[v][1] = t[1];
            lp += nc;
        }
*/
    }

    // read faces: #N v[1] v[2] ... v[n-1]
    std::vector<typename gsSurfMesh<Scalar>::Vertex> vertices;
    for (i=0; i<nF; ++i)
    {
        // read line
        //lp = fgets(line, 200, in);
        //lp = line;

        // #vertices
        is >> nV;
        vertices.resize(nV);

        // indices
        for (j=0; j<nV; ++j)
        {
            is >> idx;
            vertices[j] = typename gsSurfMesh<Scalar>::Vertex(idx);
        }
        mesh.add_face(vertices);
    }

    return true;
}


template <class Scalar>
bool read_off_ascii(gsSurfMesh<Scalar>& mesh,
                    FILE* in,
                    const bool has_normals,
                    const bool has_texcoords,
                    const bool has_colors)
{
    typedef typename gsSurfMesh<Scalar>::Point Point;
    typedef typename gsSurfMesh<Scalar>::Normal Normal;
//    std::fstream fs;
//    fs.getline(in, sizeof buffer );

    char                 line[200], *lp;
    int                  nc;
    unsigned int         i, j, items, idx;
    unsigned int         nV, nF, nE;
    Point                p, n, c;
    gsVector<real_t,2>  t;
    typename gsSurfMesh<Scalar>::Vertex v;

    // properties
    typename gsSurfMesh<Scalar>::template Vertex_property<Normal>              normals;
    typename gsSurfMesh<Scalar>::template Vertex_property<Point>  texcoords;
    typename gsSurfMesh<Scalar>::template Vertex_property<Point>               colors;
    if (has_normals)   normals   = mesh.template vertex_property<Normal>("v:normal",Point(0,0,0));
    if (has_texcoords) texcoords = mesh.template vertex_property<Point>("v:texcoord",Point(0,0,0));
    if (has_colors)    colors    = mesh.template vertex_property<Point>("v:color",Point(0,0,0));


    // #Vertice, #Faces, #Edges
    items = fscanf(in, "%d %d %d\n", (int*)&nV, (int*)&nF, (int*)&nE);
    (void)items;
    mesh.clear();
    mesh.reserve(nV, std::max(3*nV, nE), nF);


    // read vertices: pos [normal] [color] [texcoord]
    for (i=0; i<nV && !feof(in); ++i)
    {
        // read line
        lp = fgets(line, 200, in);
        lp = line;

        // position
        // note: scan into double temporaries. Writing "%f" through a
        // (float*) cast of a Scalar (=double) component overwrites only half
        // of it and leaves garbage -- it must not be re-introduced here.
        double px(0), py(0), pz(0);
        items = sscanf(lp, "%lf %lf %lf%n", &px, &py, &pz, &nc);
        assert(items==3);
        p[0] = static_cast<Scalar>(px);
        p[1] = static_cast<Scalar>(py);
        p[2] = static_cast<Scalar>(pz);
        v = mesh.add_vertex(p.template cast<Scalar>());
        lp += nc;

        // normal
        if (has_normals)
        {
            double nx(0), ny(0), nz(0);
            if (sscanf(lp, "%lf %lf %lf%n", &nx, &ny, &nz, &nc) == 3)
            {
                n[0] = static_cast<Scalar>(nx);
                n[1] = static_cast<Scalar>(ny);
                n[2] = static_cast<Scalar>(nz);
                normals[v] = n;
            }
            lp += nc;
        }

        // color
        if (has_colors)
        {
            double cx(0), cy(0), cz(0);
            if (sscanf(lp, "%lf %lf %lf%n", &cx, &cy, &cz, &nc) == 3)
            {
                c[0] = static_cast<Scalar>(cx);
                c[1] = static_cast<Scalar>(cy);
                c[2] = static_cast<Scalar>(cz);
                if (c[0]>1.0f || c[1]>1.0f || c[2]>1.0f) c *= (1.0/255.0);
                colors[v] = c;
            }
            lp += nc;
        }

        // tex coord
        if (has_texcoords)
        {
            double tx(0), ty(0);
            items = sscanf(lp, "%lf %lf%n", &tx, &ty, &nc);
            assert(items == 2);
            t[0] = static_cast<real_t>(tx);
            t[1] = static_cast<real_t>(ty);
            texcoords[v][0] = t[0];
            texcoords[v][1] = t[1];
            lp += nc;
        }
    }



    // read faces: #N v[1] v[2] ... v[n-1]
    std::vector<typename gsSurfMesh<Scalar>::Vertex> vertices;
    for (i=0; i<nF; ++i)
    {
        // read line
        lp = fgets(line, 200, in);
        lp = line;

        // #vertices
        items = sscanf(lp, "%d%n", (int*)&nV, &nc);
        assert(items == 1);
        vertices.resize(nV);
        lp += nc;

        // indices
        for (j=0; j<nV; ++j)
        {
            items = sscanf(lp, "%d%n", (int*)&idx, &nc);
            assert(items == 1);
            vertices[j] = typename gsSurfMesh<Scalar>::Vertex(idx);
            lp += nc;
        }
        mesh.add_face(vertices);
    }


    return true;
}

template <class Scalar>
bool read_off_binary(gsSurfMesh<Scalar>& mesh,
                     FILE* in,
                     const bool has_normals,
                     const bool has_texcoords,
                     const bool has_colors)
{
    typedef typename gsSurfMesh<Scalar>::Point  Point;
    typedef typename gsSurfMesh<Scalar>::Normal Normal;
    unsigned int       i, j, idx;
    unsigned int       nV, nF, nE;
    Point               p, n, c;
    gsVector<real_t,2>  t;
    typename gsSurfMesh<Scalar>::Vertex  v;


    // binary cannot (yet) read colors
    if (has_colors) return false;


    // properties
    typename gsSurfMesh<Scalar>::template Vertex_property<Normal>              normals;
    typename gsSurfMesh<Scalar>::template Vertex_property<Point>  texcoords;
    if (has_normals)   normals   = mesh.template vertex_property<Normal>("v:normal",Point(0,0,0));
    if (has_texcoords) texcoords = mesh.template vertex_property<Point>("v:texcoord", Point(0,0,0));


    // #Vertice, #Faces, #Edges
    read(in, nV);
    read(in, nF);
    read(in, nE);
    mesh.clear();
    mesh.reserve(nV, std::max(3*nV, nE), nF);


    // read vertices: pos [normal] [color] [texcoord]
    for (i=0; i<nV && !feof(in); ++i)
    {
        // position
        read(in, p);
        v = mesh.add_vertex(p.template cast<Scalar>());

        // normal
        if (has_normals)
        {
            read(in, n);
            normals[v] = n;
        }

        // tex coord
        if (has_texcoords)
        {
            read(in, t);
            texcoords[v][0] = t[0];
            texcoords[v][1] = t[1];
        }
    }


    // read faces: #N v[1] v[2] ... v[n-1]
    std::vector<typename gsSurfMesh<Scalar>::Vertex> vertices;
    for (i=0; i<nF; ++i)
    {
        read(in, nV);
        vertices.resize(nV);
        for (j=0; j<nV; ++j)
        {
            read(in, idx);
            vertices[j] = typename gsSurfMesh<Scalar>::Vertex(idx);
        }
        mesh.add_face(vertices);
    }


    return true;
}

template <class Scalar>
bool read_off(gsSurfMesh<Scalar>& mesh, const std::string& filename)
{
    char  line[200];
    bool  has_texcoords = false;
    bool  has_normals   = false;
    bool  has_colors    = false;
    bool  has_hcoords   = false;
    bool  has_dim       = false;
    bool  is_binary     = false;


    // open file (in ASCII mode)
    FILE* in = fopen(filename.c_str(), "r");
    if (!in) return false;


    // read header: [ST][C][N][4][n]OFF BINARY
    char *c = fgets(line, 200, in);
    assert(c != NULL);
    c = line;
    if (c[0] == 'S' && c[1] == 'T') { has_texcoords = true; c += 2; }
    if (c[0] == 'C') { has_colors  = true; ++c; }
    if (c[0] == 'N') { has_normals = true; ++c; }
    if (c[0] == '4') { has_hcoords = true; ++c; }
    if (c[0] == 'n') { has_dim     = true; ++c; }
    if (strncmp(c, "OFF", 3) != 0) { fclose(in); return false; } // no OFF
    if (strncmp(c+4, "BINARY", 6) == 0) is_binary = true;


    // homogeneous coords, and vertex dimension != 3 are not supported
    if (has_hcoords || has_dim)
    {
        fclose(in);
        return false;
    }


    // if binary: reopen file in binary mode
    if (is_binary)
    {
        fclose(in);
        in = fopen(filename.c_str(), "rb");
        c = fgets(line, 200, in);
        assert(c != NULL);
    }


    // read as ASCII or binary
    bool ok = (is_binary ?
               read_off_binary<Scalar>(mesh, in, has_normals, has_texcoords, has_colors) :
               read_off_ascii<Scalar>(mesh, in, has_normals, has_texcoords, has_colors));


    fclose(in);
    return ok;
}

template <class Scalar>
bool write_off(const gsSurfMesh<Scalar>& mesh, const std::string& filename)
{
    typedef typename gsSurfMesh<Scalar>::Point Point;
    typedef typename gsSurfMesh<Scalar>::Normal Normal;
    FILE* out = fopen(filename.c_str(), "w");
    if (!out)
        return false;


    bool  has_normals   = false;
    bool  has_texcoords = false;
    bool  has_colors = false;
    bool  has_sharp = false;
    typename gsSurfMesh<Scalar>::template Vertex_property<Normal> normals = mesh.template get_vertex_property<Normal>("v:normal");
    typename gsSurfMesh<Scalar>::template Vertex_property<Point>  texcoords = mesh.template get_vertex_property<Point>("v:texcoord");
    typename gsSurfMesh<Scalar>::template Vertex_property<Point> colors = mesh.template get_vertex_property<Point>("v:color");
    typename gsSurfMesh<Scalar>::template Halfedge_property<bool> sharp = mesh.template get_halfedge_property<bool>("h:sharp");
    if (normals)   has_normals = true;
    if (texcoords) has_texcoords = true;
    if (colors) has_colors = true;
    if (sharp) has_sharp = true;


    // header
    if(has_texcoords)
        fprintf(out, "ST");
    if(has_colors)
        fprintf(out, "C");
    if(has_normals)
        fprintf(out, "N");
    fprintf(out, "OFF\n%d %d 0\n", mesh.n_vertices(), mesh.n_faces());


    // vertices, and optionally normals and texture coordinates
    typename gsSurfMesh<Scalar>::template Vertex_property<Point> points = mesh.template get_vertex_property<Point>("v:point");
    for (typename gsSurfMesh<Scalar>::Vertex_iterator vit=mesh.vertices_begin(); vit!=mesh.vertices_end(); ++vit)
    {
        const Point& p = points[*vit];
        fprintf(out, "%.10f %.10f %.10f", cast<real_t,double>(p[0]), cast<real_t,double>(p[1]), cast<real_t,double>(p[2]));

        if (has_normals)
        {
            const Normal& n = normals[*vit];
            fprintf(out, " %.10f %.10f %.10f", cast<real_t,double>(n[0]), cast<real_t,double>(n[1]), cast<real_t,double>(n[2]));
        }

        if (has_colors)
        {
            const Point& c = colors[*vit];
            fprintf(out, " %.10f %.10f %.10f", cast<real_t,double>(c[0]), cast<real_t,double>(c[1]), cast<real_t,double>(c[2]));
        }

        if (has_texcoords)
        {
            const Point& t = texcoords[*vit];
            fprintf(out, " %.10f %.10f", cast<real_t,double>(t[0]), cast<real_t,double>(t[1]));
        }

        fprintf(out, "\n");
    }


    // faces
    for (typename gsSurfMesh<Scalar>::Face_iterator fit=mesh.faces_begin(); fit!=mesh.faces_end(); ++fit)
    {
        int nV = mesh.valence(*fit);
        fprintf(out, "%d", nV);
        typename gsSurfMesh<Scalar>::Vertex_around_face_circulator fvit=mesh.vertices(*fit), fvend=fvit;
        do
        {
            fprintf(out, " %d", (*fvit).idx());
        }
        while (++fvit != fvend);
        fprintf(out, "\n");
    }

    // sharp edges 
    if (has_sharp){
        for (typename gsSurfMesh<Scalar>::Edge_iterator it = mesh.edges_begin();
            it != mesh.edges_end();
            ++it)
        {
            auto e = *it;
            auto v1 = mesh.vertex(e, 0);
            auto v2 = mesh.vertex(e, 1);

            if (sharp[mesh.halfedge(e, 0)]){
                fprintf(out, "%d %d\n", v1.idx(), v2.idx());
            }
        }
    }
    
    fclose(out);
    return true;
}


//=============================================================================
} // namespace gismo
//=============================================================================
