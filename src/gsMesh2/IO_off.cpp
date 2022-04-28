
#include <gsMesh2/IO.h>

#include <cstdio>


//== NAMESPACE ================================================================


namespace gismo {


//== IMPLEMENTATION ===========================================================


// helper function
template <typename T> int read(FILE* in, T& t)
{
    int err = 0;
    err = fread(&t, 1, sizeof(t), in);
    return err;
}


// avoid copying in stringstreams !!
struct membuf : std::streambuf
{
    membuf(const char* begin)
    {
        char* p(const_cast<char*>(begin));
        auto sz = sizeof(begin)/(strlen(begin)*sizeof(char));
        this->setg(p, p, p+sz);
    }
};

struct imemstream: virtual membuf, std::istream
{
    imemstream(char const* base)
    : membuf(base), std::istream(static_cast<std::streambuf*>(this)) { }
};

//-----------------------------------------------------------------------------


bool read_off_ascii(gsSurfMesh& mesh,
                    char * node)
{
//    std::fstream fs;
//    fs.getline(in, sizeof buffer );

    char line[200];
    //int                  nc;
    unsigned int         i, j, idx;
    unsigned int         nV, nF, nE;
    Point                p, n, c;
    Vec2f                t;
    gsSurfMesh::Vertex v;

    gsDebugVar( strlen(node) );
    imemstream is(node);
    gsDebugVar( is.eof() );
    
/*
    // properties
    gsSurfMesh::Vertex_property<Normal>              normals;
    gsSurfMesh::Vertex_property<Texture_coordinate>  texcoords;
    gsSurfMesh::Vertex_property<Color>               colors;
    if (has_normals)   normals   = mesh.vertex_property<Normal>("v:normal",Point(0,0,0));
    if (has_texcoords) texcoords = mesh.vertex_property<Texture_coordinate>("v:texcoord",Point(0,0,0));
    if (has_colors)    colors    = mesh.vertex_property<Color>("v:color",Color(0,0,0));
*/

    if ( is.getline (line,200) ) gsDebugVar(std::string(line)); else std::cout<< "error\n";

    // #Vertice, #Faces, #Edges
//    items = fscanf(in, "%d %d %d\n", (int*)&nV, (int*)&nF, (int*)&nE);
    is >> nV >> nF >> nE;
    gsDebugVar(nV);
    gsDebugVar(nF);
    gsDebugVar(nE);
    //(void)items;
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
        v = mesh.add_vertex(p.cast<gsSurfMesh::Scalar>());
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
    std::vector<gsSurfMesh::Vertex> vertices;
    for (i=0; i<nF; ++i)
    {
        // read line
        //lp = fgets(line, 200, in);
        //lp = line;

        // #vertices
        is >> nV;
        gsDebugVar(nV);
        //items = sscanf(lp, "%d%n", (int*)&nV, &nc);
        //assert(items == 1);
        vertices.resize(nV);
        //lp += nc;

        // indices
        for (j=0; j<nV; ++j)
        {
            is >> idx;
            gsDebugVar(idx);
            //items = sscanf(lp, "%d%n", (int*)&idx, &nc);
            //assert(items == 1);
            vertices[j] = gsSurfMesh::Vertex(idx);
            //lp += nc;
        }
        mesh.add_face(vertices);
    }

    return true;
}


bool read_off_ascii(gsSurfMesh& mesh,
                    FILE* in,
                    const bool has_normals,
                    const bool has_texcoords,
                    const bool has_colors)
{
//    std::fstream fs;
//    fs.getline(in, sizeof buffer );

    char                 line[200], *lp;
    int                  nc;
    unsigned int         i, j, items, idx;
    unsigned int         nV, nF, nE;
    Point                p, n, c;
    Vec2f                t;
    gsSurfMesh::Vertex v;

    // properties
    gsSurfMesh::Vertex_property<Normal>              normals;
    gsSurfMesh::Vertex_property<Texture_coordinate>  texcoords;
    gsSurfMesh::Vertex_property<Color>               colors;
    if (has_normals)   normals   = mesh.vertex_property<Normal>("v:normal",Point(0,0,0));
    if (has_texcoords) texcoords = mesh.vertex_property<Texture_coordinate>("v:texcoord",Point(0,0,0));
    if (has_colors)    colors    = mesh.vertex_property<Color>("v:color",Color(0,0,0));


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
        items = sscanf(lp, "%f %f %f%n", (float*)&p[0], (float*)&p[1], (float*)&p[2], &nc);
        assert(items==3);
        v = mesh.add_vertex(p.cast<gsSurfMesh::Scalar>());
        lp += nc;

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
    }



    // read faces: #N v[1] v[2] ... v[n-1]
    std::vector<gsSurfMesh::Vertex> vertices;
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
            vertices[j] = gsSurfMesh::Vertex(idx);
            lp += nc;
        }
        mesh.add_face(vertices);
    }


    return true;
}


//-----------------------------------------------------------------------------


bool read_off_binary(gsSurfMesh& mesh,
                     FILE* in,
                     const bool has_normals,
                     const bool has_texcoords,
                     const bool has_colors)
{
    unsigned int       i, j, idx;
    unsigned int       nV, nF, nE;
    Point              p, n, c;
    Vec2f              t;
    gsSurfMesh::Vertex  v;


    // binary cannot (yet) read colors
    if (has_colors) return false;


    // properties
    gsSurfMesh::Vertex_property<Normal>              normals;
    gsSurfMesh::Vertex_property<Texture_coordinate>  texcoords;
    if (has_normals)   normals   = mesh.vertex_property<Normal>("v:normal",Point(0,0,0));
    if (has_texcoords) texcoords = mesh.vertex_property<Texture_coordinate>("v:texcoord", Texture_coordinate(0,0,0));


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
        v = mesh.add_vertex(p.cast<gsSurfMesh::Scalar>());

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
    std::vector<gsSurfMesh::Vertex> vertices;
    for (i=0; i<nF; ++i)
    {
        read(in, nV);
        vertices.resize(nV);
        for (j=0; j<nV; ++j)
        {
            read(in, idx);
            vertices[j] = gsSurfMesh::Vertex(idx);
        }
        mesh.add_face(vertices);
    }


    return true;
}


//-----------------------------------------------------------------------------


bool read_off(gsSurfMesh& mesh, const std::string& filename)
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
               read_off_binary(mesh, in, has_normals, has_texcoords, has_colors) :
               read_off_ascii(mesh, in, has_normals, has_texcoords, has_colors));


    fclose(in);
    return ok;
}


//-----------------------------------------------------------------------------


bool write_off(const gsSurfMesh& mesh, const std::string& filename)
{
    FILE* out = fopen(filename.c_str(), "w");
    if (!out)
        return false;


    bool  has_normals   = false;
    bool  has_texcoords = false;
    bool  has_colors = false;
    gsSurfMesh::Vertex_property<Normal> normals = mesh.get_vertex_property<Normal>("v:normal");
    gsSurfMesh::Vertex_property<Texture_coordinate>  texcoords = mesh.get_vertex_property<Texture_coordinate>("v:texcoord");
    gsSurfMesh::Vertex_property<Color> colors = mesh.get_vertex_property<Color>("v:color");
    if (normals)   has_normals = true;
    if (texcoords) has_texcoords = true;
    if (colors) has_colors = true;


    // header
    if(has_texcoords)
        fprintf(out, "ST");
    if(has_colors)
        fprintf(out, "C");
    if(has_normals)
        fprintf(out, "N");
    fprintf(out, "OFF\n%d %d 0\n", mesh.n_vertices(), mesh.n_faces());


    // vertices, and optionally normals and texture coordinates
    gsSurfMesh::Vertex_property<Point> points = mesh.get_vertex_property<Point>("v:point");
    for (gsSurfMesh::Vertex_iterator vit=mesh.vertices_begin(); vit!=mesh.vertices_end(); ++vit)
    {
        const Point& p = points[*vit];
        fprintf(out, "%.10f %.10f %.10f", (double)p[0], (double)p[1], (double)p[2]);

        if (has_normals)
        {
            const Normal& n = normals[*vit];
            fprintf(out, " %.10f %.10f %.10f", (double)n[0], (double)n[1], (double)n[2]);
        }

        if (has_colors)
        {
            const Color& c = colors[*vit];
            fprintf(out, " %.10f %.10f %.10f", (double)c[0], (double)c[1], (double)c[2]);
        }

        if (has_texcoords)
        {
            const Texture_coordinate& t = texcoords[*vit];
            fprintf(out, " %.10f %.10f", (double)t[0], (double)t[1]);
        }

        fprintf(out, "\n");
    }


    // faces
    for (gsSurfMesh::Face_iterator fit=mesh.faces_begin(); fit!=mesh.faces_end(); ++fit)
    {
        int nV = mesh.valence(*fit);
        fprintf(out, "%d", nV);
        gsSurfMesh::Vertex_around_face_circulator fvit=mesh.vertices(*fit), fvend=fvit;
        do
        {
            fprintf(out, " %d", (*fvit).idx());
        }
        while (++fvit != fvend);
        fprintf(out, "\n");
    }

    fclose(out);
    return true;
}


//=============================================================================
} // namespace gismo
//=============================================================================
