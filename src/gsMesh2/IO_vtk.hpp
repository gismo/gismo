
#include <gsMesh2/IO.h>

#include <cstdio>


//== NAMESPACE ================================================================


namespace gismo {


//== IMPLEMENTATION ===========================================================


// helper function
static inline int io_vtk_read_file_type(FILE* in, unsigned char *buf, size_t size)
{
    return (int)fread(buf, 1, size, in);
}

namespace {
template <typename T>
static inline int io_vtk_read(FILE* in, T& t)
{
    int err = 0;
    err = fread(&t, 1, sizeof(t), in);
    return err;
}
}


//-----------------------------------------------------------------------------

template <class T>
bool read_vtk(gsSurfMesh<T>& mesh,
              FILE* in,
              const bool has_normals,
              const bool has_texcoords,
              const bool has_colors)
{
    char                 line[200], *lp;
    int                  nc;
    unsigned int         i, j, items, idx;
    unsigned int         nV, nF, nE;
    Point                p, n, c;
    Vec2f                t;
    typename gsSurfMesh<T>::Vertex v;


    // properties
    typename gsSurfMesh<T>::Vertex_property<Normal>              normals;
    typename gsSurfMesh<T>::Vertex_property<Texture_coordinate>  texcoords;
    typename gsSurfMesh<T>::Vertex_property<Color>               colors;
    if (has_normals)   normals   = mesh.template vertex_property<Normal>("v:normal",Point(0,0,0));
    if (has_texcoords) texcoords = mesh.template vertex_property<Texture_coordinate>("v:texcoord",Point(0,0,0));
    if (has_colors)    colors    = mesh.template vertex_property<Color>("v:color",Color(0,0,0));


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
        v = mesh.add_vertex(p.cast<typename gsSurfMesh<T>::Scalar>());
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
    std::vector<typename gsSurfMesh<T>::Vertex> vertices;
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
            vertices[j] = typename gsSurfMesh<T>::Vertex(idx);
            lp += nc;
        }
        mesh.add_face(vertices);
    }


    return true;
}

//-----------------------------------------------------------------------------

template <class T>
bool write_vtk(const gsSurfMesh<T>& mesh, const std::string& filename)
{
    FILE* out = fopen(filename.c_str(), "w");
    if (!out)
        return false;


    bool  has_normals   = false;
    bool  has_texcoords = false;
    bool  has_colors = false;
    typename gsSurfMesh<T>::Vertex_property<Normal> normals = mesh.template get_vertex_property<Normal>("v:normal");
    typename gsSurfMesh<T>::Vertex_property<Texture_coordinate>  texcoords = mesh.template get_vertex_property<Texture_coordinate>("v:texcoord");
    typename gsSurfMesh<T>::Vertex_property<Color> colors = mesh.template get_vertex_property<Color>("v:color");
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
    typename gsSurfMesh<T>::Vertex_property<Point> points = mesh.template get_vertex_property<Point>("v:point");
    for (typename gsSurfMesh<T>::Vertex_iterator vit=mesh.vertices_begin(); vit!=mesh.vertices_end(); ++vit)
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
            const Color& c = colors[*vit];
            fprintf(out, " %.10f %.10f %.10f", cast<real_t,double>(c[0]), cast<real_t,double>(c[1]), cast<real_t,double>(c[2]));
        }

        if (has_texcoords)
        {
            const Texture_coordinate& t = texcoords[*vit];
            fprintf(out, " %.10f %.10f", cast<real_t,double>(t[0]), cast<real_t,double>(t[1]));
        }

        fprintf(out, "\n");
    }


    // faces
    for (typename gsSurfMesh<T>::Face_iterator fit=mesh.faces_begin(); fit!=mesh.faces_end(); ++fit)
    {
        int nV = mesh.valence(*fit);
        fprintf(out, "%d", nV);
        typename gsSurfMesh<T>::Vertex_around_face_circulator fvit=mesh.vertices(*fit), fvend=fvit;
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

// wrapper that takes a filename and forwards to the FILE* based reader
template <class T>
bool read_vtk(gsSurfMesh<T>& mesh, const std::string & filename)
{
    FILE* in = fopen(filename.c_str(), "r");
    if (!in) return false;
    bool ok = read_vtk(mesh, in, false, false, false);
    fclose(in);
    return ok;
}


//=============================================================================
} // namespace gismo
//=============================================================================
