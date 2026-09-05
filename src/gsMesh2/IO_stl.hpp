
//== INCLUDES =================================================================

#include <gsMesh2/IO.h>

#include <cstdio>
#include <cfloat>
#include <map>
#include <fstream>


//== NAMESPACES ===============================================================


namespace gismo {


//== IMPLEMENTATION ===========================================================


namespace {
// helper function
template <typename T> int read_stl_data(FILE* in, T& t)
{
    size_t n_items(0);
    (void)n_items;
    n_items = fread((char*)&t, 1, sizeof(t), in);
    GISMO_ASSERT(n_items > 0, "Failed to read from file.");
    return n_items;
}
} // anonymous namespace


//-----------------------------------------------------------------------------


// helper class for STL reader
template <class PointT>
class CmpVec
{
public:

    CmpVec(float _eps=FLT_MIN) : eps_(_eps) {}

    bool operator()(const PointT& v0, const PointT& v1) const
    {
        if (math::abs(v0[0] - v1[0]) <= eps_)
        {
            if (math::abs(v0[1] - v1[1]) <= eps_)
            {
                return (v0[2] < v1[2] - eps_);
            }
            else return (v0[1] < v1[1] - eps_);
        }
        return (v0[0] < v1[0] - eps_);
    }

private:
    float eps_;
};


//-----------------------------------------------------------------------------


template <class Scalar>
bool read_stl(gsSurfMesh<Scalar>& mesh, const std::string& filename)
{
    typedef typename gsSurfMesh<Scalar>::Point Point;
    
    char                             line[100], *c;
    unsigned int                     i, nT;
    typename gsSurfMesh<Scalar>::Point                 p;
    typename gsSurfMesh<Scalar>::Vertex               v;
    std::vector<typename gsSurfMesh<Scalar>::Vertex>  vertices(3);
    size_t n_items(0);
    (void)n_items;
    
    CmpVec<Point> comp(FLT_MIN);
    std::map<Point, typename gsSurfMesh<Scalar>::Vertex, CmpVec<Point> >            vMap(comp);
    typename std::map<Point, typename gsSurfMesh<Scalar>::Vertex, CmpVec<Point> >::iterator  vMapIt;


    // clear mesh
    mesh.clear();


    // open file (in ASCII mode)
    FILE* in = fopen(filename.c_str(), "r");
    if (!in) return false;


    // ASCII or binary STL?
    c = fgets(line, 6, in);
    GISMO_ASSERT(c != NULL, "Failed to read file header.");
    const bool binary = ((strncmp(line, "SOLID", 5) != 0) &&
                         (strncmp(line, "solid", 5) != 0));


    // parse binary STL
    if (binary)
    {
        // re-open file in binary mode
        fclose(in);
        in = fopen(filename.c_str(), "rb");
        if (!in) return false;

        // skip dummy header
        n_items = fread(line, 1, 80, in);
        GISMO_ASSERT(n_items > 0, "Failed to read from file.");

        // read number of triangles
        read_stl_data(in, nT);

        // read triangles
        while (nT)
        {
            // skip triangle normal
            n_items = fread(line, 1, 12, in);
            GISMO_ASSERT(n_items > 0, "Failed to read from file.");
            // triangle's vertices
            for (i=0; i<3; ++i)
            {
                read_stl_data(in, p);

                // has vector been referenced before?
                if ((vMapIt=vMap.find(p)) == vMap.end())
                {
                    // No : add vertex and remember idx/vector mapping
                    v = mesh.add_vertex(p.template cast<Scalar>());
                    vertices[i] = v;
                    vMap[p] = v;
                }
                else
                {
                    // Yes : get index from map
                    vertices[i] = vMapIt->second;
                }
            }

            // Add face only if it is not degenerated
            if ((vertices[0] != vertices[1]) &&
                (vertices[0] != vertices[2]) &&
                (vertices[1] != vertices[2]))
                mesh.add_face(vertices);

            n_items = fread(line, 1, 2, in);
            GISMO_ASSERT(n_items > 0, "Failed to read from file.");
            --nT;
        }
    }


    // parse ASCII STL
    else
    {
        // parse line by line
        while (in && !feof(in) && fgets(line, 100, in))
        {
            // skip white-space
            for (c=line; isspace(*c) && *c!='\0'; ++c) {};

            // face begins
            if ((strncmp(c, "outer", 5) == 0) ||
                (strncmp(c, "OUTER", 5) == 0))
            {
                // read three vertices
                for (i=0; i<3; ++i)
                {
                    // read line
                    c = fgets(line, 100, in);
                    GISMO_ASSERT(c != NULL, "Failed to read file header.");

                    // skip white-space
                    for (c=line; isspace(*c) && *c!='\0'; ++c) {};

                    // read x, y, z
                    // note: scan into double temporaries. Writing "%f" through
                    // a (float*) cast of a Scalar (=double) component fills
                    // only half of it and leaves the rest garbage.
                    double px(0), py(0), pz(0);
                    sscanf(c+6, "%lf %lf %lf", &px, &py, &pz);
                    p[0] = static_cast<Scalar>(px);
                    p[1] = static_cast<Scalar>(py);
                    p[2] = static_cast<Scalar>(pz);

                    // has vector been referenced before?
                    if ((vMapIt=vMap.find(p)) == vMap.end())
                    {
                        // No : add vertex and remember idx/vector mapping
                        v = mesh.add_vertex(p.template cast<Scalar>());
                        vertices[i] = v;
                        vMap[p] = v;
                    }
                    else
                    {
                        // Yes : get index from map
                        vertices[i] = vMapIt->second;
                    }
                }

                // Add face only if it is not degenerated
                if ((vertices[0] != vertices[1]) &&
                    (vertices[0] != vertices[2]) &&
                    (vertices[1] != vertices[2]))
                    mesh.add_face(vertices);
            }
        }
    }


    fclose(in);
    return true;
}


//-----------------------------------------------------------------------------


template <class Scalar>
bool write_stl(const gsSurfMesh<Scalar>& mesh, const std::string& filename)
{
    typedef typename gsSurfMesh<Scalar>::Point Point;
    typedef typename gsSurfMesh<Scalar>::Normal Normal;
    
    if (!mesh.is_triangle_mesh())
    {
        std::cerr << "write_stl: not a triangle mesh!" << std::endl;
        return false;
    }

    auto fnormals = mesh.template get_face_property<Normal>("f:normal");
    if (!fnormals)
    {
        std::cerr << "write_stl: no a face normals present!" << std::endl;
        return false;
    }

    std::ofstream ofs(filename.c_str());
    auto points = mesh.template get_vertex_property<Point>("v:point");

    ofs << "solid stl" << std::endl;
    Normal n;
    Point p;

    for (auto f : mesh.faces())
    {
        n = fnormals[f];
        ofs << "  facet normal ";
        ofs << n[0] << " " << n[1] << " " << n[2] << std::endl;
        ofs << "    outer loop" << std::endl;
        for (auto v : mesh.vertices(f))
        {
            p = points[v];
            ofs << "      vertex ";
            ofs << p[0] << " " << p[1] << " " << p[2] << std::endl;
        }
        ofs << "    endloop" << std::endl;
        ofs << "  endfacet" << std::endl;
    }
    ofs << "endsolid" << std::endl;
    ofs.close();
    return true;
}


//=============================================================================
} // namespace gismo
//=============================================================================
