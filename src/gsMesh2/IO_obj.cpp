
#include <gsMesh2/IO.h>

#include <cstdio>
#include <istream>
#include <sstream>
#include <string>
#include <vector>

namespace gismo {

bool read_obj(gsSurfMesh& mesh, std::istream& in)
{
    float x, y, z;

    std::vector<gsSurfMesh::Vertex> vertices;
    std::vector<gsSurfMesh::Point> all_tex_coords;
    std::vector<int> halfedge_tex_idx;

    auto tex_coords =
        mesh.halfedge_property<gsSurfMesh::Point>(
            "h:texcoord",
            gsSurfMesh::Point(0,0,0));

    bool with_tex_coord = false;

    mesh.clear();

    std::string line;

    while (std::getline(in, line))
    {
        if (line.empty())
            continue;

        // comment
        if (line[0] == '#')
            continue;

        std::istringstream ls(line);

        std::string tag;
        ls >> tag;

        // vertex
        if (tag == "v")
        {
            if (ls >> x >> y >> z)
            {
                mesh.add_vertex(gsSurfMesh::Point(x,y,z));
            }
        }

        // normal (currently ignored)
        else if (tag == "vn")
        {
            ls >> x >> y >> z;
        }

        // texture coordinate
        else if (tag == "vt")
        {
            if (ls >> x >> y)
            {
                all_tex_coords.emplace_back(x, y, 1);
            }
        }

        // face
        else if (tag == "f")
        {
            vertices.clear();
            halfedge_tex_idx.clear();
            with_tex_coord = false;

            std::string vertexSpec;

            while (ls >> vertexSpec)
            {
                int component = 0;
                std::stringstream vs(vertexSpec);
                std::string token;

                while (std::getline(vs, token, '/'))
                {
                    if (!token.empty())
                    {
                        switch (component)
                        {
                            case 0: // vertex index
                                vertices.push_back(
                                    gsSurfMesh::Vertex(std::stoi(token) - 1));
                                break;

                            case 1: // texture index
                                halfedge_tex_idx.push_back(
                                    std::stoi(token) - 1);
                                with_tex_coord = true;
                                break;

                            case 2: // normal index
                                break;
                        }
                    }

                    ++component;
                }
            }

            gsSurfMesh::Face f = mesh.add_face(vertices);

            if (with_tex_coord)
            {
                auto h_fit = mesh.halfedges(f);
                auto h_end = h_fit;

                unsigned v_idx = 0;

                do
                {
                    tex_coords[*h_fit] =
                        all_tex_coords.at(
                            halfedge_tex_idx.at(v_idx));

                    ++v_idx;
                    ++h_fit;
                }
                while (h_fit != h_end);
            }
        }
    }

    return true;
}

//-----------------------------------------------------------------------------


bool write_obj(const gsSurfMesh& mesh, const std::string& filename)
{
    typedef gsSurfMesh::Point Point;
    FILE* out = fopen(filename.c_str(), "w");
    if (!out)
        return false;

    // comment
    fprintf(out, "# OBJ export from G+Smo\n");

    //vertices
    gsSurfMesh::Vertex_property<Point> points =
        mesh.get_vertex_property<Point>("v:point");
    for (gsSurfMesh::Vertex_iterator vit=mesh.vertices_begin(); vit!=mesh.vertices_end(); ++vit)
    {
        const Point& p = points[*vit];
        fprintf(out, "v %.10f %.10f %.10f\n", cast<real_t,real_t>(p[0]), cast<real_t,real_t>(p[1]), cast<real_t,real_t>(p[2]) );
    }

    //normals
    gsSurfMesh::Vertex_property<Point> normals = mesh.get_vertex_property<Point>("v:normal");
    if(normals)
    {
        for (gsSurfMesh::Vertex_iterator vit=mesh.vertices_begin(); vit!=mesh.vertices_end(); ++vit)
        {
            const Point& p = normals[*vit];
            fprintf(out, "vn %.10f %.10f %.10f\n", cast<real_t,real_t>(p[0]), cast<real_t,real_t>(p[1]), cast<real_t,real_t>(p[2]) );
        }
    }

    //optionally texture coordinates
    // do we have them?
    std::vector<std::string> h_props= mesh.halfedge_properties();
    bool with_tex_coord = false;
    std::vector<std::string>::iterator h_prop_end = h_props.end();
    std::vector<std::string>::iterator h_prop_start= h_props.begin();
    while(h_prop_start!=h_prop_end)
    {
        if(0==(*h_prop_start).compare("h:texcoord"))
        {
            with_tex_coord=true;
        }
        ++h_prop_start;
    }

    //if so then add
    if(with_tex_coord)
    {
        gsSurfMesh::Halfedge_property<Point> tex_coord = mesh.get_halfedge_property<Point>("h:texcoord");
        for (gsSurfMesh::Halfedge_iterator hit=mesh.halfedges_begin(); hit!=mesh.halfedges_end(); ++hit)
        {
            const Point& pt = tex_coord[*hit];
            fprintf(out, "vt %.10f %.10f %.10f\n", cast<real_t,real_t>(pt[0]), cast<real_t,real_t>(pt[1]), cast<real_t,real_t>(pt[2]) );
        }
    }

    //faces
    for (gsSurfMesh::Face_iterator fit=mesh.faces_begin(); fit!=mesh.faces_end(); ++fit)
    {
        fprintf(out, "f");
        gsSurfMesh::Vertex_around_face_circulator fvit=mesh.vertices(*fit), fvend=fvit;
        gsSurfMesh::Halfedge_around_face_circulator fhit=mesh.halfedges(*fit);
        do
        {
            if(with_tex_coord)
            {
                // write vertex index, tex_coord index and normal index
                fprintf(out, " %d/%d/%d", (*fvit).idx()+1, (*fhit).idx()+1, (*fvit).idx()+1);
                ++fhit;
            }
            else
            {
                // write vertex index and normal index
                fprintf(out, " %d//%d", (*fvit).idx()+1, (*fvit).idx()+1);
            }
        }
        while (++fvit != fvend);
        fprintf(out, "\n");
    }

    fclose(out);
    return true;
}

} // namespace gismo
