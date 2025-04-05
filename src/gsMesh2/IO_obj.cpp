
#include <gsMesh2/IO.h>

#include <cstdio>


//== NAMESPACES ===============================================================


namespace gismo {


//== IMPLEMENTATION ===========================================================


bool read_obj(gsSurfMesh& mesh, const std::string& filename)
{
    char   s[200];
    float  x, y, z;
    std::vector<gsSurfMesh::Vertex>  vertices;
    std::vector<Texture_coordinate> all_tex_coords;   //individual texture coordinates
    std::vector<int> halfedge_tex_idx; //texture coordinates sorted for halfedges
    gsSurfMesh::Halfedge_property <Texture_coordinate> tex_coords = mesh.halfedge_property<Texture_coordinate>("h:texcoord", Texture_coordinate(0,0,0));
    bool with_tex_coord=false;

    // clear mesh
    mesh.clear();


    // open file (in ASCII mode)
    FILE* in = fopen(filename.c_str(), "r");
    if (!in) return false;


    // clear line once
    memset(&s, 0, 200);


    // parse line by line (currently only supports vertex positions & faces
    while(in && !feof(in) && fgets(s, 200, in))
    {
        // comment
        if (s[0] == '#' || isspace(s[0])) continue;

        // vertex
        else if (strncmp(s, "v ", 2) == 0)
        {
            if (sscanf(s, "v %f %f %f", &x, &y, &z))
            {
                mesh.add_vertex(Point(x,y,z));
            }
        }
        // normal
        else if (strncmp(s, "vn ", 3) == 0)
        {
          if (sscanf(s, "vn %f %f %f", &x, &y, &z))
          {
            // problematic as it can be either a vertex property when interpolated
            // or a halfedge property for hard edges
          }
        }

        // texture coordinate
        else if (strncmp(s, "vt ", 3) == 0)
        {
          if (sscanf(s, "vt %f %f", &x, &y))
          {
            z=1;
            all_tex_coords.push_back(Texture_coordinate(x,y,z));
          }
        }

        // face
        else if (strncmp(s, "f ", 2) == 0)
        {
          int component(0), nV(0);
          bool endOfVertex(false);
          char *p0, *p1(s+1);

          vertices.clear();
          halfedge_tex_idx.clear();

          // skip white-spaces
          while (*p1==' ') ++p1;

          while (p1)
          {
            p0 = p1;

            // overwrite next separator

            // skip '/', '\n', ' ', '\0', '\r' <-- don't forget Windows
            while (*p1!='/' && *p1!='\r' && *p1!='\n' && *p1!=' ' && *p1!='\0') ++p1;

            // detect end of vertex
            if (*p1 != '/')
            {
              endOfVertex = true;
            }

            // replace separator by '\0'
            if (*p1 != '\0')
            {
              *p1 = '\0';
              p1++; // point to next token
            }

            // detect end of line and break
            if (*p1 == '\0' || *p1 == '\n')
            {
              p1 = 0;
            }

            // read next vertex component
            if (*p0 != '\0')
            {
              switch (component)
              {
                case 0: // vertex
                {
                  vertices.push_back( gsSurfMesh::Vertex(atoi(p0) - 1) );
                  break;
                }
                case 1: // texture coord
                {
                  int idx = atoi(p0)-1;
                  halfedge_tex_idx.push_back(idx);
                  with_tex_coord=true;
                  break;
                }
                case 2: // normal
                  break;
              }
            }

            ++component;

            if (endOfVertex)
            {
              component = 0;
              nV++;
              endOfVertex = false;
            }
          }

          gsSurfMesh::Face f=mesh.add_face(vertices);


          // add texture coordinates
          if(with_tex_coord)
          {
              gsSurfMesh::Halfedge_around_face_circulator h_fit = mesh.halfedges(f);
              gsSurfMesh::Halfedge_around_face_circulator h_end = h_fit;
              unsigned v_idx =0;
              do
              {
                  tex_coords[*h_fit]=all_tex_coords.at(halfedge_tex_idx.at(v_idx));
                  ++v_idx;
                  ++h_fit;
              }
              while(h_fit!=h_end);
          }
        }
        // clear line
        memset(&s, 0, 200);
    }

    fclose(in);
    return true;
}


//-----------------------------------------------------------------------------


bool write_obj(const gsSurfMesh& mesh, const std::string& filename)
{
    FILE* out = fopen(filename.c_str(), "w");
    if (!out)
        return false;

    // comment
    fprintf(out, "# OBJ export from gsSurfMesh\n");

    //vertices
    gsSurfMesh::Vertex_property<Point> points = mesh.get_vertex_property<Point>("v:point");
    for (gsSurfMesh::Vertex_iterator vit=mesh.vertices_begin(); vit!=mesh.vertices_end(); ++vit)
    {
        const Point& p = points[*vit];
        fprintf(out, "v %.10f %.10f %.10f\n", cast<real_t,double>(p[0]), cast<real_t,double>(p[1]), cast<real_t,double>(p[2]) );
    }

    //normals
    gsSurfMesh::Vertex_property<Point> normals = mesh.get_vertex_property<Point>("v:normal");
    if(normals)
    {
        for (gsSurfMesh::Vertex_iterator vit=mesh.vertices_begin(); vit!=mesh.vertices_end(); ++vit)
        {
            const Point& p = normals[*vit];
            fprintf(out, "vn %.10f %.10f %.10f\n", cast<real_t,double>(p[0]), cast<real_t,double>(p[1]), cast<real_t,double>(p[2]) );
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
        gsSurfMesh::Halfedge_property<Texture_coordinate> tex_coord = mesh.get_halfedge_property<Texture_coordinate>("h:texcoord");
        for (gsSurfMesh::Halfedge_iterator hit=mesh.halfedges_begin(); hit!=mesh.halfedges_end(); ++hit)
        {
            const Texture_coordinate& pt = tex_coord[*hit];
            fprintf(out, "vt %.10f %.10f %.10f\n", cast<real_t,double>(pt[0]), cast<real_t,double>(pt[1]), cast<real_t,double>(pt[2]) );
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


//=============================================================================
} // namespace gismo
//=============================================================================
