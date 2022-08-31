#include <gismo.h>
#include <math.h>
#include <gsParasolid/gsWriteParasolid.h>

using namespace gismo;


void sinus2D(real_t u_min, real_t u_max, real_t v_min, real_t v_max, int n)
{
    gsInfo << "Start example 2D Paper 2, Fast Formation of LS Matrices (but without hole!): \n"
              "Sampling n points from sinus with a hole in [u_min,u_max]." << std::endl;

    // Initializing
    std::vector<real_t> u;
    std::vector<real_t> v;
    std::vector<real_t> x;
    std::vector<real_t> y;
    std::vector<real_t> z;

    for(index_t i=-n; i<=n; i++)
    {
        for(index_t j=-n; j<=n; j++)
        {
            real_t par_u = (1.0*i)/(1.0*n);
            real_t par_v = (1.0*j)/(1.0*n);
            if (!(u_min < par_u && par_u < u_max && v_min < par_v && par_v < v_max ))
            {
                u.push_back((1.0*i)/(1.0*n));
                v.push_back((1.0*j)/(1.0*n));
                x.push_back((1.0*i)/(1.0*n));
                y.push_back((1.0*j)/(1.0*n));
                z.push_back(1.0/3.0*math::sin(4*M_PI*((1.0*i)/(1.0*n)+0.1/8.0))*math::sin(4*M_PI*((1.0*j)/(1.0*n)+0.1/8.0)));
            }
        }
    }

    gsMatrix<real_t> uv(2,u.size());
    gsMatrix<real_t> xyz(3,u.size());

    for(size_t i=0; i<u.size();i++)
    {
        uv(0,i) =  u[i];
        uv(1,i) =  v[i];
        xyz(0,i) = x[i];
        xyz(1,i) = y[i];
        xyz(2,i) = z[i];
    }

    gsFileData<> fd;
    fd << uv;
    fd << xyz;
    fd.dump("example_2DSinus");

    gsWriteParaviewPoints(xyz, "xyz");
}

void holecutter(real_t u_min, real_t u_max, real_t v_min, real_t v_max, std::string fn)
{
    gsInfo << "Start Holecutter" << std::endl;

    gsFileData<> fd_in(fn);
    std::cerr << fd_in << std::endl;
    gsMatrix<> uv, xyz;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, xyz);
    //! [Read data]

    // Cutting
    std::vector<real_t> chosen_u;
    std::vector<real_t> chosen_v;
    std::vector<real_t> chosen_x;
    std::vector<real_t> chosen_y;
    std::vector<real_t> chosen_z;
    for(int col=0; col<uv.cols(); col++)
    {
        if ( !(u_min <= uv(0,col) && uv(0,col) <= u_max && v_min <= uv(1,col) && uv(1,col) <= v_max ) )
        {
            chosen_u.push_back( uv(0,col));
            chosen_v.push_back( uv(1,col));
            chosen_x.push_back(xyz(0,col));
            chosen_y.push_back(xyz(1,col));
            chosen_z.push_back(xyz(2,col));
        }
    }

    gsMatrix<>  chosen_uv(2,chosen_u.size());
    gsMatrix<> chosen_xyz(3,chosen_x.size());

    for(size_t i=0; i<chosen_u.size();i++)
    {
        chosen_uv(0,i) =  chosen_u[i];
        chosen_uv(1,i) =  chosen_v[i];
        chosen_xyz(0,i) = chosen_x[i];
        chosen_xyz(1,i) = chosen_y[i];
        chosen_xyz(2,i) = chosen_z[i];
    }

     gsFileData<> fd;
     fd << chosen_uv;
     fd << chosen_xyz ;

     fd.dump("surface_with_hole");

     gsWriteParaviewPoints(chosen_xyz, "xyz_with_hole");
}

void surface_intersecting(index_t n = 100)
{
    gsInfo << "Start Example Intersecting from Paper 1, support--guided method" << std::endl;

    std::string file = "fitting/examplethomas1.xml";

    // Read Surface
    gsFileData<> fd_in(file);
    gsTensorBSpline<2> geo;
    fd_in.getId<gsTensorBSpline<2> >(0,geo);
    gsInfo << geo;

    std::vector<real_t> u;
    std::vector<real_t> v;

    // parameter samples
    for(int i=0; i<=n; i++)
    {
        for(int j=0; j<=n; j++)
        {
            real_t curr_u = (1.0*i)/n;
            real_t curr_v = (1.0*j)/n;
            if (true)//(curr_u - origin)*(curr_u - origin)+(curr_v - origin)*(curr_v - origin)>=radius*radius)
            //if(math::abs(curr_u - origin) >= radius  || math::abs(curr_v - origin) >= radius )
            {
                u.push_back(curr_u);
                v.push_back(curr_v);
            }
        }
    }

    // -----------------------------Saving u,v and x,y,z into xml file:
    gsMatrix<>  uv(2,u.size());
    gsMatrix<> xyz(3,u.size());

    for(size_t i=0; i<u.size();i++)
    {
        uv(0,i) =  u[i];
        uv(1,i) =  v[i];
    }

    geo.eval_into(uv, xyz);

     gsFileData<> fd;
     fd << uv;
     fd << xyz;
     fd.dump("example_intersecting");

     gsWriteParaviewPoints(xyz, "xyz_intersecting");

     holecutter(0.06,0.94,0.07,0.94,"example_intersecting");
}

void split_in_3(index_t dir,real_t s1, real_t s2, const gsTensorBSpline<2> &input, gsTensorBSpline<2> &out1, gsTensorBSpline<2> &out2, gsTensorBSpline<2> &out3)
{
    // Prüfe s1 < s2 und in domain
    gsTensorBSpline<2> temp;

    input.splitAt( dir, s1, out1, temp);
    temp.splitAt( dir, s2, out2, out3);
}

void surface_subdivision_9(std::string fn)
{
    // Inititalization
    real_t u_min = 0.05;
    real_t u_max = 0.95;
    real_t v_min = 0.05;
    real_t v_max = 0.95;

    gsFileData<> fd_in(fn);
    gsTensorBSpline<2> geo;
    fd_in.getId<gsTensorBSpline<2> >(0,geo);
    gsInfo << geo;

    gsTensorBSpline<2> uleft, umiddle, uright;

    split_in_3(0, u_min, u_max, geo, uleft, umiddle, uright);

    gsTensorBSpline<2> vbottom, vmiddle, vtop;

    split_in_3(1, v_min, v_max, uleft, vbottom, vmiddle, vtop);

    gsInfo << "Writing NX File" << std::endl;
    extensions::gsWritePK_SHEET(vbottom, "uleftvbottom");
    extensions::gsWritePK_SHEET(vmiddle, "uleftvmiddle");
    extensions::gsWritePK_SHEET(vtop, "uleftvtop");

    split_in_3(1, v_min, v_max, umiddle, vbottom, vmiddle, vtop);

    gsInfo << "Writing NX File" << std::endl;
    extensions::gsWritePK_SHEET(vbottom, "umiddlevbottom");
    extensions::gsWritePK_SHEET(vmiddle, "umiddlevmiddle");
    extensions::gsWritePK_SHEET(vtop, "umiddlevtop");

    split_in_3(1, v_min, v_max, uright, vbottom, vmiddle, vtop);

    gsInfo << "Writing NX File" << std::endl;
    extensions::gsWritePK_SHEET(vbottom, "urightvbottom");
    extensions::gsWritePK_SHEET(vmiddle, "urightvmiddle");
    extensions::gsWritePK_SHEET(vtop, "urightvtop");
}

int main(int argc, char *argv[])
{
    // Parameter, Input
    index_t example = 1;
    std::string fn = "/ya/ya135/ya13515/x/build-gismo-Desktop-Debug/bin/result.xml";

    gsCmdLine cmd("Creating samples data from a surface:\n"
                  "Example 1: Intersecting surface with hole from Paper 1. \n"
                  "Example 2: Sinus times Sinus example with hole. \n"
                  "Example 3: Subdivision 9 of fn for Example 1.");
    cmd.addInt("e", "example", "Switch for example", example);
    cmd.addString("d", "data", "Input surface", fn);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (example == 1)
        surface_intersecting();
    else if (example == 2)
        sinus2D(-0.5,0.5,-0.5,0.5,2000);
    else if (example == 3)
        surface_subdivision_9(fn);

    return 0;
}
