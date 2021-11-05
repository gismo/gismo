// Fast Fitting idea Bert example surfaces
#include <gismo.h>
#include <math.h>
#include <gsIO/gsFileData.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Initializing
    int example    = 1;
    int n          = 30;
    bool random    = false;

    // String to the file in case the parameters should be read from the file
    std::string fn = "";
    // /ya/ya135/ya13515/x/gismo/filedata/fitting/deepdrawingC.xml
    // /ya/ya135/ya13515/x/Testing/Examples_MTU/Blade/blade.xml

    gsCmdLine cmd("Example for parametrized sample data. Output file is an XML "
            "file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. "
            "Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a "
            "3 x N matrix. Every column represents a point (x,y,z) in space.");
    cmd.addInt("e", "example", "number of example", example);
    cmd.addInt("n", "number", "number of samples in one direction", n);
    cmd.addString("d", "data", "Input sample data", fn);
    cmd.addSwitch("random", "uses random parameter points", random);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    std::vector<real_t> u;
    std::vector<real_t> v;
    std::vector<real_t> x;
    std::vector<real_t> y;
    std::vector<real_t> z;

    gsMatrix<> uv;

    if (fn != "")
    {
        gsFileData<> fd_in(fn);
        fd_in.getId<gsMatrix<> >(0, uv );
    }
    else
    {
        if (random)
        {
            for (index_t i=0; i<n*n; i++)
            {
                u.push_back((double)rand() / (double)RAND_MAX);
                v.push_back((double)rand() / (double)RAND_MAX);
            }
        }
        else
        {
            for (index_t i=0; i<n; i++)
            {
                for (index_t j=0; j<n; j++)
                {
                    u.push_back((1.0*i)/n);
                    v.push_back((1.0*j)/n);
                }
            }
        }
    }

    if (fn == "")
    {
        uv.setZero(2,u.size());
        for(size_t i=0; i<u.size();i++)
        {
            uv(0,i) =  u[i];
            uv(1,i) =  v[i];
        }
    }

    if (example == 1) // Sinus times Sinus function
    {
        for (index_t i=0; i<uv.cols(); i++)
        {
            x.push_back(uv(0,i));
            y.push_back(uv(1,i));
            z.push_back(1.0/3*sin(4*M_PI*(uv(0,i)+0.1/8.0))*sin(4*M_PI*(uv(1,i)+0.1/8.0)));
        }
    }

    if (example == 2) // 3 peaks function
    {
        for (index_t i=0; i<uv.cols(); i++)
        {
            x.push_back(uv(0,i));
            y.push_back(uv(1,i));
            z.push_back(2.0/(3.0*std::exp(std::sqrt((10.0*uv(0,i)-3)*(10.0*uv(0,i)-3)+(10.0*uv(1,i)-3)*(10.0*uv(1,i)-3))))+2.0/(3.0*std::exp(std::sqrt((10.0*uv(0,i)+3)*(10.0*uv(0,i)+3)+(10.0*uv(1,i)+3)*(10.0*uv(1,i)+3))))+2.0/(3.0*std::exp(std::sqrt((10.0*uv(0,i))*(10.0*uv(0,i))+(10.0*uv(1,i))*(10.0*uv(1,i))))));
        }
    }

    if (example == 3) // Noras function
    {
        for (index_t i=0; i<uv.cols(); i++)
        {
            x.push_back(uv(0,i));
            y.push_back(uv(1,i));
            z.push_back(0.1*uv(0,i)*uv(0,i)*uv(0,i)*uv(0,i)*uv(0,i)*uv(0,i)*uv(0,i)*(2.0-2.0*(1.0+0.4*std::sin(60*uv(1,i)))*std::abs(std::cos(2*M_PI*uv(1,i)))));
        }
    }

    if (example == 4) // Sinus times Sinus function, smaller, less sinuses
    {
        for (index_t i=0; i<uv.cols(); i++)
        {
            x.push_back(uv(0,i));
            y.push_back(uv(1,i));
            z.push_back(1.0/10*sin(4*M_PI*(uv(0,i)/2.0))*sin(4*M_PI*(uv(1,i)/2.0)));
        }
    }

    // putting all values in matrix xyz
    gsMatrix<> xyz(3,uv.cols());
    for(index_t i=0; i<uv.cols();i++)
    {
        xyz(0,i) = x[i];
        xyz(1,i) = y[i];
        xyz(2,i) = z[i];
    }

    // Writing XML file:
     gsFileData<> fd;
     fd << uv;
     fd << xyz;
     fd.dump("example_fastfitting");

     // Writing points as a paraview file:
     gsWriteParaviewPoints(xyz, "xyz");

     return 0;
}
