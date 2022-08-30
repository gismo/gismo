#include <gismo.h>
#include <math.h>

using namespace gismo;

void sinus1D(real_t u_min, real_t u_max, int n)
{
    gsInfo << "Start example 1D Paper 1, support--guided: \n"
              "Sampling points with mesh size 1/n from sinus with a hole in [u_min,u_max]." << std::endl;

    // Initializing
    std::vector<real_t> u;
    std::vector<real_t> x;
    real_t par;

    for(int i=-n; i<=n; i++)
    {
        par = (1.0*i)/n;
        if(par <= u_min || par >= u_max)
        {
            u.push_back(par);
            x.push_back(sin(4*M_PI*(par+0.1/8.0)));
        }
    }

    gsMatrix<>  chosen_u(1,u.size());
    gsMatrix<> chosen_x(1,x.size());

    for(size_t i=0; i< u.size();i++)
    {
        chosen_u(0,i) = u[i];
        chosen_x(0,i) = x[i];
    }

    gsFileData<> fd;
    fd << chosen_u;
    fd << chosen_x;
    fd.dump("example_1D");

    // Plot u and x as data.
    gsWriteParaviewPoints<>(chosen_u, chosen_x, "u-x_1D");
}

void linecosinus1D(real_t n1, real_t n2, real_t n)
{
    gsInfo << "Start example 1D Paper 1, error--guided: \n"
              "4 points on line, then damped sinus." << std::endl;

    // Initialization
    gsMatrix<>  u(1,n1+n2);
    gsMatrix<> xy(2,n1+n2);

    // n1 points with parameters in [0,2/3[ on x-axis, with x=u.
    for (int i=0; i<n1; i++)
    {
        u(i,0) = 2.0/3.0* i/n1;
        xy(0,i) = u(i,0);
        xy(1,i) = 0.0;
    }

    // n2 points with parameters in [2/3,1], with x=u and y a damped cosinus.
    for(int i=0; i<n2; i++)
    {
        u(i+n1,0) = (1.5*(i+1))/(n2)+1.5;

        // Intervall [0,1]
        xy(0,i+n1) = (u(i+n1,0)+1.5)/4.5;
        xy(1,i+n1) = (u(i+n1,0)-1.5) * std::cos(5*(u(i+n1,0)-1.5))/4.5;
        // Random noise in the y component
        real_t noise = n*(((double)rand() / (double)RAND_MAX)-0.5)*2.0 ;
        xy(1,i+n1) = xy(1,i+n1) + noise;

        // Bringing intervall to [0,1]
        u(i+n1,0) = (u(i+n1,0)+1.5)/4.5;
    }

    gsFileData<> fd;
    fd << u;
    fd << xy;
    fd.dump("example_v_curve");

    gsWriteParaviewPoints(xy, "xy");
}

int main(int argc, char *argv[])
{
    // Initializing
    real_t u_min = -1.0/8.0;    // left boundary hole for Example 1
    real_t u_max =  1.0/8.0;    // right boundary hole for Example 1
    int n=200;                  // sampling with mesh size 1/n for Example 1

    index_t n1 = 3;             // number samples on line (without 1.point of cosinus)
    index_t n2 = 40;            // number samples on damped cosinus, 1.point (2/3,0) on line.
    real_t noise = 0.001;       // absolute value of maximal noise

    // Example 1: sinus with hole, example Paper 1, convergence univariate
    // Example 2: straight line + sinus: example Paper 1, error-guided univariate
    index_t example = 1;

    gsCmdLine cmd("Specify Example: \n"
                  "1: sinus with hole, example Paper 1, convergence univariate.\n"
                  "   u: left boundary hole, v: right boundary hole, 1/n sampling mesh size."
                  "2: straight line + sinus: example Paper 1, error-guided univariate.\n"
                  "   n: #samples on line, m: #samples on damped cosinus, r: maximal level of random noise added.");
    cmd.addReal("u","u_min","Example 1: start of hole",u_min);
    cmd.addReal("v","u_max","Example 1: end of hole",u_max);
    cmd.addInt("e", "example", "Switch between given examples", example);
    cmd.addInt("k", "n_samples", "Example 1: 1/n = mesh size for sampling", n);
    cmd.addInt("n", "NumberSamples1", "Example 2: number of samples straight part curve", n1);
    cmd.addInt("m", "NumberSamples2", "Example 2: number of samples damped cosinus part curve", n2);
    cmd.addReal("r", "rvalue" , "Maximal value of random noise added for Example 2", noise);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (example == 1)
        sinus1D(u_min, u_max, n);
    else if (example == 2)
        linecosinus1D(n1,n2,noise);

    return 0;
}
