/// gsFitting_general.cpp
/// Author:Gabor Kiss
/// for testing the class gsFitting

#include <iostream>
#include <time.h>

#include <gismo.h>




using namespace gismo;


int main(int argc, char* argv[])
{
  
    int nd = 100;
    int function = 1;
    real_t eps = 0.0;
    std::string filename = "fitting/face.xml";


    gsCmdLine cmd("Hi, give me a file (.xml) with some multipatch geometry and basis.");
    cmd.addString("g", "filename" , "File containing mp geometry and basis (.xml).", filename);
    cmd.addInt("f", "function", "Number of the function", function);
    cmd.addInt("p", "nd", "Number of fitting points", nd);
    
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // The file data
    gsFileData<> data(filename);
    gsFunctionExpr<> f;
    gsMultiPatch<> mp;
    

    // Load data as multipatch structure

      data.getFirst(mp);
   

    switch (function) {
    case 1:
        f = gsFunctionExpr<>("(exp(52*sqrt((10*x-2)^2+(10*y-3)^2)))^(-1)", 3);//one peak
        gsInfo << "Source function: " << f << "\n";
        break;
    case 2:
        f = gsFunctionExpr<>("y/2*(cos(4*(x^2+y-1)))^4", 3);//waves
        gsInfo << "Source function: " << f << "\n";
        break;
    case 3:
        f = gsFunctionExpr<>("3/4*exp(-((9*x-2)^2 + (9*y-2)^2)/4)+3/4*exp(-((9*x+1)^2)/49 - (9*y+1)/10)+1/2*exp(-((9*x-7)^2 + (9*y-3)^2)/4)-1/5*exp(-(9*x-4)^2 - (9*y-7)^2)", 3);//hills
        gsInfo << "Source function: " << f << "\n";
        break;
    case 4:
        f = gsFunctionExpr<>("sin((x*x+y*y+2/(5*pi))^(-1))", 3);//crater
        gsInfo << "Source function: " << f << "\n";
        break;
    case 5:
        f = gsFunctionExpr<>("(exp(2*sqrt((x)^2+(y)^2)))^(-1)", 3);//cusp
        gsInfo << "Source function: " << f << "\n";
        break;
    case 6:
        f = gsFunctionExpr<>("(exp(2*sqrt((10*x+3)^2+(10*y-3)^2)))^(-1) + (exp(2*sqrt((10*x-3)^2+(10*y+3)^2)))^(-1)", 3);//2 peaks
        gsInfo << "Source function: " << f << "\n";
        break;
    case 7:
        f = gsFunctionExpr<>("(1.5*exp(sqrt((10*x-3)^2+(10*y-3)^2)))^(-1)+ (1.5*exp(sqrt((10*x+3)^2+(10*y+3)^2)))^(-1) + (1.5*exp(sqrt((10*x)^2+(10*y)^2)))^(-1)", 3);//3 peaks
        gsInfo << "Source function: " << f << "\n";
        break;
    case 8:
        f = gsFunctionExpr<>("(x+y)/2+sqrt(((x-y)/2)^2)", 3);//Rvachev functions (R-functions)
        gsInfo << "Source function: " << f << "\n";
        break;
    case 9:
        f = gsFunctionExpr<>("x", 3);//constant
        gsInfo << "Source function: " << f << "\n";
        break;
    case 10:
        f = gsFunctionExpr<>("x^2+y^2", 3);//paraboloid
        gsInfo << "Source function: " << f << "\n";
        break;
    case 11:
        f = gsFunctionExpr<>("(tanh(9*y-9*x)+1)/9", 3);//sigmoid
        gsInfo << "Source function: " << f << "\n";
        break;
    default:
        gsInfo << "Unknown function, please pick one of the functions 1 - 8" << "\n";
        return 0;
        break;
    }

    int nbp = mp.nPatches();


    gsVector<index_t> offset(nbp+1);
    offset[0] = 0;
    for (int i = 1; i<nbp; ++i)
    {
        offset[i] = offset[i - 1] + nd;
    }
    offset[nbp] = offset[nbp - 1] + nd;


    gsMatrix<real_t> Mpar (2 , nbp * nd);
    gsMatrix<real_t> fval (3, nbp * nd);


    // loop on the nb of patches
    gsMatrix<> para, pts, mp_eval;
    gsVector<> c0, c1;
    for (int i = 0; i < nbp ; i++)
    {
        para = mp.patch(i).parameterRange();

        c0 = para.col(0);
        c1 = para.col(1);
        c0.array() += eps; // avoid double points on boundary
        c1.array() -= eps;
        //the parameter values for the fitting
        pts = uniformPointGrid(c0, c1, nd);
        //gsInfo << "Parameter values used for fitting: " << "\n" << pts << "\n";
        Mpar.middleCols(offset[i], nd) = pts;

        //the evaluated values of the original surface
        mp_eval = mp.patch(i).eval(pts);

        fval.middleCols(offset[i], nd).topRows(2) = mp_eval.topRows(2);
        fval.middleCols(offset[i], nd).row(2)     = f.eval(mp_eval);

    }


    gsFileData<> out;
    out << Mpar;
    out << fval;
    out << gsMatrix<index_t>(offset); //conversion
    out.save("function_pointcloud");

   // end loop

    return EXIT_SUCCESS;
}

   