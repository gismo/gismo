/// gsFitting_general.cpp
/// Author:Gabor Kiss
/// for testing the class gsFitting

#include <iostream>
#include <time.h>

#include <gismo.h>




using namespace gismo;


int main(int argc, char* argv[])
{
    bool plot = false;
    //measuring the computational time
    //int clo=clock();
    std::string filename = "face.xml";

    int np = 700;
    //int nd = 100;
    //int err_type = 1;
    //int function = 1;
    real_t lambda = 0;
    //real_t eps = 0.01;
    std::string fn = "fitting/deepdrawingC.xml";

    /*gsCmdLine cmd("Hi, give me a file (.xml) with some multipatch geometry and basis.");
    cmd.addPlainString("filename", "File containing mp geometry and basis (.xml).", filename);
    cmd.addInt("f", "function", "Number of the function", function);
    cmd.addReal ("l", "lambda", "smoothing parameter", lambda);
    cmd.addInt("p", "nd", "Number of fitting points", nd);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    */

    gsCmdLine cmd("Fit parametrized sample data with a surface patch. Expected input file is an XML "
        "file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. "
        "Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a "
        "3 x N matrix. Every column represents a point (x,y,z) in space.");
    cmd.addString("g", "filename", "File containing mp geometry and basis (.xml).", filename);
    cmd.addString("d", "data", "Input sample data", fn);
    cmd.addReal("l", "lambda", "smoothing parameter", lambda);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFileData<> fd_in(fn);
    gsMatrix<> Mpar, fval;
    gsMatrix<> offset;
    fd_in.getId<gsMatrix<> >(0, Mpar);
    fd_in.getId<gsMatrix<> >(1, fval);
    fd_in.getId<gsMatrix<> >(2, offset);

    // The file data
    gsFileData<> data( filename );
    gsFunctionExpr<> f;
    gsMultiPatch<> mp;
    gsMultiBasis<> mb;
    gsSparseMatrix<> cf;
    gsMappedBasis<2, real_t> mbasis;

    // Load data as multipatch structure

    data.getFirst(mp);
    data.getFirst(mb);
    data.getFirst(cf);

    mbasis.init(mb, cf);

    gsInfo << "Number of patches: " << mp.nPatches() << "\n";
    //gsInfo << "Number of points per patch: " << nd << "\n";
    gsInfo << "Number of basis functions: " << mbasis.size() << "\n";
    //gsInfo << "Basis" << mb << "\n";


    /*switch (function) {
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
    out.save("fitting_test");

   // end loop

   */

   


    //create the fitting object
    gsInfo << "//////////////////////////////////////////////////////////" << "\n";
    gsInfo << "Creating the multipatch fitting object" << "\n";
    gsInfo << "//////////////////////////////////////////////////////////" << "\n";


    gsFitting<> fitting(Mpar, fval, offset, mbasis);


    gsInfo << "Fit class created" << "\n";

    fitting.compute(lambda); // smoothing parameter lambda

    gsInfo << "I computed the fitting" << "\n";

    gsMappedSpline<2, real_t>  test;
    test = fitting.mresult();

    std::vector<real_t> errors;
    fitting.get_Error(errors, 0);

    real_t min_error = *std::min_element(errors.begin(), errors.end());
    real_t max_error = *std::max_element(errors.begin(), errors.end());

    gsInfo << "Min error: " << min_error << "\n";
    gsInfo << "L_inf: " << max_error << "\n";

    real_t error;

    fitting.computeApproxError(error, 0);

    gsInfo << "L_2: " << error << "\n";

    gsFileData<> newdata;
    gsMultiPatch<> surf = test.exportToPatches();
    newdata << surf;
    newdata.save("fitting_result");

    if (plot)
    {
        gsInfo << "Plotting in Paraview..." << "\n";
        gsWriteParaviewPoints(fval, "point_data");
        gsWriteParaview(surf, "multipatch_spline", np);
        gsFileManager::open("multipatch_spline.pvd");
    }

    return EXIT_SUCCESS;
}
