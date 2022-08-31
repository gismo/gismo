#include <gismo.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsHSplines/gsHFitting.h>
#include <gsModeling/gsFittingRWF.h>
#include <gsParasolid/gsWriteParasolid.h>

#include <iostream>
#include <fstream>



using namespace gismo;

int main(int argc, char *argv[])
{
    // Options with default values
    real_t alpha      = 1e-10;
    real_t eps        = 1e-06;
    real_t epsmin     = 0.0;
    int dreg          = 2;
    //bool iter_lambda  = false;
    //bool iter_refine  = false;
    bool save         = false;
    bool condcheck    = false;
    bool errguided    = false;
    bool suppguided   = false;
    bool saveLogLambda = false;
    index_t numURef   = 0;
    index_t numUknots = 0;
    index_t numVknots = 0;
    index_t numLknots = 0;
    // iter: 0 -> only a fitting with given settings, >0 -> number of iterations with uniform refinement
    index_t iter      = 2;
    index_t deg_x     = 2;
    index_t deg_y     = 2;
    index_t deg_l     = -1;
    index_t numLRef   = 6;
    real_t tolerance  = 1e-02;
    real_t toll       = 1e-02;
    real_t u_min      = 1;
    real_t u_max      = -1;
    real_t v_min      = 1;
    real_t v_max      = -1;
    std::string fn    = "fitting/deepdrawingC.xml";
    std::string lambda_file = "";
    //real_t ltemp      = 1.0;

    // Reading options from the command line
    gsCmdLine cmd("Fit parametrized sample data with a surface patch. Expected input file is an XML "
            "file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. "
            "Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a "
            "3 x N matrix. Every column represents a point (x,y,z) in space.");
    cmd.addSwitch("save", "Save result in XML format", save);
    cmd.addSwitch("loglambda", "Save logarithm of lambda in PVD format", saveLogLambda);
    cmd.addSwitch("cc", "switches on condition number check of assembled System-Matrix", condcheck);
    cmd.addSwitch("sg", "Support guided method", suppguided);
    cmd.addSwitch("eg", "Error guided method", errguided);
    cmd.addInt("i", "iter", "number of iterations", iter);
    cmd.addInt("x", "deg_x", "degree in x direction", deg_x);
    cmd.addInt("y", "deg_y", "degree in y direction", deg_y);
    cmd.addInt("q", "deg_l","degree lambda", deg_l);
    cmd.addReal("f", "eps","for support--guided, maximum lambda in hole", eps);
    cmd.addReal("k", "epsmin","for support--guided, minimum lambda in hole", epsmin);
    cmd.addReal("j", "alpha","for error--guided, maximal step alpha^h", alpha);

    cmd.addInt("n", "numderiv", "number of which derivativeis used for regularization", dreg);
    cmd.addInt("r", "urefine", "initial uniform refinement steps of the fitting basis", numURef);
    cmd.addInt("a", "lrefine", "initial uniform refinement steps of lambda", numLRef);
    cmd.addInt("b", "uknots", "initial uniform refinement (number of inserted knots) of fitting basis in u knot direction", numUknots);
    cmd.addInt("c", "vknots", "initial uniform refinement (number of inserted knots) of fitting basis in v knot direction", numVknots);
    cmd.addInt("m", "lknots", "initial uniform refinement (number of inserted knots) of lambda basis in u and v direction", numLknots);
    cmd.addReal("e", "tolerance", "error tolerance (desired upper bound for pointwise error)", tolerance);

    cmd.addReal("g","tolerancelambda","error tolerance for lambda reduction/adaptation",toll);
    cmd.addString("d", "data", "Input sample data", fn);
    cmd.addReal("u", "u_min", "u_min in uv box", u_min);
    cmd.addReal("v", "v_min", "v_min in uv box", v_min);
    cmd.addReal("w", "u_max", "u_max in uv box", u_max);
    cmd.addReal("z", "v_max", "v_max in uv box", v_max);
    cmd.addString("l", "lambda", "regularization function lambda", lambda_file);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (deg_x < 1)
    { gsInfo << "Degree x must be positive.\n";  return 0;}
    if (deg_y < 1)
    { gsInfo << "Degree y must be positive.\n"; return 0;}

    if ( tolerance < 0 )
    {
        gsInfo << "Error tolerance cannot be negative, setting it to default value.\n";
        tolerance = 1e-02;
    }

    //! [Read data]
    // Surface fitting
    // Expected input is a file with matrices with:
    // id 0:  u,v   -- parametric coordinates, size 2 x N
    // id 1:  x,y,z -- corresponding mapped values, size 3 x N
    gsFileData<> fd_in(fn);
    gsMatrix<> uv, xyz;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, xyz);
    //! [Read data]

    // This is for outputing an XML file, if requested
    gsFileData<> fd;
    gsFileData<> fe;

    // Check if matrix sizes are OK
    GISMO_ASSERT( uv.cols() == xyz.cols() && uv.rows() == 2 && xyz.rows() == 3,
                  "Wrong input");

    // Determine the parameter domain by mi/max of parameter values
    u_min = std::min( u_min, uv.row(0).minCoeff());
    u_max = std::max( u_max, uv.row(0).maxCoeff());
    v_min = std::min( v_min, uv.row(1).minCoeff());
    v_max = std::max( v_max, uv.row(1).maxCoeff());

    // Create knot-vectors without interior knots
    gsKnotVector<> u_knots (u_min, u_max, numUknots, deg_x+1 ) ;
    gsKnotVector<> v_knots (v_min, v_max, numVknots, deg_y+1 ) ;

    // Regularization function lambda for surface fitting
    gsTensorBSpline<2> lambda;
    if (lambda_file=="")
    {
        gsKnotVector<> uknots_l (u_min,u_max,numLknots,1);
        gsKnotVector<> vknots_l (v_min,v_max,numLknots,1);
        gsTensorBSplineBasis<2> l_basis(uknots_l,vknots_l);
        gsMatrix<> coefs ((numLknots+1)*(numLknots+1),1);
        for(index_t i=0; i<coefs.rows(); i++)
            coefs(i, 0) = eps;
        lambda = gsTensorBSpline<2>(l_basis,coefs);
    }
    else
    {
        gsInfo << "Reading input lambda tensor-product B-spline." << std::endl;
        gsFileData<> lbd_in(lambda_file);
        lbd_in.getId<gsTensorBSpline<2> >(0, lambda );
    }

    // Setting the correct theoretical degree, in case deg_l is not a input
    if(suppguided)  // For strategy 6 we need q=0 such that degree is not elevated
    {
        if (deg_l == -1)
        {
            deg_l = 2*(std::max(deg_x,deg_y)) +2*dreg + 2;  // Optimal degree theory: 2q + 2 dreg + d, in Paper dreg = 2
        }
        lambda.degreeElevate(deg_l-lambda.degree(0));    // Could be extended to both directions
        lambda.uniformRefine((1<<numLRef)-1);
    }
    else if (errguided)
        lambda.uniformRefine((1<<numLRef)-1);

    // Create a tensor-basis and apply initial uniform refinement
    gsTensorBSplineBasis<2> T_tbasis( u_knots, v_knots );
    T_tbasis.uniformRefine( (1<<numURef)-1 );

    // Create hierarchical refinement object
    gsFittingRWFErrorGuided<2, real_t> ref( uv, xyz, T_tbasis);

    gsInfo << "Fitting Basis: " << ref.getBasis() << std::endl;
    gsInfo << "Lambda Basis: " << lambda << std::endl;
    // Support-guided method for the choice of lambda
    if (suppguided)
    {
        ref.findLambda(lambda, epsmin , eps);
    }

    // Adapt alpha for the lambda adaptation -for 2D is it alpha^h or alpha^h^2?
    gsInfo << "Alpha: " << alpha << std::endl;
    real_t h_rel  = lambda.basis().getMaxCellLength();
    gsInfo << "Relative h lambda" << h_rel << std::endl;

    const std::vector<real_t> & errors = ref.pointWiseErrors();

    // Print settings summary
    gsInfo<<"Fitting "<< xyz.cols() <<" samples.\n";
    gsInfo<<"----------------\n";
    gsInfo<<"Initial uniform refinement: "<< numURef<<".\n";
    gsInfo<<"Degree u/v direction      : "<< deg_x << " " << deg_y <<".\n";
    gsInfo<<"Error tolerance           : "<< tolerance<<".\n";
    gsInfo<<"Smoothing fct Refinement  : "<< numLRef  << std::endl;
    gsInfo<<"Smoothing function        : "<< lambda<<".\n";

    gsStopwatch time;

    if (errguided)
    {
        alpha = std::pow(alpha, h_rel);
        iter = std::floor(std::log(1e-02)/std::log(alpha))+1;   // 1e-02 should be determined from l_min, l_max - TODO
        gsInfo << "----------------\n"
               << "Iteration number calculated:" << iter << std::endl
               << "Alpha as alpha^h new calculated: " << alpha << std::endl
               << "Tolerance for lambda reduction: " << toll << std::endl;
    }
    else
        iter = iter + 1;

    gsMatrix<real_t> table(6,iter);
    table.setZero();

    for(int i = 0; i < iter; i++)
    {
        gsInfo<<"----------------\n";
        gsInfo<<"Iteration "<<i<<".."<<"\n";
        time.restart();

        if (suppguided)
            ref.gsFittingRWF<2,real_t>::nextIteration(alpha,toll,tolerance, lambda, true, condcheck, dreg, saveLogLambda);
        else if (errguided)
            ref.nextIteration(alpha,toll,tolerance, lambda, false, true, condcheck, dreg, saveLogLambda);
        else // given (input) lambda
            ref.gsFittingRWF<2,real_t>::nextIteration(alpha,toll,tolerance, lambda, true, condcheck, dreg, saveLogLambda);

        time.stop();
        gsInfo<<"Fitting time: "<< time <<"\n";
        gsInfo<<"Fitted with "<< ref.result()->basis() <<"\n";
        gsInfo<<"Min distance : "<< ref.minPointError() <<" / ";
        gsInfo<<"Max distance : "<< ref.maxPointError() <<"\n";
        gsInfo<<"Points below tolerance: "<< 100.0 * ref.numPointsBelow(tolerance)/errors.size()<<"%.\n";

        table(0,i) = ref.result()->basis().getMinCellLength();
        table(1,i) = ref.maxPointError();
        table(2,i) = ref.getL2ApproxErrorMidpointUniform(errors);
        table(3,i) = 1.0* ref.numPointsBelow(tolerance)/errors.size();
        table(4,i) = ref.result()->basis().size();
        table(5,i) = ref.getRMSE(errors);

        if ( ref.maxPointError() < tolerance )
            break;
    }

    gsInfo<<"----------------\n";
    gsInfo<<"Table: Span length and errors for lambda function with degree "<<lambda.degree(0)<<" and Span length "<<lambda.basis().getMaxCellLength()<<".\n";
    gsInfo<<" Line 0: Span length.\n Line 1: Maximum error.\n Line 2: L^2 error.\n Line 3: Pts below tol.\n Line 4: DoFs.\n Line 5: RMSE.\n";
    gsInfo<<table<<"\n";
    gsInfo<<"----------------\nFinished.\n";

    // Plot Errors in parameter domain
    const std::vector<real_t>& eval_field = ref.pointWiseErrors();
    gsMatrix<real_t> bigMatrix(4, eval_field.size());
    for(size_t i=0; i<eval_field.size(); i++)
    {
        bigMatrix(0,i) = uv(0,i);
        bigMatrix(1,i) = uv(1,i);
        bigMatrix(2,i) = 0;
        bigMatrix(3,i) = eval_field[i];
    }
    gsWriteParaviewPoints(bigMatrix, "errorplot");

    // New 3D writing error:
    // Writing table for errorplot and saving table to errorplot.xml
    gsInfo<<"Writing errorplot.xml\n";
    fe << table;
    fe.dump("error");

    // Save result in pvd and xml format
    gsWriteParaview(*(ref.result()),"result", 50000, false, true);
    if ( save )
    {
        gsInfo<<"Writing fitting_sandra.xml"<<"\n";

        // Writing xml as a TensorBSpline not as a THB:
        gsTensorBSpline<2>* TP_spline = static_cast<gsTensorBSpline<2>*>(ref.result());
        fd << *TP_spline;
        fd.dump("result");
    }

    bool nx = true;
    if(nx)
    {
        gsInfo << "Writing NX File" << std::endl;
        gsTensorBSpline<2>* TP_spline = static_cast<gsTensorBSpline<2>*>(ref.result());
        extensions::gsWritePK_SHEET(*TP_spline, "result");
    }

    // TODO
//    if (saveLogLambda)
//        ref.writeParaviewLog(lambda,"lambda", 50000, false);


    // Write out lambda min and max
    real_t lambdamax(0);
    real_t lambdamin(1);
    for (index_t i=0 ; i < lambda.basis().size() ; i++)
    {
        if(lambda.coef(i,0)>lambdamax)
            lambdamax = lambda.coef(i,0);
        if(lambda.coef(i,0)<lambdamin)
            lambdamin = lambda.coef(i,0);
    }
    gsInfo << "Lambda Max: " << lambdamax << std::endl;
    gsInfo << "Lambda Min: " << lambdamin << std::endl;

    return 0;
}
