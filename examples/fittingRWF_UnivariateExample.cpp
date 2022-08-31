#include <gismo.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsHSplines/gsHFitting.h>
#include <gsModeling/gsFittingRWF.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    // Options with default values
    bool save     = false;
    bool condcheck = false;
    bool errguided    = false;
    bool suppguided   = false;
    bool saveLogLambda = false;
    index_t iter = 0;
    index_t numURef   = 0;
    index_t numUknots = 0;
    index_t numLRef   = 0;
    index_t numLknots = 0;
    index_t deg_x     = 2;
    index_t deg_l     = 0;
    index_t dreg      = 2;
    real_t tolerance = 1e-05;
    real_t toll = 1e-05;
    real_t alpha = std::sqrt(1e-5);
    real_t tol_slope = 18;

    real_t l_min = 1e-12;
    real_t l_max = 1e-08;
    real_t u_min = 1;
    real_t u_max = -1;
    std::string fn = "fitting/example_1D.xml";
    std::string lambda_file = "";

    // Reading options from the command line
    gsCmdLine cmd("Fit parametrized sample data with a BSpline Curve. Expected input file is an XML "
            "file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 1 x N matrix. "
            "The column represents a u parametric coordinate\nMatrix id 1 : contains a "
            "1 x N matrix. The column represents a point-value x.");
    cmd.addSwitch("save", "Save result in XML format", save);
    cmd.addSwitch("cc", "switches on condition number check of assembled System-Matrix", condcheck);
    cmd.addSwitch("sg", "Support guided method", suppguided);
    cmd.addSwitch("eg", "Error guided method", errguided);
    cmd.addInt("x", "deg_x", "degree in x direction", deg_x);
    cmd.addInt("i", "iter", "number of iterations", iter);
    cmd.addInt("y", "deg_l", "degree of lambda function", deg_l);
    cmd.addInt("r", "urefine", "initial uniform refinement steps", numURef);
    cmd.addInt("b", "uknots", "initial uniform refinement (number of inserted knots) of fitting basis in u knot direction", numUknots);
    cmd.addReal("e", "tolerance", "error tolerance (desired upper bound for pointwise error)", tolerance);
    cmd.addReal("f", "toleranceslope", "slope tolerance", tol_slope);

    cmd.addReal("g","tolerancelambda","error tolerance for lambda reduction/adaptation",toll);
    cmd.addReal("p","FactorReductionLambda", "in each lambda iteration coefficients with errors above toll will be reduced by alpha",alpha);

    cmd.addString("d", "data", "Input sample data", fn);
    cmd.addReal("u", "u_min", "u_min in uv box", u_min);
    cmd.addReal("w", "u_max", "u_max in uv box", u_max);
    cmd.addReal("z", "l_min", "for support--guided, minimum lambda in hole", l_min);
    cmd.addReal("o", "l_max", "for support--guided, maximum lambda in hole", l_max);
    cmd.addString("l", "lambda", "regularization function lambda", lambda_file);
    cmd.addInt("a", "lrefine", "initial uniform refinement steps of lambda", numLRef);
    cmd.addInt("c", "lknots", "initial uniform refinement (number of inserted knots) of lambda basis in u knot direction", numLknots);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (deg_x < 1)
    { gsInfo << "Degree must be positive.\n";  return 0;}

    if ( tolerance < 0 )
    {
        gsInfo << "Error tolerance cannot be negative, setting it to default value.\n";
        tolerance = 1e-02;
    }

    //! [Read data]
    gsFileData<> fd_in(fn);
    gsMatrix<> params, points;
    fd_in.getId<gsMatrix<> >(0, params );
    fd_in.getId<gsMatrix<> >(1, points);
    //! [Read data]

    // This is for outputing an XML file, if requested
    gsFileData<> fd;
    gsFileData<> fe;

    // Check if matrix sizes are OK
    GISMO_ASSERT( params.cols() == points.cols() && params.rows() == 1 && points.rows() <= 2,
                  "Wrong input");

    // Determine the parameter domain by mi/max of parameter values
    u_min = std::min( u_min, params.row(0).minCoeff());
    u_max = std::max( u_max, params.row(0).maxCoeff());

    // Create knot-vectors without interior knots
    gsKnotVector<> knots (u_min, u_max, numUknots, deg_x+1 ) ;

    // Create a tensor-basis nad apply initial uniform refinement
    gsBSplineBasis<> basis( knots );
    basis.uniformRefine( (1<<numURef)-1 );
    gsFittingRWFErrorGuided<1, real_t> ref( params, points, basis);

    gsTensorBSpline<1,real_t> lambda;
    // Regularization function lambda for surface fitting
    if (!(lambda_file==""))
    {
        gsInfo << "Reading input lambda tensor-product B-spline." << std::endl;
        gsFileData<> lbd_in(lambda_file);
        lbd_in.getId<gsTensorBSpline<1,real_t> >(0, lambda );
    }
    else
    {
        gsKnotVector<> knots_l (u_min,u_max,numLknots,1);
        gsBSplineBasis<> l_basis(knots_l);
        gsMatrix<> coefs (numLknots+1,1);
        for(index_t i=0; i<coefs.rows(); i++)
            coefs(i, 0) = l_max;
        lambda = gsTensorBSpline<1>(l_basis,coefs);
    }

    // Setting the correct theoretical degree, in case deg_l is not a input
    if(suppguided)  // For strategy 6 we need q=0 such that degree is not elevated
    {
        if (deg_l == -1)
        {
            deg_l = 2*deg_x +2*dreg + 1;  // Optimal degree theory: 2q + 2 dreg + d, in Paper dreg = 2
        }
        lambda.degreeElevate(deg_l-lambda.basis().degree(0));    // Could be extended to both directions
        lambda.uniformRefine((1<<numLRef)-1);
        ref.findLambda(lambda,l_min,l_max);
    }
    else if (errguided)
        lambda.uniformRefine((1<<numLRef)-1);

    const std::vector<real_t> & errors = ref.pointWiseErrors();

    // Print settings summary
    gsInfo<<"Fitting "<< points.cols() <<" samples.\n";
    gsInfo<<"----------------\n";
    gsInfo<<"Initial uniform refinement: "<< numURef<<".\n";
    gsInfo<<"Degree                    : "<< deg_x << ".\n";
    gsInfo<<"Error tolerance           : "<< tolerance<<".\n";
    gsInfo<<"Smoothing function        : "<< "Created BSpline curve of degree "<<lambda.basis().degree(0)<<".\n";

    gsStopwatch time;

    if (errguided)
    {
        real_t h_rel  = lambda.basis().getMaxCellLength();
        alpha = std::pow(alpha, h_rel);
        iter = std::floor(std::log(l_min/l_max)/std::log(alpha))+1;   // 1e-02 should be determined from l_min, l_max - TODO
        gsInfo << "----------------\n"
               << "Iteration number calculated:" << iter << std::endl
               << "Alpha as alpha^h new calculated: " << alpha << std::endl
               << "Tolerance for lambda reduction: " << toll << std::endl;
    }
    else
        iter = iter + 1;

    gsMatrix<real_t> table(7,iter);
    table.setZero();

    for (int i=0; i < iter; i++)
    {
        gsInfo<<"----------------\n";
        gsInfo<<"Iteration "<<i<<".."<<"\n";

        time.restart();
        if (suppguided)
            ref.gsFittingRWF<1,real_t>::nextIteration(alpha,toll,tolerance, lambda, true, condcheck, dreg, saveLogLambda);
        else if (errguided)
            ref.nextIteration(alpha,toll,tolerance, lambda, false, true, condcheck, dreg, saveLogLambda);
        else // given (input) lambda
            ref.gsFittingRWF<1,real_t>::nextIteration(alpha,toll,tolerance, lambda, true, condcheck, dreg, saveLogLambda);

        time.stop();

        gsInfo<<"Fitting time: "<< time <<"\n";
        gsInfo<<"Fitted with "<< ref.result()->basis() <<"\n";
        gsInfo<<"Min distance : "<< ref.minPointError() <<" / ";
        gsInfo<<"Max distance : "<< ref.maxPointError() <<"\n";
        gsInfo<<"Points below tolerance: "<< 100.0 * ref.numPointsBelow(tolerance)/errors.size()<<"%.\n";
        gsInfo<<"Span length: " << ref.result()->basis().getMinCellLength() <<".\n";
        gsInfo<<"Regularisation function: " << lambda << ".\n";

        table(0,i) = ref.result()->basis().getMinCellLength();
        table(1,i) = ref.maxPointError();
        table(2,i) = ref.getL2ApproxErrorMidpointUniform(errors);
        table(3,i) = ref.getL2ApproxErrorMidpoint(errors);
        table(4,i) = 1.0* ref.numPointsBelow(tolerance)/errors.size();
        table(5,i) = ref.result()->basis().size();
        table(6,i) = ref.getRMSE(errors);

        if ( ref.maxPointError() < tolerance )
            break;
    }

    gsInfo<<"----------------\n";
    gsInfo<<"Writing result_sandra.\n";
    gsWriteParaview(*(ref.result()),"result_sandra_deg"+std::to_string(deg_x)+"_deg"+std::to_string(deg_l), 1000, false, true);

    //ref.writeParaviewDeriv2Log(*(ref.result()),"deriv2log",10000);
    //gsInfo<<"Writing lambda_result.\n";
    //gsWriteParaview(lambda,"lambda_result", 1000, false, true);

    //gsFileData<> lambdaout;
    //lambdaout << lambda;
    //lambdaout.dump("lambda_result");

    // Use log for CP
    //gsBSpline<real_t> lambda_log (lambda);
    //for (index_t i = 0; i<lambda_log.basis().size(); i++)
    //{
    //   real_t tmp = lambda_log.coef(i)(0);
    //    //lambda_log.coef(i)(0)=mpfr::log10(tmp);
    //    lambda_log.coef(i)(0)=math::log10(tmp);
    //}

    //gsWriteParaview(lambda_log,"lambda_log_result", 1000, false, true);

    // Writing table for errorplot and saving table to errorplot.xml
    gsMatrix<real_t> bigMatrix(3, iter);
    for(index_t i=0; i<iter; i++)
    {
        bigMatrix(0,i) = table(0,i);
        bigMatrix(1,i) = table(1,i);
        bigMatrix(2,i) = table(2,i);
    }
    gsInfo<<"Writing errorplot.xml\n";
    fe << bigMatrix;
    fe.dump("L2_errorplot_withHole_deg"+std::to_string(deg_x)+"_deg"+std::to_string(deg_l));

    // Plot Errors
    //const std::vector<real_t>& eval_field = ref.pointWiseErrors();
    //gsMatrix<real_t> bigM(3, eval_field.size());
    //for(size_t i=0; i<eval_field.size(); i++)
    //{
    //    bigM(0,i) = params(0,i);
    //    bigM(1,i) = 0;
    //    bigM(2,i) = math::log10(eval_field[i]);
    //}
    //gsWriteParaviewPoints(bigM, "errorplot");

    if ( save )
    {
        gsInfo<<"Writing fitting_sandra.xml"<<"\n";
        fd << *ref.result() ;
        fd.dump("fitting_sandra");
    }

    gsInfo << "----------------\n";
    gsInfo << "----------------\n"
           << "Iteration number calculated:" << iter << std::endl
           << "Alpha as alpha^h new calculated: " << alpha << std::endl;
    gsInfo<<"Tolerance for lambda reduction: " << toll << std::endl;
    gsInfo<<"Factor tau for lambda reduction: " << alpha << std::endl;

    gsInfo<<"----------------\n";
    gsInfo<<"Table: Span length and errors for lambda function with degree "<<lambda.basis().degree()<<" and Span length "<<lambda.basis().getMaxCellLength()<<".\n";
    gsInfo<<" Line 0: Span length.\n Line 1: Maximum error.\n Line 2: L² error approx uniform.\n"
            " Line 3: L² error approx midpoint.\n Line 4: % below tol.\n Line 5: DOFs.\n"
            " Line 6: RMSE.\n";
    gsInfo<<table<<"\n";
    gsInfo<<"----------------\n Finished.\n";

    real_t lambdamax = 0e-0;
    for (int i=0 ; i < lambda.basis().size() ; i++)
    {
        if(lambda.coef(i,0)>lambdamax)
            lambdamax = lambda.coef(i,0);
    }
    real_t lambdamin = 1.0;
    for (int i=0 ; i < lambda.basis().size() ; i++)
    {
        if(lambda.coef(i,0)<lambdamin)
            lambdamin = lambda.coef(i,0);
    }

    gsInfo << "Lambda Max: " << lambdamax << std::endl;
    gsInfo << "Lambda Min: " << lambdamin << std::endl;

    return 0;
}
