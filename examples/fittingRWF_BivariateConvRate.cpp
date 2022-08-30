#include <gismo.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsHSplines/gsHFitting.h>
#include <gsModeling/gsFittingRWF.h>
#include <gsParasolid/gsWriteParasolid.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Options with default values
    real_t alpha      = 1e-10;
    real_t eps        = 1e-06;
    real_t epsmin     = 0.0;
    int dreg          = 2;
    bool iter_lambda  = false;
    bool iter_refine  = false;
    bool save         = false;
    bool condcheck    = false;
    bool saveLogLambda = false;
    index_t numURef   = 0;
    index_t numUknots = 0;
    index_t numVknots = 0;
    index_t iter      = 2;
    index_t deg_x     = 2;
    index_t deg_y     = 2;
    index_t deg_l     = -1;
    index_t numLRef   = 6;
    real_t tolerance  = 1e-10;
    real_t toll       = 1e-02;
    real_t u_min      = 1;
    real_t u_max      = -1;
    real_t v_min      = 1;
    real_t v_max      = -1;

    // Reading options from the command line
    gsCmdLine cmd("Fit parametrized sample data with a surface patch. Expected input file is an XML "
            "file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. "
            "Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a "
            "3 x N matrix. Every column represents a point (x,y,z) in space.");
    cmd.addSwitch("save", "Save result in XML format", save);
    cmd.addSwitch("iter_l", "Iterate lambda", iter_lambda);
    cmd.addSwitch("iter_r", "Refine solution", iter_refine);
    cmd.addSwitch("cc", "switches on condition number check of assembled System-Matrix", condcheck);
    cmd.addSwitch("loglambda", "Save logarithm of lambda in PVD format", saveLogLambda);
    cmd.addInt("i", "iter", "number of iterations", iter);
    cmd.addInt("x", "deg_x", "degree in x direction", deg_x);
    cmd.addInt("y", "deg_y", "degree in y direction", deg_y);
    cmd.addInt("q", "deg_l","degree lambda", deg_l);
    cmd.addReal("f", "eps","for strategy 7, maximum lambda in hole", eps);
    cmd.addReal("k", "epsmin","for strategy 7, minimum lambda in hole", epsmin);

    cmd.addInt("n", "numderiv", "number of which derivativeis used for regularization", dreg);
    cmd.addInt("r", "urefine", "initial uniform refinement steps of the fitting basis", numURef);
    cmd.addInt("a", "lrefine", "initial uniform refinement steps of lambda", numLRef);
    cmd.addInt("b", "uknots", "initial uniform refinement (number of inserted knots) of fitting basis in u knot direction", numUknots);
    cmd.addInt("c", "vknots", "initial uniform refinement (number of inserted knots) of fitting basis in v knot direction", numVknots);
    cmd.addReal("e", "tolerance", "error tolerance (desired upper bound for pointwise error)", tolerance);

    cmd.addReal("g","tolerancelambda","error tolerance for lambda reduction/adaptation",toll);
    cmd.addReal("u", "u_min", "u_min in uv box", u_min);
    cmd.addReal("v", "v_min", "v_min in uv box", v_min);
    cmd.addReal("w", "u_max", "u_max in uv box", u_max);
    cmd.addReal("z", "v_max", "v_max in uv box", v_max);

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

    // This is for outputing an XML file, if requested
    gsFileData<> fd;
    gsFileData<> fe;

    // Creating the data by sampling from the function graph f(u,v)=1/3*sin(4pi*(u+1/8))sin(4pi*(v+1/8))
    // mesh size 1/nn, domain [-1,1]^2, hole ]-0.5,0.5[^2
    // Initializing boundary hole and number of samples
    real_t hv_min = -0.5;
    real_t hv_max = 0.5;
    real_t hu_min = -0.5;
    real_t hu_max = 0.5;
    index_t nn = 2000;
    gsMatrix<real_t> uv(2,(2*nn+1)*(2*nn+1)-(nn-1)*(nn-1));
    gsMatrix<real_t> xyz(3,(2*nn+1)*(2*nn+1)-(nn-1)*(nn-1));

    gsInfo << "Creating " << (nn+1)*(nn+1) <<" data points" << std::endl;
    index_t ii = 0;
    for(index_t i=-nn; i<=nn; i++)
    {
        for(index_t j=-nn; j<=nn; j++)
        {
            if (!(hu_min < (1.0*i)/(1.0*nn) && (1.0*i)/(1.0*nn) < hu_max && hv_min < (1.0*j)/(1.0*nn) && (1.0*j)/(1.0*nn) < hv_max ))
            {
                uv(0,ii) = (1.0*i)/(1.0*nn);
                uv(1,ii) = (1.0*j)/(1.0*nn);
                xyz(0,ii) = uv(0,ii);
                xyz(1,ii) = uv(1,ii);
                xyz(2,ii) = 1.0/3.0*math::sin(4*M_PI*(uv(0,ii)+0.1/8.0))*math::sin(4*M_PI*(uv(1,ii)+0.1/8.0));
                ii += 1;
            }
        }
    }

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

    //! [Read data]
    // Regularization function lambda for surface fitting, Initialization by input lambda degree
    gsKnotVector<> u_lknots (u_min, u_max, 0, deg_l+1 ) ;
    gsKnotVector<> v_lknots (v_min, v_max, 0, deg_l+1 ) ;
    gsMatrix<real_t> lambda_coeffs ((deg_l+1)*(deg_l+1),1);
    lambda_coeffs.setZero();
    gsTensorBSpline<2> lambda (u_lknots, v_lknots,lambda_coeffs);
    lambda.uniformRefine((1<<numLRef)-1);

    // Create a tensor-basis and apply initial uniform refinement
    gsTensorBSplineBasis<2> T_tbasis( u_knots, v_knots );
    T_tbasis.uniformRefine( (1<<numURef)-1 );

    // Create hierarchical refinement object
    gsFittingRWFErrorGuided<2, real_t> ref( uv, xyz, T_tbasis);
    gsInfo << "Fitting Basis: " << ref.getBasis() << std::endl;

    ref.findLambda(lambda, epsmin , eps);

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
    gsMatrix<real_t> table(6,iter+1);
    table.setZero();

    for(int i = 0; i <= iter; i++)
    {
        gsInfo<<"----------------\n";
        gsInfo<<"Iteration "<<i<<".."<<"\n";

        time.restart();
        ref.gsFittingRWF<2,real_t>::nextIteration(alpha,toll,tolerance, lambda, false, condcheck, dreg, saveLogLambda);//nextLambdaIteration(alpha,toll,tolerance, threshold, lambda, true, false, false, dreg);
        time.stop();
        gsInfo<<"Fitting time: "<< time <<"\n";
        gsInfo<<"Fitted with "<< ref.result()->basis() <<"\n";
        gsInfo<<"Min distance : "<< ref.minPointError() <<" / ";
        gsInfo<<"Max distance : "<< ref.maxPointError() <<"\n";
        gsInfo<<"Points below tolerance: "<< 100.0 * ref.numPointsBelow(tolerance)/errors.size()<<"%."<< std::endl;

        table(0,i) = ref.result()->basis().getMinCellLength();
        table(1,i) = ref.maxPointError();
        table(2,i) = ref.getL2ApproxErrorMidpointUniform(errors);
        table(3,i) = 1.0* ref.numPointsBelow(tolerance)/errors.size();
        table(4,i) = ref.result()->basis().size();
        table(5,i) = ref.getRMSE(errors);

        if ( ref.maxPointError() < tolerance )
            break;
    }

    // Writing table for errorplot and saving table to errorplot.xml for Matlab convergence plot
    gsInfo<<"Writing error*.xml\n";
    fe << table;
    fe.dump("error_deg"+std::to_string(deg_x)+"_deg"+std::to_string(deg_l));

    //gsInfo<<"----------------\n";
    //gsInfo<<"Strategy "<<strategy<<" \n";
    gsInfo<<"----------------\n";
    gsInfo<<"Table: Span length and errors for lambda function with degree "<<lambda.degree(0)<<" and Span length "<<lambda.basis().getMaxCellLength()<<".\n";
    gsInfo<<" Line 0: Span length.\n Line 1: Maximum error.\n Line 2: L^2 error.\n Line 3: Pts below tol.\n Line 4: DoFs.\n Line 5: RMSE.\n";
    gsInfo<<table<<"\n";
    gsInfo<<"----------------\nFinished.\n";

    // Plot Errors
    /*const std::vector<real_t>& eval_field = ref.pointWiseErrors();
    gsMatrix<real_t> bigMatrix(4, eval_field.size());
    for(size_t i=0; i<eval_field.size(); i++)
    {
        bigMatrix(0,i) = uv(0,i);
        bigMatrix(1,i) = uv(1,i);
        bigMatrix(2,i) = 0;
        bigMatrix(3,i) = eval_field[i];
    }
    gsWriteParaviewPoints(bigMatrix, "errorplot");*/

    // Save result in pvd and xml format
    gsWriteParaview(*(ref.result()),"result_deg"+std::to_string(deg_x)+"_deg"+std::to_string(deg_l), 50000, false, true);
    if ( save )
    {
        gsInfo<<"Writing fitting_sandra.xml"<<"\n";

        // Writing xml as a TensorBSpline not as a THB:
        gsTensorBSpline<2>* TP_spline = static_cast<gsTensorBSpline<2>*>(ref.result());
        fd << *TP_spline;
        fd.dump("result");
    }

    // Saving last lambda:
    ref.writeParaviewLog(lambda,"lambda", 50000, false);

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
