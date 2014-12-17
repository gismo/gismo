
#include <iostream>
#include <algorithm>
#include <time.h>

#include <gismo.h>


using std::cout;
using namespace gismo;

int main(int argc, char *argv[])
{
    // Options with default values
    bool save     = false;
    int numURef   = 3;
    int iter      = 2;
    int deg_x     = 2;
    int deg_y     = 2;
    double lambda = 1e-07;
    double threshold = 1e-02;
    double tolerance = 1e-02;
    int extension = 2;
    double refPercent = 0.1;
    std::string fn = GISMO_DATA_DIR "/fitting/deepdrawingC.xml";

    try
    {   // Reading options from the command line
        gsCmdLine cmd("Fit parametrized sample data with a surface patch. Expected input file is an XML file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a 3 x N matrix. Every column represents a point (x,y,z) in space.");
        gsArgSwitch arg_save("", "save", "Save result in XML format", cmd); 
        gsArgVal<int> arg_iter("i", "iter", "number of iterations", false, iter, "int", cmd);
        gsArgVal<int> arg_deg_x("x", "deg_x", "degree in x direction",false,deg_x,"int",cmd);
        gsArgVal<int> arg_deg_y("y", "deg_y", "degree in y direction",false,deg_y,"int",cmd);
        gsArgVal<> arg_lambda("s", "lambda", "smoothing coefficient",false,lambda,"double",cmd);
        gsArgVal<> arg_err_threshold("t", "threshold", "error threshold (special valule -1)",false, threshold, "double", cmd);
        gsArgVal<> arg_perc("p", "refPercent", "percentage of points to refine in each iteration", false, refPercent, "double", cmd);
        gsArgVal<int> arg_extension("q", "extension", "extension size", false, extension, "int", cmd);
        gsArgVal<int> arg_uref("r", "urefine", "initial uniform refinement steps",
                               false, numURef, "int", cmd);
        gsArgVal<> arg_tolerance("e", "tolerance",  "error tolerance (desired upper bound for pointwise error)", false, tolerance, "double", cmd);
        gsArgVal<std::string> arg_data("d", "data", "Input sample data", false, fn, "double", cmd);

        cmd.parse(argc,argv);
        iter =  arg_iter.getValue();
        deg_x  =  arg_deg_x.getValue();
        deg_y = arg_deg_y.getValue();
        numURef = arg_uref.getValue();
        lambda = arg_lambda.getValue();
        threshold = arg_err_threshold.getValue();
        fn  = arg_data.getValue();
        extension = arg_extension.getValue();
        refPercent = arg_perc.getValue();
        //if ( extension == -1 )
        // {
        //    extension = deg_x +1;
        // }
        tolerance = arg_tolerance.getValue();
        save = arg_save.getValue();

        if (deg_x < 1)
        { cout << "Degree x must be positive.\n";  return 0;} //throw TCLAP::ExitException(0);
        if (deg_y < 1)
        { cout << "Degree y must be positive.\n"; return 0;}
        if (extension < 0)
        { cout << "Extension must be non negative.\n"; return 0;}

        if ( tolerance < 0 )
        { 
            cout << "Error tolerance cannot be negative, setting it to default value.\n";
            tolerance = 1e-02;
        }

        if (threshold > 0 && threshold > tolerance )
        { 
            cout << "Refinement threshold is over tolerance, setting it the same as tolerance.\n";
            threshold = tolerance;
        }

    } catch ( gsArgException& e )
    { std::cout << "Error: " << e.error() << " " << e.argId() <<"\n"; return -1; }

    // Surface fitting
    // Expected input is a file with matrices with:
    // id 0:  u,v   -- parametric coordinates, size 2 x N
    // id 1:  x,y,z -- corresponding mapped values, size 3 x N
    gsFileData<> fd_in(fn);
    gsMatrix<> uv      = safe( fd_in.getId<gsMatrix<> >(0) );
    gsMatrix<> xyz     = safe( fd_in.getId<gsMatrix<> >(1) );

    // This is for outputing an XML file, if requested
    gsFileData<> fd;

    // Check if matrix sizes are OK
    GISMO_ASSERT( uv.cols() == xyz.cols() && uv.rows() == 2 && xyz.rows() == 3,
                  "Wrong input");

    // Determine the parameter domain by mi/max of parameter values
    double u_min = uv.row(0).minCoeff(),
        u_max = uv.row(0).maxCoeff(),
        v_min = uv.row(1).minCoeff(),
        v_max = uv.row(1).maxCoeff();

    // Create knot-vectors without interior knots
    gsKnotVector<> u_knots (u_min, u_max, 0, deg_x+1 ) ;
    gsKnotVector<> v_knots (v_min, v_max, 0, deg_y+1 ) ;

    // Create a tensor-basis nad apply initial uniform refinement    
    gsTensorBSplineBasis<2> T_tbasis( u_knots, v_knots );
    T_tbasis.uniformRefine( (1<<numURef)-1 );

    // Create Initial hierarchical basis
    
    gsTHBSplineBasis<2>  THB ( T_tbasis , iter+1) ;
    //gsHBSplineBasis<2>  THB ( T_tbasis , iter+1) ;

    // Specify extension size in u and v cells
    std::vector<int> ext;
    ext.push_back(extension);
    ext.push_back(extension);

    // Create hierarchical refinement object
    gsHFitting<real_t> ref( uv, xyz, THB, refPercent, ext, lambda);
    
    const std::vector<real_t> & errors = ref.pointWiseErrors();

    // Print settings summary
    cout<<"Fitting "<< xyz.cols() <<" samples.\n";
    cout<<"----------------\n";
    cout<<"Cell extension     : "<< ext[0]<<" "<<ext[1]<<".\n";
    if ( threshold >= 0.0 )
        cout<<"Ref. threshold     : "<< threshold<<".\n";
    else
        cout<<"Cell refinement    : "<< 100*refPercent<<"%%.\n";
    cout<<"Error tolerance    : "<< tolerance<<".\n";
    cout<<"Smoothing parameter: "<< lambda<<".\n";

    gsStopwatch time;

    for(int i = 0; i <= iter; i++)
    {
        cout<<"----------------\n";
        cout<<"Iteration "<<i<<".."<<"\n";
        
        time.restart();
        ref.nextIteration(tolerance, threshold);
        const double clock = time.stop();
        cout<<"Fitting time: "<< clock <<"\n";

        cout<<"Fitted with "<< ref.result()->basis() <<"\n";
        cout<<"Min distance : "<< ref.minPointError() <<" / ";
        cout<<"Max distance : "<< ref.maxPointError() <<"\n";
        cout<<"Points below tolerance: "<< 100.0 * ref.numPointsBelow(tolerance)/errors.size()<<"%.\n";

        if ( ref.maxPointError() < tolerance )
        {
            cout<<"Error tolerance achieved after "<<i<<" iterations.\n";
            break;
        }
    }

    cout<<"----------------\nFinished.\n";

    if ( save )
    {
        cout<<"Writing fitting_out.xml"<<"\n";
        fd << *ref.result() ;

        fd.dump("fitting_out");
    }

    return 0;
}

