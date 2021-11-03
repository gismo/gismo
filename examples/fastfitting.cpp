// Fast Fitting idea Bert
#include <gismo.h>
#include <gsModeling/gsFastFitting.h>
#include <gsModeling/gsFitting.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    index_t numURef   = 0;
    index_t numUknots = 0;
    index_t numVknots = 0;
    index_t deg_x     = 3;
    index_t deg_y     = 3;
    real_t u_min      = 1;
    real_t u_max      = -1;
    real_t v_min      = 1;
    real_t v_max      = -1;
    index_t n_ugrid   = 11;
    index_t n_vgrid   = 11;
    bool borders      = false;

    std::string fn    = "fitting/fasttest.xml";

    // Reading options from the command line
    gsCmdLine cmd("Fit parametrized sample data with a surface patch. Expected input file is an XML "
            "file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. "
            "Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a "
            "3 x N matrix. Every column represents a point (x,y,z) in space.");
    cmd.addSwitch("borders", "For grid use also borders", borders);
    cmd.addInt("x", "deg_x", "degree in x direction", deg_x);
    cmd.addInt("y", "deg_y", "degree in y direction", deg_y);
    cmd.addInt("r", "urefine", "initial uniform refinement steps of the fitting basis", numURef);
    cmd.addInt("b", "uknots", "initial uniform refinement (number of inserted knots) of fitting basis in u knot direction", numUknots);
    cmd.addInt("c", "vknots", "initial uniform refinement (number of inserted knots) of fitting basis in v knot direction", numVknots);
    cmd.addString("d", "data", "Input sample data", fn);
    cmd.addReal("u", "u_min", "u_min in uv box", u_min);
    cmd.addReal("v", "v_min", "v_min in uv box", v_min);
    cmd.addReal("w", "u_max", "u_max in uv box", u_max);
    cmd.addReal("z", "v_max", "v_max in uv box", v_max);
    cmd.addInt("n", "n_ugrid", "number points in ugrid", n_ugrid);
    cmd.addInt("m", "n_vgrid", "number points in vgrid", n_vgrid);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (deg_x < 1)
    { gsInfo << "Degree x must be positive.\n";  return 0;}
    if (deg_y < 1)
    { gsInfo << "Degree y must be positive.\n"; return 0;}

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
    gsFileData<> fl;
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

    // Create a tensor-basis and apply initial uniform refinement
    gsTensorBSplineBasis<2> T_tbasis( u_knots, v_knots );
    T_tbasis.uniformRefine( (1<<numURef)-1 );

    // Create fast fitting object
    gsFastFitting<real_t> ref( uv, xyz, T_tbasis);

    // Create grid which is used to compute solution:
    gsMatrix<real_t> ugrid(1,n_ugrid);
    gsMatrix<real_t> vgrid(1,n_vgrid);
    for (index_t i=0; i<n_ugrid; i++)
    {
        if (borders)
            ugrid(0,i)=i/(n_ugrid-1.0);                     // Grid including borders
        else
            ugrid(0,i)= i/(n_ugrid*1.0) + 1/(2.0*n_ugrid);        // Grid without borders
    }
    for (index_t i=0; i<n_vgrid; i++)
    {
        if (borders)
            vgrid(0,i)=i/(n_vgrid-1.0);                     // Grid including borders
        else
            vgrid(0,i)= i/(n_vgrid*1.0) + 1/(2.0*n_vgrid);        // Grid without borders
    }

    // Print settings summary
    gsInfo<<"--------------------------------\n";
    gsInfo<<"Fast Fitting vs. Standard Fitting \n";
    gsInfo<<"--------------------------------\n";
    gsInfo<<"D: Number of samples       : "<< xyz.cols() <<".\n";
    gsInfo<<"Initial uniform refinement : "<< numURef<<".\n";
    gsInfo<<"Degree u/v direction       : "<< deg_x << " " << deg_y <<".\n";
    gsInfo<<"DoFs in u/v direction      : "<< T_tbasis.size(0) <<  " " << T_tbasis.size(1) <<".\n";
    gsInfo<<"Grid with/without borders  : "<< n_ugrid << ", " << n_vgrid << ". " << borders << ".\n";
    gsInfo<<"Delta (should be <1)       : "<< T_tbasis.size(0)*T_tbasis.size(1)/(1.0*xyz.cols()) <<".\n";
    gsInfo<<"delta (should be <1)       : "<< (1.0*xyz.cols())/(n_ugrid*n_vgrid)                 <<".\n";
    gsInfo<<"--------------------------------\n";

    gsStopwatch time;
    time.restart();
    ref.compute(ugrid, vgrid, false);
    time.stop();
    gsInfo<<"Fitting time                      : "<< time <<"\n";
    gsWriteParaview(*(ref.result()),"resultfast", 50000, false, true);

    gsFitting<real_t> ref2(uv, xyz, T_tbasis);
    time.restart();
    ref2.compute();
    time.stop();
    gsInfo<<"Fitting time                      : "<< time <<"\n";
    ref2.computeErrors();
    gsWriteParaview(*(ref2.result()),"resultslow", 50000, false, true);

    ref.plotErrors("errorplot_gk_average");

    gsInfo<<"--------------------------------\n";
    gsInfo<<"Max error fast fitting gk-average : "<< ref.maxPointError() <<".\n";
    ref.computeAllProjectedErrors(uv,xyz,ugrid,vgrid);

    ref.plotErrors("errorplot_gk_all");

    gsInfo<<"Max error fast fitting gk-all     : "<< ref.maxPointError() <<".\n";
    ref.changeParam(uv,xyz);
    ref.computeErrors();

    ref.plotErrors("errorplot_uk");

    gsInfo<<"Max error fast fitting uk         : "<< ref.maxPointError() <<".\n";
    gsInfo<<"Max error slow fitting            : "<< ref2.maxPointError() <<".\n";
    return 0;
}
