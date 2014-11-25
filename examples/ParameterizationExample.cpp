// Example used for testing parameterization applying BEM and
//predictor corrector segmentation technique
#include <iostream>
#include <math.h>

#include <gismo.h>



using namespace std;
using namespace gismo;

int main(int argc, char *argv[])
{
    gsPlanarDomain<> * Pdomain = NULL;
    int n_points(20);
    real_t tolerance = 1e-4;
    unsigned smpl(60);
    bool plot = false;

    try
    {
        gsCmdLine cmd("Segmentation in quadrangular patches of a planar domain ");
        std::string fn;
        gsArgVal<std::string> a5("g","geometry","File containing Geometry (.axl, .txt)",
                                 false,"", "geometry file", cmd );
        gsArgSwitch a4("", "plot", "Plot result with ParaView ", cmd);
        gsArgVal<int> a2("s","samples",
                         "Number of samples",
                         false,60, "samples", cmd );
        gsArgVal<int> a1("p","n_points",
                         "Number of points traced per curve",
                         false,20, "points traced", cmd );
        gsArgVal<> a0("t","tolerance",
                            "Required accuracy",
                            false,0.0001, "tolerance", cmd );

        cmd.parse(argc,argv);


        fn = a5.getValue();

        if (fn.empty() )
        {
            fn = GISMO_SOURCE_DIR;
            fn+="/filedata/planar/planarDomainPuzzle1.xml";
        }

        // Pdomain =  gsReadFile<>( fn ) ;
        gsFileData<>  filedata(fn);

        if( filedata.has<gsPlanarDomain<> >() )
            Pdomain = filedata.getFirst< gsPlanarDomain<> >();

        else
        {

            if(filedata.has<gsCurve <> >() )
            {
                gsCurve<> * geo = filedata.getAnyFirst< gsCurve<> >();
                gsCurveLoop<> * cp = new gsCurveLoop<>( geo );
                Pdomain = new gsPlanarDomain<>( cp );
            }


            else
            {
                gsWarn<< "Did not find any planar domain or geometry in "<< fn<<", quitting.\n";
                return false;
            }
        }

        plot = a4.getValue();
        smpl = a2.getValue();
        n_points = a1.getValue();
        tolerance = a0.getValue();
    } catch ( gsArgException& e )
    { cout << "Error: " << e.error() << " " << e.argId() << endl; }



//    gsMatrix<> corners(2,4);
//    corners<< 0, 1, 0, 1,
//              0, 1, 1, 1;
   // cout<<"\n corner matrix : \n"<<corners; v
//corners.conservativeResize(2,4);
//cout<<"\n corner matrix : \n"<<corners;

    //constructing a square
//    char q;
//    gsTemplate<> square(q,1);
 bool cc = false;

 gsMesh<> par;


 Pdomain->segment( n_points, tolerance, cc, par);


    int exitCommand = 0;
    if(plot)
    {

        gsWriteParaview(*Pdomain, "quadJaka",smpl);

        gsWriteParaview<>( par, "Parameterization");

        //run: paraview
        exitCommand = system("paraview Parameterization.vtp &");
    }

    delete Pdomain;

    return exitCommand;
}


