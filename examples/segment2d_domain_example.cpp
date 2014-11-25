// example segmentation

#include <iostream>
#include <math.h>

#include <gismo.h>



using namespace std;
using namespace gismo;

int main(int argc, char *argv[])
{
    gsPlanarDomain<> * Pdomain = NULL;
    int n_points(10);
    double tolerance = 1e-4;
    unsigned smpl(30);
    bool plot = false;
    std::vector<real_t> start_v;
    std::string fn("");

    try
    {
        gsCmdLine cmd("Segmentation in quadrangular patches of a planar domain ");

        gsArgVal<std::string> a5("g","geometry","File containing Geometry (.axl, .txt)",
                                 false,"", "geometry file", cmd );
        gsArgSwitch a4("", "plot", "Plot result with ParaView ", cmd);
        gsArgMultiVal<real_t> a3("i","starting_point",
                                 "point inside the computational domain from which looking for the starting point of the tracing curves", 
                                 false,"fixed starting point", cmd );
        gsArgVal<int> a2("s","samples", 
                         "Number of samples", 
                         false,smpl, "samples", cmd );
        gsArgVal<int> a1("p","n_points", 
                         "Number of points traced per curve", 
                         false,n_points, "points traced", cmd );
        gsArgVal<> a0("t","tolerance", 
                            "Required accuracy", 
                            false,tolerance, "tolerance", cmd );
        
        cmd.parse(argc,argv);
        
        
        fn = a5.getValue();
        
        if (fn.empty() )
        {
            fn = GISMO_SOURCE_DIR;
            fn+="/filedata/planar/amoeba0_pdomain.xml";
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
        start_v = a3.getValue();
        if(start_v.empty())
        {
            start_v.push_back(1.5);
            start_v.push_back(2);
        }
        
        smpl = a2.getValue();
        n_points = a1.getValue();
        tolerance = a0.getValue();
    } catch ( gsArgException& e )
    { cout << "Error: " << e.error() << " " << e.argId() << endl; }
    
    
    std::cout << "------------------------------------------------------------"
                 "\nInput Arguments: \n\n"
                 "input geometry: " << fn <<
                 "\n how many samples:  "<<smpl<<
                 "\n number of points traced per curve: "<<n_points<<
                 "\n tolerance: "<<tolerance<<"\n"
                 "------------------------------------------------------------"
                 "\n\n";


    gsAsMatrix<> start (start_v,2,1);
    
    bool c= true;
    gsMesh<> res;
    
    Pdomain->segment( n_points, tolerance, c, res );
    //  cout<<"\n res \n "<<res << "\n";

    int exitCommand = 0;
    if(plot)
    {
        
        gsWriteParaview<>( *Pdomain, "Pdomain", smpl) ;
        
        gsWriteParaview<>( res, "Segmentation");

        //run: paraview        
        exitCommand = system("paraview Segmentation.vtp &");        
    }

    delete Pdomain;
    
    return exitCommand;
}
