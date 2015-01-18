/** @file gsView.cpp

    @brief Produce Paraview file output from XML input, fo Visualizing  G+Smo objects

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    std::string fn;
    unsigned np(1000);
    int choice(0);
    bool plot_mesh = false;
    bool plot_net = false;
    try 
    {
        gsCmdLine cmd("Hi, give me a file (.txt, .axl) and I will try to draw it!");    
        gsArgSwitch a5("g", "geometry", "Try to find and plot a geometry contained in the file", cmd);
        gsArgSwitch a4("m", "mesh", "Try to find and plot a mesh contained in the file", cmd);
        gsArgSwitch a3("b", "basis", "Try to find and plot a basis contained in the file", cmd);
        gsArgVal<int> a2("s","samples", 
                         "Number of samples to use for viewing", 
                         false,1000, "samples", cmd );
        gsArgSwitch arg_plot_mesh("e", "element", "Plot the element mesh (when applicable)", cmd);
        gsArgSwitch arg_plot_net("c", "controlNet", "Plot the control net (when applicable)", cmd);

        gsArgValPlain<std::string> a1("filename","File containing data to draw (.xml, .axl, .txt)", 
                                      false, "", "string",cmd );
        cmd.parse(argc,argv);
        fn = a1.getValue();
        np = a2.getValue();
        plot_mesh = arg_plot_mesh.getValue();
        plot_net  = arg_plot_net.getValue();

        if (fn.empty() )
        {
            gsInfo<< cmd.getMessage();
            gsInfo<<"\nType "<< argv[0]<< " -h, to get the list of command line options.\n";
            return 0;
        }

        if (a3.getValue() )
            choice= 3;
        else if (a4.getValue() )
            choice= 4;
        else if (a5.getValue() )
            choice= 5;
        
    } catch ( gsArgException& e )
    { gsInfo << "Error: " << e.error() << " " << e.argId() << "\n"; }
    
    gsFileData<>  filedata(fn);
    
    switch ( choice)
    {
    case 3:
    {
        gsBasis<> * bb = filedata.getAnyFirst< gsBasis<> >();
        if (bb)
            gsInfo<< "Got "<< *bb <<"\n";
        else
        {
            gsInfo<< "Did not find any basis to plot in "<<fn<<", quitting."<<"\n";
            return 0;
        }
        
        gsWriteParaview( *bb , "gsview", np, true);
        
        delete bb;

        return system("paraview gsview.pvd &");
    }
    case 4:
    {
        gsMesh<> * msh = filedata.getAnyFirst< gsMesh<> >();
        if (msh)
            gsInfo<< "Got "<< *msh <<"\n";
        else
        {
            gsInfo<< "Did not find any mesh to plot in "<<fn<<", quitting."<<"\n";
            return 0;
        }
        gsWriteParaview( *msh, "gsview");
        delete msh;

        return system("paraview gsview.pvd &");
    }
    case 5:{
        gsGeometry<> * geo = filedata.getAnyFirst< gsGeometry<> >();
        if (geo)
            gsInfo<< "Got "<< *geo <<"\n";
        else
        {
            gsInfo<< "Did not find any geometry to plot in "<<fn<<", quitting."<<"\n";
            return 0;
        }

        gsWriteParaview( *geo , "gsview", np, plot_mesh, plot_net);
        delete geo;

        return system("paraview gsview.pvd &");
    }
    default:
        if ( filedata.has< gsGeometry<> >() )
        {
            std::vector<gsGeometry<>* > geo = filedata.getAll< gsGeometry<> >();
            if ( ! geo.empty() )
                gsInfo<< "Got "<< geo.size() <<" patch"<<(geo.size() == 1 ? "." : "es.") <<"\n";
            else
            {
                gsInfo<< "Problem encountered in file "<<fn<<", quitting." <<"\n";
                return 0;
            }
            
            gsWriteParaview( geo , "gsview", np, plot_mesh, plot_net);
            
            freeAll( geo );
            
            return system("paraview gsview.pvd &");            
        }
        
        if ( filedata.has< gsMesh<> >() )
        {
            gsMesh<> * msh = filedata.getFirst< gsMesh<> >();
            if (msh)
                gsInfo<< "Got "<< *msh <<"\n";
            else
            {
                gsInfo<< "Problem encountered in file "<<fn<<", quitting." <<"\n";
                return 0;
            }
            
            gsWriteParaview( *msh, "gsview");
            delete msh;
            
        return system("paraview gsview.pvd &");
        }

    if ( filedata.has< gsBasis<> >() )
    {
        gsBasis<> * bb = filedata.getFirst< gsBasis<> >();
        //bb->uniformRefine(3);
        
        if (bb)
            gsInfo<< "Got "<< *bb <<"\n";
        else
        {
            gsInfo<< "Problem encountered in file "<<fn<<", quitting." <<"\n";
            return 0;
        }
        
        gsWriteParaview( *bb , "gsview", np, plot_mesh);
        
        delete bb;

        return system("paraview gsview.pvd &");        
    }

    if ( filedata.has< gsSolid<> >() )
    {
        gsSolid<> * bb = filedata.getFirst< gsSolid<> >();

        if (bb)
            gsInfo<< "Got "<< *bb <<"\n";
        else
        {
            gsInfo<< "Problem encountered in file "<<fn<<", quitting." <<"\n";
            return 0;
        }
        
        gsWriteParaviewSolid( *bb, "gsview", np);
        //gsWriteParaview( *bb, "gsview", np, 0, 0.02);
        
        delete bb;

        return system("paraview gsview.pvd &");        
    }

    if ( filedata.has< gsTrimSurface<> >() )
    {
        gsTrimSurface<> * bb = filedata.getFirst< gsTrimSurface<> >();
        //bb->uniformRefine(3);

        if (bb)
            gsInfo<< "Got "<< *bb <<"\n";
        else
        {
            gsInfo<< "Problem encountered in file "<<fn<<", quitting." <<"\n";
            return 0;
        }

        gsWriteParaview( *bb, "gsview", np);
        
        delete bb;

        return system("paraview gsview.pvd &");
    }

    if ( filedata.has< gsPlanarDomain<> >() )
    {
        gsPlanarDomain<> * bb = filedata.getFirst< gsPlanarDomain<> >();
        //bb->uniformRefine(3);

        if (bb)
            gsInfo<< "Got "<< *bb <<"\n";
        else
        {
            gsInfo<< "Problem encountered in file "<<fn<<", quitting." <<"\n";
            return 0;
        }
        
        gsMesh<> * msh = bb->toMesh(np);

        gsWriteParaview( *msh , "gsview");
        
        delete bb;
        delete msh;

        return system("paraview gsview.pvd &");        
    }
        
    if ( filedata.has< gsMatrix<> >() )
    {
        gsMatrix<> * bb = filedata.getFirst< gsMatrix<> >();
        //bb->uniformRefine(3);
        
        if (bb)
            gsInfo<< "Got Matrix with "<< bb->cols() <<" points.\n";
        else
        {
            gsInfo<< "Problem encountered in file "<<fn<<", quitting." <<"\n";
            return 0;
        }

        gsInfo<< "Plot "<< bb->rows() <<"D points.\n";
        if (bb->rows() == 2 )    
        {         
            gsWriteParaviewPoints<real_t>( bb->row(0),
                             bb->row(1), "gsview");
        }
        else if (bb->rows() == 3 )    
            gsWriteParaviewPoints<real_t>( bb->row(0).eval(),
                             bb->row(1).eval() , 
                             bb->row(2).eval(), "gsview");
        else
            gsInfo<< "In trouble...\n";

        delete bb;

        return system("paraview gsview.pvd &");
    }

    gsInfo<< "Did not find anything to plot in "<<fn<<", quitting."<<"\n";
    }

    return 0;
}
