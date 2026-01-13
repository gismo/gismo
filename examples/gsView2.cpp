/** @file gsView2.cpp

    @brief Produce Paraview file output from XML input using gsParaview class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>
#include <gsIO/gsParaview.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    std::string pname("gsview"), fn("");
    index_t numSamples(1000);
    index_t choice(0);
    bool plot_mesh = false;
    bool plot_net = false;
    bool plot_patchid = false;
    bool get_basis = false;
    bool get_mesh = false;
    bool get_geo = false;
    bool show = true;

    //! [Parse Command line]
    gsCmdLine cmd("Hi, give me a file (eg: .xml) and I will try to draw it!");

    cmd.addSwitch("geometry", "Try to find and plot a geometry contained in the file", get_geo);
    cmd.addSwitch("mesh"    , "Try to find and plot a mesh contained in the file", get_mesh);
    cmd.addSwitch("basis"   , "Try to find and plot a basis contained in the file", get_basis);
    cmd.addInt   ("s", "samples", "Number of samples to use for viewing", numSamples);
    cmd.addSwitch("element"   , "Plot the element mesh (when applicable)", plot_mesh);
    cmd.addSwitch("controlNet", "Plot the control net (when applicable)", plot_net);
    cmd.addSwitch("pid"  , "Plot the ID of each patch and boundaries as color", plot_patchid);
    cmd.addSwitch("noshow", "Do not open Paraview after writing", show);
    cmd.addPlainString("filename", "File containing data to draw (.xml or third-party)", fn);
    cmd.addString("o", "oname", "Filename to use for the ParaView output", pname);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse Command line]

    if ( fn.empty() )
    {
        gsInfo<< cmd.getMessage();
        gsInfo<<"\nType "<< argv[0]<< " -h, to get the list of command line options.\n";
        return 0;
    }

    if (get_basis)
        choice = 3;
    else if (get_mesh)
        choice = 4;
    else if (get_geo)
        choice = 5;

    // Create gsParaview object and configure options
    gsParaview<real_t> pv;
    pv.options().setInt("numPoints", numSamples);
    pv.options().setSwitch("plotElements", plot_mesh);
    pv.options().setSwitch("plotControlNet", plot_net);
    pv.options().setSwitch("show", show);

    gsFileData<>  filedata(fn);

    switch ( choice )
    {
    case 3:
    {
        gsBasis<>::uPtr bb = filedata.getAnyFirst< gsBasis<> >();
        if (bb)
            gsInfo<< "Got "<< *bb <<"\n";
        else
        {
            gsInfo<< "Did not find any basis to plot in "<<fn<<", quitting."<<"\n";
            return 0;
        }

        pv.options().setSwitch("plotElements", true);
        pv.write(*bb, pname);
        break;
    }
    case 4:
    {
        gsMesh<>::uPtr msh = filedata.getAnyFirst< gsMesh<> >();
        if (msh)
            gsInfo<< "Got "<< *msh <<"\n";
        else
        {
            gsInfo<< "Did not find any mesh to plot in "<<fn<<", quitting."<<"\n";
            return 0;
        }
        pv.write(*msh, pname);
        break;
    }
    case 5:
    {
        gsGeometry<>::uPtr geo = filedata.getAnyFirst< gsGeometry<> >();
        if (geo)
            gsInfo<< "Got "<< *geo <<"\n";
        else
        {
            gsInfo<< "Did not find any geometry to plot in "<<fn<<", quitting."<<"\n";
            return 0;
        }

        pv.write(*geo, pname);
        break;
    }
    default:
        if ( filedata.has< gsMultiPatch<> >() )
        {
            gsMultiPatch<> mp;
            filedata.getFirst(mp);
            gsInfo<< "Got "<< mp <<"\n";

            if (plot_patchid)
            {
                gsField<> nfield = gsFieldCreator<>::patchIds(mp);
                pv.write(nfield, pname);
            }
            else
            {
                pv.write(mp, pname);
            }

            break;
        }
        if ( filedata.has< gsGeometry<> >() )
        {
            std::vector<gsGeometry<>::uPtr> geo = filedata.getAll< gsGeometry<> >();
            if ( ! geo.empty() )
                gsInfo<< "Got "<< geo.size() <<" patch"<<(geo.size() == 1 ? "." : "es.") <<"\n";
            else
            {
                gsInfo<< "Problem encountered in file "<<fn<<", quitting." <<"\n";
                return 0;
            }

            pv.write(memory::get_raw(geo), pname);
            break;
        }

        if ( filedata.has< gsSurfMesh >() )
        {
            auto msh = filedata.getFirst< gsSurfMesh >();
            if (msh)
                gsInfo<< "Got "<< *msh <<"\n";
            else
            {
                gsInfo<< "Problem encountered in file "<<fn<<", quitting." <<"\n";
                return 0;
            }

            // gsSurfMesh uses the old free function directly
            gsWriteParaview( *msh, pname);
            if (show)
                gsFileManager::open(pname+".vtk");
            return EXIT_SUCCESS;
        }

        if ( filedata.has< gsBasis<> >() )
        {
            gsBasis<>::uPtr bb = filedata.getFirst< gsBasis<> >();

            if (bb)
                gsInfo<< "Got "<< *bb <<"\n";
            else
            {
                gsInfo<< "Problem encountered in file "<<fn<<", quitting." <<"\n";
                return 0;
            }

            pv.write(*bb, pname);
            break;
        }

        if ( filedata.has< gsSolid<> >() )
        {
            gsSolid<>::uPtr bb = filedata.getFirst< gsSolid<> >();

            if (bb)
                gsInfo<< "Got "<< *bb <<"\n";
            else
            {
                gsInfo<< "Problem encountered in file "<<fn<<", quitting." <<"\n";
                return 0;
            }

            pv.write(*bb, pname);
            break;
        }

        if ( filedata.has< gsTrimSurface<> >() )
        {
            gsTrimSurface<>::uPtr bb = filedata.getFirst< gsTrimSurface<> >();

            if (bb)
                gsInfo<< "Got "<< *bb <<"\n";
            else
            {
                gsInfo<< "Problem encountered in file "<<fn<<", quitting." <<"\n";
                return 0;
            }

            pv.write(*bb, pname);
            break;
        }

        if ( filedata.has< gsPlanarDomain<> >() )
        {
            gsPlanarDomain<>::uPtr bb = filedata.getFirst< gsPlanarDomain<> >();

            if (bb)
                gsInfo<< "Got "<< *bb <<"\n";
            else
            {
                gsInfo<< "Problem encountered in file "<<fn<<", quitting." <<"\n";
                return 0;
            }

            pv.write(*bb, pname);
            break;
        }

        if ( filedata.has< gsMatrix<> >() )
        {
            gsMatrix<> bb;
            filedata.getFirst(bb);
            gsInfo<< "Got Matrix with "<< bb.cols() <<" points.\n";
            gsInfo<< "Plot "<< bb.rows() <<"D points.\n";
            pv.writePoints(bb, pname);
            break;
        }
        gsInfo<< "Did not find anything to plot in "<<fn<<", quitting."<<"\n";
        return 0;
    }

    return EXIT_SUCCESS;
}
