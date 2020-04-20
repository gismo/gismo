/** @file thbRefinement_example.cpp

    @brief Demonstates THB refinement and provides info on the resulting basis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

using namespace gismo;

void refineMode(int rf, int lvl, unsigned meshSize,
                unsigned extent, std::vector<index_t> & boxes);

int main(int argc, char *argv[])
{
    index_t refLevels = 5;
    index_t refmode   = 0;
    index_t numknots  = 3;
    index_t degree    = 2;
    bool plot     = false;

    gsCmdLine cmd("Create standard refined TH>B meshes.");
    cmd.addInt("l","levels",
               "Number of refinement levels", refLevels);
    cmd.addInt("m","mode",
               "Refinement mode (0, 1, 2 3 4)", refmode);
    cmd.addInt("p","degree",
               "Spline degree", degree);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Create a tensot-producte basis
    gsKnotVector<> KV (0, 1, numknots, degree+1, 1);
    gsTensorBSplineBasis<2> tp(KV,KV);
    gsInfo<< "Coarsest level: "<< tp <<"\n";

    // Get refinement boxes based on pattern requested
    std::vector<index_t> boxes;
    refineMode(refmode, refLevels, numknots+1, degree, boxes);
    //gsInfo<< "boxes: "<< boxes.size()/5 <<"\n";

    // Construct the hierarchical basis
    //gsHBSplineBasis<2> hb(tp);
    gsTHBSplineBasis<2> thb(tp, boxes);
    gsInfo<< "THB-spline basis: "<< thb <<"\n";

    // Count the knots used
    gsInfo <<"Nested grid dimensions: ";
    unsigned kcount = 0;
    for( unsigned l = 0; l <= thb.maxLevel(); ++l)
    {
        kcount += thb.tensorLevel(l).knots(0).uSize()
                + thb.tensorLevel(l).knots(1).uSize();
        gsInfo<< thb.tensorLevel(l).size() <<", ";
    }
    gsInfo<<"\n";

    // Count the number of truncated functions
    const int numTr = thb.numTruncated();
    gsVector<int> trcount;
    trcount.setZero(thb.maxLevel()+1);
    if (numTr < 1000)
        gsInfo <<"\nCoefficient count for each truncated function: \n";
    unsigned ccount = 0;
    typedef std::map<index_t, gsSparseVector<> >::const_iterator trIter;
    for( trIter it = thb.truncatedBegin(); it != thb.truncatedEnd(); ++it)
    {
        const int lvl = thb.levelOf(it->first);
        if (numTr < 1000)
            gsInfo << it->second.nonZeros() <<", ";
        trcount[lvl]++;
        ccount += it->second.nonZeros();
    }
    gsInfo<<"\n\n";

    // Print statistics
    gsInfo <<"Truncated functions: "<< numTr <<"\n";
    gsInfo <<"          per level: "<< trcount.transpose() <<"\n";
    gsInfo <<"Total coeffs stored: "<< ccount
         <<" ("<< (sizeof(real_t)*ccount >> 20) <<"MB)\n";
    gsInfo <<"Total knots stored : "<< kcount
         <<" ("<< ( (sizeof(real_t)+sizeof(index_t))*kcount >> 20) <<"MB)\n";

    // Output paraview plot of the basis
    if ( plot )
        gsWriteParaview( thb , "thb_refined", 1000, true);
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return 0;
}

// Provides boxes for refinement
void refineMode(int rf, int lvl, unsigned meshSize,
                unsigned extent, std::vector<index_t> & boxes)
{
    boxes.clear();
    const unsigned nb  = meshSize*(1<<lvl);
    const unsigned mid = nb / 2;
    const unsigned bpd = 5*nb;

    switch( rf )
    {
    case 0: // one corner
        boxes.resize(5);
        boxes[0] = lvl;
        boxes[1] = boxes[2] = 0;
        boxes[3] = boxes[4] = 1;
        gsInfo << "Corner refinement.\n";
        break;
    case 1: // diagonal
        boxes.resize(5*nb);
        for ( unsigned i = 0; i!= nb; ++i)
        {
            boxes[5*i  ] = lvl;
            boxes[5*i+1] = boxes[5*i+2] = i-extent/2;
            boxes[5*i+3] = boxes[5*i+4] = i+extent/2;
        }
        gsInfo << "Diagonal refinement.\n";
        break;
    case 2: // cross
        boxes.resize(2 *  bpd);
        for ( unsigned i = 0; i!= nb; ++i)
        {
            for ( unsigned k = 0; k!= 2; ++k)
            {
                boxes[k*bpd+5*i  ] = lvl;
                boxes[k*bpd+5*i+1] = boxes[k*bpd+5*i+2] = i-extent/2;
                boxes[k*bpd+5*i+3] = boxes[k*bpd+5*i+4] = i+extent/2;
                boxes[k*bpd+5*i+1+k] = mid-1;
                boxes[k*bpd+5*i+3+k] = mid+1;
            }
        }
        gsInfo << "Cross refinement.\n";
        break;
    case 3: // central
        boxes.resize(5);
        boxes[0] = lvl;
        boxes[1] = boxes[2] = mid - extent/2;
        boxes[3] = boxes[4] = mid + extent/2;
        gsInfo << "Central refinement.\n";
        break;
    case 4: // stripes
        boxes.resize(5*lvl);
        for ( int i = 0; i<lvl; ++i)
        {
            boxes[5*i  ] = i+1;
            boxes[5*i+1] = boxes[5*i+2] = 0 ;
            boxes[5*i+3] = ( meshSize <<(i+1) );
            boxes[5*i+4] = meshSize;
        }
        gsInfo << "Stripes refinement.\n";
        break;

    default:
        gsInfo << "No refinement.\n";
        break;
    }
}
