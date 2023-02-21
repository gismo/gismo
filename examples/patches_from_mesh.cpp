/** @file refit_patches.cpp

    @brief Computes patches from structured (tensor-product) data samples by fitting.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    std::string filename;
    std::string out;
    bool plot = false;

    gsCmdLine cmd("Computes patches from structured (tensor-product) data samples by fitting. Give a file path to an XML or 3dm file to refit the patches!");
    cmd.addPlainString("filename", "File containing multipatch input (.xml).", filename);
    cmd.addSwitch("plot", "plot results", plot);
    cmd.addString("w","write", "write results", out);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(!filename.empty(),"Filename is empty");
    gsStopwatch time;
    gsFileData<> fd(filename);

    gsSurfMesh HEmesh;
    fd.getFirst<gsSurfMesh>(HEmesh);

    gsInfo<<"Reading mesh:\t"<<time.stop()<<" seconds\n";
    if (plot) gsWriteParaview(HEmesh,"HEmesh");

    time.restart();

    gsMultiPatch<> mp;

    // Counts for each vertex the number of passes
    auto vpassed = HEmesh.vertex_property<gsSurfMesh::Scalar>("v:passed", 0 );
    // Index of the curve loop (negative) or of the patch (positive)
    auto hindex = HEmesh.halfedge_property<int>("h:index", 0 );
    // Patch index of each face
    auto findex = HEmesh.face_property<int>("f:index", -1 );

    // Create a stack of EVs and boundary EVs
    std::list<gsSurfMesh::Vertex> EVs;
    for (auto v : HEmesh.vertices())
    {
        if (HEmesh.valence(v) != 4 && !HEmesh.is_boundary(v))//interior
            EVs.push_back(v);
        else if (HEmesh.valence(v) > 3 && HEmesh.is_boundary(v)) //boundary
            EVs.push_back(v);
        // else if (HEmesh.valence(v) == 2 && HEmesh.is_boundary(v)) //boundary corner
        //     EVs.push_back(v);
    }

    // For all EVs, find the curve-loop over all adjacent edges
    index_t curveloopIdx = -1;
    for (auto EV = EVs.begin(); EV!=EVs.end(); EV++)
    {
        // std::list<gsSurfMesh::Vertex>::iterator EV = EVs.begin();
        for ( auto he : HEmesh.halfedges(*EV) ) // iterate over all HE that come from a EV
        {
            auto h = he;
            // Mark the vertex from which we departed as passed
            auto v = HEmesh.from_vertex(h);
            auto vold = HEmesh.from_vertex(h);
            vpassed[v] += 1;

            // stopping conditions:
            // - boundary is hit
            // - HE is already assigned to another loop
            // - hit another V that has been crossed,
            while (true)
            {
                v = HEmesh.to_vertex(h);
                vold = HEmesh.from_vertex(h);
                // If h is already assigned to a curve loop, we stop
                if (hindex[h]!=0)
                    break;

                hindex[h] =
                hindex[HEmesh.opposite_halfedge(h)] = curveloopIdx;
                if (!HEmesh.is_boundary(v))
                {
                    if (HEmesh.valence(v)==4) // interior ordinary vertex
                    {
                        h = HEmesh.next_halfedge(h);
                        h = HEmesh.opposite_halfedge(h);
                        h = HEmesh.next_halfedge(h);
                    }
                    else
                        break; // EV, thus stop this curveloop
                }
                else // is boundary
                {
                    if (HEmesh.valence(v) > 3)
                        break; // EV, thus stop this curveloop
                }
                // Mark the to-vertex as passed, if it is not an EV
                vpassed[v] += 1;

                // Check if the to-vertex is a boundary vertex. If yes, we stop
                if (HEmesh.is_boundary(v))
                    break;
            }
            curveloopIdx--;
        }
    }

    // gsWriteParaview(HEmesh,"HEmesh",{"v:passed"});
    // Collect intersections and EVs (points that have been passed more than once)
    EVs.clear();
    EVs.resize(0);
    for (auto v : HEmesh.vertices())
        if (vpassed[v]>1) // Intersection or EV
            EVs.push_back(v);
        // Probably not needed:
        // else if (vpassed[v]==1 && HEmesh.is_boundary(v)) // Point that splits the boundary
        //     EVs.push_back(v);

    // Patch index which will be assigned to half-edges, starting from 1 now, because 0 is reserved for non-assigned half edges
    index_t patchIdx = 1;
    // Stores the number of element in both directions for each patch
    std::vector<std::pair<index_t,index_t>> dirSizes;
    // From each starting point, we take each half edge again and we assign patch indices to the incident faces
    for (auto EV = EVs.begin(); EV!=EVs.end(); EV++)
    {
        for ( auto he : HEmesh.halfedges(*EV) ) // iterate over all HE that come from a EV
        {
            if (!HEmesh.is_valid(he))
                continue;
            // Check if the half-edge is already assigned to a patch
            if (hindex[he] > 0)
                continue;

            auto h = he;
            if (HEmesh.is_boundary(h))
                continue;
            // Store other patch direction sizes for a check
            std::pair<index_t,index_t> dirSize(1,1);
            bool dir = 0;
            // Stopping conditions:
            // 1. Startpoint reached
            // 2. The half-edge h has a positive index, meaning it is assigned to a patch.
            while (true)
            {
                if (hindex[h] > 0)
                {
                    gsWarn<<"Half-edge has already been passed\n";
                    break;
                }
                // ---- Assign the patch index to the half-edge
                hindex[h] = patchIdx;
                // ---- Count on the direction
                if (dir)
                    dirSize.second++;
                else
                    dirSize.first++;

                // If h goes to the original EV, we stop
                if (HEmesh.to_vertex(h)==*EV)
                    break;

                // Move on (goes around the corner on the same face)
                h = HEmesh.next_halfedge(h);

                // Check if the next half-edge has a negative index OR is a boundary OR if the opposite has a face with patchIdx.
                // If so, it is a patch boundary and we ended up in a corner
                // ---- Go around the corner
                if (hindex[h]<0
                    || ( ( HEmesh.is_boundary(h) || HEmesh.is_boundary(HEmesh.opposite_halfedge(h)) ) && HEmesh.is_boundary(HEmesh.to_vertex(h)) )
                    )
                {
                    // if dir was 1, pushback dirSize and empty temporary one
                    if (dir==1)
                    {
                        dirSizes.push_back(dirSize);
                        dirSize.first = 1;
                        dirSize.second= 1;
                    }
                    dir = !dir;
                    continue;
                }
                // ---- Stopping condition 2
                else if (hindex[h] > 0)
                {
                    gsWarn<<"Something went wrong. Half-edge is already assigned\n";
                    break;
                }
                // check if the next half-edge has an index of 0. If so, it is an interior half-edge
                // ---- Go straight
                else if (hindex[h]==0)
                {
                    h = HEmesh.opposite_halfedge(h);
                    h = HEmesh.next_halfedge(h);
                    continue;
                }
                else
                    gsWarn<<"Something went wrong?\n";
            }

            // Make a linear tensor basis of size dir.first,dir.second
            gsKnotVector<real_t> KV0(0,dirSize.first -1,dirSize.first -2,2,1,1);
            gsKnotVector<real_t> KV1(0,dirSize.second-1,dirSize.second-2,2,1,1);
            gsTensorBSplineBasis<2,real_t> tbasis(KV0,KV1);
            gsMatrix<> coefs(dirSize.first*dirSize.second,3);
            coefs.setZero();

            // Start from the EV again. We start in dirsize
            h = he;
            if (HEmesh.is_boundary(h))
                continue;
            if (!HEmesh.is_valid(he))
                continue;

            auto h0 = h;
            if (HEmesh.is_boundary(h0))
                gsWarn<<"Edge is already a boundary...\n";
            bool flip = false;
            gsMatrix<> tmpcoefs(dirSize.first,3);
            for (index_t i=0; i!=dirSize.second; i++)
            {
                for (index_t j=0; j!=dirSize.first-2; j++) //
                {
                    tmpcoefs.row(j) = HEmesh.position(HEmesh.from_vertex(h)).transpose();

                    h = HEmesh.next_halfedge(h);
                    if (!HEmesh.is_boundary(h))
                    {
                        h = HEmesh.opposite_halfedge(h);
                        h = HEmesh.next_halfedge(h);
                    }
                }
                tmpcoefs.row(dirSize.first-2) = HEmesh.position(HEmesh.from_vertex(h)).transpose();
                tmpcoefs.row(dirSize.first-1) = HEmesh.position(HEmesh.to_vertex(h)).transpose();

                if (flip)
                    h = HEmesh.opposite_halfedge(h);

                h = HEmesh.next_halfedge(h);
                h = HEmesh.next_halfedge(h);

                if (flip)
                {
                    h = HEmesh.opposite_halfedge(h);
                    tmpcoefs = tmpcoefs.colwise().reverse().eval();
                }
                coefs.block(i*dirSize.first,0,dirSize.first,3) = tmpcoefs;

                flip = !flip;
            }
            mp.addPatch(tbasis.makeGeometry(give(coefs)));
            patchIdx++;
        }
    }

    gsInfo<<"Making multipatch:\t"<<time.stop()<<" seconds\n";
    time.restart();
    if (plot)
    {
        gsWriteParaview(mp,"mp",1000,true);
        gsInfo<<"Plotting multipatch:\t"<<time.stop()<<" seconds\n";
    }
    time.restart();
    if (!out.empty())
    {
        gsWrite(mp,out);
        gsInfo<<"Writing multipatch:\t"<<time.stop()<<" seconds\n";
    }

    for (index_t k = 0; k!=dirSizes.size(); k++)
        gsDebug<<"Patch "<<k<<": Dir 0: "<<dirSizes[k].first<<"; Dir 1: "<<dirSizes[k].second<<"\n";


    return EXIT_SUCCESS;
}
