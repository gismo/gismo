/** @file geometry_example.cpp

    @brief Tutorial on gsGeometry abstract class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <iostream>
#include <gismo.h>
#include <gsCore/gsGeometry.h>

using namespace gismo;


// Returns the string with the size of a matrix.
template <typename T>
std::string size(const gsMatrix<T>& matrix)
{
    std::string result = "(" + util::to_string(matrix.rows()) + " x " +
        util::to_string(matrix.cols()) + ")";

    return result;
}


template <typename T>
void rotate2D(gsGeometry<T> & geom, const T turndeg = 0.0, const T Tx = 0.0, const T Ty = 0.0)
{
    const T pi = 3.1415926535897932384626433832795;
    T r = turndeg / 180 * pi;

    T tx, ty;
    for(index_t i =0; i < geom.coefs().rows(); i++)
    {
        tx = geom.coefs()(i,0);
        ty = geom.coefs()(i,1);
        geom.coefs()(i,0) = math::cos(r) * (tx-Tx) - math::sin(r) * (ty-Ty) + Tx;
        geom.coefs()(i,1) = math::sin(r) * (tx-Tx) + math::cos(r) * (ty-Ty) + Ty;
    }
}

template <typename T>
void shift2D(gsGeometry<T> & geom, const T dx = 0.0, const T dy = 0.0)
{
    geom.coefs().col(0) += gsVector<T>::Ones(geom.coefs().rows())*dx;
    geom.coefs().col(1) += gsVector<T>::Ones(geom.coefs().rows())*dy;
}

template <typename T>
void shift2D(gsMultiPatch<T> & mp, const T dx = 0.0, const T dy = 0.0)
{
    for (size_t k = 0; k!=mp.nPatches(); k++)
    {
        mp.patch(k).coefs().col(0) += gsVector<T>::Ones(mp.patch(k).coefs().rows())*dx;
        mp.patch(k).coefs().col(1) += gsVector<T>::Ones(mp.patch(k).coefs().rows())*dy;
    }
}

template <typename T>
void mirror(gsGeometry<T> & geom, bool axis)
{
    T mid = geom.coefs().col(!axis).maxCoeff() - geom.coefs().col(!axis).minCoeff();
    geom.coefs().col(!axis) -= gsVector<T>::Ones(geom.coefs().rows())*mid;
    geom.coefs().col(!axis) *= -1;

}

template <typename T>
void mirror(gsMultiPatch<T> & mp, bool axis)
{
    gsMatrix<T> bbox;
    mp.boundingBox(bbox);

    T mid = (bbox(!axis,1)+bbox(!axis,0))/2;
    for (size_t p = 0; p!=mp.nPatches(); p++)
    {
        mp.patch(p).coefs().col(!axis) -= gsVector<T>::Ones(mp.patch(p).coefs().rows())*mid;
        mp.patch(p).coefs().col(!axis) *= -1;
    }
}

template <typename T>
void scale2D(gsGeometry<T> & geom, T factor = 1.0)
{

    // gsMatrix<T> bbox;
    // mp.boundingBox(bbox);

    // T mid = (bbox(!axis,1)+bbox(!axis,0))/2;
    for (index_t k = 0; k!= geom.coefs().cols(); k++)
    {
        geom.coefs().col(k) *= factor;
    }
}

template <typename T>
void scale2D(gsGeometry<T> & geom, std::vector<T> factors)
{
    GISMO_ENSURE(factors.size()==geom.coefs().cols(),"Number of scaling factors must be the same as the number of dimensions");
    // gsMatrix<T> bbox;
    // mp.boundingBox(bbox);

    // T mid = (bbox(!axis,1)+bbox(!axis,0))/2;
    for (index_t k = 0; k!= geom.coefs().cols(); k++)
    {
        geom.coefs().col(k) *= factors.at(k);
    }
}

template <typename T>
void scale2D(gsMultiPatch<T> & mp,  T factor = 1.0)
{

    // gsMatrix<T> bbox;
    // mp.boundingBox(bbox);

    // T mid = (bbox(!axis,1)+bbox(!axis,0))/2;
    for (size_t p = 0; p!=mp.nPatches(); p++)
        scale2D(mp.patch(p),factor);
}

template <typename T>
void scale2D(gsMultiPatch<T> & mp, std::vector<T> factors)
{

    // gsMatrix<T> bbox;
    // mp.boundingBox(bbox);

    // T mid = (bbox(!axis,1)+bbox(!axis,0))/2;
    for (size_t p = 0; p!=mp.nPatches(); p++)
        scale2D(mp.patch(p),factors);
}

template <typename T>
void makeGrid(gsMultiPatch<T> & mp, const index_t M=0, const index_t N=0)
{
    gsMultiPatch<T> mp_ori(mp);
    gsMatrix<T> bbox;
    mp.boundingBox(bbox);

    T L = bbox(0,1)-bbox(0,0);
    T H = bbox(1,1)-bbox(1,0);

    mp.clear();
    for (index_t m = 0; m!=M; m++)
        for (index_t n = 0; n!=N; n++)
        {
            gsMultiPatch<T> mp_tmp(mp_ori);
            shift2D(mp_tmp,(m)*L,(n)*H);
            for (size_t k=0; k!=mp_tmp.nPatches(); k++)
                mp.addPatch(mp_tmp.patch(k));
        }
    mp.computeTopology();
}

template <typename T>
gsMultiPatch<T> makeGrid(std::vector<gsMultiPatch<T>> & mps, const index_t M=0, const index_t N=0)
{
    gsMultiPatch<T> mp;

    std::vector<gsMultiPatch<T>> mps_ori(mps);
    std::vector<T> Hs,Ls;
    Hs.reserve(mps.size());
    Ls.reserve(mps.size());
    gsMatrix<T> bbox;
    for (typename std::vector<gsMultiPatch<T>>::iterator it = mps.begin(); it!=mps.end(); it++)
    {
        it->boundingBox(bbox);
        Ls.push_back(bbox(0,1)-bbox(0,0));
        Hs.push_back(bbox(1,1)-bbox(1,0));
    }

    typename std::vector<gsMultiPatch<T>>::iterator mp_it = mps.begin();
    typename std::vector<T>::iterator L_it = Ls.begin();
    typename std::vector<T>::iterator H_it = Hs.begin();
    for (index_t m = 0; m!=M; m++)
    {
        for (index_t n = 0; n!=N; n++)
        {
            gsMultiPatch<T> mp_tmp(*mp_it);
            shift2D(mp_tmp,(m)*(*L_it),(n)*(*H_it));
            for (size_t k=0; k!=mp_tmp.nPatches(); k++)
                mp.addPatch(mp_tmp.patch(k));

            mp_it==mps.end()-1 ? mp_it = mps.begin() : mp_it++;
            L_it==Ls.end()-1 ? L_it = Ls.begin() : L_it++;
            H_it==Hs.end()-1 ? H_it = Hs.begin() : H_it++;
        }

        // if
        if (mps.size() % 2 == 0 && N % 2 == 0)
        {
            mp_it==mps.end()-1 ? mp_it = mps.begin() : mp_it++;
            L_it==Ls.end()-1 ? L_it = Ls.begin() : L_it++;
            H_it==Hs.end()-1 ? H_it = Hs.begin() : H_it++;
        }
    }
    mp.computeTopology();
    return mp;
}

// template <typename T>
// gsMultiPatch<T> makeGrid(gsGeometry<T> & g, const index_t m = 0.0, const index_t n = 0.0)
// {
//     return makeGrid(g.clone(),m,n);
// }


int main(int argc, char* argv[])
{

    std::string input("surfaces/simple.xml");
    std::string output("");

    gsCmdLine cmd("Tutorial on gsGeometry class.");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addString("o", "output", "Name of the output file", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======================================================================
    // reading the geometry
    // ======================================================================

    gsFileData<> fileData(input);



    // ======================================================================
    // writing to paraview
    // ======================================================================

    if (!output.empty())
    {
        gsMultiPatch<> mp, mp_copy;
        std::vector<gsMultiPatch<real_t>> container(2);

        // mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,1,1,45));
        // gsWriteParaview(mp,output + "_a",1000,true,true);

        // mp.clear();
        // mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,1,1,0));
        // rotate2D(mp.patch(0),45.0,0.0,0.0);
        // gsWriteParaview(mp,output + "_b",1000,true,true);


        // mp.clear();
        // mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,1,1,0));
        // mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,1,1,0));
        // shift2D(mp.patch(1),1.0,0.);
        // shift2D(mp,1.0,1.0);
        // makeGrid(mp,3,3);
        // gsWriteParaview(mp,output + "_c",1000,false);

        mp.clear();
        mp.addPatch(gsNurbsCreator<>::BSplineTrapezium(1,0.5,1.0,0.0,0.0));
        mp.addPatch(gsNurbsCreator<>::BSplineTrapezium(1,0.5,1.0,0.0,90.));
        mp.addPatch(gsNurbsCreator<>::BSplineTrapezium(1,0.5,1.0,0.0,180));
        mp.addPatch(gsNurbsCreator<>::BSplineTrapezium(1,0.5,1.0,0.0,-90));
        mp_copy = mp;

        mp.clear();
        mp = mp_copy;
        mirror(mp,1);

        container[0] = mp_copy;
        container[1] = mp;
        mp = makeGrid(container,6,4);
        mp.uniformRefine(2);
        gsWriteParaview(mp,output + "_diamond",1000,true);


        mp.clear();
        mp.addPatch(gsNurbsCreator<>::NurbsArcTrapezium(1,0.5,1.0,0.0,0.0));
        mp.addPatch(gsNurbsCreator<>::NurbsArcTrapezium(1,0.5,1.0,0.0,90.));
        mp.addPatch(gsNurbsCreator<>::NurbsArcTrapezium(1,0.5,1.0,0.0,180));
        mp.addPatch(gsNurbsCreator<>::NurbsArcTrapezium(1,0.5,1.0,0.0,-90));
        mp_copy = mp;

        mp.clear();
        mp = mp_copy;
        mirror(mp,1);

        container[0] = mp_copy;
        container[1] = mp;
        mp = makeGrid(container,6,4);
        mp.uniformRefine(2);
        gsWriteParaview(mp,output + "_ellipse",1000,true);

        // mp.clear();
        // mp.addPatch(gsNurbsCreator<>::NurbsArcTrapezium(1,0.5,1.0,0.0,0.0));
        // scale2D(mp,0.5);
        // mp.uniformRefine(2);
        // gsWriteParaview(mp,output + "_g",1000,true,false);
    }
    else
        gsInfo << "Done. No output created, re-run with --output <fn> to get a ParaView "
                  "file containing the solution.\n";




    return 0;
}


