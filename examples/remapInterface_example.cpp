/** @file remapInterface_example.cpp

    @brief Some tests for gsRemapInterface

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <gismo.h>

#include <gsAssembler/gsRemapInterface.h>

using namespace gismo;


void showCorners( const gsGeometry<>& geo )
{
    gsMatrix<> in(4,2);
    in << 0,0,    1,0,    0,1,    1,1;
    gsMatrix<> out;
    geo.eval_into(in.transpose(),out);
    gsInfo << out.transpose() << "\n";

}

void showCorners3D( const gsGeometry<>& geo )
{
    gsMatrix<> in(8,3);
    in << 0,0,0,    1,0,0,    0,1,0,    1,1,0,   0,0,1,    1,0,1,    0,1,1,    1,1,1;
    gsMatrix<> out;
    geo.eval_into(in.transpose(),out);
    gsInfo << out.transpose() << "\n";

}

int main(int argc, char* argv[])
{

    gsCmdLine cmd("remapInterface_example");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    {
        gsInfo << "************* Test 1 *************\n";
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,0,1,1).release() );
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,1,1,2).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.computeTopology();

        gsInfo << "Patch 0:\n"; showCorners(mp[0]);
        gsInfo << "Patch 1:\n"; showCorners(mp[1]);


        GISMO_ENSURE ( mp.nInterfaces() == 1, "mp.nInterfaces() == "<<mp.nInterfaces());
        const boundaryInterface &bi = *(mp.iBegin());
        gsInfo << bi.first() << "\n";
        gsInfo << bi.second() << "\n";

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi);

        gsInfo << ri << "\n";

        gsMatrix<> in(4,2);

        in <<
              0,    0,
              1,    0,
              0.5,  0,
              0.25, 0;

        gsInfo << "In:\n" << in << "\n\n";

        gsMatrix<> expected(4,2);

        expected <<
              0,    1,
              1,    1,
              0.5,  1,
              0.25, 1;

        gsInfo << "Expected:\n" << expected << "\n\n";

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        gsInfo << "Out:\n" << out.transpose() << "\n\n";

        GISMO_ENSURE ( (expected - out.transpose()).norm() < 1.e-4, "");
    }
    {
        gsInfo << "************* Test 2 *************\n";
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,0,1,1).release() );
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,0,1,1,90).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.computeTopology();

        gsInfo << "Patch 0:\n"; showCorners(mp[0]);
        gsInfo << "Patch 1:\n"; showCorners(mp[1]);


        GISMO_ENSURE ( mp.nInterfaces() == 1, "mp.nInterfaces() == "<<mp.nInterfaces());
        const boundaryInterface &bi = *(mp.iBegin());
        gsInfo << bi.first() << "\n";
        gsInfo << bi.second() << "\n";

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi);

        gsInfo << ri << "\n";

        gsMatrix<> in(4,2);

        in <<
              0,    0,
              1,    0,
              0.5,  0,
              0.25, 0;


        gsInfo << "In:\n" << in << "\n\n";

        gsMatrix<> expected(4,2);

        expected <<
              0,    0,
              0,    1,
              0,    0.5,
              0,    0.25;

        gsInfo << "Expected:\n" << expected << "\n\n";



        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        gsInfo << "Out:\n" << out.transpose() << "\n\n";

        GISMO_ENSURE ( (expected - out.transpose()).norm() < 1.e-4, "");
    }
    {
        gsInfo << "************* Test 3 *************\n";
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,-1,1,0).release() );
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,0,1,1,180).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.computeTopology();

        gsInfo << "Patch 0:\n"; showCorners(mp[0]);
        gsInfo << "Patch 1:\n"; showCorners(mp[1]);


        GISMO_ENSURE ( mp.nInterfaces() == 1, "mp.nInterfaces() == "<<mp.nInterfaces());
        const boundaryInterface &bi = *(mp.iBegin());
        gsInfo << bi.first() << "\n";
        gsInfo << bi.second() << "\n";

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi);

        gsInfo << ri << "\n";

        gsMatrix<> in(4,2);

        in <<
              0,    0,
              0,    1,
              0,    0.5,
              0,    0.25;


        gsInfo << "In:\n" << in << "\n\n";

        gsMatrix<> expected(4,2);

        expected <<
              0,    1,
              0,    0,
              0,    0.5,
              0,    0.75;

        gsInfo << "Expected:\n" << expected << "\n\n";



        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        gsInfo << "Out:\n" << out.transpose() << "\n\n";

        GISMO_ENSURE ( (expected - out.transpose()).norm() < 1.e-4, "");
    }

    {
        gsInfo << "************* Test 4 *************\n";
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::lift3D(*gsNurbsCreator<>::BSplineRectangle(0,0,1,1)).release() );
        pc.push_back( gsNurbsCreator<>::lift3D(*gsNurbsCreator<>::BSplineRectangle(0,1,1,2)).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.computeTopology();

        gsInfo << "Patch 0:\n"; showCorners3D(mp[0]);
        gsInfo << "Patch 1:\n"; showCorners3D(mp[1]);


        GISMO_ENSURE ( mp.nInterfaces() == 1, "mp.nInterfaces() == "<<mp.nInterfaces());
        const boundaryInterface &bi = *(mp.iBegin());
        gsInfo << bi.first() << "\n";
        gsInfo << bi.second() << "\n";

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi);

        gsInfo << ri << "\n";

        gsMatrix<> in(8,3);

        in <<
              0,    0, 0.3,
              1,    0, 0.3,
              0.5,  0, 0.3,
              0.25, 0, 0.3,
              0,    0, 1,
              1,    0, 1,
              0.5,  0, 1,
              0.25, 0, 1;

        gsInfo << "In:\n" << in << "\n\n";

        gsMatrix<> expected(8,3);

        expected <<
              0,    1, 0.3,
              1,    1, 0.3,
              0.5,  1, 0.3,
              0.25, 1, 0.3,
              0,    1, 1,
              1,    1, 1,
              0.5,  1, 1,
              0.25, 1, 1;

        gsInfo << "Expected:\n" << expected << "\n\n";

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        gsInfo << "Out:\n" << out.transpose() << "\n\n";

        GISMO_ENSURE ( (expected - out.transpose()).norm() < 1.e-4, "");
    }
    {
        gsInfo << "************* Test 5 *************\n";
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::lift3D(*gsNurbsCreator<>::BSplineRectangle(0,0,1,1)).release() );
        pc.push_back( gsNurbsCreator<>::lift3D(*gsNurbsCreator<>::BSplineRectangle(0,0,1,1,90)).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.computeTopology();

        gsInfo << "Patch 0:\n"; showCorners3D(mp[0]);
        gsInfo << "Patch 1:\n"; showCorners3D(mp[1]);


        GISMO_ENSURE ( mp.nInterfaces() == 1, "mp.nInterfaces() == "<<mp.nInterfaces());
        const boundaryInterface &bi = *(mp.iBegin());
        gsInfo << bi.first() << "\n";
        gsInfo << bi.second() << "\n";

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi);

        gsInfo << ri << "\n";

        gsMatrix<> in(8,3);

        in <<
              0,    0, 0.3,
              1,    0, 0.3,
              0.5,  0, 0.3,
              0.25, 0, 0.3,
              0,    0, 1,
              1,    0, 1,
              0.5,  0, 1,
              0.25, 0, 1;


        gsInfo << "In:\n" << in << "\n\n";

        gsMatrix<> expected(8,3);

        expected <<
              0,    0,   0.3,
              0,    1,   0.3,
              0,    0.5, 0.3,
              0,    0.25,0.3,
              0,    0,   1,
              0,    1,   1,
              0,    0.5, 1,
              0,    0.25,1;

        gsInfo << "Expected:\n" << expected << "\n\n";



        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        gsInfo << "Out:\n" << out.transpose() << "\n\n";

        GISMO_ENSURE ( (expected - out.transpose()).norm() < 1.e-4, "");
    }
    {
        gsInfo << "************* Test 6 *************\n";
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::lift3D(*gsNurbsCreator<>::BSplineRectangle(0,-1,1,0)).release() );
        pc.push_back( gsNurbsCreator<>::lift3D(*gsNurbsCreator<>::BSplineRectangle(0,0,1,1,180)).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.computeTopology();

        gsInfo << "Patch 0:\n"; showCorners3D(mp[0]);
        gsInfo << "Patch 1:\n"; showCorners3D(mp[1]);


        GISMO_ENSURE ( mp.nInterfaces() == 1, "mp.nInterfaces() == "<<mp.nInterfaces());
        const boundaryInterface &bi = *(mp.iBegin());
        gsInfo << bi.first() << "\n";
        gsInfo << bi.second() << "\n";

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi);

        gsInfo << ri << "\n";

        gsMatrix<> in(8,3);

        in <<
              0,    0,    0.3,
              0,    1,    0.3,
              0,    0.5,  0.3,
              0,    0.25, 0.3,
              0,    0,    1,
              0,    1,    1,
              0,    0.5,  1,
              0,    0.25, 1;


        gsInfo << "In:\n" << in << "\n\n";

        gsMatrix<> expected(8,3);

        expected <<
              0,    1,    0.3,
              0,    0,    0.3,
              0,    0.5,  0.3,
              0,    0.75, 0.3,
              0,    1,    1,
              0,    0,    1,
              0,    0.5,  1,
              0,    0.75, 1;

        gsInfo << "Expected:\n" << expected << "\n\n";



        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        gsInfo << "Out:\n" << out.transpose() << "\n\n";

        GISMO_ENSURE ( (expected - out.transpose()).norm() < 1.e-4, "");
    }
    {
        gsInfo << "************* Test 7 *************\n";
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,0,2,1).release() );
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0.5,1,1.5,2).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.addInterface( 1, boxSide(boundary::south), 0, boxSide(boundary::north));

        gsInfo << "Patch 0:\n"; showCorners(mp[0]);
        gsInfo << "Patch 1:\n"; showCorners(mp[1]);


        GISMO_ENSURE ( mp.nInterfaces() == 1, "mp.nInterfaces() == "<<mp.nInterfaces());
        const boundaryInterface &bi = *(mp.iBegin());
        gsInfo << bi.first() << "\n";
        gsInfo << bi.second() << "\n";

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi);

        gsInfo << ri << "\n";

        gsMatrix<> in(3,2);

        in <<
              0.0,  0,
              0.5,  0,
              1.0,  0;

        gsInfo << "In:\n" << in << "\n\n";

        gsMatrix<> expected(3,2);

        expected <<
              0.25,  1,
              0.5,   1,
              0.75,  1;

        gsInfo << "Expected:\n" << expected << "\n\n";

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        gsInfo << "Out:\n" << out.transpose() << "\n\n";

        GISMO_ENSURE ( (expected - out.transpose()).norm() < 1.e-4, "");
    }
    {
        gsInfo << "************* Test 8 *************\n";
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,0,2,1).release() );
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0.5,1,1.5,2).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.addInterface( 0, boxSide(boundary::north), 1, boxSide(boundary::south));

        gsInfo << "Patch 0:\n"; showCorners(mp[0]);
        gsInfo << "Patch 1:\n"; showCorners(mp[1]);


        GISMO_ENSURE ( mp.nInterfaces() == 1, "mp.nInterfaces() == "<<mp.nInterfaces());
        const boundaryInterface &bi = *(mp.iBegin());
        gsInfo << bi.first() << "\n";
        gsInfo << bi.second() << "\n";

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi);

        gsInfo << ri << "\n";

        gsMatrix<> in(3,2);

        in <<
              0.25,  1,
              0.5,   1,
              0.75,  1;

        gsInfo << "In:\n" << in << "\n\n";

        gsMatrix<> expected(3,2);

        expected <<
              0.0,  0,
              0.5,  0,
              1.0,  0;

        gsInfo << "Expected:\n" << expected << "\n\n";

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        gsInfo << "Out:\n" << out.transpose() << "\n\n";

        GISMO_ENSURE ( (expected - out.transpose()).norm() < 1.e-4, "");
    }
    {
        gsInfo << "************* Test 9 *************\n";
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,0,1,1).release() );
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0.5,1,1.5,2).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.addInterface( 0, boxSide(boundary::north), 1, boxSide(boundary::south));

        gsInfo << "Patch 0:\n"; showCorners(mp[0]);
        gsInfo << "Patch 1:\n"; showCorners(mp[1]);


        GISMO_ENSURE ( mp.nInterfaces() == 1, "mp.nInterfaces() == "<<mp.nInterfaces());
        const boundaryInterface &bi = *(mp.iBegin());
        gsInfo << bi.first() << "\n";
        gsInfo << bi.second() << "\n";

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi);

        gsInfo << ri << "\n";

        gsMatrix<> in(3,2);

        in <<
              0.5,  1,
              0.6,  1,
              1.0,  1;

        gsInfo << "In:\n" << in << "\n\n";

        gsMatrix<> expected(3,2);

        expected <<
              0.0,  0,
              0.1,  0,
              0.5,  0;

        gsInfo << "Expected:\n" << expected << "\n\n";

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        gsInfo << "Out:\n" << out.transpose() << "\n\n";

        GISMO_ENSURE ( (expected - out.transpose()).norm() < 1.e-4, "");
    }
    {
        gsInfo << "************* Test 10 *************\n";
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::BSplineCube(1, 0,0,0).release() );
        pc.push_back( gsNurbsCreator<>::BSplineCube(1, 1,0.3,0.5 ).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.addInterface( 0, boxSide(boundary::east), 1, boxSide(boundary::west));

        gsInfo << "Patch 0:\n"; showCorners3D(mp[0]);
        gsInfo << "Patch 1:\n"; showCorners3D(mp[1]);


        GISMO_ENSURE ( mp.nInterfaces() == 1, "mp.nInterfaces() == "<<mp.nInterfaces());
        const boundaryInterface &bi = *(mp.iBegin());
        gsInfo << bi.first() << "\n";
        gsInfo << bi.second() << "\n";

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi);

        gsInfo << ri << "\n";

        gsMatrix<> in(5,3);

        in <<
              1, 0.3, 0.5,
              1, 0.6, 0.5,
              1, 0.6, 1.0,
              1, 1.0, 1.0,
              1, 1.0, 0.5;

        gsInfo << "In:\n" << in << "\n\n";

        gsMatrix<> expected(5,3);

        expected <<
              0, 0.0, 0.0,
              0, 0.3, 0.0,
              0, 0.3, 0.5,
              0, 0.7, 0.5,
              0, 0.7, 0.0;

        gsInfo << "Expected:\n" << expected << "\n\n";

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        gsInfo << "Out:\n" << out.transpose() << "\n\n";

        GISMO_ENSURE ( (expected - out.transpose()).norm() < 1.e-4, "");
    }

    return 0;
}
