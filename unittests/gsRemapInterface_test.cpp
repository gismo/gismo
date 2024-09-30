/** @file gsRemapInterface_test.cpp

    @brief Some tests for gsRemapInterface

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include "gismo_unittest.h"
#include <gsAssembler/gsRemapInterface.h>


using namespace gismo;

namespace {
gsOptionList optHelper(index_t checkAffine)
{
    gsOptionList opt = gsRemapInterface<real_t>::defaultOptions();
    opt.setInt("CheckAffine", checkAffine);
    return opt;
}
}

SUITE(gsRemapInterface_test)
{
    TEST(MatchingAffineTest1)
    {
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,-1,1,0).release() );
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,0,1,1,180).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.computeTopology();

        CHECK ( mp.nInterfaces() == 1 );
        const boundaryInterface &bi = *(mp.iBegin());

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi,optHelper(5));

        CHECK ( ri.isAffine() );
        CHECK ( ri.isMatching() );

        gsMatrix<> in(4,2);

        in <<
              0,    0,
              0,    1,
              0,    0.5,
              0,    0.25;

        gsMatrix<> expected(4,2);

        expected <<
              0,    1,
              0,    0,
              0,    0.5,
              0,    0.75;

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        CHECK ( (expected - out.transpose()).norm() < 1.e-4 );
    }
    TEST(MatchingPseudoNonAffineTest)
    {
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,-1,1,0).release() );
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,0,1,1,180).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.computeTopology();

        CHECK ( mp.nInterfaces() == 1 );
        const boundaryInterface &bi = *(mp.iBegin());

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi,optHelper(gsRemapInterface<real_t>::neverAffine));

        CHECK ( !ri.isAffine() );
        CHECK ( ri.isMatching() );

        gsMatrix<> in(4,2);

        in <<
              0,    0,
              0,    1,
              0,    0.5,
              0,    0.25;

        gsMatrix<> expected(4,2);

        expected <<
              0,    1,
              0,    0,
              0,    0.5,
              0,    0.75;

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        CHECK ( (expected - out.transpose()).norm() < 1.e-4 );
    }
    TEST(MatchingAffineTest3D)
    {
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::lift3D(*gsNurbsCreator<>::BSplineRectangle(0,0,1,1)).release() );
        pc.push_back( gsNurbsCreator<>::lift3D(*gsNurbsCreator<>::BSplineRectangle(0,1,1,2)).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.computeTopology();

        CHECK ( mp.nInterfaces() == 1 );
        const boundaryInterface bi = mp.iBegin()->getInverse();
        
        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi,optHelper(5));

        CHECK ( ri.isAffine() );
        CHECK ( ri.isMatching() );

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

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        CHECK ( (expected - out.transpose()).norm() < 1.e-4 );
    }
    TEST(NonMatchingAffine)
    {
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,0,2,1).release() );
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0.5,1,1.5,2).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.addInterface( 1, boxSide(boundary::south), 0, boxSide(boundary::north));

        CHECK ( mp.nInterfaces() == 1 );
        const boundaryInterface &bi = *(mp.iBegin());

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi,optHelper(5));

        CHECK ( ri.isAffine() );
        CHECK ( !ri.isMatching() );

        gsMatrix<> in(3,2);

        in <<
              0.0,  0,
              0.5,  0,
              1.0,  0;

        gsMatrix<> expected(3,2);

        expected <<
              0.25,  1,
              0.5,   1,
              0.75,  1;

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        CHECK ( (expected - out.transpose()).norm() < 1.e-4 );
    }
    TEST(NonMatchingAffine2)
    {
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,0,2,1).release() );
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0.5,1,1.5,2).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.addInterface( 0, boxSide(boundary::north), 1, boxSide(boundary::south));

        CHECK ( mp.nInterfaces() == 1 );
        const boundaryInterface &bi = *(mp.iBegin());

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi,optHelper(5));

        CHECK ( ri.isAffine() );
        CHECK ( !ri.isMatching() );

        gsMatrix<> in(3,2);

        in <<
              0.25,  1,
              0.5,   1,
              0.75,  1;

        gsMatrix<> expected(3,2);

        expected <<
              0.0,  0,
              0.5,  0,
              1.0,  0;

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        CHECK ( (expected - out.transpose()).norm() < 1.e-4 );
    }
    TEST(NonMatchingAffine3)
    {
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,0,1,1).release() );
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0.5,1,1.5,2).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.addInterface( 0, boxSide(boundary::north), 1, boxSide(boundary::south));

        CHECK ( mp.nInterfaces() == 1 );
        const boundaryInterface &bi = *(mp.iBegin());

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi,optHelper(5));

        CHECK ( ri.isAffine() );
        CHECK ( !ri.isMatching() );

        gsMatrix<> in(3,2);

        in <<
              0.5,  1,
              0.6,  1,
              1.0,  1;

        gsMatrix<> expected(3,2);

        expected <<
              0.0,  0,
              0.1,  0,
              0.5,  0;

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        CHECK ( (expected - out.transpose()).norm() < 1.e-4 );
    }
    TEST(NonMatchingAffine3D)
    {
        // xlow, ylow, xup, yup, rotate
        std::vector< gsGeometry<>* > pc;
        pc.push_back( gsNurbsCreator<>::BSplineCube(1, 0,0,0).release() );
        pc.push_back( gsNurbsCreator<>::BSplineCube(1, 1,0.3,0.5 ).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.addInterface( 0, boxSide(boundary::east), 1, boxSide(boundary::west));

        CHECK ( mp.nInterfaces() == 1 );
        const boundaryInterface &bi = *(mp.iBegin());

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi,optHelper(5));

        CHECK ( ri.isAffine() );
        CHECK ( !ri.isMatching() );

        gsMatrix<> in(5,3);

        in <<
              1, 0.3, 0.5,
              1, 0.6, 0.5,
              1, 0.6, 1.0,
              1, 1.0, 1.0,
              1, 1.0, 0.5;

        gsMatrix<> expected(5,3);

        expected <<
              0, 0.0, 0.0,
              0, 0.3, 0.0,
              0, 0.3, 0.5,
              0, 0.7, 0.5,
              0, 0.7, 0.0;

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        CHECK ( (expected - out.transpose()).norm() < 1.e-4 );
    }
    TEST(NonSimpleAffineTest)
    {
        gsMultiPatch<> mp(*gsNurbsCreator<>::NurbsQuarterAnnulus(1,2));
        mp = mp.uniformSplit();

        CHECK ( mp.nInterfaces() == 4 );
        const boundaryInterface &bi = *(--mp.iEnd());

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi,optHelper(5));

        CHECK ( ri.isAffine() );
        CHECK ( ri.isMatching() );

        gsMatrix<> in(3,2);

        in <<
              0.8,  .5,
              1,    .5,
              0.5,  .5;

        gsMatrix<> expected(3,2);

        expected <<
              0.8,  .5,
              1,    .5,
              0.5,  .5;

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        CHECK ( (expected - out.transpose()).norm() < 1.e-4 );
    }
    TEST(NonAffineTest)
    {
        std::vector< gsGeometry<>* > pc;
        gsKnotVector<real_t> KV (0,1,0,3);
        gsMatrix<> C(9,2);
        C <<  0.00,0.00,  0.00,0.25,  0.00,1.00,
              0.25,0.00,  0.25,0.25,  0.25,1.00,
              1.00,0.00,  1.00,0.25,  1.00,1.00;

        pc.push_back( new gsTensorBSpline<2,real_t>(KV,KV, give(C)) );
        pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,1,1,2).release() );
        gsMultiPatch<> mp(pc); // consumes ptrs
        mp.addInterface( 1, boxSide(boundary::south), 0, boxSide(boundary::east));

        CHECK ( mp.nInterfaces() == 1 );

        const boundaryInterface &bi = *(mp.iBegin());

        gsMultiBasis<> mb(mp); // extract basis

        gsRemapInterface<real_t> ri(mp,mb,bi,optHelper(5));

        CHECK ( !ri.isAffine() );
        CHECK ( ri.isMatching() );

        gsMatrix<> in(5,2);

        in <<
              0,    0,
              0.25, 0,
              0.5,  0,
              0.75, 0,
              1,    0;

        gsMatrix<> expected(5,2);

        expected <<
              1,    0,
              1,    0.366025,
              1,    0.618034,
              1,    0.822876,
              1,    1;

        gsMatrix<> out;
        ri.eval_into(in.transpose(),out);

        CHECK ( (expected - out.transpose()).norm() < 1.e-4 );
    }
}
