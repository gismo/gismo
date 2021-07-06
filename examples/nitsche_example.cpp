/** @file nitsche_example.cpp

    @brief Some tests for gsVisitorNitsche

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <gismo.h>

using namespace gismo;


int main(int argc, char* argv[])
{

    real_t y = 1, x = 1;
    real_t alpha = 1, beta  = 1, penalty = -1;

    gsCmdLine cmd("dg_example");
    cmd.addReal("y","Domain.Height","",y);
    cmd.addReal("x","Domain.Width","",x);

    cmd.addReal("a","Nitsche.Alpha","",alpha);
    cmd.addReal("b","Nitsche.Beta","",beta);
    cmd.addReal("d","Nitsche.Penalty","",penalty);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << cmd.getOptionList() << "\n";

    // xlow, ylow, xup, yup, rotate
    std::vector< gsGeometry<>* > pc;
    pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,0,y,x).release() );
    gsMultiPatch<> mp(pc); // consumes ptrs
    mp.computeTopology();

    gsMultiBasis<> mb(mp); // extract basis

    {
        const boundary_condition<real_t>  bc(
                        /*patch*/ 0,
                        /*boxSide*/ boundary::north,
                        /*f_shptr*/ gsConstantFunction<>::make(1,2),
                        /*label*/ "Dirichlet",
                        /*unknown*/ 0,
                        /*unkcomp*/ 0, 
                        /*parametric*/ false
                        );


        gsOptionList opt = gsGenericAssembler<>::defaultOptions();
        opt.setReal( "Nitsche.Alpha", alpha );
        opt.setReal( "Nitsche.Beta", beta );
        opt.setReal( "Nitsche.Penalty", penalty );
        opt.setInt( "DirichletStrategy", dirichlet::nitsche );
        gsGenericAssembler<> ass(mp,mb,opt);
        gsSparseMatrix<> dgmat = ass.assembleNitsche(bc);

        gsInfo << std::fixed << std::setprecision(1) << dgmat.toDense() << "\n\n";

        gsInfo << ass.rhs().transpose() << "\n\n";

    }



    return 0;
}
