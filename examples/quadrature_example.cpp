/** @file quadrature_example.cpp

    @brief Playing with quadrature!

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

// Look also in tuturialBasis for more functionality of gsBSplineBasis.

#include <gismo.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    // ======================================================================
    // different construction of a knot vector
    // ======================================================================


    real_t a = 0; // starting knot
    real_t b = 1; // ending knot
    index_t interior = 4; // number of interior knots
    index_t multEnd = 3; // multiplicity at the two end knots
    index_t multInt = 1; // multiplicity at the interior knots
    bool plot = false;

    gsCmdLine cmd("Quadratire rules in G+Smo.");
    cmd.addReal("","starting","Starting knot",a);
    cmd.addReal("","ending","Ending knot",b);
    cmd.addInt("n","interior","Number of interior knots",interior);
    cmd.addInt("m","multI","Multiplicity at the interior knots",multInt);
    cmd.addInt("M","multE","Multiplicity at the two end knots",multEnd);
    cmd.addSwitch("plot","Plot with paraview",plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "------------- Constructions -----------------------------\n";

    gsKnotVector<> kv1(a, b, interior, multEnd, multInt);
    gsKnotVector<> kv2 = kv1;

    gsBSplineBasis<> bsb0(kv1);
    gsTensorBSplineBasis<2,real_t> tbsb(kv1,kv2);
    gsInfo<<tbsb<<"\n";


    gsWriteParaview(tbsb,"basis",1000,true);

    // ======================================================================
    // some properties
    // ======================================================================


    gsInfo << "------------- Some properties    -----------------------\n\n";

    gsInfo << "bsb0.size(): " << bsb0.size() << "\n\n"
              << "bsb0.numElements(): " << bsb0.numElements() << "\n\n"
              << "bsb0.degree(): " << bsb0.degree() << "\n\n";

    // printing some properties of the basis
    gsInfo << "Dimension of the parameter space: " << tbsb.dim() << "\n"
         << "Number of basis functions: " << tbsb.size() << "\n"
         << "Number of elements: " << tbsb.numElements() << "\n"
         << "Max degree of the basis: " << tbsb.maxDegree() << "\n"
         << "Min degree of the basis: " << tbsb.minDegree() << "\n"
         << "\n";

    // ======================================================================
    // some operations
    // ======================================================================

    // quadrature
    gsQuadRule<real_t> QuadRule;
    gsMatrix<> points;
    gsVector<> weights;

    gsOptionList options;
    options.addReal("quA", "Number of quadrature points: quA*deg + quB", 0  );
    options.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );

    QuadRule = gsQuadrature::get(tbsb, options);

    typename gsBasis<real_t>::domainIter domIt =  // add patchInd to domainiter ?
                tbsb.makeDomainIterator();

    gsMatrix<> allPoints(tbsb.dim(),tbsb.numElements()*QuadRule.numNodes());
    index_t k=0;
    for (; domIt->good(); domIt->next() )
    {

        gsInfo<<"Element corners:\n"<<domIt->lowerCorner().transpose()<<"\n"<<domIt->upperCorner().transpose()<<"\n";

        // Map the Quadrature rule to the element
        QuadRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                        points, weights);

        allPoints.block(0,k*QuadRule.numNodes(),tbsb.dim(),QuadRule.numNodes()) = points;
        k++;

        gsDebugVar(domIt->isBoundaryElement());

        gsDebugVar(points);
    }
        gsDebugVar(allPoints);

    gsWriteParaviewPoints(allPoints,"quadPoints");

    return 0;


    // typename gsBasis<real_t>::domainIter domIt =  // add patchInd to domainiter ?
    //             bsb0.makeDomainIterator();

    // gsOptionList options;
    // options.addReal("quA", "Number of quadrature points: quA*deg + quB", 1  );
    // options.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );

    // QuadRule = gsQuadrature::get(tbsb, options);
    // QuadRule.mapTo( a, b,points, weights);

    // gsDebugVar(points);

    // gsVector<> exact(bsb0.size());
    // exact.setZero();
    // gsMatrix<index_t> actives;
    // gsMatrix<real_t> values;
    // for (; domIt->good(); domIt->next() )
    // {

    //     gsInfo<<"Element corners:\n"<<domIt->lowerCorner().transpose()<<"\n"<<domIt->upperCorner().transpose()<<"\n";

    //     gsInfo<<domIt->side()<<"\n";

    //     // Map the Quadrature rule to the element
    //     QuadRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
    //                     points, weights);

    //     bsb0.active_into(points,actives);
    //     bsb0.eval_into(points,values);

    //     for (index_t p = 0; p!=actives.cols(); p++)
    //     {
    //         for (index_t r=0; r!=actives.rows(); r++)
    //         {
    //             exact.at(actives(r,p)) += weights.at(p) * values(r,p);
    //         }
    //     }
    // }
    // gsDebugVar(exact);

    // index_t ndof = bsb0.size();
    // index_t nquad = std::ceil(ndof/2.);

    // gsMatrix<> allPoints;
    // for (; domIt->good(); domIt->next() )
    // {

    //     gsInfo<<"Element corners:\n"<<domIt->lowerCorner().transpose()<<"\n"<<domIt->upperCorner().transpose()<<"\n";

    //     // Map the Quadrature rule to the element
    //     QuadRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
    //                     points, weights);

    //     bsb0.active_into(points,actives);
    //     bsb0.eval_into(points,values);

    //     for (index_t p = 0; p!=actives.cols(); p++)
    //     {
    //         for (index_t r=0; r!=actives.rows(); r++)
    //         {
    //             exact.at(actives(r,p)) += weights.at(p) * values(r,p);
    //         }
    //     }
    // }


    // points.resize(nquad,1);
    // points.col(0).setLinSpaced(nquad,a,b);
    // gsDebugVar(points);

    // bsb0.active_into(points.transpose(),actives);
    // bsb0.eval_into(points.transpose(),values);

    // gsMatrix<> shape(ndof,nquad);
    // shape.setZero();
    // for (index_t p = 0; p!=actives.cols(); p++)
    // {
    //     for (index_t r=0; r!=actives.rows(); r++)
    //     {
    //         shape(actives(r,p),p) = values(r,p);
    //     }
    // }
    // gsDebugVar(shape);

    // // gsVector<> resVec(2*nquad);
    // // res.segment(0,ndof) = shape * wq



    return 0;


    /*
        Other way of integration


        gsExprAssembler<> A;
        gsMultiBasis<> basis;
        basis.addBasis(&bsb0);
        A.setIntegrationElements(basis);
        typedef gsExprAssembler<>::space       space;
        gsExprEvaluator<> ev(A);
        space u = A.getSpace(basis);
        A.initSystem(true);
        A.assemble( u );
        gsDebugVar(A.rhs());
    */

}