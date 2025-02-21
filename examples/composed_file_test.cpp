/** @file composed_domain_test.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsNurbs/gsMobiusDomain.h>

using namespace gismo;
//! [Include namespace]

template <short_t _DIM, class T>
gsTensorBSplineBasis<_DIM,T> integrationBasis(const gsTensorBSplineBasis<_DIM,T> & basis1,
                                           const gsTensorBSplineBasis<_DIM,T> & basis2)
{
    gsTensorBSplineBasis<_DIM,T> ibasis(basis1);
    // Integration basis: parent basis with knots of composition basis inserted, and the degree is the sum of the two degrees (?)
    index_t targetDegree;
    for (size_t d = 0; d!=_DIM; d++)
    {
        // 1. Insert interior knots of composition basis
        for (typename gsKnotVector<T>::uiterator it = std::next(basis2.knots(d).ubegin());
                                                    it!= std::prev(basis2.knots(d).uend());
                                                    ++it)
            {
                if (ibasis.knots(d).has(*it))
                    continue;
                ibasis.insertKnot(*it,d);
            }
        // 2. Increase the degree
        targetDegree = ibasis.degree(d) * basis2.degree(d);
        ibasis.degreeIncrease(targetDegree-ibasis.degree(d),d);

    }
    return ibasis;
}

int main(int argc, char *argv[])
{
    //! [Parse command line]
    std::string fn = "../filedata/monitor_results/monitor_example_face_nn_r0_e1_R3_E1/composedGeometry.xml";

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addString( "f", "file", "Input XML file", fn );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    GISMO_ASSERT(!fn.empty(),"Please provide an input file.");

    //! [Read input file]
    gsFileData<> fd(fn);
    GISMO_ASSERT(fd.hasAny<gsGeometry<real_t>>(),"The input file must contain a geometry.");
    gsGeometry<>::uPtr geom = fd.getFirst<gsGeometry<>>();

    // Extract the composition
    GISMO_ASSERT((dynamic_cast<gsComposedGeometry<real_t>*>(geom.get())),"The geometry must be a composed geometry.");
    gsComposedGeometry<>::uPtr cgeom = memory::make_unique(static_cast<gsComposedGeometry<real_t>*>(geom.release()));
    gsComposedBasis<>::uPtr    cbasis= memory::make_unique(cgeom->basis().clone().release());

    GISMO_ASSERT((dynamic_cast<gsGeometry<real_t>*>(&cgeom->composition())),"The composition must be a geometry.");
    gsGeometry<>::uPtr composition = memory::make_unique(static_cast<gsGeometry<real_t>*>(&cgeom->composition())->clone().release());
    // gsGeometry<>::uPtr geometry    = memory::make_unique(cgeom->geometry().clone().release());
    gsBasis<>::uPtr    basis       = memory::make_unique(cbasis->basis().clone().release());
    gsGeometry<>::uPtr geometry    = basis->makeGeometry(cgeom->coefs());

    gsInfo<<"Composition:\n"<<*composition<<"\n";
    gsInfo<<"Geometry:\n"<<*geometry<<"\n";
    gsInfo<<"Basis:\n"<<*basis<<"\n";


    gsWriteParaview(*geometry,"geometry");
    gsWriteParaview(*cgeom,"cgeometry");

    GISMO_ASSERT((dynamic_cast<gsTensorBSplineBasis<2>*>(basis.get())),"The basis must be a tensor basis.");
    gsTensorBSplineBasis<2> * tbasis = static_cast<gsTensorBSplineBasis<2>*>(basis.get());
    GISMO_ASSERT((dynamic_cast<gsTensorBSplineBasis<2>*>(&composition->basis())),"The composition must be a tensor basis.");
    gsTensorBSplineBasis<2> * composition_basis = static_cast<gsTensorBSplineBasis<2>*>(&composition->basis());

    gsTensorBSplineBasis<2> ibasis = integrationBasis(*tbasis,*composition_basis);

    gsExprEvaluator<> ev;
    gsMultiBasis<> mb(ibasis);
    ev.setIntegrationElements(mb);
    auto G = ev.getMap(*geometry);
    auto cG = ev.getMap(*cgeom);

    auto GArea = ev.integral(meas(G));
    auto cGArea = ev.integral(meas(cG));

    gsInfo<<"Area = "<< GArea  <<"\n";
    gsInfo<<"Area = "<< cGArea <<"\n";

    /*
    gsExprEvaluator<T> evaluator;

    evaluator.setIntegrationElements(m_mb); // does not work when in constructor
    m_comp->setControls(u);

    // Penalty constant
    gsConstantFunction<T> pen(m_options.getReal("Penalty"), m_cgeom.domainDim());
    geometryMap G = evaluator.getMap(m_mp);
    */

    //jacobian determinant for a surface, i.e. the measure
    auto fform = jac(G).tr()*jac(G);
    auto detG = pow(fform.det().val(),0.5); 
    auto G_frob = jac(G).norm();

    

    auto cfform = jac(cG).tr()*jac(G);
    auto detcG = pow(cfform.det().val(),0.5); 
    auto cG_frob = jac(cG).norm();

    gsInfo << " ----------------- G ---------- cG ---------- \n";
    
    
    //ev.eval(expr,p)
    gsInfo << "Area distortion: " << ev.integral(detG*meas(G))/GArea << ' ----- '<< ev.integral(detcG*meas(cG))/cGArea << "\n";
    //gsInfo << "Angular distortion:" << detG/GArea << ' ----- '<< detcG/cGArea << "\n";

    return EXIT_SUCCESS;

}// end main
