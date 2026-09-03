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


template <class T>
class funDerivative : public gsFunction<T>
{
public:
    funDerivative(const gsFunction<T> & fun)
    :
    m_fun(fun)
    {}

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        m_fun.eval_into(u,result);
    }

    void exact_deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        m_fun.deriv_into(u,result);
    }

    void exact_deriv2_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        m_fun.deriv2_into(u,result);
    }

    short_t domainDim() const { return m_fun.domainDim(); }
    short_t targetDim() const { return m_fun.targetDim(); }

protected:
    const gsFunction<T> & m_fun;

};

template <class T>
class basisDerivative : public gsFunction<T>
{
public:
    basisDerivative(const gsBasis<T> & basis, const gsMatrix<T> & points)
    :
    m_basis(basis)
    {
      m_actives = m_basis.active(points).rows();
    }

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        m_basis.eval_into(u,result);
    }

    void exact_deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        m_basis.deriv_into(u,result);
    }

    void exact_deriv2_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        m_basis.deriv2_into(u,result);
    }

    short_t domainDim() const { return m_basis.domainDim(); }
    short_t targetDim() const { return m_basis.targetDim()*m_actives; }

protected:
    const gsBasis<T> & m_basis;
    index_t m_actives;

};

int main(int argc, char *argv[])
{
    // Input options
    int numElevate    = 0;
    int numElevateMap = 0;
    int numRefine     = 1;
    int numRefineMap  = 1;
    bool plot         = false;

    // Arc length method options
    std::string bvp;

    index_t nmodes = 1;

    std::string wn("data.csv");

    gsCmdLine cmd("Shell modal solver.");

    cmd.addInt("r","ref", "Number of dyadic h-refinement (bisection) steps to perform before solving", numRefine);
    cmd.addInt("e","elev", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addInt("R","refM", "Number of dyadic h-refinement (bisection) steps to perform before solving", numRefineMap);
    cmd.addInt("E","elevM", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevateMap);

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    cmd.addString("i","inputFile", "Input file", bvp);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp;

    // ![Read data]
    gsFileData<real_t> fd(bvp);
    gsInfo<<"Reading geometry (ID=0) ...";
    fd.getId(0,mp);
    gsInfo<<"Finished\n";

    gsComposedGeometry<real_t> & cgeom = static_cast<gsComposedGeometry<real_t> &>(mp.patch(0));

    // cgeom.degreeElevate(numElevate);
    for (index_t k = 0; k < numRefine; k++)
        cgeom.uniformRefine();

    gsDebugVar(cgeom.basis());

    gsTensorBSpline<2,real_t> & composition = static_cast<gsTensorBSpline<2,real_t> &>(cgeom.composition());

    // composition.degreeElevate(numElevateMap);
    for (index_t k = 0; k < numRefineMap; k++)
        composition.uniformRefine();

    gsDebugVar(composition.basis());

    gsMatrix<> pts(2,3);
    pts.col(0).setConstant(0.25);
    pts.col(1).setConstant(0.50);
    pts.col(2).setConstant(0.75);

    // gsMatrix<> pts(2,1);
    // pts.col(0).setConstant(0.50);

    gsInfo<<"-------------------------------------------------------------------------\n";

    funDerivative<real_t> testG(cgeom);
    gsInfo<<"-------------------------------------------------------------------------\n";

    gsMatrix<> val, der_ex, der_num, der2_ex, der2_num;
    testG.eval_into(pts,val);
    gsDebugVar(val);
    testG.exact_deriv_into(pts,der_ex);
    gsDebugVar(der_ex);
    testG.deriv_into(pts,der_num);
    gsDebugVar(der_num);
    testG.exact_deriv2_into(pts,der2_ex);
    gsDebugVar(der2_ex);
    testG.deriv2_into(pts,der2_num);
    gsDebugVar(der2_num);

    gsInfo<<"-------------------------------------------------------------------------\n";

    basisDerivative<real_t> testB(mp.basis(0),pts);
    gsInfo<<"-------------------------------------------------------------------------\n";

    testB.eval_into(pts,val);
    gsDebugVar(val);
    testB.exact_deriv_into(pts,der_ex);
    gsDebugVar(der_ex);
    testB.deriv_into(pts,der_num);
    gsDebugVar(der_num);
    testB.exact_deriv2_into(pts,der2_ex);
    gsDebugVar(der2_ex);
    testB.deriv2_into(pts,der2_num);
    gsDebugVar(der2_num);


    gsDebugVar(composition.basis().active(pts));
    gsDebugVar(cgeom.basis().active(pts));

    gsMatrix<index_t> actives = cgeom.basis().active(pts);
    short_t d=cgeom.targetDim();
    gsMatrix<> gval(d,pts.cols()), gder(2*d,pts.cols()), gder2(3*d,pts.cols());
    gsMatrix<> bval(1,pts.cols()), bder(2,pts.cols()), bder2(3,pts.cols());
    gval.setZero();
    gder.setZero();
    gder2.setZero();
    for (index_t i = 0; i < actives.cols(); i++)
        for (index_t j = 0; j < actives.rows(); j++)
        {
            gsDebug<<"Basis function "<<actives(j,i)<<" on point "<<i<<":\n";
            bval = cgeom.basis().evalSingle(actives(j,i),pts.col(i));
            bder = cgeom.basis().derivSingle(actives(j,i),pts.col(i));
            bder2 = cgeom.basis().deriv2Single(actives(j,i),pts.col(i));
            gsDebug<<"* Value: \t"<<bval.transpose()<<"\n";
            gsDebug<<"* Deriv: \t"<<bder.transpose()<<"\n";
            gsDebug<<"* Deriv2:\t"<<bder2.transpose()<<"\n";

            for (index_t d=0; d!=cgeom.coefs().cols(); d++)
            {
                gval.block(d*1,i,1,1) += bval*cgeom.coefs()(actives(j,i),d);
                gder.block(d*2,i,2,1) += bder*cgeom.coefs()(actives(j,i),d);
                gder2.block(d*3,i,3,1) += bder2*cgeom.coefs()(actives(j,i),d);

            }
        }

    gsDebugVar(gval);
    gsDebugVar(gder);
    gsDebugVar(gder2);

    gsWriteParaview(mp,"mp",1000,true);

    return EXIT_SUCCESS;

}// end main
