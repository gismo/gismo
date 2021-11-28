#/** @file gsMappedBasis_test.cpp

    @brief Tests using the gsMappedBasis and gsMappedSpline class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, P. Weinmueller, H. Verhelst
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]


template <class T>
class gsSingleBasis : public gismo::gsFunction<T>
{

protected:
    gsBasis<T> & _basis;
    mutable gsMapData<T> _tmp;
    index_t m_bfID;


public:
    /// Shared pointer for gsSingleBasis
    typedef memory::shared_ptr< gsSingleBasis > Ptr;

    /// Unique pointer for gsSingleBasis
    typedef memory::unique_ptr< gsSingleBasis > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsSingleBasis(gsBasis<T> & basis, index_t bfID) :
            _basis(basis), m_bfID(bfID), _basis_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN;
    }

    ~gsSingleBasis() { delete _basis_piece; }

GISMO_CLONE_FUNCTION(gsSingleBasis)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 1;}

    mutable gsSingleBasis<T> * _basis_piece; // why do we need this?

    const gsFunction<T> & piece(const index_t k) const
    {
        //delete _basis_piece;
        _basis_piece = new gsSingleBasis(_basis, m_bfID);
        return *_basis_piece;
    }

    // Input is parametric coordinates of 1-D \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( targetDim() , u.cols() );
        result = _basis.evalSingle(m_bfID, u);
    }
};



int main(int argc, char *argv[]) {
    //! [Parse command line]
    bool plot = false;
    index_t discreteDegree = 3;
    index_t numRefine = 1;

    std::string fn("msplines/spline_test");

    gsCmdLine cmd("Example using mapped spline.");
    cmd.addString("f", "file", "Input XML file prefix", fn);
    cmd.addInt("p", "degree", "Set the degree for the basis", discreteDegree);
    cmd.addInt("l", "loop", "Set the uniform refinement for the basis", numRefine);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.getValues(argc, argv);
    //! [Parse command line]

    //! [Initialization]
    gsMultiPatch<real_t> mp;
    gsMultiBasis<real_t> mb;
    gsSparseMatrix<real_t> cf;

    gsExprAssembler<> A(1, 1);
    gsExprEvaluator<> ev(A);
    A.setIntegrationElements(mb);
    //! [Initialization]

    //! [Read the input]
    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1, 1)); // patch 0
    // mp.embed(3);
    mp.computeTopology();

    mb = gsMultiBasis<>(mp);
    mb.setDegree(discreteDegree);
    mb.uniformRefine(numRefine);
    //! [Read the input]


    gsMultiBasis<> mb_mspline = mb;
    mb_mspline.uniformRefine();

    A.setIntegrationElements(mb_mspline);

    // Set the discretization space
    auto u = A.getSpace(mb_mspline);

    gsMatrix<> mat_coarser, mat_finer;
    mat_coarser.setZero(mb.basis(0).size(), mb_mspline.basis(0).size());
    mat_finer.setIdentity(mb.basis(0).size(), mb.basis(0).size());
    for (index_t bfID = 0; bfID < mb.basis(0).size(); bfID++) {
        gsSingleBasis<real_t> sb(mb.basis(0), bfID);
        auto aa = A.getCoeff(sb);

        gsBoundaryConditions<> bc_empty;
        u.setup(bc_empty, dirichlet::homogeneous, 0);
        A.initSystem();

        A.assemble(u * u.tr(), u * aa);

        gsSparseSolver<>::CGDiagonal solver;
        solver.compute(A.matrix());
        gsMatrix<> solVector = solver.solve(A.rhs());

        auto u_sol = A.getSolution(u, solVector);
        gsMatrix<> sol;
        u_sol.extract(sol);
        mat_coarser.row(bfID) = sol.transpose();
    }
    mat_coarser = mat_coarser.transpose();
    cf = mat_coarser.sparseView(1,1e-10);


    //! [Setup the Mapped Basis]
    gsMappedBasis<2, real_t> mbasis(mb_mspline, cf);

    gsInfo << "The MappedBasis has " << mbasis.size() << " basis functions for all patches! \n";
    gsInfo << "The SplineBasis has " << mb.size() << " basis functions for all patches! \n";

    //! [Setup the Mapped Basis]

    //! [Some computation on the basis]
    // Works only for Bspline and identity
    gsMatrix<> points(2, 4), result_mspline, result_bspline;
    points << 0, 0, 0.2, 0, 0.5, 0.5, 1, 1;

    result_mspline = mbasis.basis(0).eval(points);
    result_bspline = mb.basis(0).eval(points);
    gsInfo << ".eval(points) " << ((result_mspline - result_bspline).norm() < 1e-10 ? "passed" : "failed") << "\n";

    result_mspline = mbasis.basis(0).deriv(points);
    result_bspline = mb.basis(0).deriv(points);
    gsInfo << ".deriv(points) " << ((result_mspline - result_bspline).norm() < 1e-10 ? "passed" : "failed") << "\n";

    result_mspline = mbasis.basis(0).deriv2(points);
    result_bspline = mb.basis(0).deriv2(points);
    gsInfo << ".deriv2(points) " << ((result_mspline - result_bspline).norm() < 1e-10 ? "passed" : "failed") << "\n";

    bool failed = false;
    for (index_t bfID = 0; bfID < mb.basis(0).size(); bfID++)
    {
        result_mspline = mbasis.basis(0).evalSingle(bfID, points);
        result_bspline = mb.basis(0).evalSingle(bfID, points);
        if ((result_mspline - result_bspline).norm() > 1e-10)
            failed = true;
    }
    gsInfo << ".evalSingle(bfID, points) " << ( failed ? "failed" : "passed") << "\n";

    failed = false;
    for (index_t bfID = 0; bfID < mb.basis(0).size(); bfID++)
    {
        result_mspline = mbasis.basis(0).support(bfID);
        result_bspline = mb.basis(0).support(bfID);
        if ((result_mspline - result_bspline).norm() > 1e-10)
            failed = true;
    }
    gsInfo << ".support(bfID) " << ( failed ? "failed" : "passed") << "\n";

    real_t real_mspline, real_bspline;
    real_mspline = mbasis.basis(0).numElements();
    real_bspline = mb.basis(0).numElements();
    gsInfo << ".numElements() " << (real_mspline - real_bspline < 1e-10 ? "passed" : "failed") << "\n";

    // TODO anchors_into() not implemented
    gsInfo << "anchors_into() not implemented \n";
    //result_mspline = mbasis.basis(0).anchors();
    //result_bspline = mbasis_exact.basis(0).anchors();
    //gsInfo<<( (result_mspline-result_bspline).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    failed = false;
    gsMatrix<index_t> act, act_exact;
    for (index_t i = 0; i < points.cols(); i++) {
        act = mbasis.basis(0).active(points.col(i));
        act_exact = mb.basis(0).active(points.col(i));
        if ((act - act_exact).norm() > 1e-10)
            failed = true;
    }
    gsInfo << ".active(points.col(i)) " << (failed ? "failed" : "passed") << "\n";

    failed = false;
    for (index_t i = 1; i < 5; i++) {
        act = mbasis.basis(0).boundaryOffset(boxSide(i), 0);
        act_exact = mb.basis(0).boundaryOffset(boxSide(i), 0);
        if (act.rows() != act_exact.rows())
        {
            gsInfo << ".boundaryOffset(boxSide("<< i <<"), 0) " << "failed: actives have different size \n";
            failed = true;
        }
        else if ((act - act_exact).norm() > 1e-10)
            failed = true;
    }
    gsInfo << ".boundaryOffset(boxSide(i), 0) " << (failed ? "failed" : "passed") << "\n";

    failed = false;
    for (index_t i = 1; i < 5; i++) {
        act = mbasis.basis(0).boundaryOffset(boxSide(i), 1);
        act_exact = mb.basis(0).boundaryOffset(boxSide(i), 1);
        if (act.rows() != act_exact.rows())
        {
            //gsDebugVar(act);
            //gsDebugVar(act_exact);
            gsInfo << ".boundaryOffset(boxSide("<< i <<"), 1) " << "failed: actives have different size \n";
            failed = true;
        }
        else if ((act - act_exact).norm() > 1e-10)
            failed = true;
    }
    gsInfo << ".boundaryOffset(boxSide(i), 1) " << (failed ? "failed" : "passed") << "\n";

    result_mspline = mbasis.basis(0).component(0).eval(points);
    result_bspline = mb.basis(0).component(0).eval(points);
    gsInfo << ".component(0).eval(points) " << ((result_mspline - result_bspline).norm() < 1e-10 ? "passed" : "failed") << "\n";

    for (index_t i = 0; i < points.cols(); i++) {
        act = mbasis.basis(0).component(1).active(points.col(i));
        act_exact = mb.basis(0).component(1).active(points.col(i));
        if ((act - act_exact).norm() > 1e-10)
            failed = true;
    }
    gsInfo << ".component(1).active(points.col(i)) " << (failed ? "failed" : "passed") << "\n";

    failed = false;
    for (index_t bfID = 0; bfID < mb.basis(0).size(); bfID++)
    {
        result_mspline = mbasis.basis(0).function(bfID).eval(points);
        result_bspline = mb.basis(0).function(bfID).eval(points);
        if ((result_mspline - result_bspline).norm() > 1e-10)
            failed = true;
    }
    gsInfo << ".function(bfID).eval(points) " << ( failed ? "failed" : "passed") << "\n";
    //! [Some computation on the basis]

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( mbasis.basis(0), "MappedBasis", 1000, false);

        gsWriteParaview<>( mb.basis(0), "SplineBasis", 1000, false);
    }

    return EXIT_SUCCESS;
}