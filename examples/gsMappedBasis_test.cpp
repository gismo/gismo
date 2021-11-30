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
    //mb_mspline.degreeElevate(1);
    mb_mspline.reduceContinuity(1);

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
    //gsDebugVar(cf.toDense());

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
            gsInfo << ".boundaryOffset(boxSide("<< i <<"), 1) " << "failed: actives have different size \n";
            failed = true;
        }
        else if ((act - act_exact).norm() > 1e-10)
            failed = true;
    }
    gsInfo << ".boundaryOffset(boxSide(i), 1) " << (failed ? "failed" : "passed") << "\n";

    gsInfo << "component(i) not working \n";
    //result_mspline = mbasis.basis(0).component(0).eval(points);
    //result_bspline = mb.basis(0).component(0).eval(points);
    //gsInfo << ".component(0).eval(points) " << ((result_mspline - result_bspline).norm() < 1e-10 ? "passed" : "failed") << "\n";
/*
    for (index_t i = 0; i < points.cols(); i++) {
        act = mbasis.basis(0).component(1).active(points.col(i));
        act_exact = mb.basis(0).component(1).active(points.col(i));
        if ((act - act_exact).norm() > 1e-10)
            failed = true;
    }
    gsInfo << ".component(1).active(points.col(i)) " << (failed ? "failed" : "passed") << "\n";
*/
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

    //! [Computing the kernel for offset(s,1)]
    gsMatrix<index_t> act2;

    boxSide side = boxSide(3);
    //act = mbasis.basis(0).boundaryOffset(side, 0); // since offset 0 is deleted from offset 1.
    act = mbasis.basis(0).boundaryOffset(side, 1);

    gsInfo << "act: " << act << "\n";

    // Four 2D-points
    index_t numPoints = 7;
    gsMatrix<> points2D(2,numPoints), normal(2,1), res_deriv, res_normalDer;
    points2D.setZero();
    gsVector<> linspace(numPoints);
    linspace.setLinSpaced(numPoints,0.0,1.0);
    points2D.row(0) = linspace.transpose();

    normal <<  0, -1;
    res_deriv = mbasis.basis(0).deriv(points2D);
    for (index_t i = 0; i < act.size(); i++) {
        index_t ii = act.at(i);
        result_mspline = mbasis.basis(0).function(ii).eval(points2D);
        res_normalDer.setZero(1, numPoints);
        for (index_t j = 0; j < numPoints; j++)
            res_normalDer.col(j) = mbasis.basis(0).function(ii).deriv(points2D.col(j)).transpose() * normal;

        gsInfo << "Bf " << ii << " Value:\n";
        gsInfo << result_mspline << "\n";
        gsInfo << "Bf " << ii << " Normal Der:\n";
        gsInfo << res_normalDer << "\n";

    }
    // Only working if the basis functions are linear dependence!
    bool linearDepedence = false;
    if (linearDepedence) {
        gsMatrix<index_t> act2;

        boxSide side = boxSide(3);
        act = mbasis.basis(0).boundaryOffset(side, 0); // since offset 0 is deleted from offset 1.
        act2 = mbasis.basis(0).boundaryOffset(side, 1);

        act.conservativeResize(act.rows() + act2.rows(), 1);
        act.bottomRows(act2.rows()) = act2;

        index_t dir = side < 3 ? 1 : 0;
        index_t dim_u = mbasis.getBase(0).component(0).size(); // Size of the underlying basis functions
        index_t dim_v = mbasis.getBase(0).component(1).size(); // Size of the underlying basis functions

        gsMatrix<> coefs_mat(dir == 0 ? 2 * dim_u : 2 * dim_v, act.size()); // the first two rows/cols
        for (index_t i = 0; i < act.size(); i++) {
            index_t ii = act.at(i);
            if (dir == 0) {
                coefs_mat.block(0, i, dim_u, 1) = cf.block(0, ii, dim_u, 1); // Todo replace with mbasis.matrix() ?
                coefs_mat.block(dim_u, i, dim_u, 1) = cf.block(dim_u, ii, dim_u,
                                                               1); // Todo replace with mbasis.matrix() ?
            } else if (dir == 1) {
                // Not implemented yet.
            }
        }

        real_t threshold = 1e-10;
        Eigen::FullPivLU<gsMatrix<>> KernelCorner(coefs_mat);
        KernelCorner.setThreshold(threshold);

        gsMatrix<> basis_trafo;
        basis_trafo.setIdentity(act.size(), act.size());

        gsMatrix<> kernel = KernelCorner.kernel();

        size_t count = 0;
        while (kernel.cols() < act.size()) {
            kernel.conservativeResize(kernel.rows(), kernel.cols() + 1);
            kernel.col(kernel.cols() - 1) = basis_trafo.col(count);

            Eigen::FullPivLU<gsMatrix<>> ker_temp(kernel);
            ker_temp.setThreshold(threshold);
            if (ker_temp.dimensionOfKernel() != 0) {
                kernel = kernel.block(0, 0, kernel.rows(), kernel.cols() - 1);
            }
            count++;
        }

        gsMatrix<> coef_new = cf.toDense();
        for (index_t j = 0; j < act.size(); ++j) {
            coef_new.col(act.at(j)).setZero();
            for (index_t i = 0; i < act.size(); ++i)
                if (kernel(i, j) * kernel(i, j) > 1e-25)
                    coef_new.col(act.at(j)) += cf.col(act.at(i)) * kernel(i, j);
        }

        gsMappedBasis<2, real_t> mbasis2(mb_mspline, cf);
        gsWriteParaview<>(mbasis2.basis(0), "MappedBasis2", 1000, false);
    }
    //! [Computing the kernel for offset(s,1)]

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( mbasis.basis(0), "MappedBasis", 1000, false);

        gsWriteParaview<>( mb.basis(0), "SplineBasis", 1000, false);
    }

    return EXIT_SUCCESS;
}