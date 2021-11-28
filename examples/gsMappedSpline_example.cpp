#/** @file gsMappedSpline_example.cpp

    @brief Example using the gsMappedBasis and gsMappedSpline class

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


void createHBSpline(gsMultiPatch<> mp_coarser, gsMultiBasis<> mb_coarser, std::string name)
{
    gsMultiBasis<> mb_finer = mb_coarser;
    mb_finer.uniformRefine();

    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);

    A.setIntegrationElements(mb_finer);

    // Set the discretization space
    auto u = A.getSpace(mb_finer);

    gsMatrix<> mat_coarser, mat_finer;
    mat_coarser.setZero(mb_coarser.basis(0).size(),mb_finer.basis(0).size());
    mat_finer.setIdentity(mb_finer.basis(0).size(),mb_finer.basis(0).size());
    for (index_t bfID = 0; bfID < mb_coarser.basis(0).size(); bfID ++) {
        gsSingleBasis<real_t> sb(mb_coarser.basis(0), bfID);
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
    gsMatrix<real_t> point(2,1);
    point << 0.1, 0.1;
    gsMatrix<index_t> act = mb_finer.basis(0).active(point);
    gsMatrix<index_t> act2 = mb_coarser.basis(0).active(point);

    index_t rows = act.size() + mb_coarser.basis(0).size() - act2.size();
    gsMatrix<real_t> mat_hb(rows,mb_finer.basis(0).size());

    for (index_t j = 0; j < act.size(); j++) {
        index_t jj = act.at(j);
        mat_hb.row(j) = mat_finer.row(jj);
    }
    index_t row_shift = act.size();
    for (index_t i = 0; i < mb_coarser.basis(0).size(); i++) {

        bool not_active = true;
        for (index_t j = 0; j < act2.size(); j++) {
            index_t jj = act2.at(j);
            if (jj == i)
                not_active = false;
        }
        if (not_active)
        {
            mat_hb.row(row_shift) = mat_coarser.row(i);
            ++row_shift;
        }
    }
    //gsInfo << mat_hb << "\n";
    //gsInfo << mat_coarser << "\n";
    mat_hb = mat_hb.transpose();
    gsSparseMatrix<> sparse_sol = mat_hb.sparseView(1,1e-10);
    //gsInfo << sparse_sol << "\n";

    gsFileData<> fd;
    fd << mb_finer;
    fd << sparse_sol;
    fd.save(name + ".xml");

    mat_coarser = mat_coarser.transpose();
    gsSparseMatrix<> sparse_coarser = mat_coarser.sparseView(1,1e-10);
    //gsInfo << sparse_sol << "\n";

    fd.clear();
    fd << mp_coarser;
    fd << mb_finer;
    fd << sparse_coarser;
    fd.save( "spline_test.xml");
}



int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;

    std::string fn("msplines/spline_test");

    gsCmdLine cmd("Example using mapped spline.");
    cmd.addString( "f", "file", "Input XML file prefix", fn );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.getValues(argc,argv);
    //! [Parse command line]

    //! [Initialization]
    gsFileData<> fd;
    gsMultiPatch<real_t> mp;
    gsMultiBasis<real_t> mb;
    gsSparseMatrix<real_t> cf;

    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);
    A.setIntegrationElements(mb);
    //! [Initialization]

    //! [Read the input]
    fd.read(fn + ".xml");
    fd.getFirst(mp);
    fd.getFirst(mb);
    fd.getFirst(cf);
    //! [Read the input]


// TEMPORARLY TODO DELETE MAYBE

    //mb.uniformRefine();
    //createHBSpline(mp,mb, "hb_basis_2");

// TEMPORARLY TODO DELETE MAYBE



    //! [Setup the Mapped Basis]
    gsMappedBasis<2,real_t> mbasis(mb,cf);

    gsMultiBasis<> mbasis_exact(mp);

    gsInfo << "The MappedBasis has " << mbasis.size() << " basis functions for all patches! \n";
    gsInfo << "The SplineBasis has " << mbasis_exact.size() << " basis functions for all patches! \n";

    //! [Setup the Mapped Basis]

    //! [Some computation on the basis]
    // Works only for Bspline and identity
    gsMatrix<> points(2,3), result_mspline, result_bspline;
    points << 0,0, 0.2,0, 0.5,0.5;

    result_mspline = mbasis.basis(0).eval(points);
    result_bspline = mbasis_exact.basis(0).eval(points);
    gsInfo<<( (result_mspline-result_bspline).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    result_mspline = mbasis.basis(0).deriv(points);
    result_bspline = mbasis_exact.basis(0).deriv(points);
    gsInfo<<( (result_mspline-result_bspline).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    result_mspline = mbasis.basis(0).deriv2(points);
    result_bspline = mbasis_exact.basis(0).deriv2(points);
    gsInfo<<( (result_mspline-result_bspline).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    index_t bfID = 0;
    result_mspline = mbasis.basis(0).evalSingle(bfID, points);
    result_bspline = mbasis_exact.basis(0).evalSingle(bfID, points);
    gsInfo<<( (result_mspline-result_bspline).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    result_mspline = mbasis.basis(0).support(bfID);
    result_bspline = mbasis_exact.basis(0).support(bfID);
    gsInfo << result_mspline << "\n";
    gsInfo << result_bspline << "\n";
    gsInfo<<( (result_mspline-result_bspline).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    // TODO anchors_into() not implemented
    //result_mspline = mbasis.basis(0).anchors();
    //result_bspline = mbasis_exact.basis(0).anchors();
    //gsInfo<<( (result_mspline-result_bspline).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    bool failed = false;
    gsMatrix<index_t> act, act_exact;
    for (index_t i = 0; i < points.cols(); i++)
    {
        act = mbasis.basis(0).active(points.col(i));
        act_exact = mbasis_exact.basis(0).active(points.col(i));
        if ((act-act_exact).norm() > 1e-10)
            failed = true;
    }
    gsInfo<<( failed ? "failed" : "passed" )<<"\n";


    act = mbasis.basis(0).boundaryOffset(boxSide(boundary::west), 0);
    act_exact = mbasis_exact.basis(0).boundaryOffset(boxSide(boundary::west), 0);
    if (act.rows() != act_exact.rows())
        gsInfo << "failed: actives are different size \n";
    else
        gsInfo<<( (act-act_exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    act = mbasis.basis(0).boundaryOffset(boxSide(boundary::south), 1);
    act_exact = mbasis_exact.basis(0).boundaryOffset(boxSide(boundary::south), 1);
    if (act.rows() != act_exact.rows())
        gsInfo << "failed: actives are different size \n";
    else
        gsInfo<<( (act-act_exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    //! [Some computation on the basis]

    //! [Setup the Mapped Spline]
    gsMatrix<> coefs;
    coefs.setOnes(mbasis.size(),1);
    gsMappedSpline<2,real_t> mspline(mbasis,coefs);
    //! [Setup the Mapped Spline]

    //! [Simple L2 Projection example]
    if (cf.toDense().isIdentity(1e-10))
    {
        gsFunctionExpr<> functionExpr("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
        gsMatrix<> coefs2, coefs3;
        gsQuasiInterpolate<real_t>::localIntpl(mbasis.basis(0), functionExpr, coefs2);
        gsMappedSpline<2,real_t> mspline2(mbasis,coefs2);
        gsField<> solField2(mp, mspline2,true);
        gsWriteParaview<>( solField2, "MappedSpline_sol", 1000, false);

        gsQuasiInterpolate<real_t>::localIntpl(mb.basis(0), functionExpr, coefs3);
        gsMultiPatch<> mp2;
        mp2.addPatch(mb.basis(0).makeGeometry(coefs3));
        gsField<> solField3(mp, mp2,true);
        gsWriteParaview<>( solField3, "BSpline_sol", 1000, false);

        gsInfo<<( (coefs2-coefs3).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    }
    //! [Simple L2 Projection example]

    //! [Export visualization in ParaView]
    gsInfo<<"Plotting in Paraview...\n";
    gsWriteParaview<>( mbasis.basis(0), "MappedBasis", 1000, false);

    gsWriteParaview<>( mbasis_exact.basis(0), "SplineBasis", 1000, false);
    gsField<> solField(mp, mspline,true);
    gsWriteParaview<>( solField, "MappedSpline", 1000, false);
    //! [Export visualization in ParaView]
}
