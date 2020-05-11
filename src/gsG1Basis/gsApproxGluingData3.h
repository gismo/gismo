/** @file gsGluingData.h

    @brief Compute the gluing data for one interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

# include <gsG1Basis/gsG1OptionList.h>
# include <gsG1Basis/gsGlobalGDAssembler3.h>
# include <gsG1Basis/gsGlobalGDAssembler3Lambda.h>
# include <gsG1Basis/gsGlobalGDNorm3.h>

namespace gismo
{

template<class T>
class gsApproxGluingData3
{
public:
    gsApproxGluingData3()
    { }

    gsApproxGluingData3(gsMultiPatch<T> const & mp,
                        gsMultiBasis<T> const & mb,
                        gsG1OptionList const & optionList)
        : m_mp(mp), m_mb(mb), m_optionList(optionList)
    {
        p_tilde = m_optionList.getInt("p_tilde");
        r_tilde = m_optionList.getInt("r_tilde");
    }

    // Computed the gluing data globally
    void setGlobalGluingData();

    void setGlobalGluingDataWithLambda();

    const gsBSpline<T> get_alpha_tilde(index_t i) const {return alpha_tilde[i]; }
    const gsBSpline<T> get_beta_tilde(index_t i) const {return beta_tilde[i]; }

protected:

    // Spline space for the gluing data (p_tilde,r_tilde,k)
    index_t p_tilde, r_tilde;

    gsMultiPatch<> m_mp;
    gsMultiBasis<> m_mb;

    gsG1OptionList m_optionList;

    std::vector<gsBSpline<T>> alpha_tilde, beta_tilde;

}; // class gsGluingData


template<class T>
void gsApproxGluingData3<T>::setGlobalGluingData()
{
    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsKnotVector<T> kv(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD(kv);

    // First patch: the interface always at v
    gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_mb.basis(0).component(1));

    index_t degree = temp_basis_first.maxDegree();

    gsKnotVector<T> kv_beta(0,1,0,2*p_tilde+1,2*(p_tilde - r_tilde)); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD_beta(kv_beta);

    for (size_t i = degree+1; i < temp_basis_first.knots().size() - (degree+1); i += temp_basis_first.knots().multiplicityIndex(i))
    {
        bsp_gD.insertKnot(temp_basis_first.knot(i),p_tilde-r_tilde);
        bsp_gD_beta.insertKnot(temp_basis_first.knot(i),2*(p_tilde - r_tilde));
    }


    gsGlobalGDAssembler3<T> globalGdAssembler3(m_mp, bsp_gD, bsp_gD_beta, m_optionList);
    globalGdAssembler3.assemble();

    gsSparseSolver<real_t>::CGDiagonal solver; // Matrix is symmetric (?)
    solver.analyzePattern(globalGdAssembler3.matrix());
    solver.factorize(globalGdAssembler3.matrix());
    solver.compute(globalGdAssembler3.matrix());
    //gsInfo << "DET: " << solver.determinant() << "\n";
    if(solver.info()!=Eigen::Success)
        gsInfo << "Solver failed \n";
    gsMatrix<> sol = solver.solve(globalGdAssembler3.rhs());
    if(solver.info()!=Eigen::Success)
        gsInfo << "Solver failed 2\n";


    //gsInfo << "matrix : " << globalGdAssembler3.matrix() << "\n";
    //gsInfo << "rhs : " << globalGdAssembler3.rhs() << "\n";
    gsInfo << "sol : " << sol.transpose() << "\n";

    gsGlobalGDNorm3<T> globalGdNorm3(m_mp, bsp_gD, bsp_gD_beta, sol);
    globalGdNorm3.compute();
    gsInfo << "Termination value: " << globalGdNorm3.value() << "\n";


    gsGeometry<>::uPtr tilde_temp;

    tilde_temp = bsp_gD.makeGeometry(sol.block(0,0,bsp_gD.size(),1));
    gsBSpline<T> alpha_L = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    tilde_temp = bsp_gD.makeGeometry(sol.block(bsp_gD.size(),0,bsp_gD.size(),1));
    gsBSpline<T> alpha_R = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    tilde_temp = bsp_gD_beta.makeGeometry(sol.block(2*bsp_gD.size(),0,bsp_gD_beta.size(),1));
    gsBSpline<T> beta = dynamic_cast<gsBSpline<T> &> (*tilde_temp);


    gsWriteParaview(alpha_L,"alpha_L",2000);
    gsWriteParaview(alpha_R,"alpha_R",2000);
    gsWriteParaview(beta,"beta",2000);


} // setGlobalGluingData


template<class T>
void gsApproxGluingData3<T>::setGlobalGluingDataWithLambda()
{
    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsKnotVector<T> kv(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD(kv);

    gsKnotVector<T> kv_beta(0,1,0,2*p_tilde+1,2*(p_tilde-r_tilde)); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD_beta(kv_beta);

    // First patch: the interface always at v
    gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_mb.basis(0).component(1));

    index_t degree = temp_basis_first.maxDegree();
    for (size_t i = degree+1; i < temp_basis_first.knots().size() - (degree+1); i += temp_basis_first.knots().multiplicityIndex(i))
    {
        bsp_gD.insertKnot(temp_basis_first.knot(i),p_tilde-r_tilde);
        bsp_gD_beta.insertKnot(temp_basis_first.knot(i),2*(p_tilde-r_tilde));
    }

    gsMatrix<> new_sol(bsp_gD.size() * 2 + bsp_gD_beta.size() + 1,1);
    new_sol = new_sol.setOnes(); // startvector

    gsMatrix<> old_sol(new_sol.rows(),1);
    old_sol.setZero(); // different to startvector

    new_sol.block(bsp_gD.size()*2,0,bsp_gD_beta.size(),1) = old_sol.bottomRows(bsp_gD_beta.size()); // zero value for beta

    gsGlobalGDNorm3<T> globalGdNorm3(m_mp, bsp_gD, bsp_gD_beta, new_sol );
    globalGdNorm3.compute();
    real_t termination_value = globalGdNorm3.value();

    new_sol(bsp_gD.size() * 2 + bsp_gD_beta.size(),0) = 1;
    //gsInfo << "Initial sol: " << termination_value << "\n";
    //gsInfo << "Initial lambda: " << new_sol << "\n";

    //while ( (old_sol - new_sol).norm() > 1e-6 )
    while ( termination_value  > m_optionList.getReal("lambda") )
    {
        old_sol = new_sol;

        gsGlobalGDAssembler3Lambda<T> globalGDAssembler3Lambda(m_mp, bsp_gD, bsp_gD_beta, m_optionList);

        globalGDAssembler3Lambda.assemble(old_sol);

        //gsSparseSolver<real_t>::CGDiagonal solver; // Matrix is non-symmetric
        // alpha^S
        //solver.compute(globalGDAssembler3Lambda.matrix());
        //new_sol = solver.solve(globalGDAssembler3Lambda.rhs());


        gsSparseSolver<real_t>::CGDiagonal solver; // Matrix is symmetric (?)
        solver.analyzePattern(globalGDAssembler3Lambda.matrix());
        solver.factorize(globalGDAssembler3Lambda.matrix());
        solver.compute(globalGDAssembler3Lambda.matrix());
        if(solver.info()!=Eigen::Success)
            gsInfo << "Solver failed \n";
        new_sol = solver.solve(globalGDAssembler3Lambda.rhs());
        if(solver.info()!=Eigen::Success)
        {
            gsInfo << "Solver failed 2\n";
            gsInfo << new_sol.transpose() << "\n";
        }

/*
        gsSparseSolver<real_t>::BiCGSTABILUT solver;
        solver.compute(globalGDAssembler3Lambda.matrix());
        new_sol = solver.solveWithGuess(globalGDAssembler3Lambda.rhs(), old_sol);

        std::cout << "#iterations:     " << solver.iterations() << std::endl;
        std::cout << "estimated error: " << solver.error()      << std::endl;
*/
        new_sol.block(0,0, bsp_gD.size() * 2 + bsp_gD_beta.size(),1) += old_sol.block(0,0, bsp_gD.size() * 2 + bsp_gD_beta.size(),1);

        //gsInfo << "matrix : " << globalGDAssembler3Lambda.matrix() << "\n";
        //gsInfo << "rhs : " << globalGDAssembler3Lambda.rhs() << "\n";
        //gsInfo << "sol : " << new_sol << "\n";

        gsGlobalGDNorm3<T> globalGdNorm3(m_mp, bsp_gD, bsp_gD_beta, new_sol);
        globalGdNorm3.compute();
        termination_value = globalGdNorm3.value();
        gsInfo << "Termination value: " << termination_value << "\n";
        gsInfo << "Lambda value: " << new_sol.bottomRows(1) << "\n";
    }

    gsInfo << "sol : " << new_sol.topRows(10) << "\n";

    gsGeometry<>::uPtr tilde_temp;

    tilde_temp = bsp_gD.makeGeometry(new_sol.block(0,0,bsp_gD.size(),1));
    gsBSpline<T> alpha_L = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    tilde_temp = bsp_gD.makeGeometry(new_sol.block(bsp_gD.size(),0,bsp_gD.size(),1));
    gsBSpline<T> alpha_R = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    tilde_temp = bsp_gD_beta.makeGeometry(new_sol.block(2*bsp_gD.size(),0,bsp_gD_beta.size(),1));
    gsBSpline<T> beta = dynamic_cast<gsBSpline<T> &> (*tilde_temp);


    gsWriteParaview(alpha_L,"alpha_L",2000);
    gsWriteParaview(alpha_R,"alpha_R",2000);
    gsWriteParaview(beta,"beta",2000);


} // setGlobalGluingData

} // namespace gismo

