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
# include <gsG1Basis/gsGlobalGDAssembler2.h>
# include <gsG1Basis/gsGlobalGDNorm.h>

namespace gismo
{

template<class T>
class gsApproxGluingData2
{
public:
    gsApproxGluingData2()
    { }

    gsApproxGluingData2(gsMultiPatch<T> const & mp,
                       gsMultiBasis<T> const & mb,
                       gsG1OptionList const & optionList)
        : m_mp(mp), m_mb(mb), m_optionList(optionList)
    {
        p_tilde = m_optionList.getInt("p_tilde");
        r_tilde = m_optionList.getInt("r_tilde");
    }

    // Computed the gluing data globally
    void setGlobalGluingData();

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
void gsApproxGluingData2<T>::setGlobalGluingData()
{
    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsKnotVector<T> kv(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD(kv);

    // First patch: the interface always at v
    gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_mb.basis(0).component(1));

    index_t degree = temp_basis_first.maxDegree();
    for (size_t i = degree+1; i < temp_basis_first.knots().size() - (degree+1); i += temp_basis_first.knots().multiplicityIndex(i))
        bsp_gD.insertKnot(temp_basis_first.knot(i),p_tilde-r_tilde);


    gsMatrix<> new_sol(bsp_gD.size() * 4,1);
    new_sol = new_sol.setOnes(); // startvector

    gsMatrix<> old_sol(bsp_gD.size() * 4,1);
    old_sol.setZero(); // different to startvector

    new_sol.block(bsp_gD.size()*2,0,bsp_gD.size()*2,1) = old_sol.bottomRows(bsp_gD.size()*2);

    gsGlobalGDNorm<T> globalGdNorm(m_mp, bsp_gD, new_sol);
    globalGdNorm.compute();
    real_t termination_value = globalGdNorm.value();

    gsInfo << "Initial sol: " << termination_value << "\n";

    while ( termination_value > 1e-03) // Different here?
    //while ( (old_sol - new_sol).norm() > 1e-6 )
    {
        old_sol = new_sol;

        gsGlobalGDAssembler2<T> globalGdAssembler2(m_mp, bsp_gD,m_optionList);
        globalGdAssembler2.assemble(old_sol);

        gsSparseSolver<real_t>::LU solver; // Matrix is non-symmetric
        solver.analyzePattern(globalGdAssembler2.matrix());
        solver.factorize(globalGdAssembler2.matrix());
        solver.compute(globalGdAssembler2.matrix());

        if(solver.info()!=Eigen::Success)
            gsInfo << "Solver failed \n";
        gsMatrix<> sol = solver.solve(globalGdAssembler2.rhs());
        if(solver.info()!=Eigen::Success)
            gsInfo << "Solver failed 2\n";
        // alpha^S
        solver.compute(globalGdAssembler2.matrix());
        new_sol = old_sol + solver.solve(globalGdAssembler2.rhs());

        //gsInfo << "matrix : " << globalGdAssembler2.matrix() << "\n";
        //gsInfo << "rhs : " << globalGdAssembler2.rhs() << "\n";
        //gsInfo << "sol : " << new_sol << "\n";

        gsGlobalGDNorm<T> globalGdNorm(m_mp, bsp_gD, new_sol );
        globalGdNorm.compute();
        termination_value = globalGdNorm.value();
        gsInfo << "Termination value: " << termination_value << "\n";

    }

    gsInfo << "sol : " << new_sol << "\n";

    gsGeometry<>::uPtr tilde_temp;

    tilde_temp = bsp_gD.makeGeometry(new_sol.block(0,0,bsp_gD.size(),1));
    gsBSpline<T> alpha_L = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    tilde_temp = bsp_gD.makeGeometry(new_sol.block(bsp_gD.size(),0,bsp_gD.size(),1));
    gsBSpline<T> alpha_R = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    tilde_temp = bsp_gD.makeGeometry(new_sol.block(2*bsp_gD.size(),0,bsp_gD.size(),1));
    gsBSpline<T> beta_L = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    tilde_temp = bsp_gD.makeGeometry(new_sol.block(3*bsp_gD.size(),0,bsp_gD.size(),1));
    gsBSpline<T> beta_R = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    gsWriteParaview(alpha_L,"alpha_L",2000);
    gsWriteParaview(alpha_R,"alpha_R",2000);
    gsWriteParaview(beta_R,"beta_R",2000);
    gsWriteParaview(beta_L,"beta_L",2000);


} // setGlobalGluingData



} // namespace gismo

