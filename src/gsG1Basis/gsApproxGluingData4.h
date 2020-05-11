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

# include <gsG1Basis/gsGlobalGDAssembler4.h>
# include <gsG1Basis/gsGlobalGDAssemberBeta4.h>
# include <gsG1Basis/gsGlobalGDNorm4.h>

namespace gismo
{

template<class T>
class gsApproxGluingData4
{
public:
    gsApproxGluingData4()
    { }

    gsApproxGluingData4(gsMultiPatch<T> const & mp,
                        gsMultiBasis<T> const & mb,
                        gsG1OptionList const & optionList)
        : m_mp(mp), m_mb(mb), m_optionList(optionList)
    {
        p_tilde = m_optionList.getInt("p_tilde");
        r_tilde = m_optionList.getInt("r_tilde");
    }

    // Computed the gluing data globally
    void setGlobalGluingData();

    const gsBSpline<T> get_alpha_tilde(index_t i) const {return alpha_S[i]; }
    const gsBSpline<T> get_beta_tilde(index_t i) const {return beta_S[i]; }

protected:

    // Spline space for the gluing data (p_tilde,r_tilde,k)
    index_t p_tilde, r_tilde;

    gsMultiPatch<> m_mp;
    gsMultiBasis<> m_mb;

    gsG1OptionList m_optionList;

    std::vector<gsBSpline<T>> alpha_S, beta_S;

}; // class gsGluingData


template<class T>
void gsApproxGluingData4<T>::setGlobalGluingData()
{
    // ALPHA^S

    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsKnotVector<T> kv(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD(kv);

    // First patch: the interface always at v
    gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_mb.basis(0).component(1));
    index_t degree = temp_basis_first.maxDegree();

    for (size_t i = degree+1; i < temp_basis_first.knots().size() - (degree+1); i += temp_basis_first.knots().multiplicityIndex(i))
    {
        bsp_gD.insertKnot(temp_basis_first.knot(i),p_tilde-r_tilde);
    }

    gsGlobalGDAssembler4<T> globalGdAssembler4(m_mp, bsp_gD, m_optionList);
    globalGdAssembler4.assemble();

    gsSparseSolver<real_t>::CGDiagonal solver; // Matrix is symmetric (?)
    solver.analyzePattern(globalGdAssembler4.matrix());
    solver.factorize(globalGdAssembler4.matrix());
    solver.compute(globalGdAssembler4.matrix());
    //gsInfo << "DET: " << solver.determinant() << "\n";
    gsMatrix<> sol = solver.solve(globalGdAssembler4.rhs());

    //gsInfo << "matrix : " << globalGdAssembler3.matrix() << "\n";
    //gsInfo << "rhs : " << globalGdAssembler3.rhs() << "\n";
    gsInfo << "sol : " << sol.transpose() << "\n";


    gsGeometry<>::uPtr tilde_temp;

    tilde_temp = bsp_gD.makeGeometry(sol.block(0,0,bsp_gD.size(),1));
    gsBSpline<T> alpha_L = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    tilde_temp = bsp_gD.makeGeometry(sol.block(bsp_gD.size(),0,bsp_gD.size(),1));
    gsBSpline<T> alpha_R = dynamic_cast<gsBSpline<T> &> (*tilde_temp);


    gsWriteParaview(alpha_L,"alpha_L",2000);
    gsWriteParaview(alpha_R,"alpha_R",2000);

    alpha_S.push_back(alpha_L);
    alpha_S.push_back(alpha_R);

    gsGlobalGDAssemblerBeta4<T> globalGdAssemblerBeta4(m_mp, bsp_gD, alpha_S, m_optionList);
    globalGdAssemblerBeta4.assemble();

    gsSparseSolver<real_t>::CGDiagonal solver_beta; // Matrix is symmetric (?)
    solver_beta.analyzePattern(globalGdAssemblerBeta4.matrix());
    solver_beta.factorize(globalGdAssemblerBeta4.matrix());
    solver_beta.compute(globalGdAssemblerBeta4.matrix());
    //gsInfo << "DET: " << solver.determinant() << "\n";
    sol = solver_beta.solve(globalGdAssemblerBeta4.rhs());

    //gsInfo << "matrix : " << globalGdAssemblerBeta4.matrix() << "\n";
    //gsInfo << "rhs : " << globalGdAssemblerBeta4.rhs() << "\n";
    gsInfo << "sol : " << sol.transpose() << "\n";


    tilde_temp = bsp_gD.makeGeometry(sol.block(0,0,bsp_gD.size(),1));
    gsBSpline<T> beta_L = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    tilde_temp = bsp_gD.makeGeometry(sol.block(bsp_gD.size(),0,bsp_gD.size(),1));
    gsBSpline<T> beta_R = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    gsWriteParaview(beta_L,"beta_L",2000);
    gsWriteParaview(beta_R,"beta_R",2000);

    beta_S.push_back(beta_L);
    beta_S.push_back(beta_R);


    gsGlobalGDNorm4<T> globalGdNorm4(m_mp, alpha_S, beta_S);
    globalGdNorm4.compute();
    gsInfo << "Termination value: " << globalGdNorm4.value() << "\n";


} // setGlobalGluingData


} // namespace gismo

