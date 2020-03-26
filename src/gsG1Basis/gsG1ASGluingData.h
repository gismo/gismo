//
// Created by afarahat on 3/23/20.
//


#pragma once

#include <gsG1Basis/gsG1ASGluingDataAssembler.h>
#include <gsG1Basis/gsGluingData.h>
# include <gsG1Basis/gsG1OptionList.h>



namespace gismo
{

template<class T>
class gsG1ASGluingData : public gsGluingData<T>
{
public:
    gsG1ASGluingData()
    { }

    gsG1ASGluingData(gsMultiPatch<T> const & mp,
                 gsMultiBasis<T> const & mb,
                 index_t uv,
                 bool isBoundary,
                     gsG1OptionList const & optionList)
        : gsGluingData<T>(mp, mb, uv, isBoundary, optionList)
    {
        p_tilde = this->m_optionList.getInt("p_tilde");
        r_tilde = this->m_optionList.getInt("r_tilde");

        m_r = this->m_optionList.getInt("regularity");

        if (this->m_optionList.getSwitch("local"))
            gsInfo << "Is not yet implemented \n";
        else
            setGlobalGluingData();
    }


    // Computed the gluing data globally
    void setGlobalGluingData();



protected:

    // Spline space for the gluing data (p_tilde,r_tilde,k)
    index_t p_tilde, r_tilde;

    // Regularity of the geometry
    index_t m_r;

}; // class gsGluingData


    template<class T>
    void gsG1ASGluingData<T>::setGlobalGluingData()
    {
        // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
        gsKnotVector<T> kv(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
        gsBSplineBasis<T> bsp_gD(kv);


        gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(this->m_mb.basis(0).component(this->m_uv)); // u


        index_t degree = temp_basis_first.maxDegree();

        bsp_gD.insertKnot(temp_basis_first.knot(0),1); // Increase multiplicity of the first knot by one
        bsp_gD.insertKnot(temp_basis_first.knot(temp_basis_first.knots().size() - 1 ),1); // Increase multiplicity of the last knot by one

        if(temp_basis_first.knots().size() != (2 * (degree + 1)))  // If we have inner knots
        {
            for (size_t i = degree+1; i < temp_basis_first.knots().size() - (degree+1); i = i+(degree-m_r))
                bsp_gD.insertKnot(temp_basis_first.knot(i),2); // Increase the multiplicity of the inner knots by two
        }

        gsG1ASGluingDataAssembler<T> globalGdAssembler(bsp_gD, this->m_uv, this->m_mp, this->m_gamma, this->m_isBoundary);
        globalGdAssembler.assemble();

        gsSparseSolver<real_t>::CGDiagonal solver;
        gsVector<> sol_a, sol_b;

        // alpha^S
        solver.compute(globalGdAssembler.matrix_alpha());
        sol_a = solver.solve(globalGdAssembler.rhs_alpha());

        gsGeometry<>::uPtr tilde_temp;
        tilde_temp = bsp_gD.makeGeometry(sol_a);
        gsBSpline<T> alpha_tilde_L_2 = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
        this->alpha_tilde = alpha_tilde_L_2;

        // beta^S
        solver.compute(globalGdAssembler.matrix_beta());
        sol_b = solver.solve(globalGdAssembler.rhs_beta());

        tilde_temp = bsp_gD.makeGeometry(sol_b);
        gsBSpline<T> beta_tilde_L_2 = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
        this->beta_tilde = beta_tilde_L_2;

    } // setGlobalGluingData

} // namespace gismo


