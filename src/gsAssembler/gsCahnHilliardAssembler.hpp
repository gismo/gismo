/** @file gsCahnHilliardAssembler.hpp

    @brief Provides assemblers for the Cahn-Hilliard equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        m_assembler. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsAssembler/gsCahnHilliardAssembler.h>
#include <gsMSplines/gsMappedBasis.h>
#include <gsMSplines/gsMappedSpline.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsCore/gsFunctionExpr.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsPiecewiseFunction.h>

namespace gismo
{

template <class T>
gsCahnHilliardAssembler<T>::gsCahnHilliardAssembler(const gsMultiPatch<T> & mp,
                                                    const gsMultiBasis<T> & mb,
                                                    const gsBoundaryConditions<T> & bcs
                                                    )
:
m_patches(mp),
m_integrationBasis(&mb),
m_spaceBasis(&mb),
m_bcs(bcs),
m_initialized(false)
{
    this->_defaultOptions();
};

template <class T>
gsCahnHilliardAssembler<T>& gsCahnHilliardAssembler<T>::operator=( const gsCahnHilliardAssembler& other )
{
    if (this!=&other)
    {
        m_penalty=other.m_penalty;
        m_lambda=other.m_lambda;
        // m_continuity=other.m_continuity;

        m_patches=other.m_patches;
        m_integrationBasis=other.m_integrationBasis;
        m_spaceBasis=other.m_spaceBasis;
        m_bcs=other.m_bcs;

        m_options=other.m_options;

        // To do: make copy constructor for the gsExprAssembler
        m_assembler.setIntegrationElements(*m_integrationBasis);
        m_assembler.setOptions(m_options);
    }
    return *this;
}

template <class T>
gsCahnHilliardAssembler<T>& gsCahnHilliardAssembler<T>::operator=( gsCahnHilliardAssembler&& other )
{
    m_penalty=give(other.m_penalty);
    m_lambda=give(other.m_lambda);
    // m_continuity=give(other.m_continuity);

    m_patches=give(other.m_patches);
    m_integrationBasis=give(other.m_integrationBasis);
    m_spaceBasis=give(other.m_spaceBasis);
    m_bcs=give(other.m_bcs);

    m_options=give(other.m_options);

    // To do: make copy constructor for the gsExprAssembler
    m_assembler.setIntegrationElements(*m_integrationBasis);
    m_assembler.setOptions(m_options);
    return *this;
}

template <class T>
void gsCahnHilliardAssembler<T>::_defaultOptions()
{
    m_options.addReal("Lambda","Lambda parameter",1e0);
    m_options.addInt("Continuity","Continuity between patches: C^{-1} (-1) or C^0 (0, default)",0);
    m_options.addReal("Penalty","Penalty parameter for Nitsche boundary conditions (default: 1e4)",1e4);
    m_options.addSwitch("AssembleWeakBCs","Assemble Nitsche boundary conditions in every iteration",false);

    /* UNUSED:
    m_options.addInt("Mobility","Mobility function: 0 for constant, 1 for double well",0);
    m_options.addReal("M0","M0 parameter",1e0);
     */
    // Assembler options
    gsOptionList assemblerOptions = m_assembler.defaultOptions().wrapIntoGroup("ExprAssembler");
    m_options.update(assemblerOptions,gsOptionList::addIfUnknown);
}

template <class T>
void gsCahnHilliardAssembler<T>::_getOptions()
{
    m_lambda = m_options.getReal("Lambda");
    m_continuity = m_options.getInt("Continuity");
    m_penalty = m_options.getReal("Penalty");

    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));
}

template<class T>
void gsCahnHilliardAssembler<T>::setOptions(gsOptionList & options)
{
    m_options.update(options,gsOptionList::ignoreIfUnknown);
}

template <class T>
void gsCahnHilliardAssembler<T>::initialize()
{
    this->_getOptions();

    // Elements used for numerical integration
    m_assembler.setIntegrationElements(*m_integrationBasis);
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    GISMO_ASSERT(m_bcs.hasGeoMap(),"No geometry map was assigned to the boundary conditions. Use bc.setGeoMap to assign one!");

    typename gsExprAssembler<T>::space u = m_assembler.getSpace(*m_spaceBasis, 1, 0); // last argument is the space ID

    // Setup the system
    u.setup(m_bcs,
            m_assembler.options().askInt("DirichletValues",dirichlet::l2Projection),
            m_options.askInt("Continuity",-1));
    m_assembler.initSystem();

    // Compute sparsity patter: this is done automatically - but
    // is needed if assemble(.) is called twice
    m_assembler.computePattern( igrad(u) * igrad(u).tr() );

    m_initialized = true;
}

template <class T>
void gsCahnHilliardAssembler<T>::assemble(const gsFunctionSet<T> & C)
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    // m_assembler.clearRhs(); // Resets to zero the values of the already allocated to residual (RHS)
    m_assembler.initVector();
    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Get the solution and its derivative
    auto c  = m_assembler.getCoeff(C);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    // Mobility
    auto M_c  = 1.0 + 0.0*c.val(); // replace with const_expr(1.0) instead of using 0*c
    // auto dM_c = 0.0 * igrad(c,G); // replace with const_expr(1.0) instead of using 0*c!!

    // Psi(c) = 1/4*(c^2-1)^2 and its derivatives
    // Psi'(c) = c^3 - c
    // Psi''(c) = 3*c^2 - 1
    auto psi_prime = (c*c*c).val() - c.val();    // c^3 - c
    auto psi_double_prime = 3.0*(c*c).val() - 1.0; // 3c^2 - 1
    auto psi_triple_prime = 6.0*c.val();           // 6c

    // Weak form residual
    // mu = psi'(c) - lambda*laplacian(c)
    // term1: M*grad(w)*psi''(c)*grad(c)
    // term2: M*lapl(w)*lambda*lapl(c)
    auto residual = M_c.val() * igrad(w,G) * psi_double_prime * igrad(c,G).tr() + // F_bar
                    M_c.val() * ilapl(w,G) * m_lambda * ilapl(c,G).val(); // K_laplacian

    auto jacobian = M_c.val() * psi_double_prime * igrad(w,G) * igrad(w,G).tr() + // K_f1
                    M_c.val() * psi_triple_prime * igrad(w,G) * igrad(c,G).tr() * w.tr() + // K_f2
                    M_c.val() * m_lambda * ilapl(w,G) * ilapl(w,G).tr(); // K_laplacian

    m_assembler.assemble(jacobian * meas(G), residual * meas(G));

    // ASSEMBLE NITSCHE
    if (m_options.getSwitch("AssembleWeakBCs"))
    {
        if (m_bcs.get("Weak Clamped").size()==0)
            gsWarn<<"Nitsche boundary assembly is requested, but no boundaries are marked 'Weak Clamped'";

        T hmax = 0;
        for (size_t p=0; p!=m_patches.nPatches(); p++)
            hmax = math::max(hmax, m_patches.basis(p).getMaxCellLength());

        auto weak_bc_residual = - igrad(w,G) * nv(G) * m_lambda * ilapl(c,G).val()
                                +(igrad(w,G) * nv(G).normalized()) * hmax * m_penalty * (igrad(c,G) * nv(G))
                                - m_lambda * ilapl(w,G) * igrad(c,G) * nv(G);
        auto weak_bc_jacobian = - m_lambda * igrad(w,G) *  nv(G)  * ilapl(w,G).tr() + // consistency term
                                m_penalty * (igrad(w,G) * nv(G).normalized()) * hmax * (igrad(w,G) * nv(G)).tr() - // penalty (stabilizing) term
                                m_lambda * ilapl(w,G) * (igrad(w,G)  * nv(G)).tr(); // symmetry term

        m_assembler.assembleBdr(m_bcs.get("Weak Clamped"), weak_bc_jacobian, weak_bc_residual);
    }
}

template <class T>
void gsCahnHilliardAssembler<T>::assembleResidual(const gsFunctionSet<T> & C)
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.initVector();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Get the solution
    auto c  = m_assembler.getCoeff(C);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    // Mobility
    auto M_c  = 1.0 + 0.0*c.val(); // replace with const_expr(1.0) instead of using 0*c
    // auto dM_c = 0.0 * igrad(c,G); // replace with const_expr(1.0) instead of using 0*c!!

    // Psi(c) = 1/4*(c^2-1)^2 and its derivatives
    // Psi'(c) = c^3 - c
    // Psi''(c) = 3*c^2 - 1
    auto psi_prime = (c*c*c).val() - c.val();    // c^3 - c
    auto psi_double_prime = 3.0*(c*c).val() - 1.0; // 3c^2 - 1

    // Weak form residual
    // mu = psi'(c) - lambda*laplacian(c)
    // term1: M*grad(w)*psi''(c)*grad(c)
    // term2: M*lapl(w)*lambda*lapl(c)
    auto residual = M_c.val() * igrad(w,G) * psi_double_prime * igrad(c,G).tr() + // F_bar
                    M_c.val() * ilapl(w,G) * m_lambda * ilapl(c,G).val(); // K_laplacian

    m_assembler.assemble(residual * meas(G));

    // ASSEMBLE NITSCHE
    if (m_options.getSwitch("AssembleWeakBCs"))
    {
        if (m_bcs.get("Weak Clamped").size()==0)
            gsWarn<<"Nitsche boundary assembly is requested, but no boundaries are marked 'Weak Clamped'";

        T hmax = 0;
        for (size_t p=0; p!=m_patches.nPatches(); p++)
            hmax = math::max(hmax, m_patches.basis(p).getMaxCellLength());

        m_assembler.assembleBdr(m_bcs.get("Weak Clamped"),
            - igrad(w,G) * nv(G) * m_lambda * ilapl(c,G).val()
            +(igrad(w,G) * nv(G).normalized()) * hmax * m_penalty * (igrad(c,G) * nv(G))
            - m_lambda * ilapl(w,G) * igrad(c,G) * nv(G));
    }
    // NOTE: no m_assembler.cleanUp() here — parse() handles cleanup at the start of the next call
}

template <class T>
void gsCahnHilliardAssembler<T>::assembleMassMatrix()
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    // Initialize the system
    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.assemble(w*w.tr()*meas(G));// K_m
}

template <class T>
void gsCahnHilliardAssembler<T>::assembleJacobian(const gsFunctionSet<T> & C)
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

        // Get the solution and its derivative
    auto c  = m_assembler.getCoeff(C);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    // Mobility
    auto M_c  = 1.0 + 0.0*c.val(); // replace with const_expr(1.0) instead of using 0*c
    // auto dM_c = 0.0 * igrad(c,G); // replace with const_expr(1.0) instead of using 0*c!!

    // Psi(c) = 1/4*(c^2-1)^2 and its derivatives
    // Psi'(c) = c^3 - c
    // Psi''(c) = 3*c^2 - 1
    auto psi_prime = (c*c*c).val() - c.val();    // c^3 - c
    auto psi_double_prime = 3.0*(c*c).val() - 1.0; // 3c^2 - 1
    auto psi_triple_prime = 6.0*c.val();           // 6c

    // Weak form jacobian
    auto jacobian = M_c.val() * psi_double_prime * igrad(w,G) * igrad(w,G).tr() + // K_f1
                    M_c.val() * psi_triple_prime * igrad(w,G) * igrad(c,G).tr() * w.tr() + // K_f2
                    M_c.val() * m_lambda * ilapl(w,G) * ilapl(w,G).tr(); // K_laplacian
    m_assembler.assemble(jacobian * meas(G));

    // ASSEMBLE NITSCHE
    if (m_options.getSwitch("AssembleWeakBCs"))
    {
        if (m_bcs.get("Weak Clamped").size()==0)
            gsWarn<<"Nitsche boundary assembly is requested, but no boundaries are marked 'Weak Clamped'";

        // Determine maximum mesh size
        // gsWarn<<"Mesh size computation needs to be checked!\n";
        T hmax = 0;
        for (size_t p=0; p!=m_patches.nPatches(); p++)
            hmax = math::max(hmax, m_patches.basis(p).getMaxCellLength());


        m_assembler.assembleBdr(m_bcs.get("Weak Clamped"), - m_lambda * igrad(w,G) *  nv(G)  * ilapl(w,G).tr() + // consistency term
                  m_penalty * (igrad(w,G) * nv(G).normalized()) * hmax * (igrad(w,G) * nv(G)).tr() - // penalty (stabilizing) term
                  m_lambda * ilapl(w,G) * (igrad(w,G)  * nv(G)).tr()); // symmetry term
    }

    // NOTE: no m_assembler.cleanUp() here — parse() handles cleanup at the start of the next call
}

template <class T>
void gsCahnHilliardAssembler<T>::assembleNitscheVector(const gsFunctionSet<T> & C)
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    if (m_bcs.get("Weak Clamped").size()==0)
    {
        gsWarn<<"Nitsche boundary assembly is requested, but no boundaries are marked 'Weak Clamped'";
        return;
    }

    // m_assembler.clearRhs(); // Resets to zero the values of the already allocated to residual (RHS)
    m_assembler.initVector();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Get the solution and its derivative
    auto c  = m_assembler.getCoeff(C);
    // auto dc = m_assembler.getCoeff(DC);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    // Determine maximum mesh size
    // gsWarn<<"Mesh size computation needs to be checked!\n";
    T hmax = 0;
    for (size_t p=0; p!=m_patches.nPatches(); p++)
        hmax = math::max(hmax, m_patches.basis(p).getMaxCellLength());

    m_assembler.assembleBdr(m_bcs.get("Weak Clamped"), - igrad(w,G) * nv(G) * m_lambda * ilapl(c,G).val() // consistency term
                                                       +(igrad(w,G) * nv(G).normalized()) * hmax * m_penalty * (igrad(c,G) * nv(G)) // penalty term
                                                       - m_lambda * ilapl(w,G) * igrad(c,G) * nv(G)); // symmetry term

    // NOTE: no m_assembler.cleanUp() here — parse() handles cleanup at the start of the next call
}

template <class T>
void gsCahnHilliardAssembler<T>::assembleNitscheMatrix()
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearMatrix();
    if (m_bcs.get("Weak Clamped").size()==0)
    {
        gsWarn<<"Nitsche boundary assembly is requested, but no boundaries are marked 'Weak Clamped'";
        m_assembler.makeMatrix();
        return;
    }


    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    // Determine maximum mesh size
    // gsWarn<<"Mesh size computation needs to be checked!\n";
    T hmax = 0;
    for (size_t p=0; p!=m_patches.nPatches(); p++)
        hmax = math::max(hmax, m_patches.basis(p).getMaxCellLength());

    m_assembler.assembleBdr(m_bcs.get("Weak Clamped"), - m_lambda * igrad(w,G) *  nv(G)  * ilapl(w,G).tr() + // consistency term
                  m_penalty * (igrad(w,G) * nv(G).normalized()) * hmax * (igrad(w,G) * nv(G)).tr() - // penalty (stabilizing) term
                  m_lambda * ilapl(w,G) * (igrad(w,G)  * nv(G)).tr()); // symmetry term
}

template <class T>
void gsCahnHilliardAssembler<T>::constructSolution(const gsMatrix<T> & Cvec,
                                                   gsMultiPatch<T>   & C) const
{
    auto w  = m_assembler.trialSpace(0);
    w.fixedPart() = m_ddofs;
    C.clear();
    for (size_t p = 0; p != m_patches.nPatches(); ++p)
    {
        gsMatrix<T> cf;
        w.getCoeffs(Cvec, cf, p);
        typename gsGeometry<T>::uPtr patch = m_spaceBasis->basis(p).makeGeometry(give(cf));
        C.addPatch(give(patch));
    }
}

template <class T>
void gsCahnHilliardAssembler<T>::constructSolution(const gsMatrix<T>   & Cvec,
                                                   gsMappedSpline<2,T> & C) const
{
    auto w  = m_assembler.trialSpace(0);
    w.fixedPart() = m_ddofs;
    gsMatrix<T> Cvec_nc = Cvec; // getSolution needs non-const ref
    auto c  = m_assembler.getSolution(w, Cvec_nc);
    c.extract(C);
}

template <class T>
void gsCahnHilliardAssembler<T>::constructSolution(const gsMultiPatch<T> & C,
                                                         gsMatrix<T>     & Cvec) const
{
    GISMO_ASSERT(C.geoDim()==1,"C must be a scalar function");
    const gsMultiBasis<T> * basis = dynamic_cast<const gsMultiBasis<T>*>(m_spaceBasis);
    GISMO_ASSERT(basis,"The space basis must be a multi-basis to construct the solution from a multi-patch");
    GISMO_ASSERT(C.nPatches()==basis->nBases(),"Number of patches in C must be equal to the number of bases in the assembler");
    auto w  = m_assembler.trialSpace(0);

    Cvec.setZero(this->numDofs(),1);
    for (size_t b=0; b!=basis->nBases(); b++)
    {
        for (index_t i = 0; i < basis->basis(b).size(); i++)
        {
            GISMO_ASSERT(C.basis(b).size()==basis->basis(b).size(),"Number of basis functions in C must be equal to the number of basis functions in the assembler");
            if (w.mapper().is_free(i,b))
                Cvec(w.mapper().index(i,b)) = C.patch(b).coefs()(i,0);
        }
    }
}

template <class T>
void gsCahnHilliardAssembler<T>::assemble(const gsMatrix<T> & Cvec)
{
    gsMultiPatch<T> C, DC;
    constructSolution(Cvec,  C);
    assemble(C);
}

template <class T>
void gsCahnHilliardAssembler<T>::assembleResidual(const gsMatrix<T> & Cvec)
{
    gsMultiPatch<T> C;
    constructSolution(Cvec,  C);
    assembleResidual(C);
}

template <class T>
void gsCahnHilliardAssembler<T>::assembleJacobian(const gsMatrix<T> & Cvec)
{
    gsMultiPatch<T> C;
    constructSolution(Cvec,  C);
    assembleJacobian(C);
}

template <class T>
void gsCahnHilliardAssembler<T>::assembleNitscheVector(const gsMatrix<T> & Cvec)
{
    gsMultiPatch<T> C;
    constructSolution(Cvec,  C);
    assembleNitscheVector(C);
}

template <class T>
void gsCahnHilliardAssembler<T>::setSpaceBasis(const gsFunctionSet<T> & spaceBasis)
{
    m_spaceBasis = &spaceBasis;
}

template <class T>
void gsCahnHilliardAssembler<T>::setIntegrationBasis(const gsMultiBasis<T> & integrationBasis)
{
    m_integrationBasis = &integrationBasis;
}

template <class T>
T gsCahnHilliardAssembler<T>::computeMass(const gsFunctionSet<T> & C)
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before computing mass.");
    gsExprEvaluator<T> evaluator(m_assembler);

    // Set the geometry map
    geometryMap G = evaluator.getMap(m_patches);
    
    // Get the solution
    auto c = evaluator.getVariable(C);
    
    // Compute M(c) = ∫ c dx
    return evaluator.integral(c * meas(G));
}

template <class T>
T gsCahnHilliardAssembler<T>::computeMass(const gsMatrix<T> & Cvec)
{
    gsMultiPatch<T> C;
    constructSolution(Cvec,  C);
    return computeMass(C);
}

template <class T>
T gsCahnHilliardAssembler<T>::computeDissipation(const gsFunctionSet<T> & C, const gsFunctionSet<T> & mu)
{
    // NOTE: this implementation requires a field for the chemical potential mu, avoiding the requirement for igrad(ilapl(c)) which requires a C^2 space.
    // NOTE2: C is not really used now
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet.");
    gsExprEvaluator<T> evaluator(m_assembler);
    geometryMap G = evaluator.getMap(m_patches);
    auto c = evaluator.getVariable(C, G);
    auto Mu= evaluator.getVariable(mu, G);

    // 1. Define Mobility (M)
    auto M_c = 1.0 * c.val(); // replace with const_expr(1.0) instead of using 0*c
    // 2. Define grad(mu)
    auto grad_mu = igrad(Mu, G);

    // Dissipation is integral of M * |grad mu|^2
    return evaluator.integral(M_c * grad_mu.sqNorm() * meas(G));
}

template <class T>
T gsCahnHilliardAssembler<T>::computeDissipation(const gsMatrix<T> & Cvec, const gsMatrix<T> & muVec)
{
    gsMultiPatch<T> C, mu;
    constructSolution(Cvec,  C);
    constructSolution(muVec, mu);
    return computeDissipation(C, mu);
}

template <class T>
T gsCahnHilliardAssembler<T>::computeDissipation(const gsFunctionSet<T> & C)
{
    // NOTE: Below is an implementation that requires igrad(ilapl(c)) which requires a C^2 space. 
    // igrad(ilapl) is not currently supported by the assembler, but it can be implemented in the future if needed. 
    /* 
        GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet.");
        gsExprEvaluator<T> evaluator(m_assembler);
        geometryMap G = evaluator.getMap(m_patches);
        auto c = evaluator.getVariable(C, G);

        // 1. Define Mobility (M)
        auto M_c = 1.0 * c.val(); // replace with const_expr(1.0) instead of using 0*c

        // 2. Define grad(mu) = Psi''(c)*grad(c) - lambda*grad(lapl(c))
        // Psi(c) = 1/4 (c^2 - 1)^2 
        // Psi'(c) = c^3 - c
        // Psi''(c) = 3*c^2 - 1
        // mu(c) = Psi'(c) - lambda*lapl(c) = c^3 - c - lambda*lapl(c)
        // grad(mu) = Psi''(c)*grad(c) - lambda*grad(lapl(c)) = (3*c^2 - 1)*grad(c) - lambda*grad(lapl(c))
        auto psi_double_prime = 3.0 * (c * c).val() - 1.0;
        auto grad_mu = psi_double_prime * igrad(c, G) - m_lambda * igrad(ilapl(c, G), G);
        
        // Dissipation is integral of M * |grad mu|^2
        return evaluator.integral(M_c * (grad_mu * grad_mu).sum() * meas(G));
    */
    GISMO_UNUSED(C);
    GISMO_NO_IMPLEMENTATION;
}

template <class T>
T gsCahnHilliardAssembler<T>::computeDissipation(const gsMatrix<T> & Cvec)
{
    gsMultiPatch<T> C;
    constructSolution(Cvec,  C);
    return computeDissipation(C);
}

template <class T>
T gsCahnHilliardAssembler<T>::computeEnergy(const gsFunctionSet<T> & C)
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet.");
    gsExprEvaluator<T> evaluator(m_assembler);
    geometryMap G = evaluator.getMap(m_patches);
    auto c = evaluator.getVariable(C, G);
    
    // Bulk energy: 1/4 * (c^2 - 1)^2
    auto c2 = (c * c).val();
    auto bulk = 0.25 * (c2 - 1.0) * (c2 - 1.0);
    
    // Interface energy: lambda/2 * |grad c|^2
    auto grad_c = igrad(c, G);
    auto interface = 0.5 * m_lambda *  grad_c.sqNorm();
    
    return evaluator.integral((bulk + interface) * meas(G));
}

template <class T>
T gsCahnHilliardAssembler<T>::computeEnergy(const gsMatrix<T> & Cvec)
{
    gsMultiPatch<T> C;
    constructSolution(Cvec,  C);
    return computeEnergy(C);
}

}// namespace gismo
