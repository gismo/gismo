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
m_basis(mb),
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
        m_basis=other.m_basis;
        m_spaceBasis=other.m_spaceBasis;
        m_bcs=other.m_bcs;

        m_options=other.m_options;

        // To do: make copy constructor for the gsExprAssembler
        m_assembler.setIntegrationElements(m_basis);
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
    m_basis=give(other.m_basis);
    m_spaceBasis=give(other.m_spaceBasis);
    m_bcs=give(other.m_bcs);

    m_options=give(other.m_options);

    // To do: make copy constructor for the gsExprAssembler
    m_assembler.setIntegrationElements(m_basis);
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
    m_assembler.setIntegrationElements(m_basis);
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
void gsCahnHilliardAssembler<T>::assembleResidual(const gsMultiPatch<T> & C, const gsMultiPatch<T> & DC)
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    // m_assembler.clearRhs(); // Resets to zero the values of the already allocated to residual (RHS)
    m_assembler.initVector();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Get the solution and its derivative
    auto c  = m_assembler.getCoeff(C);
    auto dc = m_assembler.getCoeff(DC);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    // Derivatives of the double well potential (Gomez et al., 2008)
    // @lventavinuela where to use the double well potential option?
    auto dmu_c = - 1.0 + 3.0 * (c*c).val(); // f_2 (second derivative of double well)
    // auto ddmu_c = 6*c.val(); // f_3 (third derivative of double well)

    // // Mobility
    // T m0 = m_options.getReal("M0");
    auto M_c  = 1.0 + 0.0*c.val(); // replace with const_expr(1.0) instead of using 0*c
    // auto dM_c = 0.0 * igrad(c,G); // replace with const_expr(1.0) instead of using 0*c!!

    auto residual = w*dc + // M
                    M_c.val() * igrad(w,G)  * dmu_c * igrad(c,G).tr() + // F_bar
                    M_c.val() * ilapl(w,G)*m_lambda*ilapl(c,G).val(); // K_laplacian
                    // lambda*ilapl(c,G).val()*igrad(w,G)*dM_c.tr() + // term gradient mobility!
    m_assembler.assemble(residual * meas(G));

    // ASSEMBLE NITSCHE
    if (m_options.getSwitch("AssembleWeakBCs"))
    {
        if (m_bcs.get("Weak Clamped").size()==0)
            gsWarn<<"Nitsche boundary assembly is requested, but no boundaries are marked 'Weak Clamped'";
        // Determine maximum mesh size
        gsWarn<<"Mesh size computation needs to be checked!\n";
        T hmax = 0;
        for (size_t p=0; p!=m_patches.nPatches(); p++)
            hmax = math::max(hmax, m_patches.basis(p).getMaxCellLength());

        m_assembler.assembleBdr(m_bcs.get("Weak Clamped"), - igrad(w,G) * nv(G) * m_lambda * ilapl(c,G).val() // consistency term
                                                        +(igrad(w,G) * nv(G).normalized()) * hmax * m_penalty * (igrad(c,G) * nv(G)) // penalty term
                                                        - m_lambda * ilapl(w,G) * igrad(c,G) * nv(G)); // symmetry term
    }

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
void gsCahnHilliardAssembler<T>::assembleJacobian(const gsMultiPatch<T> & C, const gsMultiPatch<T> & DC)
{
    GISMO_UNUSED(DC);

    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Get the solution and its derivative
    auto c  = m_assembler.getCoeff(C);
    // auto dc = m_assembler.getCoeff(DC);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    // Derivatives of the double well potential (Gomez et al., 2008)
    auto dmu_c = - 1.0 + 3.0 * (c*c).val(); // f_2 (second derivative of double well)
    auto ddmu_c = 6*c.val(); // f_3 (third derivative of double well)

    // // Mobility
    // auto M_c  = 1.0 + 0.0*c.val(); // replace with const_expr(1.0) instead of using 0*c
    // auto dM_c = 0.0 * igrad(c,G); // replace with const_expr(1.0) instead of using 0*c!!

    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)

    m_assembler.assemble((
                dmu_c *igrad(w,G) * igrad(w,G).tr() + // K_f1
                ddmu_c * igrad(w,G) * igrad(c,G).tr() * w.tr() + // K_f2
                m_lambda * ilapl(w,G) * ilapl(w,G).tr()
                ) * meas(G)); // K_laplacian

    // ASSEMBLE NITSCHE
    if (m_options.getSwitch("AssembleWeakBCs"))
    {
        if (m_bcs.get("Weak Clamped").size()==0)
            gsWarn<<"Nitsche boundary assembly is requested, but no boundaries are marked 'Weak Clamped'";

        // Determine maximum mesh size
        gsWarn<<"Mesh size computation needs to be checked!\n";
        T hmax = 0;
        for (size_t p=0; p!=m_patches.nPatches(); p++)
            hmax = math::max(hmax, m_patches.basis(p).getMaxCellLength());


        m_assembler.assembleBdr(m_bcs.get("Weak Clamped"), - m_lambda * igrad(w,G) *  nv(G)  * ilapl(w,G).tr() + // consistency term
                  m_penalty * (igrad(w,G) * nv(G).normalized()) * hmax * (igrad(w,G) * nv(G)).tr() - // penalty (stabilizing) term
                  m_lambda * ilapl(w,G) * (igrad(w,G)  * nv(G)).tr()); // symmetry term
    }
}

template <class T>
void gsCahnHilliardAssembler<T>::assembleNitscheVector(const gsMultiPatch<T> & C, const gsMultiPatch<T> & DC)
{
    GISMO_UNUSED(DC);

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
    gsWarn<<"Mesh size computation needs to be checked!\n";
    T hmax = 0;
    for (size_t p=0; p!=m_patches.nPatches(); p++)
        hmax = math::max(hmax, m_patches.basis(p).getMaxCellLength());

    m_assembler.assembleBdr(m_bcs.get("Weak Clamped"), - igrad(w,G) * nv(G) * m_lambda * ilapl(c,G).val() // consistency term
                                                       +(igrad(w,G) * nv(G).normalized()) * hmax * m_penalty * (igrad(c,G) * nv(G)) // penalty term
                                                       - m_lambda * ilapl(w,G) * igrad(c,G) * nv(G)); // symmetry term
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
    gsWarn<<"Mesh size computation needs to be checked!\n";
    T hmax = 0;
    for (size_t p=0; p!=m_patches.nPatches(); p++)
        hmax = math::max(hmax, m_patches.basis(p).getMaxCellLength());

    m_assembler.assembleBdr(m_bcs.get("Weak Clamped"), - m_lambda * igrad(w,G) *  nv(G)  * ilapl(w,G).tr() + // consistency term
                  m_penalty * (igrad(w,G) * nv(G).normalized()) * hmax * (igrad(w,G) * nv(G)).tr() - // penalty (stabilizing) term
                  m_lambda * ilapl(w,G) * (igrad(w,G)  * nv(G)).tr()); // symmetry term
}

template <class T>
void gsCahnHilliardAssembler<T>::constructSolution(gsMatrix<T>     & Cvec,
                                                   gsMultiPatch<T> & C)
{
    auto w  = m_assembler.trialSpace(0);
    auto c  = m_assembler.getSolution(w,  Cvec);
    c.extract(C);
}

template <class T>
void gsCahnHilliardAssembler<T>::constructSolution(gsMatrix<T>         & Cvec,
                                                   gsMappedSpline<2,T> & C)
{
    auto w  = m_assembler.trialSpace(0);
    auto c  = m_assembler.getSolution(w,  Cvec);
    c.extract(C);
}

template <class T>
void gsCahnHilliardAssembler<T>::constructSolution(const gsMultiPatch<T> & C,
                                                         gsMatrix<T>     & Cvec)
{
    GISMO_ASSERT(C.geoDim()==1,"C must be a scalar function");
    GISMO_ASSERT(C.nPatches()==m_basis.nBases(),"Number of patches in C must be equal to the number of bases in the assembler");
    auto w  = m_assembler.trialSpace(0);

    Cvec.setZero(this->numDofs(),1);
    for (size_t b=0; b!=m_basis.nBases(); b++)
    {
        for (index_t i = 0; i < m_basis.basis(b).size(); i++)
        {
            GISMO_ASSERT(C.basis(b).size()==m_basis.basis(b).size(),"Number of basis functions in C must be equal to the number of basis functions in the assembler");
            if (w.mapper().is_free(i,b))
                Cvec(w.mapper().index(i,b)) = C.patch(b).coefs()(i,0);
        }
    }
}

template <class T>
void gsCahnHilliardAssembler<T>::setSpaceBasis(const gsFunctionSet<T> & spaceBasis)
{
    m_spaceBasis = &spaceBasis;
}

}// namespace gismo
