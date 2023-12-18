/** @file gsThinShellAssembler.hpp

    @brief Provides linear and nonlinear assemblers for thin shells

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsAssembler/gsBiharmonicExprAssembler.h>
#include <gsMSplines/gsMappedBasis.h>
#include <gsMSplines/gsMappedSpline.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsCore/gsFunctionExpr.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsPiecewiseFunction.h>

#ifdef gsSpectra_ENABLED
#include <gsSpectra/gsSpectra.h>
#endif

namespace gismo
{

template <class T>
gsBiharmonicExprAssembler<T>::gsBiharmonicExprAssembler(const gsMultiPatch<T> & mp,
                                                        const gsMultiBasis<T> & mb,
                                                        const gsFunctionSet<T> & force,
                                                        const gsBoundaryConditions<T> & bcs
                                                        )
:
m_patches(mp),
m_basis(mb),
m_spaceBasis(&mb),
m_force(&force),
m_bcs(bcs)
{
    this->_defaultOptions();
    this->_getOptions();
    this->_initialize();
};

template <class T>
gsBiharmonicExprAssembler<T>& gsBiharmonicExprAssembler<T>::operator=( const gsBiharmonicExprAssembler& other )
{
    if (this!=&other)
    {
        m_penalty=other.m_penalty;
        m_lambda=other.m_lambda;
        m_second=other.m_second;
        m_continuity=other.m_continuity;

        m_patches=other.m_patches;
        m_basis=other.m_basis;
        m_spaceBasis=other.m_spaceBasis;
        m_bcs=other.m_bcs;
        m_force=other.m_force;

        m_options=other.m_options;

        // To do: make copy constructor for the gsExprAssembler
        m_assembler.setIntegrationElements(m_basis);
        m_assembler.setOptions(m_options);
    }
    return *this;
}

template <class T>
gsBiharmonicExprAssembler<T>& gsBiharmonicExprAssembler<T>::operator=( gsBiharmonicExprAssembler&& other )
{
    m_penalty=give(other.m_penalty);
    m_lambda=give(other.m_lambda);
    m_second=give(other.m_second);
    m_continuity=give(other.m_continuity);

    m_patches=give(other.m_patches);
    m_basis=give(other.m_basis);
    m_spaceBasis=give(other.m_spaceBasis);
    m_bcs=give(other.m_bcs);
    m_force=give(other.m_force);

    m_options=give(other.m_options);

    // To do: make copy constructor for the gsExprAssembler
    m_assembler.setIntegrationElements(m_basis);
    m_assembler.setOptions(m_options);
    return *this;
}

template <class T>
void gsBiharmonicExprAssembler<T>::_defaultOptions()
{
    m_options.addReal("PenaltyIfc","Parameter Nitsche's method",-1);
    m_options.addReal("Lambda","Parameter for BC projection",1e-5);
    m_options.addSwitch("Second","Assemble the second biharmonic equation",false);
    m_options.addInt("Continuity","Set the continuity for the space",-1);
    // Assembler options
    gsOptionList assemblerOptions = m_assembler.defaultOptions().wrapIntoGroup("ExprAssembler");
    m_options.update(assemblerOptions,gsOptionList::addIfUnknown);
}

template <class T>
void gsBiharmonicExprAssembler<T>::_getOptions()
{
    m_penalty = m_options.getReal("PenaltyIfc");
    m_lambda = m_options.getReal("Lambda");
    m_second = m_options.getSwitch("Second");
    m_continuity = m_options.getInt("Continuity");

    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));
}

template<class T>
void gsBiharmonicExprAssembler<T>::setOptions(gsOptionList & options)
{
    m_options.update(options,gsOptionList::ignoreIfUnknown);
    this->_getOptions();
    this->_initialize();
}

template <class T>
void gsBiharmonicExprAssembler<T>::_initialize()
{
    // Elements used for numerical integration
    m_assembler.setIntegrationElements(m_basis);
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    GISMO_ASSERT(m_bcs.hasGeoMap(),"No geometry map was assigned to the boundary conditions. Use bc.setGeoMap to assign one!");

    m_assembler.getSpace(*m_spaceBasis, 1, 0); // last argument is the space ID
}

template <class T>
void gsBiharmonicExprAssembler<T>::_setup(const expr::gsFeSpace<T> & u)
{
    // Setup the mapper
    gsDofMapper map;
    _setMapperForBiharmonic(m_bcs, *m_spaceBasis, map);

    // Setup the system
    u.setupMapper(map);
    _getDirichletNeumannValuesL2Projection(m_patches, m_basis, m_bcs, *m_spaceBasis, u);
}

template <class T>
void gsBiharmonicExprAssembler<T>::assembleMass()
{
    m_assembler.cleanUp();
    this->_getOptions();
    this->_initialize();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Set the discretization space
    auto u = m_assembler.trialSpace(0);

    // Setup the mapper
    this->_setup(u);

    // Initialize the system
    m_assembler.initSystem();

    // Compute the system matrix and right-hand side
    m_assembler.assemble(u * u.tr() * meas(G));
}

template <class T>
void gsBiharmonicExprAssembler<T>::assemble()
{
    m_assembler.cleanUp();
    this->_getOptions();
    this->_initialize();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Set the source term
    auto ff = m_assembler.getCoeff(*m_force, G);

    // Set the discretization space
    auto u = m_assembler.trialSpace(0);

    // Setup the mapper
    this->_setup(u);

    // Initialize the system
    m_assembler.initSystem();

    // Compute the system matrix and right-hand side
    m_assembler.assemble(ilapl(u, G) * ilapl(u, G).tr() * meas(G),u * ff * meas(G));

    // Enforce Laplace conditions to right-hand side
    auto g_L = m_assembler.getBdrFunction(G); // Set the laplace bdy value
    //auto g_L = A.getCoeff(laplace, G);
    m_assembler.assembleBdr(m_bcs.get("Laplace"), (igrad(u, G) * nv(G)) * g_L.tr() );

    // Optionally, assemble Nitsche
    if (m_continuity>0)
    {
        gsMatrix<T> mu_interfaces;
        if (m_penalty == -1)
            _computeStabilityParameter(m_patches, m_basis, mu_interfaces);

        index_t i = 0;
        for ( typename gsMultiPatch<T>::const_iiterator it = m_patches.iBegin(); it != m_patches.iEnd(); ++it, ++i)
        {
            T stab     = 4 * ( m_basis.maxCwiseDegree() + m_basis.dim() ) * ( m_basis.maxCwiseDegree() + 1 );
            T m_h      = m_basis.basis(0).getMinCellLength(); // m_basis.basis(0).getMinCellLength();
            T mu       = 2 * stab / m_h;
            T alpha = 1;

            //mu = penalty_init == -1.0 ? mu : penalty_init / m_h;
            if (m_penalty == -1)
                mu = mu_interfaces(i,0) / m_h;
            else
                mu = m_penalty / m_h;

            std::vector<boundaryInterface> iFace;
            iFace.push_back(*it);
            m_assembler.assembleIfc(iFace,
                    //B11
                          -alpha * 0.5 * igrad(u.left(), G) * nv(G.left()).normalized() *
                          (ilapl(u.left(), G)).tr() * nv(G.left()).norm(),
                          -alpha * 0.5 *
                          (igrad(u.left(), G) * nv(G.left()).normalized() * (ilapl(u.left(), G)).tr()).tr() *
                          nv(G.left()).norm(),
                    //B12
                          -alpha * 0.5 * igrad(u.left(), G.left()) * nv(G.left()).normalized() *
                          (ilapl(u.right(), G.right())).tr() * nv(G.left()).norm(),
                          -alpha * 0.5 * (igrad(u.left(), G.left()) * nv(G.left()).normalized() *
                                          (ilapl(u.right(), G.right())).tr()).tr() * nv(G.left()).norm(),
                    //B21
                          alpha * 0.5 * igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                          (ilapl(u.left(), G.left())).tr() * nv(G.left()).norm(),
                          alpha * 0.5 * (igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                                         (ilapl(u.left(), G.left())).tr()).tr() * nv(G.left()).norm(),
                    //B22
                          alpha * 0.5 * igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                          (ilapl(u.right(), G.right())).tr() * nv(G.left()).norm(),
                          alpha * 0.5 * (igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                                         (ilapl(u.right(), G.right())).tr()).tr() * nv(G.left()).norm(),

                    // E11
                          mu * igrad(u.left(), G.left()) * nv(G.left()).normalized() *
                          (igrad(u.left(), G.left()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //-E12
                          -mu * (igrad(u.left(), G.left()) * nv(G.left()).normalized()) *
                          (igrad(u.right(), G.right()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //-E21
                          -mu * (igrad(u.right(), G.right()) * nv(G.left()).normalized()) *
                          (igrad(u.left(), G.left()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    // E22
                          mu * igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                          (igrad(u.right(), G.right()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm()
            );

            if (m_continuity>1)
                gsDebug<<"Put here the terms for continuity 1\n";
        }
    }
}


template <class T>
void gsBiharmonicExprAssembler<T>::assembleRHS()
{
    m_assembler.cleanUp();
    this->_getOptions();
    this->_initialize();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Set the source term
    auto ff = m_assembler.getCoeff(*m_force, G);

    // Set the discretization space
    auto u = m_assembler.trialSpace(0);

    // Setup the mapper
    this->_setup(u);

    // Initialize the system
    m_assembler.initVector(1);

    // Compute the system matrix and right-hand side
    m_assembler.assemble(u * ff * meas(G));

    // Enforce Laplace conditions to right-hand side
    auto g_L = m_assembler.getBdrFunction(G); // Set the laplace bdy value
    //auto g_L = A.getCoeff(laplace, G);
    m_assembler.assembleBdr(m_bcs.get("Laplace"), (igrad(u, G) * nv(G)) * g_L.tr() );
}

template <class T>
void gsBiharmonicExprAssembler<T>::assembleLHS()
{
    m_assembler.cleanUp();
    this->_getOptions();
    this->_initialize();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Set the discretization space
    auto u = m_assembler.trialSpace(0);

    // Setup the mapper
    this->_setup(u);

    // Initialize the system
    m_assembler.initMatrix();

    // Compute the system matrix and right-hand side
    m_assembler.assemble(ilapl(u, G) * ilapl(u, G).tr() * meas(G));

    // Optionally, assemble Nitsche
    if (m_continuity>0)
    {
        gsMatrix<T> mu_interfaces;
        if (m_penalty == -1)
            _computeStabilityParameter(m_patches, m_basis, mu_interfaces);

        index_t i = 0;
        for ( typename gsMultiPatch<T>::const_iiterator it = m_patches.iBegin(); it != m_patches.iEnd(); ++it, ++i)
        {
            T stab     = 4 * ( m_basis.maxCwiseDegree() + m_basis.dim() ) * ( m_basis.maxCwiseDegree() + 1 );
            T m_h      = m_basis.basis(0).getMinCellLength(); // m_basis.basis(0).getMinCellLength();
            T mu       = 2 * stab / m_h;
            T alpha = 1;

            //mu = penalty_init == -1.0 ? mu : penalty_init / m_h;
            if (m_penalty == -1)
                mu = mu_interfaces(i,0) / m_h;
            else
                mu = m_penalty / m_h;

            std::vector<boundaryInterface> iFace;
            iFace.push_back(*it);
            m_assembler.assembleIfc(iFace,
                    //B11
                          -alpha * 0.5 * igrad(u.left(), G) * nv(G.left()).normalized() *
                          (ilapl(u.left(), G)).tr() * nv(G.left()).norm(),
                          -alpha * 0.5 *
                          (igrad(u.left(), G) * nv(G.left()).normalized() * (ilapl(u.left(), G)).tr()).tr() *
                          nv(G.left()).norm(),
                    //B12
                          -alpha * 0.5 * igrad(u.left(), G.left()) * nv(G.left()).normalized() *
                          (ilapl(u.right(), G.right())).tr() * nv(G.left()).norm(),
                          -alpha * 0.5 * (igrad(u.left(), G.left()) * nv(G.left()).normalized() *
                                          (ilapl(u.right(), G.right())).tr()).tr() * nv(G.left()).norm(),
                    //B21
                          alpha * 0.5 * igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                          (ilapl(u.left(), G.left())).tr() * nv(G.left()).norm(),
                          alpha * 0.5 * (igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                                         (ilapl(u.left(), G.left())).tr()).tr() * nv(G.left()).norm(),
                    //B22
                          alpha * 0.5 * igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                          (ilapl(u.right(), G.right())).tr() * nv(G.left()).norm(),
                          alpha * 0.5 * (igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                                         (ilapl(u.right(), G.right())).tr()).tr() * nv(G.left()).norm(),

                    // E11
                          mu * igrad(u.left(), G.left()) * nv(G.left()).normalized() *
                          (igrad(u.left(), G.left()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //-E12
                          -mu * (igrad(u.left(), G.left()) * nv(G.left()).normalized()) *
                          (igrad(u.right(), G.right()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //-E21
                          -mu * (igrad(u.right(), G.right()) * nv(G.left()).normalized()) *
                          (igrad(u.left(), G.left()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    // E22
                          mu * igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                          (igrad(u.right(), G.right()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm()
            );

            if (m_continuity>1)
                gsDebug<<"Put here the terms for continuity 1\n";
        }
    }
}

template <class T>
void gsBiharmonicExprAssembler<T>::constructSolution(gsMatrix<T> & solVector)
{
    auto u = m_assembler.trialSpace(0);
    auto u_sol = m_assembler.getSolution(u, solVector);

    if (const gsMappedBasis<2,T> * bb2 = dynamic_cast<const gsMappedBasis<2,T> *>(m_spaceBasis))
    {
        u_sol.extract(m_mspline);
    }
    else
    {
        u_sol.extract(m_sol);
    }

}

template <class T>
typename gsFunctionSet<T>::Ptr gsBiharmonicExprAssembler<T>::getSolution() const
{
    if (const gsMappedBasis<2,T> * bb2 = dynamic_cast<const gsMappedBasis<2,T> *>(m_spaceBasis))
    {
        return memory::make_shared_not_owned(static_cast<gsFunctionSet<T>*>(&m_mspline));
    }
    else
    {
        return memory::make_shared_not_owned(static_cast<gsFunctionSet<T>*>(&m_sol));
    }
}

template <class T>
T gsBiharmonicExprAssembler<T>::l2error(gsMatrix<T> & solVector, const gsFunctionSet<T> & exact)
{
    auto u = m_assembler.trialSpace(0);
    auto u_sol = m_assembler.getSolution(u, solVector);
    geometryMap G = m_assembler.getMap(m_patches);

    gsExprEvaluator<T> ev(m_assembler);
    auto u_ex = ev.getVariable(exact, G);

    T l2err= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
    return l2err;
}

template <class T>
T gsBiharmonicExprAssembler<T>::h1error(gsMatrix<T> & solVector, const gsFunctionSet<T> & exact)
{
    geometryMap G = m_assembler.getMap(m_patches);
    auto u = m_assembler.trialSpace(0);
    auto u_sol = m_assembler.getSolution(u, solVector);

    gsExprEvaluator<T> ev(m_assembler);
    auto u_ex = ev.getVariable(exact, G);

    T l2err = math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) ); // / ev.integral(f.sqNorm()*meas(G)) );
    T h1err = l2err +
            math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
    return h1err;
}

template <class T>
T gsBiharmonicExprAssembler<T>::h2error(gsMatrix<T> & solVector, const gsFunctionSet<T> & exact)
{
    geometryMap G = m_assembler.getMap(m_patches);
    auto u = m_assembler.trialSpace(0);
    auto u_sol = m_assembler.getSolution(u, solVector);

    gsExprEvaluator<T> ev(m_assembler);
    auto u_ex = ev.getVariable(exact, G);

    T l2err = math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) ); // / ev.integral(f.sqNorm()*meas(G)) );
    T h1err = l2err +
            math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
    T h2err = h1err +
            math::sqrt(ev.integral( ( ihess(u_ex) - ihess(u_sol,G) ).sqNorm() * meas(G) )); // /ev.integral( ihess(f).sqNorm()*meas(G) )
    return h2err;
}


template <class T>
std::tuple<T,T,T> gsBiharmonicExprAssembler<T>::errors(gsMatrix<T> & solVector, const gsFunctionSet<T> & exact)
{
    geometryMap G = m_assembler.getMap(m_patches);
    auto u = m_assembler.trialSpace(0);
    auto u_sol = m_assembler.getSolution(u, solVector);

    gsExprEvaluator<T> ev(m_assembler);
    auto u_ex = ev.getVariable(exact, G);

    T l2err = math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) ); // / ev.integral(f.sqNorm()*meas(G)) );
    T h1err = l2err +
            math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
    T h2err = h1err +
            math::sqrt(ev.integral( ( ihess(u_ex) - ihess(u_sol,G) ).sqNorm() * meas(G) )); // /ev.integral( ihess(f).sqNorm()*meas(G) )
    return std::make_tuple(l2err,h1err,h2err);
}

template <class T>
T gsBiharmonicExprAssembler<T>::interfaceError(gsMatrix<T> & solVector, const gsFunctionSet<T> & exact)
{
    if (const gsMappedBasis<2,T> * bb2 = dynamic_cast<const gsMappedBasis<2,T> *>(m_spaceBasis))
    {
        geometryMap G = m_assembler.getMap(m_patches);
        auto u = m_assembler.trialSpace(0);
        auto u_sol = m_assembler.getSolution(u, solVector);

        gsMatrix<T> solFull;
        u_sol.extractFull(solFull);
        gsMappedSpline<2, T> mappedSpline(*bb2, solFull);

        auto ms_sol = m_assembler.getCoeff(mappedSpline);

        gsExprEvaluator<T> ev(m_assembler);
        T IFaceErr = math::sqrt(ev.integralInterface(((igrad(ms_sol.left(), G.left()) -
                                                       igrad(ms_sol.right(), G.right())) *
                                                       nv(G).normalized()).sqNorm() * meas(G)));

        return IFaceErr;
    }
    else
    {
        geometryMap G = m_assembler.getMap(m_patches);
        auto u = m_assembler.trialSpace(0);
        auto u_sol = m_assembler.getSolution(u, solVector);

        gsMultiPatch<T> sol;
        u_sol.extract(sol);

        auto ms_sol = m_assembler.getCoeff(sol);

        gsExprEvaluator<T> ev(m_assembler);
        T IFaceErr = math::sqrt(ev.integralInterface(((igrad(ms_sol.left(), G.left()) -
                                                       igrad(ms_sol.right(), G.right())) *
                                                       nv(G).normalized()).sqNorm() * meas(G)));

        return IFaceErr;
    }
}

template <class T>
void gsBiharmonicExprAssembler<T>::setSpaceBasis(const gsFunctionSet<T> & spaceBasis)
{
    m_spaceBasis = &spaceBasis;
    this->_getOptions();
    this->_initialize();
}

template <class T>
void gsBiharmonicExprAssembler<T>::_setMapperForBiharmonic(  const gsBoundaryConditions<T> & bc,
                                                            const gsFunctionSet<T> & spaceBasis,
                                                            gsDofMapper & mapper)
{
    if (const gsMappedBasis<2,T> * bb2 = dynamic_cast<const gsMappedBasis<2,T> *>(&spaceBasis))
    {
        mapper.setIdentity(bb2->nPatches(), bb2->size(), 1);

        gsMatrix<index_t> bnd;
        for (typename gsBoundaryConditions<T>::const_iterator it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it)
        {
          bnd = bb2->basis(it->ps.patch).boundary(it->ps.side());
          mapper.markBoundary(it->ps.patch, bnd, 0);
        }

        for (typename gsBoundaryConditions<T>::const_iterator it = bc.begin("Neumann"); it != bc.end("Neumann"); ++it)
        {
          bnd = bb2->basis(it->ps.patch).boundaryOffset(it->ps.side(),1);
          mapper.markBoundary(it->ps.patch, bnd, 0);
        }
        mapper.finalize();
    }
    else if (const gsMultiBasis<T> * dbasis = dynamic_cast<const gsMultiBasis<T> *>(&spaceBasis))
    {
        mapper.init(*dbasis);

        for (gsBoxTopology::const_iiterator it = dbasis->topology().iBegin();
             it != dbasis->topology().iEnd(); ++it) // C^0 at the interface
        {
            dbasis->matchInterface(*it, mapper);
        }

        gsMatrix<index_t> bnd;
        for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it)
        {
            bnd = dbasis->basis(it->ps.patch).boundary(it->ps.side());
            mapper.markBoundary(it->ps.patch, bnd, 0);
        }

        for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Neumann"); it != bc.end("Neumann"); ++it)
        {
            bnd = dbasis->basis(it->ps.patch).boundaryOffset(it->ps.side(),1);
            mapper.markBoundary(it->ps.patch, bnd, 0);
        }
        mapper.finalize();
    }
    else
    {
        GISMO_ERROR("bb2 is not a gsMappedBasis");
    }
};


template <class T>
void gsBiharmonicExprAssembler<T>::_getDirichletNeumannValuesL2Projection(
                                                                        const gsMultiPatch<T> & mp,
                                                                        const gsMultiBasis<T> & dbasis,
                                                                        const gsBoundaryConditions<T> & bc,
                                                                        const gsFunctionSet<T> & spaceBasis,
                                                                        const expr::gsFeSpace<T> & u
                                                                        )
{
    gsDebugVar(bc.dirichletSides().size());
    gsDebugVar(bc.neumannSides().size());

    if (bc.dirichletSides().size()==0 && bc.neumannSides().size()==0)
        return;
    if (const gsMappedBasis<2,T> * bb2 = dynamic_cast<const gsMappedBasis<2,T> *>(&spaceBasis))
    {
        const gsDofMapper & mapper = u.mapper();

        gsMatrix<index_t> bnd = mapper.findFree(mapper.numPatches()-1);
        gsDofMapper mapperBdy;
        mapperBdy.setIdentity(bb2->nPatches(), bb2->size(), 1);  // bb2->nPatches() == 1
        mapperBdy.markBoundary(0, bnd, 0);
        mapperBdy.finalize();

        gsExprAssembler<T> A(1,1);
        A.setIntegrationElements(dbasis);

        auto G = A.getMap(mp);
        auto uu = A.getSpace(*bb2);
        auto g_bdy = A.getBdrFunction(G);

        uu.setupMapper(mapperBdy);
        gsMatrix<T> & fixedDofs_A = const_cast<expr::gsFeSpace<T>&>(uu).fixedPart();
        fixedDofs_A.setZero( uu.mapper().boundarySize(), 1 );

        A.initSystem();
        A.assembleBdr(bc.get("Dirichlet"), uu * uu.tr() * meas(G));
        A.assembleBdr(bc.get("Dirichlet"), uu * g_bdy * meas(G));
        A.assembleBdr(bc.get("Neumann"),m_lambda * (igrad(uu, G) * nv(G).normalized()) * (igrad(uu, G) * nv(G).normalized()).tr() * meas(G));
        A.assembleBdr(bc.get("Neumann"),m_lambda *  (igrad(uu, G) * nv(G).normalized()) * (g_bdy.tr() * nv(G).normalized()) * meas(G));

        typename gsSparseSolver<T>::SimplicialLDLT solver;
        solver.compute( A.matrix() );
        gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>& >(u).fixedPart();
        fixedDofs = solver.solve(A.rhs());
    }
    else if (dynamic_cast<const gsMultiBasis<T> *>(&spaceBasis)) // assumes spacebasis==dbasis
    {
        gsDofMapper mapper = u.mapper();
        gsDofMapper mapperBdy(dbasis, u.dim());
        for (gsBoxTopology::const_iiterator it = dbasis.topology().iBegin();
             it != dbasis.topology().iEnd(); ++it) // C^0 at the interface
        {
            dbasis.matchInterface(*it, mapperBdy);
        }
        for (index_t np = 0; np < mp.nPieces(); np++)
        {
            gsMatrix<index_t> bnd = mapper.findFree(np);
            mapperBdy.markBoundary(np, bnd, 0);
        }
        mapperBdy.finalize();

        gsExprAssembler<T> A(1,1);
        A.setIntegrationElements(dbasis);

        auto G = A.getMap(mp);
        auto uu = A.getSpace(dbasis);
        auto g_bdy = A.getBdrFunction(G);

        uu.setupMapper(mapperBdy);
        gsMatrix<T> & fixedDofs_A = const_cast<expr::gsFeSpace<T>&>(uu).fixedPart();
        fixedDofs_A.setZero( uu.mapper().boundarySize(), 1 );

        A.initSystem();
        A.assembleBdr(bc.get("Dirichlet"), uu * uu.tr() * meas(G));
        A.assembleBdr(bc.get("Dirichlet"), uu * g_bdy * meas(G));
        A.assembleBdr(bc.get("Neumann"),
                      m_lambda * (igrad(uu, G) * nv(G).normalized()) * (igrad(uu, G) * nv(G).normalized()).tr() * meas(G));
        A.assembleBdr(bc.get("Neumann"),
                      m_lambda *  (igrad(uu, G) * nv(G).normalized()) * (g_bdy.tr() * nv(G).normalized()) * meas(G));

        typename gsSparseSolver<T>::SimplicialLDLT solver;
        solver.compute( A.matrix() );
        gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>& >(u).fixedPart();
        gsMatrix<T> fixedDofs_temp = solver.solve(A.rhs());

        // Reordering the dofs of the boundary
        fixedDofs.setZero(mapper.boundarySize(),1);
        index_t sz = 0;
        for (index_t np = 0; np < mp.nPieces(); np++)
        {
            gsMatrix<index_t> bnd = mapperBdy.findFree(np);
            bnd.array() += sz;
            for (index_t i = 0; i < bnd.rows(); i++)
            {
                index_t ii = mapperBdy.asVector()(bnd(i,0));
                fixedDofs(mapper.global_to_bindex(mapper.asVector()(bnd(i,0))),0) = fixedDofs_temp(ii,0);
            }
            sz += mapperBdy.patchSize(np,0);
        }
    }
    else
    {
        GISMO_ERROR("bb2 is not a gsMappedBasis");
    }
};

template <class T>
void gsBiharmonicExprAssembler<T>::_computeStabilityParameter(
                                                                const gsMultiPatch<T> & mp,
                                                                const gsMultiBasis<T> & dbasis,
                                                                gsMatrix<T> & mu_interfaces)
{
    mu_interfaces.setZero(mp.nInterfaces(),1);

    index_t i = 0;
    for ( typename gsMultiPatch<T>::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it, ++i)
    {
        gsMultiPatch<T> mp_temp;
        mp_temp.addPatch(mp.patch(it->first().patch));
        mp_temp.addPatch(mp.patch(it->second().patch));
        mp_temp.computeTopology();

        gsMultiBasis<T> dbasis_temp;
        dbasis_temp.addBasis(dbasis.basis(it->first().patch).clone().release());
        dbasis_temp.addBasis(dbasis.basis(it->second().patch).clone().release());

        gsBoundaryConditions<T> bc;

//        patchSide pS1 = mp_temp.interfaces()[0].first();
//        patchSide pS2 = mp_temp.interfaces()[0].second();
//
//
//        index_t side = pS1.index() < 3 ? (pS1.index() == 1 ? 2 : 1) : (pS1.index() == 3 ? 4 : 3);
//        bc.addCondition(patchSide(pS1.patchIndex(), side), condition_type::dirichlet, 0);
//
//        side = pS2.index() < 3 ? (pS2.index() == 1 ? 2 : 1) : (pS2.index() == 3 ? 4 : 3);
//        bc.addCondition(patchSide(pS2.patchIndex(), side), condition_type::dirichlet, 0);

        // Make the Eigenvalue problem to a homogeneous one
        for (typename gsMultiPatch<T>::const_biterator bit = mp_temp.bBegin(); bit != mp_temp.bEnd(); ++bit)
            bc.addCondition(*bit, condition_type::dirichlet, 0);

        gsExprAssembler<T> A2(1, 1), B2(1, 1);

        // Elements used for numerical integration
        A2.setIntegrationElements(dbasis_temp);
        B2.setIntegrationElements(dbasis_temp);

        // Set the geometry map
        auto GA = A2.getMap(mp_temp);
        auto GB = B2.getMap(mp_temp);

        // Set the discretization space
        auto uA = A2.getSpace(dbasis_temp);
        auto uB = B2.getSpace(dbasis_temp);

        uA.setup(bc, dirichlet::homogeneous, 0);
        uB.setup(bc, dirichlet::homogeneous,0);
        //uA.setup(0);
        //uB.setup(0);

        A2.initSystem();
        B2.initSystem();

        T c = 0.25;
        A2.assembleIfc(mp_temp.interfaces(),
                       c * ilapl(uA.left(), GA.left()) * ilapl(uA.left(), GA.left()).tr() * nv(GA.left()).norm(),
                       c * ilapl(uA.left(), GA.left()) * ilapl(uA.right(), GA.right()).tr() * nv(GA.left()).norm(),
                       c * ilapl(uA.right(), GA.right()) * ilapl(uA.left(), GA.left()).tr() * nv(GA.left()).norm(),
                       c * ilapl(uA.right(), GA.right()) * ilapl(uA.right(), GA.right()).tr() * nv(GA.left()).norm());

        B2.assemble(ilapl(uB, GB) * ilapl(uB, GB).tr() * meas(GB));

        // TODO INSTABLE && SLOW
        gsMatrix<T> AA = A2.matrix().toDense();
        gsMatrix<T> BB = B2.matrix().toDense();

#ifdef gsSpectra_ENABLED
        gsSpectraGenSymSolver<gsSparseMatrix<T>,Spectra::GEigsMode::Cholesky> ges(A2.matrix(),B2.matrix(),1,10);
        ges.compute(Spectra::SortRule::LargestMagn,1000,1e-6,Spectra::SortRule::LargestMagn);
        T maxEV = ges.eigenvalues()(0,0);
#else
        gsEigen::GeneralizedSelfAdjointEigenSolver<gsMatrix<T>> ges(AA, BB);
        T maxEV = ges.eigenvalues().array().maxCoeff();
#endif
        T m_h      = dbasis_temp.basis(0).getMinCellLength(); // dbasis.basis(0).getMinCellLength();
        mu_interfaces(i,0) = 16.0 * m_h * maxEV;
/*
        gsSparseSolver<T>::SimplicialLDLT sol;
        sol.compute(B2.matrix());
        gsSparseMatrix<T> R = sol.matrixU();
        gsSparseMatrix<T> RT = sol.matrixL();
        gsMatrix<T> AAA = RT.toDense().inverse() * AA * R.toDense().inverse();

        gsConjugateGradient<T> cg(AAA);

        cg.setCalcEigenvalues(true);
        cg.setTolerance(1e-15);
        cg.setMaxIterations(100000);

        gsMatrix<T> rhs, result;
        rhs.setRandom( AAA.rows(), 1 );
        result.setRandom( AAA.rows(), 1 );

        cg.solve(rhs,result);

        gsInfo << "Tol: " << cg.error() << "\n";
        gsInfo << "Max it: " << cg.iterations() << "\n";

        gsMatrix<T> eigenvalues;
        cg.getEigenvalues(eigenvalues);

        gsInfo << "Cond Number: " << eigenvalues.bottomRows(1)(0,0) << "\n";
*/
    }
}


}// namespace gismo
