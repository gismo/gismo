/** @file gsG1BiharmonicAssembler.h
    @brief Provides assembler for a homogenius Biharmonic equation.
    This file is part of the G+Smo library.
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    Author(s): P. Weinmueller
*/

#pragma once

#include <gsAssembler/gsAssembler.h>

#include <gsPde/gsBiharmonicPde.h>

#include <gsAssembler/gsVisitorBiharmonic.h>
#include <gsAssembler/gsVisitorNeumann.h>
#include <gsAssembler/gsVisitorNeumannBiharmonic.h>
//#include <gsAssembler/gsVisitorNitscheBiharmonic.h>

namespace gismo
{

/** @brief
    Implementation of a homogeneous Biharmonic Assembler.
    It sets up an assembler and assembles the system patch wise and
    combines the patch-local stiffness matrices into a global system.
    Dirichlet boundary can only be enforced strongly (i.e Nitsche is
    not implemented).
*/
template <class T, class bhVisitor = gsVisitorBiharmonic<T> >
class gsG1BiharmonicAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
/** @brief
    Constructor of the assembler object.
    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] bconditions  is a gsBoundaryConditions object that holds boundary conditions on the form:
    \f[ \text{Dirichlet: } u = g \text{ on } \Gamma, \text{ and Neumann: } \nabla \Delta u \cdot \mathbf{n} = h \text{ on } \Gamma\f]
    \param[in] bconditions2 is a gsBoundaryConditions object that holds Neumann boundary conditions on the form:
    \f[\text{Neumann: } \nabla \Delta u \cdot \mathbf{n} = g\, \rightarrow \,(g,\nabla v \cdot \mathbf{n})_\Gamma, \f] where \f$ g \f$ is the Neumann data,
    \f$ v \f$ is the test function and \f$ \Gamma\f$ is the boundary side.
    \param[in] rhs is the right-hand side of the Biharmonic equation, \f$\mathbf{f}\f$.
    \param[in] dirStrategy option for the treatment of Dirichlet boundary in the \em bconditions object.
    \param[in] intStrategy option for the treatment of patch interfaces
*/
    gsG1BiharmonicAssembler(gsMultiPatch<T> const        & patches,
                            gsMultiBasis<T> const         & bases,
                            gsBoundaryConditions<T> const & bconditions,
                            gsBoundaryConditions<T> const & bconditions2,
                            const gsPiecewiseFunction<T>          & rhs)
        : m_ppde(patches,bconditions,bconditions2,rhs)
    {
        iFace::strategy intStrategy = iFace::none;
        dirichlet::strategy dirStrategy = dirichlet::none;

        m_options.setInt("DirichletStrategy", dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);

        Base::initialize(m_ppde, bases, m_options);

        const gsMultiBasis<T> & mbasis = m_bases[m_system.colBasis(0)];
        const gsDofMapper & mapper =
            dirichlet::elimination == m_options.getInt("DirichletStrategy") ?
            m_system.colMapper(0) :
            mbasis.getMapper(dirichlet::elimination,
                             static_cast<iFace::strategy>(m_options.getInt("InterfaceStrategy")),
                             m_pde_ptr->bc(), 0);
        map_boundary = mapper;

        // For interface:
        // 1. Obtain a map from basis functions to matrix columns and rows
        gsVector<index_t> sz(m_ppde.patches().nPatches()); // n Patches
        for (index_t i = 0; i < sz.size(); i++)
            sz[i] = m_bases[0].basis(i).size();
        gsDofMapper mapper_interface(sz);

        gsMatrix<unsigned> act1, act2;
        for (unsigned i = 0 ; i < m_ppde.patches().interfaces().size(); i++)  // Iterate over the interfaces
        {
            const boundaryInterface & iFace = m_ppde.patches().interfaces()[i];// assume one single interface

            gsBasis<> & B1 = m_bases[0].basis(iFace.first().patch);
            gsBasis<> & B2 = m_bases[0].basis(iFace.second().patch);

            // glue interface
            //B1.matchWith(iFace, B2, act1, act2);

            // mark dofs
            act1 = B1.boundaryOffset(iFace.first().side(), 0);
            mapper_interface.markBoundary(iFace.first().patch, act1);//interface
            act2 = B2.boundaryOffset(iFace.second().side(), 0);
            mapper_interface.markBoundary(iFace.second().patch, act2);//interface

            act1 = B1.boundaryOffset(iFace.first().side(), 1);
            mapper_interface.markBoundary(iFace.first().patch, act1); //first adj. face
            act2 = B2.boundaryOffset(iFace.second().side(), 1);
            mapper_interface.markBoundary(iFace.second().patch, act2);//second adj. face
        }
        mapper_interface.finalize();
        //mapper_interface.print();
        map_interface = mapper_interface;

        // END
    }

    void refresh();

    void assemble();

    void constructSolution(const gsMatrix<T>& solVector,
                           gsMultiPatch<T>& result, int unk = 0);

    void writeParaview(const gsField<T> & field,
                       std::string const & fn, gsMultiPatch<T> mp_L,
                       gsMultiPatch<T> mp_R,
                       unsigned npts = 1000, bool mesh = false);

    void writeParaview(const gsField<T> & field,
                       std::string const & fn,
                       std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>> mp_g1,
                       unsigned npts = 1000, bool mesh = false);


    void computeDirichletDofsL2Proj(gsMultiPatch<T> & basisG1_L,
                                    gsMultiPatch<T> & basisG1_R,
                                    index_t n_tilde, index_t n_bar,
                                    size_t unk = 0);

    void computeDirichletDofsL2Proj(std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>> & basisG1,
                                    std::vector<index_t> n_tilde, std::vector<index_t> n_bar,
                                    size_t unk = 0);

    gsDofMapper get_mapper_boundary() { return map_boundary; };
    gsDofMapper get_mapper_interface() { return map_interface; };
    gsDofMapper get_mapper() { return m_system.colMapper(0); };

    gsMatrix<T> get_g1dofs() { return m_g1_ddof; };


protected:

    // fixme: add constructor and remove this
    gsBiharmonicPde<T> m_ppde;

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

protected:

    gsMatrix<T> m_g1_ddof;

    gsDofMapper map_boundary, map_interface;

};

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::constructSolution(const gsMatrix<T> & solVector,
                                                             gsMultiPatch<T> & result,
                                                             int unk)
{
    const gsDofMapper & mapper = m_system.colMapper(unk); // unknown = 0
    result.clear();

    const index_t dim = ( 0!=solVector.cols() ? solVector.cols() :  m_ddof[unk].cols() );

    gsMatrix<T> coeffs;
    for (size_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p)
    {
        // Reconstruct solution coefficients on patch p
        const int sz  = m_bases[m_system.colBasis(unk)][p].size();
        coeffs.resize(sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if (map_boundary.is_free(i, p)) // DoF value is in the solVector
            {
                //gsInfo << "mapper index: " << i << " : " << mapper.index(i, p) << " : " << solVector.row( mapper.index(i, p) ) << "\n";
                coeffs.row(i) = solVector.row( mapper.index(i, p) );
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                if (map_interface.is_free(i,p))
                {
                    //gsInfo << "mapper index dirichlet: " << i << " : " << map_boundary.bindex(i, p) << " : " << m_g1_ddof.row( map_boundary.bindex(i, p) ) << "\n";
                    //coeffs.row(i) = m_ddof[unk].row( mapper.bindex(i, p) ).head(dim);
                    coeffs.row(i) = m_g1_ddof.row( map_boundary.bindex(i, p) );
                }
                else // For interface boundary
                {
                    //gsInfo << "mapper index dirichlet: " << i << " : " << mapper.index(i, p) << " : " << solVector.row( mapper.index(i, p) ) << "\n";
                    coeffs.row(i) = solVector.row( mapper.index(i, p) ); // == 0
                }
            }
        }

        result.addPatch( m_bases[m_system.colBasis(unk)][p].makeGeometry( give(coeffs) ) );
    }
}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::refresh()
{
    // We use predefined helper which initializes the system matrix
    // rows and columns using the same test and trial space
    Base::scalarProblemGalerkinRefresh();

}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::assemble()
{
    GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
    m_system.reserve(nz, this->pde().numRhs());

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    m_ddof.resize(m_system.numUnknowns());
    m_ddof[0].setZero(m_system.colMapper(0).boundarySize(), m_system.unkSize(0)*m_system.rhs().cols() );

    // Assemble volume integrals
    Base::template push<bhVisitor >();

    // Neuman conditions of first kind
    //Base::template push<gsVisitorNeumann<T> >(
    //    m_ppde.bcFirstKind().neumannSides() );

    // Neuman conditions of second kind
    Base::template push<gsVisitorNeumannBiharmonic<T> >(
        m_ppde.bcSecondKind().neumannSides() );

    if ( m_options.getInt("InterfaceStrategy") == iFace::dg )
        gsWarn <<"DG option ignored.\n";

    /*
    // If requested, force Dirichlet boundary conditions by Nitsche's method
    this->template push<gsVisitorNitscheBiharmonic<T> >(
    m_ppde.bcSecondKind().dirichletSides() );
    */

    // Assembly is done, compress the matrix
    Base::finalize();
}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::computeDirichletDofsL2Proj(gsMultiPatch<T> & basisG1_L,
                                                                      gsMultiPatch<T> & basisG1_R,
                                                                      index_t n_tilde,
                                                                      index_t n_bar,
                                                                      size_t unk_)
{
    int num_boundary = 6;
    std::vector<bool> bInt(2);
    bInt.at(0) = false;
    bInt.at(1) = false;



    m_g1_ddof.resize( map_boundary.boundarySize(), m_system.unkSize(unk_)*m_system.rhs().cols());  //m_pde_ptr->numRhs() );

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>        globProjRhs;
    globProjRhs.setZero( map_boundary.boundarySize() + num_boundary, m_system.unkSize(unk_)*m_system.rhs().cols() );

    // Temporaries
    gsVector<T> quWeights;

    gsMatrix<T> rhsVals;
    gsMatrix<unsigned> globIdxAct, locIdxAct;
    gsMatrix<T> basisVals;

    gsMapData<T> md(NEED_MEASURE);



    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
              iter = m_pde_ptr->bc().dirichletBegin();
          iter != m_pde_ptr->bc().dirichletEnd(); ++iter )
    {
        if (iter->isHomogeneous() )
            continue;

        GISMO_ASSERT(iter->function()->targetDim() == m_system.unkSize(unk_)*m_system.rhs().cols(),
                     "Given Dirichlet boundary function does not match problem dimension."
                         <<iter->function()->targetDim()<<" != "<<m_system.unkSize(unk_)<<"\n");

        const size_t unk = iter->unknown();
        if(unk!=unk_)
            continue;
        const int patchIdx   = iter->patch();
        const gsBasis<T> & basis = (m_bases[unk])[patchIdx];

        const gsGeometry<T> & patch = m_pde_ptr->patches()[patchIdx];

        // Set up quadrature to degree+1 Gauss points per direction,
        // all lying on iter->side() except from the direction which
        // is NOT along the element

        gsGaussRule<T> bdQuRule(basis, 1.0, 1, iter->side().direction());

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(iter->side());

        for(; bdryIter->good(); bdryIter->next() )
        {
            bdQuRule.mapTo( bdryIter->lowerCorner(), bdryIter->upperCorner(),
                            md.points, quWeights);

            //geoEval->evaluateAt( md.points );
            patch.computeMap(md);

            // the values of the boundary condition are stored
            // to rhsVals. Here, "rhs" refers to the right-hand-side
            // of the L2-projection, not of the PDE.
            rhsVals = iter->function()->eval( m_pde_ptr->domain()[patchIdx].eval( md.points ) );

            basis.eval_into( md.points, basisVals);

            index_t temp_i = basisVals.rows();
            basisVals.conservativeResize(basisVals.rows() + num_boundary,basisVals.cols());

            if (patchIdx == 0)
            {
                basisVals.row(temp_i) = basisG1_L.patch(0).eval(md.points);
                basisVals.row(temp_i+1) = basisG1_L.patch(n_tilde-1).eval(md.points);
                basisVals.row(temp_i+2) = basisG1_L.patch(n_tilde).eval(md.points);
                basisVals.row(temp_i+3) = basisG1_L.patch(n_tilde+n_bar-1).eval(md.points);

                basisVals.row(temp_i+4) = basisG1_L.patch(1).eval(md.points);
                basisVals.row(temp_i+5) = basisG1_L.patch(n_tilde-2).eval(md.points);
            }
            if (patchIdx == 1)
            {
                basisVals.row(temp_i) = basisG1_R.patch(0).eval(md.points);
                basisVals.row(temp_i+1) = basisG1_R.patch(n_tilde-1).eval(md.points);
                basisVals.row(temp_i+2) = basisG1_R.patch(n_tilde).eval(md.points);
                basisVals.row(temp_i+3) = basisG1_R.patch(n_tilde+n_bar-1).eval(md.points);

                basisVals.row(temp_i+4) = basisG1_R.patch(1).eval(md.points);
                basisVals.row(temp_i+5) = basisG1_R.patch(n_tilde-2).eval(md.points);
            }


            // Indices involved here:
            // --- Local index:
            // Index of the basis function/DOF on the patch.
            // Does not take into account any boundary or interface conditions.
            // --- Global Index:
            // Each DOF has a unique global index that runs over all patches.
            // This global index includes a re-ordering such that all eliminated
            // DOFs come at the end.
            // The global index also takes care of glued interface, i.e., corresponding
            // DOFs on different patches will have the same global index, if they are
            // glued together.
            // --- Boundary Index (actually, it's a "Dirichlet Boundary Index"):
            // The eliminated DOFs, which come last in the global indexing,
            // have their own numbering starting from zero.

            // Get the global indices (second line) of the local
            // active basis (first line) functions/DOFs:
            basis.active_into(md.points.col(0), locIdxAct );

            map_boundary.localToGlobal( locIdxAct, patchIdx, globIdxAct);
            // Out of the active functions/DOFs on this element, collect all those
            // which correspond to a boundary DOF.
            // This is checked by calling mapper.is_boundary_index( global Index )

            // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
            // something like a "element-wise index"
            std::vector<index_t> eltBdryFcts;
            eltBdryFcts.reserve(map_boundary.boundarySize());
            for( index_t i=0; i < globIdxAct.rows(); i++)
                if( map_boundary.is_boundary_index( globIdxAct(i,0)) )
                    if (!(map_interface.is_boundary(locIdxAct(i), patchIdx)))
                        eltBdryFcts.push_back( i );

            if (patchIdx == 0)
            {
                eltBdryFcts.push_back( globIdxAct.rows());
                eltBdryFcts.push_back( globIdxAct.rows()+1);
                eltBdryFcts.push_back( globIdxAct.rows()+2);
                eltBdryFcts.push_back( globIdxAct.rows()+3);

                if (bInt.at(0))
                    eltBdryFcts.push_back( globIdxAct.rows()+4);
                if (bInt.at(1))
                    eltBdryFcts.push_back( globIdxAct.rows()+5);
            }
            if (patchIdx == 1)
            {
                eltBdryFcts.push_back( globIdxAct.rows());
                eltBdryFcts.push_back( globIdxAct.rows()+1);
                eltBdryFcts.push_back( globIdxAct.rows()+2);
                eltBdryFcts.push_back( globIdxAct.rows()+3);

                if (bInt.at(0))
                    eltBdryFcts.push_back( globIdxAct.rows()+4);
                if (bInt.at(1))
                    eltBdryFcts.push_back( globIdxAct.rows()+5);
            }


            // Do the actual assembly:
            for( index_t k=0; k < md.points.cols(); k++ )
            {
                const T weight_k = quWeights[k] * md.measure(k);

                // Only run through the active boundary functions on the element:
                for( size_t i0=0; i0 < eltBdryFcts.size(); i0++ )
                {
                    unsigned ii, jj;
                    // Each active boundary function/DOF in eltBdryFcts has...
                    // ...the above-mentioned "element-wise index"
                    const index_t i = eltBdryFcts[i0];
                    // ...the boundary index.
                    if ( i < globIdxAct.rows())
                        ii = map_boundary.global_to_bindex( globIdxAct( i ));
                    else
                        ii = map_boundary.boundarySize() + i - globIdxAct.rows();

                    for( size_t j0=0; j0 < eltBdryFcts.size(); j0++ )
                    {
                        const index_t j = eltBdryFcts[j0];
                        if ( j < globIdxAct.rows())
                            jj = map_boundary.global_to_bindex( globIdxAct( j ));
                        else
                            jj = map_boundary.boundarySize() + j - globIdxAct.rows();

                        // Use the "element-wise index" to get the needed
                        // function value.
                        // Use the boundary index to put the value in the proper
                        // place in the global projection matrix.
                        projMatEntries.add(ii, jj, weight_k * basisVals(i,k) * basisVals(j,k));
                    } // for j
                    globProjRhs.row(ii) += weight_k *  basisVals(i,k) * rhsVals.col(k).transpose();
                } // for i
            } // for k
        } // bdryIter
    } // boundaryConditions-Iterator

    gsSparseMatrix<T> globProjMat( map_boundary.boundarySize()+ num_boundary,
                                   map_boundary.boundarySize()+ num_boundary);
    globProjMat.setFrom( projMatEntries );
    globProjMat.makeCompressed();

    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    typename gsSparseSolver<T>::CGDiagonal solver;
    m_g1_ddof = solver.compute( globProjMat ).solve ( globProjRhs );
}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::computeDirichletDofsL2Proj(std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>> & basisG1,
                                                                      std::vector<index_t> n_tilde,
                                                                      std::vector<index_t> n_bar,
                                                                      size_t unk_)
{
    int num_boundary = 4 * m_pde_ptr->patches().interfaces().size();

    m_g1_ddof.resize( map_boundary.boundarySize(), m_system.unkSize(unk_)*m_system.rhs().cols());  //m_pde_ptr->numRhs() );

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>        globProjRhs;
    globProjRhs.setZero( map_boundary.boundarySize() + num_boundary, m_system.unkSize(unk_)*m_system.rhs().cols() );

    // Temporaries
    gsVector<T> quWeights;

    gsMatrix<T> rhsVals;
    gsMatrix<unsigned> globIdxAct, locIdxAct;
    gsMatrix<T> basisVals;

    gsMapData<T> md(NEED_MEASURE);



    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
              iter = m_pde_ptr->bc().dirichletBegin();
          iter != m_pde_ptr->bc().dirichletEnd(); ++iter )
    {
        if (iter->isHomogeneous() )
            continue;

        GISMO_ASSERT(iter->function()->targetDim() == m_system.unkSize(unk_)*m_system.rhs().cols(),
                     "Given Dirichlet boundary function does not match problem dimension."
                         <<iter->function()->targetDim()<<" != "<<m_system.unkSize(unk_)<<"\n");

        const size_t unk = iter->unknown();
        if(unk!=unk_)
            continue;
        const int patchIdx   = iter->patch();
        const gsBasis<T> & basis = (m_bases[unk])[patchIdx];

        const gsGeometry<T> & patch = m_pde_ptr->patches()[patchIdx];

        // Set up quadrature to degree+1 Gauss points per direction,
        // all lying on iter->side() except from the direction which
        // is NOT along the element

        gsGaussRule<T> bdQuRule(basis, 1.0, 1, iter->side().direction());

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(iter->side());

        for(; bdryIter->good(); bdryIter->next() )
        {
            bdQuRule.mapTo( bdryIter->lowerCorner(), bdryIter->upperCorner(),
                            md.points, quWeights);

            //geoEval->evaluateAt( md.points );
            patch.computeMap(md);

            // the values of the boundary condition are stored
            // to rhsVals. Here, "rhs" refers to the right-hand-side
            // of the L2-projection, not of the PDE.
            rhsVals = iter->function()->eval( m_pde_ptr->domain()[patchIdx].eval( md.points ) );

            basis.eval_into( md.points, basisVals);

            // Indices involved here:
            // --- Local index:
            // Index of the basis function/DOF on the patch.
            // Does not take into account any boundary or interface conditions.
            // --- Global Index:
            // Each DOF has a unique global index that runs over all patches.
            // This global index includes a re-ordering such that all eliminated
            // DOFs come at the end.
            // The global index also takes care of glued interface, i.e., corresponding
            // DOFs on different patches will have the same global index, if they are
            // glued together.
            // --- Boundary Index (actually, it's a "Dirichlet Boundary Index"):
            // The eliminated DOFs, which come last in the global indexing,
            // have their own numbering starting from zero.

            // Get the global indices (second line) of the local
            // active basis (first line) functions/DOFs:
            basis.active_into(md.points.col(0), locIdxAct );

            map_boundary.localToGlobal( locIdxAct, patchIdx, globIdxAct);
            // Out of the active functions/DOFs on this element, collect all those
            // which correspond to a boundary DOF.
            // This is checked by calling mapper.is_boundary_index( global Index )

            // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
            // something like a "element-wise index"
            std::vector<index_t> eltBdryFcts;
            eltBdryFcts.reserve(map_boundary.boundarySize());
            for( index_t i=0; i < globIdxAct.rows(); i++)
                if( map_boundary.is_boundary_index( globIdxAct(i,0)) )
                    if (!(map_interface.is_boundary(locIdxAct(i), patchIdx)))
                        eltBdryFcts.push_back( i );

            index_t temp_i = basisVals.rows();
            basisVals.conservativeResize(basisVals.rows() + num_boundary,basisVals.cols());

            std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>>::iterator i_face;
            for (i_face=basisG1.equal_range(patchIdx).first; i_face!=basisG1.equal_range(patchIdx).second; ++i_face)
            {
                std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>::iterator i_side;
                for (i_side=i_face->second.begin(); i_side!=i_face->second.end(); ++i_side)
                {
                    for (std::map<index_t, gsMultiPatch<real_t>>::iterator i_mp=i_side->second.begin();
                         i_mp!=i_side->second.end(); ++i_mp)
                    {
                        gsMultiPatch<real_t> mp_side = i_mp->second;
                        index_t mp_index = i_mp->first;
                        if (iter->side().index() >= 3) // North or south
                        {
                            if (i_side->first <= 2) // West or east
                            {
                                if (iter->side().index() == 3)
                                {
                                    basisVals.row(temp_i + mp_index * 4) = mp_side.patch(0).eval(md.points);
                                    basisVals.row(temp_i + mp_index * 4 + 2) =
                                        mp_side.patch(n_tilde.at(mp_index)).eval(md.points);

                                    eltBdryFcts.push_back(globIdxAct.rows() + mp_index * 4);
                                    eltBdryFcts.push_back(globIdxAct.rows() + mp_index * 4 + 2);
                                }
                                else if (iter->side().index() == 4)
                                {
                                    basisVals.row(temp_i + mp_index * 4 + 1) =
                                        mp_side.patch(n_tilde.at(mp_index) - 1).eval(md.points);
                                    basisVals.row(temp_i + mp_index * 4 + 3) =
                                        mp_side.patch(n_tilde.at(mp_index) + n_bar.at(mp_index) - 1).eval(md.points);

                                    eltBdryFcts.push_back(globIdxAct.rows() + mp_index * 4 + 1);
                                    eltBdryFcts.push_back(globIdxAct.rows() + mp_index * 4 + 3);
                                }

                            }
                        }
                        if (iter->side().index() <= 2) // West or east
                        {
                            if (i_side->first >= 3) // North or south
                            {
                                if (iter->side().index() == 1)
                                {
                                    basisVals.row(temp_i + mp_index * 4) = mp_side.patch(0).eval(md.points);
                                    basisVals.row(temp_i + mp_index * 4 + 2) =
                                        mp_side.patch(n_tilde.at(mp_index)).eval(md.points);

                                    eltBdryFcts.push_back(globIdxAct.rows() + mp_index * 4);
                                    eltBdryFcts.push_back(globIdxAct.rows() + mp_index * 4 + 2);
                                }
                                else if (iter->side().index() == 2)
                                {
                                    basisVals.row(temp_i + mp_index * 4 + 1) =
                                        mp_side.patch(n_tilde.at(mp_index) - 1).eval(md.points);
                                    basisVals.row(temp_i + mp_index * 4 + 3) =
                                        mp_side.patch(n_tilde.at(mp_index) + n_bar.at(mp_index) - 1).eval(md.points);

                                    eltBdryFcts.push_back(globIdxAct.rows() + mp_index * 4 + 1);
                                    eltBdryFcts.push_back(globIdxAct.rows() + mp_index * 4 + 3);
                                }
                            }
                        }
                    }
                }
            }

            // Do the actual assembly:
            for( index_t k=0; k < md.points.cols(); k++ )
            {
                const T weight_k = quWeights[k] * md.measure(k);

                // Only run through the active boundary functions on the element:
                for( size_t i0=0; i0 < eltBdryFcts.size(); i0++ )
                {
                    unsigned ii, jj;
                    // Each active boundary function/DOF in eltBdryFcts has...
                    // ...the above-mentioned "element-wise index"
                    const index_t i = eltBdryFcts[i0];
                    // ...the boundary index.
                    if ( i < globIdxAct.rows())
                        ii = map_boundary.global_to_bindex( globIdxAct( i ));
                    else
                        ii = map_boundary.boundarySize() + i - globIdxAct.rows();

                    for( size_t j0=0; j0 < eltBdryFcts.size(); j0++ )
                    {
                        const index_t j = eltBdryFcts[j0];
                        if ( j < globIdxAct.rows())
                            jj = map_boundary.global_to_bindex( globIdxAct( j ));
                        else
                            jj = map_boundary.boundarySize() + j - globIdxAct.rows();

                        // Use the "element-wise index" to get the needed
                        // function value.
                        // Use the boundary index to put the value in the proper
                        // place in the global projection matrix.
                        projMatEntries.add(ii, jj, weight_k * basisVals(i,k) * basisVals(j,k));
                    } // for j
                    globProjRhs.row(ii) += weight_k *  basisVals(i,k) * rhsVals.col(k).transpose();
                } // for i
            } // for k
        } // bdryIter
    } // boundaryConditions-Iterator

    gsSparseMatrix<T> globProjMat( map_boundary.boundarySize()+ num_boundary,
                                   map_boundary.boundarySize()+ num_boundary);
    globProjMat.setFrom( projMatEntries );
    globProjMat.makeCompressed();

    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    typename gsSparseSolver<T>::CGDiagonal solver;
    m_g1_ddof = solver.compute( globProjMat ).solve ( globProjRhs );

}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::writeParaview(const gsField<T> & field,
                                                         std::string const & fn,
                                                         gsMultiPatch<T> mp_L,
                                                         gsMultiPatch<T> mp_R,
                                                         unsigned npts,
                                                         bool mesh)
{
    const unsigned n = field.nPieces();
    gsParaviewCollection collection(fn);
    std::string fileName;

    for ( unsigned i=0; i < n; ++i ) // Patches
    {
        //const gsBasis<T> & dom = field.isParametrized() ?
        //                         field.igaFunction(i).basis() : field.patch(i).basis();

        fileName = fn + util::to_string(i);
        //writeSinglePatchField( field, i, fileName, npts );

        const gsFunction<T> & geometry = field.patch(i);
        const gsFunction<T> & parField = field.function(i);

        const int n = geometry.targetDim();
        const int d = geometry.domainDim();

        gsMatrix<T> ab = geometry.support();
        gsVector<T> a = ab.col(0);
        gsVector<T> b = ab.col(1);

        gsVector<unsigned> np = uniformSampleCount(a, b, npts);
        gsMatrix<T> pts = gsPointGrid(a, b, np);

        gsMatrix<T> eval_geo = geometry.eval(pts);//pts
        gsMatrix<T> eval_field = field.isParametric() ? parField.eval(pts) : parField.eval(eval_geo);

        // Hier g1 basis dazu addieren!!!
        if (i == 0)
        {
            for (unsigned j = 0; j < mp_L.nPatches(); j++)
            {
                gsField<> solField(m_ppde.patches().patch(i),mp_L.patch(j));
                eval_field += solField.value(pts);
            }

        }
        if (i == 1)
        {
            for (unsigned j = 0; j < mp_R.nPatches(); j++)
            {
                gsField<> solField(m_ppde.patches().patch(i),mp_R.patch(j));
                eval_field += solField.value(pts);
            }

        }

        //gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
        //gsField<> exact( field.patch(i), solVal, false );
        //eval_field -= exact.value(pts);

        if ( 3 - d > 0 )
        {
            np.conservativeResize(3);
            np.bottomRows(3-d).setOnes();
        }
        else if (d > 3)
        {
            gsWarn<< "Cannot plot 4D data.\n";
            return;
        }

        if ( 3 - n > 0 )
        {
            eval_geo.conservativeResize(3,eval_geo.cols() );
            eval_geo.bottomRows(3-n).setZero();
        }
        else if (n > 3)
        {
            gsWarn<< "Data is more than 3 dimensions.\n";
        }

        if ( eval_field.rows() == 2)
        {
            eval_field.conservativeResize(3,eval_geo.cols() );
            eval_field.bottomRows(1).setZero(); // 3-field.dim()
        }

        gsWriteParaviewTPgrid(eval_geo, eval_field, np.template cast<index_t>(), fileName);


        collection.addPart(fileName, ".vts");
        if ( mesh )
        {
            fileName+= "_mesh";
            //writeSingleCompMesh(dom, field.patch(i), fileName);
            collection.addPart(fileName, ".vtp");
        }

    }
    collection.save();
}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::writeParaview(const gsField<T> & field,
                                                         std::string const & fn,
                                                         std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>> mp_g1,
                                                         unsigned npts,
                                                         bool mesh)
{
    const unsigned n = field.nPieces();
    gsParaviewCollection collection(fn);
    std::string fileName;

    for ( unsigned i=0; i < n; ++i ) // Patches
    {
        //const gsBasis<T> & dom = field.isParametrized() ?
        //                         field.igaFunction(i).basis() : field.patch(i).basis();

        fileName = fn + util::to_string(i);
        //writeSinglePatchField( field, i, fileName, npts );

        const gsFunction<T> & geometry = field.patch(i);
        const gsFunction<T> & parField = field.function(i);

        const int n = geometry.targetDim();
        const int d = geometry.domainDim();

        gsMatrix<T> ab = geometry.support();
        gsVector<T> a = ab.col(0);
        gsVector<T> b = ab.col(1);

        gsVector<unsigned> np = uniformSampleCount(a, b, npts);
        gsMatrix<T> pts = gsPointGrid(a, b, np);

        gsMatrix<T> eval_geo = geometry.eval(pts);//pts
        gsMatrix<T> eval_field = field.isParametric() ? parField.eval(pts) : parField.eval(eval_geo);

        // Hier g1 basis dazu addieren!!!
        std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>>::iterator i_face;
        for (i_face=mp_g1.equal_range(i).first; i_face!=mp_g1.equal_range(i).second; ++i_face)
        {
            std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>::iterator i_side;
            for (i_side = i_face->second.begin(); i_side != i_face->second.end(); ++i_side)
            {
                for (std::map<index_t, gsMultiPatch<real_t>>::iterator i_mp = i_side->second.begin();
                     i_mp != i_side->second.end(); ++i_mp)
                {
                    gsMultiPatch<real_t> mp_side = i_mp->second;
                    for (unsigned j = 0; j < mp_side.nPatches(); j++)
                    {
                        gsField<> solField(m_ppde.patches().patch(i),mp_side.patch(j));
                        eval_field += solField.value(pts);
                    }
                }
            }
        }
        //gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
        //gsField<> exact( field.patch(i), solVal, false );
        //eval_field -= exact.value(pts);

        if ( 3 - d > 0 )
        {
            np.conservativeResize(3);
            np.bottomRows(3-d).setOnes();
        }
        else if (d > 3)
        {
            gsWarn<< "Cannot plot 4D data.\n";
            return;
        }

        if ( 3 - n > 0 )
        {
            eval_geo.conservativeResize(3,eval_geo.cols() );
            eval_geo.bottomRows(3-n).setZero();
        }
        else if (n > 3)
        {
            gsWarn<< "Data is more than 3 dimensions.\n";
        }

        if ( eval_field.rows() == 2)
        {
            eval_field.conservativeResize(3,eval_geo.cols() );
            eval_field.bottomRows(1).setZero(); // 3-field.dim()
        }

        gsWriteParaviewTPgrid(eval_geo, eval_field, np.template cast<index_t>(), fileName);


        collection.addPart(fileName, ".vts");
        if ( mesh )
        {
            fileName+= "_mesh";
            //writeSingleCompMesh(dom, field.patch(i), fileName);
            //collection.addPart(fileName, ".vtp");
        }

    }
    collection.save();
}

} // namespace gismo
