/** @file gsPdeAssembler.hpp

    @brief Provides base interface for linear PDE's.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, J. Sogn
*/

#pragma once

#include <gsCore/gsDebug.h>
#include <gsAssembler/gsAssemblerUtils.h>

#include <gsCore/gsBasisEvaluator.hpp>

//JS2: What can we put in here?
namespace gismo {



template<class T>
void gsPdeAssembler<T>::initDofMapper(iFace::strategy interfaceStrategy, 
                                      dirichlet::strategy dirStrategy, bool hasMatrixRhs)
{
    // Setup the degrees of freedom, interface matching etc

    // Clean up
    freeAll(m_dofMapper);

    // Finding the total number of unknowns (including each
    // components) and resize m_dofMapper
    signed numOfUnknowns_tmp = m_unknownDim.sum();
    m_dofMapper.resize(numOfUnknowns_tmp);
    int bases_sz = m_bases.size();

    // For each unknown initialize dofMapper with correct basis
    const bool match     = (interfaceStrategy == iFace::glue);
    const bool eliminate = ( dirStrategy == dirichlet::elimination );
    int counter_tmp ;

    // If one basis for each unknown (and not for each component)
    if (bases_sz == m_unknownDim.size())
    {
        counter_tmp = 0;
        for (index_t k=0; k< m_unknownDim.size(); ++k )
        {
            m_bases[k].setTopology(m_patches);

            for(int j=0; j< m_unknownDim[k]; ++j )
            {
                GISMO_ASSERT ( counter_tmp <= numOfUnknowns_tmp,
                               "Counter is lager then number of unknown " );

                if ( eliminate && m_bconditions[k]) // Last argument is to check for NULL pointers
                    m_dofMapper[counter_tmp] =  m_bases[k].makeMapper(match, *m_bconditions[k]);
                else
                    m_dofMapper[counter_tmp] =  m_bases[k].makeMapper(match);

                ++counter_tmp;
            }
        }
    }
    // If one basis for each component (i.e. compatible discretisation)
    else if (bases_sz == numOfUnknowns_tmp)
    {
        counter_tmp = 0;
        for (index_t k=0; k< m_unknownDim.size(); ++k )
        {
            for(int j=0; j< m_unknownDim[k]; ++j )
            {
                GISMO_ASSERT ( counter_tmp <= numOfUnknowns_tmp,
                               "Counter is lager then number of unknown " );
                m_bases[counter_tmp].setTopology(m_patches);
                if ( eliminate && m_bconditions[k]) // Last argument is to check for NULL pointers
                    m_dofMapper[counter_tmp] =  m_bases[counter_tmp].makeMapper(match, *m_bconditions[k]);
                else
                    m_dofMapper[counter_tmp] =  m_bases[counter_tmp].makeMapper(match);

                ++counter_tmp;
            }
        }
    }
    else
    {
        GISMO_ERROR("Number of basis does not match to number of unknowns!");
    }


    // Set internal and total number of DoFs
    m_idofs = 0;
    m_dofs  = 0;

    //Vector valued Poisson case (Matrix equation)
    if (m_unknownDim.size() == 1 && hasMatrixRhs)
    {
        m_idofs += m_dofMapper[0]->freeSize();
        m_dofs  += m_dofMapper[0]->size();
    }
    else
    {
        for(signed k=0; k< numOfUnknowns_tmp; ++k)
        {
            m_idofs += m_dofMapper[k]->freeSize();
            m_dofs  += m_dofMapper[k]->size();
        }
    }
}


template<class T>
void gsPdeAssembler<T>::reconstructSolution(bool hasMatrixRhs)
{
    if (this->m_patches.nPatches() == 1)     // single-patch case
    {
        for (signed k=0; k< m_unknownDim.size(); ++k )
        {
            this->m_solutions.push_back(
                new gsField<T>(this->m_patches, reconstructPatchSolution(k,0, hasMatrixRhs)));
        }
    }
    else        // multipatch case
    {
        const gsMultiPatch<T> & mp = this->m_patches;
        for (signed k=0; k< m_unknownDim.size(); ++k ) // For each unknown
        {
            std::vector<gsFunction<T> * > sols;
            for (size_t np=0; np < mp.nPatches(); ++np ) // For each patch
            {
                sols.push_back( reconstructPatchSolution(k,np, hasMatrixRhs) );
            }
            this->m_solutions.push_back(new gsField<T>(mp, sols));
            sols.clear();
        }
    }
}

// Solution field(s)
template<class T>
gsFunction<T> * gsPdeAssembler<T>::reconstructPatchSolution(index_t unk, int p, bool hasMatrixRhs) const
{
    const gsMatrix<T> & data = m_sysSolution;

    //Check if each velocity component has the same basis (i.e => A1 = A2 = A3)
    bool equalBasis;
    if ((int) m_bases.size() == m_unknownDim.size())
        equalBasis = true;
    else if ((int) m_bases.size() == m_unknownDim.sum())
        equalBasis = false;
    else
        GISMO_ERROR("Number of basis does not match to number of unknowns!");

    //Find the basis index
    index_t basis_ind = 0;
    if (equalBasis)
        basis_ind = unk;
    else
    {
        for (index_t k = 0 ; k < unk; ++k)
            basis_ind += m_unknownDim[k];
    }


    const int sz  = m_bases[basis_ind][p].size();

    // Target dimension for unknown unk
    const int targetDim = m_unknownDim[unk];
    // The position of the unknown in the m_sysSolution vector/matrix
    int UnkPosSysSol = 0;
    // The (local) position of the component in the m_sysSolution vector/matrix starting from UnkPosSysSol
    int CompPosSysSol = 0;


    // Counter to help find correct index in m_dofMapper with
    // respect to the unknown unk.
    int compNr = 0; const int unkint = unk;
    for (int k = 0; k < unkint; ++k)
    {
        // Finding position in solution vector (for unknowns indexes less then unk)
        for (int k1 = 0; k1 < m_unknownDim[k]; ++k1)
            UnkPosSysSol += m_dofMapper[k1+compNr]->freeSize();

        compNr += m_unknownDim[k];
    }

    //Coefficient matrix
    gsMatrix<T> coeffs(sz, targetDim);

    for (index_t i = 0; i < sz; ++i)
    {
        for (int ki = 0; ki < targetDim; ++ki) // For each component in unknown unk
        {
            if ( m_dofMapper[ki+compNr]->is_free(i, p) ) // internal or interface
            {
                if (hasMatrixRhs) //Specially handle the if using matrix right hand side
                    coeffs(i,ki) = data(m_dofMapper[ki+compNr]->index(i, p),ki);
                else
                {
                    coeffs(i,ki) = data(UnkPosSysSol + CompPosSysSol + m_dofMapper[ki+compNr]->index(i, p),0);
                }
            }
            else // eliminated DoFs: fill with Dirichlet data
            {
                coeffs(i,ki) = m_fixedDofs[unk]( m_dofMapper[ki+compNr]->bindex(i, p),ki);
            }
            // Increase position with the size of the component
            CompPosSysSol += m_dofMapper[ki+compNr]->freeSize();
        }
        // Resett local position
        CompPosSysSol = 0;
    }

    return m_bases[basis_ind][p].makeGeometry( give(coeffs) );
}


// SKleiss: Note that this implementation is not useable for (T)HB-Splines!
//
// 1. Computation of the Dirichlet values explicitly uses gsPointGrid(rr)
// Computing a grid of evaluation points does not make sense for any locally
// refined basis.
// Also, "component(i)" is used.
// I'm afraid this makes sense ONLY FOR TENSOR-PRODUCT bases.
//
// 2. As of now (16.May 2014), the boundaryBasis of (T)HB-spline basis is not
// implemented, as far as I know.
//
// 3. gsInterpolate uses the anchors of the boundary basis.
// With truncated hierarchical B-splines, the use of classical anchors does
// not work, because functions might be truncated to zero at these points.
template<class T>
void gsPdeAssembler<T>::computeDirichletDofs(dirichlet::strategy dirStrategy)
{
    int component_count = 0;
    int bases_sz = m_bases.size();

    // Initialize the matrices in m_fixedDofs
    for (unsigned k=0; k< m_fixedDofs.size(); ++k)
    {
        int basisBmaxsz = 0;
        for (signed j = 0; j < m_unknownDim[k]; ++j)
        {
            if (basisBmaxsz < m_dofMapper[j + component_count]->boundarySize())
                basisBmaxsz = m_dofMapper[j + component_count]->boundarySize();
        }
        m_fixedDofs[k].resize( basisBmaxsz, m_unknownDim[k]);
        m_fixedDofs[k].setZero();
        component_count += m_unknownDim[k];
    }

    component_count = 0;
    for (unsigned k=0; k< m_fixedDofs.size(); ++k)
    {
        if (!m_bconditions[k]) //if bconditions[k] is an Null pointer
            continue;
        for ( typename gsBoundaryConditions<T>::const_iterator
                  it = m_bconditions[k]->dirichletBegin(); it != m_bconditions[k]->dirichletEnd(); ++it )
        {
            // If one basis for each unknown (and not for each component)
            if (bases_sz == m_unknownDim.size())
                {
                // Get DoFs on this boundary
                gsMatrix<unsigned> * boundary = m_bases[k][it->patch()].boundary(it->side()) ;

                // If the condition is homogeneous then fill with zeros
                if ( it->isHomogeneous())
                {//JS2: I think this can be removed: Start
                    for (index_t j=0; j!= boundary->size(); ++j)
                    {
                        const int ii= m_dofMapper[component_count]->bindex( (*boundary)(j) , it->patch() );
                        m_fixedDofs[k].row(ii).setZero();
                    }
                    delete boundary;
                    component_count += m_unknownDim[k];//JS2: I think this can be removed: End
                    continue;
                }

                // Get the side information
                int dir = it->side().direction();
                index_t param = ( it->side().parameter() ? 1 : 0);

                // Compute grid of points on the face ("face anchors")
                std::vector< gsVector<T> > rr;
                rr.reserve( m_patches.parDim() );

                for( int i=0; i < m_patches.parDim(); ++i)
                {
                    if ( i==dir )
                    {
                        gsVector<T> b(1);
                        b[0] = ( m_bases[k][it->patch()].component(i).support() ) (0, param);
                        rr.push_back(b);
                    }
                    else
                    {
                        rr.push_back( m_bases[k][it->patch()].component(i).anchors()->transpose() );
                    }
                }

                // Compute Dirichlet values
                gsMatrix<T> fpts =
                        it->function()->eval(m_patches.patch(it->patch()).eval(  gsPointGrid( rr ) ) );

                // Interpolate Dirichlet boundary
                gsBasis<T> * h = m_bases[k][it->patch()].boundaryBasis(it->side());

                gsGeometry<T> * geo = gsInterpolate( *h, *h->anchors(), fpts );
                const gsMatrix<T> & dVals =  geo->coefs();

                // Save corresponding boundary DoFs
                for (index_t j=0; j!= boundary->size(); ++j)
                {
                    const int ii= m_dofMapper[component_count]->bindex( (*boundary)(j) , it->patch() );
                    m_fixedDofs[k].row(ii) = dVals.row(j);
                }
                // REMARK: Note that correcting the right-hand-side by the contributions
                // from the fixed DoFs is done/has to be done in the assembling process.
                delete h;
                delete geo;
                delete boundary;
            }
            // If one basis for each component (i.e. compatible discretisation)
            else if (bases_sz == m_unknownDim.sum())
            {
                for (int comp = 0; comp < m_unknownDim[k]; ++comp)
                {
                    // Get DoFs on this boundary
                    gsMatrix<unsigned> * boundary = m_bases[comp+component_count][it->patch()].boundary(it->side()) ;

                    // If the condition is homogeneous then fill with zeros
                    if ( it->isHomogeneous() )
                    {//JS2: I think this can be removed: Start
                        for (index_t j=0; j!= boundary->size(); ++j)
                        {
                            const int ii= m_dofMapper[component_count]->bindex( (*boundary)(j) , it->patch() );
                            //m_fixedDofs[k].row(ii).setZero();
                            m_fixedDofs[k](ii,comp) = 0.0;
                        }
                        delete boundary;
                        component_count += m_unknownDim[k];//JS2: I think this can be removed: end
                        continue;
                    }
                    // Get the side information
                    int dir =  it->side().direction();
                    index_t param = ( it->side().parameter() ? 1 : 0);

                    // Compute grid of points on the face ("face anchors")
                    std::vector< gsVector<T> > rr;
                    rr.reserve( m_patches.parDim() );

                    for( int i=0; i < m_patches.parDim(); ++i)
                    {
                        if ( i==dir )
                        {
                            gsVector<T> b(1);
                            b[0] = ( m_bases[comp+component_count][it->patch()].component(i).support() ) (0, param);
                            rr.push_back(b);
                        }
                        else
                        {
                            rr.push_back( m_bases[comp+component_count][it->patch()].component(i).anchors()->transpose() );
                        }
                    }

                    // Compute Dirichlet values
                    gsMatrix<T> fpts =
                            it->function()->eval(m_patches.patch(it->patch()).eval(  gsPointGrid( rr ) ) );
                    gsMatrix<T> fpts_row = fpts.row(comp);

                    // Interpolate Dirichlet boundary
                    gsBasis<T> * h = m_bases[comp+component_count][it->patch()].boundaryBasis(it->side());

                    gsGeometry<T> * geo = gsInterpolate( *h, *h->anchors(), fpts_row );
                    const gsMatrix<T> & dVals =  geo->coefs();

                    // Save corresponding boundary DoFs
                    for (index_t j=0; j!= boundary->size(); ++j)
                    {
                        const int ii= m_dofMapper[component_count]->bindex( (*boundary)(j) , it->patch() );
                        m_fixedDofs[k](ii,comp) = dVals(j,0);
                    }
                    // REMARK: Note that correcting the right-hand-side by the contributions
                    // from the fixed DoFs is done/has to be done in the assembling process.
                    delete h;
                    delete geo;
                    delete boundary;
                }
            }
            else
            {
                GISMO_ERROR("Number of basis does not match to number of unknowns!");
            }
        }
        component_count += m_unknownDim[k];
    }
}


template<class T> void
gsPdeAssembler<T>::boundaryNeumann( const gsBasis<T> & B,
                                 const int patchIndex,
                                 const gsFunction<T> & f,
                                 const boxSide s,
                                 const gsDofMapper& mapper,
                                 gsVector<index_t> blockMapper)
{
    index_t tarDim = f.targetDim();
    std::vector< gsBasis<T> *> B_vec;
    std::vector< gsDofMapper *> mapper_vec;
    gsDofMapper mapper_copy = mapper;
    for (index_t k = 0; k< tarDim; ++k)
    {
        B_vec.push_back(B.clone());
        mapper_vec.push_back(&mapper_copy);
    }
    boundaryNeumann( B_vec, patchIndex, f, s ,mapper_vec, blockMapper);
    freeAll( B_vec);
    //freeAll( mapper_vec);
}

template<class T> void
gsPdeAssembler<T>::boundaryNeumann( std::vector< gsBasis<T> *>  const & B,
                                 const int patchIndex,
                                 const gsFunction<T> & f,
                                 const boxSide s,
                                 std::vector< gsDofMapper *> const & mapper,
                                 gsVector<index_t> blockMapper)
{
    //gsDebug<<"Vector Neumann boundary: side="<< s<<", patch="<<patchIndex<<"\n";
    index_t tarDim = f.targetDim();
    GISMO_ASSERT ((index_t) B.size() == tarDim, "Number of basis and target dimention does not match!");
    for (index_t comp = 0; comp < tarDim; ++comp)
    {
        const int d   = B[comp]->dim();
        int numBlock = blockMapper.size();

        GISMO_ASSERT ( numBlock <= f.targetDim(), "Miss match between number of block and target dimension of f" );

        //gsVector<index_t>;
        // Quadrature for boundary integral: we fix coordinates for
        // direction = dir to be the fixed coordinate on the edge/face
        gsVector<int> bd_intNodes = gsAssemblerUtils<T>::getNumIntNodesForSide( *B[comp], s.direction() );

        std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( m_patches.patch(patchIndex).evaluator(NEED_VALUE |NEED_JACOBIAN) );

        // Temporaries
        gsMatrix<T> fev;
        gsMatrix<T> localRhs;
        gsVector<T>  unormal(d);


        // iterate over all boundary grid cells
        for (typename gsDomainIterator<T>::uPtr domIter = B[comp]->makeDomainIterator(s);
             domIter->good(); domIter->next())
        {
            // Compute the quadrature rule (nodes and weights)
            domIter->computeQuadratureRule(bd_intNodes);

            // Evaluate the geometry on the Gauss points
            geoEval->evaluateAt(domIter->quNodes);

            // Evaluate the basis functions
            const index_t numActive = domIter->computeActiveDofs(*mapper[comp], patchIndex).rows();
            domIter->evaluateBasis();


            // Evaluate the Neumann data
            f.eval_component_into(geoEval->values(), comp, fev);

            localRhs.setZero(numActive, 1);

            for (index_t k=0; k!= domIter->numQuNodes(); ++k) // For all quadrature points
            {
                // Compute the outer normal vector on the side
                // The length of the outer normal (unormal.norm())
                // is the local boundary area segment.

                geoEval->outerNormal(k, s, unormal);


                // Sum up quadrature evaluations
                const T fff = domIter->quWeights[k] * fev(0,k) *unormal.norm();
                localRhs.col(0).noalias() += fff * domIter->basisValues().col(k);
            }

            // Push element contribution to the global load vector
            //
            // if right hand side does not have any special block structure
            // we use the matrix right hand side
            if (numBlock == 1 && blockMapper[0] ==  0)
            {
                for (index_t j=0; j!=numActive; ++j)
                {
                    // convert local DoFs index to global DoFs index
                    const unsigned jj = domIter->activeDofs[j];
                    if (mapper[comp]->is_free_index(jj))
                    {
                        if (m_rhs.cols() == tarDim) //If the right hand side is a matrix
                            m_rhs(jj,comp) += localRhs(j,0);
                        else
                            m_rhs(jj, 0) += localRhs(j,0);
                    }
                }
            }
            // if right hand side has block structure
            else
            {
                for (index_t j=0; j!=numActive; ++j)
                {
                    // convert local dof index to global dof index
                    const unsigned jj = domIter->activeDofs[j];
                    if (mapper[comp]->is_free_index(jj))
                    {
                        m_rhs(blockMapper[comp] + jj,0) += localRhs(j,0);
                    }
                }
            }
        }
    }
}

template<class T> void
gsPdeAssembler<T>::boundaryNitsche( const gsBasis<T> & B,
                                 const int patchIndex,
                                 const gsFunction<T> & f,
                                 const boxSide s,
                                 const gsDofMapper& mapper,
                                 gsMatrix<index_t> blockMapper,
                                 const T kappa)
{
    index_t tarDim = f.targetDim();
    std::vector< gsBasis<T> *> B_vec;
    std::vector< gsDofMapper *> mapper_vec;
    gsDofMapper mapper_copy = mapper;
    for (index_t k = 0; k< tarDim; ++k)
    {
        B_vec.push_back(B.clone());
        mapper_vec.push_back(&mapper_copy);
    }
    boundaryNitsche( B_vec, patchIndex, f, s ,mapper_vec, blockMapper, kappa);
    freeAll( B_vec);
    //freeAll( mapper_vec);
}

template<class T> void
gsPdeAssembler<T>::boundaryNitsche( std::vector< gsBasis<T> *>  const & B,
                                 const int patchIndex,
                                 const gsFunction<T> & f,
                                 const boxSide s,
                                 std::vector< gsDofMapper *> const & mapper,
                                 gsMatrix<index_t> blockMapper,
                                 const T kappa)
{
    index_t tDim = f.targetDim();
    GISMO_ASSERT ((index_t) B.size() == tDim, "Number of basis and target dimention does not match!");
    for (index_t comp = 0; comp < tDim; ++comp)
    {
        const int d   = B[comp]->dim();
        int numBlock = blockMapper.col(0).size();
        const T mu = gsAssemblerUtils<T>::getMu(*B[comp]);


        // Quadrature for boundary integral: we fix coordinates for
        // direction = dir to be the fixed coordinate on the edge/face
        gsVector<int> bd_intNodes = gsAssemblerUtils<T>::getNumIntNodesForSide( *B[comp], s.direction() );
        std::auto_ptr< gsGeometryEvaluator<T> >
            geoEval ( m_patches.patch(patchIndex).evaluator(NEED_VALUE |
                                    NEED_JACOBIAN | NEED_GRAD_TRANSFORM) );

        // Temporaries
        gsMatrix<T> fev, grads_k;
        gsVector<T> unormal(d);

        // basisData contains stacked the values and the gradients of all
        // basis functions at one quadrature node
        gsMatrix<T> basisData;

        // active basis functions at one quadrature node
        gsMatrix<unsigned> actives;

        // Local matrix and load vector
        gsMatrix<T> LM;
        gsMatrix<T> LB;


        // Make domain element iterator
        typename gsBasis<T>::domainIter domIter = B[comp]->makeDomainIterator(s);

        // Compute the quadrature rule (nodes and weights)
        domIter->computeQuadratureRule(bd_intNodes);

        // Iterate over all boundary grid cells
        for (; domIter->good(); domIter->next())
        {
            // Evaluate the geometry on the Gauss points
            geoEval->evaluateAt(domIter->quNodes);

            // Compute the active basis functions
            B[comp]->active_into(domIter->center, actives);
            const index_t numActive = actives.rows();

            // Map the local (indices) to global dofs
            mapper[comp]->localToGlobal(actives, patchIndex, actives);

            // Evaluate basis functions and their first derivatives
            B[comp]->evalAllDers_into(domIter->quNodes, 1, basisData);
            const typename gsMatrix<T>::Block ev          = basisData.topRows(numActive);
            const typename gsMatrix<T>::Block basisGrads  = basisData.middleRows(numActive, numActive*d);

            // Evaluate the Dirichlet data
            f.eval_component_into(geoEval->values(), comp, fev);

            // Initialize element matrices to zero
            LM.setZero(numActive, numActive);
            LB.setZero(numActive, 1);


            for (index_t k=0; k!= domIter->numQuNodes(); ++k) // For all quadrature points
            {
                // Compute the outer normal vector
                geoEval->outerNormal(k, s, unormal);

                // Integral transformation and quadrature weight
                const T fff = kappa * domIter->quWeights[k] * unormal.norm();

                // Compute the unit normal vector
                unormal.normalize();

                // Transform the basis gradients
                geoEval->transformGradients(k, basisGrads, grads_k);

                // Sum up quadrature point evaluations
                // volElement * F_k * ( dN_i * n - mu * N_i )
                LB.col(0).noalias() += fff * fev(0,k) * ( grads_k.transpose() * unormal - mu * ev.col(k) );

                LM.noalias() += fff * ( ev.col(k) * unormal.transpose() * grads_k
                                    +  (ev.col(k) * unormal.transpose() * grads_k).transpose()
                                    -  mu * ev.col(k) * ev.col(k).transpose() );
            }

            // Push element contribution to the global matrix and load vector
            //
            // if right hand side does not have any special block structure
            // we use the matrix right hand side
            if (numBlock == 1 && blockMapper(0,0) ==  0)
            {
                GISMO_ASSERT (m_rhs.cols() == tDim, "Target dimention and size of right hand size does not match");
                for (index_t i=0; i!=numActive; ++i)
                {
                    const unsigned ii = actives(i);
                    m_rhs(ii,comp) -= LB(i,0);

                    // Global matrix contribution should only be added ones
                    if (comp == 0)
                    {
                        for (index_t j=0; j!=numActive; ++j)
                        {
                            const unsigned jj = actives(j);
                            if ( !m_isSymmetric || jj <= ii )
                                m_matrix( ii, jj ) -= LM(i,j);
                        }
                    }
                }
            }
            // if right hand side has block structure
            else
            {
                for (index_t i=0; i!=numActive; ++i)
                {
                    // convert local dof index to global dof index
                    const unsigned ii = actives(i);
                    m_rhs(blockMapper(comp,0) + ii,0) -= LB(i,0);

                    for (index_t j=0; j!=numActive; ++j)
                    {
                        const unsigned jj = actives(j);
                        m_matrix(blockMapper(comp,0) + ii,blockMapper(comp,1) + jj) -=  LM(i,j);
                    }

                }
            }
        }
    }
}

// TO BE MOVED INSIDE INIT DOF MAPPERS WHEN ALL CODE IN PLACE
template<typename T>
void gsPdeAssembler<T>::remapToFullSystem (std::vector<gsDofMapper*> &input)
{
    std::vector<index_t> free_shifts(input.size());
    index_t              total_free=0;

//    std::vector<index_t> fixed_shifts(input.size());
//    index_t              total_fixed=0;

    for (size_t i=0; i<input.size();++i)
    {
        free_shifts[i]=total_free;
        total_free+=input[i]->freeSize();
//        fixed_shifts[i]=total_free;
//        total_fixed+=input[i]->boundarySize();
    }
    for (size_t i=0; i<input.size();++i)
    {
        input[i]->setShift(free_shifts[i]);
    }
}

template <typename T>
void gsPdeAssembler<T>::reconstructPatchCoefficients (
    const gsDofMapper &mapper,
          index_t      patch_id,
    const gsMatrix<T> &solution,
    const gsMatrix<T> &eliminated_values,
          gsMatrix<T> &result
    )
{
    result.resize(mapper.freeSize(), solution.cols());

    for (index_t row= 0; row<mapper.freeSize(); ++row)
    {
        if (mapper.is_free(row, patch_id))
        {
            index_t index = mapper.index(row,patch_id);
            result.row(row)=solution.row(index);
        }
        else
        {
            index_t index = mapper.bindex(row,patch_id);
            result.row(row).array()=eliminated_values.row(index);
        }
    }
    return result;
}


}// namespace gismo


