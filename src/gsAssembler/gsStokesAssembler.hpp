/** @file gsStokesAssembler.hpp

    @brief Assembler and solver for the Stokes problem

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#pragma once

#include <gsCore/gsDebug.h>
#include <gsUtils/gsStopwatch.h>

#include <gsAssembler/gsAssemblerUtils.h>
#include <gsAssembler/gsLocalToGlobal.h>
#include <gsAssembler/gsGenericAssembler.h> //Needed for preconditinor
#include <gsMatrix/gsMatrixPreconditioner.h> //Needed for preconditinor
#include <gsSolverDev/gsCoarsening.hpp>
#include <gsSolverDev/gsMultiGrid.h>
#include <gsSolver/gsBlockPreconditioner.h>
#include <gsSolver/gsMinimalResidual.hpp> //Solver
//#include <gismo.h>


namespace gismo {

template<class T>
void gsStokesAssembler<T>::initialize()
{
    //Tells the user that elimination strategy does not work in combination with dg.
    GISMO_ASSERT(!(m_interfaceStrategy == iFace::dg && 
                   m_dirStrategy == dirichlet::elimination),
                 "Must use Nitsche method for Dirichlet BC when using discontinuous Galerkin for patch interface");

    m_fixedDofs.resize(2); //One unknown u (velocity) and unknown p (pressure)
    m_fixedDofs[0].resize(0, m_unknownDim[0]);
    index_t tarDim = m_unknownDim[0];

    //Check if each velocity component has the same basis (i.e => A1 = A2 = A3)
    if ((int) m_bases.size() == m_unknownDim.size())
    {
        //Save and remove the bases for the pressure (added at the end)
        gsMultiBasis<T> pressure_bases = m_bases[m_bases.size()-1];
        m_bases.pop_back();
        //Add basis for velocity
        for (index_t k = 0; k < tarDim-1; ++k)
        {
            m_bases.push_back(gsMultiBasis<T>(m_bases[0]));
        }
        m_bases.push_back(pressure_bases);
    }
    else if ((int) m_bases.size() == m_unknownDim.sum())
    {}
    else
        GISMO_ERROR("Number of basis does not match to number of unknowns!");

    //initDofMapper needs to be called before initAssembler
    this->gsPdeAssembler<T>::initDofMapper(m_interfaceStrategy, m_dirStrategy);

    // Initialize the assembler by resize the global matrix and the right-hand side.
    m_matrix.resize(m_idofs,m_idofs);
    m_rhs.resize(m_idofs,1);
    m_matrix.setZero();
    m_rhs.setZero();

    // Estimate max non-zeros per row (only valid for tensor product bases right now)
    // Assumes that degree is equal for all velocity components and patches
    unsigned nzRowsPerColAii = 1;
    unsigned nzRowsPerColB = 1;
    //const gsBasis<T> & basis       = *m_bases[0][patchIndex];
    for (int i = 0; i < m_bases[0][0].dim(); ++i)
    {
        nzRowsPerColAii *= 2 * m_bases[0][0].degree(i) + 1;
        nzRowsPerColB   *= 2 * m_bases[tarDim][0].degree(i) + 1; //This is probably not right
    }

    gsVector<int> vectorA, vectorB;
    if (m_geoTrans == DIV_CONFORMING)
    {
        vectorA = gsVector<int>::Constant(m_idofs - m_dofMapper[tarDim]->freeSize(),
                                          nzRowsPerColAii*tarDim + nzRowsPerColB);
    }
    else
    {
        vectorA = gsVector<int>::Constant(m_idofs - m_dofMapper[tarDim]->freeSize(),
                                          nzRowsPerColAii + nzRowsPerColB);
    }
    vectorB = gsVector<int>::Constant(m_dofMapper[tarDim]->freeSize(), tarDim*nzRowsPerColB);
    gsVector<int> vectorAB(m_idofs);
    vectorAB << vectorA, vectorB;
    // Reserve space
    // \todo get the fragment of dofs in patch numbered patchIndex
    m_matrix.reserve( vectorAB );


    // Get a block view of the matrix and right hand side:
    //
    // (A11 A12 A13 B1^T)      (rhs_u1)
    // (A21 A22 A23 B2^T)  and (rhs_u2)
    // (A31 A32 A33 B3^T)      (rhs_u3)
    // (B1  B2  B3     0)      (rhs_p )

    //Initialize block structure of m_matrix and m_rhs
    gsVector<index_t> blockPositions(m_dofMapper.size());
    for (unsigned k = 0; k < m_dofMapper.size(); ++k)
    {
        blockPositions[k] = m_dofMapper[k]->freeSize();
    }
    gsVector<index_t> singleCol(1);
    singleCol << 1;

    m_matrixBlocks  = m_matrix.blockView(blockPositions, blockPositions);
    m_rhsBlocks= m_rhs.blockView(blockPositions, singleCol);
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
void gsStokesAssembler<T>::computeDirichletDofs(dirichlet::strategy dirStrategy)
{
    int component_count = 0;

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

            for (int comp = 0; comp < m_unknownDim[k]; ++comp)
            {
                // Get DoFs on this boundary
                gsMatrix<unsigned> * boundary = m_bases[comp+component_count][it->patch()].boundary(it->side()) ;

                // If the condition is homogeneous then fill with zeros
                if ( it->isHomogeneous() || dirStrategy == dirichlet::homogeneous)
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
                int dir = direction( it->side() );
                index_t param = (parameter( it->side() ) ? 1 : 0);

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

                //Do invers Piola tarnsform on Dirichlet data
                gsMatrix<T> fpts_row;
                gsMatrix<T> fptsPiolaInv;
                gsMatrix<T> fpts =
                        it->function()->eval(m_patches.patch(it->patch()).eval(  gsPointGrid( rr ) ) );
                if (k == 0 && m_geoTrans == DIV_CONFORMING)
                {
                    // Evaluate the geometry
                    std::auto_ptr< gsGeometryEvaluator<T> > geoEval (
                                m_patches.patch(it->patch()).evaluator(NEED_VALUE|NEED_GRAD_TRANSFORM|NEED_MEASURE));
                    geoEval->evaluateAt(gsPointGrid( rr ));

                    gsMatrix<T> paraPts = gsPointGrid( rr );
                    index_t numPts = paraPts.cols();
                    fptsPiolaInv.setZero(fpts.rows(), numPts);
                    for (index_t kp = 0; kp < numPts; ++kp)
                    {
                        fptsPiolaInv.col(kp) = (geoEval->gradTransform(kp).transpose())*fpts.col(kp)*geoEval->measure(kp);
                    }
                    fpts_row = fptsPiolaInv.row(comp);
                }
                else
                {
                    fpts_row = fpts.row(comp);
                }


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
        component_count += m_unknownDim[k];
    }
}

// Assembles the final system with all boundary conditions contained
template<class T>
void gsStokesAssembler<T>::assemble()
{

    // If the Dirichlet strategy is elimination then precompute
    // Dirichlet dofs (m_dofMapper excludes these from the system)
    if ( m_dirStrategy == dirichlet::elimination || 
         m_dirStrategy == dirichlet::homogeneous)
    {
        computeDirichletDofs(m_dirStrategy);
    }


    if (m_idofs == 0 ) // Are there any interior dofs ?
    {
        gsWarn << " No internal DOFs. Computed dirichlet boundary only.\n" <<"\n" ;
        return;
    }
    // Assemble the system matrix and right-hand side
    else
    {
        for (unsigned patchIndex=0; patchIndex < m_patches.nPatches(); ++patchIndex )
        {
            // Assemble stiffness matrix and rhs for the
            // local patch and add to m_matrix and m_rhs
            assemblePatch(patchIndex);
        }
    }

    // Enforce Dirichlet boundary conditions by Nitsche's method
    if ( m_dirStrategy == dirichlet::nitsche )
    {
        for ( typename gsBVProblem<T>::const_iterator
              BCiterator = m_bconditions[0]->dirichletBegin();
              BCiterator != m_bconditions[0]->dirichletEnd(); ++BCiterator )
        {
            gsStokesAssembler<T>::boundaryNitsche(BCiterator->patch(),
                                               *BCiterator->function(), BCiterator->side());
            //boundaryNitschePressure(BCiterator->patch(), *BCiterator->function(), BCiterator->side());
            //freeAll( basis_vector);
        }
    }

    // Enforce Neumann boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
          BCiterator = m_bconditions[0]->neumannBegin();
          BCiterator != m_bconditions[0]->neumannEnd(); ++BCiterator )
    {
        gsStokesAssembler<T>::boundaryNeumann(BCiterator->patch(),
                                           *BCiterator->function(),BCiterator->side() ); 
        //freeAll( B_vec);
    }


    // Enforce Pressure boundary conditions
    if (m_bconditions[1]) // If pressure BC
    {
        GISMO_ERROR("Dirichlet condition for pressure not implemented as it alone is not sufficient to get coersivity. Use instead grad(u)*n + pn = h(x) as a Neumann condition on u.");
        // Changeset 2820 as some code trying to add Dirichlet pressure BC.
    }
    m_matrix.makeCompressed();

}


template<class T>
void gsStokesAssembler<T>::assemblePatch(int patchIndex)
{
    //Comment explinations:
    //phi_j1 are the basis functions for the velocity field's first component
    //phi_j2 are the basis functions for the velocity field's second component
    //phi_j3 are the basis functions for the velocity field's thired component
    //psi_j  are the basis functions for the pressue

    //Target dimention for the velocity
    const int targetDim = m_unknownDim[0];
    const int tar2      = targetDim*targetDim;
    const T nu = m_nu;

    // Copy data
    std::vector<gsBasis<T> *> basis_u_vec;
    for (int k = 0; k < targetDim; ++k)
    {
        basis_u_vec.push_back(&m_bases[k][patchIndex]);
    }
    const gsBasis<T> & basis_p = m_bases[m_bases.size()-1][patchIndex];
    //const std::vector<gsDofMapper  *> mapper = m_dofMapper;
    const std::vector<gsMatrix<T> > ddof = m_fixedDofs;

    // Quadrature  nodes
    gsVector<int> numNodes(targetDim);
    for (index_t comp = 0 ; comp < targetDim; ++comp)
    {
        numNodes(comp) = basis_u_vec[comp]->maxDegree() + 1;
    }
    gsGaussRule<T> quad(numNodes);
    gsMatrix<T> quNodes  ; // Mapped nodes
    gsVector<T> quWeights; // Mapped weights

    // values of the right-hand side
    gsMatrix<T> rhsVals;

    // basis(Values/Grad/Div) contains the stacked (Values/Gradients/Divergence)
    // of all basis functions at one quadrature node for each collum
    gsMatrix<T> basisValues_u, basisGrad_u, basisValue_p, basisDiv_u;

    // Active basis functions at one quadrature node
    std::vector<gsMatrix<unsigned> > actives_vec(targetDim+1);

    // (dense) stiffness matrice within one grid cell
    // Contains the local A_{i,j} and B_i for i,j=1,...,targetDim
    std::vector<gsMatrix<T> > localBlockMatrice;
    localBlockMatrice.resize(targetDim + tar2);

    //Define some counters (used in inner loops)
    index_t actCount_i, actCount_j;

    //size of momentum matrix(A)
    index_t pShift = 0;
    for (index_t comp = 0 ; comp < targetDim; ++comp)
        pShift += m_dofMapper[comp]->freeSize();

    // local load vector
    gsVector<T> localRhs_p;
    std::vector<gsMatrix<T> > localRhs_u_vec(targetDim);

    // Evaluate the geometry
    std::auto_ptr< gsGeometryEvaluator<T> > geoEval (
                m_patches.patch(patchIndex).evaluator(NEED_VALUE   |
                                                      NEED_MEASURE |
                                                      NEED_GRAD_TRANSFORM |
                                                      NEED_2ND_DER) );


    // Make domain element iterator
    typename gsBasis<T>::domainIter domIter = basis_u_vec[0]->makeDomainIterator();

    //Initilize basis evaluator
    gsBasisEvaluator<T>     &basis_eval_u=*makeBasisEvaluator( basis_u_vec,
                                                               NEED_VALUE | NEED_GRAD | NEED_DIV,
                                                               &m_patches.patch(patchIndex),
                                                               m_geoTrans);
    gsBasisEvaluator<T>     &basis_eval_p=*makeBasisEvaluator( basis_p,
                                                               NEED_VALUE | NEED_GRAD,
                                                               &m_patches.patch(patchIndex),
                                                               INVERSE_COMPOSITION);


    //Initialize the local to global methods
    gsLocalToGlobalMapper<T> **vel_mappers = new gsLocalToGlobalMapper<T>*[targetDim*targetDim];
    gsLocalToGlobalMapper<T> **div_mappers = new gsLocalToGlobalMapper<T>*[targetDim];
    gsLocalToGlobalMapper<T> **rhs_mappers = new gsLocalToGlobalMapper<T>*[targetDim];

    gsSparseMatrix<T> *rhs_mod = new gsSparseMatrix<T>[targetDim] ;

    typedef gsShiftWriter<gsSparseMatrix<> > SW;
    typedef gsShiftWriter<gsMatrix<T> >      SWR;

    unsigned tShift=0;
    unsigned uShift=0;
    for (index_t comp_i=0; comp_i<targetDim; ++comp_i)
    {
        index_t comp_j = comp_i;
        index_t stopComp_j = comp_i + 1;
        if (m_geoTrans == DIV_CONFORMING)
        {
            comp_j = 0;
            stopComp_j = targetDim;
        }
        for (; comp_j < stopComp_j; ++comp_j) //For eash component of the parametric velocity
        {
            SW shift_m(m_matrix, tShift, uShift);
            SW shift_rhs(rhs_mod[comp_i], tShift, 0);
            vel_mappers[comp_i*targetDim+comp_j]= new gsL2GMapped<SW,SW>(
                        shift_m, shift_rhs,
                        *m_dofMapper[comp_i], *m_dofMapper[comp_j],
                        patchIndex);
            uShift+=m_dofMapper[comp_j]->freeSize();
        }
        // this would more elegantly be in another loop
        rhs_mod[comp_i].resize(m_matrix.rows(),ddof[0].rows());

        div_mappers[comp_i] = new gsL2GMappedMultiplier<gsSparseMatrix<T> >(
                    m_matrix, rhs_mod[comp_i],
                    *m_dofMapper[comp_i], *m_dofMapper[targetDim],
                    patchIndex, tShift, pShift);//ddof[0].col(comp_i)
        SWR shift_rhs(m_rhs, tShift, 0);
        rhs_mappers[comp_i] = new gsL2GMappedRhs<SWR> (shift_rhs, *m_dofMapper[comp_i], patchIndex);
        tShift+=m_dofMapper[comp_i]->freeSize();
        if (m_geoTrans == DIV_CONFORMING)
            uShift = 0;
    }


    // Start iteration over elements
    for (; domIter->good(); domIter->next())
    {
        gsVector<index_t> numActs(targetDim+1); //Number of Active components for velocity components and pressure

        // Map the Quadrature rule to the element
        quad.mapTo( domIter->lowerCorner(),
                    domIter->upperCorner(), quNodes, quWeights);

        for (index_t comp = 0 ; comp < targetDim+1; ++comp) //For each velocity component
        {
            // Compute the active basis functions
            if (comp != targetDim)
                basis_u_vec[comp]->active_into(domIter->center, actives_vec[comp]);
            else
                basis_p.active_into(domIter->center, actives_vec[comp]);
            numActs[comp]= actives_vec[comp].rows();
        }


        // Evaluate basis functions on element
        geoEval->evaluateAt(quNodes);

        basis_eval_u.evaluateAt (quNodes, *geoEval);
        basis_eval_p.evaluateAt (quNodes, *geoEval);

        basisValues_u= basis_eval_u.values(); //Transformed basis values
        basisGrad_u  = basis_eval_u.derivs(); //Transformed gardients values
        basisDiv_u   = basis_eval_u.divs();   //Transformed divergence values
        basisValue_p = basis_eval_p.values(); //Transformed values

        // Evaluate right-hand side at the geometry points
        if (m_rhs_function)
        {
            m_rhs_function->eval_into( geoEval->values(), rhsVals );
        }

        // Initialize local linear system to 0
        for (index_t comp_i = 0 ; comp_i < targetDim; ++comp_i) //For each velocity component
        {
            //The first targetDim**2 possitions are for the A_ij block matrix
            //the remaning targetdim possitions are the the B_i block matrix
            //Numbering of A_ij: A_00 A_01 A_02 A_10 A_11 ... A_21 A_22 B_0 B_1 B_2
            for (index_t comp_j = 0 ; comp_j < targetDim; ++comp_j) //For each velocity component
            {
                localBlockMatrice[comp_i*(targetDim)+comp_j].setZero(numActs[comp_i], numActs[comp_j]);//local_A_ij
            }
            localBlockMatrice[comp_i+tar2].setZero(numActs[targetDim], numActs[comp_i]);//localB_i
            localRhs_u_vec[comp_i].setZero(numActs[comp_i], targetDim);
        }
        localRhs_p.setZero(numActs[targetDim]);



        // Loop over quadrature nodes for velocity
        for (index_t node = 0; node < quNodes.cols(); ++node)
        {
            // weight * abs(det J), where J is geometry Jacobian.
            T weight = quWeights(node) * geoEval->measure(node);;

            actCount_i = 0;
            actCount_j = 0;
            for (index_t comp_i = 0 ; comp_i < targetDim; ++comp_i) //For each velocity component
            {

                //The first targetDim**2 possitions are for the A_ij block matrix
                //the remaning targetdim possitions are the the B_i block matrix
                //Numbering of A_ij: A_00 A_01 A_02 A_10 A_11 ... A_21 A_22 B_0 B_1 B_2

                //If using divergence presering transformation calculate off-diagonal
                //lacal matrix like A_ij, where i != j. Else do not calculate
                index_t comp_j = comp_i;
                index_t stopComp_j = comp_i + 1;
                if (m_geoTrans == DIV_CONFORMING)
                {
                    comp_j = 0;
                    stopComp_j = targetDim;
                    actCount_j = 0;
                }
                for ( ; comp_j < stopComp_j; ++comp_j) //For each velocity component
                {
                    for (index_t ai = 0; ai< numActs[comp_i]; ai++)
                    {
                        for (index_t aj = 0; aj< numActs[comp_j]; aj++)
                        {
                            //Siffness matrix for A_ij: nu (D phi, D phi)
                            //Frobenius mutiplication
                            localBlockMatrice[comp_i*targetDim + comp_j](ai,aj) += weight * nu *
                                    (basisGrad_u.block(tar2*(ai+actCount_i) ,node, tar2,1).transpose() *
                                     basisGrad_u.block(tar2*(aj+actCount_j) ,node, tar2,1) ).value();
                        }
                    }
                    if (m_geoTrans == DIV_CONFORMING)
                        actCount_j += numActs[comp_j];
                }
                //Calculateing the B_i local Matrices
                for (index_t ai = 0; ai< numActs[comp_i]; ai++)
                {
                    // B_i: (div(phi), psi)
                    for (index_t ap = 0; ap< numActs[targetDim]; ap++)
                    {
                        localBlockMatrice[tar2+comp_i](ap,ai) +=
                                weight * basisDiv_u(ai+actCount_i,node) * basisValue_p(ap,node);
                    }
                    //Calculateing the b_i local load vector: (phi,f)
                    localRhs_u_vec[comp_i](ai,comp_i) += weight *
                            (basisValues_u.block(targetDim*(ai+actCount_i) ,node, targetDim,1).transpose() *
                             rhsVals.col(node)).value();
                }
                actCount_i += numActs[comp_i];
                actCount_j += numActs[comp_j];
            }
        }  // end loop Gauss nodes for velocity

        //----- ADD LOCAL TO GLOBAL -----//
        for (index_t comp_i=0; comp_i<targetDim; ++comp_i)
        {
            index_t comp_j = comp_i;
            index_t stopComp_j = comp_i + 1;
            if (m_geoTrans == DIV_CONFORMING)
            {
                comp_j = 0;
                stopComp_j = targetDim;
            }
            for (; comp_j < stopComp_j; ++comp_j) //For eash component of the parametric velocity
            {
                vel_mappers[comp_i*targetDim+comp_j]->store(actives_vec[comp_i], actives_vec[comp_j], localBlockMatrice[comp_j + comp_i*targetDim]);
            }
            div_mappers[comp_i]->store(actives_vec[comp_i], actives_vec[targetDim], localBlockMatrice[comp_i+tar2].transpose());
            rhs_mappers[comp_i]->store(actives_vec[comp_i], actives_vec[comp_i],localRhs_u_vec[comp_i].col(comp_i));
        }
        // Use if div(u) = g(x)
        // Adding pressure rhs (compression term)
        //gsL2GMappedRhs<T> L2Gg(m_rhs, *m_dofMapper[targetDim], patchIndex,pShift);
        //L2Gg.store(actives_vec[targetDim], localRhs_p );


    } // end loop over all domain elements

    // add rhs modifications
    for (index_t comp_i=0; comp_i<targetDim; ++comp_i)
    {
        m_rhs-=rhs_mod[comp_i]*ddof[0].col(comp_i);
    }
    // free mappers
    for (index_t comp_i=0; comp_i<targetDim; ++comp_i)
    {
        index_t comp_j = comp_i;
        index_t stopComp_j = comp_i + 1;
        if (m_geoTrans == DIV_CONFORMING)
        {
            comp_j = 0;
            stopComp_j = targetDim;
        }
        for (; comp_j < stopComp_j; ++comp_j) //For eash component of the parametric velocity
        {
            delete vel_mappers[comp_i*targetDim+comp_j];
        }
        delete div_mappers[comp_i];
        delete rhs_mappers[comp_i];
    }
    delete[] vel_mappers;
    delete[] div_mappers;
    delete[] rhs_mappers;
    delete[] rhs_mod;

    delete &basis_eval_u;
    delete &basis_eval_p;
}

template<class T>  void
gsStokesAssembler<T>::boundaryNeumann( const int patchIndex,
                                    const gsFunction<T> & f,
                                    const boundary::side s)
{
    // Copy data
    const index_t tarDim = m_unknownDim[0];
    std::vector<gsBasis<T> *> basis_u_vec;
    for (int k = 0; k < tarDim; ++k)
    {
        basis_u_vec.push_back(&m_bases[k][patchIndex]);
    }
    const std::vector<gsDofMapper  *> & mapper = m_dofMapper;

    // Quadrature  nodes
    gsVector<int> numNodes(tarDim);
    for (index_t comp = 0 ; comp < tarDim; ++comp)
    {
        if (comp == direction(s))
            numNodes(comp) =  1;
        else
            numNodes(comp) = basis_u_vec[comp]->maxDegree() + 1;

    }
    gsGaussRule<T> quad(numNodes);

    gsMatrix<T> quNodes  ; // Mapped nodes
    gsVector<T> quWeights; // Mapped weights

    std::auto_ptr< gsGeometryEvaluator<T> >
            geoEval ( m_patches.patch(patchIndex).evaluator(NEED_VALUE | NEED_MEASURE) );

    // Temporaries
    gsMatrix<T> fev, basisValues;
    gsVector<T> unormal(tarDim);


    // Active basis functions at one quadrature node
    std::vector<gsMatrix<unsigned> > actives_vec(tarDim);

    // Local load vactor
    std::vector<gsMatrix<T> > localBlockRhs(tarDim);

    //Define some counters (used in inner loops)
    index_t actCount;

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIter = basis_u_vec[0]->makeDomainIterator(s);

    //Initilize basis evaluator
    gsBasisEvaluator<T>     &basis_eval_u=*makeBasisEvaluator( basis_u_vec,NEED_VALUE, &m_patches.patch(patchIndex), m_geoTrans);

    // Iterate over all boundary grid cells
    for (; domIter->good(); domIter->next())
    {
        gsVector<index_t> numActs(tarDim); //Number of Active components for basis components

        // Map the Quadrature rule to the element
        quad.mapTo( domIter->lowerCorner(),
                    domIter->upperCorner(), quNodes, quWeights);

        for (index_t comp = 0 ; comp < tarDim; ++comp) //For each basis component
        {
            // Compute the active basis functions
            basis_u_vec[comp]->active_into(domIter->center, actives_vec[comp]);
            numActs[comp]= actives_vec[comp].rows();

            // Local Dofs to global dofs
            mapper[comp]->localToGlobal(actives_vec[comp], patchIndex, actives_vec[comp]);
        }

        // Evaluate the geometry and basis functions on element
        geoEval->evaluateAt(quNodes);
        basis_eval_u.evaluateAt (quNodes, *geoEval);

        //ev and basisGrads
        basisValues = basis_eval_u.values(); //Transformed velocity values

        // Evaluate the Dirichlet data
        //JS2: NB This might need to change for PIOLA TRANSFORM
        f.eval_into(geoEval->values(), fev);

        // Initialize local linear system to 0
        for (index_t comp = 0 ; comp < tarDim; ++comp) //For each basis component
        {
            localBlockRhs[comp].setZero(numActs[comp], tarDim);
        }

        // Loop over quadrature nodes
        for (index_t node = 0; node < quNodes.cols(); ++node)
        {
            geoEval->outerNormal(node, s, unormal);

            T weight = quWeights(node) * unormal.norm();

            actCount = 0;

            for (index_t comp = 0 ; comp < tarDim; ++comp)
            {
                //Calculateing the right hand side
                for (index_t ai = 0; ai< numActs[comp]; ai++)
                {
                    localBlockRhs[comp](ai,comp) += weight *
                            (fev.col(node).transpose() *
                             basisValues.block(tarDim*(ai+actCount) ,node, tarDim,1)).value();
                }
                actCount += numActs[comp];
            }
        }

        //----- ADD LOCAL TO GLOBAL -----//
        for (index_t comp = 0 ; comp < tarDim; ++comp)
        {
            gsVector<int> loc2glob_v = actives_vec[comp].cast<int>();
            m_rhsBlocks(comp);
            for (index_t j=0; j!=numActs[comp]; ++j)
            {
                // convert local dof index to global dof index
                const index_t jj = loc2glob_v[j]; //actives_vec[comp]
                if (mapper[comp]->is_free_index(jj))
                {
                    m_rhsBlocks(comp)(jj,0) += localBlockRhs[comp](j,comp);
                }
            }
        }

    } // end loop over all domain elements

    // clean up other stuff
    delete &basis_eval_u;
}

template<class T>  void
gsStokesAssembler<T>::boundaryNitsche( const int patchIndex,
                                       const gsFunction<T> & f,
                                       const boundary::side s)
{
    // Copy data
    const index_t tarDim = m_unknownDim[0];
    const index_t tar2 = tarDim*tarDim;
    std::vector<gsBasis<T> *> basis_u_vec;
    for (int k = 0; k < tarDim; ++k)
    {
        basis_u_vec.push_back(&m_bases[k][patchIndex]);
    }
    const gsBasis<T> & basis_p = m_bases.back()[patchIndex];
    const T kappa = m_nu;

    //const T mu = gsAssemblerUtils<T>::getMu(*basis_u_vec[0]);
    T mu = 5*(basis_u_vec[1][0].degree(0) + 1);//Evans


    //const std::vector<gsDofMapper  *> & mapper = m_dofMapper;

    // Quadrature  nodes
    gsVector<int> numNodes(tarDim);
    for (index_t comp = 0 ; comp < tarDim; ++comp)
    {
        if (comp == direction(s))
            numNodes(comp) =  1;
        else
            numNodes(comp) = basis_u_vec[comp]->maxDegree() + 1;

    }
    gsGaussRule<T> quad(numNodes);
    gsMatrix<T> quNodes  ; // Mapped nodes
    gsVector<T> quWeights; // Mapped weights

    unsigned flags;
    if (m_geoTrans == DIV_CONFORMING)
        flags = NEED_VALUE | NEED_JACOBIAN | NEED_GRAD_TRANSFORM | NEED_2ND_DER | NEED_MEASURE;
    else
        flags = NEED_VALUE | NEED_JACOBIAN | NEED_GRAD_TRANSFORM | NEED_MEASURE;
    std::auto_ptr< gsGeometryEvaluator<T> >
            geoEval ( m_patches.patch(patchIndex).evaluator(flags) );

    // Temporaries
    gsMatrix<T> fev, tmp_gradTj, tmp_gradTi;
    gsVector<T> unormal(tarDim), tmp_gradNormi, tmp_gradNormj;
    // Zero matrix, needed for the LocalToGlobal function
    gsMatrix<T> ddofZero = gsMatrix<T>::Zero(m_fixedDofs[0].rows(),m_fixedDofs[0].cols());

    // basis(Values/Grad) contains the stacked (Values/Gradients)
    // of all basis functions at one quadrature node for each collum
    gsMatrix<T> basisValues, basisGrad, pressureValues;

    // Active basis functions at one quadrature node
    std::vector<gsMatrix<unsigned> > actives_vec(tarDim+1);

    // (dense) local stiffness matrice within one grid cell
    // Contains the local A_{i,j}
    std::vector<gsMatrix<T> > localBlockMatrice;
    localBlockMatrice.resize(tar2+tarDim);

    // Local load vactor
    std::vector<gsMatrix<T> > localBlockRhs(tarDim);
    gsVector<T> localRhs_p;

    //Define some counters (used in inner loops)
    index_t actCount_i, actCount_j;

    //size of momentum matrix(A)
    index_t pShift = 0;
    for (index_t comp = 0 ; comp < tarDim; ++comp)
        pShift += m_dofMapper[comp]->freeSize();

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIter = basis_u_vec[0]->makeDomainIterator(s);

    //Initilize basis evaluator
    gsBasisEvaluator<T>     &basis_eval_u=*makeBasisEvaluator( basis_u_vec,
                                                               NEED_VALUE | NEED_GRAD | NEED_DIV,
                                                               &m_patches.patch(patchIndex),
                                                               m_geoTrans);
    gsBasisEvaluator<T>     &basis_eval_p=*makeBasisEvaluator( basis_p,
                                                               NEED_VALUE | NEED_GRAD,
                                                               &m_patches.patch(patchIndex),
                                                               INVERSE_COMPOSITION);


    // Iterate over all boundary grid cells
    for (; domIter->good(); domIter->next())
    {
        gsVector<index_t> numActs(tarDim+1); //Number of Active components for basis components

        // Map the Quadrature rule to the element
        quad.mapTo( domIter->lowerCorner(),
                                      domIter->upperCorner(), quNodes, quWeights);

        for (index_t comp = 0 ; comp < tarDim; ++comp) //For each basis component
        {
            // Compute the active basis functions
            basis_u_vec[comp]->active_into(domIter->center, actives_vec[comp]);
            numActs[comp]= actives_vec[comp].rows();

        }
        basis_p.active_into(domIter->center, actives_vec[tarDim]);
        numActs[tarDim] = actives_vec[tarDim].rows();
        //mapper[tarDim]->localToGlobal(actives_vec[tarDim], patchIndex, actives_vec[tarDim]);

        // Evaluate the geometry and basis functions on element
        geoEval->evaluateAt(quNodes);

        //Find cell size.
        T h_Q = domIter->getPerpendicularCellSize();
        T jac_inf = geoEval->jacobians().maxCoeff();
        if (fabs(geoEval->jacobians().minCoeff()) > jac_inf)
            jac_inf = fabs(geoEval->jacobians().minCoeff());
        jac_inf = 1.0;
        T h_K = h_Q*jac_inf;
        //gsDebug<< "jac_inf: " << jac_inf << "  h_Q: "<< h_Q<< "  h_K: "<< h_K<< "  nu: "<< mu << "  mu: "<< mu/h_K<< " side: "<< s <<"\n";


        basis_eval_u.evaluateAt (quNodes, *geoEval);
        basis_eval_p.evaluateAt (quNodes, *geoEval);

        //ev and basisGrads
        basisValues = basis_eval_u.values(); //Transformed velocity values
        basisGrad   = basis_eval_u.derivs(); //Transformed gardients values
        pressureValues = basis_eval_p.values();

        // Evaluate the Dirichlet data
        //JS2: NB This might need to change for PIOLA TRANSFORM
        f.eval_into(geoEval->values(), fev);

        // Initialize local linear system to 0
        for (index_t comp_i = 0 ; comp_i < tarDim; ++comp_i) //For each basis component
        {
            //Numbering of A_ij: A_00 A_01 A_02 A_10 A_11 ... A_21 A_22 B_0 B_1 B_2
            for (index_t comp_j = 0 ; comp_j < tarDim; ++comp_j) //For each basis component
            {
                localBlockMatrice[comp_i*tarDim + comp_j].setZero(numActs[comp_i], numActs[comp_j]);//local_A_ij
            }
            localBlockMatrice[comp_i+tar2].setZero(numActs[tarDim], numActs[comp_i]);//localB_i
            localBlockRhs[comp_i].setZero(numActs[comp_i], tarDim);
        }
        localRhs_p.setZero(numActs[tarDim]);

        // Loop over quadrature nodes
        for (index_t node = 0; node < quNodes.cols(); ++node)
        {
            //gsDebug <<"node: "<<node<<"\n";
            // Compute the outer normal vector
            geoEval->outerNormal(node, s, unormal);

            T weight = kappa * quWeights(node) * unormal.norm();

            // Compute the unit normal vector
            unormal.normalize();
            actCount_i = 0;
            actCount_j = 0;
            for (index_t comp_i = 0 ; comp_i < tarDim; ++comp_i)
            {
                //Calculateing the right hand side
                for (index_t ai = 0; ai< numActs[comp_i]; ai++)
                {
                    tmp_gradTi.noalias() = basisGrad.block(tar2*(ai+actCount_i) ,node, tar2,1);
                    tmp_gradTi.resize(tarDim,tarDim);
                    tmp_gradNormi.noalias() = unormal.transpose()*tmp_gradTi;
                    localBlockRhs[comp_i](ai,comp_i) += weight *
                            ((mu/h_K)*(fev.col(node).transpose() *
                                       basisValues.block(tarDim*(ai+actCount_i) ,node, tarDim,1)).value()
                             - fev.col(node).dot((tmp_gradNormi)));
                    // B_i: (div(phi), psi)
                    for (index_t ap = 0; ap< numActs[tarDim]; ap++)
                    {
                        localBlockMatrice[tar2+comp_i](ap,ai) -= weight/kappa * pressureValues(ap,node) *
                                (unormal.transpose() * basisValues.block(tarDim*(ai+actCount_i) ,node, tarDim,1)).value();
                    }
                }

                //If using divergence presering transformation calculate off-diagonal
                //lacal matrix like A_ij, where i != j. Else do not calculate
                index_t comp_j = comp_i;
                index_t stopComp_j = comp_i + 1;
                if (m_geoTrans == DIV_CONFORMING)
                {
                    comp_j = 0;
                    stopComp_j = tarDim;
                    actCount_j = 0;
                }
                for ( ; comp_j < stopComp_j; ++comp_j) //For each velocity component
                {
                    for (index_t ai = 0; ai< numActs[comp_i]; ai++)
                    {
                        for (index_t aj = 0; aj< numActs[comp_j]; aj++)
                        {   //JS2: This can be optimised because of symmertri properties!
                            tmp_gradTi = basisGrad.block(tar2*(ai+actCount_i) ,node, tar2,1); ;
                            tmp_gradTj = basisGrad.block(tar2*(aj+actCount_j) ,node, tar2,1);
                            tmp_gradTi.resize(tarDim,tarDim);
                            tmp_gradTj.resize(tarDim,tarDim);

                            tmp_gradNormi.noalias() = unormal.transpose()*tmp_gradTi;
                            tmp_gradNormj.noalias() = unormal.transpose()*tmp_gradTj;
                            localBlockMatrice[comp_i*(tarDim)+comp_j](ai,aj) += weight *
                                    ((mu/h_K)*(basisValues.block(tarDim*(aj+actCount_j) ,node, tarDim,1).transpose() *
                                               basisValues.block(tarDim*(ai+actCount_i) ,node, tarDim,1)).value() -
                                     (basisValues.block(tarDim*(ai+actCount_i) ,node, tarDim,1).transpose() * tmp_gradNormj).value() -
                                     (basisValues.block(tarDim*(aj+actCount_j) ,node, tarDim,1).transpose() * tmp_gradNormi).value());
                        }
                    }
                    if (m_geoTrans == DIV_CONFORMING)
                        actCount_j += numActs[comp_j];
                }
                actCount_i += numActs[comp_i];
                actCount_j += numActs[comp_j];
            }
            for (index_t ap = 0; ap< numActs[tarDim]; ap++)
            {
                localRhs_p(ap) -= weight/kappa * unormal.dot(fev.col(node)) * pressureValues(ap,node);
            }
        }

        //----- ADD LOCAL TO GLOBAL -----//


        index_t tShift = 0;
        index_t uShift = 0;
        for (index_t comp_i = 0 ; comp_i < tarDim; ++comp_i) //For eash component of the parametric velocity
        {

            //If using divergence presering transformation calculate off-diagonal
            //lacal matrix like A_ij, where i != j. Else do not calculate
            index_t comp_j = comp_i;
            index_t stopComp_j = comp_i + 1;
            if (m_geoTrans == DIV_CONFORMING)
            {
                comp_j = 0;
                stopComp_j = tarDim;
            }
            for (; comp_j < stopComp_j; ++comp_j) //For eash component of the parametric velocity
            {
                //Adding momentum terms A_ij
                typedef gsShiftWriter<gsSparseMatrix<T> > SW;
                gsSparseMatrix<T> rhs_mod;
                rhs_mod.resize(m_matrix.rows(),ddofZero.rows() );

                SW myshift_view(m_matrix, tShift, uShift);
                SW myshift_view_rhs(rhs_mod, tShift, 0);
                gsL2GMapped<SW,SW> L2GA(myshift_view, myshift_view_rhs, *m_dofMapper[comp_i], *m_dofMapper[comp_j], patchIndex);
                L2GA.store(actives_vec[comp_i], actives_vec[comp_j], localBlockMatrice[comp_j + comp_i*tarDim]);
                m_rhs-=rhs_mod*ddofZero.col(comp_j);
                uShift += m_dofMapper[comp_j]->freeSize();
            }
            //Adding divergence terms B and BT  -(v*n,q)_\Gamma
            gsSparseMatrix<T> rhs_mod;
            rhs_mod.resize(m_matrix.rows(),ddofZero.rows() );
            gsL2GMappedMultiplier<gsSparseMatrix<T> > L2GB(m_matrix,rhs_mod, *m_dofMapper[comp_i], *m_dofMapper[tarDim],
                                                           patchIndex,tShift, pShift);
            L2GB.store(actives_vec[comp_i], actives_vec[tarDim], localBlockMatrice[comp_i+tar2].transpose());
            m_rhs-=rhs_mod*ddofZero.col(comp_i);

            //Adding rhs terms
            typedef gsShiftWriter<gsMatrix<T> > SWR;
            SWR shift_rhs(m_rhs, tShift, 0);
            gsL2GMappedRhs<SWR > L2Gf(shift_rhs, *m_dofMapper[comp_i], patchIndex);
            L2Gf.store(actives_vec[comp_i], actives_vec[comp_i], localBlockRhs[comp_i].col(comp_i));

            tShift += m_dofMapper[comp_i]->freeSize();
            if (m_geoTrans == DIV_CONFORMING)
                uShift = 0;
        }
        //Add -(u*n,q)_\Gamma terms
        typedef gsShiftWriter<gsMatrix<T> > SWR;
        SWR shift_rhs(m_rhs, pShift, 0);
        gsL2GMappedRhs<SWR > L2Gp(shift_rhs, *m_dofMapper[tarDim], patchIndex);
        L2Gp.store(actives_vec[tarDim], actives_vec[tarDim], localRhs_p );

    } // end loop over all domain elements

    // clean up other stuff
    delete &basis_eval_u;
    delete &basis_eval_p;
}



// Solves the linear system and fills up \a m_sysSolution
template<class T>
void gsStokesAssembler<T>::solveSystem()
{
    std::cout << "Solve linear system of size: " << m_idofs << "\n";

    //Direct solvers
    //Eigen::SparseLU< gsSparseMatrix<T> > solver;
    //Eigen::BiCGSTAB< gsSparseMatrix<T> > solver;
    Eigen::SparseQR< gsSparseMatrix<T>, Eigen::COLAMDOrdering<int> > solver; //COLAMDOrdering<int>
    gsInfo << "Using direct solver" << std::endl;

    //solver.compute( m_matrix );
    solver.analyzePattern(m_matrix);
    solver.factorize(m_matrix);
    gsDebug << "   Eigen  lastErrorMessage: " << solver.lastErrorMessage () << "\n";
    m_sysSolution = solver.solve(m_rhs);
    //_sysSolution = solver.compute(m_matrix).solve (m_rhs);

    gsInfo << "Solved with SparceLU\n" ;
    //gsInfo << "Solved with BiCGSTAB \n" ;
}

// Solution field(s)
template<class T>
gsFunction<T> * gsStokesAssembler<T>::reconstructPatchSolution(index_t unk, int p, bool hasMatrixRhs) const
{
    gsDebug << "Im now using the derived class reconstruction\n";
    const gsMatrix<T> & data = m_sysSolution;

    GISMO_ASSERT ( (signed) m_bases.size() == m_unknownDim.sum(), "Implementation assumes one basis for each unknown component");

    // Target dimension for unknown unk
    const index_t targetDim = m_unknownDim[unk];
    // The position of the unknown in the m_sysSolution vector/matrix
    int UnkPosSysSol = 0;
    // The (local) position of the component in the m_sysSolution vector/matrix starting from UnkPosSysSol
    int CompPosSysSol = 0;

    //Find the basis index
    index_t basis_ind = 0;
    for (index_t k = 0 ; k < unk; ++k)
        basis_ind += m_unknownDim[k];

    int szMax  = m_bases[basis_ind][p].size();
    for (index_t k = 1; k < targetDim; ++k)
    {
        if (szMax < m_bases[basis_ind + k][p].size())
            szMax = m_bases[basis_ind + k][p].size();
    }

    if (unk==0)
    {
        gsDebug << "size of Basis 0: " << m_bases[0][0].size()<< "\n";
        gsDebug << "Degree of Basis 0: " << m_bases[0][0].degree(0)<< " and "<< m_bases[0][0].degree(1)<<"\n";
        gsDebug << "size of Basis 1: " << m_bases[1][0].size()<< "\n";
        gsDebug << "Degree of Basis 1: " << m_bases[1][0].degree(0)<< " and "<< m_bases[1][0].degree(1)<<"\n";
        gsDebug << "size of Basis 2: " << m_bases[2][0].size()<< "\n";
        gsDebug << "Degree of Basis 2: " << m_bases[2][0].degree(0)<< " and "<< m_bases[2][0].degree(1)<<"\n";
    }

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
    gsMatrix<T> coeffs(szMax, targetDim);

    for (index_t k = 0; k < targetDim; ++k)
    {
        for (index_t i = 0; i < m_bases[basis_ind + k][p].size(); ++i)
        {
            if ( m_dofMapper[k+compNr]->is_free(i, p) ) // internal or interface
            {
                if (hasMatrixRhs) //Specially handle the if using matrix right hand side
                    coeffs(i,k) = data(m_dofMapper[k+compNr]->index(i, p),k);
                else
                {
                    coeffs(i,k) = data(UnkPosSysSol + CompPosSysSol + m_dofMapper[k+compNr]->index(i, p),0);
                }
            }
            else // eliminated DoFs: fill with Dirichlet data
            {
                coeffs(i,k) = m_fixedDofs[unk]( m_dofMapper[k+compNr]->bindex(i, p),k);
            }
            // Increase position with the size of the component

        }
        CompPosSysSol += m_dofMapper[k+compNr]->freeSize();
    }

    //If reconstructing the pressure:
    if (unk == 1)
        return m_bases[basis_ind][p].makeGeometry( give(coeffs) );
    //If reconstruction the velocity:
    else if (unk == 0)
    {
        if (m_geoTrans == DIV_CONFORMING)
        {
            std::vector<gsBasis<T> *> bases_vec;
            for (index_t k = 0; k < targetDim; ++k)
            {
                bases_vec.push_back( const_cast<gsBasis<T> *>(&m_bases[basis_ind + k][p]) );
            }
            return new gsDivConSolution<T, gsBSplineBasis<T,gsKnotVector<T> > > (coeffs,m_patches[p],bases_vec);
        }
        else
            return m_bases[0][p].makeGeometry( give(coeffs) );
    }
    else
        GISMO_ERROR("Number of unknowns is wrong");
}

} // namespace gismo

