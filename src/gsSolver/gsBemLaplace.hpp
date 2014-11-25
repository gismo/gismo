/** @file gsBemLaplace.hpp

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Falini, A. Mantzaflaris
*/


#pragma once

#include <gsCore/gsBasis.h>
#include <gsCore/gsGeometry.h>
#include <gsModeling/gsCurveLoop.h>
#include <gsModeling/gsPlanarDomain.h>
#include <gsCore/gsFunctionExpr.h>
#include <gsSolver/gsBemUtils.h>
#include <gsAssembler/gsGaussAssembler.h>
#include <gsUtils/gsQuadrature.h>
#include <gsSolver/gsBemUtils.h>

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsBSplineBasis.h>

namespace gismo
{

template<class T>
gsBemLaplace<T>::gsBemLaplace( gsPlanarDomain<T> * pl,
                               std::vector<gsFunction<T> *> boundaries,
                               std::vector< gsBasis<T> *> basis )
{
    //  GISMO_ASSERT( ) ;
    m_boundary_fun = boundaries;
    m_pdomain= pl;
    m_basis=basis;
}


template<class T>
gsBemSolution<T> *
gsBemLaplace<T>::solve(bool parametric_bc)
{
    const int nloops = m_pdomain->numLoops();

    std::vector<gsGeometry<T> *>  Fluxf     ;
    std::vector<gsGeometry<T>*>   all_loops ;
    std::vector< std::vector<T> > all_breaks;
    Fluxf     .reserve(nloops);
    all_loops .reserve(nloops);
    all_breaks.reserve(nloops);


    int n_v;
    gsMatrix<T> gr_points, col_points;
    gsMatrix<unsigned> act;

    gsVector<T> unormal;


    // Generate quadrature grids
    int numDofs_s = 0, numDofs_v = 0;
    for(int v=0; v != nloops; v++)
    {
        // Get a single curve representing every loop
        all_loops.push_back(  m_pdomain->loop(v).singleCurve() );


        // Periodic
        dynamic_cast<gsBSplineBasis<T>*>(m_basis[v])->setPeriodic();
        

        numDofs_v +=  m_basis[v]->size();

        // Add the collocation points to the breaks, to avoid singuarities
        std::vector<T> breaks;
        m_basis[v]->anchors_into(gr_points);
        breaks = m_basis[v]->domain()->breaks();

        for ( index_t i= 1; i<gr_points.size()-1; ++i )
            breaks.push_back( gr_points(0,i) );

        std::sort(breaks.begin(), breaks.end() ) ;

        typename std::vector<T>::iterator itr =
            unique( breaks.begin(), breaks.end(), math::template almostEqual<T> );

        breaks.resize( itr - breaks.begin() ) ;

        all_breaks.push_back(breaks);

//<<<<<<< .mine
//        gsMatrix<T> * ngrid = new gsMatrix<T>;
//        gsVector<T> * wgrid = new gsVector<T>;
//        iteratedGaussRule(*ngrid, *wgrid,3* m_basis[v]->minDegree(), breaks ) ;
//        all_ngrid.push_back(ngrid);
//        all_wgrid.push_back(wgrid);
//=======

    }


    gsGaussRule<T> QuRule   ;
    gsMatrix<T>    quNodes  ;
    gsVector<T>    quWeights;
    gsVector<index_t> numNodes(1);


    // Initialize matrix and right-hand side
    gsMatrix<T> globalA(numDofs_v, numDofs_v ) ;
    gsMatrix<T> globalB( numDofs_v, 1 ) ;
    globalA.setZero();
    globalB.setZero();

    gsMatrix<T> ev, green_val, green_grad, bd_val;

    // For all loops
    numDofs_v = 0;
    for(int v = 0; v < nloops; v++)
    {
        // Greville points
        m_basis[v]->anchors_into(gr_points) ;

        // Collocation points
        all_loops[v]->eval_into( gr_points, col_points ) ;

        n_v = m_basis[v]->size();

        numDofs_s = 0;
        for(int s=0; s < nloops; s++)
        {
            // Point to the current blocks
            typename gsMatrix<T>::Block A = globalA.block(numDofs_v, numDofs_s,
                                            n_v, m_basis[s]->size() );

            typename gsMatrix<T>::Block B = globalB.block(numDofs_v, 0, n_v, 1 );

//<<<<<<< .mine
//                // Evaluator on the current loop
//                typename gsGeometry<T>::Evaluator loop_s_eval (
//                    all_loops[s]->evaluator(GEO_NEED_VALUE   |
//                                            GEO_NEED_JACOBIAN) );
//=======
            // Evaluator on the current loop
            typename gsGeometry<T>::Evaluator loop_s_eval (
                all_loops[s]->evaluator(NEED_VALUE   |
                                        NEED_JACOBIAN) );

            numNodes(0) = m_basis[s]->minDegree()+1;
            QuRule.setNodes(numNodes);


            // Point to the current sub-elements
            const std::vector<T> & breaks = all_breaks[s];
            
            // Fill in block (v,s)
            //for ( int i = 0; i != n_v -1 ; ++i)
            for ( int i = 0; i != n_v; ++i)
            {
                // Set collocation point in the Green function
                m_green_fun.setSourcePoint(col_points.col(i) ) ;                


                // Last row expresses the condition that first and last
                // Greville points coincide

                    if (v==s)
                    {
                        m_boundary_fun[v]->eval_into( parametric_bc ?
                        gr_points.col(i) : col_points.col(i), bd_val );

                        B( i, 0 ) -= bd_val(0,0) / T(2.0);
                    }

                    // Start iteration over sub-elements
                    for (std::size_t e=1; e!= breaks.size(); ++e)
                    {
                        // Map the Quadrature rule to the element
                        QuRule.mapTo(breaks[e-1], breaks[e], quNodes, quWeights );

                        // Evaluate loop values and derivatives
                        loop_s_eval->evaluateAt(quNodes);

                        // Evaluate the boundary function
                        m_boundary_fun[s]->eval_into( parametric_bc ?
                        quNodes : loop_s_eval->values(), bd_val );

                        // Evaluate the Green function
                        m_green_fun.eval_into( loop_s_eval->values(), green_val );

                        // Evaluate the Green normal derivative
                        m_green_fun.grad_into( loop_s_eval->values(), green_grad );


                        // evaluate basis functions
                        m_basis[s]->eval_into  (quNodes, ev );

                        // indices of evaluated basis functions
                        m_basis[s]->active_into(quNodes, act);

                        // For all quadrature points in this interval
                        for (int k=0; k!= numNodes[0]; ++k)
                        {
                            // Compute the outer normal
                            loop_s_eval->normal(k, unormal);

                            const T intElement = unormal.norm();
                            unormal.normalize();

                            // Evaluate the Green normal derivative
                            const T greenNDeriv = ( green_grad.col(k).transpose() * unormal ).value();

                            // Update right hand side
                            B( i, 0 )    += quWeights[k] * bd_val(0,k) * greenNDeriv * intElement;

                            const T weight = quWeights[k] *  green_val(0,k) * intElement;
                            // Update BEM collocation matrix
                            for (index_t j=0; j!=act.rows(); ++j)
                                A( i , act(j,k) ) +=   weight * ev(j,k);
                        }
                     }

                }
                numDofs_s += m_basis[s]->size();
            }
            numDofs_v += n_v;
        }

        // Cleanup loops
        freeAll( all_loops );

        globalB = globalA.colPivHouseholderQr().solve( globalB );

        numDofs_v= 0;
        for(int  v=0; v< nloops; v++)
        {
            n_v = m_basis[v]->size();
            gsMatrix<T> B = globalB.block(numDofs_v, 0, n_v, 1 );
            Fluxf.push_back( m_basis[v]->makeGeometry(give(B)) );
            numDofs_v += n_v;
        }

        // gsDomain: should have a "sample(int)" and return points in the domain
        // also boundary? or not!
        //( knotvector, planardomain, tensordomain, hierarchical domain )

        //solvers: return void, and provide an eval() -- ie inherit from gsFunction ?
        //assertions
        return new gsBemSolution<T>(m_pdomain, Fluxf, m_boundary_fun,
                                    all_breaks, parametric_bc);

    }


    };// namespace gismo
