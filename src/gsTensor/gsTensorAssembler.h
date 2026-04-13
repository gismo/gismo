/** @file gsTHBAssembler.h

    @brief Implements matrix assembly of common matrices by use of low
    tensor-rank approximation and Hadamar or Kronecker-product sum representation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, M. Matucci
*/

#pragma once

#include <gsTensor/gsTensorTools.h>
#include <gsTensor/gsKroneckerMatrix.h>
#include <gsTensor/gsTensorFunction.h>
#include <gsMatrix/gsFiberMatrix.h>


namespace gismo
{


enum GalerkinMatrix
{
    MASS      = 1U << 1,
    STIFFNESS = 1U << 2
};


template <class T = real_t>
class gsTensorAssembler
{
    const gsBasis<T> * m_basis;
    typedef gsFiberMatrix<T,ColMajor> FiberMatrix;
    typedef gsKroneckerMatrix<T,-1> Kronecker_t;
    typedef typename Kronecker_t::uPtr Kronecker_ptr;
    typedef typename Kronecker_t::cwiseMat cwiseMat;
    //typedef typename Kronecker_t::blockRowColumnRep SkeletonMat;

    Kronecker_ptr m_kprod;

    gsOptionList m_option;

public:

    gsTensorAssembler() : m_basis(nullptr), m_option(defaultOptions()) {}

    //explicit gsTensorAssembler(const gsBasis<T> & basis) : m_kprod(&basis), m_option(defaultOptions()) {}

    const Kronecker_t & kronecker() const { return *m_kprod; }

    static gsOptionList defaultOptions();

    gsOptionList & options() { return m_option; }

public:

    void compute(const gsBasis<T> & bb, GalerkinMatrix mode,
                 const gsTensorFunction<T> & kernel )
    {
        if (const gsTensorBSplineBasis<2,T>* hb2 = dynamic_cast<const gsTensorBSplineBasis<2,T>*>(&bb))
            compute(*hb2, mode, kernel);
        else if (const gsTensorBSplineBasis<3,T>* hb3 = dynamic_cast<const gsTensorBSplineBasis<3,T>*>(&bb))
            compute(*hb3, mode, kernel);
        else
            GISMO_ERROR("Unable to process "<< bb);
    }

    void compute() { GISMO_ASSERT(m_basis!=nullptr, "No basis given."); compute(*m_basis); }

    template<short_t d>
    void compute(const gsTensorBasis<d,T> & hb, GalerkinMatrix mode,
                 const gsTensorFunction<T> & kernel = gsTensorFunction<T>(d))
    {
        m_basis = & hb;

        cwiseMat  skeleton_mat;
        if (mode==MASS)
            assembleMass<d>(kernel, skeleton_mat);
        if (mode==STIFFNESS)
            assembleStiffness<d>(kernel, skeleton_mat);

        const index_t rank = kernel.totalRank();
        m_kprod.reset( new gsKroneckerMatrix<T,d>(give(skeleton_mat)) );
    }

    template<short_t d>
    void assembleMass(const gsTensorFunction<T> & kernel,
                      cwiseMat & skeleton)
    {
#ifdef DO_TIMING
        timer.restart();
#endif
        const index_t r = kernel.rank(0);
        const gsTensorBSplineBasis<d,T> & hb = *dynamic_cast<const gsTensorBSplineBasis<d,T>*>(m_basis);
        gsExprAssembler<> A(1, 1);

        skeleton.reserve(d);

        for (short_t s=0; s!=d; ++s)
        {
            A.cleanUp();
            const gsBSplineBasis<T> & hbs = hb.component(s);
            A.setIntegrationDomain( hbs.domain() );

            auto ul = A.getSpace(hbs, r);
            auto vk = A.getTestSpace(hbs, 1);
            auto w  = A.getCoeff(kernel.component(s));
            A.initMatrix();
            A.assemble( vk * w.tr() * ul.tr() );

            skeleton.push_back( A.giveFiberMatrix() );
        }

#ifdef DO_TIMING
        assembly += timer.stop();
#endif
    }

    template<short_t d>
    void assembleStiffness(const gsTensorFunction<T> & tf,
                           cwiseMat & skeleton)
    {
#ifdef DO_TIMING
        timer.restart();
#endif
        const gsTensorBSplineBasis<d,T> & basis = *dynamic_cast<const gsTensorBSplineBasis<d,T>*>(m_basis);
        const gsVector<int> & ranks  = tf.ranks(); // Ranks of stiffness

        std::vector<gsMatrix<T> > ev_b;
        gsMatrix<T>  qNodes, ev_gf, localMat;
        gsVector<T> qWeights;
        gsMatrix<index_t> actives;

        skeleton.resize(d);
        const int Grank = ranks.sum();

        for (int k = 0; k!=d; ++k)
        {
            const gsBSplineBasis<T> & b = basis.component(k);
            const gsFunction<T> & Gk = tf.component(k);

            const int sz = b.size();
            gsFiberMatrix<T> & mat_k = skeleton[k];
            const index_t msz = Grank * sz;
            mat_k.resize(sz,msz);
            mat_k.setZero();
            mat_k.reserve( gsVector<size_t>::Constant(msz, 2*b.degree(0)+1) );

            const index_t numNodes = b.degree(0) + 1;
            gsGaussRule<T> QuRule(numNodes);

            // Make domain element iterator
            typename gsBasis<T>::domainIter domIt = b.domain()->beginAll();
            typename gsBasis<T>::domainIter domItEnd = b.domain()->endAll();

            // Start iteration over elements
            for (; domIt!=domItEnd; ++domIt)
            {
                b.active_into(domIt.centerPoint() , actives);
                const index_t numActive = actives.rows();

                // Compute the quadrature rule on the element
                QuRule.mapTo(domIt.lowerCorner().value(),
                            domIt.upperCorner().value(),
                            qNodes, qWeights );

                b.evalAllDers_into(qNodes, 1, ev_b );
                const gsMatrix<T> & ev = ev_b[0];
                const gsMatrix<T> & grad = ev_b[1];

                Gk.eval_into(qNodes, ev_gf); // d^2 ?
                int pos = 0;
                for (int r = 0; r < d; ++r)
                    for (int s = 0; s < d; ++s)
                    //for (int s = 0; s <= r; ++s)
                    {
                        const int rs = r*d+s;
                        for (index_t t = 0; t < ranks[rs]; ++t) // for each rank-1 part
                        {
                            localMat.setZero(numActive, numActive);

                            for (index_t j = 0; j < numNodes; ++j) // for all quadrature nodes
                            {
                                localMat.noalias() += qWeights[j] *
                                    ev_gf(pos,j)* // Redundant computations: no symmetry of bilinear form
                                    // dd, de, ed, dd
                                    (
                                    ( r==k ? grad.col(j)             : ev.col(j)             ) *
                                    ( s==k ? grad.col(j).transpose() : ev.col(j).transpose() )
                                    );
                            }

                            // Local to global
                            const index_t qq = sz*pos;
                            for (index_t i = 0; i < numActive; ++i)
                            {
                                const int ii = actives(i,0);
                                for (index_t j = 0; j < numActive; ++j)
                                {
                                    const int jj = actives(j,0);
                                    mat_k.coeffRef(ii, qq+jj) += localMat(i, j);
                                }
                            }
                            ++pos;
                        }//t
                }
            }
        }
#ifdef DO_TIMING
        assembly += timer.stop();
#endif
    }

/*
    template<short_t d>
    void assembleStiffness( const gsTensorFunction<T> & kernel,
                              const std::vector<gsCompressedSparseFiber<d,index_t> > & CM,
                              cwiseMat & skeleton)
    {
#ifdef DO_TIMING
        timer.restart();
#endif
        const gsTensorBasis<d,T> & hb = *dynamic_cast<const gsTensorBasis<d,T>*>(m_basis);
        gsExprAssembler<> A(1, d);//total rank

        std::array<typename gsExprAssembler<>::space,d> u;
        for (short_t s=0; s!=d; ++s)
        {
            auto u_q  = A.getSpace(hbs, s, r);
            auto v_k  = A.getTestSpace(hbs, 1);
            auto w_kq = A.getCoeff(kernel.component(s));
        }
        
        skeleton.resize(d);
        index_t r;
        for (short_t s=0; s!=d; ++s)
        {
            gsHBasisComponent<T> hbs(hb, CM, s);
            A.setIntegrationDomain( hbs.domain() );
            for (short_t k=0; k!=d; ++k)
                for (short_t q=0; q!=d; ++q)
                {
                    r = kernel.rank(k*d + q);
                    A.cleanUp();

                    auto u_q = A.getSpace(hbs, s, r); //r,0
                    auto v_k = A.getTestSpace(hbs, 1);
                    auto w_kq  = A.getCoeff(kernel.component(s));
                    A.initMatrix();
            
                    switch ( (k==s) + 2 *(q==s) )
                    {
                    default://0
                        A.assemble(v_k * w_kq.tr() * u_q.tr());
                        break;
                    case 1:
                        A.assemble(grad(v_k)* w_kq.tr() * u_q.tr());
                        break;
                    case 2:
                        A.assemble(v_k * w_kq.tr() * grad(u_q).tr());
                        break;
                    case 3:
                        A.assemble(grad(v_k)* w_kq.tr() * grad(u_q).tr());
                        break;
                    }
                }
        }
            skeleton.push_back( A.giveFiberMatrix() );


#ifdef DO_TIMING
            setup += timer.stop();
            timer.restart();
#endif
    }
*/

};

template<class T>
gsOptionList gsTensorAssembler<T>::defaultOptions()
{
    gsOptionList opt;
    opt.addSwitch("highestIntegrationLevel", "Set integration domain to the highest level", false);
    return opt;
}


} // namespace gismo
