/** @file gsTensorFunction.h

    @brief Provides declaration of TensorFunction class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsTensor/gsTensorUtils.h>

#pragma once


namespace gismo
{

/**
    \brief This class provides a tensor function represented by
    skeleton (component) functions with certain rank per component.

    f = sum

    \ingroup function
    \ingroup Tensor
*/
template<class T = real_t>
class gsTensorFunction : public gsFunction<T>
{
    std::vector< const gsFunction<T>*> m_skeleton;
    gsVector<int>                      m_ranks;

public:

    /*
    gsTensorFunction(const gsFunction<T> & f1,
                     const gsFunction<T> & f2,
                     const gsVector<int> & ranks)
    : m_ranks(ranks)
    {
        m_skeleton.push_back(f1.clone());
        m_skeleton.push_back(f2.clone());
    }
    */
    
    gsTensorFunction(short_t _dim = 2, int type = 0) : m_skeleton(_dim)
    {
        // /*
        const index_t sz = (0==type ? 1 : _dim);
        gsMatrix<T> mm = gsMatrix<T>::Identity(sz,sz);
        mm.resize(sz*sz,1);
        m_ranks = mm.template cast<int>();
        mm.setOnes(sz,1);
        for (size_t i=0; i!= m_skeleton.size(); ++i)
            m_skeleton[i] = new gsConstantFunction<T>(mm,1);
        //*/
        
        /*
        int rr = 3;
        gsVector<T> cc(rr);
        cc.setConstant( math::pow((real_t)rr, -1.0/_dim) );
        m_ranks[0] = rr;
        for (size_t i=0; i!= m_skeleton.size(); ++i)
            m_skeleton[i] = new gsConstantFunction<T>(cc, 1);
        //*/
    }

    gsTensorFunction(const gsTensorFunction & other)
    {
        this->operator=(other);
    }

    /// Assignment operator
    gsTensorFunction& operator= (const gsTensorFunction& other )
    {
        m_ranks = other.m_ranks;
        cloneAll(other.m_skeleton.begin(), other.m_skeleton.end(), m_skeleton.begin() );
        return *this;
    }

    /// Move assignment operator
    gsTensorFunction& operator= ( gsTensorFunction&& other )
    {
        m_ranks = give(other.m_ranks);
        m_skeleton = give(other.m_skeleton);
        other.m_skeleton.clear();
        return *this;
    }


    /// \a funcs are consumed
    explicit gsTensorFunction(const std::vector<gsFunction<T>*> & funcs)
    {
        m_ranks.setZero();// ! ++
        m_skeleton.swap(funcs);
    }

    // Decomposes Afap as sum of univariate products
    gsTensorFunction(const gsGeometry<T> & Afap, T tol, int maxrank = -1);

    // todo: Projection + decomposition
    /*
    gsTensorFunction(const gsFunction<T> & Afap, const gsBasis<T> & prBasis, T tol)
    :
    m_ranks(Afap.targetDim())
    {
    // Evaluate on anchors

    // Decompose

    // Interpolate skeleton
    }
    */

    // gsTensorFunction(const gsBasis<T> & basis, const std::vector<gsMatrix<T> > & skeleton)
    // : // Assumes single scalar function
    // m_ranks(1)
    // {
    //     const index_t d = skeleton.size();
    //     GISMO_ASSERT( basis.dim() == d, "Inconsistent skeleton");
    //
    //     m_ranks[0] = skeleton.front().cols();
    //
    //     m_skeleton.resize(d);
    //     for (index_t i=0; i!= d; ++i)
    //         m_skeleton[i] = basis.component(i).makeGeometry( skeleton[i] );
    //
    //
    // }

    gsTensorFunction(const gsBasis<T> & basis, const std::vector<gsMatrix<T> > & skeleton)
    : m_ranks(1)
        {
            const index_t d = skeleton.size();
            GISMO_ASSERT(basis.dim() == d, "Inconsistent skeleton");
            GISMO_ASSERT(d > 0, "Skeleton must not be empty.");

            const index_t R = skeleton.front().cols();
            GISMO_ASSERT(R > 0, "Skeleton rank must be positive.");

            for (index_t i = 0; i != d; ++i)
            {
                GISMO_ASSERT(skeleton[i].cols() == R,
                    "All skeleton components must have the same number of columns.");
            }

            m_ranks[0] = R;

            m_skeleton.resize(d);
            for (index_t i = 0; i != d; ++i)
                m_skeleton[i] = basis.component(i).makeGeometry(skeleton[i]).release();
        }

    gsTensorFunction(const gsBasis<T> & basis,
                     const std::vector<gsMatrix<T> > & skeleton,
                     gsVector<index_t> ranks)
    : m_ranks(ranks)
    {
        const index_t d = skeleton.size();
        GISMO_ASSERT( basis.dim() == d, "Inconsistent skeleton");
        GISMO_ASSERT( ranks.sum() == skeleton.front().cols(), "Inconsistent rank/skeleton");

        m_skeleton.resize(d);
        for (index_t i=0; i!= d; ++i)
            m_skeleton[i] = basis.component(i).makeGeometry( give(skeleton[i]) ).release();
    }

    ~gsTensorFunction()
    {
        freeAll(m_skeleton);
    }

    std::ostream &print(std::ostream &os) const
    {
        os << "Tensor function, rank="<<m_ranks.transpose()<<".\n";
        for ( size_t i = 0; i!= m_skeleton.size(); ++i)
            os << "Skeleton function "<<i<<" ("<<m_skeleton[i]->domainDim()<<"):\n"<<*m_skeleton[i]<<"\n";

        return os;
    }


public:

    int domainDim() const
    {
        // Assuming components are 1D
        //return m_skeleton[0]->domainDim() + m_skeleton[1]->domainDim();
        return m_skeleton.size();
    }

    int targetDim() const
    {
        return m_ranks.size();
    }

    int rank( int i) const { return m_ranks[i]; }

    const gsVector<int> & ranks() const { return m_ranks; }

    int totalRank() const { return m_ranks.sum(); }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        const index_t n = m_ranks.size();
        const index_t d = u.rows(); //Dimension of source space
        result.resize(n, u.cols() );

        std::vector<gsMatrix<T> > ev(d);
        for (index_t j=0; j!= u.cols(); ++j) //eval points
        {
            for (index_t k=0; k!= d; ++k)    //Dimension of source space/Number of univariate functions
                m_skeleton[k]->eval_into(u.col(j).row(k), ev[k]);

            // Hadamard product
            gsMatrix<T> & Hprod  = ev.front();
            for (index_t i=1; i!=d; ++i)
                Hprod.array() *= ev[i].array();

            index_t pos = 0;
            for (index_t k=0; k!=n; ++k)
            {
                result(k,j) = Hprod.middleRows(pos, m_ranks[k]).sum();
                pos += m_ranks[k];
            }
        }
    }

    // Part in variable \a s of Function coordinate \a k
    void evalComponent_into(int k, int s, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(m_ranks[k], u.cols() );
        GISMO_ASSERT( u.rows() ==1, "expect 1d pts");

        gsMatrix<T> ev;
        const int pos = m_ranks.segment(1,k).sum();// coord. position

        for (index_t j=0; j!= u.cols(); ++j)
        {
            m_skeleton[s]->eval_into(u.col(j).row(s), ev);
            result.col(j) = ev.middleRows(pos, m_ranks[k]);
        }
    }

    // Part in variable \a k of Function coordinate \a i
    const gsFunction<T> & component(int k) const
    {
        return *m_skeleton[k];
    }
};



template<class T>
gsTensorFunction<T>::gsTensorFunction(const gsGeometry<T> & Afap, T tol, int maxrank)
{
    const gsBasis<T> & bb = Afap.basis();
    //gsDebugVar( bb);
    gsVector<index_t> sz;

    if ( const gsTensorBasis<2,T>* tb2 = // to do: improve access
         dynamic_cast<const gsTensorBasis<2,T>*>( & Afap.basis()))
        tb2->size_cwise(sz);
    else if ( const gsTensorBasis<3,T>* tb3 = // to do: improve access
         dynamic_cast<const gsTensorBasis<3,T>*>( & Afap.basis()))
        tb3->size_cwise(sz);
    else if ( const gsTensorBasis<4,T>* tb4 = // to do: improve access
         dynamic_cast<const gsTensorBasis<4,T>*>( & Afap.basis()))
        tb4->size_cwise(sz);

    GISMO_ASSERT( is_finite  (Afap.coefs()), "Has infs" );
    GISMO_ASSERT( is_nan_free(Afap.coefs()), "Has nans" );

    std::vector< gsMatrix<T> > skeleton;

    const index_t d = sz.size();
    if ( d==2)
    {
        //--Test
        //gsCpDec<0>(Afap.coefs(), sz, skeleton, m_ranks, tol, 50);
        //--Test
        gsSvdDecomposition(Afap.coefs(), sz, m_ranks, skeleton, tol, maxrank);
    }
    else
    {
        // ALS with incremental rank-1 steps
        //gsCpDec<0>(Afap.coefs(), sz, skeleton, m_ranks, tol, 50);

        // ALS with incremental rank-1 steps and Tucker compression
        gsCpDec<3>(Afap.coefs(), sz, skeleton, m_ranks, tol, 100, maxrank); //Maximal rank for each component
    }

    GISMO_ASSERT( 0 != totalRank(), "The approximation is NULL for tolerance "<< tol);

    // Create the skeleton functions
    m_skeleton.reserve(d);
    for (index_t i = 0; i!= d; ++i)
    {
        // gsDebugVar(skeleton[i].transpose());
        //GISMO_ASSERT( is_finite(skeleton[i]), "Has infs: "<< i );
        m_skeleton.push_back( bb.component(i).makeGeometry( give(skeleton[i])).release() );
    }
}


}// namespace gismo
