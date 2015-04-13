/** @file gsTensorFunction.h

    @brief Provides declaration of TensorFunction class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once


namespace gismo
{

/** 
    \brief This class provides a tensor function represented by
    skeleton (component) functions with certain rank per component.

    \ingroup function
    \ingroup Tensor
*/
template<class T>
class gsTensorFunction : public gsFunction<T>
{
private:
    gsTensorFunction( ){ } 

public:

    gsTensorFunction(const gsFunction<T> & f1, 
                     const gsFunction<T> & f2,
                     const gsVector<int> & ranks)
    : m_ranks(ranks)
    { 
        m_skeleton.push_back(f1.clone());
        m_skeleton.push_back(f2.clone());
    } 
    
    explicit gsTensorFunction(const std::vector<gsFunction<T>*> & funcs)
    : m_ranks(ranks)
    {
        cloneAll(funcs.start(), funcs.end(), m_skeleton.begin() );
    }

    gsTensorFunction(const gsGeometry<T> & Afap, T tol);

    ~gsTensorFunction( )
    { 
        freeAll(m_skeleton);
    } 

public:

    int domainDim() const                     
    {
        return m_skeleton[0]->domainDim() + m_skeleton[1]->domainDim();
    }
    
    int targetDim() const                     
    {
        return m_ranks.size();
    }
    
    int rank( int i) const { return m_ranks[i]; }

    const gsVector<int> & ranks() const { return m_ranks; }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT( u.rows() == 2, "wrong dimension" );
        const index_t n = m_ranks.size();
        result.resize(n, u.cols() );

        gsMatrix<T> ev1, ev2;
        for (index_t j=0; j!= u.cols(); ++j)
        {
            m_skeleton[0]->eval_into(u.col(j).row(0), ev1);
            m_skeleton[1]->eval_into(u.col(j).row(1), ev2);            

            int pos = 0;
            for (index_t k=0; k!=n; ++k)
            {
                result(k,j) = (ev1.middleRows(pos, m_ranks[k]).transpose() * 
                               ev2.middleRows(pos, m_ranks[k])
                              ).sum();
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
            m_skeleton[s]->eval_into(u.col(j), ev);
            result.col(j) = ev.middleRows(pos, m_ranks[k]);
        }
    }

    // Part in variable \a k of Function coordinate \a i
    const gsFunction<T> & component(int k) const
    {
        return *m_skeleton[k];
    }
    
protected:

    void initBySVD(const gsTensorBasis<2,T> & tb, const gsMatrix<T> & AA, T tol);

protected:

    std::vector< const gsFunction<T>*> m_skeleton;

    gsVector<int>                      m_ranks;
};



template<class T>
gsTensorFunction<T>::gsTensorFunction(const gsGeometry<T> & Afap, T tol)
{
    if ( const gsTensorBasis<2,T>* tb2 = 
         dynamic_cast<const gsTensorBasis<2,T>*>( & Afap.basis() ) 
       )
    {
        initBySVD( *tb2, Afap.coefs(), tol );
    }
    else
        gsWarn<<"Not implemented yet.\n";

}


template<class T>
void gsTensorFunction<T>::initBySVD(const gsTensorBasis<2,T> & tb, const gsMatrix<T> & AA, T tol)
{
    gsVector<int,2> sz;
    tb.size_cwise(sz);

    GISMO_ASSERT( sz.prod() == AA.rows(), "Invalid tensor structure");

    const int n = sz.prod();
    //gsDebugVar( sz.transpose() );

    std::vector< gsMatrix<T> > TC(2);

    m_ranks.resize(AA.cols());

    gsMatrix<T> tmp;
    // Permutations ?

    Eigen::JacobiSVD< typename gsMatrix<T>::Base > svd;
    const unsigned svdOpts = Eigen::ComputeThinU | Eigen::ComputeThinV;

    //gsDebugVar(AA.cols() );
    //gsDebugVar(AA.rows() );
    for ( index_t i = 0; i<AA.cols(); ++i )
    {
        gsAsConstMatrix<T> Amat( AA.col(i).data(), sz[0], n / sz[0]);
        svd.compute(Amat, svdOpts);
        const index_t numSV = svd.singularValues().size();
        // gsDebugVar( svd.matrixU().rows());
        // gsDebugVar( svd.matrixU().cols());
        // gsDebugVar( svd.matrixV().rows());
        // gsDebugVar( svd.matrixV().cols());
        // gsDebugVar( svd.singularValues().size() );
        //gsDebug<<"SV:" <<  svd.singularValues().transpose() <<"\n";
        //gsDebug<<"SVsqrt:" <<  svd.singularValues().array().sqrt().transpose()<<"\n";;

        const double * svalptr = svd.singularValues().data();
        const double * pos = std::lower_bound( svalptr ,svalptr + numSV, tol, 
                                               std::greater<const T>() );
        int rank = pos - svalptr;
        //gsDebug << "Rank: "<< rank <<" (full was "<<svd.singularValues().size()<<")\n";
        gsVector<T> sqsv = svd.singularValues().topRows(rank).array().sqrt();
        gsDebug<< "Singular values: "<< svd.singularValues().topRows(rank).transpose() 
               <<" ("<< ( rank < numSV ? svd.singularValues()[rank] : 0)<<",..)\n";

        // special case of zero rank
        if (rank == 0)
        {
            rank++;
            sqsv.setZero(1);
        }
        
        // Append to matrix
        if ( i == 0)
        {
            TC[0] = svd.matrixU().leftCols(rank) * sqsv.asDiagonal();
            TC[1] = svd.matrixV().leftCols(rank) * sqsv.asDiagonal();
        }
        else// Append to matrix
        {
           tmp.resize( TC[0].rows(), TC[0].cols()+rank);
           tmp << TC[0], // 2D
               svd.matrixU().leftCols(rank)*sqsv.asDiagonal();
           tmp.swap(TC[0]);
            // TC[0].conservativeResize(Eigen::NoChange, TC[0].cols()+rank);
            // TC[0].rightCols(rank)= svd.matrixU().leftCols(rank)*sqsv.asDiagonal();

            tmp.resize( TC[1].rows(), TC[1].cols()+rank);
            tmp << TC[1], // 2D
                svd.matrixV().leftCols(rank)*sqsv.asDiagonal() ;
            tmp.swap(TC[1]);
            // TC[1].conservativeResize(Eigen::NoChange, TC[1].cols()+rank);
            // TC[1].rightCols(rank)= svd.matrixV().leftCols(rank)*sqsv.asDiagonal();
        }
        m_ranks[i] = rank;

        gsDebug<< "Tolerance:"<< tol <<"\n";
        gsDebug<< "Low-rank Approx. error (rank "<<rank
               <<", out of "<<svd.singularValues().size()<<"): "
        << ( TC[0].rightCols(rank) * TC[1].rightCols(rank).transpose() - Amat ).norm() <<"\n";
    }

    // Skeleton functions
    m_skeleton.push_back( tb.component(0).makeGeometry(TC[0]) );
    m_skeleton.push_back( tb.component(1).makeGeometry(TC[1]) );
}



}// namespace gismo

