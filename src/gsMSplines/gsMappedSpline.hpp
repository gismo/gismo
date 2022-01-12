/** @file gsMappedSpline.hpp

    @brief implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsMSplines/gsMappedSpline.h>

namespace gismo
{

// template<short_t d,class T>
// gsMappedSpline<d,T>::gsMappedSpline( gsMultiPatch<T> const & mp,std::string pathToMap ) : Base()
// {
//     short_t geoDim = mp.geoDim();
//     std::vector<gsMatrix<T> * > coefs;
//     for(int i = 0;i<mp.size();++i)
//         coefs.push_back( new gsMatrix<T>(mp.patch(i).coefs() ) );
//     m_compBasis = new gsMappedBasis<d,T>(mp,pathToMap);
//     Base::m_basis=new gsMappedSingleBasis<d,T>(m_compBasis->getMappedSingleBasis(0));
//     int start = 0, end = -1;
//     gsMatrix<T> localCoefs;
//     localCoefs.resize(m_compBasis->localSize(),geoDim);
//     for(size_t i = 0;i<mp.nPatches();i++)
//     {
//         start=end+1;
//         end+=m_compBasis->localSize(i);
//         localCoefs.block(start,0,end-start+1,geoDim) << *(coefs[i]);
//     }
//     m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
//     freeAll(coefs);
// }

// template<short_t d,class T>
// gsMappedSpline<d,T>::gsMappedSpline(const gsMultiPatch<T> & mp, const gsSparseMatrix<T> & bmap ) : Base()
// {
//     m_compBasis = new gsMappedBasis<d,T>(gsMultiBasis<T>(mp),bmap);
//     Base::m_basis=new gsMappedSingleBasis<d,T>(m_compBasis->getMappedSingleBasis(0));
//     int start = 0, end = -1;
//     short_t geoDim = mp.geoDim();
//     gsMatrix<T> localCoefs(m_compBasis->localSize(),geoDim);
//     for(size_t i = 0;i<mp.nPatches();i++)
//     {
//         start=end+1;
//         end+=m_compBasis->localSize(i);
//         localCoefs.block(start,0,end-start+1,geoDim) << mp.patch(i).coefs();
//     }
//     m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
// }

template<short_t d,class T>
gsMappedSpline<d,T>::gsMappedSpline( const gsMultiPatch<T> & mp, const gsSparseMatrix<T> & m )
{
    GISMO_ASSERT(mp.nPatches()>0,"MultiPatch is empty?");
    gsMultiBasis<T> mb(mp);
    m_mbases = new gsMappedBasis<d,T>(mb,m);

    // TRANSFORM COEFFICIENTS
    index_t rows = 0;
    index_t cols = mp.patch(0).coefs().cols();
    for (size_t p=0; p!=mp.nPatches(); ++p)
        rows += mp.patch(p).coefs().rows();

    gsMatrix<T> local;
    local.resize(rows,cols);

    index_t offset = 0;
    for (size_t p=0; p!=mp.nPatches(); ++p)
    {
        local.block(offset,0,mp.patch(p).coefs().rows(),cols) = mp.patch(p).coefs();
        offset += mp.patch(p).coefs().rows();
    }

    m_mbases->local_coef_to_global_coef(local,m_global);

    init(*m_mbases);
}

template<short_t d,class T>
gsMappedSpline<d,T>::gsMappedSpline( const gsMappedBasis<d,T> & mbases, const gsMatrix<T> & coefs )
:
m_global(coefs)
{
    m_mbases=mbases.clone().release();
    init(mbases);
}

template<short_t d,class T>
gsMappedSpline<d,T>::gsMappedSpline( const gsMappedSpline& other )
:
m_global(other.m_global)
{
    m_mbases=other.m_mbases->clone().release();
}

template<short_t d,class T>
gsMappedSpline<d,T> & gsMappedSpline<d,T>::operator=( const gsMappedSpline& other )
{
    delete m_mbases;
    m_mbases=other.m_mbases->clone().release();
    return *this;
}

template<short_t d,class T>
void gsMappedSpline<d,T>::eval_into(const unsigned patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    gsMatrix<index_t> actives;
    gsMatrix<T> evals;

    index_t N = targetDim();
    // gsC1Basis<d,t> * basis = dynamic_cast<gsC1Basis<d,t> *>(this->getBase(patch));
    // if (basis==NULL)
    // {
        // m_mbases->active_into(patch,u,actives);
        // m_mbases->eval_into(patch,u,evals);
        // m_mbases->getBase(patch).linearCombination_into(m_coefs,actives,evals,result);
    // }
    // else
    // {
        gsMatrix<T> tmp;
        result.resize( N,u.cols());
        // This loop enables that the number of actives can be different for each column in u
        for (index_t k = 0; k!=u.cols(); k++)
        {
            m_mbases->active_into(patch,u.col(k),actives);
            m_mbases->eval_into(patch,u.col(k),evals);
            m_mbases->getBase(patch).linearCombination_into(m_global,actives,evals,tmp);
            result.col(k) = tmp;
        }
    // }
}

template<short_t d,class T>
void gsMappedSpline<d,T>::deriv_into(const unsigned patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    gsMatrix<index_t> actives;
    gsMatrix<T> evals;

    index_t N = targetDim();
    index_t M = domainDim();
    result.resize( N * M,u.cols());

    gsMatrix<T> tmp;
    // This loop enables that the number of actives can be different for each column in u
    for (index_t k = 0; k!=u.cols(); k++)
    {
        m_mbases->active_into(patch,u.col(k),actives);
        m_mbases->deriv_into(patch,u.col(k),evals);
        m_mbases->getBase(patch).linearCombination_into(m_global,actives,evals,tmp);
        result.col(k) = tmp;
    }
}

template<short_t d,class T>
void gsMappedSpline<d,T>::deriv2_into(const unsigned patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    gsMatrix<index_t> actives;
    gsMatrix<T> evals;

    index_t N = targetDim();
    index_t M = domainDim();
    index_t S = M*(M+1)/2;
    result.resize( S*N,u.cols());
    gsMatrix<T> tmp;
    // This loop enables that the number of actives can be different for each column in u
    for (index_t k = 0; k!=u.cols(); k++)
    {
        m_mbases->active_into(patch,u.col(k),actives);
        m_mbases->deriv2_into(patch,u.col(k),evals);
        m_mbases->getBase(patch).linearCombination_into(m_global,actives,evals,tmp);
        result.col(k) = tmp;
    }
}

template<short_t d,class T>
void gsMappedSpline<d,T>::evalAllDers_into(const unsigned patch, const gsMatrix<T> & u,
                                             const int n, std::vector<gsMatrix<T> >& result) const
{
    result.resize(n+1);

    gsMatrix<index_t> actives;
    std::vector< gsMatrix<T> > evals;

    index_t N = targetDim();
    index_t m = domainDim();
    index_t S = m*(m+1)/2;

    std::vector<index_t> blocksizes(3);
    blocksizes[0] = N;
    blocksizes[1] = N*m;
    blocksizes[2] = N*S;

    gsMatrix<T> tmp;
    // todo: change the loop over i and the loop over k
    for( int i = 0; i <= n; i++)
    {
        result[i].resize(blocksizes[i],u.cols());
        // This loop enables that the number of actives can be different for each column in u
        for (index_t k = 0; k!=u.cols(); k++)
        {
            m_mbases->active_into(patch,u.col(k),actives);
            m_mbases->evalAllDers_into(patch,u.col(k),n,evals);
            m_mbases->getBase(patch).linearCombination_into(m_global,actives,evals[i],tmp);
            result[i].col(k) = tmp;

        }
    }

    // for( int i = 0; i <= n; i++)
    //     result[i].resize(blocksizes[i],u.cols());

    // // todo: change the loop over i and the loop over k
    // for (index_t k = 0; k!=u.cols(); k++)
    // {
    //     // This loop enables that the number of actives can be different for each column in u
    //     m_mbases->active_into(patch,u.col(k),actives);
    //     m_mbases->evalAllDers_into(patch,u.col(k),n,evals);
    //     for( int i = 0; i <= n; i++)
    //     {
    //         m_mbases->getBase(patch).linearCombination_into(m_global,actives,evals[i],tmp);
    //         result[i].col(k) = tmp;
    //     }
    // }
}

template<short_t d,class T>
gsMultiPatch<T> gsMappedSpline<d,T>::exportToPatches() const
{
    gsMatrix<T> local;
    m_mbases->global_coef_to_local_coef(m_global,local);
    return m_mbases->exportToPatches(local);

    // gsMultiPatch<T> mp;
    // gsFunction<T> * msinglesplinefun;
    // for (size_t p=0; p!=this->nPatches(); ++p)
    // {
    //     msinglesplinefun = const_cast<gsFunction<real_t> *>(dynamic_cast<const gsFunction<real_t> * >(&(m_ss[p])));
    //     typename gsGeometry<T>::uPtr geom = dynamic_cast<typename gsGeometry<T>::uPtr>(msinglesplinefun);
    //     // mp.addPatch((gsGeometry<T>::uPtr) msinglesplinefun.clone());
    // }

    // return mp;

    // GISMO_NO_IMPLEMENTATION;
    //
    // gsMatrix<T> localCoef;
    // m_compBasis->global_coef_to_local_coef(Base::m_coefs,localCoef);
    // return m_compBasis->exportToPatches(localCoef);

    /*
     template<short_t d,class T>
     gsMultiPatch<T> gsMappedBasis<d,T>::exportToPatches(gsMatrix<T> const & localCoef) const
     {
         std::vector<gsGeometry<T> *> patches(nPatches());
         for(size_t i = 0; i<nPatches() ; ++i)
             patches[i]= exportPatch(i,localCoef);
         return gsMultiPatch<T>(patches,m_topol.boundaries(),m_topol.interfaces());
     }
     */
}

template<short_t d,class T>
gsGeometry<T> * gsMappedSpline<d,T>::exportPatch(int i,gsMatrix<T> const & localCoef) const
{
    return m_mbases->exportPatch(i,localCoef);

    /*

template<short_t d,class T>
gsGeometry<T>* gsMappedBasis<d,T>::exportPatch(const int i,gsMatrix<T> const & localCoef) const
{
    const short_t geoDim=localCoef.cols();
    const int start = _getFirstLocalIndex(i);
    const int end   = _getLastLocalIndex(i);
    gsMatrix<T> coefs = localCoef.block(start,0,end-start+1,geoDim);
    return getBase(i).makeGeometry( give(coefs) ).release();
}


     */

}

}
