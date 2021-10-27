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
gsMappedSpline<d,T>::gsMappedSpline( const gsMappedBasis<d,T> & mbases, const gsMatrix<T> & coefs )
:
m_coefs(coefs)
{
    m_mbases=mbases.clone().release();
    init(mbases);
}

template<short_t d,class T>
gsMappedSpline<d,T>::gsMappedSpline( const gsMappedSpline& other )
:
m_coefs(other.m_coefs)
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

    index_t n = targetDim();
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
        result.resize( n,u.cols());
        // This loop enables that the number of actives can be different for each column in u
        for (index_t k = 0; k!=u.cols(); k++)
        {
            m_mbases->active_into(patch,u.col(k),actives);
            m_mbases->eval_into(patch,u.col(k),evals);
            m_mbases->getBase(patch).linearCombination_into(m_coefs,actives,evals,tmp);
            result.col(k) = tmp;
        }
    // }
}

template<short_t d,class T>
void gsMappedSpline<d,T>::deriv_into(const unsigned patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    gsMatrix<index_t> actives;
    gsMatrix<T> evals;

    index_t n = targetDim();
    index_t m = domainDim();
    result.resize( n * m,u.cols());

    gsMatrix<T> tmp;
    // This loop enables that the number of actives can be different for each column in u
    for (index_t k = 0; k!=u.cols(); k++)
    {
        m_mbases->active_into(patch,u.col(k),actives);
        m_mbases->deriv_into(patch,u.col(k),evals);
        m_mbases->getBase(patch).linearCombination_into(m_coefs,actives,evals,tmp);
        result.col(k) = tmp;
    }
}

template<short_t d,class T>
void gsMappedSpline<d,T>::deriv2_into(const unsigned patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    gsMatrix<index_t> actives;
    gsMatrix<T> evals;

    index_t n = targetDim();
    index_t m = domainDim();
    index_t S = n*(n+1)/2;
    result.resize( S * m,u.cols());
    gsMatrix<T> tmp;
    // This loop enables that the number of actives can be different for each column in u
    for (index_t k = 0; k!=u.cols(); k++)
    {
        m_mbases->active_into(patch,u.col(k),actives);
        m_mbases->deriv_into(patch,u.col(k),evals);
        m_mbases->getBase(patch).linearCombination_into(m_coefs,actives,evals,tmp);
        result.col(k) = tmp;
    }

    // gsMatrix<index_t> actives;
    // m_mbases->active_into(patch,u,actives);
    // gsMatrix<T> evals;
    // m_mbases->deriv2_into(patch,u,evals);
    // m_mbases->getBase(patch).linearCombination_into(m_coefs,actives,evals,result);
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
    index_t S = N*(N+1)/2;

    std::vector<index_t> blocksizes(3);
    blocksizes[0] = 1;
    blocksizes[1] = N;
    blocksizes[2] = S;

    gsMatrix<T> tmp;
    for( int i = 0; i <= n; i++)
    {
        result[i].resize(blocksizes[i] * m,u.cols());
        // This loop enables that the number of actives can be different for each column in u
        for (index_t k = 0; k!=u.cols(); k++)
        {
            m_mbases->active_into(patch,u.col(k),actives);
            m_mbases->evalAllDers_into(patch,u.col(k),n,evals);
            m_mbases->getBase(patch).linearCombination_into(m_coefs,actives,evals[i],tmp);
            result[i].col(k) = tmp;

        }
    }
}

template<short_t d,class T>
gsMultiPatch<T> gsMappedSpline<d,T>::exportToPatches() const
{
    GISMO_NO_IMPLEMENTATION;
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
    GISMO_NO_IMPLEMENTATION;

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
