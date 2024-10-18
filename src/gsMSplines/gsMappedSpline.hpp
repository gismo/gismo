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

template<short_t d,class T>
gsMappedSpline<d,T>::gsMappedSpline( const gsMultiPatch<T> & mp, const gsSparseMatrix<T> & m )
{
    GISMO_ASSERT(mp.nPatches()>0,"MultiPatch is empty?");
    m_mbases = new gsMappedBasis<d,T>(gsMultiBasis<T>(mp),m);

    // collect and transform the coefficients
    gsMatrix<T> local = mp.coefs();
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
: gsFunctionSet<T>(), m_global(other.m_global)
{
    m_mbases=other.m_mbases->clone().release();
}

template<short_t d,class T>
gsMappedSpline<d,T> & gsMappedSpline<d,T>::operator=( const gsMappedSpline& other )
{
    delete m_mbases;
    m_mbases=other.m_mbases->clone().release();
    m_global = other.m_global;
    m_ss = other.m_ss;
    for (auto & s : m_ss) s.setSource(*this);
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
                                           const int n, std::vector<gsMatrix<T> >& result,
                                           bool sameElement) const
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
            m_mbases->active_into(patch,u.col(k),actives);//..
            m_mbases->evalAllDers_into(patch,u.col(k),n,evals,sameElement);
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
}

template<short_t d,class T>
gsGeometry<T> * gsMappedSpline<d,T>::exportPatch(int i,gsMatrix<T> const & localCoef) const
{
    return m_mbases->exportPatch(i,localCoef);
}

template<short_t d,class T>
std::map<index_t, internal::ElementBlock> gsMappedSpline<d,T>::BezierOperator() const
{
    GISMO_ENSURE( 2==domainDim(), "Anything other than bivariate splines is not yet supported!");

    // Loop over all the elements of the given multipatch and collect all relevant
    // information in ElementBlocks. These will be grouped in a std::map
    // with respect to the number of active basis functions ( = NN/NCV )
    // of each Bezier element
    std::map<index_t, internal::ElementBlock> ElementBlocks;

    index_t NN; // Number of control points of the Bezier element
    gsMatrix<index_t> actives; // Active basis functions

    gsMatrix<T> tmp, center_point(2,1);
    center_point.setConstant( (T)(0.5) );

    for (size_t p=0; p<nPatches(); ++p)
    {
        const gsBasis<T> * basis = & this->getBase(p);
        actives = basis->active( center_point );
        NN = actives.size();
        ElementBlocks[NN].numElements += 1;    // Increment the Number of Elements contained in the ElementBlock
        ElementBlocks[NN].actives.push_back(actives);  // Append the active basis functions ( = the Node IDs ) for this element.
        ElementBlocks[NN].PR = basis->degree(0);
        ElementBlocks[NN].PS = basis->degree(1);
        ElementBlocks[NN].PT = 0; // TODO: if implemented for trivariates fix this
        
        this->exportPatch(p, tmp);
        ElementBlocks[NN].coefVectors.push_back( tmp ); //(!) If it is not  bezier
    }

    return ElementBlocks;
}


} //namespace gismo
