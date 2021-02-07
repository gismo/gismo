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
gsMappedSpline<d,T>::gsMappedSpline( gsMultiPatch<T> const & mp,std::string pathToMap ) : Base()
{
    short_t geoDim = mp.geoDim();
    std::vector<gsMatrix<T> * > coefs;
    for(int i = 0;i<mp.size();++i)
        coefs.push_back( new gsMatrix<T>(mp.patch(i).coefs() ) );
    m_compBasis = new gsMappedBasis<d,T>(mp,pathToMap);
    Base::m_basis=new gsMappedSingleBasis<d,T>(m_compBasis->getMappedSingleBasis(0));
    int start = 0, end = -1;
    gsMatrix<T> localCoefs;
    localCoefs.resize(m_compBasis->localSize(),geoDim);
    for(size_t i = 0;i<mp.nPatches();i++)
    {
        start=end+1;
        end+=m_compBasis->localSize(i);
        localCoefs.block(start,0,end-start+1,geoDim) << *(coefs[i]);
    }
    m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
    freeAll(coefs);
}

template<short_t d,class T>
gsMappedSpline<d,T>::gsMappedSpline(const gsMultiPatch<T> & mp, const gsSparseMatrix<T> & bmap ) : Base()
{
    m_compBasis = new gsMappedBasis<d,T>(gsMultiBasis<T>(mp),bmap);
    Base::m_basis=new gsMappedSingleBasis<d,T>(m_compBasis->getMappedSingleBasis(0));
    int start = 0, end = -1;
    short_t geoDim = mp.geoDim();
    gsMatrix<T> localCoefs(m_compBasis->localSize(),geoDim);
    for(size_t i = 0;i<mp.nPatches();i++)
    {
        start=end+1;
        end+=m_compBasis->localSize(i);
        localCoefs.block(start,0,end-start+1,geoDim) << mp.patch(i).coefs();
    }
    m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
}

template<short_t d,class T>
gsMappedSpline<d,T>::gsMappedSpline( const gsMappedBasis<d,T> & basis, const gsMatrix<T> & coefs ) : Base()
{
    m_compBasis=basis.clone().release();
    Base::m_basis=new gsMappedSingleBasis<d,T>(m_compBasis->getMappedSingleBasis(0));
    Base::m_coefs=coefs;
}

template<short_t d,class T>
gsMappedSpline<d,T>::gsMappedSpline( const gsMappedSpline& other ) : Base( other )
{
    m_compBasis=other.m_compBasis->clone().release();
    delete Base::m_basis;
    Base::m_basis = new gsMappedSingleBasis<d,T>(m_compBasis->getMappedSingleBasis(0));
}

template<short_t d,class T>
gsMappedSpline<d,T> & gsMappedSpline<d,T>::operator=( const gsMappedSpline& other )
{
    Base::operator= (other);
    delete m_compBasis;
    m_compBasis=other.m_compBasis->clone().release();
    delete Base::m_basis;
    Base::m_basis=new gsMappedSingleBasis<d,T>(m_compBasis->getMappedSingleBasis(0));
    return *this;
}

template<short_t d,class T>
gsMultiPatch<T> gsMappedSpline<d,T>::exportToPatches() const
{
    gsMatrix<T> localCoef;
    m_compBasis->global_coef_to_local_coef(Base::m_coefs,localCoef);
    return m_compBasis->exportToPatches(localCoef);
}

}
