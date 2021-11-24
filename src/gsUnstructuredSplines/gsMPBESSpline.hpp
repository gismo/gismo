/** @file gsMPBESSpline.hpp

    @brief implementation file

    paper: https://hal.inria.fr/hal-02272215/

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsMSplines/gsMPBESSpline.h>
#include <gsMSplines/gsMPBESUtils.h>

#define TO_INCRSMOOTHNESS(x) static_cast<gsMPBESBasis<d,T> *>(x)

namespace gismo
{

template<short_t d,class T>
gsMPBESSpline<d,T>::gsMPBESSpline( const gsMPBESBasis<d,T> & basis, const gsMatrix<T> & coefs )
:
Base(basis,coefs)
{
}

template<short_t d,class T>
gsMPBESSpline<d,T>::gsMPBESSpline( gsMultiPatch<T> const & mp,int incrSmoothness,int minEVDistance ) : Base()
{
    short_t geoDim = mp.geoDim();
    std::vector<gsMatrix<T> * > coefs;
    for(size_t i = 0;i<mp.nPatches();++i)
        coefs.push_back( new gsMatrix<T>(mp.patch(i).coefs()) );
    m_mbases=getCompBasisFromMultiPatch_withCoefs<d>(mp,coefs,incrSmoothness,minEVDistance);
    if(m_mbases==NULL)
        GISMO_ERROR("no known basis for gsMappedGeom");
    int start = 0, end = -1;
    gsMatrix<T> localCoefs;
    localCoefs.resize(m_mbases->localSize(),geoDim);
    for(size_t i = 0;i<mp.nPatches();i++)
    {
        start=end+1;
        end+=m_mbases->localSize(i);
        localCoefs.block(start,0,end-start+1,geoDim) << *(coefs[i]);
    }
    m_mbases->local_coef_to_global_coef(localCoefs,m_global);
    freeAll(coefs);
}

template<short_t d,class T>
gsMPBESSpline<d,T>::gsMPBESSpline(gsMultiPatch<T> const  & mp,
                                                                  std::vector<patchCorner> C0List,
                                                                  int                      /*incrSmoothness*/,
                                                                  int                      /*minEVDistance*/)
                                                                  : Base()
{
    short_t geoDim = mp.geoDim();
    std::vector<gsMatrix<T> * > coefs;
    for(size_t i = 0;i<mp.nPatches();++i)
        coefs.push_back(new gsMatrix<T>(mp.patch(i).coefs()) );
    m_mbases=getCompBasisFromMultiPatch<d,T>(mp);
    if(m_mbases==NULL)
        GISMO_ERROR("no known basis for gsMappedGeom");
    for(unsigned i=0;i<C0List.size();i++)
    {
        TO_INCRSMOOTHNESS(m_mbases)->setC0(C0List[i]);
    }
    TO_INCRSMOOTHNESS(m_mbases)->updateTopol();
    int start = 0, end = -1;
    gsMatrix<T> localCoefs;
    localCoefs.resize(m_mbases->localSize(),geoDim);
    for(size_t i = 0;i<mp.nPatches();i++)
    {
        start=end+1;
        end+=m_mbases->localSize(i);
        localCoefs.block(start,0,end-start+1,geoDim) << *(coefs[i]);
    }
    m_mbases->local_coef_to_global_coef(localCoefs,m_global);
}

template<short_t d,class T>
void gsMPBESSpline<d,T>::setCornerC0(patchCorner const & pc)
{
    gsMatrix<T> localCoefs;
    m_mbases->global_coef_to_local_coef(m_global,localCoefs);
    TO_INCRSMOOTHNESS(m_mbases)->setC0(pc);
    m_mbases->local_coef_to_global_coef(localCoefs,m_global);
}

template<short_t d,class T>
void gsMPBESSpline<d,T>::smoothCornerEdge(const patchCorner& pc,const patchSide& ps,bool updateBasis)
{
    gsMatrix<T> localCoefs;
    m_mbases->global_coef_to_local_coef(m_global,localCoefs);
    TO_INCRSMOOTHNESS(m_mbases)->smoothCornerEdge(pc,ps,updateBasis);
    m_mbases->local_coef_to_global_coef(localCoefs,m_global);
}

template<short_t d,class T>
void gsMPBESSpline<d,T>::smoothEverything()
{
    gsMatrix<T> localCoefs;
    m_mbases->global_coef_to_local_coef(m_global,localCoefs);
    TO_INCRSMOOTHNESS(m_mbases)->smoothEverything();
    m_mbases->local_coef_to_global_coef(localCoefs,m_global);
}

template<short_t d,class T>
void gsMPBESSpline<d,T>::uniformRefine(int numKnots, int mul)
{
    gsMatrix<T> localCoefs;
    m_mbases->global_coef_to_local_coef(m_global,localCoefs);
    TO_INCRSMOOTHNESS(m_mbases)->uniformRefine_withCoefs(localCoefs, numKnots,mul,true);
    m_mbases->local_coef_to_global_coef(localCoefs,m_global);
}

template<short_t d,class T>
void gsMPBESSpline<d,T>::refine(const index_t patch, const gsMatrix<T> &boxes)
{
    gsMatrix<T> localCoefs;
    m_mbases->global_coef_to_local_coef(m_global,localCoefs);
    //std::cout << m_global << std::endl <<std::endl;
    TO_INCRSMOOTHNESS(m_mbases)->refine_withCoefs(localCoefs, patch, boxes,true);
    //std::cout << localCoefs << std::endl <<std::endl;
    m_mbases->local_coef_to_global_coef(localCoefs,m_global);
    //std::cout << m_global << std::endl <<std::endl;
}

template<short_t d,class T>
void gsMPBESSpline<d,T>::refineElements(const index_t patch, std::vector<index_t> const & boxes)
{
    gsMatrix<T> localCoefs;
    m_mbases->global_coef_to_local_coef(m_global,localCoefs);
    TO_INCRSMOOTHNESS(m_mbases)->refineElements_withCoefs(localCoefs, patch, boxes,true);
    m_mbases->local_coef_to_global_coef(localCoefs,m_global);
}

template<short_t d,class T>
void gsMPBESSpline<d,T>::refineElements(std::vector<std::vector<index_t> > const & boxes)
{
    GISMO_ASSERT(boxes.size()==m_mbases->nPatches(),"number of refinementvectors and number of patches must agree");
    gsMatrix<T> localCoefs;
    m_mbases->global_coef_to_local_coef(m_global,localCoefs);
    for(unsigned i = 0;i<boxes.size();++i)
    {
        TO_INCRSMOOTHNESS(m_mbases)->refineElements_withCoefs(localCoefs, i, boxes[i],false);
    }
    TO_INCRSMOOTHNESS(m_mbases)->repairPatches(localCoefs);
    TO_INCRSMOOTHNESS(m_mbases)->updateTopol();
    m_mbases->local_coef_to_global_coef(localCoefs,m_global);
}

template<short_t d,class T>
void gsMPBESSpline<d,T>::uniformRefineAndSmooth(int numKnots)
{
    gsMatrix<T> localCoefs;
    m_mbases->global_coef_to_local_coef(m_global,localCoefs);
    TO_INCRSMOOTHNESS(m_mbases)->uniformRefine_withCoefs(localCoefs, numKnots,false);
    TO_INCRSMOOTHNESS(m_mbases)->smoothEverything();
    m_mbases->local_coef_to_global_coef(localCoefs,m_global);
}

template<short_t d,class T>
void gsMPBESSpline<d,T>::refineAndSmooth(const index_t patch, const gsMatrix<T> &boxes)
{
    gsMatrix<T> localCoefs;
    m_mbases->global_coef_to_local_coef(m_global,localCoefs);
    TO_INCRSMOOTHNESS(m_mbases)->refine_withCoefs(localCoefs, patch,boxes,false);
    TO_INCRSMOOTHNESS(m_mbases)->smoothEverything();
    m_mbases->local_coef_to_global_coef(localCoefs,m_global);
}

template<short_t d,class T>
void gsMPBESSpline<d,T>::refineElementsAndSmooth(const index_t patch, std::vector<index_t> const & boxes)
{
    gsMatrix<T> localCoefs;
    m_mbases->global_coef_to_local_coef(m_global,localCoefs);
    TO_INCRSMOOTHNESS(m_mbases)->refineElements_withCoefs(localCoefs, patch,boxes,false);
    TO_INCRSMOOTHNESS(m_mbases)->smoothEverything();
    m_mbases->local_coef_to_global_coef(localCoefs,m_global);
}

} // namespace gismo

#undef TO_INCRSMOOTHNESS
