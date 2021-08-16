/** @file gsCompositeIncrSmoothnessGeom.hpp

    @brief implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsMSplines/gsCompositeIncrSmoothnessGeom.h>
#include <gsMSplines/gsCompositeUtils.h>

#define TO_INCRSMOOTHNESS(x) static_cast<gsCompositeIncrSmoothnessBasis<d,T> *>(x)

namespace gismo
{

template<short_t d,class T>
gsCompositeIncrSmoothnessGeom<d,T>::gsCompositeIncrSmoothnessGeom( const gsCompositeIncrSmoothnessBasis<d,T> & basis, const gsMatrix<T> & coefs ) : Base()
{
    m_compBasis=basis.clone().release();
    Base::m_basis=new gsMappedSingleBasis<d,T>(m_compBasis->getMappedSingleBasis(0));
    Base::m_coefs=coefs;
}

template<short_t d,class T>
gsCompositeIncrSmoothnessGeom<d,T>::gsCompositeIncrSmoothnessGeom( gsMultiPatch<T> const & mp,int incrSmoothness,int minEVDistance ) : Base()
{
    short_t geoDim = mp.geoDim();
    std::vector<gsMatrix<T> * > coefs;
    for(size_t i = 0;i<mp.nPatches();++i)
        coefs.push_back( new gsMatrix<T>(mp.patch(i).coefs()) );
    m_compBasis=getCompBasisFromMultiPatch_withCoefs<d>(mp,coefs,incrSmoothness,minEVDistance);
    if(Base::m_compBasis==NULL)
        GISMO_ERROR("no known basis for gsMappedGeom");
    Base::m_basis=new gsMappedSingleBasis<d,T>(Base::m_compBasis->getMappedSingleBasis(0));
    int start = 0, end = -1;
    gsMatrix<T> localCoefs;
    localCoefs.resize(Base::m_compBasis->localSize(),geoDim);
    for(size_t i = 0;i<mp.nPatches();i++)
    {
        start=end+1;
        end+=Base::m_compBasis->localSize(i);
        localCoefs.block(start,0,end-start+1,geoDim) << *(coefs[i]);
    }
    Base::m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
    freeAll(coefs);
}

template<short_t d,class T>
gsCompositeIncrSmoothnessGeom<d,T>::gsCompositeIncrSmoothnessGeom(gsMultiPatch<T> const  & mp,
                                                                  std::vector<patchCorner> C0List,
                                                                  int                      /*incrSmoothness*/,
                                                                  int                      /*minEVDistance*/)
                                                                  : Base()
{
    short_t geoDim = mp.geoDim();
    std::vector<gsMatrix<T> * > coefs;
    for(size_t i = 0;i<mp.nPatches();++i)
        coefs.push_back(new gsMatrix<T>(mp.patch(i).coefs()) );
    m_compBasis=getCompBasisFromMultiPatch<d,T>(mp);
    if(Base::m_compBasis==NULL)
        GISMO_ERROR("no known basis for gsMappedGeom");
    for(unsigned i=0;i<C0List.size();i++)
    {
        TO_INCRSMOOTHNESS(Base::m_compBasis)->setC0(C0List[i]);
    }
    TO_INCRSMOOTHNESS(Base::m_compBasis)->updateTopol();
    Base::m_basis=new gsMappedSingleBasis<d,T>(m_compBasis->getMappedSingleBasis(0));
    int start = 0, end = -1;
    gsMatrix<T> localCoefs;
    localCoefs.resize(Base::m_compBasis->localSize(),geoDim);
    for(size_t i = 0;i<mp.nPatches();i++)
    {
        start=end+1;
        end+=Base::m_compBasis->localSize(i);
        localCoefs.block(start,0,end-start+1,geoDim) << *(coefs[i]);
    }
    Base::m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
}

template<short_t d,class T>
void gsCompositeIncrSmoothnessGeom<d,T>::setCornerC0(patchCorner const & pc)
{
    gsMatrix<T> localCoefs;
    m_compBasis->global_coef_to_local_coef(Base::m_coefs,localCoefs);
    TO_INCRSMOOTHNESS(m_compBasis)->setC0(pc);
    m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
}

template<short_t d,class T>
void gsCompositeIncrSmoothnessGeom<d,T>::smoothCornerEdge(const patchCorner& pc,const patchSide& ps,bool updateBasis)
{
    gsMatrix<T> localCoefs;
    m_compBasis->global_coef_to_local_coef(Base::m_coefs,localCoefs);
    TO_INCRSMOOTHNESS(m_compBasis)->smoothCornerEdge(pc,ps,updateBasis);
    m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
}

template<short_t d,class T>
void gsCompositeIncrSmoothnessGeom<d,T>::smoothEverything()
{
    gsMatrix<T> localCoefs;
    m_compBasis->global_coef_to_local_coef(Base::m_coefs,localCoefs);
    TO_INCRSMOOTHNESS(m_compBasis)->smoothEverything();
    m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
}

template<short_t d,class T>
void gsCompositeIncrSmoothnessGeom<d,T>::uniformRefine(int numKnots, int mul)
{
    gsMatrix<T> localCoefs;
    m_compBasis->global_coef_to_local_coef(Base::m_coefs,localCoefs);
    TO_INCRSMOOTHNESS(m_compBasis)->uniformRefine_withCoefs(localCoefs, numKnots,mul,true);
    m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
}

template<short_t d,class T>
void gsCompositeIncrSmoothnessGeom<d,T>::refine(const index_t patch, const gsMatrix<T> &boxes)
{
    gsMatrix<T> localCoefs;
    m_compBasis->global_coef_to_local_coef(Base::m_coefs,localCoefs);
    //std::cout << Base::m_coefs << std::endl <<std::endl;
    TO_INCRSMOOTHNESS(m_compBasis)->refine_withCoefs(localCoefs, patch, boxes,true);
    //std::cout << localCoefs << std::endl <<std::endl;
    m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
    //std::cout << Base::m_coefs << std::endl <<std::endl;
}

template<short_t d,class T>
void gsCompositeIncrSmoothnessGeom<d,T>::refineElements(const index_t patch, std::vector<index_t> const & boxes)
{
    gsMatrix<T> localCoefs;
    m_compBasis->global_coef_to_local_coef(Base::m_coefs,localCoefs);
    TO_INCRSMOOTHNESS(m_compBasis)->refineElements_withCoefs(localCoefs, patch, boxes,true);
    m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
}

template<short_t d,class T>
void gsCompositeIncrSmoothnessGeom<d,T>::refineElements(std::vector<std::vector<index_t> > const & boxes)
{
    GISMO_ASSERT(boxes.size()==m_compBasis->nPatches(),"number of refinementvectors and number of patches must agree");
    gsMatrix<T> localCoefs;
    m_compBasis->global_coef_to_local_coef(Base::m_coefs,localCoefs);
    for(unsigned i = 0;i<boxes.size();++i)
    {
        TO_INCRSMOOTHNESS(m_compBasis)->refineElements_withCoefs(localCoefs, i, boxes[i],false);
    }
    TO_INCRSMOOTHNESS(m_compBasis)->repairPatches(localCoefs);
    TO_INCRSMOOTHNESS(m_compBasis)->updateTopol();
    m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
}

template<short_t d,class T>
void gsCompositeIncrSmoothnessGeom<d,T>::uniformRefineAndSmooth(int numKnots)
{
    gsMatrix<T> localCoefs;
    m_compBasis->global_coef_to_local_coef(Base::m_coefs,localCoefs);
    TO_INCRSMOOTHNESS(m_compBasis)->uniformRefine_withCoefs(localCoefs, numKnots,false);
    TO_INCRSMOOTHNESS(m_compBasis)->smoothEverything();
    m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
}

template<short_t d,class T>
void gsCompositeIncrSmoothnessGeom<d,T>::refineAndSmooth(const index_t patch, const gsMatrix<T> &boxes)
{
    gsMatrix<T> localCoefs;
    m_compBasis->global_coef_to_local_coef(Base::m_coefs,localCoefs);
    TO_INCRSMOOTHNESS(m_compBasis)->refine_withCoefs(localCoefs, patch,boxes,false);
    TO_INCRSMOOTHNESS(m_compBasis)->smoothEverything();
    m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
}

template<short_t d,class T>
void gsCompositeIncrSmoothnessGeom<d,T>::refineElementsAndSmooth(const index_t patch, std::vector<index_t> const & boxes)
{
    gsMatrix<T> localCoefs;
    m_compBasis->global_coef_to_local_coef(Base::m_coefs,localCoefs);
    TO_INCRSMOOTHNESS(m_compBasis)->refineElements_withCoefs(localCoefs, patch,boxes,false);
    TO_INCRSMOOTHNESS(m_compBasis)->smoothEverything();
    m_compBasis->local_coef_to_global_coef(localCoefs,Base::m_coefs);
}

} // namespace gismo

#undef TO_INCRSMOOTHNESS
