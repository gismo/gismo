/** @file gsMPBESHSplineBasis.hpp

    @brief implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsUnstructuredSplines/gsMPBESHSplineBasis.h>
#include <gsCore/gsMultiBasis.h>

namespace gismo
{

template<short_t d, class T>
gsMPBESHSplineBasis<d,T>::gsMPBESHSplineBasis (std::vector<BasisType *> const & bases, gsBoxTopology const & topol,
                   int increaseSmoothnessLevel, int minEVDistance)
{
    typedef typename std::vector<BasisType *>::const_iterator tempBasisIter;
    m_topol = topol;
    for(tempBasisIter it = bases.begin();it!=bases.end();++it)
        m_bases.push_back((BasisType*)(*it)->clone().release());
    if(increaseSmoothnessLevel==-1)
        m_incrSmoothnessDegree=this->maxDegree()-1;
    else
        m_incrSmoothnessDegree=increaseSmoothnessLevel;
    m_minDist=static_cast<unsigned>(std::max<int>(std::min<int>(m_incrSmoothnessDegree+1,maxDegree()),minEVDistance));
    _initVertices();
    _setDistanceOfAllVertices();

    repairPatches();

    if(!_checkTopologyWithBases())
        GISMO_ERROR("topology and bases not suitable for composite basis");
    _setMapping();
}

template<short_t d, class T>
gsMPBESHSplineBasis<d,T>::gsMPBESHSplineBasis (std::vector<BasisType *> const & bases, gsBoxTopology const & topol,
                   std::vector<gsMatrix<T> * > & coefs,int increaseSmoothnessLevel, int minEVDistance)
{
    m_topol = topol;
    for(typename BasisContainer::const_iterator it = bases.begin();it!=bases.end();++it)
        m_bases.push_back( (gsBasis<T>*)(*it)->clone().release());
    if(increaseSmoothnessLevel==-1)
        m_incrSmoothnessDegree=this->maxDegree()-1;
    else
        m_incrSmoothnessDegree=increaseSmoothnessLevel;
    m_minDist=static_cast<unsigned>(std::max<int>(std::min<int>(m_incrSmoothnessDegree+1,maxDegree()),minEVDistance));
    _initVertices();
    _setDistanceOfAllVertices();

    repairPatches(coefs);

    if(!_checkTopologyWithBases())
        GISMO_ERROR("topology and bases not suitable for composite basis");
    _setMapping();
}

template<short_t d, class T>
gsMPBESHSplineBasis<d,T>::gsMPBESHSplineBasis (BasisType const & base, gsBoxTopology const & topol)
{
    m_topol = topol;
    m_bases.push_back((BasisType *)base.clone().release());
    m_incrSmoothnessDegree=this->maxDegree()-1;
    m_minDist=static_cast<unsigned>(std::min<int>(m_incrSmoothnessDegree+1,maxDegree()));
    _initVertices();
    _setDistanceOfAllVertices();

    repairPatches();

    _checkTopologyWithBases();
    _setMapping();
}

template<short_t d, class T>
gsMPBESHSplineBasis<d,T>::gsMPBESHSplineBasis( gsMultiPatch<T> const & mp, int increaseSmoothnessLevel,
                   int minEVDistance)
{
    //topol.computeAllVertices();
    for (size_t i = 0; i < mp.nPatches(); i++)
    {
        GISMO_ASSERT(dynamic_cast<gsHTensorBasis<2> * >(& mp.basis(i))!=NULL,"Bases is not of type gsHTensorBasis<2>");
    }
    m_topol = mp;
    m_bases = mp.basesCopy();
    if(increaseSmoothnessLevel==-1)
        m_incrSmoothnessDegree=this->maxDegree()-1;
    else
        m_incrSmoothnessDegree=increaseSmoothnessLevel;
    m_minDist=static_cast<unsigned>(std::max<int>(std::min<int>(m_incrSmoothnessDegree+1,maxDegree()),minEVDistance));
    _initVertices();
    _setDistanceOfAllVertices(); // this is wrong, every vertex has distance according to degree_u or degree_V

    repairPatches();

    _checkTopologyWithBases();
    _setMapping();
}

template<short_t d, class T>
gsMPBESHSplineBasis<d,T>::gsMPBESHSplineBasis( gsMultiBasis<T> const & mb,gsBoxTopology const & topol, int increaseSmoothnessLevel,
                   int minEVDistance)
{
    //topol.computeAllVertices();
    for (size_t i = 0; i < mb.nBases(); i++)
    {
        GISMO_ASSERT(dynamic_cast<const gsHTensorBasis<2> * >(&mb.basis(i))!=NULL,"Bases is not of type gsHTensorBasis<2>");
        m_bases.push_back((BasisType *)mb.basis(i).clone().release());
    }
    m_topol = topol;
    if(increaseSmoothnessLevel==-1)
        m_incrSmoothnessDegree=this->maxDegree()-1;
    else
        m_incrSmoothnessDegree=increaseSmoothnessLevel;
    m_minDist=static_cast<unsigned>(std::max<int>(std::min<int>(m_incrSmoothnessDegree+1,maxDegree()),minEVDistance));
    _initVertices();
    _setDistanceOfAllVertices(); // this is wrong, every vertex has distance according to degree_u or degree_V

    repairPatches();

    _checkTopologyWithBases();
    _setMapping();
}

template<short_t d, class T>
gsMPBESHSplineBasis<d,T>::gsMPBESHSplineBasis( const gsMPBESHSplineBasis& other )
{
    m_topol = other.m_topol;
    m_vertices=other.m_vertices;
    m_distances=other.m_distances;
    m_minDist=other.m_minDist;
    // clone all geometries
    for(ConstBasisIter it = other.m_bases.begin();it!=other.m_bases.end();++it)
    {
        m_bases.push_back( (BasisType *)(*it)->clone().release() );
    }
    m_incrSmoothnessDegree=other.m_incrSmoothnessDegree;
    this->_checkTopologyWithBases();
    //_setMapping();
    m_mapper=new gsWeightMapper<T>(*other.m_mapper);
    //m_mapFactory= need a method to clone mapFactory
}

template<short_t d, class T>
gsMPBESHSplineBasis<d,T>& gsMPBESHSplineBasis<d,T>::operator=( const gsMPBESHSplineBasis& other )
{
    m_topol=other.m_topol;
    m_vertices=other.m_vertices;
    m_distances=other.m_distances;
    m_minDist=other.m_minDist;
    freeAll(m_bases);
    m_bases.clear();
    for(ConstBasisIter it = other.m_bases.begin();it!=other.m_bases.end();++it)
        m_bases.push_back((BasisType*)(*it)->clone().release());
    if(m_mapper)
        delete m_mapper;
    m_mapper=NULL;//other.m_mapper->clone();
    m_incrSmoothnessDegree=other.m_incrSmoothnessDegree;
    //m_mapFactory= need a method to clone mapFactory
    return *this;
}

template<short_t d, class T>
void gsMPBESHSplineBasis<d,T>::_setMapping()
{
    // * Initializer mapper
    gsMPBESMapHB2D<d,T> maker(m_incrSmoothnessDegree,& m_topol,this);
    if(m_mapper)
        delete m_mapper;
    m_mapper = maker.makeMapper();
}

template<short_t d, class T>
unsigned gsMPBESHSplineBasis<d,T>::basisFunctionsOnSide(const patchSide& ps) const
{
    unsigned nr=0;
    std::vector<bool> actives;
    for(unsigned level = 0;level<=basis(ps.patch).maxLevel();++level)
    {
        basis(ps.patch).activeBoundaryFunctionsOfLevel(level,ps.side(),actives);
        for(std::vector<bool>::const_iterator iter = actives.begin();iter!=actives.end();++iter)
            if(*iter)
                nr++;
    }
    return nr;
}

template<short_t d, class T>
bool gsMPBESHSplineBasis<d,T>::isLocallyConnected(indexType i,indexType j) const
{
    unsigned patch_i = _getPatch(i), patch_j = _getPatch(j);
    if( patch_i != patch_j )
        return false;
    indexType loc_i = _getPatchIndex(i), loc_j = _getPatchIndex(j);
    unsigned level_i = basis(patch_i).levelOf(loc_i);
    unsigned level_j = basis(patch_j).levelOf(loc_j);
    if( level_i != level_j )
        return false;
    unsigned tensorIndex_i = basis(patch_i).flatTensorIndexOf(loc_i);
    unsigned tensorIndex_j = basis(patch_j).flatTensorIndexOf(loc_j);
    gsVector<index_t,d> vec_i = basis(patch_i).getBases()[level_i]->tensorIndex(tensorIndex_i);
    gsVector<index_t,d> vec_j = basis(patch_j).getBases()[level_j]->tensorIndex(tensorIndex_j);
    unsigned distance = 0;
    for( unsigned i2=0; i2< d;++i2 )
        if( vec_i[j]-vec_j[i2] != 0 )
            distance++;
    return distance==1;
}

template<short_t d, class T>
void gsMPBESHSplineBasis<d,T>::refine(const index_t patch, gsMatrix<T> const & boxes, bool updateBasis)
{
    basis(patch).refine(boxes);
    std::vector<gsMatrix<T> *> coefs;
    for (size_t i = 0; i < nPatches(); ++i)
        coefs.push_back(NULL);
    if(updateBasis)
    {
        repairPatches(coefs,patch);
        updateTopol();
    }
}

template<short_t d, class T>
void gsMPBESHSplineBasis<d,T>::refineElements(const index_t patch, std::vector<index_t> const & boxes, bool updateBasis)
{
    basis(patch).refineElements(boxes);
    std::vector<gsMatrix<T> *> coefs;
    for (size_t i = 0; i < nPatches(); ++i)
        coefs.push_back(NULL);
    if(updateBasis)
    {
        repairPatches(coefs,patch);
        updateTopol();
    }
}

template<short_t d, class T>
void gsMPBESHSplineBasis<d,T>::refine_withCoefs(gsMatrix<T>& localCoef, const index_t patch, gsMatrix<T> const & boxes, bool updateBasis)
{
    //std::cout << localCoef << std::endl << std::endl;
    std::vector<gsMatrix<T> *> coefs;
    unsigned geoDim = localCoef.cols();
    int start, end = -1;
    for (size_t i = 0; i < nPatches(); ++i)
    {
        start=end+1;
        end+=basis(i).size();
        gsMatrix<T>* localMat = new gsMatrix<T>(end-start+1,geoDim);
        *localMat << localCoef.block(start,0,end-start+1,geoDim);
        coefs.push_back(localMat);
    }
    basis(patch).refine_withCoefs(*coefs[patch],boxes);
    if(updateBasis)
        repairPatches(coefs,patch);
    unsigned totalSize=0;
    for (size_t i = 0; i < nPatches(); ++i)
    {
        totalSize+=basis(i).size();
    }
    localCoef.resize(totalSize,geoDim);
    end = -1;
    for (size_t i = 0; i < nPatches(); i++)
    {
        start=end+1;
        end+=basis(i).size();
        localCoef.block(start,0,end-start+1,geoDim) << *coefs[i];
    }
    //std::cout << localCoef << std::endl << std::endl;
    if(updateBasis)
        updateTopol();
    freeAll(coefs);
}

template<short_t d, class T>
void gsMPBESHSplineBasis<d,T>::refineElements_withCoefs(gsMatrix<T>& localCoef, const index_t patch, std::vector<index_t> const & boxes, bool updateBasis)
{
    std::vector<gsMatrix<T> *> coefs;
    unsigned geoDim = localCoef.cols();
    int start, end = -1;
    for (size_t i = 0; i < nPatches(); ++i)
    {
        start=end+1;
        end+=basis(i).size();
        gsMatrix<T>* localMat = new gsMatrix<T>(end-start+1,geoDim);
        *localMat << localCoef.block(start,0,end-start+1,geoDim);
        coefs.push_back(localMat);
    }
    basis(patch).refineElements_withCoefs(*coefs[patch],boxes);
    if(updateBasis)
        repairPatches(coefs,patch);
    unsigned totalSize=0;
    for (size_t i = 0; i < nPatches(); ++i)
    {
        totalSize+=basis(i).size();
    }
    localCoef.resize(totalSize,geoDim);
    end = -1;
    for (size_t i = 0; i < nPatches(); i++)
    {
        start=end+1;
        end+=basis(i).size();
        localCoef.block(start,0,end-start+1,geoDim) << *coefs[i];
    }
    if(updateBasis)
        updateTopol();
    freeAll(coefs);
}

template<short_t d, class T>
void gsMPBESHSplineBasis<d,T>::refineWithExtension(const index_t patch,gsMatrix<T> const & boxes, int refExt,bool updateBasis)
{
    m_bases[patch]->refine( boxes, refExt);
    std::vector<gsMatrix<T> *> coefs;
    for (size_t i = 0; i < nPatches(); ++i)
        coefs.push_back(NULL);
    if(updateBasis)
    {
        repairPatches(coefs,patch);
        updateTopol();
    }
}

template<short_t d, class T>
void gsMPBESHSplineBasis<d,T>::repairPatches(std::vector<gsMatrix<T> *> & coefs, index_t startFromPatch)
{
    std::vector<size_t> toCheck; //set of all indizes of patches, which have to be checked
    if(startFromPatch==-1)
        for(size_t i = 0; i < m_bases.size(); ++i)
            toCheck.push_back(i);
    else
        toCheck.push_back(startFromPatch);
    size_t curElement;
    do
    {
        curElement=toCheck[0];
        toCheck.erase(toCheck.begin());
        std::vector<index_t> boxes;
        std::vector<index_t> checkPatches(0);
        bool matched = true;
        matched = matched && _innerBoxesAreSuitable(curElement,boxes);
        if(!matched)
        {
            if(coefs[curElement]!=NULL)
                basis(curElement).refineElements_withCoefs(*(coefs[curElement]),boxes);
            else
                basis(curElement).refineElements(boxes);
            toCheck.push_back(curElement);
            continue;
        }
        matched = true;
        boxes.clear();
        matched = matched && _boxesMatchNeighbours(curElement,boxes,checkPatches);
        if(!matched)
        {
            if(coefs[curElement]!=NULL)
                basis(curElement).refineElements_withCoefs(*(coefs[curElement]),boxes);
            else
                basis(curElement).refineElements(boxes);
        }
        for(size_t i = 0; i < checkPatches.size(); ++i)
            if(std::find(toCheck.begin(), toCheck.end(), checkPatches[i])==toCheck.end())
                toCheck.push_back(checkPatches[i]);
    }while(!toCheck.empty());
}

template<short_t d, class T>
bool gsMPBESHSplineBasis<d,T>::_innerBoxesAreSuitable(const index_t patch,
                            std::vector<index_t>& boxes)
{
    size_t sz = boxes.size();
    short_t dist_u = 2*std::min(m_bases[patch]->degree(0),(short_t)(m_incrSmoothnessDegree+1));
    short_t dist_v = 2*std::min(m_bases[patch]->degree(1),(short_t)(m_incrSmoothnessDegree+1));
    patchSide ps_north(patch,boundary::north);
    patchSide ps_south(patch,boundary::south);
    patchSide ps_east(patch,boundary::east);
    patchSide ps_west(patch,boundary::west);
    bool north = m_topol.isInterface(ps_north);
    bool south = m_topol.isInterface(ps_south);
    bool east  = m_topol.isInterface(ps_east);
    bool west  = m_topol.isInterface(ps_west);
    gsMatrix<index_t> b1,b2;
    gsVector<index_t> level;
    // make the interior of the patch ok
    basis(patch).tree().getBoxesInLevelIndex(b1,b2,level); // sorted by level?!?!
    for(index_t i = 0; i < level.rows(); i++)
    {
        index_t l = level(i);
        if(l==0)
            continue;
        index_t b_uMin = b1(i,0),b_vMin = b1(i,1);
        index_t b_uMax = b2(i,0),b_vMax = b2(i,1);
        index_t uMin = 0, uMax = basis(patch).getBases()[l]->knots(0).uSize() - 1;
        index_t vMin = 0, vMax = basis(patch).getBases()[l]->knots(1).uSize() - 1;
        // if the box does not touch the boundary we have to do nothing
        if( uMin < b_uMin && b_uMax < uMax && vMin < b_vMin && b_vMax < vMax )
            continue;
        // if the box touches the boundary, we have to check if goes far enough from
        // the boundary, such that all the functions needed for gluing are there in the same level
        if( west && uMin == b_uMin && b_uMax < uMin+dist_u )
            _addBox(patch,b_uMin,b_vMin,uMin+dist_u,b_vMax,l,boxes);
        if( east && uMax == b_uMax && b_uMin > uMax-dist_u )
            _addBox(patch,uMax-dist_u,b_vMin,b_uMax,b_vMax,l,boxes);
        if( south && vMin == b_vMin && b_vMax < vMin+dist_v )
            _addBox(patch,b_uMin,b_vMin,b_uMax,vMin+dist_v,l,boxes);
        if( north && vMax == b_vMax && b_vMin > vMax-dist_v )
            _addBox(patch,b_uMin,vMax-dist_v,b_uMax,b_vMax,l,boxes);
    }
    return boxes.size()==sz;
}

template<short_t d, class T>
bool gsMPBESHSplineBasis<d,T>::_boxesMatchNeighbours(const index_t patch,
                           std::vector<index_t>& boxes, std::vector<index_t>& checkPatches)
{
    unsigned sz = boxes.size();
    //check if neighbour has the same refined boxes.
    patchSide ps, ps_neigh;
    std::vector<bool> sideToCheck(4,false);
    std::vector<int> neighbours(4,-1);
    gsVector<bool> orient;
    index_t patch_max_level = basis(patch).maxLevel();
    for(unsigned side=1;side<=4;++side)
    {
        ps=patchSide(patch,side);
        if( !m_topol.getNeighbour(ps,ps_neigh) )
            continue;
        neighbours[side-1]=ps_neigh.patch; //fill with the patch numbers of the neighbours
        unsigned neighbour_max_level = basis(ps_neigh.patch).maxLevel();
        unsigned max_level = std::max<unsigned>(patch_max_level,neighbour_max_level);
        for(unsigned l=1;l<=max_level;++l)
        {
            std::vector<bool> s_patch,s_neigh,s_res;
            basis(patch).activeBoundaryFunctionsOfLevel(l,ps.side(),s_patch);
            basis(ps_neigh.patch).activeBoundaryFunctionsOfLevel(l,ps_neigh.side(),s_neigh);
            unsigned sz2 = s_patch.size();
            boundaryInterface bf;
            m_topol.getInterface(ps,bf);
            if(!bf.dirOrientation(ps,1-ps.direction()))
                std::reverse(s_neigh.begin(),s_neigh.end());
            s_res.resize(sz2);
            for(unsigned i = 0;i<sz2;i++)
            {
                s_res[i]=s_neigh[i]&&!s_patch[i];
                if(s_patch[i]&&!s_neigh[i])
                    sideToCheck[side-1]=true; //set this side to true
            }
            unsigned start = 0, end=0;
            do
            {
                if(!s_res[start])
                {
                    start++;
                    end++;
                }
                else if(end<sz2&&s_res[end])
                    end++;
                else
                {
                    _addBoundaryBox(patch,ps.side(),start,end-1,l,boxes,sideToCheck);
                    start=end;
                }
            }while(start!=sz2);
        }
    }
    //find the patches we will have to check again:
    for(unsigned i = 0;i<4;i++)
    {
        if(sideToCheck[i]&&neighbours[i]!=-1)
            checkPatches.push_back(neighbours[i]);
    }
    return boxes.size()==sz;
}

template<short_t d, class T>
void gsMPBESHSplineBasis<d,T>::_addBoundaryBox(const index_t patch,const boxSide s,const int start, const int end,const unsigned level, std::vector<index_t> & boxes, std::vector<bool> & sideToCheck)
{
    short_t u_max = basis(patch).getBases()[level]->size(0)-1;
    short_t v_max = basis(patch).getBases()[level]->size(1)-1;
    short_t dist_u = 2*std::min(m_bases[patch]->degree(0),(short_t)(m_incrSmoothnessDegree+1));
    short_t dist_v = 2*std::min(m_bases[patch]->degree(1),(short_t)(m_incrSmoothnessDegree+1));
    switch(s)
    {
    case 1:
        _addFunBox(patch,0,start,dist_u-1,end,level,boxes);
        if(dist_u>=u_max)
            sideToCheck[2-1]=true;
        break;
    case 2:
        _addFunBox(patch,u_max-(dist_u-1),start,u_max,end,level,boxes);
        if(dist_u>=u_max)
            sideToCheck[1-1]=true;
        break;
    case 3:
        _addFunBox(patch,start,0,end,dist_v-1,level,boxes);
        if(dist_v>=v_max)
            sideToCheck[4-1]=true;
        break;
    case 4:
        _addFunBox(patch,start,v_max-(dist_v-1),end,v_max,level,boxes);
        if(dist_v>=v_max)
            sideToCheck[3-1]=true;
        break;
    default:
        GISMO_ERROR("only 2D possible");
    }
    if(s==1||s==2)
    {
        if(start<=0)
            sideToCheck[3-1]=true;
        if(end>=static_cast<int>(v_max)-1)
            sideToCheck[4-1]=true;
    }
    else
    {
        if(start<=0)
            sideToCheck[1-1]=true;
        if(end>=static_cast<int>(u_max)-1)
            sideToCheck[2-1]=true;
    }
}

template<short_t d, class T>
void gsMPBESHSplineBasis<d,T>::_addFunBox(const index_t patch,const unsigned uMin,const unsigned vMin,const unsigned uMax,const unsigned vMax,const unsigned level, std::vector<index_t> & boxes)
{
    boxes.push_back(level);
    gsMatrix<index_t> supportIndex;
    basis(patch).getBases()[level]->knots(0).supportIndex_into(uMin,supportIndex);
    boxes.push_back(supportIndex(0,0));
    basis(patch).getBases()[level]->knots(1).supportIndex_into(vMin,supportIndex);
    boxes.push_back(supportIndex(0,0));
    basis(patch).getBases()[level]->knots(0).supportIndex_into(uMax,supportIndex);
    boxes.push_back(supportIndex(0,1));
    basis(patch).getBases()[level]->knots(1).supportIndex_into(vMax,supportIndex);
    boxes.push_back(supportIndex(0,1));
}

template<short_t d, class T>
void gsMPBESHSplineBasis<d,T>::_addBox(const index_t patch,const unsigned uMin,const unsigned vMin,const unsigned uMax,const unsigned vMax,const unsigned level, std::vector<index_t> & boxes)
{
    gsVector<index_t,d> lowerLeft;
    lowerLeft(0)=uMin;
    lowerLeft(1)=vMin;
    gsVector<index_t,d> upperRight;
    upperRight(0)=uMax;
    upperRight(1)=vMax;
    if( uMin<uMax && vMin<vMax && static_cast<int>(level)>basis(patch).tree().query3(lowerLeft,upperRight,level))
    {
        boxes.push_back(level);
        boxes.push_back(uMin);
        boxes.push_back(vMin);
        boxes.push_back(uMax);
        boxes.push_back(vMax);
    }
}

template<short_t d, class T>
void gsMPBESHSplineBasis<d,T>::_endpointsOfActiveBoundaryFunctions(patchSide const & ps,bool orient,std::vector<T>& endpoints) const
{
    int patch = ps.patch;
    unsigned deg = degree(patch,1-(ps.direction()));
    std::vector<bool> actives;
    for(unsigned level = 0;level<=basis(patch).maxLevel();++level)
    {
        gsKnotVector<T> knots = basis(patch).getBases()[level]->knots(1-(ps.direction()));
        basis(patch).activeBoundaryFunctionsOfLevel(level,ps.side(),actives);
        if(orient)
        {
            knots.reverse();
            std::reverse(actives.begin(), actives.end());
        }
        for(unsigned i = 0; i<actives.size();++i)
            if(actives[i])
                endpoints.push_back(knots.at(i+deg+1));
    }
    std::sort(endpoints.begin(),endpoints.end());
}

template<short_t d, class T>
T gsMPBESHSplineBasis<d,T>::findParameter(patchSide const & ps,patchCorner const & pc,unsigned nrBasisFuncs) const
{
    if(nrBasisFuncs==0)
        return 0.0;
    std::vector<T> endpoints;
    gsVector<bool> pars;
    pc.parameters_into(d,pars);
    _endpointsOfActiveBoundaryFunctions(ps,pars(1-ps.direction()),endpoints);
    return endpoints[nrBasisFuncs-1];
}

}
