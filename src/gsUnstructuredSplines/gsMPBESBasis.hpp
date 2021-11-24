/** @file gsMPBESBasis.hpp

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsUnstructuredSplines/gsMPBESBasis.h>

namespace gismo
{

// template<short_t d, class T>
// void gsMPBESBasis<d,T>::_setMapping()
// {
//     // * Initializer mapper
//     gsMapFactory * maker = _getMapFactory();
//     if(m_mapper)
//         delete m_mapper;
//     m_mapper = maker->makeMapper();
//     delete maker;
// }

template<short_t d, class T>
bool gsMPBESBasis<d,T>::_check() const
{
    bool consistent = true;
    consistent = consistent && _checkTopologyWithBases();
    return consistent;
}

template<short_t d, class T>
void gsMPBESBasis<d,T>::connectivity(const gsMatrix<T> & nodes,gsMesh<T> & mesh) const
{
    const indexType sz  = size();
    GISMO_ASSERT( nodes.rows() == sz, "Invalid input.");

    // Add vertices
    for( indexType i = 0; i< sz; ++i )
        mesh.addVertex( nodes.row(i).transpose() );

    for( indexType i = 0; i< sz; ++i )
        for( indexType j = i+1; j< sz; ++j )
            if(isConnected(i,j))
                mesh.addEdge(i,j);
}

template<short_t d, class T>
bool gsMPBESBasis<d,T>::isConnected(indexType i,indexType j) const
{
    std::vector<indexType> locals_i,locals_j;
    m_mapper->targetToSource(i,locals_i);
    m_mapper->targetToSource(j,locals_j);
    indexType sz_i = locals_i.size(), sz_j = locals_j.size();

    if( sz_i==1 && sz_j==1 )
        return isLocallyConnected(locals_i[0],locals_j[0]);

    if( sz_i==1 || sz_j==1 )
    {
        std::vector<indexType> & lookthrough=locals_j;
        indexType element = sz_i==1 ? locals_i[0] : locals_j[0];
        if( sz_i==1 )
            lookthrough=locals_j;
        else
            lookthrough=locals_i;

        bool connected = false;
        for( unsigned i2 = 0; i2<lookthrough.size(); ++i2 )
            if( isLocallyConnected(element,lookthrough[i2]) )
            {
                connected = true;
                break;
            }
        return connected;
    }

    std::vector<indexType> setIntersection(sz_i+sz_j);
    std::vector<indexType>::iterator it;
    it=std::set_intersection (locals_i.begin(), locals_i.end(), locals_j.begin(), locals_j.end(), setIntersection.begin());
    indexType sz_intersect = it-setIntersection.begin();
    setIntersection.resize(sz_intersect);

    if( sz_intersect == sz_i-1 || sz_intersect == sz_j-1 )
        return true;

    if( sz_intersect == 0 )
    {
        bool connected,allConnected=true;
        for( indexType i3 = 0; i3< sz_i; ++i3 )
        {
            connected = false;
            for( indexType j3 = 0; j3< sz_j; ++j3 )
                if( isLocallyConnected(locals_i[i3],locals_j[j3]) )
                {
                    connected = true;
                    break;
                }
            if( !connected )
            {
                allConnected=false;
                break;
            }
        }
        if ( allConnected )
            return true;
        allConnected=true;
        for( indexType j4 = 0; j4< sz_j; ++j4 )
        {
            connected = false;
            for( indexType i5 = 0; i5< sz_i; ++i5 )
                if( isLocallyConnected(locals_i[i5],locals_j[j4]) )
                {
                    connected = true;
                    break;
                }
            if( !connected )
            {
                allConnected=false;
                break;
            }
        }
        if ( allConnected )
            return true;
    }

    return false;
}

template<short_t d, class T>
void gsMPBESBasis<d,T>::uniformRefine(int numKnots, int mul,bool updateBasis)
{
    int start, end = -1;
    for (BasisIter it=m_bases.begin();it!=m_bases.end();++it)
    {
        start=end+1;
        end=start+(*it)->size()-1;
        (*it)->uniformRefine(numKnots,mul);
    }
    if(updateBasis)
        updateTopol();
}

template<short_t d, class T>
void gsMPBESBasis<d,T>::uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots, int mul, bool updateBasis)
{
    int start, end = -1, geoDim=coefs.cols(), totalLength = 0;
    std::vector<gsMatrix<T> *> newCoefs;
    for (BasisIter it=m_bases.begin();it!=m_bases.end();++it)
    {
        start=end+1;
        end=start+(*it)->size()-1;
        gsMatrix<T> * it_coef = new gsMatrix<T>(end-start+1,geoDim);
        it_coef->setZero();
        *it_coef << coefs.block(start,0,end-start+1,geoDim);
        (*it)->uniformRefine_withCoefs(*it_coef, numKnots,mul);
        newCoefs.push_back(it_coef);
        totalLength+=it_coef->rows();
    }
    coefs.resize(totalLength,geoDim);
    end = -1;
    for(ConstMatrixPtrIter it=newCoefs.begin();it!=newCoefs.end();++it)
    {
        start=end+1;
        end=start+(*it)->rows()-1;
        coefs.block(start,0,end-start+1,geoDim)<<**it;
    }
    freeAll(newCoefs);
    if(updateBasis)
        updateTopol();
}

template<short_t d, class T>
void gsMPBESBasis<d,T>::degreeElevate( index_t amount , bool updateBasis )
{
    for (size_t i = 0; i < nPatches(); ++i)
        m_bases[i]->degreeElevate(amount);  // TODO: must be something else then gsBasis::degreeElevate
    if(updateBasis)
        updateTopol();
}

template<short_t d, class T>
void gsMPBESBasis<d,T>::degreeIncrease( index_t amount , index_t dir, bool updateBasis )
{
    for (size_t i = 0; i < nPatches(); ++i)
        m_bases[i]->degreeElevate(amount, dir);
    if(updateBasis)
        updateTopol();
}

template<short_t d, class T>
void gsMPBESBasis<d,T>::repairPatches(index_t startFromPatch)
{
    std::vector<gsMatrix<T> *> coefs;
    for (size_t i = 0; i < nPatches(); ++i)
        coefs.push_back(NULL);
    repairPatches(coefs,startFromPatch);
}

template<short_t d, class T>
void gsMPBESBasis<d,T>::smoothCornerEdge(const patchCorner&pc,const patchSide& ps,bool updateBasis)
{
    boundaryInterface interface;
    bool special = isSpecialVertex(pc);
    if( special && m_topol.getInterface(ps,interface) )
    {
        for(unsigned i=0;i<m_distances.size();++i)
            if(m_distances[i].isDistancesOfInterface(interface))
                m_distances[i].setParamDist(m_minDist,pc,*this);
    }
    if(updateBasis)
        updateTopol();
}

template<short_t d, class T>
void gsMPBESBasis<d,T>::smoothEverything()
{
    std::vector<std::vector<patchCorner> >cornerLists;
    std::vector<patchSide> psides;
    m_topol.getEVs(cornerLists);
    for(unsigned i=0;i<cornerLists.size();++i)
        for(unsigned j=0;j<cornerLists[i].size();++j)
        {
            patchCorner& pc=cornerLists[i][j];
            pc.getContainingSides(d,psides);
            for(unsigned k=0;k<psides.size();++k)
                smoothCornerEdge(pc,psides[k],false);
        }
    updateTopol();
}

template<short_t d, class T>
void gsMPBESBasis<d,T>::repairPatches(gsMatrix<T> & localCoef, index_t startFromPatch)
{
    std::vector<gsMatrix<T> *> patchCoefs;
    unsigned geoDim = localCoef.cols();
    int start, end = -1;
    for (size_t i = 0; i < nPatches(); ++i)
    {
        start=end+1;
        end+=m_bases[i]->size();
        gsMatrix<T>* localMat = new gsMatrix<T>(end-start+1,geoDim);
        *localMat << localCoef.block(start,0,end-start+1,geoDim);
        patchCoefs.push_back(localMat);
    }
    repairPatches(patchCoefs,startFromPatch);
    unsigned totalSize=0;
    for (size_t i = 0; i < nPatches(); ++i)
    {
        totalSize+=m_bases[i]->size();
    }
    localCoef.resize(totalSize,geoDim);
    end = -1;
    for (size_t i = 0; i < nPatches(); i++)
    {
        start=end+1;
        end+=m_bases[i]->size();
        localCoef.block(start,0,end-start+1,geoDim) << *patchCoefs[i];
    }
}

template<short_t d, class T>
T gsMPBESBasis<d,T>::getWeight(const patchSide & ps) const
{
    typedef typename std::vector<std::pair<patchSide,T> >::const_iterator cWeightPairIter;
    for(cWeightPairIter it=m_patchSideWeights.begin();it!=m_patchSideWeights.end();++it)
        if(it->first == ps)
            return it->second;
    return 1.0;
}

template<short_t d, class T>
bool gsMPBESBasis<d,T>::setWeight(const patchSide & ps, const T weight)
{
    typedef typename std::vector<std::pair<patchSide,T> >::iterator weightPairIter;
    bool found = false;
    for(weightPairIter it=m_patchSideWeights.begin();it!=m_patchSideWeights.end();++it)
        if(it->first == ps)
        {
            it->second = weight;
            found = true;
        }
    if(!found && m_topol.isInterface(ps))
    {
        m_patchSideWeights.push_back(std::make_pair(ps,weight));
        found = true;
    }
    return found;
}

template<short_t d, class T>
void gsMPBESBasis<d,T>::_initVertices()
{
    m_vertices.clear();
    std::vector<std::vector<patchCorner> > vertexLists;
    m_topol.getEVs(vertexLists);
    for(unsigned i = 0;i<vertexLists.size();++i)
        for(unsigned j=0;j<vertexLists[i].size();++j)
            m_vertices.push_back(std::make_pair(vertexLists[i][j],true));
    m_topol.getOVs(vertexLists);
    for(unsigned i = 0;i<vertexLists.size();++i)
        for(unsigned j=0;j<vertexLists[i].size();++j)
            m_vertices.push_back(std::make_pair(vertexLists[i][j],false));
}

template<short_t d, class T>
void gsMPBESBasis<d,T>::_setDistanceOfAllVertices()
{
    m_distances.clear();
    for(std::vector<boundaryInterface>::const_iterator iter = m_topol.iBegin();
        iter!=m_topol.iEnd();++iter)
    {
        const patchSide & ps = (*iter).first();
        gsVector<bool>pars(2);
        int dir = ps.direction();
        pars(dir  ) = ps.parameter();
        pars(1-dir)=0;
        patchCorner pc1 = patchCorner(ps.patch,pars);
        pars(1-dir)=1;
        patchCorner pc2 = patchCorner(ps.patch,pars);
        distances dist(*iter,pc1,pc2,*this);
        m_distances.push_back(dist);
    }
}

template<short_t d, class T>
void gsMPBESBasis<d,T>::setC0(patchCorner pc)
{
    std::vector<patchCorner> vertexList;
    m_topol.getCornerList(pc,vertexList);
    for(unsigned i = 0;i<vertexList.size();++i)
    {
        bool found = false;
        for(unsigned j = 0;j<m_vertices.size();++j)
            if(m_vertices[j].first==vertexList[i])
            {
                m_vertices[j].second=true;
                found=true;
            }
        if(!found)
            m_vertices.push_back(std::make_pair(vertexList[i],true));
    }
    std::vector<patchSide> sides;
    pc.getContainingSides(dim(),sides);
    for(unsigned i = 0;i<sides.size();++i)
    {
        boundaryInterface interface;
        if(m_topol.getInterface(sides[i],interface))
            for(unsigned j=0;j<m_distances.size();++j)
                if(m_distances[j].isDistancesOfInterface(interface))
                    m_distances[j].setParamDist(m_minDist,pc,*this);
    }
}

template<short_t d, class T>
bool gsMPBESBasis<d,T>::isSpecialVertex(const patchCorner & pc) const
{
    typedef std::vector<std::pair<patchCorner,bool> >::const_iterator cvertices_iter;
    for(cvertices_iter iter=m_vertices.begin(); iter!=m_vertices.end();++iter)
        if(pc==(*iter).first)
            return (*iter).second;
    return false;
}

template<short_t d, class T>
T gsMPBESBasis<d,T>::getParametricDistanceOfVertex(const patchCorner& pc,const patchSide& ps) const
{
    boundaryInterface interface;
    T param = 0.0;
    bool special = isSpecialVertex(pc);
    if( special && m_topol.getInterface(ps,interface) )
    {
        for(unsigned i=0;i<m_distances.size();++i)
            if(m_distances[i].isDistancesOfInterface(interface))
                param = m_distances[i].getParamDist(pc,*this);
    }
    else if( special )
       param = findParameter(ps,pc,m_minDist);
    else
    {
       const int deg=degree(pc.patch,1-ps.side().direction());
       const unsigned degree = std::min<int>(deg,m_incrSmoothnessDegree+1);
       for(unsigned i=0;i<m_vertices.size();++i)
           if(pc==m_vertices[i].first)
               param = findParameter(ps,pc,degree);
    }
    return param;
}

template<short_t d, class T>
gsMPBESBasis<d,T>::distances::distances(const boundaryInterface&  iface, const patchCorner& pc1,
          const patchCorner& pc2,const gsMPBESBasis<d,T>& basis) :
    interface(iface),corner1(pc1),corner2(pc2)
{
    const patchSide& ps = interface.first();
    const unsigned patch = ps.patch;
    const boxSide& side = ps.side();
    const int deg=basis.degree(patch,1-side.direction());
    const int max=basis.basisFunctionsOnSide(ps);

    GISMO_ASSERT(d==2,"only works for dimension 2");
    std::vector<patchSide> psides;
    pc1.getContainingSides(2,psides);
    psides.erase(std::remove(psides.begin(), psides.end(), ps), psides.end());
    const patchSide ps1 = psides[0];
    pc2.getContainingSides(2,psides);
    psides.erase(std::remove(psides.begin(), psides.end(), ps), psides.end());
    const patchSide ps2 = psides[0];

    const unsigned degree = std::min<int>(deg,basis.getIncrSmoothnessDegree()+1);
    const int dist=std::max<int>(degree,basis.getMinDist());
    unsigned val1,val2;

    _determineValues(ps,ps1,ps2,dist,degree,max,val1,val2,basis);
    parametricDistance1=basis.findParameter(interface.first(),pc1,val1);
    parametricDistance2=basis.findParameter(interface.first(),pc1,val2);
}

template<short_t d, class T>
void gsMPBESBasis<d,T>::distances::setParamDist(unsigned absoluteVal,const patchCorner& pc,const gsMPBESBasis<d,T>& basis)
{
    std::vector<patchCorner> cornerList;
    basis.getTopol().getCornerList(corner1,cornerList);
    if( std::find(cornerList.begin(), cornerList.end(), pc)!=cornerList.end() )
    {
        if(basis.isSpecialVertex(pc))
        {
            if(corner1.patch==interface.first().patch)
                parametricDistance1=basis.findParameter(interface.first(),corner1,absoluteVal);
            else
                parametricDistance1=basis.findParameter(interface.second(),corner1,absoluteVal);
        }
    }
    else
    {
        basis.getTopol().getCornerList(corner2,cornerList);
        if( std::find(cornerList.begin(), cornerList.end(), pc)!=cornerList.end() )
            if(basis.isSpecialVertex(pc))
            {
                if(corner2.patch==interface.first().patch)
                    parametricDistance2=basis.findParameter(interface.first(),corner2,absoluteVal);
                else
                    parametricDistance2=basis.findParameter(interface.second(),corner2,absoluteVal);
            }
    }
}

template<short_t d, class T>
T gsMPBESBasis<d,T>::distances::getParamDist(const patchCorner& pc,const gsMPBESBasis<d,T>& basis) const
{
    T param = -1.0;
    std::vector<patchCorner> cornerList;
    basis.getTopol().getCornerList(corner1,cornerList);
    if( std::find(cornerList.begin(), cornerList.end(), pc)!=cornerList.end() )
        param = parametricDistance1;
    else
    {
        basis.getTopol().getCornerList(corner2,cornerList);
        if( std::find(cornerList.begin(), cornerList.end(), pc)!=cornerList.end() )
            param = parametricDistance2;
    }
    return param;
}

template<short_t d, class T>
void gsMPBESBasis<d,T>::distances::_determineValues(patchSide side,patchSide ls,patchSide rs,int dist,unsigned degree,unsigned max,
                      unsigned& left,unsigned& right,const gsMPBESBasis<d,T>& basis) const
{
    gsVector<bool>pars(2);
    int dir = side.direction();
    int par = side.parameter();
    pars(dir) = par==0 ? false : true;
    pars(1-dir)=ls.parameter();
    patchCorner pc1 = patchCorner(side.patch,pars);
    pars(1-dir)=rs.parameter();
    patchCorner pc2 = patchCorner(side.patch,pars);
    bool leftSpecialVert = basis.isSpecialVertex(pc1);
    bool rightSpecialVert = basis.isSpecialVertex(pc2);
    bool leftInter = basis.getTopol().isInterface(ls);
    bool rightInter = basis.getTopol().isInterface(rs);
    unsigned leftDist = dist;
    unsigned rightDist = dist;
    if(leftSpecialVert&&rightSpecialVert)
    {
        if(leftDist+rightDist<=max)
        {
            left=leftDist;
            right=rightDist;
        }
        else
        {
            left=max/2;
            right=max-left;
        }
    }
    else if(leftSpecialVert&&rightInter)
    {
        right=degree;
        left=max-degree<leftDist?max-degree:leftDist;
    }
    else if(rightSpecialVert&&leftInter)
    {
        left=degree;
        right=max-degree<rightDist?max-degree:rightDist;
    }
    else if(rightSpecialVert)
    {
        right=max<rightDist?max:rightDist;
        left=0;
    }
    else if(leftSpecialVert)
    {
        left=max<leftDist?max:leftDist;
        right=0;
    }
    else
    {
        left=0;
        right=0;
        if(leftInter)
            left=degree;
        if(rightInter)
            right=degree;
    }
}

} // namespace gismo
