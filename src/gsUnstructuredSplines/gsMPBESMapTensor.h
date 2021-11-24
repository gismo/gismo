/** @file gsMPBESMapTensor.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

// This class gsMPBESMapTensor has the sole purpose of creating a mapping of the type gsWeightMapper

#pragma once

#include <gsNurbs/gsBoehm.h>
#include <gsCore/gsBoxTopology.h>

namespace gismo
{
/** @brief
      A univariate Lagrange basis.

      \tparam T coefficient type

      \ingroup basis
  */

template<short_t d,class T>
class gsMPBESMapTensor
{
private:

    typedef T weightType;
    typedef index_t indexType; //indizes of gsMatrix
    typedef std::vector<std::pair<int,int> >::iterator step_iter;
    typedef gsBasis<T> BasisType;
    typedef typename std::vector<BasisType *>::const_iterator ConstBasisIter;
    typedef typename std::vector<BasisType *>::iterator BasisIter;
    typedef typename std::vector<gsMatrix<T> *>::const_iterator ConstMatrixPtrIter;
    typedef typename gsSparseMatrix<weightType,0,indexType>::InnerIterator InIterMat;
    typedef std::vector<indexType> IndexContainer;
    typedef std::vector<indexType>::const_iterator ConstIndexIter;
    //int m_dim = 2;

public:
    gsMPBESMapTensor(int incrSmoothnessDegree, gsBoxTopology * topol, gsMPBESBasis<d,T> * basis) :
        m_incrSmoothnessDegree(incrSmoothnessDegree), m_topol(topol), m_basis(basis), m_global(0)
    {
        m_mapper=NULL;
    }

    virtual ~gsMPBESMapTensor() { }

    //////////////////////////////////////////////////
    // Virtual member functions required by the base class
    //////////////////////////////////////////////////

    gsWeightMapper<T> * makeMapper() const
    {
        unsigned locals = _getNrOfLocalBasisFunktions();
        unsigned size = m_basis->nPatches();
        std::vector<index_t> offsets;
        for(unsigned i = 0;i<size;i++)
            offsets.push_back(_getFirstLocalIndex(i));
        m_mapper = new gsWeightMapper<T>(locals,locals);
        for(unsigned i = 0;i<size;i++)
            _setMappingOfPatch(i);
        _finalize();
        gsWeightMapper<T> * mapper = m_mapper;
        m_mapper = NULL;
        return mapper;
    }

protected:
    //////////////////////////////////////////////////
    // general help functions for checking, finalizing and building of the mapping
    //////////////////////////////////////////////////

    virtual bool _checkMapping() const = 0;

    virtual void _finalize() const = 0;

    virtual void _setMappingOfPatch(unsigned const patch) const = 0;

    int _getNrOfLocalBasisFunktions() const
    {
        int size = 0;
        for (size_t i = 0; i < m_basis->nPatches(); ++i)
            size += m_basis->getBase(i).size();
        return size;
    }

protected:
    //////////////////////////////////////////////////
    // function to construct a mapping of a tensorstructured patch
    //////////////////////////////////////////////////

    virtual index_t getDistanceOfVertex(const patchCorner& pc,const patchSide& ps) const = 0;

    void _setTensorMappingOfPatch(unsigned const patch) const
    {
        const gsBasis<T> & base = m_basis->getBase(patch);
        unsigned degree_u = std::min(base.degree(0),(short_t)(m_incrSmoothnessDegree+1));
        unsigned degree_v = std::min(base.degree(1),(short_t)(m_incrSmoothnessDegree+1));
        unsigned u_amount =_getParMax(patch,0)+1, v_amount=_getParMax(patch,1)+1;
        patchSide ps_north(patch,boundary::north);
        patchSide ps_south(patch,boundary::south);
        patchSide ps_east(patch,boundary::east);
        patchSide ps_west(patch,boundary::west);
        bool north = m_topol->isInterface(ps_north);
        bool south = m_topol->isInterface(ps_south);
        bool east  = m_topol->isInterface(ps_east);
        bool west  = m_topol->isInterface(ps_west);
        patchCorner pc_northeast(patch,boundary::northeast);
        patchCorner pc_northwest(patch,boundary::northwest);
        patchCorner pc_southeast(patch,boundary::southeast);
        patchCorner pc_southwest(patch,boundary::southwest);
        bool northeast = m_basis->isSpecialVertex(pc_northeast);
        bool northwest = m_basis->isSpecialVertex(pc_northwest);
        bool southeast = m_basis->isSpecialVertex(pc_southeast);
        bool southwest = m_basis->isSpecialVertex(pc_southwest);

        unsigned ne_l_u = getDistanceOfVertex(pc_northeast,ps_north);
        unsigned nw_l_u = getDistanceOfVertex(pc_northwest,ps_north);
        unsigned se_l_u = getDistanceOfVertex(pc_southeast,ps_south);
        unsigned sw_l_u = getDistanceOfVertex(pc_southwest,ps_south);
        unsigned se_l_v = getDistanceOfVertex(pc_southeast,ps_east);
        unsigned ne_l_v = getDistanceOfVertex(pc_northeast,ps_east);
        unsigned sw_l_v = getDistanceOfVertex(pc_southwest,ps_west);
        unsigned nw_l_v = getDistanceOfVertex(pc_northwest,ps_west);

        if(ne_l_u+nw_l_u>u_amount)
        {
            ne_l_u=u_amount/2;
            nw_l_u=u_amount-ne_l_u;
        }
        if(se_l_u+sw_l_u>u_amount)
        {
            se_l_u=u_amount/2;
            sw_l_u=u_amount-se_l_u;
        }
        if(se_l_v+ne_l_v>v_amount)
        {
            se_l_v=v_amount/2;
            ne_l_v=v_amount-se_l_v;
        }
        if(sw_l_v+nw_l_v>v_amount)
        {
            sw_l_v=v_amount/2;
            nw_l_v=v_amount-sw_l_v;
        }

        unsigned n_l_v = north ? degree_v : 0;
        unsigned s_l_v = south ? degree_v : 0;
        unsigned e_l_u = east ? degree_u : 0;
        unsigned w_l_u = west ? degree_u : 0;
        // added for one-patch gluing
        bool north_same=false;
        bool east_same=false;
        if(1==m_basis->nPatches())
        {
            north_same = south;
            east_same = east;
        }

        unsigned localIndex;
        // south-west
        for(unsigned j = 0;j<sw_l_v ;j++)
            for(unsigned i = 0;i<sw_l_u ;i++)
            {
                if(!_getLocalIndex_into(patch,i,j,localIndex))
                    continue;
                if(i>=degree_u&&j>=degree_v)
                    ;
                else if((southwest||degree_u==1||degree_v==1) && i == 0 && j==0)
                    _addSingularCornerToMap(localIndex);
                else if(southwest && west && i==0)
                    _addCombinedLineToMap(ps_west,localIndex,1,0);
                else if(southwest && south && j==0)
                    _addCombinedLineToMap(ps_south,localIndex,1,0);
                else if(!southwest && south && west)
                    _addCombinedBlockToMap(ps_south,ps_west,localIndex,degree_u,degree_v,j,i);
                else if(!southwest && south && !west)
                    _addCombinedLineToMap(ps_south,localIndex,degree_v,j);
                else if(!southwest && !south && west)
                    _addCombinedLineToMap(ps_west,localIndex,degree_u,i);
                else
                    _addFreeToMap(localIndex);
            }
        // south-middle
        for(unsigned j=0;j<s_l_v;j++)
            for(unsigned i=sw_l_u;i<u_amount-se_l_u;i++)
            {
                if(!_getLocalIndex_into(patch,i,j,localIndex))
                    continue;
                if(south)
                    _addCombinedLineToMap(ps_south,localIndex,degree_v,j);
                else
                    _addFreeToMap(localIndex);
            }
        // south-east
        for(unsigned j = 0;j<se_l_v ;j++)
            for(unsigned i = u_amount-se_l_u; i<u_amount ;i++)
            {
                if(!_getLocalIndex_into(patch,i,j,localIndex))
                    continue;
                if(i<u_amount-degree_u&&j>=degree_v)
                    continue;
                else if(!southeast&&east_same)
                    continue;
                else if((southeast||degree_u==1||degree_v==1) && i == u_amount-1 && j==0)
                    _addSingularCornerToMap(localIndex);
                else if(southeast && east && i==u_amount-1)
                    _addCombinedLineToMap(ps_east,localIndex,1,0);
                else if(southeast && south && j==0)
                    _addCombinedLineToMap(ps_south,localIndex,1,0);
                else if(!southeast && south && east)
                    _addCombinedBlockToMap(ps_south,ps_east,localIndex,degree_u,degree_v,j,u_amount-1-i);
                else if(!southeast && south && !east)
                    _addCombinedLineToMap(ps_south,localIndex,degree_v,j);
                else if(!southeast && !south && east)
                    _addCombinedLineToMap(ps_east,localIndex,degree_u,u_amount-1-i);
                else
                    _addFreeToMap(localIndex);
            }
        // west-middle
        for(unsigned j=sw_l_v;j<v_amount-nw_l_v;j++)
            for(unsigned i=0;i<w_l_u;i++)
            {
                if(!_getLocalIndex_into(patch,i,j,localIndex))
                    continue;
                if(west)
                    _addCombinedLineToMap(ps_west,localIndex,degree_u,i);
                else
                    _addFreeToMap(localIndex);
            }
        // middle-middle
        for(unsigned j=s_l_v;j<v_amount-n_l_v;j++)
            for(unsigned i=w_l_u;i<u_amount-e_l_u;i++)
            {
                if(j<sw_l_v&&j<degree_v&&i<sw_l_u&&i<degree_u)
                    continue;
                if(j<se_l_v&&j<degree_v&&i>=u_amount-se_l_u&&i>=u_amount-degree_u)
                    continue;
                if(j>=v_amount-ne_l_v&&j>=v_amount-degree_v&&i>=u_amount-ne_l_u&&i>=u_amount-degree_u)
                    continue;
                if(j>=v_amount-nw_l_v&&j>=v_amount-degree_v&&i<nw_l_u&&i<degree_u)
                    continue;
                if(!_getLocalIndex_into(patch,i,j,localIndex))
                    continue;
                _addFreeToMap(localIndex);
            }
        // east-middle
        for(unsigned j=se_l_v;j<v_amount-ne_l_v;j++)
            for(unsigned i=u_amount-e_l_u;i<u_amount;i++)
            {
                if(!_getLocalIndex_into(patch,i,j,localIndex))
                    continue;
                if(east&&east_same)
                    continue;
                if(east)
                    _addCombinedLineToMap(ps_east,localIndex,degree_u,u_amount-1-i);
                else
                    _addFreeToMap(localIndex);
            }
        // north-west
        for(unsigned j = v_amount-nw_l_v;j<v_amount ;j++)
            for(unsigned i = 0; i<nw_l_u ;i++)
            {
                if(!_getLocalIndex_into(patch,i,j,localIndex))
                    continue;
                if(i>=degree_u&&j<v_amount-degree_v)
                    continue;
                else if(!northwest&&north_same)
                    continue;
                else if((northwest||degree_u==1||degree_v==1) && i == 0 && j==v_amount-1)
                    _addSingularCornerToMap(localIndex);
                else if(northwest && west && i==0)
                    _addCombinedLineToMap(ps_west,localIndex,1,0);
                else if(northwest && north && j==v_amount-1)
                    _addCombinedLineToMap(ps_north,localIndex,1,0);
                else if(!northwest && north && west)
                    _addCombinedBlockToMap(ps_north,ps_west,localIndex,degree_u,degree_v,v_amount-1-j,i);
                else if(!northwest && north && !west)
                    _addCombinedLineToMap(ps_north,localIndex,degree_v,v_amount-1-j);
                else if(!northwest && !north && west)
                    _addCombinedLineToMap(ps_west,localIndex,degree_u,i);
                else
                    _addFreeToMap(localIndex);
            }
        // north-middle
        for(unsigned j=v_amount-n_l_v;j<v_amount;j++)
            for(unsigned i=nw_l_u;i<u_amount-ne_l_u;i++)
            {
                if(!_getLocalIndex_into(patch,i,j,localIndex))
                    continue;
                if(north&&north_same)
                    continue;
                if(north)
                    _addCombinedLineToMap(ps_north,localIndex,degree_v,v_amount-1-j);
                else
                    _addFreeToMap(localIndex);
            }
        // north-east
        for(unsigned j = v_amount-ne_l_v;j<v_amount ;j++)
            for(unsigned i = u_amount-ne_l_u; i<u_amount ;i++)
            {
                if(!_getLocalIndex_into(patch,i,j,localIndex))
                    continue;
                if(i<u_amount-degree_u&&j<v_amount-degree_v)
                    continue;
                else if(!northeast&&(north_same||east_same))
                    continue;
                else if((northeast||degree_u==1||degree_v==1) && i == u_amount-1 && j==v_amount-1)
                    _addSingularCornerToMap(localIndex);
                else if(northeast && east && i==u_amount-1)
                    _addCombinedLineToMap(ps_east,localIndex,1,0);
                else if(northeast && north && j==v_amount-1)
                    _addCombinedLineToMap(ps_north,localIndex,1,0);
                else if(!northeast && north && east)
                    _addCombinedBlockToMap(ps_north,ps_east,localIndex,degree_u,degree_v,v_amount-1-j,u_amount-1-i);
                else if(!northeast && north && !east)
                    _addCombinedLineToMap(ps_north,localIndex,degree_v,v_amount-1-j);
                else if(!northeast && !north && east)
                    _addCombinedLineToMap(ps_east,localIndex,degree_u,u_amount-1-i);
                else
                    _addFreeToMap(localIndex);
            }
//        gsVector<unsigned> point(m_dim);
//        std::vector<index_t> lengths(m_dim);
//        lengths.push_back(degree_u);
//        lengths.push_back(degree_v);
//        std::vector<index_t> maxima(m_dim);
//        maxima.push_back(u_amount);
//        maxima.push_back(v_amount);
//        std::vector<index_t> distances(m_dim);
//        point.setZero();
//        do
//        {
//            if(!_getLocalIndex_into(patch,point[0],point[1],localIndex))
//                continue;
//            if(_isInMiddle(point,lengths,maxima))
//                addFreeToMap();
//            _getDistancesToBoundary(point,distances);
//            for(unsigned i = 0;i<point.rows();++i)
//            {

//                //dists
//                //sides
//            }
//            _addBlockToMap(localIndex, sides, lengths, dists);

//            if(m_topol->getDistancesToBoundary(patch,point,lengths,distances))
//                addFreeToMap();
//            else if(m_topol->pointIsNearCp(patch,point))
//                addCombinedLineToMap();
//            else if(m_topol->pointIsNextToSmoothCorner(patch,point))
//                addCombinedBlockToMap();
//            else
//                addFreeToMap();
//        }while(_nextPoint(lengths,point));
    }

//    bool _isInMiddle(const gsVector<unsigned>& point,const std::vector<index_t>& lengths,
//                     const std::vector<index_t>& maxima)
//    {
//        for(unsigned i = 0;i<point.rows();++i)
//        {
//            if(point(i)<lengths[i]||maxima[i]-lengths[i]<point(i))
//                return false;
//        }
//        return true;
//    }

//    void _getDistancesToBoundary(const gsVector<unsigned>& point, std::vector<index_t>& dists)
//    {
//        for(unsigned i = 0;i<point.rows();++i)
//        {
//            dists.push_back(point((i+1)%point.rows()));
//        }
//    }
private:
    //////////////////////////////////////////////////
    // functions calculating the weights for the mapping
    //////////////////////////////////////////////////

    void _getBoehmCoefs(gsKnotVector<T> kv0i,bool dir0,T weight0,gsKnotVector<T> kv1i,bool dir1,T weight1,unsigned degree,unsigned knot_index,std::vector<T> & result) const
    {
        result.clear();
        if(degree==1)
        {
            result.push_back(1.0);
            result.push_back(1.0);
            return;
        }
        typedef typename std::vector<T>::iterator tIter;
        std::vector<T> temp0 = kv0i.get();
        std::vector<T> temp1 = kv1i.get();
        for(tIter it = temp0.begin();it!=temp0.end();++it)
            (*it)*=weight0;
        for(tIter it = temp1.begin();it!=temp1.end();++it)
            (*it)*=weight1;
        gsKnotVector<T> kv0(give(temp0), kv0i.degree());
        gsKnotVector<T> kv1(give(temp1), kv1i.degree());

        unsigned l0 = kv0.size();
        if (!dir0)
          kv0.reverse();
        if (dir1)
          kv1.reverse();
        kv1.addConstant(kv0.last());
        T inserted = kv0.last();
        kv0.remove(kv0.last(),degree);
        kv1.remove(kv1.first(),kv1.degree()+1);
        kv0.append(kv1.begin(),kv1.end());
        gsMatrix<T> boem_mat;
        boem_mat.setZero(kv0.size()-kv0.degree()-1,1);
        boem_mat.coeffRef(l0-(kv1.degree()+2)-knot_index)=1;
        gsBoehm<T,gsKnotVector<T>,gsMatrix<T> >(kv0,boem_mat,inserted,degree-1,false);
        for(unsigned i = 0;i<degree;i++)
        {
            //unsigned ind=dir0 ? boem_mat.rows()-1-i : i;
            result.push_back(boem_mat.coeff(i+l0-(kv0.degree()+2)-knot_index,0));
            if(i==knot_index)
                result.push_back(boem_mat.coeff(i+l0-(kv0.degree()+2)-knot_index,0));
        }
    }

protected:
    virtual gsKnotVector<T> _getKnotVector(unsigned const patch,unsigned const par) const = 0;

    virtual index_t _getParMax(unsigned patch,bool par) const = 0;

private:
    //////////////////////////////////////////////////
    // functions to look for one vertex in all the patches
    //////////////////////////////////////////////////

    bool _alreadyHandled(patchSide const & ps,bool flag) const
    {
        std::vector<indexType> vertices = _findVertexInPatches(ps,flag);
        for(unsigned i = 0;i<vertices.size();i++)
            if(_getPatch(vertices[i]) < ps.patch )
                return true;
        return false;
    }

    bool _alreadyHandled(std::vector<index_t> local_BFs,unsigned startPatch) const
    {
        for(size_t i = 0;i<local_BFs.size();i++)
            if((index_t)startPatch>_getPatch(local_BFs[i]))
                return true;
        return false;
    }

    std::vector<indexType> _findVertexInPatches(patchSide ps, bool flag) const
    {
        gsVector<bool> v(2);
        v(ps.direction())=ps.parameter();
        v(1-ps.direction())=flag;
        patchCorner start(ps.patch,boxCorner(v));
        std::vector<patchCorner> cornerList;
        m_topol->getCornerList(start,cornerList);
        std::vector<indexType> vertices;
        for(unsigned i = 0;i<cornerList.size();++i)
            vertices.push_back(_getLocalIndex(cornerList[i]));
        return vertices;
    }

    void _flipOverCorner(patchSide & ps,bool & flag) const
    {
        int patch = ps.patch;
        switch(ps.side())
        {
        case boundary::north : ps = patchSide(patch,flag ? boundary::east : boundary::west);
            flag = 1;
            break;
        case boundary::south : ps = patchSide(patch,flag ? boundary::east : boundary::west);
            flag = 0;
            break;
        case boundary::west : ps = patchSide(patch,flag ? boundary::north : boundary::south);
            flag = 0;
            break;
        case boundary::east : ps = patchSide(patch,flag ? boundary::north : boundary::south);
            flag = 1;
            break;
        default :
            GISMO_ERROR("patchSide has no valid side.");
        }
    }

private:
    //////////////////////////////////////////////////
    // functions for adding new map entries
    //////////////////////////////////////////////////


    void _addToMap(std::vector<indexType> indices,std::vector<weightType> weights) const
    {
        for(unsigned i = 0;i<indices.size();++i)
            m_mapper->setEntry(indices[i],m_global,weights[i]);
        m_global++;
    }

    void _addFreeToMap(unsigned localIndex) const
    {
        std::vector<indexType> indices;
        indices.push_back(localIndex);
        std::vector<weightType> weights;
        weights.push_back(1);
        _addToMap(indices,weights);
    }

    void _addSingularCornerToMap(unsigned localIndex) const
    {
        unsigned u=_getPar(localIndex,0), v=_getPar(localIndex,1), patch=_getPatch(localIndex);
        bool flag=0;
        patchSide ps;
        if(u==0&&v==0)
            ps=patchSide(patch,boundary::south);
        else if(u==0)
            ps=patchSide(patch,boundary::north);
        else if(v==0)
            ps=patchSide(patch,boundary::east);
        else
        {
            ps=patchSide(patch,boundary::north);
            flag=1;
        }
        if(_alreadyHandled(ps,flag))
            return;
        std::vector<indexType> vertices = _findVertexInPatches(ps,flag);
        std::vector<indexType> locals;
        std::vector<weightType> weights;
        for(unsigned i = 0;i<vertices.size();i++)
        {
            locals.push_back(vertices[i]);
            weights.push_back(static_cast<weightType>(1));
        }
        _addToMap(locals,weights);
    }

    void _addCombinedLineToMap(patchSide & ps,unsigned localStartIndex,unsigned length,int distToPs) const
    {
        std::vector<patchSide> sides;
        std::vector<index_t> lengths;
        std::vector<index_t> dists;
        sides.push_back(ps);
        lengths.push_back(length);
        dists.push_back(distToPs);
        _addBlockToMap(localStartIndex,sides,lengths,dists);
    }

    void _addCombinedBlockToMap(patchSide & ps_u,patchSide & ps_v,unsigned localStartIndex,unsigned length_u,unsigned length_v,int distToPs_u, int distToPs_v) const
    {
        std::vector<patchSide> sides;
        std::vector<index_t> lengths;
        std::vector<index_t> dists;
        sides.push_back(ps_u);
        sides.push_back(ps_v);
        lengths.push_back(length_u);
        lengths.push_back(length_v);
        dists.push_back(distToPs_u);
        dists.push_back(distToPs_v);
        _addBlockToMap(localStartIndex,sides,lengths,dists);
    }

    void _addBlockToMap(unsigned localStartIndex, std::vector<patchSide> const & sides, std::vector<index_t> const & lengths, std::vector<index_t> const & dists) const
    {
        patchSide ps,ps_neighbour;
        bool par,par_neigh;
        int dir,dir_neigh;
        gsKnotVector<T> kv_neigh,kv_this;
        std::vector<std::vector<T> > weights1D;
        std::vector<bool> minus_dir;
        std::vector<index_t> dirs;
        std::vector<weightType> weight;
        T w0,w1;
        for(unsigned i = 0; i<sides.size();++i)
        {
            ps=sides[i];
            if(!m_topol->getNeighbour(ps,ps_neighbour))
                GISMO_ERROR("no neighbour found on this side.");
            par = ps.parameter();
            dir = ps.direction();
            minus_dir.push_back(!par);
            dirs.push_back(dir);
            kv_this=_getKnotVector(ps.patch,dir);
            par_neigh = ps_neighbour.parameter();
            dir_neigh = ps_neighbour.direction();
            kv_neigh=_getKnotVector(ps_neighbour.patch,dir_neigh);
            w0 = m_basis->getWeight(ps);
            w1 = m_basis->getWeight(ps_neighbour);
            weight.clear();
            _getBoehmCoefs(kv_this,par,w0,kv_neigh,par_neigh,w1,lengths[i],dists[i],weight);
            weights1D.push_back(weight);
        }
        std::vector<indexType> locals;
        std::vector<weightType> weights;
        std::vector<std::pair<int,int> >steps;
        gsVector<index_t> distances(sides.size());
        distances.setZero();
        do
        {
            steps.clear();
            weightType weight2 = 1.0;
            for(unsigned i=0;i<sides.size();++i)
            {
                steps.push_back(std::make_pair(minus_dir[i] ? -distances(i) : distances(i),dirs[i]));
                weight2 *= weights1D[i][distances(i)];
            }
            unsigned localIndexTravelled = _travelUVSteps(localStartIndex,steps);
            locals.push_back(localIndexTravelled);
            weights.push_back(weight2);
        }while(_nextPoint(lengths,distances));
        if(_alreadyHandled(locals,_getPatch(localStartIndex)))
            return;
        _addToMap(locals,weights);
    }

    bool _nextPoint(std::vector<index_t> const & lengths, gsVector<index_t> & point) const
    {
        for(unsigned i = 0; i<lengths.size(); ++i)
        {
            if(++point(i)<=lengths[i])
            {
                return true;
            }
            else
                point(i)=0;
        }
        return false;
    }

private:
    //////////////////////////////////////////////////
    // functions for travelling through the basis-function indizes of patches
    //////////////////////////////////////////////////

    unsigned _transformToNeighbour(unsigned localIndex,step_iter current,step_iter end) const
    {
        index_t patch = _getPatch(localIndex);
        index_t u = _getPar(localIndex,0);
        index_t v = _getPar(localIndex,1);
        patchSide ps,ps_neighbour;
        if(current->second==0)
            ps = patchSide(patch,current->first>0 ? 2 :  1);
        else
            ps = patchSide(patch,current->first>0 ? 4 : 3);
        boundaryInterface bf;
        if(!m_topol->getNeighbour(ps,ps_neighbour)||!m_topol->getInterface(ps,bf))
            GISMO_ERROR("ps no interface");
        GISMO_ASSERT((ps.side()==1&&u==0)||(ps.side()==2&&u==_getParMax(patch,0))||
                     (ps.side()==3&&v==0)||(ps.side()==4&&v==_getParMax(patch,1)),
                     "steps does not fit to point");
        GISMO_ASSERT(current->first!=0,"stepsize is zero, which is not allowed in this position.");
        bool orient = bf.dirOrientation(ps,1-ps.direction());
        unsigned non_fixed = current->second==1 ? u : v;
        current->first>0 ? current->first-- : current->first++;
        _transformStepsToNeighbour(current,end,ps,ps_neighbour,orient);
        localIndex=_transformIndexToNeighbour(localIndex,ps_neighbour,orient,non_fixed);
        return localIndex;
    }

    unsigned _transformIndexToNeighbour(unsigned localIndex,patchSide const & ps_neighbour,bool orient, int non_fixed) const
    {
        unsigned patch = ps_neighbour.patch;
        int u_max = _getParMax(patch,0);
        int v_max = _getParMax(patch,1);
        if((ps_neighbour.direction()))
            localIndex = _getLocalIndex(patch,orient ? non_fixed : u_max-non_fixed,ps_neighbour.parameter() ? v_max : 0);
        else
            localIndex = _getLocalIndex(patch,(ps_neighbour.parameter()) ? u_max : 0,orient ? non_fixed : v_max-non_fixed);
        return localIndex;
    }
    void _transformStepsToNeighbour(step_iter current,step_iter end,patchSide const & ps,patchSide const & ps_neighbour,bool orient) const
    {
        boxSide s1=ps.side(),s2=ps_neighbour.side();
        if(s1==s2)
            for(step_iter it = current;it!=end;++it)
            {
                if(it->second==s1.direction())
                    it->first *= -1;
                else
                    if(!orient)
                        it->first *= -1;
            }
        else if((s1==1&&s2==3)||(s1==2&&s2==4)||(s1==3&&s2==1)||(s1==4&&s2==2))
            for(step_iter it = current;it!=end;++it)
            {
                if(it->second==s1.direction())
                    it->first *= -1;
                else
                    if(!orient)
                        it->first *= -1;
                it->second=(1-it->second);
            }
        else if((s1==1&&s2==4)||(s1==2&&s2==3)||(s1==3&&s2==2)||(s1==4&&s2==1))
            for(step_iter it = current;it!=end;++it)
            {
                if(!(it->second==s1.direction()))
                    if(!orient)
                        it->first *= -1;
                it->second=(1-it->second);
            }
    }

    unsigned _travelInsidePatch(unsigned localIndex,bool dir,step_iter current) const
    {
        unsigned patch = _getPatch(localIndex);
        int par = _getPar(localIndex,dir);
        int par_max = _getParMax(patch,dir);
        int fixed_par = _getPar(localIndex,!dir);
        if(0<=par+current->first&&par+current->first<=par_max)
        {
            par+=current->first;
            current->first=0;
        }
        else if(current->first>0)
        {
            current->first-=par_max-par;
            par=par_max;
        }
        else
        {
            current->first+=par;
            par=0;
        }
        localIndex = _getLocalIndex(patch,dir ? fixed_par : par,dir ? par : fixed_par);
        return localIndex;
    }

    unsigned _travelUVSteps(unsigned localIndex,std::vector<std::pair<int, int> > steps) const
    {
        step_iter current = steps.begin();
        step_iter end = steps.end();
        while(current!=end)
        {
            GISMO_ASSERT(current->second==0||current->second==1,"only 0 and 1 direction possible");
            localIndex=_travelInsidePatch(localIndex, current->second != 0, current);
            if(current->first==0)
                ++current;
            else
                localIndex=_transformToNeighbour(localIndex,current,end);
        }
        return localIndex;
    }

protected:
    //////////////////////////////////////////////////
    // functions for working with Indexes
    //////////////////////////////////////////////////

    //localIndex: this is the index going over all local basis functions of all patches
    virtual bool _getLocalIndex_into(unsigned const patch,unsigned const u,unsigned const v,unsigned & localIndex) const=0;

    index_t _getLocalIndex(unsigned const patch,unsigned const patchIndex) const
    {
        return _getFirstLocalIndex(patch)+patchIndex;
    }

    index_t _getLocalIndex(patchCorner const & pc) const
    {
        gsVector<bool> param = pc.parameters(d);
        patchSide ps = param(1) ? patchSide(pc.patch,boxSide(4)) : patchSide(pc.patch,boxSide(3));
        return _getLocalIndex(ps,param(0));
    }

    index_t _getLocalIndex(patchSide const & ps,bool const flag) const
    {
        return _getLocalIndex(ps.patch,ps.side(),flag);
    }

    index_t _getLocalIndex(unsigned const patch,boxSide const side,bool const flag) const
    {
        return _getLocalIndex(patch,_getPatchIndex(patch,side,flag));
    }

    virtual index_t _getLocalIndex(unsigned const patch,unsigned u,unsigned v) const = 0;

    index_t _getFirstLocalIndex(unsigned const patch) const
    {
        unsigned index=0;
        for(unsigned i=0;i<patch;i++)
        {
            index+=m_basis->localSize(i);
        }
        return index;
    }

    index_t _getLastLocalIndex(unsigned const patch) const
    {
        return _getFirstLocalIndex(patch)+m_basis->localSize(patch)-1;
    }

    //patchIndex: this is the index going over all local basis functions of one patch
    virtual index_t _getPatchIndex(unsigned const patch,boxSide const side,bool const flag) const = 0;

    index_t _getPatchIndex(unsigned localIndex) const
    {
        unsigned patchIndex=localIndex;
        for(index_t i = 0;i<_getPatch(localIndex);i++)
        {
            patchIndex-=m_basis->localSize(i);
        }
        return patchIndex;
    }

    index_t _getPatch(unsigned localIndex) const
    {
        size_t patch;
        for (patch = 0; patch < m_basis->nPatches(); patch++)
        {
            if (localIndex >= m_basis->localSize(patch))
                localIndex -= m_basis->localSize(patch);
            else
                break;
        }
        return patch;
    }

    virtual index_t _getPar(index_t localIndex,bool par) const = 0;

protected:

    short_t const m_incrSmoothnessDegree;
    gsBoxTopology const *m_topol;
    gsMPBESBasis<d,T> const *m_basis;
    mutable gsWeightMapper<T> *m_mapper;
    mutable unsigned m_global;
};// class gsMPBESMapTensor

}
