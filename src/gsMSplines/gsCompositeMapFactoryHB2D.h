/** @file gsHB2DMapMaker.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

// This class gsMapMaker has the sole purpose of creating a mapping of the type gsCompositeMapper

#pragma once

#include <gsCore/gsBasis.h>
#include <gsCore/gsBoxTopology.h>
#include <gsNurbs/gsKnotVector.h>

#include <gsMSplines/gsCompositeIncrSmoothnessBasis.h>
#include <gsMSplines/gsCompositeMapFactoryTensor.h>

#define TO_HTENSOR(x) static_cast<const gsHTensorBasis<d,T> *>(x)
#define TO_BSPLINE(x) static_cast<const gsTensorBSplineBasis<d,T> *>(x)

namespace gismo
{

/** @brief
      A univariate Lagrange basis.

      \tparam T coefficient type

      \ingroup basis
  */

template<short_t d,class T,class MapType>
class gsCompositeMapFactoryHB2D : public gsCompositeMapFactoryTensor<d,T,MapType>
{
    //static const int d = 2;
private:
    typedef gsBasis<T> BasisType;
    typedef gsCompositeMapFactoryTensor<d,T,MapType> Base;
public:
    gsCompositeMapFactoryHB2D(int incrSmoothnessDegree, gsBoxTopology * topol, gsCompositeIncrSmoothnessBasis<d,T> * basis) :
        Base(incrSmoothnessDegree,topol,basis)
    { }

    ~gsCompositeMapFactoryHB2D() { }

private:
    using Base::m_basis;
    using Base::m_incrSmoothnessDegree;
    using Base::m_topol;
    using Base::m_mapper;
    using Base::m_global;
    using Base::_setTensorMappingOfPatch;
    using Base::_getPatch;
    using Base::_getPatchIndex;
    using Base::_getLocalIndex;

private:
    //////////////////////////////////////////////////
    // general help functions for checking, finalizing and building of the mapping
    //////////////////////////////////////////////////

    bool _checkMapping() const
    {
        return true;
    }

    void _finalize() const
    {
        m_level = _getMaxLevel();
        gsSparseMatrix<real_t> mat=m_mapper->asMatrix();
        mat.conservativeResize(mat.rows(),m_global);
        delete m_mapper;
        m_mapper=new gsWeightMapper<T>(mat);
        m_mapper->optimize(gsWeightMapper<T>::optTargetToSource);
    }

    void _setMappingOfPatch(unsigned const patch) const
    {
        m_level=0;
        for(index_t i = 0;i<=_getMaxLevel();i++)
        {
            if(m_level<=_getMaxLevel(patch))
                _setTensorMappingOfPatch(patch);
            m_level++;
        }
    }

    index_t _getMaxLevel() const
    {
        index_t level = 0;
        for (size_t i = 0; i < m_basis->nPatches(); i++)
        {
            level = std::max(level, _getMaxLevel(i));
        }
        return level;
    }

    index_t _getMaxLevel(int patch) const
    {
        return TO_HTENSOR(&(m_basis->getBase(patch)))->maxLevel();
    }


    index_t getDistanceOfVertex(const patchCorner& pc,const patchSide& ps) const
    {
        std::vector<T> endpoints;
        T parametricDistance = m_basis->getParametricDistanceOfVertex(pc,ps);
        if(math::almostEqual<T>(parametricDistance,0.0))
            return 0;
        int patch = ps.patch;
        unsigned deg = m_basis->degree(patch,1-ps.direction());
        gsTensorBSplineBasis<d,T>* basis = TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[m_level];
        gsKnotVector<T> knots = basis->knots(1-(ps.direction()));
        gsVector<bool> pars;
        pc.parameters_into(d,pars);
        if(!pars(1-ps.direction()))
            knots.reverse();
        for(size_t i = deg+1;i<knots.size();i++)
            endpoints.push_back(knots.at(i));
        std::sort(endpoints.begin(),endpoints.end());
        unsigned nr=0;
        for(;nr<endpoints.size();nr++)
            if(math::almostEqual<T>(endpoints[nr],parametricDistance)||endpoints[nr]>=parametricDistance)
                break;
        return nr+1;
    }

private:
    //////////////////////////////////////////////////
    // functions calculating the weights for the mapping
    //////////////////////////////////////////////////

    gsKnotVector<T> _getKnotVector(unsigned const patch,unsigned const par) const
    {// todo: remove
        gsKnotVector<T> kvComp = TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[m_level]->knots(par);
        return gsKnotVector<T>(kvComp);
    }

    index_t _getParMax(unsigned patch,bool par) const
    {
        return TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[m_level]->size(par)-1;
    }

private:
    //////////////////////////////////////////////////
    // functions for working with Indexes
    //////////////////////////////////////////////////
    // localIndex = hierachical index of hierachical splines collected over all patches.

    bool _getLocalIndex_into(unsigned const patch,unsigned const u,unsigned const v,unsigned & localIndex) const
    {
        int patchIndex = _getPatchIndex(patch,u,v);
        localIndex=_getLocalIndex(patch,patchIndex);
        if(patchIndex==-1)
            return false;
        else
            return true;
    }

    index_t _getLocalIndex(unsigned const patch,unsigned u, unsigned v) const
    {
        return _getLocalIndex(patch,_getPatchIndex(patch,u,v));
    }

    index_t _getPatchIndex(unsigned const patch,boxSide const side,bool const flag) const
    {
        unsigned u,v;
        unsigned level = 0;
        gsVector<index_t,d> vec;
        int index, patchindex;
        do
        {
            if(level>TO_HTENSOR(&(m_basis->getBase(patch)))->maxLevel())
                GISMO_ERROR("did not find the patchindex");

            const unsigned u_amount=TO_HTENSOR(&(m_basis->getBase(patch)))->tensorLevel(level).size(0);
            const unsigned v_amount=TO_HTENSOR(&(m_basis->getBase(patch)))->tensorLevel(level).size(1);
            if(side.direction())
                if(side.parameter())
                {
                    u=flag ? (u_amount-1) : 0;
                    v=v_amount-1;
                }
                else
                {
                    u=flag ? (u_amount-1) :0;
                    v=0;
                }
            else
                if(side.parameter())
                {
                    u=u_amount-1;
                    v=flag ? (v_amount-1) : 0;
                }
                else
                {
                    u=0;
                    v=flag ? (v_amount-1) : 0;
                }
            vec(0)=u,vec(1)=v;
            index=TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[level]->index(vec);
            patchindex=TO_HTENSOR(&(m_basis->getBase(patch)))->flatTensorIndexToHierachicalIndex(index,level);
            level++;
        }while(-1==patchindex);
        return patchindex;
    }

    index_t _getPatchIndex(unsigned const patch,unsigned const u,unsigned const v) const
    {
        unsigned index=_getTensorIndex(patch,u,v);
        return TO_HTENSOR(&(m_basis->getBase(patch)))->flatTensorIndexToHierachicalIndex(index,m_level);
    }

    // tensorIndex = flat tensor index of one level of hierarchical splines
    index_t _getTensorIndex(unsigned const patch,unsigned const u, unsigned const v) const
    {
        gsVector<index_t,d> vec;
        vec(0)=u,vec(1)=v;
        return TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[m_level]->index(vec);
    }

    index_t _getTensorIndex(unsigned const patch,unsigned const patchIndex) const
    {
        unsigned combIndex = TO_HTENSOR(&(m_basis->getBase(patch)))->flatTensorIndexOf(patchIndex,m_level);
//        for(unsigned i = 0;i<m_level;i++)
//            combIndex-= TO_HTENSOR((*m_bases)[patch])->m_bases[i]->size();
        return combIndex;
    }

    index_t _getPar(index_t localIndex,bool par) const
    {
        unsigned patch = _getPatch(localIndex);
        unsigned patchIndex = _getPatchIndex(localIndex);
        return _getPar(patch,_getTensorIndex(patch,patchIndex),par);
    }

    index_t _getPar(unsigned patch,unsigned tensorIndex, bool par) const
    {
        gsVector<index_t,d> vec = TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[m_level]->tensorIndex(tensorIndex);
        GISMO_ASSERT(static_cast<int>(vec(par))<TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[m_level]->size(par),"wrong tensorIndex");
        GISMO_ASSERT(static_cast<int>(vec(!par))<TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[m_level]->size(!par),"wrong tensorIndex");
        return vec(par);
    }

private:
    mutable index_t m_level; // used in the construction to determine the level, after the construction it is used to say the max_level of the H-Splines
};

}

#undef TO_BSPLINE
#undef TO_HTENSOR
