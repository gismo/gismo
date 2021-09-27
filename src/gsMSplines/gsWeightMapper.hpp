/** @file gsWeightMapper.cpp

    @brief Provides implementation of gsWeightMapper class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, F. Buchegger
*/

#include <gsMSplines/gsWeightMapper.h>

namespace gismo
{

namespace
{
void increment(const index_t ** ptr){++(*ptr);}

template <class T, typename matrixT>
void getNonZeroInnerOfMatrix ( matrixT const &mat, typename gsWeightMapper<T>::IndexContainer const & outer, typename gsWeightMapper<T>::IndexContainer & inner)
{
    const index_t numSource=outer.size();
    std::vector<const index_t*>   nzIndices(numSource);
    std::vector<const index_t*>   end(numSource);
    std::vector<const index_t**>  nzMin;
    nzMin.reserve(numSource);
    inner.clear();
    inner.reserve(10*numSource); // reserve more memory then expected to be used

    // init the vectors
    for (index_t i=0; i<numSource;++i)
    {
        nzIndices[i]=mat.innerIndexPtr()+mat.outerIndexPtr()[outer[i]];
        end[i]=mat.innerIndexPtr()+mat.outerIndexPtr()[outer[i]+1];
    }

    const index_t idMax=mat.cols();
    index_t curMin;

    while (true) {
        curMin=idMax;
        nzMin.clear();
        for(index_t i=0; i<numSource;++i)
        {
            if(nzIndices[i]<end[i])
            {
                if (*nzIndices[i]<curMin)
                {
                    curMin=*nzIndices[i];
                    nzMin.clear();
                    nzMin.push_back(&nzIndices[i]);
                }
                else if (*nzIndices[i]==curMin)
                {
                    nzMin.push_back(&nzIndices[i]);
                }
            }
        }
        if (curMin>=idMax)
            break;
        inner.push_back(curMin);
        std::for_each(nzMin.begin(),nzMin.end(),increment);
    }
}

}

template<class T>
void gsWeightMapper<T>::mapToTargetCoefs(gsMatrix<weightType> const & sourceCoefs,gsMatrix<weightType> & targetCoefs) const
{
    // THIS DOES NOT WORK AS EIGEN SOLVERS ARE NOT ABLE TO SOLVE A ROW MAJOR MATRIX (last checked EIGEN::3.4)
    //    GISMO_ASSERT(m_matrix.isCompressed(),"Mapping has to be compressed for this.");
    //    gsSparseSolver<gsWeightMapper<T>::LToGMatrix>::QR solver;
    //    solver.compute(m_matrix);
    //    if(solver.info()!=Eigen::Success)
    //        GISMO_ERROR("Could not compute the solver SparseQR");
    //    targetCoefs=solver.solve(sourceCoefs);
    //    if(solver.info()!=Eigen::Success)
    //        GISMO_ERROR("Could not solve the QR system for the specific b");

    // WORKAROUND
    if(!m_optimizationMatrix)
        optimize(gsWeightMapper<T>::optTargetToSource);
    // solve system with least squares method
    typename gsSparseSolver<T>::QR solver;
    solver.compute(*m_optimizationMatrix);
    if(solver.info()!=Eigen::Success)
        GISMO_ERROR("Could not compute the solver SparseQR");
    targetCoefs=solver.solve(sourceCoefs);
    if(solver.info()!=Eigen::Success)
        GISMO_ERROR("Could not solve the QR system for the specific b");
}

template<class T>
void gsWeightMapper<T>::sourceToTarget(indexType source, IndexContainer & target, WeightContainer & weights) const
{
    target.clear();
    weights.clear();
    // todo: add reserve
    for(InIterMat it = InIterMat(m_matrix,source);it;++it)
    {
        target.push_back(it.col());
        weights.push_back(it.value());
    }
}

template<class T>
void gsWeightMapper<T>::targetToSource(indexType target,IndexContainer & indizes,WeightContainer & weights) const
{
    indizes.clear();
    weights.clear();
    // todo: add reserve
    for(indexType i=0;i<m_matrix.rows();i++)
    {
        const weightType w = m_matrix.coeff(i,target);
        if(w!=0)
        {
            indizes.push_back(i);
            weights.push_back(w);
        }
    }
}

template<class T>
void gsWeightMapper<T>::fastSourceToTarget(IndexContainer const & source,IndexContainer & target) const
{
    GISMO_ASSERT(m_matrix.isCompressed(),"optimize() must be called on the mapper with fastSourceToTarget flag before using this function.");
    getNonZeroInnerOfMatrix<T>(m_matrix,source,target);
}

template<class T>
void gsWeightMapper<T>::fastTargetToSource(const IndexContainer &target, IndexContainer &source) const
{
    GISMO_ASSERT(m_optimizationMatrix,"optimize() must be called on the mapper with fastTargetToSource flag before using this function.");
    getNonZeroInnerOfMatrix<T>(*m_optimizationMatrix,target,source);
}


template<class T>
void gsWeightMapper<T>::getLocalMap (IndexContainer const & source, IndexContainer const & target, gsMatrix<T> &map) const
{
    GISMO_ASSERT(m_matrix.isCompressed(),"optimize() must be called on the mapper with fastSourceToTarget flag before using this function.");

    const index_t numRow=source.size();
    const index_t numCol=target.size();

    map.resize(numRow,numCol);

    for (index_t r=0;r<numRow;++r)
    {
        index_t c=0;
        Iterator coef = fastSourceToTarget(source[r]);

        while (c<numCol && coef)
        {
            if (target[c]< coef.index())
            {
                map(r,c)=0;
                ++c;
                continue;
            }
            if (target[c]== coef.index())
            {
                map(r,c)=coef.weight();
                ++c;
            }
            ++coef;
        }
        while(c<numCol)
        {
            map(r,c)=0;
            ++c;
        }
    }
}

template<class T>
void gsWeightMapper<T>::getLocalMap (IndexContainer const & source, gsMatrix<T> &map) const
{
    GISMO_ASSERT(m_matrix.isCompressed(),"optimize() must be called on the mapper with fastSourceToTarget flag before using this function.");

    const index_t numRow=source.size();
    const index_t numCol=getNrOfTargets();

    map.resize(numRow,numCol);

    for (index_t r=0;r<numRow;++r)
    {
        index_t c=0;
        Iterator coef = fastSourceToTarget(source[r]);

        while (c<numCol && coef)
        {
            if (c< coef.index())
            {
                map(r,c)=0;
                ++c;
                continue;
            }
            if (c== coef.index())
            {
                map(r,c)=coef.weight();
                ++c;
            }
            ++coef;
        }
        while(c<numCol)
        {
            map(r,c)=0;
            ++c;
        }
    }
}


}
