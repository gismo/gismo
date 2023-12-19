/** @file gsMapperUtils.cpp

    @brief implementation of the utility functions declared in gsMapperUtils.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include <gsMSplines/gsWeightMapperUtils.h>
#include <vector>
#include <algorithm>

namespace gismo {

template <class T>
index_t reorderMapperTarget(gsWeightMapper<T> &mapper, const std::vector<index_t>& perm, gsPermutationMatrix* permMatrix)
{
    const bool matSupplied = (permMatrix != NULL);
    const index_t size=mapper.getNrOfTargets();
    const index_t start=size-perm.size();
    const index_t flag=size+1; // invalid index

    std::vector<index_t> reordered (size,flag);
    for (size_t i=0; i<perm.size(); ++i)
    {
        reordered[perm[i]]=start+i;
    }
    index_t pos=0;
    for (size_t i=0; i<reordered.size(); ++i)
    {
        if(reordered[i]==flag)
        {
            reordered[i]=pos;
            ++pos;
        }
    }
    gsAsVector<index_t> permVec(reordered);
    gsPermutationMatrix m_perm(permVec);
    if(matSupplied)
        *permMatrix=m_perm;
    mapper*=m_perm.inverse();
    return start;
}


// assumes that the indices are not preset in the destination
// assumes row mayor
template <class T>
void copyToBlock (const gsWeightMapper<T> &src, typename gsWeightMapper<T>::LToGMatrix &dst, index_t rowShift, index_t colShift)
{
    src.optimize();
    for (index_t r=0; r<src.getNrOfSources();  ++r)
        for (typename gsWeightMapper<T>::Iterator entry=src.fastSourceToTarget(r); entry;  ++entry)
            dst.insert(r+rowShift,entry.index()+colShift)=entry.weight();
}


template <class T>
gsWeightMapper<T>* combineMappers (const std::vector<gsWeightMapper<T>*> &mappers, std::vector<index_t> &shifts, bool needShifting)
{
    gsWeightMapper<T>* result=new gsWeightMapper<T>();
    combineMappers(mappers,needShifting,shifts,*result);
    return result;
}

template <class T>
void combineMappers (const std::vector<gsWeightMapper<T>*> &mappers, bool needShifting, std::vector<index_t> &shifts, gsWeightMapper<T>& result)
{
    typename gsWeightMapper<T>::LToGMatrix &resultMatrix=result.asMatrix();

    shifts.clear();
    shifts.resize(mappers.size()+1);
    shifts[0]=0;

    index_t totTarSize=0;
    index_t totNNZ=0;
    index_t totShift=0;

    for (size_t m=0;m<mappers.size();++m)
    {
        totShift+=mappers[m]->getNrOfSources();
        shifts[m+1]=totShift;
        if (needShifting)
            totTarSize+=mappers[m]->getNrOfTargets();
        else
            totTarSize=std::max(totTarSize,mappers[m]->getNrOfTargets());
        totNNZ+=mappers[m]->asMatrix().nonZeros();
    }

    resultMatrix.resize(shifts.back(),totTarSize);
    resultMatrix.reserve(totNNZ);

    index_t startingGlobal=0;
    for (size_t m=0;m<mappers.size();++m)
    {
        copyToBlock( mappers[m]->asMatrix(), resultMatrix, shifts[m], startingGlobal);
        if (needShifting)
           startingGlobal+= mappers[m]->getNrOfTargets();
    }
}





}


