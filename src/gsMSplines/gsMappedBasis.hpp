/** @file gsMappedBasis.hpp

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsMSplines/gsMappedBasis.h>
#include <gsIO/gsFileData.h>

namespace gismo
{

template<short_t d,class T>
gsMappedBasis<d,T>::gsMappedBasis( const gsMappedBasis& other )
{
    m_topol = other.m_topol;
    // clone all geometries
    for(ConstBasisIter it = other.m_bases.begin();it!=other.m_bases.end();++it)
    {
        m_bases.push_back( (BasisType*)(*it)->clone().release() );
    }
    m_mapper=new gsWeightMapper<T>(*other.m_mapper);

    //m_sb = other.m_sb; //no: other.m_sb refers to other
}

template<short_t d,class T>
gsMappedBasis<d,T>::~gsMappedBasis()
{
    freeAll(m_bases);
    delete m_mapper;
}

template<short_t d,class T>
const std::vector<gsBasis<T>*> gsMappedBasis<d,T>::getBases() const
{
    std::vector<gsBasis<T>*> result;
    for (size_t i=0; i<m_bases.size();++i)
        result.push_back(m_bases[i]);
    return result;
}

template<short_t d,class T>
index_t gsMappedBasis<d,T>::size(const index_t index) const
{
    if(index == static_cast<int>(nPatches())-1)
        return size()-(size()/nPatches())*(nPatches()-1);
    else
        return  size()/nPatches();
}

template<short_t d,class T>
short_t gsMappedBasis<d,T>::maxDegree() const
{
    index_t deg = degree(0,0);
    for(size_t i=0;i<nPatches();i++)
        for(short_t j=0;j<m_bases[i]->dim();++j)
            if(degree(i,j)>deg)
                deg=degree(i,j);
    return deg;
}

template<short_t d,class T>
void gsMappedBasis<d,T>::addLocalIndicesOfPatchSide(const patchSide& ps,index_t offset,std::vector<index_t>& locals) const
{
    index_t patch = ps.patch;
    index_t side  = ps.side();
    index_t localOffset = _getFirstLocalIndex(patch);
    gsMatrix<index_t> indizes;
    //for(index_t i = 0;i<=offset;++i)
    {
        indizes=m_bases[patch]->boundaryOffset(side,offset);
        for(index_t j=0;j<indizes.rows();++j)
            locals.push_back(indizes.at(j)+localOffset);
    }
}

// HV[05/10/2023]: check for correctness
template<short_t d,class T>
void gsMappedBasis<d,T>::boundary(std::vector<index_t> & indices,index_t offset) const
{
    std::vector<index_t> locals;
    locals.reserve(this->size());
    typedef std::vector< patchSide >::const_iterator b_const_iter;
    for(b_const_iter iter = m_topol.bBegin();iter!=m_topol.bEnd();++iter)
        addLocalIndicesOfPatchSide(*iter,offset,locals);
    sort( locals.begin(), locals.end() );
    locals.erase( unique( locals.begin(), locals.end() ), locals.end() );
    m_mapper->sourceToTarget(locals,indices);
}

// HV[05/10/2023]: check for correctness
template<short_t d,class T>
void gsMappedBasis<d,T>::innerBoundaries(std::vector<index_t> & indices,index_t offset) const
{
    std::vector<index_t> locals;
    locals.reserve(this->size());
    typedef std::vector< gismo::boundaryInterface >::const_iterator i_const_iter;
    for(i_const_iter iter = m_topol.iBegin();iter!=m_topol.iEnd();++iter)
    {
        patchSide firstPs = iter->first();
        addLocalIndicesOfPatchSide(firstPs,offset,locals);
        patchSide secondPs = iter->second();
        addLocalIndicesOfPatchSide(secondPs,offset,locals);
    }
    sort( locals.begin(), locals.end() );
    locals.erase( unique( locals.begin(), locals.end() ), locals.end() );
    m_mapper->sourceToTarget(locals,indices);
}

template<short_t d,class T>
gsGeometry<T>* gsMappedBasis<d,T>::exportPatch(const index_t i,gsMatrix<T> const & localCoef) const
{
    const short_t geoDim=localCoef.cols();
    const index_t start = _getFirstLocalIndex(i);
    const index_t end   = _getLastLocalIndex(i);
    gsMatrix<T> coefs = localCoef.block(start,0,end-start+1,geoDim);
    return getBase(i).makeGeometry( give(coefs) ).release();
}

template<short_t d,class T>
gsMultiPatch<T> gsMappedBasis<d,T>::exportToPatches(gsMatrix<T> const & localCoef) const
{
    std::vector<gsGeometry<T> *> patches(nPatches());
    for(size_t i = 0; i<nPatches() ; ++i)
        patches[i]= exportPatch(i,localCoef);
    return gsMultiPatch<T>(patches,m_topol.boundaries(),m_topol.interfaces());
}

template<short_t d,class T>
void gsMappedBasis<d,T>::active_into(const index_t patch, const gsMatrix<T> & u,
                 gsMatrix<index_t>& result) const //global BF active on patch at point
{
    const index_t start = _getFirstLocalIndex(patch);
    gsMatrix<index_t> pActive;
    m_bases[patch]->active_into(u, pActive);
    // Shift actives by the offset of the current patch
    pActive.array() += start;
    const index_t numact  = pActive.rows();
    IndexContainer temp;
    result.resizeLike( pActive );
    std::vector<std::vector<index_t> > temp_output;//collects the outputs
    temp_output.resize( pActive.cols() );
    size_t max = 0;
    for(index_t i = 0; i< pActive.cols();i++)
    {
        std::vector<index_t> act_i( pActive.col(i).data(), pActive.col(i).data() + numact); // to be removed
        m_mapper->sourceToTarget(act_i, temp);
        temp_output[i]=temp;
        if(temp.size()>max)
            max=temp.size();
                //result.col(i) = gsAsConstMatrix<index_t>(temp, temp.size(), 1 ).cast<index_t>();
    }
    result.resize(max,u.cols());
    for(index_t i = 0; i < result.cols(); i++)
        for (index_t j = 0; j < result.rows();j++)
            if (size_t(j) < temp_output[i].size())
                result(j,i) = temp_output[i][j];
            else
                result(j,i) = 0 ;
}

template<short_t d,class T>
void gsMappedBasis<d,T>::eval_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    gsMatrix<index_t> bact;
    std::vector<index_t>  act, act0;
    gsMatrix<T> beval, map;//r:B,c:C
    const index_t shift=_getFirstLocalIndex(patch);

    gsVector<index_t> numAct;
    std::vector<gsMatrix<T>> result_tmp;
    result_tmp.resize(u.cols());
    numAct.resize(u.cols());
    for (index_t i = 0; i!=u.cols(); ++i)
    {
        m_bases[patch]->active_into(u.col(i), bact);
        act0 = std::vector<index_t>(bact.data(), bact.data()+bact.rows());
        m_bases[patch]->eval_into(u.col(i), beval);
        std::transform(act0.begin(), act0.end(), act0.begin(),
                       GS_BIND2ND(std::plus<index_t>(), shift));

        m_mapper->fastSourceToTarget(act0,act);
        m_mapper->getLocalMap(act0, act, map);
        result_tmp[i] = map.transpose() * beval; // todo: remove transpose()
        numAct[i] = result_tmp[i].rows();
    }

    result.setZero(numAct.maxCoeff(), u.cols());
    for (index_t i = 0; i!=u.cols(); ++i)
        for(index_t j = 0; j != result_tmp[i].rows(); j++) // result_tmp[i] == dim(rows,1)
            result(j,i) = result_tmp[i](j,0);

}

template<short_t d,class T>
void gsMappedBasis<d,T>::deriv_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    gsMatrix<index_t> actives;
    gsMatrix<T> curr_res;
    active_into(patch,u,actives);
    index_t rows = actives.rows();
    index_t cols = actives.cols();
    result.setZero(d*rows,cols);
    for(index_t j = 0; j<cols; j++) // For all points u.col(j)
    {
        const gsMatrix<T> & curr_u =  u.col(j);
        for(index_t i = 0; i<rows; i++) // For all actives at the point
        {
            if(actives(i,j)==0&&i>0)
                continue;
            derivSingle_into(patch, actives(i,j), curr_u, curr_res); // Evaluate N_j(u_i)
            result.block(i*d,j,d,1).noalias() = curr_res;
        }
    }
}

template<short_t d,class T>
void gsMappedBasis<d,T>::deriv2_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    gsMatrix<index_t> bact;
    m_bases[patch]->active_into(u, bact);
    std::vector<index_t>  act, act0(bact.data(), bact.data()+bact.rows());
    const index_t shift=_getFirstLocalIndex(patch);
    std::transform(act0.begin(), act0.end(), act0.begin(),
                   GS_BIND2ND(std::plus<index_t>(), shift));
    m_mapper->fastSourceToTarget(act0,act);
    gsMatrix<T> map;//r:B,c:C
    m_mapper->getLocalMap(act0, act, map);

    m_bases[patch]->deriv2_into(u, result);
    index_t       nr = bact.rows();
    const index_t nc = bact.cols();

    static const index_t sd = d*(d+1)/2;

    const index_t mr = map.cols();
    gsMatrix<T> tmp(sd*mr, nc);

    for (unsigned i = 0; i!=sd; ++i)
    {
        gsEigen::Map<typename gsMatrix<T>::Base, 0, gsEigen::Stride<-1,sd> >
            s(result.data()+i, nr, nc, gsEigen::Stride<-1,sd>(sd*nr,sd) );
        gsEigen::Map<typename gsMatrix<T>::Base, 0, gsEigen::Stride<-1,sd> >
            t(tmp.data()+i, mr, nc, gsEigen::Stride<-1,sd>(sd*mr,sd) );
        t = map.transpose() * s; //.noalias() bug
    }

    result.swap(tmp);
}

template<short_t d,class T>
void gsMappedBasis<d,T>::evalSingle_into(const index_t patch, const index_t global_BF, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    BasisType * this_patch = m_bases[patch];
    index_t start = _getFirstLocalIndex(patch), end = _getLastLocalIndex(patch);
    result.setZero(1,u.cols());
    if ( m_mapper->targetIsId(global_BF) )
    {
        IndexContainer indlist;
        m_mapper->targetToSource(global_BF, indlist);
        if(start <= indlist[0] && indlist[0] <= end)
            this_patch->evalSingle_into(_getPatchIndex(indlist[0]), u, result);
        else
            result.setZero(1,u.cols());
    }
    else
    {
        gsMatrix<T> allLocals;
        gsMatrix<index_t> pActive;
        for (index_t k = 0; k != u.cols(); k++)
        {
            this_patch->eval_into(u.col(k), allLocals);
            this_patch->active_into(u.col(k), pActive);
            gsSparseMatrix<T> L(allLocals.cols(),localSize()), Coefs(size(),1);
            const index_t offset = _getFirstLocalIndex(patch);
            for(index_t p = 0;p<allLocals.cols();++p)
                for(index_t j=0;j<pActive.rows();++j)
                    L(p,pActive(j,p)+offset)=allLocals(j,p);
            Coefs(global_BF,0)=1;
            gsSparseMatrix<T> temp = L*(m_mapper->asMatrix())*Coefs;
            result.col(k) = temp.transpose().toDense();
        }
    }
}

template<short_t d,class T>
void gsMappedBasis<d,T>::derivSingle_into(const index_t patch, const index_t global_BF, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    BasisType * this_patch = m_bases[patch];
    index_t start = _getFirstLocalIndex(patch), end = _getLastLocalIndex(patch);
    if ( m_mapper->targetIsId(global_BF) )
    {
        IndexContainer indlist;
        m_mapper->targetToSource(global_BF, indlist);
        if(start <= indlist[0] && indlist[0] <= end)
            this_patch->derivSingle_into(_getPatchIndex(indlist[0]), u, result);
        else
            result.setZero(d,u.cols());
    }
    else
    {
        gsMatrix<T> allLocals;
        gsMatrix<index_t> pActive;
        this_patch->deriv_into(u, allLocals);
        this_patch->active_into(u, pActive);
        gsSparseMatrix<T> L(2*allLocals.cols(),localSize()),Coefs(size(),1);
        const index_t offset = _getFirstLocalIndex(patch);
        for(index_t p = 0;p<allLocals.cols();++p)
            for(index_t j=0;j<pActive.rows();++j)
            {
                L(2*p,pActive(j,p)+offset)=allLocals(2*j,p);
                L(2*p+1,pActive(j,p)+offset)=allLocals(2*j+1,p); 
            }
        Coefs(global_BF,0)=1;
        gsSparseMatrix<T> temp = L*(m_mapper->asMatrix())*Coefs;
        result = temp.toDense();
    }
}

template<short_t d,class T>
void gsMappedBasis<d,T>::deriv2Single_into(const index_t patch, const index_t global_BF, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    const index_t blocksize = d*(d + 1)/2;
    BasisType * this_patch = m_bases[patch];
    index_t start = _getFirstLocalIndex(patch), end = _getLastLocalIndex(patch);
    if ( m_mapper->targetIsId(global_BF) )
    {
        IndexContainer indlist;
        m_mapper->targetToSource(global_BF, indlist);
        if(start <= indlist[0] && indlist[0] <= end)
            this_patch->deriv2Single_into(_getPatchIndex(indlist[0]), u, result);
        else
            result.setZero(blocksize,u.cols());
    }
    else
    {
        gsMatrix<T> allLocals;
        gsMatrix<index_t> pActive;
        this_patch->deriv2_into(u, allLocals);
        this_patch->active_into(u, pActive);
        gsSparseMatrix<T> L(3*allLocals.cols(),localSize()),Coefs(size(),1);
        const index_t offset = _getFirstLocalIndex(patch);
        for(index_t p = 0;p<allLocals.cols();++p)
            for(index_t j=0;j<pActive.rows();++j)
            {
                L(3*p,pActive(j,p)+offset)=allLocals(3*j,p);
                L(3*p+1,pActive(j,p)+offset)=allLocals(3*j+1,p);
                L(3*p+2,pActive(j,p)+offset)=allLocals(3*j+2,p);
            }
        Coefs(global_BF,0)=1;
        gsSparseMatrix<T> temp = L*(m_mapper->asMatrix())*Coefs;
        result = temp.toDense();
    }
}

template<short_t d,class T>
void gsMappedBasis<d,T>::evalAllDers_into(const index_t patch, const gsMatrix<T> & u,
                                             const index_t n, std::vector<gsMatrix<T> >& result) const
{
    gsMatrix<index_t> bact;
    m_bases[patch]->active_into(u, bact);
    std::vector<index_t>  act, act0(bact.data(), bact.data()+bact.rows());
    const index_t shift=_getFirstLocalIndex(patch);
    std::transform(act0.begin(), act0.end(), act0.begin(),
                   GS_BIND2ND(std::plus<index_t>(), shift));
    m_mapper->fastSourceToTarget(act0,act);
    gsMatrix<T> map;//r:B,c:C
    m_mapper->getLocalMap(act0, act, map);

    m_bases[patch]->evalAllDers_into(u, n, result);
    index_t       nr = result.front().rows();
    const index_t nc = result.front().cols();
    result.front() = map.transpose() * result.front();

    if ( n>0 )
    {
        const index_t mr = map.cols();
        gsMatrix<T> tmp(d*mr, nc);

        for (index_t i = 0; i!=d; ++i)
        {
            gsEigen::Map<typename gsMatrix<T>::Base, 0, gsEigen::Stride<-1,d> >
                s(result[1].data()+i, nr, nc, gsEigen::Stride<-1,d>(d*nr,d) );
            gsEigen::Map<typename gsMatrix<T>::Base, 0, gsEigen::Stride<-1,d> >
                t(tmp.data()+i, mr, nc, gsEigen::Stride<-1,d>(d*mr,d) );
            t = map.transpose() * s; //.noalias() bug
        }
        result[1].swap(tmp);

        if ( n>1 )
        {
            static const index_t sd = d*(d+1)/2;
            tmp.resize(sd*mr, nc);

            for (index_t i = 0; i!=sd; ++i)
            {
                gsEigen::Map<typename gsMatrix<T>::Base, 0, gsEigen::Stride<-1,sd> >
                    s(result[2].data()+i, nr, nc, gsEigen::Stride<-1,sd>(sd*nr,sd) );
                gsEigen::Map<typename gsMatrix<T>::Base, 0, gsEigen::Stride<-1,sd> >
                    t(tmp.data()+i, mr, nc, gsEigen::Stride<-1,sd>(sd*mr,sd) );
                t = map.transpose() * s; //.noalias() bug
            }
            result[2].swap(tmp);
        }
    }
    GISMO_ASSERT( n < 3, "gsMappedBasis::evalAllDers() not implemented for 2<n." );
}

template<short_t d,class T>
void gsMappedBasis<d,T>::evalAllDersSingle_into(const index_t patch, const index_t global_BF, const gsMatrix<T> & u,const index_t n,gsMatrix<T> & result ) const
{
    GISMO_ASSERT( n<2, "gsTensorBasis::evalAllDers() not implemented for 1<n." );
    result.resize(( 2*n + 1 ), u.cols());
    BasisType * this_patch = m_bases[patch];
    result.setZero();
    IndexContainer allLocalIndizes,localIndizes;
    WeightContainer allWeights,weights;
    gsMatrix<T> res;
    m_mapper->targetToSource(global_BF,allLocalIndizes,allWeights);
    const index_t start = _getFirstLocalIndex(patch);
    const index_t end = _getLastLocalIndex(patch);
    for(size_t i = 0;i<allLocalIndizes.size();++i)
        if(allLocalIndizes[i]>=start && allLocalIndizes[i]<=end)
        {
            weights.push_back(allWeights[i]);
            localIndizes.push_back(allLocalIndizes[i]);
        }
    for(index_t i=0;i<=n;i++)
    {
        ConstWeightIter wIter = weights.begin();
        for(ConstIndexIter it=localIndizes.begin();it!=localIndizes.end();++it)
        {
            if(i==0)
                this_patch->evalSingle_into(_getPatchIndex(*it),u,res);
            else if(i==1)
                this_patch->derivSingle_into(_getPatchIndex(*it),u,res);
            else
                GISMO_ERROR("only for n<2");
            for(index_t j=0;j<u.cols();j++)
            {
                result(i,j)+=res(0,j)*(*wIter);
                if(i==1)
                    result(i+1,j)+=res(1,j)*(*wIter);
            }
            ++wIter;
        }
    }
}

template<short_t d,class T>
std::ostream & gsMappedBasis<d,T>::print(std::ostream &os) const
{
    os << "gsMappedBasis:\n";
    os << "\t Local size: "<<this->localSize()<<"\n";
    os << "\t Global size: "<<this->globalSize()<<"\n";
    return os;
}

template<short_t d,class T>
index_t gsMappedBasis<d,T>::_getPatch(const index_t localIndex) const
{
    size_t patch;
    index_t patchIndex=localIndex;
    for(patch=0;patch<m_bases.size();patch++)
    {
        if(patchIndex>=static_cast<index_t>(m_bases[patch]->size()))
            patchIndex-=m_bases[patch]->size();
        else
            break;
    }
    return patch;
}

template<short_t d,class T>
index_t gsMappedBasis<d,T>::_getPatchIndex(const index_t localIndex) const
{
    index_t patchIndex=localIndex;
    for(index_t i = 0;i<_getPatch(localIndex);i++)
        patchIndex-=m_bases[i]->size();
    return patchIndex;
}

template<short_t d,class T>
index_t gsMappedBasis<d,T>::_getFirstLocalIndex(index_t const patch) const
{
    index_t index=0;
    for(index_t i=0;i<patch;i++)
        index+=m_bases[i]->size();
    return index;
}

} // namespace gismo
