/** @file gsCompositeBSplineBasis.hpp

    @brief Provides implementation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsMSplines/gsCompositeBSplineBasis.h>
#include <gsMSplines/gsCompositeIncrSmoothnessBasis.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

namespace gismo
{

template<short_t d, class T>
gsCompositeBSplineBasis<d,T>::gsCompositeBSplineBasis (const BasisContainer & bases, gsBoxTopology const & topol,
                                                       int increaseSmoothnessLevel, int minEVDistance)
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

    repairPatches();

    if(!_checkTopologyWithBases())
        GISMO_ERROR("topology and bases not suitable for composite basis");
    _setMapping();
}

template<short_t d, class T>
gsCompositeBSplineBasis<d,T>::gsCompositeBSplineBasis (const BasisContainer & bases, gsBoxTopology const & topol,std::vector<gsMatrix<T> * > coefs,
                                                       int increaseSmoothnessLevel, int minEVDistance)
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
gsCompositeBSplineBasis<d,T>::gsCompositeBSplineBasis( gsMultiPatch<T> const & mp, int increaseSmoothnessLevel,
                                                       int minEVDistance)
{
    //topol.computeAllVertices();
    for (size_t i = 0; i < mp.nPatches(); i++)
    {
        GISMO_ASSERT(dynamic_cast<gsTensorBSplineBasis<2> * >(& mp.basis(i))!=NULL,"Bases is not of type gsTensorBSplineBasis<2>");
    }
    m_topol = mp;
    m_bases = mp.basesCopy();
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
gsCompositeBSplineBasis<d,T>::gsCompositeBSplineBasis( const gsCompositeBSplineBasis& other )
{
    m_topol = other.m_topol;
    m_vertices = other.m_vertices;
    m_distances=other.m_distances;
    m_minDist=other.m_minDist;
    // clone all geometries
    for(ConstBasisIter it = other.m_bases.begin();it!=other.m_bases.end();++it)
    {
        m_bases.push_back( (gsBasis<T>*)(*it)->clone().release() );
    }
    m_incrSmoothnessDegree=other.m_incrSmoothnessDegree;
    if(!this->_checkTopologyWithBases())
        GISMO_ERROR("topology and bases not suitable for composite basis");
    //_setMapping();
    m_mapper=new gsWeightMapper<T>(*other.m_mapper);
    //m_mapFactory= need a method to clone mapFactory
}

template<short_t d, class T>
gsCompositeBSplineBasis<d,T> &gsCompositeBSplineBasis<d,T>::operator=( const gsCompositeBSplineBasis& other )
{
    m_topol=other.m_topol;
    m_vertices=other.m_vertices;
    m_distances=other.m_distances;
    m_minDist=other.m_minDist;
    freeAll(m_bases);
    m_bases.clear();
    for(ConstBasisIter it = other.m_bases.begin();it!=other.m_bases.end();++it)
        m_bases.push_back( (gsBasis<T>*)(*it)->clone().release() );
    if(m_mapper)
        delete m_mapper;
    m_mapper=NULL;//other.m_mapper->clone();
    m_incrSmoothnessDegree=other.m_incrSmoothnessDegree;
    //m_mapFactory= need a method to clone mapFactory
    return *this;
}

template<short_t d, class T>
unsigned gsCompositeBSplineBasis<d,T>::basisFunctionsOnSide(const patchSide& ps) const
{
    if(ps.side()==1||ps.side()==2)
        return basis(ps.patch).size(1);
    else
        return basis(ps.patch).size(0);
}

template<short_t d, class T>
bool gsCompositeBSplineBasis<d,T>::isLocallyConnected(indexType i,indexType j) const
{
    unsigned patch_i = _getPatch(i), patch_j = _getPatch(j);
    if( patch_i != patch_j )
        return false;
    indexType loc_i = _getPatchIndex(i), loc_j = _getPatchIndex(j);
    gsVector<index_t,d> vec_i = basis(patch_i).tensorIndex(loc_i);
    gsVector<index_t,d> vec_j = basis(patch_j).tensorIndex(loc_j);
    unsigned distance = 0;
    for( unsigned i2=0; i2< d;++i2 )
    {
        unsigned dist = std::abs(static_cast<int>(vec_i[i2])-static_cast<int>(vec_j[i2]));
        if( dist != 0 )
            distance+=dist;
    }
    return distance==1;
}

template<short_t d, class T>
void gsCompositeBSplineBasis<d,T>::numActive_into(const index_t patch,const gsMatrix<T> & u, gsVector<index_t>& result) const
{
    result.resize(u.cols());
    for(int i = 0; i<u.cols();++i)
    {
        result(i)=static_cast<unsigned>((degree(patch,0)+1)*(degree(patch,1)+1));
    }
}

template<short_t d, class T>
void gsCompositeBSplineBasis<d,T>::refine(const index_t patch,const std::vector<T>& knots_u, const std::vector<T>& knots_v,bool updateBasis)
{
    std::vector<std::vector<T> >refineKnots(2);
    refineKnots[0]=knots_u;
    refineKnots[1]=knots_v;
    basis(patch).insertKnots(refineKnots);
    if(updateBasis)
    {
        repairPatches();
        updateTopol();
    }
}

template<short_t d, class T>
void gsCompositeBSplineBasis<d,T>::refine(const index_t patch, gsMatrix<T> const & boxes, bool updateBasis)
{
    std::vector<std::vector<T> > refineKnots;
    _boxToRefineKnots(patch,boxes,refineKnots);
    refine(patch,refineKnots[0],refineKnots[1],updateBasis);
}

template<short_t d, class T>
void gsCompositeBSplineBasis<d,T>::refineElements(const index_t patch, std::vector<index_t> const & boxes, bool updateBasis)
{
    gsMatrix<T> mat_boxes;
    _boxesVectorToMatrix(boxes,mat_boxes);
    refine(patch,mat_boxes,updateBasis);
}

template<short_t d, class T>
void gsCompositeBSplineBasis<d,T>::refine_withCoefs(gsMatrix<T>& localCoef, const index_t patch,const std::vector<T>& knots_u, const std::vector<T>& knots_v,
                                                    bool updateBasis)
{
    std::vector<std::vector<T> >refineKnots(2);
    refineKnots[0]=knots_u;
    refineKnots[1]=knots_v;
    std::vector<gsMatrix<T> *> coefs;
    unsigned geoDim = localCoef.cols();
    int start, end = -1;
    for (size_t i = 0; i<nPatches(); ++i)
    {
        start=end+1;
        end+=basis(i).size();
        gsMatrix<T>* localMat = new gsMatrix<T>(end-start+1,geoDim);
        *localMat << localCoef.block(start,0,end-start+1,geoDim);
        coefs.push_back(localMat);
    }
    basis(patch).refine_withCoefs(*coefs[patch],refineKnots);
    if(updateBasis)
        repairPatches(coefs,patch);
    unsigned totalSize=0;
    for (size_t i = 0; i<nPatches(); ++i)
    {
        totalSize+=basis(i).size();
    }
    localCoef.resize(totalSize,geoDim);
    end = -1;
    for (size_t i = 0;i<nPatches();i++)
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
void gsCompositeBSplineBasis<d,T>::refine_withCoefs(gsMatrix<T>& coefs, const index_t patch, gsMatrix<T> const & boxes,
                                                    bool updateBasis)
{
    std::vector<std::vector<T> > refineKnots;
    _boxToRefineKnots(patch,boxes,refineKnots);
    refine_withCoefs(coefs,patch,refineKnots[0],refineKnots[1],updateBasis);
}

template<short_t d, class T>
void gsCompositeBSplineBasis<d,T>::refineElements_withCoefs(gsMatrix<T>& coefs, const index_t patch, std::vector<index_t> const & boxes,
                                                            bool updateBasis)
{
    gsMatrix<T> mat_boxes;
    _boxesVectorToMatrix(boxes,mat_boxes);
    refine_withCoefs(coefs,patch,mat_boxes,updateBasis);
}

template<short_t d, class T>
void gsCompositeBSplineBasis<d,T>::_boxesVectorToMatrix(const std::vector<index_t> & boxes, gsMatrix<T> & mat_boxes)
{
    GISMO_ASSERT( boxes.size() % (2 * d + 1 ) == 0, "The points did not define boxes properly. The boxes were not added to the basis.");
    mat_boxes.resize(2,2*(boxes.size()/(2*d + 1)));
    for( unsigned int i = 0; i < (boxes.size())/(2*d+1); i++)
    {
        mat_boxes(0,2*i)=boxes[(i*(2*d+1))+1];
        mat_boxes(1,2*i)=boxes[(i*(2*d+1))+2];
        mat_boxes(0,2*i+1)=boxes[(i*(2*d+1))+3];
        mat_boxes(1,2*i+1)=boxes[(i*(2*d+1))+4];
    }
}

template<short_t d, class T>
void gsCompositeBSplineBasis<d,T>::_boxToRefineKnots(const index_t patch,gsMatrix<T> const & boxes,std::vector<std::vector<T> > & refineKnots)
{
    gsTensorBSplineBasis<d,T> & this_base=basis(patch);
    GISMO_ASSERT( boxes.rows() == this_base.dim() , "Number of rows of refinement boxes must equal dimension of parameter space.");
    GISMO_ASSERT( boxes.cols() % 2 == 0, "Refinement boxes must have even number of columns.");

    const T tol = 0.000000001;
    gsKnotVector<T> knots_u = this_base.component(0).knots();
    gsKnotVector<T> knots_v = this_base.component(1).knots();
    std::vector<T> ins_knots_u;
    std::vector<T> ins_knots_v;
    for(size_t i = 1;i < knots_u.size();i++)
        if( knots_u[i]-knots_u[i-1] > tol)
        {
            T midpt = knots_u[i-1] + (knots_u[i]-knots_u[i-1])/2;
            for( index_t j=0; j < boxes.cols(); j+=2 ) // loop over all boxes
            {
                if( boxes(0,j) < midpt && midpt < boxes(0,j+1) )
                    ins_knots_u.push_back(midpt);
            }
        }
    for(size_t i = 1;i < knots_v.size();i++)
        if( knots_v[i]-knots_v[i-1] > tol)
        {
            T midpt = knots_v[i-1] + (knots_v[i]-knots_v[i-1])/2;
            for( index_t j=0; j < boxes.cols(); j+=2 ) // loop over all boxes
            {
                if( boxes(1,j) < midpt && midpt < boxes(1,j+1) )
                    ins_knots_v.push_back(midpt);
            }
        }
    refineKnots.resize(2);
    refineKnots[0]=ins_knots_u;
    refineKnots[1]=ins_knots_v;
}

template<short_t d, class T>
unsigned gsCompositeBSplineBasis<d,T>::getNrOfSpecialKnots(const gsKnotVector<T> kv,const std::vector<T>& new_knots,bool par,int distance)
{
    int specialKnots=0, index = 0, kv_index, nk_index;
    for(int i=0;i<=distance+kv.degree();++i)
    {
        kv_index = par ? kv.size()-i-1 : i;
        nk_index = par ? new_knots.size()-index-1 : index;
        if(nk_index>=static_cast<int>(new_knots.size())||nk_index<0)
            break;
        if((!par&&new_knots[nk_index]<=kv.at(kv_index))||
                (par&&new_knots[nk_index]>=kv.at(kv_index)))
        {
            specialKnots++;
            index++;
        }
    }
    return specialKnots;
}

template<short_t d, class T>
void gsCompositeBSplineBasis<d,T>::repairPatches(std::vector<gsMatrix<T> *> & coefs,
                                                 index_t startFromPatch)
{
    std::set<index_t> toCheck; //set of all indizes of patches, which have to be checked
    if(startFromPatch==-1)
        for(unsigned i = 0;i<m_bases.size();++i)
            toCheck.insert(i);
    else
        toCheck.insert(startFromPatch);
    unsigned curElement;
    do
    {
        curElement=*(toCheck.begin());
        toCheck.erase(curElement);
        std::vector<std::vector<T> >refineKnots;
        std::vector<index_t> checkPatches;
        bool matched = _knotsMatchNeighbours(curElement,refineKnots,checkPatches);
        if(!matched)
        {
            if(coefs[curElement]!=NULL)
                basis(curElement).refine_withCoefs(*(coefs[curElement]),refineKnots);
            else
                basis(curElement).insertKnots(refineKnots);
        }
        for(unsigned i = 0;i<checkPatches.size();++i)
        {
            toCheck.insert(checkPatches[i]);
        }
    }while(!toCheck.empty());
}

template<short_t d, class T>
bool gsCompositeBSplineBasis<d,T>::_knotsMatchNeighbours(index_t patch,std::vector<std::vector<T> >& knotsToInsert,
                                                         std::vector<index_t>& checkPatches,
                                                         T eps)
{
    typedef typename gsKnotVector<T>::iterator knotIter;
    bool matched = true;
    gsTensorBSplineBasis<d,T> & base = basis(patch);
    patchSide ps, ps_neigh;
    knotsToInsert.resize(2);
    checkPatches.clear();
    checkPatches.reserve(nPatches());
    std::vector<bool> sideToCheck(4,false);
    std::vector<int> neighbours(4,-1);
    gsKnotVector<T> knots_neigh;
    knotIter kv,kv_end,kvN,kvN_end;
    for(boxSide side=boxSide::getFirst(2);side<boxSide::getEnd(2);++side)
    {
        ps=patchSide(patch,side);
        if( !m_topol.getNeighbour(ps,ps_neigh) )
            continue;
        knots_neigh = basis(ps_neigh.patch).knots(1-(ps_neigh.direction()));
        boundaryInterface bf;
        m_topol.getInterface(ps,bf);
        if(!bf.dirOrientation(ps,1-ps.direction()))
            knots_neigh.reverse();
        kv = base.knots(1-side.direction()).begin();
        kv_end = base.knots(1-side.direction()).end();
        kvN = knots_neigh.begin();
        kvN_end = knots_neigh.end();
        neighbours[side-1]=ps_neigh.patch; //fill with the patch numbers of the neighbours
        while(kv != kv_end || kvN != kvN_end)
        {
            if(math::abs(*kv-*kvN)<eps)
            {
                ++kv;
                ++kvN;
            }
            else if(*kv < *kvN)
            {
                ++kv;
                sideToCheck[side-1]=true; //set this side to true
            }
            else
            {
                knotsToInsert[1-side.direction()].push_back(*kvN);
                matched = false;
                ++kvN;
                sideToCheck[side.opposite()-1]=true; //set the opposite side to true
            }
        }
    }
    //remove duplicates:
    sort( knotsToInsert[0].begin(), knotsToInsert[0].end() );
    knotsToInsert[0].erase( unique( knotsToInsert[0].begin(), knotsToInsert[0].end() ), knotsToInsert[0].end() );
    sort( knotsToInsert[1].begin(), knotsToInsert[1].end() );
    knotsToInsert[1].erase( unique( knotsToInsert[1].begin(), knotsToInsert[1].end() ), knotsToInsert[1].end() );
    //find the patches we will have to check again:
    for(unsigned i = 0;i<4;i++)
    {
        if(sideToCheck[i]&&neighbours[i]!=-1)
            checkPatches.push_back(neighbours[i]);
    }
    return matched;
}

template<short_t d, class T>
T gsCompositeBSplineBasis<d,T>::findParameter(patchSide const & ps,patchCorner const & pc,unsigned nrBasisFuncs) const
{
    if(nrBasisFuncs==0)
        return 0.0;
    gsKnotVector<T> knots = basis(ps.patch).knots(1-(ps.direction()));
    gsVector<bool> pars;
    pc.parameters_into(d,pars);
    typename gsKnotVector<T>::uiterator iter;
    if(!pars(1-ps.direction()))
        knots.reverse();
    iter = knots.ubegin();
    iter+=nrBasisFuncs;
    return *iter;
}

}
