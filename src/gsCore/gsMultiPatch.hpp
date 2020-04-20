/** @file gsMultiPatch.hpp

    @brief Provides declaration of the MultiPatch class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsCore/gsGeometry.h>
#include <gsCore/gsDofMapper.h>
#include <gsCore/gsAffineFunction.h>

#include <gsUtils/gsCombinatorics.h>

namespace gismo
{

template<class T>
gsMultiPatch<T>::gsMultiPatch(const gsGeometry<T> & geo )
    : BaseA( geo.parDim() )
{
    m_patches.push_back(geo.clone().release());
    //m_patches[0]->setId(0); // Note: for the single-patch constructor the id remains unchanged
    addBox();
    this->addAutoBoundaries();
}

template<class T>
gsMultiPatch<T>::gsMultiPatch( const gsMultiPatch& other )
    : BaseA( other ), BaseB( other ), m_patches( other.m_patches.size() )
{
    // clone all geometries
    cloneAll( other.m_patches.begin(), other.m_patches.end(),
              this->m_patches.begin());
}

#if EIGEN_HAS_RVALUE_REFERENCES

template<class T>
gsMultiPatch<T>& gsMultiPatch<T>::operator=( const gsMultiPatch& other )
{
    if (this!=&other)
    {
        freeAll(m_patches);
        BaseA::operator=(other);
        m_patches.resize(other.m_patches.size());
        // clone all geometries
        cloneAll( other.m_patches.begin(), other.m_patches.end(),
                this->m_patches.begin());
    }
    return *this;
}

template<class T>
gsMultiPatch<T>& gsMultiPatch<T>::operator=( gsMultiPatch&& other )
{
    freeAll(m_patches);
    BaseA::operator=(give(other));
    m_patches = give(other.m_patches);
    return *this;
}


#endif

template<class T>
gsMultiPatch<T>::gsMultiPatch(PatchContainer & patches )
    : BaseA( patches[0]->parDim(), patches.size() )
{
    m_patches.swap(patches); // patches are consumed
    setIds();
    this->addAutoBoundaries();
}

template<class T>
gsMultiPatch<T>::gsMultiPatch( PatchContainer& patches,
                               const std::vector<patchSide>& boundary,
                               const std::vector<boundaryInterface>& interfaces )
    : BaseA( patches[0]->parDim(), patches.size(), boundary, interfaces )
{
    m_patches.swap(patches); // patches are consumed
    setIds();
}

template<class T>
gsMultiPatch<T>::~gsMultiPatch()
{
    freeAll(m_patches);
}

template<class T>
void gsMultiPatch<T>::setIds()
{
    size_t id = 0;
    for ( iterator it = m_patches.begin(); it != m_patches.end(); ++it )
    {
        ( *it )->setId( id++ );
    }
}

template<class T>
std::ostream& gsMultiPatch<T>::print(std::ostream& os) const
{
    if ( !this->empty() ) {
        os << "gsMultiPatch (" << this->nPatches() << "): ";
        os << "#Boundaries= " << nBoundary() << ", ";
        os << "#Interfaces= " << nInterfaces() << ".\n";
    } else {
        os << "gsMultiPatch ( empty! ).\n";
    }
    return os;
}

template<class T>
std::string gsMultiPatch<T>::detail() const
{
    std::ostringstream os;
    print(os);
    if ( nPatches() > 0 )
    {
        BaseA::print( os );
    }
    return os.str();
}

template<class T>
short_t gsMultiPatch<T>::geoDim() const
{
    GISMO_ASSERT( m_patches.size() > 0 , "Empty multipatch object.");
    return m_patches[0]->geoDim();
}

template<class T>
short_t gsMultiPatch<T>::coDim() const
{
    GISMO_ASSERT( m_patches.size() > 0 , "Empty multipatch object.");
    return m_patches[0]->geoDim() - m_dim;
}

template<class T>
gsMatrix<T>
gsMultiPatch<T>::parameterRange(int i) const
{
    return m_patches[i]->basis().support();
}

template<class T>
gsBasis<T> &
gsMultiPatch<T>::basis( size_t i ) const
{
    GISMO_ASSERT( i < m_patches.size(), "Invalid patch index requested from gsMultiPatch" );
    return m_patches[i]->basis();
}

template<class T>
std::vector<gsBasis<T> *> gsMultiPatch<T>::basesCopy(bool NoRational) const
{
    std::vector<gsBasis<T> *> bb;
    bb.reserve(m_patches.size());

    if (NoRational)
    {
        for ( const_iterator it = m_patches.begin();
              it != m_patches.end(); ++it )
        {
            bb.push_back( (*it)->basis().source().clone().release() );
        }
    }
    else
    {
        for ( const_iterator it = m_patches.begin();
              it != m_patches.end(); ++it )
        {
            bb.push_back( (*it)->basis().clone().release() );
        }
    }
    return bb ;
}

template<class T>
void gsMultiPatch<T>::permute(const std::vector<short_t> & perm)
{
    gsAsVector<gsGeometry<T>*> a (m_patches);
    a = Eigen::PermutationMatrix<-1,-1,short_t>(gsAsConstVector<short_t>(perm)) * a;
}

template<class T>
void gsMultiPatch<T>::addPatch(typename gsGeometry<T>::uPtr g)
{
    if ( m_dim == -1 )
    {
        m_dim = g->parDim();
    } else
    {
        GISMO_ASSERT( m_dim == g->parDim(),
                      "Tried to add a patch of different dimension in a multipatch." );
    }
    g->setId( m_patches.size() );
    m_patches.push_back( g.release() ) ;
    addBox();
}

template<class T>
inline void gsMultiPatch<T>::addPatch(const gsGeometry<T> & g)
{
    addPatch(g.clone());
}

template<class T>
size_t gsMultiPatch<T>::findPatchIndex( gsGeometry<T>* g ) const
{
    const_iterator it
        = std::find( m_patches.begin(), m_patches.end(), g );
    GISMO_ASSERT( it != m_patches.end(), "Did not find the patch index." );
    // note: should return g->patchId();
    return it - m_patches.begin();
}

template<class T>
void gsMultiPatch<T>::addInterface( gsGeometry<T>* g1, boxSide s1,
                                    gsGeometry<T>* g2, boxSide s2 )
{
    int p1 = findPatchIndex( g1 );
    int p2 = findPatchIndex( g2 );
    BaseA::addInterface( p1, s1, p2, s2 );
}

template<class T>
gsMatrix<T> gsMultiPatch<T>::pointOn( const patchCorner& pc )
{
    gsMatrix<T> coordinates;
    m_patches[pc.patch]->eval_into(m_patches[pc.patch]->parameterCenter(pc),coordinates);
    return coordinates;
}


template<class T>
gsMatrix<T> gsMultiPatch<T>::pointOn( const patchSide& ps )
{
    gsMatrix<T> coordinates;
    m_patches[ps.patch]->eval_into(m_patches[ps.patch]->parameterCenter(ps),coordinates);
    return coordinates;
}

template<class T>
void gsMultiPatch<T>::uniformRefine(int numKnots, int mul)
{
    for ( typename PatchContainer::const_iterator it = m_patches.begin();
          it != m_patches.end(); ++it )
    {
        ( *it )->uniformRefine(numKnots, mul);
    }
}

template<class T>
void gsMultiPatch<T>::degreeElevate(int elevationSteps)
{
    for ( typename PatchContainer::const_iterator it = m_patches.begin();
          it != m_patches.end(); ++it )
    {
        ( *it )->degreeElevate(elevationSteps, -1);
    }
}


template<class T>
void gsMultiPatch<T>::boundingBox(gsMatrix<T> & result) const
{
    result.setZero(geoDim(),2);
    if ( m_patches.size() == 0 )
        return;

    // to do: this is the bounding box of a gsGeometry, make a member function
    result.col(0) = patch(0).coefs().colwise().minCoeff().transpose();
    result.col(1) = patch(0).coefs().colwise().maxCoeff().transpose();

    for (const_iterator it = begin()+1; it != end(); ++it)
    {
        const gsMatrix<T>  & cc = (*it)->coefs();
        result.col(0) = result.col(0).cwiseMin( cc.colwise().minCoeff().transpose() ) ;
        result.col(1) = result.col(1).cwiseMax( cc.colwise().maxCoeff().transpose() ) ;
    }
}

template<class T>
gsMultiPatch<T> gsMultiPatch<T>::uniformSplit() const
{
    int n = math::exp2(parDim());
    std::vector<gsGeometry<T>*> result;
    result.reserve(nPatches() * n);

    for (size_t np = 0; np < nPatches(); ++np)
    {
        std::vector<gsGeometry<T>*> result_temp = m_patches[np]->uniformSplit();
        result.insert(result.end(), result_temp.begin(), result_temp.end());
    }
    gsMultiPatch<T> mp(result);
    mp.computeTopology();
    return mp;
}


/*
  This is based on comparing a set of reference points of the patch
  side and thus it implicitly assumes that the patch faces match
*/
template<class T>
bool gsMultiPatch<T>::computeTopology( T tol, bool cornersOnly )
{
    BaseA::clearTopology();

    const size_t   np    = m_patches.size();
    const index_t  nCorP = 1 << m_dim;     // corners per patch
    const index_t  nCorS = 1 << (m_dim-1); // corners per side

    gsMatrix<T> supp,
    // Parametric coordinates of the reference points. These points
    // are used to decide if two sides match.
    // Currently these are the corner points and the side-centers
    coor;
    if (cornersOnly)
        coor.resize(m_dim,nCorP);
    else
        coor.resize(m_dim,nCorP + 2*m_dim);

    gsVector<bool> boxPar(m_dim);

    // each matrix contains the physical coordinates of the reference points
    std::vector<gsMatrix<T> > pCorners(np);

    std::vector<patchSide> pSide; // list of all candidate patchSides to compare
    pSide.reserve(np * 2 * m_dim);

    for (size_t p=0; p<np; ++p)
    {
        supp = m_patches[p]->parameterRange(); // the parameter domain of patch i

        // Corners' parametric coordinates
        for (boxCorner c=boxCorner::getFirst(m_dim); c<boxCorner::getEnd(m_dim); ++c)
        {
            boxPar   = c.parameters(m_dim);
            for (index_t i=0; i<m_dim;++i)
                coor(i,c-1) = boxPar(i) ? supp(i,1) : supp(i,0);
        }

        if (!cornersOnly)
        {
            // Sides' centers parametric coordinates
            index_t l = nCorP;
            for (boxSide c=boxSide::getFirst(m_dim); c<boxSide::getEnd(m_dim); ++c)
            {
                const index_t dir = c.direction();
                const index_t s   = static_cast<index_t>(c.parameter());// 0 or 1

                for (index_t i=0; i<m_dim;++i)
                    coor(i,l) = ( dir==i ?  supp(i,s) :
                                  (supp(i,1)+supp(i,0))/2.0 );
                l++;
            }
        }

        // Evaluate the patch on the reference points
        m_patches[p]->eval_into(coor,pCorners[p]);

        // Add all the patchSides of this patch to the candidate list
        for (boxSide bs=boxSide::getFirst(m_dim); bs<boxSide::getEnd(m_dim); ++bs)
            pSide.push_back(patchSide(p,bs));
    }

    gsVector<index_t>      dirMap(m_dim);
    gsVector<bool>         matched(nCorS), dirOr(m_dim);
    std::vector<boxCorner> cId1, cId2;
    cId1.reserve(nCorS);
    cId2.reserve(nCorS);

    while ( pSide.size() != 0 )
    {
        bool done = false;
        const patchSide side = pSide.back();
        pSide.pop_back();
        for (size_t other=0; other<pSide.size(); ++other)
        {
            side        .getContainedCorners(m_dim,cId1);
            pSide[other].getContainedCorners(m_dim,cId2);
            matched.setConstant(false);

            // Check whether the side center matches
            if (!cornersOnly)
                if ( ( pCorners[side.patch        ].col(nCorP+side-1        ) -
                       pCorners[pSide[other].patch].col(nCorP+pSide[other]-1)
                         ).norm() >= tol )
                    continue;

            // Check whether the vertices match and compute direction map and orientation
            if ( matchVerticesOnSide( pCorners[side.patch]        , cId1, 0,
                                      pCorners[pSide[other].patch], cId2,
                                      matched, dirMap, dirOr, tol ) )
            {
                dirMap(side.direction()) = pSide[other].direction();
                dirOr (side.direction()) = !( side.parameter() == pSide[other].parameter() );
                BaseA::addInterface( boundaryInterface(side, pSide[other], dirMap, dirOr));
                // done with pSide[other], remove it from candidate list
                std::swap( pSide[other], pSide.back() );
                pSide.pop_back();
                done=true;
                break;//for (size_t other=0..)
            }
        }
        if (!done) // not an interface ?
            BaseA::addBoundary( side );
    }

    return true;
}


template <class T>
bool gsMultiPatch<T>::matchVerticesOnSide (
    const gsMatrix<T> &cc1, const std::vector<boxCorner> &ci1, index_t start,
    const gsMatrix<T> &cc2, const std::vector<boxCorner> &ci2, const gsVector<bool> &matched,
    gsVector<index_t> &dirMap, gsVector<bool>    &dirO,
    T tol,
    index_t reference)
{
    const bool computeOrientation = !(start&(start-1)) && (start != 0); // true if start is a power of 2
    const bool setReference       = start==0;          // if we search for the first point then we set the reference

    const short_t dim = static_cast<short_t>(cc1.rows());

    index_t o_dir = 0, d_dir = 0;

    gsVector<bool> refPar, newPar, newMatched;

    if (computeOrientation)
    {
        // o_dir is the only direction in which the parameters differ between start and 0
        const gsVector<bool> parStart = ci1[start].parameters(dim);
        const gsVector<bool> parRef   = ci1[0].parameters(dim);
        for (; o_dir<dim && parStart(o_dir)==parRef(o_dir)  ;) ++o_dir;
    }

    if (!setReference)
        refPar = ci2[reference].parameters(dim);

    for (size_t j=0;j<ci2.size();++j)
    {
        if( !matched(j) && (cc1.col(ci1[start]-1)-cc2.col(ci2[j]-1)).norm() < tol )
        {
            index_t newRef =   (setReference) ? j : reference;
            if (computeOrientation)
            {
                index_t count=0;
                d_dir = 0;
                newPar =  ci2[j].parameters(dim);
                for (index_t i=0; i< newPar.rows();++i)
                {
                    if ( newPar(i)!=refPar(i) )
                    {
                        d_dir=i;
                        ++count;
                    }
                }
                if (count != 1)
                {
                    // the match is wrong: we are mapping an edge to a diagonal
                    continue;
                }
                dirMap(o_dir) = d_dir;
                dirO  (o_dir) = (static_cast<index_t>(j) > reference);
            }
            if ( start + 1 == static_cast<index_t>( ci1.size() ) )
            {
                // we matched the last vertex, we are done
                return true;
            }
            newMatched = matched;
            newMatched(j) = true;
            if (matchVerticesOnSide( cc1,ci1,start+1,cc2,ci2,newMatched,dirMap,dirO, tol,newRef))
                return true;
        }
    }

    return false;
}


template<class T>
void gsMultiPatch<T>::closeGaps(T tol)
{
    const T tol2 = tol*tol;
    gsMatrix<index_t> bdr1, bdr2; // indices of the boundary control points

    // Create a map which assigns to all meeting patch-local indices a
    // unique global id
    const index_t sz = m_patches.size();
    gsVector<index_t> patchSizes(sz);
    for (index_t i=0; i!= sz; ++i)
        patchSizes[i] = m_patches[i]->coefsSize();

    gsDofMapper mapper(patchSizes);

    for ( iiterator it = iBegin(); it != iEnd(); ++it ) // for all interfaces
    {
        const gsGeometry<T> & p1 = *m_patches[it->first() .patch];
        const gsGeometry<T> & p2 = *m_patches[it->second().patch];

        // Grab boundary control points in matching configuration
        p1.basis().matchWith(*it, p2.basis(), bdr1, bdr2);

        //mapper.matchDofs(it->first().patch, bdr1, it->second().patch, bdr2);
        for (index_t i = 0; i!= bdr1.size(); ++i )
        {
            if ( ( p1.coef(bdr1(i)) - p2.coef(bdr2(i)) ).squaredNorm() > tol2 )
                gsWarn<<"Big gap detected between patches "<< it->first().patch
                      <<" and "<<it->second().patch <<"\n";

            // Match the dofs on the interface
            mapper.matchDof(it->first().patch, bdr1(i,0), it->second().patch, bdr2(i,0) );
        }
    }

    // Finalize the mapper. At this point all patch-local dofs are
    // mapped to unique global indices
    mapper.finalize();

    gsMatrix<T> meanVal;
    std::vector<std::pair<index_t,index_t> > dof;
    const index_t start = mapper.freeSize() - mapper.coupledSize();
    const index_t end   = mapper.freeSize();

    for (index_t i = start; i!= end; ++i) // For all coupled DoFs
    {
        // Get the preimages of this global dof (as pairs (patch,index) )
        mapper.preImage(i, dof);

        // Compute the mean value
        meanVal = m_patches[dof.front().first]->coef(dof.front().second);
        for (size_t k = 1; k!=dof.size(); ++k)
            meanVal += m_patches[dof[k].first]->coef(dof[k].second);
        meanVal.array() /= dof.size();

        // Set involved control points equal to their average value
        for (size_t k = 0; k!=dof.size(); ++k)
            m_patches[dof[k].first]->coef(dof[k].second) = meanVal;
    }
}


template<class T> // to do: move to boundaryInterface
gsAffineFunction<T> gsMultiPatch<T>::getMapForInterface(const boundaryInterface &bi, T scaling) const
{
    if (scaling==0)
        scaling=1;
    gsMatrix<T> box1=m_patches[bi.first().patch ]->support();
    gsMatrix<T> box2=m_patches[bi.second().patch]->support();
    const index_t oDir1 = bi.first() .direction();
    const index_t oDir2 = bi.second().direction();
    const T len1=box1(oDir1,1)-box1(oDir1,0);

    if (bi.second().parameter())
    {
        box2(oDir2,0)=box2(oDir2,1);
        box2(oDir2,1)+=scaling*len1;
    }
    else
    {
        box2(oDir2,1)=box2(oDir2,0);
        box2(oDir2,0)-=scaling*len1;
    }

    return gsAffineFunction<T>(bi.dirMap(bi.first()),bi.dirOrientation(bi.first()) ,box1,box2);
}

template<class T>
void gsMultiPatch<T>::repairInterfaces()
{
    // A bit crude to create a new gsMultiBasis, but this is the fastest
    // way to re-use the already existing functions of gsMultiBasis
    gsMultiBasis<T> multiBasis(*this);
    std::vector< boundaryInterface > bivec = interfaces();

    size_t kmax = 2*bivec.size();
    size_t k = 0;
    bool sthChanged = false;
    bool changed = false;

    do
    {
        sthChanged = false;
        //...keep repairing until nothing changes any more

        // loop over all interfaces
        for( size_t i = 0; i < size_t( bivec.size() ); i++ )
        {
            changed = false;

            std::vector<index_t> refEltsFirst;
            std::vector<index_t> refEltsSecond;

            // For each interface, find the areas/elements that do not match...
            switch( this->dim() )
            {
            case 2:
                changed = multiBasis.template repairInterfaceFindElements<2>( bivec[i], refEltsFirst, refEltsSecond );
                break;
            case 3:
                changed = multiBasis.template repairInterfaceFindElements<3>( bivec[i], refEltsFirst, refEltsSecond );
                break;
            default:
                GISMO_ASSERT(false,"wrong dimension");
            }

            // ...and if there are any found, refine the bases accordingly
            if( changed )
            {
                if( refEltsFirst.size() > 0 )
                {
                    int pi( bivec[i].first().patch );
                    patch(pi).basis().refineElements_withCoefs( patch(pi).coefs(), refEltsFirst );
                }
                if( refEltsSecond.size() > 0 )
                {
                    int pi( bivec[i].second().patch );
                    patch(pi).basis().refineElements_withCoefs( patch(pi).coefs(), refEltsSecond );
                }
            }

            sthChanged = sthChanged || changed;
        }
        k++; // just to be sure this loop cannot go on infinitely
    }
    while( sthChanged && k <= kmax );
}




template<class T>
void gsMultiPatch<T>::locatePoints(const gsMatrix<T> & points,
                                   gsVector<index_t> & pids,
                                   gsMatrix<T> & preim) const
{
    pids.resize(points.cols());
    pids.setConstant(-1); // -1 implies not in the domain
    preim.resize(parDim(), points.cols());//uninitialized by default
    gsMatrix<T> pt, pr, tmp;

    for (index_t i = 0; i!=pids.size(); ++i)
    {
        pt = points.col(i);

        for (size_t k = 0; k!= m_patches.size(); ++k)
        {
            pr = m_patches[k]->parameterRange();
            m_patches[k]->invertPoints(pt, tmp);
            if ( (tmp.array() >= pr.col(0).array()).all()
                 && (tmp.array() <= pr.col(1).array()).all() )
            {
                pids[i] = k;
                preim.col(i) = tmp;
                break;
            }
        }
    }
}

template<class T>
void gsMultiPatch<T>::locatePoints(const gsMatrix<T> & points, index_t pid1,
                                   gsVector<index_t> & pid2, gsMatrix<T> & preim) const
{
    // Assumes points are found on pid1 and possibly on one more patch
    pid2.resize(points.cols());
    pid2.setConstant(-1); // -1 implies not in the domain
    preim.resize(parDim(), points.cols());//uninitialized by default
    gsMatrix<T> pt, pr, tmp;

    for (index_t i = 0; i!=pid2.size(); ++i)
    {
        pt = points.col(i);

        for (size_t k = 0; k!= m_patches.size(); ++k)
        {
            if (pid1==(index_t)k) continue; // skip pid1

            pr = m_patches[k]->parameterRange();
            m_patches[k]->invertPoints(pt, tmp);
            if ( (tmp.array() >= pr.col(0).array()).all()
                 && (tmp.array() <= pr.col(1).array()).all() )
            {
                pid2[i] = k;
                preim.col(i) = tmp;
                break;
            }
        }
    }
}


} // namespace gismo
