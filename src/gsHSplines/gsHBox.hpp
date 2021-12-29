/** @file gsHTensorBasis.hpp

    @brief Provides implementation of HTensorBasis common operations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#pragma once

namespace gismo
{

template <short_t d, class T>
gsHBox<d, T>::gsHBox(const gsHDomainIterator<T,d> * domHIt)
{
    m_basis = static_cast<const gsHTensorBasis<d,T> *>(domHIt->m_basis);
    GISMO_ASSERT(m_basis!=nullptr,"basis is not a gsHTensorBasis");

    m_coords.resize(d,2);
    m_coords.col(0) = domHIt->lowerCorner();
    m_coords.col(1) = domHIt->upperCorner();
    m_center = (m_coords.col(1) + m_coords.col(0))/2;
    _computeIndices();
}

template <short_t d, class T>
gsHBox<d, T>::gsHBox(const typename gsHBox<d,T>::point & low,const typename gsHBox<d,T>::point & upp, index_t level)
:
m_indices(low,upp,level),
m_basis(nullptr)
{}

template <short_t d, class T>
gsHBox<d, T>::gsHBox(const typename gsHBox<d,T>::point & low,const typename gsHBox<d,T>::point & upp, index_t level, const gsHTensorBasis<d,T> * basis)
:
gsHBox(low,upp,level)
{
    m_basis = basis;
}

template <short_t d, class T>
gsHBox<d, T>::gsHBox(const gsAabb<d,index_t> & box)
:
m_indices(box),
m_basis(nullptr)
{

}

template <short_t d, class T>
gsHBox<d, T>::gsHBox(const gsAabb<d,index_t> & box, const gsHTensorBasis<d,T> * basis)
:
gsHBox(box)
{
    m_basis = basis;
}

template <short_t d, class T>
gsHBox<d, T>::gsHBox(const std::vector<index_t> & indices)
:
m_basis(nullptr)
{
    GISMO_ENSURE(indices.size()==2*d+1,"Index size is wrong");
    typename gsHBox<d,T>::point low, upp;
    for (index_t k=0; k!=d; k++)
    {
        low[k] = indices[k+1];
        upp[k] = indices[k+d+1];
    }

    m_indices = gsAabb<d,index_t>(low,upp,indices[0]);
    //  = gsAsVector<index_t,d>()
    // upp;
}

template <short_t d, class T>
gsHBox<d, T>::gsHBox( const gsHBox<d,T> & other )
{
    operator=(other);
}

template <short_t d, class T>
gsHBox<d, T>::gsHBox( gsHBox<d,T> && other )
{
    operator=(give(other));
}

template <short_t d, class T>
gsHBox<d,T> & gsHBox<d, T>::operator= ( const gsHBox<d,T> & other )
{
    if (this!=&other)
    {
        m_indices = other.m_indices;
        m_coords  = other.m_coords;
        m_center  = other.m_center;
        m_basis   = other.m_basis;
    }
    return *this;
}

template <short_t d, class T>
gsHBox<d,T> & gsHBox<d, T>::operator= ( gsHBox<d,T> && other )
{
    m_indices = give(other.m_indices);
    m_coords  = give(other.m_coords);
    m_center  = give(other.m_center);
    m_basis   = give(other.m_basis);
    return *this;
}

template <short_t d, class T>
bool gsHBox<d, T>::isContained(const gsHBox<d,T> & other) const
{
    bool res = true;
    res &= this->level() == other.level();
    for (index_t i=0; i!=d && res; i++)
    {
        res &= this->lowerIndex().at(i) <= other.lowerIndex().at(i);
        res &= this->upperIndex().at(i) >= other.upperIndex().at(i);
    }
    return res;
}

template <short_t d, class T>
bool gsHBox<d, T>::contains(const gsHBox<d,T> & other) const
{
    bool res = true;
    res &= this->level() == other.level();
    for (index_t i=0; i!=d && res; i++)
    {
        res &= this->lowerIndex().at(i) >= other.lowerIndex().at(i);
        res &= this->upperIndex().at(i) <= other.upperIndex().at(i);
    }
    return res;
}

template <short_t d, class T>
bool gsHBox<d, T>::isSame(const gsHBox<d,T> & other) const
{
    bool res = true;
    res &= this->level() == other.level();
    for (index_t i=0; i!=d && res; i++)
    {
        res &= this->lowerIndex().at(i) == other.lowerIndex().at(i);
        res &= this->upperIndex().at(i) == other.upperIndex().at(i);
    }
    return res;
}

template <short_t d, class T>
bool gsHBox<d, T>::isActive()
{
    return m_basis->getLevelAtPoint(this->getCenter()) == this->level();
}

template <short_t d, class T>
const gsMatrix<T> & gsHBox<d, T>::getCoordinates()
{
    if (m_coords.size()==0)
        _computeCoordinates();
    return m_coords;
}

template <short_t d, class T>
const gsMatrix<T> & gsHBox<d, T>::getCenter()
{
    if (m_center.size()==0)
        _computeCoordinates();
    return m_center;
}

template <short_t d, class T>
gsVector<T,d>  gsHBox<d, T>::lowerCorner() const { return m_coords.col(0); }
template <short_t d, class T>
gsVector<T,d>  gsHBox<d, T>::upperCorner() const { return m_coords.col(1); }

template <short_t d, class T>
typename gsHBox<d,T>::point gsHBox<d, T>::lowerIndex() const { return m_indices.first;  }
template <short_t d, class T>
typename gsHBox<d,T>::point gsHBox<d, T>::upperIndex() const { return m_indices.second; }

template <short_t d, class T>
index_t gsHBox<d, T>::level() const { return m_indices.level; }

template <short_t d, class T>
gsHBox<d,T> gsHBox<d, T>::getParent() const
{
    GISMO_ENSURE(this->level()>0,"Box is at ground level and has no parent");
    gsAabb<d,index_t> box = _elevateBox(m_indices);
    return gsHBox<d,T>(box,m_basis);
}

template <short_t d, class T>
gsHBox<d,T> gsHBox<d, T>::getAncestor(index_t k) const
{
    index_t lvl = this->level();
    GISMO_ASSERT(lvl > k,"Current level should be larger than requested level, l = "<<lvl<<", k = "<<k);
    GISMO_ASSERT(k > 0,"Current level should be larger than requested level, l = "<<lvl<<", k = "<<k);
    gsHBox<d,T> parent = this->getParent();
    gsHBox<d,T> ancestor;
    if (k < lvl - 1)
    {
        ancestor = parent.getAncestor(k);
        return ancestor;
    }
    else
        return parent;
}

template <short_t d, class T>
typename gsHBox<d,T>::HContainer gsHBox<d, T>::toContainer()
{
    HContainer container;
    container.resize(this->level()+1);
    container.at(this->level()).push_back(*this);
    return container;
}

template <short_t d, class T>
typename gsHBox<d, T>::HContainer gsHBox<d, T>::getSupportExtension()
{
    GISMO_ENSURE(m_basis!=nullptr,"Basis is undefined");
    this->_computeCoordinates();
    index_t lvl = this->level();
    // Compute actives
    gsMatrix<index_t> acts;
    m_basis->tensorLevel(lvl).active_into(m_center,acts);

    // Support extension
    gsMatrix<T> cells(d,2*acts.rows());
    gsMatrix<T> support;
    gsAabb<d,index_t> aabb;
    HContainer container;
    container.resize(lvl+1);
    for (index_t act = 0; act!=acts.rows(); act++)
    {
        support = m_basis->tensorLevel(lvl).support(acts(act,0));
        aabb    = _computeIndices(support,lvl);
        container[lvl].push_back(gsHBox<d,T>(aabb,m_basis));
    }
    return container;
}

template <short_t d, class T>
typename gsHBox<d,T>::HContainer gsHBox<d, T>::getMultiLevelSupportExtension(index_t k)
{
    index_t lvl = this->level();
    GISMO_ASSERT(k <= lvl,"Ancestor must be requested from a level lower than the current level (k <= l), but k="<<k<<" and lvl="<<lvl);
    if (k==lvl)
    {
        return this->getSupportExtension();
    }
    else
    {
        gsHBox<d,T> ancestor = this->getAncestor(k);
        return ancestor.getSupportExtension();
    }
}

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBox<d, T>::getHneighborhood(index_t m)
{
    HContainer extension;
    Container neighborhood;

    index_t lvl = this->level();
    index_t k = lvl - m + 1;
    if (k>=0)
    {
        // Get multi level support extension on level k
        extension = this->getMultiLevelSupportExtension(k);
        // Eliminate elements which are too low
        neighborhood = extension.at(k);
    }
    return neighborhood;
}

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBox<d, T>::getTneighborhood(index_t m)
{
    Container   neighborhood;
    HContainer  parents, extension;

    index_t lvl = this->level();
    index_t k = lvl - m + 2;
    if (k>=0)
    {
        // Get multi level support extension on level k
        extension = this->getMultiLevelSupportExtension(k);
        // Eliminate elements which are too low
        parents = _getParents(extension);
        neighborhood = parents.at(k-1);
    }
    return neighborhood;
}

template <short_t d, class T>
std::ostream& gsHBox<d, T>::print( std::ostream& os ) const
{
    os  <<"Cell with dimension "<<d
        <<" on level "<<m_indices.level<<"; "
        <<"("<<m_indices.first.transpose()<<")"
        <<" -- "
        <<"("<<m_indices.second.transpose()<<")";
    return os;
}

template <short_t d, class T>
void gsHBox<d, T>::_computeCoordinates()
{
    GISMO_ENSURE(m_basis!=nullptr,"Basis is not provided");
    m_coords.resize(d,2);
    gsVector<T> low(d), upp(d);

    for (index_t i=0; i!=d; i++)
    {
        const gsKnotVector<T> & kv = m_basis->tensorLevel(m_indices.level).knots(i);
        low[i] = kv.uValue(m_indices.first.at(i));
        upp[i] = kv.uValue(m_indices.second.at(i));
    }
    m_coords.col(0) = low;
    m_coords.col(1) = upp;

    m_center = (m_coords.col(0) + m_coords.col(1))/2;
}

template <short_t d, class T>
void gsHBox<d, T>::_computeIndices()
{
    GISMO_ENSURE(m_basis!=nullptr,"Basis is not provided");
    m_indices = _computeIndices(m_coords);
}

template <short_t d, class T>
gsAabb<d,index_t> gsHBox<d, T>::_computeIndices(const gsMatrix<T> & coords, index_t level)
{
    GISMO_ENSURE(m_basis!=nullptr,"Basis is not provided");
    typename gsHBox<d,T>::point low,upp;
    for(index_t j = 0; j < d;j++)
    {
        // Convert the parameter coordinates to (unique) knot indices
        const gsKnotVector<T> & kv = m_basis->tensorLevel(level).knots(j);
        int k1 = (std::upper_bound(kv.domainUBegin(), kv.domainUEnd(),
                                   coords(j,0) ) - 1).uIndex();
        int k2 = (std::upper_bound(kv.domainUBegin(), kv.domainUEnd()+1,
                                   coords(j,1) ) - 1).uIndex();

        // Trivial cells trigger some refinement
        if ( k1 == k2)
        {
            if (0!=k1) {--k1;}
            ++k2;
        }

        // Store the data...
        low.at(j)  = k1;
        upp.at(j) = k2;
    }
    return gsAabb<d,index_t>(low,upp,level);
}

template <short_t d, class T>
gsAabb<d,index_t> gsHBox<d, T>::_computeIndices(const gsMatrix<T> & coords)
{
    GISMO_ENSURE(m_basis!=nullptr,"Basis is not provided");
    gsMatrix<T> center = (coords.col(0) + coords.col(1))/2;
    index_t level = m_basis->getLevelAtPoint(center);
    return _computeIndices(coords,level);
}

template <short_t d, class T>
gsAabb<d,index_t> gsHBox<d, T>::_computeIndices(const gsMatrix<T> & coords, const gsMatrix<T> & center)
{
    GISMO_ENSURE(m_basis!=nullptr,"Basis is not provided");
    index_t level = m_basis->getLevelAtPoint(center);
    return _computeIndices(coords,level);
}

template <short_t d, class T>
gsAabb<d,index_t> gsHBox<d, T>::_elevateBox(const gsAabb<d,index_t> & box) const
{
    typename gsHBox<d,T>::point low,upp;
    for (index_t i = 0; i!=d; i++)
    {
        low.at(i) = box.first .at(i) / 2;
        upp.at(i) = box.second.at(i) / 2 + (index_t)(box.second.at(i) % 2 != 0);
    }
    return gsAabb<d,index_t>(low,upp,box.level-1);
}

template <short_t d, class T>
typename gsHBox<d,T>::HContainer gsHBox<d, T>::_getParents(HContainer & container)
{
    HContainer result;
    result.resize(container.size()-1);

    // Handle level 0 separately (operation cannot be performed on level 0)
    GISMO_ASSERT(container[0].size()==0,"Boxes at level 0 cannot have a parent. Did something go wrong? You can run check() to see if the boxes are allocated coorectly");

    HIterator resIt = result.begin();
    for (HIterator hit = std::next(container.begin()); hit!=container.end(); hit++, resIt++)
        for (Iterator it=hit->begin(); it!=hit->end(); it++)
            resIt->push_back(*it);

    return result;
}

} // namespace gismo
