/** @file gsHBox.hpp

    @brief Provides implementation of HTensorBasis common operations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (TU Delft 2019-...)
*/

#pragma once

#include <gsUtils/gsCombinatorics.h>
#include <gsHSplines/gsHBSplineBasis.h>
#include <gsHSplines/gsTHBSplineBasis.h>

namespace gismo
{

template <short_t d, class T>
gsHBox<d, T>::gsHBox(const gsHDomainIterator<T,d> * domHIt)
:
gsHBox(domHIt,-1)
{ }

template <short_t d, class T>
gsHBox<d, T>::gsHBox(const gsHDomainIterator<T,d> * domHIt, const index_t pid)
:
m_pid(pid),
m_error(0),
m_index(-1),
m_marked(false)
{
    m_basis = nullptr;
    m_basis = static_cast<const gsHTensorBasis<d,T> *>(domHIt->m_basis);
    GISMO_ASSERT(m_basis!=nullptr,"basis is not a gsHTensorBasis");

    m_coords.resize(d,2);
    m_coords.col(0) = domHIt->lowerCorner();
    m_coords.col(1) = domHIt->upperCorner();
    m_center = (m_coords.col(1) + m_coords.col(0))/2;
    m_indices = _computeIndices(m_coords,m_center);
}

template <short_t d, class T>
gsHBox<d, T>::gsHBox(const typename gsHBox<d,T>::point & low,const typename gsHBox<d,T>::point & upp, index_t level, const gsHTensorBasis<d,T> * basis, const index_t pid)
:
m_indices(low,upp,level),
m_pid(pid),
m_error(0),
m_index(-1),
m_marked(false)
{
    m_basis = basis;
}

template <short_t d, class T>
gsHBox<d, T>::gsHBox(const gsAabb<d,index_t> & box, const gsHTensorBasis<d,T> * basis, const index_t pid)
:
m_indices(box),
m_pid(pid),
m_error(0),
m_index(-1),
m_marked(false)
{
    m_basis = basis;
}

template <short_t d, class T>
gsHBox<d, T>::gsHBox(const std::vector<index_t> & indices, const gsHTensorBasis<d,T> * basis, const index_t pid)
:
m_pid(pid),
m_error(0),
m_index(-1),
m_marked(false)
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

    m_basis = basis;
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
        m_pid     = other.m_pid;
        m_coords  = other.m_coords;
        m_center  = other.m_center;
        m_basis   = other.m_basis;
        m_error   = other.m_error;
        m_marked  = other.m_marked;
        m_index   = other.m_index;
    }
    return *this;
}

template <short_t d, class T>
gsHBox<d,T> & gsHBox<d, T>::operator= ( gsHBox<d,T> && other )
{
    m_indices = give(other.m_indices);
    m_pid     = give(other.m_pid);
    m_coords  = give(other.m_coords);
    m_center  = give(other.m_center);
    m_basis   = give(other.m_basis);
    m_error   = give(other.m_error);
    m_marked  = give(other.m_marked);
    m_index   = give(other.m_index);
    return *this;
}

template <short_t d, class T>
bool gsHBox<d, T>::isContained(const gsHBox<d,T> & other) const
{
    // bool res = true;
    // res &= this->level() == other.level();
    // for (index_t i=0; i!=d && res; i++)
    // {
    //     res &= this->lowerIndex().at(i) <= other.lowerIndex().at(i);
    //     res &= this->upperIndex().at(i) >= other.upperIndex().at(i);
    // }
    // return res;


    return gsHBoxIsContained<d,T>()(*this,other);

}

template <short_t d, class T>
bool gsHBox<d, T>::contains(const gsHBox<d,T> & other) const
{
    // bool res = true;
    // res &= this->level() == other.level();
    // for (index_t i=0; i!=d && res; i++)
    // {
    //     res &= this->lowerIndex().at(i) >= other.lowerIndex().at(i);
    //     res &= this->upperIndex().at(i) <= other.upperIndex().at(i);
    // }
    // return res;

    return gsHBoxContains<d,T>()(*this,other);
}

template <short_t d, class T>
bool gsHBox<d, T>::isSame(const gsHBox<d,T> & other) const
{
    bool res = true;
    res &= this->patch() == other.patch();
    res &= this->level() == other.level();
    for (index_t i=0; i!=d && res; i++)
    {
        res &= this->lowerIndex().at(i) == other.lowerIndex().at(i);
        res &= this->upperIndex().at(i) == other.upperIndex().at(i);
    }
    return res;
}

template <short_t d, class T>
bool gsHBox<d, T>::isActive() const
{
    return m_basis->getLevelAtPoint(this->getCenter()) == this->level();
}

template <short_t d, class T>
bool gsHBox<d, T>::isActiveOrContained() const
{
    return m_basis->getLevelAtPoint(this->getCenter()) >= this->level();
}

template <short_t d, class T>
index_t gsHBox<d, T>::levelInCenter() const
{
    return m_basis->getLevelAtPoint(this->getCenter());
}

template <short_t d, class T>
const gsMatrix<T> & gsHBox<d, T>::getCoordinates() const
{
    GISMO_ASSERT(m_coords.size()!=0,"Coordinates have not been computed! Please call computeCoordinates()");
    return m_coords;
}

template <short_t d, class T>
const gsMatrix<T> & gsHBox<d, T>::getCenter() const
{
    // GISMO_ASSERT(m_center.size()!=0,"Center has not been computed! Please call computeCenter()");
    if (m_center.size()==0)
        this->computeCenter();
    return m_center;
}

template <short_t d, class T>
gsVector<T,d>  gsHBox<d, T>::lowerCorner() const
{
    GISMO_ASSERT(m_coords.size()!=0,"Coordinates have not been computed! Please call computeCoordinates()");
    return m_coords.col(0);
}

template <short_t d, class T>
gsVector<T,d>  gsHBox<d, T>::upperCorner() const
{
    GISMO_ASSERT(m_coords.size()!=0,"Coordinates have not been computed! Please call computeCoordinates()");
    return m_coords.col(1);
}

template <short_t d, class T>
T gsHBox<d, T>::getCellSize() const
{
    return (upperCorner() - lowerCorner()).norm();
}

template <short_t d, class T>
T gsHBox<d, T>::getMinCellLength() const
{
    return (upperCorner() - lowerCorner()).minCoeff();
}

template <short_t d, class T>
T gsHBox<d, T>::getMaxCellLength() const
{
    return (upperCorner() - lowerCorner()).maxCoeff();
}

template <short_t d, class T>
const typename gsHBox<d,T>::point & gsHBox<d, T>::lowerIndex() const { return m_indices.first;  }
template <short_t d, class T>
const typename gsHBox<d,T>::point & gsHBox<d, T>::upperIndex() const { return m_indices.second; }

template <short_t d, class T>
index_t gsHBox<d, T>::patch() const { return m_pid; }

template <short_t d, class T>
index_t gsHBox<d, T>::level() const { return m_indices.level; }

template <short_t d, class T>
void gsHBox<d, T>::setError(T error) { m_error = error; }

template <short_t d, class T>
T gsHBox<d, T>::error() const { return m_error; }

template <short_t d, class T>
void gsHBox<d, T>::setIndex(index_t index) { m_index = index; }

template <short_t d, class T>
index_t gsHBox<d, T>::index() const { return m_index; }

template <short_t d, class T>
void gsHBox<d, T>::mark() { m_marked = true; }

template <short_t d, class T>
void gsHBox<d, T>::unmark() { m_marked = false; }

template <short_t d, class T>
bool gsHBox<d, T>::marked() const { return m_marked; }

template <short_t d, class T>
void gsHBox<d, T>::setMark(bool mark) { m_marked = mark; }

template <short_t d, class T>
gsHBox<d,T> gsHBox<d, T>::getParent() const
{
    GISMO_ENSURE(this->level()>0,"Box is at ground level and has no parent");
    gsAabb<d,index_t> box = _elevateBox(m_indices);
    return gsHBox<d,T>(box,m_basis,m_pid);
}

template <short_t d, class T>
gsHBox<d,T> gsHBox<d, T>::getAncestor(index_t k) const
{
    index_t lvl = this->level();
    // GISMO_ASSERT(lvl > k,"Current level should be larger than requested level, l = "<<lvl<<", k = "<<k);
    GISMO_ASSERT(k >= 0,"Requested ancestor k = "<<k<<" should be greater than or equal to 0");
    GISMO_ASSERT(lvl >= 0,"Level lvl = "<<lvl<<" should be larger than 0");
    // GISMO_ASSERT(k > 0,"Current level should be larger than requested level, l = "<<lvl<<", k = "<<k);
    gsHBox<d,T> parent = this->getParent();
    gsHBox<d,T> ancestor;
    if (k < lvl - 1)
    {
        ancestor = parent.getAncestor(k);
        return ancestor;
    }
    else if (k==lvl-1)
        return parent;
    else if (k==lvl)
        return *this;
    else
        GISMO_ERROR("Cannot find ancestor with index k="<<k<<" on level l="<<lvl<<". Something went wrong?");
}

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBox<d, T>::getChildren() const
{
    gsAabb<d,index_t> box = _lowerBox(m_indices);
    gsHBox<d,T> childRegion(box,m_basis,m_pid);
    return childRegion.toUnitBoxes();
}

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBox<d, T>::getDescendants(index_t k) const
{
    index_t lvl = this->level();
    // GISMO_ASSERT(lvl > k,"Current level should be larger than requested level, l = "<<lvl<<", k = "<<k);
    GISMO_ASSERT(k >= 0,"Requested descendant k = "<<k<<" should be larger than or equal to 0");
    GISMO_ASSERT(lvl >= 0,"Level lvl = "<<lvl<<" should be larger than 0");
    Container descendants(std::pow(std::pow(2,d),k-lvl));
    if (k==lvl)
        descendants.push_back(*this);
    if (k==lvl+1)
        descendants = this->getChildren();
    else
    {
        Container tmp = this->getChildren(); // children have level lvl+1
        for (index_t k_tmp = lvl+1; k_tmp!=k; k_tmp++)
        {
            descendants.clear();
            // descendants.reserve()
            for (Iterator it=tmp.begin(); it!=tmp.end(); it++)
            {
                Container tmp2 = it->getChildren();
                for (Iterator it2=tmp2.begin(); it2!=tmp2.end(); it2++)
                    descendants.push_back(*it2);
            }
            tmp = descendants;
        }
    }
    return descendants;
}

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBox<d, T>::getSiblings() const
{
    auto isthis = [this](const gsHBox<d,T> & other) { return gsHBoxEqual<d,T>()(other,*this); };
    gsHBox<d,T> parent = this->getParent();
    typename gsHBox<d,T>::Container children = parent.getChildren();
    Iterator toErase = std::find_if(children.begin(),children.end(),isthis);
    if (toErase!=children.end())
        children.erase(toErase);
    else
        GISMO_ERROR("Something went wrong");
    return children;
}

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBox<d, T>::toContainer()
{
    Container container;
    container.push_back(*this);
    return container;
}

template <short_t d, class T>
typename gsHBox<d,T>::HContainer gsHBox<d, T>::toHContainer()
{
    HContainer container;
    container.resize(this->level()+1);
    container.at(this->level()).push_back(*this);
    return container;
}

template <short_t d, class T>
typename gsHBox<d, T>::Container gsHBox<d, T>::getSupportExtension()
{
    this->computeCenter();
    index_t lvl = this->level();
    // Compute actives
    gsMatrix<index_t> acts;
    m_basis->tensorLevel(lvl).active_into(m_center,acts);

    // Support extension
    gsMatrix<T> cells(d,2*acts.rows());
    gsMatrix<T> support;
    gsHBox<d,T> supportBox;
    gsAabb<d,index_t> aabb;
    Container container;
    Container tmpContainer;
    for (index_t act = 0; act!=acts.rows(); act++)
    {
        support = m_basis->tensorLevel(lvl).support(acts(act,0));
        aabb    = _computeIndices(support,lvl);

        supportBox = gsHBox<d,T>(aabb,m_basis,m_pid);

        // Split the boxes into index interval with coordinate delta 1
        tmpContainer = supportBox.toUnitBoxes();
        for (cIterator it = tmpContainer.begin(); it!=tmpContainer.end(); it++)
            container.push_back(*it);
    }
    // Remove duplicates
    Container container2 = _makeUnique(container);
    return container2;
}

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBox<d, T>::getMultiLevelSupportExtension(index_t k)
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
    Container extension;
    Container neighborhood;

    index_t lvl = this->level();
    index_t k = lvl - m + 1;
    if (k>=0)
    {
        // Get multi level support extension on level k
        extension = this->getMultiLevelSupportExtension(k);

        // Eliminate elements which are too low
        for (Iterator it = extension.begin(); it!=extension.end(); it++)
        {
            it->computeCenter(); // needed to check active
            if (it->isActive())
                neighborhood.push_back(*it);
        }
    }
    return neighborhood;
}

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBox<d, T>::getTneighborhood(index_t m)
{
    // Everything is in the same level so we can use normal Container
    Container   neighborhood;
    Container  parents, extension;

    index_t lvl = this->level();
    index_t k = lvl - m + 2;
    if (k-1>=0)
    {
        // Get multi level support extension on level k
        extension = this->getMultiLevelSupportExtension(k);
        // Eliminate elements which are too low
        parents = _getParents(extension);

        for (Iterator it = parents.begin(); it!=parents.end(); it++)
        {
            it->computeCenter(); // needed to check active
            if (it->isActive())
                neighborhood.push_back(*it);
        }
        // neighborhood = parents;
    }
    return neighborhood;
}

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBox<d, T>::getNeighborhood(index_t m)
{
    Container result;
    if (dynamic_cast<const gsTHBSplineBasis<d,T>*>(m_basis))
        result = this->getNeighborhood<Neighborhood::T>(m);
    else if (dynamic_cast<const gsHBSplineBasis<d,T>*>(m_basis))
        result = this->getNeighborhood<Neighborhood::H>(m);
    else
        GISMO_ERROR("Basis type should be gsTHBSplineBasis or gsHBSplineBasis");
    return result;
}

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBox<d, T>::getCextension(index_t m)
{
    // See Definition 3.5 from [1]
    // [1] Carraturo, M., Giannelli, C., Reali, A. & Vázquez, R.
    // Suitably graded THB-spline refinement and coarsening: Towards an adaptive isogeometric analysis of additive manufacturing processes.
    // Comput. Methods Appl. Mech. Eng. 348, 660–679 (2019).
    //
    // This function returns the coarsening neighborhood of \a this, which is an element to be reactivated
    Container extension, children, childExtension, descendants, result;
    index_t targetLvl = this->level() + m;

    children = this->getChildren();
    for (typename gsHBox<d,T>::Iterator child=children.begin(); child!=children.end(); child++)
    {
        childExtension = child->getSupportExtension();
        extension = gsHBoxUtils<d,T>::Union(extension,childExtension);
    }

    // All elements in the support extension of the children (so at level l+1)
    extension = gsHBoxUtils<d,T>::Union(extension,children);
    // question: needed MINUS the siblings??
    // extension = gsHBoxUtils<d,T>::Difference(extension,children);

    // result = extension;
    // gsDebugVar(targetLvl);
    for (typename gsHBox<d,T>::Iterator eIt = extension.begin(); eIt!=extension.end(); eIt++)
    {
        descendants = eIt->getDescendants(targetLvl);
        result.insert(result.end(), descendants.begin(), descendants.end());
    }
    return result;
}

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBox<d, T>::getCneighborhood(index_t m)
{
    // See Definition 3.5 from [1]
    // [1] Carraturo, M., Giannelli, C., Reali, A. & Vázquez, R.
    // Suitably graded THB-spline refinement and coarsening: Towards an adaptive isogeometric analysis of additive manufacturing processes.
    // Comput. Methods Appl. Mech. Eng. 348, 660–679 (2019).
    //
    // This function returns the coarsening neighborhood of \a this, which is an element to be reactivated
    Container descendants, result;
    descendants = this->getCextension(m);

    for (typename gsHBox<d,T>::Iterator dIt = descendants.begin(); dIt!=descendants.end(); dIt++)
    {
        dIt->computeCenter(); // needed to check active
        if (dIt->levelInCenter()>=dIt->level()) // the level is even larger (i.e. even higher decendant)!!
            result.push_back(*dIt);
    }
    return result;
}

template <short_t d, class T>
std::ostream& gsHBox<d, T>::print( std::ostream& os ) const
{
    os  <<"Cell on patch "<<m_pid
        <<" on level "<<m_indices.level<<". "
        <<"\nIndices:\n"
        <<"("<<m_indices.first.transpose()<<")"
        <<" -- "
        <<"("<<m_indices.second.transpose()<<")";
    if (m_coords.cols()!=0)
    {
        os  <<"\nKnot values:\n"
            <<"("<<m_coords.col(0).transpose()<<")"
            <<" -- "
            <<"("<<m_coords.col(1).transpose()<<")";
    }
    os  <<"\nmarked = "<<m_marked<<"";
    os  <<"\nerror = "<<m_error<<"";
    return os;
}

template <short_t d, class T>
void gsHBox<d, T>::computeCoordinates() const
{
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
}

template <short_t d, class T>
void gsHBox<d, T>::computeCenter() const
{
    this->computeCoordinates();
    m_center = (m_coords.col(0) + m_coords.col(1))/2;
}

template <short_t d, class T>
void gsHBox<d, T>::_computeIndices()
{
    m_indices = _computeIndices(m_coords);
}

template <short_t d, class T>
gsAabb<d,index_t> gsHBox<d, T>::_computeIndices(const gsMatrix<T> & coords, index_t level)
{
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
    gsMatrix<T> center = (coords.col(0) + coords.col(1))/2;
    index_t level = m_basis->getLevelAtPoint(center);
    return _computeIndices(coords,level);
}

template <short_t d, class T>
gsAabb<d,index_t> gsHBox<d, T>::_computeIndices(const gsMatrix<T> & coords, const gsMatrix<T> & center)
{
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
gsAabb<d,index_t> gsHBox<d, T>::_lowerBox(const gsAabb<d,index_t> & box) const
{
    typename gsHBox<d,T>::point low,upp;
    for (index_t i = 0; i!=d; i++)
    {
        low.at(i) = box.first .at(i) * 2;
        upp.at(i) = box.second.at(i) * 2;
    }
    return gsAabb<d,index_t>(low,upp,box.level+1);
}

template <short_t d, class T>
typename gsHBox<d,T>::HContainer gsHBox<d, T>::_getParents(HContainer & container) const
{
    HContainer result;
    result.resize(container.size()-1);

    // Handle level 0 separately (operation cannot be performed on level 0)
    GISMO_ASSERT(container[0].size()==0,"Boxes at level 0 cannot have a parent. Did something go wrong? You can run check() to see if the boxes are allocated coorectly");

    HIterator resIt = result.begin();
    for (HIterator hit = std::next(container.begin()); hit!=container.end(); hit++, resIt++)
        for (Iterator it=hit->begin(); it!=hit->end(); it++)
            resIt->push_back(it->getParent());

    return result;
}

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBox<d, T>::_getParents(Container & container) const
{
    Container result;
    for (Iterator it=container.begin(); it!=container.end(); it++)
        result.push_back(it->getParent());

    return result;
}


template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBox<d, T>::toUnitBoxes() const
{
    point low, upp, cur, curupp, ones;
    ones = gsVector<index_t,d>::Ones();
    cur = low = this->lowerIndex();
    upp       = this->upperIndex();
    // index_t nboxes = (upp-low).prod();
    Container result;

    bool next = true;
    while (next)
    {
        curupp = cur + ones;
        result.push_back(gsHBox<d,T>(cur,curupp,this->level(),m_basis,m_pid));
        next = nextLexicographic(cur,low,upp);
    }
    return result;
}

template <short_t d, class T>
typename gsHBox<d,T>::RefBox gsHBox<d, T>::toBox() const
{
    std::vector<index_t> result(5);
    result[0] = this->level();
    result[1] = this->lowerIndex()[0];
    result[2] = this->lowerIndex()[1];
    result[3] = this->upperIndex()[0];
    result[4] = this->upperIndex()[1];
    return result;
}

template <short_t d, class T>
typename gsHBox<d,T>::RefBox gsHBox<d, T>::toRefBox(index_t targetLevel) const
{
    GISMO_ASSERT(targetLevel > this->level(),"Target level must be larger than the current level, but "<<targetLevel<<"=<"<<this->level());
    std::vector<index_t> result(2*d+1);
    index_t diff = targetLevel - this->level();
    result[0] = this->level() + diff;
    index_t lowerIndex, upperIndex, degree; //, maxKtIndex
    for (index_t i = 0; i!=d; i++)
    {
        degree = m_basis->degree(i);
        // maxKtIndex = m_basis->tensorLevel(this->level()+1).knots(i).size();

        if (degree % 2 == 1 && degree>1)
        {
            lowerIndex = this->lowerIndex()[i]*std::pow(2,diff);
            ( lowerIndex < (degree-1)/2-1 ? lowerIndex=0 : lowerIndex-=(degree-1)/2-1 );
            result[i+1] = lowerIndex;
            upperIndex = this->upperIndex()[i]*std::pow(2,diff);
            // ( upperIndex + (degree)/2+1 >= maxKtIndex ? upperIndex=maxKtIndex-1 : upperIndex+=(degree)/2+1);
            result[d+i+1] = upperIndex;
        }
        else
        {
            lowerIndex = this->lowerIndex()[i]*std::pow(2,diff);
            ( lowerIndex < (degree-1)/2 ? lowerIndex=0 : lowerIndex-=(degree-1)/2 );
            result[i+1] = lowerIndex;
            upperIndex = this->upperIndex()[i]*std::pow(2,diff);
            // ( upperIndex + (degree)/2 >= maxKtIndex ? upperIndex=maxKtIndex-1 : upperIndex+=(degree)/2);
            result[d+i+1] = upperIndex;
        }
    }
    return result;
}

template <short_t d, class T>
typename gsHBox<d,T>::RefBox gsHBox<d, T>::toCrsBox(index_t targetLevel) const
{
    GISMO_ASSERT(targetLevel < this->level(),"Target level must be smaller than the current level, but "<<targetLevel<<">="<<this->level());
    std::vector<index_t> result(2*d+1);
    index_t diff = this->level() - targetLevel;
    result[0] = this->level() - diff;
    for (index_t i = 0; i!=d; i++)
    {
        result[i+1] = this->lowerIndex()[i]/std::pow(2,diff);
        result[d+i+1] = this->upperIndex()[i]/std::pow(2,diff);
    }
    return result;
}


template <short_t d, class T>
typename gsHBox<d,T>::HContainer gsHBox<d, T>::boxUnion(const HContainer & container1, const HContainer & container2) const
{
    HContainer result, region1, region2;

    region1 = container1;
    region2 = container2;

    index_t lmax = std::max(region1.size(),region2.size());
    region1.resize(lmax);
    region2.resize(lmax);
    result.resize(lmax);

    for (index_t l = 0; l!=lmax; l++)
        result[l] = _boxUnion(region1[l],region2[l]);

    return result;
}

template <short_t d, class T>
typename gsHBox<d, T>::Container gsHBox<d, T>::_boxUnion(const Container & container1, const Container & container2) const
{
    // SortedContainer sortedResult;

    // SortedContainer scontainer1 = gsHBoxSort<d,T>(container1);
    // SortedContainer scontainer2 = gsHBoxSort<d,T>(container2);


    // struct
    // {
    //     bool operator()(const gsHBox<d,T> & a, const gsHBox<d,T> & b) const
    //     {
    //         return
    //          (a.patch() < b.patch())
    //         ||
    //         ((a.patch() == b.patch()) &&
    //          (a.level() < b.level())     )
    //         ||
    //         ((a.patch() == b.patch()) &&
    //          (a.level() == b.level()) &&
    //          std::lexicographical_compare(  a.lowerIndex().begin(), a.lowerIndex().end(),
    //                                     b.lowerIndex().begin(), b.lowerIndex().end())   )
    //         ||
    //         ((a.patch() == b.patch()) &&
    //          (a.level() == b.level()) &&
    //          (a.lowerIndex() == b.lowerIndex()) &&
    //          std::lexicographical_compare(  a.upperIndex().begin(), a.upperIndex().end(),
    //                                     b.upperIndex().begin(), b.upperIndex().end())    );
    //     };
    // }
    // comp;

    // sortedResult.reserve(scontainer1.size() + scontainer2.size());
    // if (scontainer1.size()!=0 && scontainer2.size()!=0)
    // {
    //     // First sort (otherwise union is wrong)
    //     std::sort(scontainer1.begin(),scontainer1.end(),comp);
    //     std::sort(scontainer2.begin(),scontainer2.end(),comp);

    //     std::set_union( scontainer1.begin(),scontainer1.end(),
    //                     scontainer2.begin(),scontainer2.end(),
    //                     std::inserter(sortedResult,sortedResult.begin()),
    //                     comp);
    // }
    // else if (scontainer1.size()!=0 && container2.size()==0)
    // {
    //     scontainer1 = SortedContainer (container1.begin(), container1.end());
    //     scontainer2 = SortedContainer (container2.begin(), container2.end());
    //     sortedResult.insert(sortedResult.end(),scontainer1.begin(),scontainer1.end());
    // }
    // else if (scontainer1.size()==0 && container2.size()!=0)
    // {

    //     sortedResult.insert(sortedResult.end(),scontainer2.begin(),scontainer2.end());
    // }
    // else    { /* Do nothing */ }

    return gsHBoxUtils<d,T>::Union(container1,container2);
}

template <short_t d, class T>
typename gsHBox<d, T>::Container gsHBox<d, T>::_makeUnique(const Container & container) const
{
    // SortedContainer scontainer = gsHBoxSort<d,T>(container);

    // struct
    // {
    //     bool operator()(const gsHBox<d,T> & a, const gsHBox<d,T> & b) const
    //     {
    //         return a.isSame(b);
    //     };
    // }
    // pred;

    // // Get unique entries
    // typename SortedContainer::iterator it = std::unique(scontainer.begin(),scontainer.end(),pred);
    // scontainer.resize(distance(scontainer.begin(), it));
    // Container result(scontainer.begin(),scontainer.end());
    // return result;

    return gsHBoxUtils<d,T>::Unique(container);
}

template <short_t d, class T>
bool gsHBox<d, T>::good() const
{
    return (m_indices.first.array() >=0).any() && (m_indices.second.array() >=0).any();
}

template <short_t d, class T>
void gsHBox<d, T>::clean(Container & container) const
{
    std::function<bool(const gsHBox<d,T> &)> pred = [](const gsHBox<d,T> & box) { return !(box.good()); };
    cIterator beg = std::remove_if(container.begin(),container.end(),pred);
    container.erase(beg,container.end());
}

} // namespace gismo
