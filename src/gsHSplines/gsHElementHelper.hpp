/** @file gsHElementHelper.hpp

    @brief Provides implementations for operations on hierarchical elements in a mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsHSplines/gsHElement.h>
#include <gsHSplines/gsTHBSplineBasis.h>

namespace gismo
{

    template <short_t d, class T>
    gsHElementHelper<d,T>::gsHElementHelper(const gsBasis<T> & basis)
    :
    gsHElementHelper(dynamic_cast<const gsHTensorBasis<d,T> &>(basis))
    {
    }

    template <short_t d, class T>
    gsHElementHelper<d,T>::gsHElementHelper(const gsHTensorBasis<d,T> & basis)
    :
    m_basis(basis)
    {
        // We make sure that the basis is not using manual levels, since we don't use lines like
        // `_knotIndexToDiadicIndex` in this code.
        // TODO: fix that
        GISMO_ASSERT(!basis.manualLevels(),"The basis must not have manual levels for this helper to work.");
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::element_t gsHElementHelper<d,T>::toElement(const gsVector<T,d> & low, const gsVector<T,d> & upp, typename gsHElementHelper<d,T>::level_t level, index_t patch) const
    {
        GISMO_ASSERT((low.array() <= upp.array()).all(), "Lower corner must be less than or equal to upper corner.");
        gsVector<index_t, d> lowIdx, uppIdx;
        for(index_t j = 0; j < d; j++)
        {
            // Convert the parameter coordinates to (unique) knot indices
            const gsKnotVector<T> & kv = m_basis.tensorLevel(level).knots(j);
            lowIdx(j) = (std::upper_bound(kv.domainUBegin(), kv.domainUEnd(),low[j] ) - 1).uIndex();
            uppIdx(j) = (std::upper_bound(kv.domainUBegin(), kv.domainUEnd()+1,upp[j] ) - 1).uIndex();
        }
        return element_t(lowIdx, uppIdx, level, patch);
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::HElementContainer gsHElementHelper<d,T>::toElements(index_t patch) const
    {
        HElementContainer elements;
        level_t level;
        for (const auto & elem : m_basis.domain()->allElements())
        {
            level = static_cast<const gsHDomainIterator<T,d> *>(&elem)->getLevel();
            // Convert knot coordinates to indices
            elements.insert(this->toElement(elem.lowerCorner(), elem.upperCorner(), level, patch));
        }
        return elements;
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::HElementContainer gsHElementHelper<d,T>::toElements(const ElementContainer & elements, level_t level) const
    {
        HElementContainer hElements;
        for (const auto & elem : elements)
        {
            hElements.emplace(elem, level);
        }
        return hElements;
    }

    template <short_t d, class T>
    std::pair<gsMatrix<T>,typename gsHElementHelper<d,T>::level_t> gsHElementHelper<d,T>::toBoxAndLevel(const typename gsHElementHelper<d,T>:: element_t & element) const
    {
        gsMatrix<T> box(d, 2);
        const point_t & low = element.lowerCorner();
        const point_t & upp = element.upperCorner();
        for(index_t j = 0; j < d; j++)
        {
            // Convert the knot indices back to parameter coordinates
            const gsKnotVector<T> & kv = m_basis.tensorLevel(element.level()).knots(j);
            box(j,0) = kv.uValue(low(j));
            box(j,1) = kv.uValue(upp(j));
        }
        return {box, element.level()};
    }

    template <short_t d, class T>
    gsMatrix<T> gsHElementHelper<d,T>::toBox(const element_t & element) const
    {
        return this->toBoxAndLevel(element).first;
    }

    template <short_t d, class T>
    std::pair<gsMatrix<T>,gsVector<typename gsHElementHelper<d,T>::level_t>> gsHElementHelper<d,T>::toBoxesAndLevels(const HElementContainer & elements) const
    {
        gsMatrix<T> boxes(d, 2*elements.size());
        gsVector<level_t> levels(elements.size());
        index_t i = 0;
        for (const auto & elem : elements)
        {
            std::pair<gsMatrix<T>, level_t> boxAndLevel = this->toBoxAndLevel(elem);
            boxes.middleCols(2*i,2) = boxAndLevel.first;
            levels(i) = boxAndLevel.second;
            i++;
        }
        return {boxes, levels};
    }

    template <short_t d, class T>
    gsMatrix<T> gsHElementHelper<d,T>::toBoxes(const HElementContainer & elements) const
    {
        return this->toBoxesAndLevels(elements).first;
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::level_t gsHElementHelper<d,T>::levelInBasis(const element_t & element) const
    {
        // // Get the level of the basis at the center of the element
        // gsVector<T,d> center = element.lowerCorner() + 0.5 * (element.upperCorner() - element.lowerCorner());
        // return m_basis.getLevelAtPoint(center);

        // Get the maximum level in the basis
        const index_t maxLevel = m_basis.maxLevel();
        // Convert the corners to the maximum level in the basis
        point_t low = element.lowerCorner() * (index_t)math::pow(2, maxLevel - element.level());
        point_t upp = element.upperCorner() * (index_t)math::pow(2, maxLevel - element.level());
        point_t center = (low + upp) / 2;
        return m_basis.getLevelAtIndex(center);
    }

    template <short_t d, class T>
    bool gsHElementHelper<d,T>::isActive(const element_t & element) const
    {
        // An element is considered active if its level in the basis matches its level in the hierarchy
        return this->levelInBasis(element) == element.level();
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::element_t gsHElementHelper<d,T>::getParent(const element_t & element) const
    {
        // GISMO_ASSERT(element.level() > 0, "Cannot get parent of an element at level 0.");
        const level_t parentLevel = element.level() - 1;
        point_t low, upp;
        for (index_t i = 0; i!=d; i++)
        {
            low.at(i) = element.lowerCorner().at(i) / 2;
            upp.at(i) = element.upperCorner().at(i) / 2 + (index_t)(element.upperCorner().at(i) % 2 != 0);
        }
        return element_t(low, upp, parentLevel, element.patch());
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::element_t gsHElementHelper<d,T>::getAncestor(const element_t & element, level_t jump) const
    {
        const level_t lvl = element.level();
        GISMO_ASSERT(jump >= 0, "Jump must be non-negative.");
        GISMO_ASSERT(lvl >= 0, "Cannot get ancestor of an element at level 0.");
        element_t parent = this->getParent(element);
        if (jump < lvl - 1)
        {
            element_t ancestor = this->getAncestor(parent, jump);
            return ancestor;
        }
        else if (jump == lvl - 1)
        {
            return parent;
        }
        else if (jump == lvl)
        {
            return element; // the element itself
        }
        else
        {
            GISMO_ERROR("Jump is larger than the level of the element.");
        }
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::HElementContainer gsHElementHelper<d,T>::getChildren(const element_t & element) const
    {
        const level_t childLevel = element.level() + 1;
        point_t low = element.lowerCorner();
        point_t upp = element.upperCorner();
        low.array() *= 2;
        upp.array() *= 2;
        HElementContainer children = this->explode(low, upp, childLevel, element.patch());
        return children;
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::HElementContainer gsHElementHelper<d,T>::getDescendants(const element_t & element, level_t jump) const
    {
        level_t lvl = element.level();
        // GISMO_ASSERT(lvl > k,"Current level should be larger than requested level, l = "<<lvl<<", k = "<<k);
        GISMO_ASSERT(lvl >= 0,"Level lvl = "<<lvl<<" should be larger than 0");
        HElementContainer descendants;
        if (jump == lvl)
        {
            descendants.insert(element);
        }
        else if (jump == lvl + 1)
        {
            descendants = this->getChildren(element);
        }
        else
        {
            HElementContainer tmp = this->getChildren(element); // children have level lvl+1
            for (level_t k_tmp = lvl + 1; k_tmp != jump; k_tmp++)
            {
                HElementContainer tmp2;
                for (const auto & elem : tmp)
                {
                    HElementContainer children = this->getChildren(elem);
                    tmp2.insert(children.begin(), children.end());
                }
                tmp = tmp2;
            }
            descendants = tmp;
        }

        return descendants;
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::HElementContainer gsHElementHelper<d,T>::getSiblings(const element_t & element, bool erase_original) const
    {
        GISMO_ASSERT(element.level() > 0, "Cannot get siblings of an element at level 0.");
        // Get the parent of the element
        element_t parent = this->getParent(element);
        // Get the children of the parent
        HElementContainer siblings = this->getChildren(parent);
        // Remove the element itself from the siblings
        if (erase_original)
            siblings.erase(element);
        return siblings;
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::box_t gsHElementHelper<d,T>::getSupportExtension(const element_t & element) const
    {
        gsElementHelper<d,T> helper(m_basis.tensorLevel(element.level()));
        return helper.getSupportExtension(element);
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::box_t gsHElementHelper<d,T>::getMultiLevelSupportExtension(const element_t & element, level_t jump) const
    {
        level_t lvl = element.level();
        GISMO_ASSERT(jump >= 0, "Jump must be non-negative.");
        GISMO_ASSERT(lvl >= 0, "Cannot get support extension of an element at level 0.");
        GISMO_ASSERT(jump <= lvl, "Jump cannot be larger than the level of the element.");

        if (jump == lvl)
        {
            return this->getSupportExtension(element);
        }
        else
        {
            // Get the ancestor at the specified jump level
            element_t ancestor = this->getAncestor(element, jump);
            // Get the support extension of the ancestor
            return this->getSupportExtension(ancestor);
        }
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::HElementContainer gsHElementHelper<d,T>::getHNeighborhood(const element_t & element, level_t m) const
    {
        HElementContainer extension;
        HElementContainer neighborhood;

        level_t lvl = element.level();
        level_t k = lvl - m + 1;
        if (k >= 0)
        {
            // Get multi level support extension on level k
            extension = this->explode(this->getMultiLevelSupportExtension(element, k), k);
            // Eliminate elements which are too low
            for (const auto & elem : extension)
            {
                if (this->isActive(elem))
                {
                    neighborhood.insert(elem);
                }
            }
        }
        return neighborhood;
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::HElementContainer gsHElementHelper<d,T>::getTNeighborhood(const element_t & element, level_t m) const
    {
        // Everything is in the same level so we can use normal HElementContainer
        HElementContainer neighborhood;
        HElementContainer extension;

        level_t lvl = element.level();
        level_t k = lvl - m + 2;
        if (k - 1 >= 0)
        {
            // Get multi-level support extension on level k
            extension = this->explode(this->getMultiLevelSupportExtension(element, k),k);
            for (const auto & elem : extension)
            {
                element_t parent = this->getParent(elem);
                if (this->isActive(parent))
                    neighborhood.insert(parent);
            }
        }
        return neighborhood;
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::HElementContainer gsHElementHelper<d,T>::getNeighborhood(const element_t & element, level_t m) const
    {
        if (dynamic_cast<const gsTHBSplineBasis<d,T> *>(&m_basis))
        {
            // If the basis is a tensor basis, use the TNeighborhood
            return this->getTNeighborhood(element, m);
        }
        else if (dynamic_cast<const gsHBSplineBasis<d,T> *>(&m_basis))
        {
            // If the basis is a HBSpline basis, use the TNeighborhood
            return this->getHNeighborhood(element, m);
        }
        else
        {
            GISMO_ERROR("The basis is neither a tensor basis nor a HBSpline basis. Cannot get neighborhood.");
        }
        return HElementContainer(); // Return an empty container if the basis is neither
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::HElementContainer gsHElementHelper<d,T>::getCextension(const element_t & element, level_t m) const
    {
        // See Definition 3.5 from [1]
        // [1] Carraturo, M., Giannelli, C., Reali, A. & Vázquez, R.
        // Suitably graded THB-spline refinement and coarsening: Towards an adaptive isogeometric analysis of additive manufacturing processes.
        // Comput. Methods Appl. Mech. Eng. 348, 660–679 (2019).
        //
        // This function returns the coarsening neighborhood of \a this, which is an element to be reactivated
        HElementContainer extension, children, childExtension, descendants, result;
        level_t targetLvl = element.level() + m;

        children = this->getChildren(element);
        for (const auto & child : children)
        {
            // Get the support extension of the child
            childExtension = this->explode(this->getSupportExtension(child),child.level());
            // Add the child extension to the extension
            extension.insert(childExtension.begin(), childExtension.end());
        }

        // All elements in the support extension of the children (so at level l+1)
        extension.insert(children.begin(), children.end());
        // question: needed MINUS the siblings??
        // extension = gsHBoxUtils<d,T>::Difference(extension,children);

        // result = extension;
        // gsDebugVar(targetLvl);
        for (const auto & eIt : extension)
        {
            descendants = this->getDescendants(eIt, targetLvl);
            // Add the descendants to the result
            result.insert(descendants.begin(), descendants.end());
        }
        return result;
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::HElementContainer gsHElementHelper<d,T>::getCneighborhood(const element_t & element, level_t m) const
    {
        // See Definition 3.5 from [1]
        // [1] Carraturo, M., Giannelli, C., Reali, A. & Vázquez, R.
        // Suitably graded THB-spline refinement and coarsening: Towards an adaptive isogeometric analysis of additive manufacturing processes.
        // Comput. Methods Appl. Mech. Eng. 348, 660–679 (2019).
        //
        // This function returns the coarsening neighborhood of \a this, which is an element to be reactivated

        HElementContainer descendants, result;
        descendants = this->getCextension(element, m);

        for (const auto & elem : descendants)
        {
            if (this->levelInBasis(elem) >= elem.level()) // the level is even larger (i.e. even higher decendant)!!
            {
                result.insert(elem);
            }
        }
        return result;
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::HElementContainer gsHElementHelper<d,T>::markRecursive(const HElementContainer & elements, level_t level, level_t m) const
    {
        HElementContainer marked_elements = elements;
        HElementContainer marked_level;
        HElementContainer marked_k;

        // Collect the elements at the current level
        for (const auto & elem : elements)
        {
            if (elem.level() == level)
                marked_level.insert(elem);
        }

        HElementContainer neighbors, tmp_neighbors;
        for (const auto & elem : marked_level)
        {
            tmp_neighbors = this->getNeighborhood(elem, m);
            neighbors.insert(tmp_neighbors.begin(), tmp_neighbors.end());
        }

        level_t k = level - m + 1;
        if (!neighbors.empty() && k >= 0)
        {
            // Collect the elements at level k
            for (const auto & elem : marked_elements)
            {
                if (elem.level() == k)
                    marked_k.insert(elem);
            }

            // Add the elements of the neighborhood which are at level k
            for (const auto & neighbor : neighbors)
            {
                if (neighbor.level() == k)
                {
                    marked_k.insert(neighbor);
                }
            }

            // Add the marked elements at level k to the marked elements
            marked_elements.insert(marked_k.begin(), marked_k.end());
            marked_elements = this->markRecursive(marked_elements, k, m);
        }
        return marked_elements;
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::HElementContainer gsHElementHelper<d,T>::markAdmissible(const HElementContainer & marked, level_t m) const
    {
        // Get the maximum level of the marked elements
        level_t maxLevel = 0;
        for (const auto & elem : marked)
        {
            if (elem.level() > maxLevel)
                maxLevel = elem.level();
        }
        // Mark recursively
        HElementContainer marked_elements = marked;
        for (level_t l = 0; l!=maxLevel+1; l++)
        {
            marked_elements = this->markRecursive(marked_elements, l, m);
        }
        return marked_elements;
    }

    template <short_t d, class T>
    std::vector<index_t> gsHElementHelper<d,T>::toRefBox(const element_t & element, level_t targetLevel) const
    {
        GISMO_ASSERT(targetLevel > element.level(),
                     "Target level must be larger than the element level. "
                     "Element level: " << element.level() << ", target level: " << targetLevel);
        std::vector<index_t> result(2*d + 1);
        result[0] = targetLevel; // Store the target level at the first position
        level_t diff = targetLevel - element.level();
        index_t lowerIndex, upperIndex, degree;

        // Get the unique knot index for the lower and upper corners of the element box
        gsElementHelper<d,T> helper(m_basis.tensorLevel(targetLevel));
        for (index_t i = 0; i < d; ++i)
        {
            // Get the degree
            degree = m_basis.degree(i);
            lowerIndex = element.lowerCorner()(i)*math::pow(2, diff);
            upperIndex = element.upperCorner()(i)*math::pow(2, diff);
            if (degree % 2 == 1 && degree > 1)
                ( (lowerIndex < (degree-1)/2-1) ? lowerIndex = 0 : lowerIndex -= (degree-1)/2-1);
            else
                ( (lowerIndex < (degree-1)/2)   ? lowerIndex = 0 : lowerIndex -= (degree-1)/2  );

            result[i+1] = lowerIndex;
            result[d+i+1] = upperIndex;
        }

        return result;
    }

    template <short_t d, class T>
    std::vector<index_t> gsHElementHelper<d,T>::toRefBox(const element_t & element) const
    {
        return this->toRefBox(element, element.level() + 1);
    }

    template <short_t d, class T>
    std::vector<index_t> gsHElementHelper<d,T>::toRefBoxes(const HElementContainer & elements) const
    {
        std::vector<index_t> refBoxes;
        refBoxes.reserve(elements.size() * (2 * d + 1));
        for (const auto & elem : elements)
        {
            std::vector<index_t> refBox = this->toRefBox(elem);
            refBoxes.insert(refBoxes.end(), refBox.begin(), refBox.end());
        }
        return refBoxes;
    }

    template <short_t d, class T>
    std::vector<index_t> gsHElementHelper<d,T>::toCrsBox(const element_t & element, level_t targetLevel) const
    {
        GISMO_ASSERT(targetLevel < element.level(),
                     "Target level must be smaller than the element level. "
                     "Element level: " << element.level() << ", target level: " << targetLevel);
        std::vector<index_t> result(2*d + 1);
        result[0] = targetLevel; // Store the target level at the first position
        level_t diff = element.level() - targetLevel;
        for (index_t i = 0; i < d; ++i)
        {
            // floor division for lower corner
            result[i+1] = element.lowerCorner()(i)/math::pow(2, diff); // floor division
            // ceil division for upper corner: (x + (2^diff-1)) / 2^diff
            result[d+i+1] = (element.upperCorner()(i) + math::pow(2, diff) - 1) / math::pow(2, diff);
        }
        return result;
    }

    template <short_t d, class T>
    std::vector<index_t> gsHElementHelper<d,T>::toCrsBox(const element_t & element) const
    {
        return this->toCrsBox(element, element.level() - 1);
    }

    template <short_t d, class T>
    std::vector<index_t> gsHElementHelper<d,T>::toCrsBoxes(const HElementContainer & elements) const
    {
        std::vector<index_t> refBoxes;
        refBoxes.reserve(elements.size() * (2 * d + 1));
        for (const auto & elem : elements)
        {
            std::vector<index_t> refBox = this->toCrsBox(elem);
            refBoxes.insert(refBoxes.end(), refBox.begin(), refBox.end());
        }
        return refBoxes;
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::HElementContainer gsHElementHelper<d,T>::explode(const point_t & low, const point_t & upp, level_t level, index_t patch) const
    {
        GISMO_ASSERT(level >= 0, "Level must be non-negative.");
        gsElementHelper<d,T> helper(m_basis.tensorLevel(level));
        ElementContainer exploded = helper.explode(low, upp, patch);
        return this->toElements(exploded, level);
    }

    template <short_t d, class T>
    typename gsHElementHelper<d,T>::HElementContainer gsHElementHelper<d,T>::explode(const box_t & box, level_t level, index_t patch) const
    {
        return this->explode(box.col(0), box.col(1), level, patch);
    }

    template <short_t d, class T>
    bool gsHElementHelper<d,T>::contains(const element_t & element1, const element_t & element2) const
    {
        // Check if element1 contains element2
        bool result = true;
        result &= element1.level() <= element2.level(); // element1 must be at the same or coarser level
        if (result)
        {
            element_t ancestor = this->getAncestor(element2, element1.level());
            // Check if the indices of the ancestor are contained in element1
            result &= (ancestor.lowerCorner().array() >= element1.lowerCorner().array()).all();
            result &= (ancestor.upperCorner().array() <= element1.upperCorner().array()).all();
        }
        return result;
    }

}; //namespace gismo
