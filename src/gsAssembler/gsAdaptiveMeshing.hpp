/** @file gsAdaptiveMeshing.hpp

    @brief Provides class for adaptive refinement.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-)
*/

#pragma once


#include <gsHSplines/gsHBSplineBasis.h>
#include <gsHSplines/gsTHBSplineBasis.h>

namespace gismo
{

template <class T>
gsAdaptiveMeshing<T>::gsAdaptiveMeshing()
{
    defaultOptions();
}

template <class T>
gsAdaptiveMeshing<T>::gsAdaptiveMeshing(gsFunctionSet<T> & input)
:
m_input(&input)
{
    defaultOptions();
    rebuild();
}

template <class T>
void gsAdaptiveMeshing<T>::rebuild()
{
    getOptions();

    m_indices.clear();
    m_boxes.clear();
    this->_makeMap(m_input,m_indices,m_boxes);

    bool check = true;
    for (typename indexMapType::iterator it=m_indices.begin(); it!=m_indices.end(); it++)
        check &= gsHBoxEqual<2,T>()(it->first,*m_boxes[it->second]);

    for (typename boxMapType::iterator it=m_boxes.begin(); it!=m_boxes.end(); it++)
        check &= it->first==m_indices[*it->second];

    GISMO_ASSERT(check,"Something went wrong in the construction of the mappers");
}

template <class T>
void gsAdaptiveMeshing<T>::_makeMap(const gsFunctionSet<T> * input, typename gsAdaptiveMeshing<T>::indexMapType & indexMap, typename gsAdaptiveMeshing<T>::boxMapType & boxMap)
{
    // typename gsAdaptiveMeshing<T>::indexMapType indexMap;
    // typename gsAdaptiveMeshing<T>::boxMapType   boxMap;

    const gsBasis<T> * basis = nullptr;

    // Make the container of boxes
    // #pragma omp parallel
    // {
// #ifdef _OPENMP
//         const int tid = omp_get_thread_num();
//         const int nt  = omp_get_num_threads();
//         index_t patch_cnt = 0;
// #endif

        index_t c = 0;
        for (index_t patchInd=0; patchInd < input->nPieces(); ++patchInd)
        {
            // Initialize domain element iterator
            if ( const gsMultiPatch<T> * mp = dynamic_cast<const gsMultiPatch<T>*>(input) ) basis = &(mp->basis(patchInd));
            if ( const gsMultiBasis<T> * mb = dynamic_cast<const gsMultiBasis<T>*>(input) ) basis = &(mb->basis(patchInd));
            GISMO_ENSURE(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");
            // for all elements in patch pn
            typename gsBasis<T>::domainIter domIt = basis->makeDomainIterator();
            gsHDomainIterator<T,2> * domHIt = nullptr;
            domHIt = dynamic_cast<gsHDomainIterator<T,2> *>(domIt.get());
            GISMO_ENSURE(domHIt!=nullptr,"Domain not loaded");

// #ifdef _OPENMP
//             c = patch_cnt + tid;
//             patch_cnt += domHIt->numElements();// a bit costy
//             for ( domHIt->next(tid); domHIt->good(); domHIt->next(nt) )
// #else
            for (; domHIt->good(); domHIt->next())
// #endif
            {
                // #pragma omp critical (gsAdaptiveMeshingmakeBoxesinsert1)
                {
                    HBox box(domHIt,patchInd);
                    std::pair<typename gsAdaptiveMeshing<T>::indexMapType::iterator,bool> mapIt = indexMap.insert({box,c});
                    if (mapIt.second)
                    {
                        // std::pair<typename gsAdaptiveMeshing<T>::boxMapType::iterator,bool> indexIt =
                        boxMap.insert({c,const_cast<gsHBox<2,real_t> *>(&(mapIt.first->first))});
    // #                   ifdef _OPENMP
    //                     c += nt;
    // #                   else
                        c++;
    // #                   endif

                    }
                }
            }
        }
    // }

    // return std::make_pair(give(indexMap),give(boxMap));
}



template <class T>
void gsAdaptiveMeshing<T>::_assignErrors(boxMapType & container, const std::vector<T> & elError)
{
    GISMO_ASSERT(elError.size()==container.size(),"The number of errors must be the same as the number of elements, but "<<elError.size()<<"!="<<container.size());
    index_t k=0;

    if (m_refRule == PBULK || m_crsRule == PBULK)
        for (typename boxMapType::iterator it = container.begin(); it!=container.end(); it++, k++)
            it->second->setAndProjectError(elError[k],m_alpha,m_beta);
    else
        for (typename boxMapType::iterator it = container.begin(); it!=container.end(); it++, k++)
            it->second->setError(elError[k]);

    m_totalError = _totalError(m_boxes);
    m_maxError = _maxError(m_boxes);
    if (m_alpha!=-1 && m_beta!=-1)
    {
        const gsBasis<T> * basis = nullptr;
        index_t patchInd=0; // Assymes now that the minDegree is lowest on patch 0, or that the degree there is representative
        if ( const gsMultiPatch<T> * mp = dynamic_cast<const gsMultiPatch<T>*>(m_input) ) basis = &(mp->basis(patchInd));
        if ( const gsMultiBasis<T> * mb = dynamic_cast<const gsMultiBasis<T>*>(m_input) ) basis = &(mb->basis(patchInd));
        GISMO_ASSERT(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");

        m_uniformRefError = math::pow(1/2.,m_alpha * basis->minDegree() + m_beta) * m_totalError;
        m_uniformCrsError = math::pow(  2.,m_alpha * basis->minDegree() + m_beta) * m_totalError;
    }
    else
    {
        m_uniformRefError = 0;
        m_uniformCrsError = 0;
    }
}

/** \brief Marks elements/cells for refinement.
 *
 * Let the global error/error estimate \f$\eta\f$ be a sum of element/cell-wise
 * local contributions:
 * \f[ \eta = \sum_{K} \eta_k \quad \mathrm{or} \quad \eta^2 = \sum_K \eta_K^2 \f]
 *
 * Computes a threshold \f$\Theta\f$ and marks all elements \f$K\f$ for refinement,
 * for which
 * \f[ \eta_K \geq \Theta \f]
 * holds.
 * Three criteria for computing \f$\Theta\f$ are currently (26.Nov.2014) implemented:
 *
 * Let \f$\rho\f$ denote the m_input parameter \em parameter.
 *
 * <b>refCriterion = 1 = treshold, GARU-criterion</b> (greatest appearing eRror utilization):\n
 * Threshold computed based on the largest of all appearing local errors:
 * \f[ \Theta = \rho \cdot \max_K \{ \eta_K \} \f]
 * The actual number of marked elements can vary in each refinement step,
 * depending on the distribution of the error.
 *
 * <b>refCriterion = 2 = cellPercentage, PUCA-criterion</b> (percentile-utilizing cutoff ascertainment):\n
 * In each step, a certain percentage of all elements are marked.
 * \f[ \Theta = (1-\rho)\cdot 100\ \textrm{-percentile of}\ \{ \eta_K \}_K \f]
 * For example, if \f$\rho = 0.8\f$, those 20% of all elements which have the
 * largest local errors are marked for refinement.
 *
 * <b>refCriterion = 3 = errorFraction, BULK-criterion</b> ("Doerfler-marking"):\n
 * The threshold is chosen in such a manner that the local
 * errors on the marked cells sum up to a certain fraction of the
 * global error:
 * \f[ \sum_{ K:\ \eta_K \geq \Theta } \eta_K \geq (1-\rho) \cdot \eta \f]
 *
 * \param elError std::vector of local errors on some elements.
 * \param refCriterion selects the criterion (see above) for marking elements.
 * \param parameter parameter \f$ \rho \f$ for refinement criterion (see above).\n
 * \f$\rho = 0\f$ corresponds to global refinement,\n
 * \f$ \rho=1\f$ corresponds to (almost) no refinement.
 * \param[out] elMarked std::vector of Booleans indicating whether the corresponding element is marked or not.
 *
 * \ingroup Assembler
 */

template <class T>
template<bool _coarsen,bool _admissible>
void gsAdaptiveMeshing<T>::_markElements(  const std::vector<T> & elError, const index_t refCriterion, const std::vector<gsHBoxCheck<2,T> *> & predicates, HBoxContainer & elMarked) const
{
    // Mark using different rules
    switch (refCriterion)
    {
    case GARU:
        _markThreshold<_coarsen,_admissible>(m_boxes,predicates,elMarked);
        break;
    case PUCA:
        _markPercentage<_coarsen,_admissible>(m_boxes,predicates,elMarked);
        break;
    case BULK:
        _markFraction<_coarsen,_admissible>(m_boxes,predicates,elMarked);
        break;
    case PBULK:
        _markProjectedFraction<_coarsen,_admissible>(m_boxes,predicates,elMarked);
        break;
    default:
        GISMO_ERROR("unknown marking strategy");
    }
}

template <class T>
void gsAdaptiveMeshing<T>::_crsPredicates_into(std::vector<gsHBoxCheck<2,T> *> & predicates)
{
    HBoxContainer empty;
    predicates.push_back(new gsMinLvlCompare<2,T>(0));
    predicates.push_back(new gsOverlapCompare<2,T>(empty,m_m));
}

template <class T>
void gsAdaptiveMeshing<T>::_crsPredicates_into(const HBoxContainer & markedRef, std::vector<gsHBoxCheck<2,T> *> & predicates)
{
    predicates.push_back(new gsMinLvlCompare<2,T>(0));
    predicates.push_back(new gsOverlapCompare<2,T>(markedRef,m_m));
}

template <class T>
void gsAdaptiveMeshing<T>::_refPredicates_into(std::vector<gsHBoxCheck<2,T> *> & predicates)
{
    predicates.push_back(new gsMaxLvlCompare<2,T>(m_maxLvl));
}


template <class T>
std::vector<index_t> gsAdaptiveMeshing<T>::_sortPermutation( const boxMapType & container)
{
    std::vector<index_t> idx(container.size());
    std::iota(idx.begin(),idx.end(),0);
    std::stable_sort(idx.begin(), idx.end(),
           [&container](size_t i1, size_t i2) {
            return container.at(i1)->error() < container.at(i2)->error();});

    return idx;
}

template <class T>
std::vector<index_t> gsAdaptiveMeshing<T>::_sortPermutationProjectedRef( const boxMapType & container)
{
    std::vector<index_t> idx(container.size());
    std::iota(idx.begin(),idx.end(),0);
    std::stable_sort(idx.begin(), idx.end(),
           [&container](size_t i1, size_t i2) {
            return container.at(i1)->projectedImprovement() < container.at(i2)->projectedImprovement();});

    return idx;
}

template <class T>
std::vector<index_t> gsAdaptiveMeshing<T>::_sortPermutationProjectedCrs( const boxMapType & container)
{
    std::vector<index_t> idx(container.size());
    std::iota(idx.begin(),idx.end(),0);
    std::stable_sort(idx.begin(), idx.end(),
           [&container](size_t i1, size_t i2) {
            return container.at(i1)->projectedSetBack() < container.at(i2)->projectedSetBack();});

    return idx;
}

// template <class T>
// void gsAdaptiveMeshing<T>::_sortPermutated( const std::vector<index_t> & permutation, boxContainer & container)
// {
//     GISMO_ASSERT(permutation.size()==container.size(),"Permutation and vector should have the same size, but "<<permutation.size()<<"!="<<container.size());
//     boxContainer sorted(container.size());
//     for (size_t k=0; k!=container.size(); k++)
//         sorted.at(k) = container.at(permutation.at(k));

//     container.swap(sorted);
// }

template <class T>
typename gsAdaptiveMeshing<T>::HBox * gsAdaptiveMeshing<T>::_boxPtr(const HBox & box) const
{
    // Fails if box is not in m_boxes
    HBox * boxPtr = m_boxes.at(m_indices.at(box));
    return boxPtr;
    // return m_boxes[m_indices[box]];
}

template <class T>
bool gsAdaptiveMeshing<T>::_checkBox( const HBox & box, const std::vector<gsHBoxCheck<2,T> *> predicates) const
{
    bool check = true;
    for (typename std::vector<gsHBoxCheck<2,T>*>::const_iterator errIt = predicates.begin(); errIt!=predicates.end(); errIt++)
        check &= (*errIt)->check(box);
    return check;
}

template <class T>
bool gsAdaptiveMeshing<T>::_checkBoxes( const typename HBox::Container & boxes, const std::vector<gsHBoxCheck<2,T> *> predicates) const
{
    bool check = true;
    for (typename HBox::cIterator it = boxes.begin(); it!=boxes.end(); it++)
        check &= _checkBox(*it,predicates);

    return check;
}

// Use the mapTypes here!!


template <class T>
void gsAdaptiveMeshing<T>::_addAndMark( HBox & box, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    _boxPtr(box)->mark();
    elMarked.add(*_boxPtr(box));
}

template <class T>
void gsAdaptiveMeshing<T>::_addAndMark( typename HBox::Container & boxes, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    for (typename HBox::Iterator it = boxes.begin(); it!=boxes.end(); it++)
        _addAndMark(*it,elMarked);
}

template <class T>
void gsAdaptiveMeshing<T>::_setContainerProperties( typename HBox::Container & boxes ) const
{
    for (typename HBox::Iterator it = boxes.begin(); it!=boxes.end(); it++)
    {
        // HBox * box = _boxPtr(*it);
        // it->setError(box->error());
        // it->setMark(box->marked());
        it->setAndProjectError(_boxPtr(*it)->error(),m_alpha,m_beta);
    }
}

template <class T>
T gsAdaptiveMeshing<T>::_totalError(const boxMapType & elements)
{
    // get total error
    // Accumulation operator for boxMapType
    auto accumulate_error_ptr = [](const T & val, const typename boxMapType::value_type & b)
    { return val + b.second->error(); };
    T totalError = std::accumulate(elements.begin(),elements.end(),0.0,accumulate_error_ptr);
    return totalError;
}

// Refinement: parameter % contributions to the total error of the highest cells are marked
// Coarsening: parameter % contributions to the total error of the lowest cells are marked
// In both cases, the total contribution is always lower than the threshold, i.e. if an element causes an exceed of the error, it is not taken into account

template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if< _coarsen &&  _admissible, void>::type
gsAdaptiveMeshing<T>::_markFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Fraction marking for coarsening...\n";
    T cummulErrMarked = T(0);
    T errorMarkSum =  m_crsParam * m_totalError;

    // Accumulation operator for neighborhoods (boxContainer)
    auto accumulate_error = [](const T & val, const typename boxContainer::value_type & b)
    { return val + b.error(); };

    auto loop_action = [this,&elements,&predicates,&cummulErrMarked,&errorMarkSum,&elMarked,&accumulate_error]
                    (const index_t & index)
    {
        // Get the box
        HBox_ptr box = elements.at(index);
        // Continue if the box is already marked or if it does not satsfy all checks
        if (_boxPtr(*box)->marked() || !_checkBox(*box,predicates))
            return false;

        // Get the neighborhoods
        typename HBox::Container sibs = box->getSiblings();
        _setContainerProperties(sibs);
        if (!gsHBoxUtils<2,T>::allActive(sibs))
            return false;

        HBoxContainer siblings(sibs);
        T siblingError = std::accumulate(sibs.begin(),sibs.end(),0.0,accumulate_error);

        // Check all siblings if they satisfy the predicates
        if (_checkBoxes(sibs,predicates))
        {
            cummulErrMarked += siblingError;
            _addAndMark(sibs,elMarked);
            cummulErrMarked += box->error();
            _addAndMark(*box,elMarked);
        }
        // return false;
        return (cummulErrMarked > errorMarkSum);
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_crsPermutation.cbegin(),m_crsPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark fraction] Marked "<<elMarked.totalSize()<<" elements with marked error = "<<cummulErrMarked<<" and threshold = "<<errorMarkSum<<((it==m_crsPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}


template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if< _coarsen && !_admissible, void>::type
gsAdaptiveMeshing<T>::_markFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Fraction marking for coarsening...\n";
    T cummulErrMarked = T(0);
    T errorMarkSum = m_crsParam * m_totalError;

    auto loop_action = [this,&elements,&predicates,&cummulErrMarked,&errorMarkSum,&elMarked]
                    (const index_t & index)
    {
        // Get the box
        HBox_ptr box = elements.at(index);

        // Check box if it satisfy the predicates
        if (_checkBox(*box,predicates))
        {
            cummulErrMarked += box->error();
            _addAndMark(*box,elMarked);
        }
        // return false;
        return (cummulErrMarked > errorMarkSum);
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_crsPermutation.cbegin(),m_crsPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark fraction] Marked "<<elMarked.totalSize()<<" elements with marked error = "<<cummulErrMarked<<" and threshold = "<<errorMarkSum<<((it==m_crsPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}

template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if<!_coarsen &&  _admissible, void>::type
gsAdaptiveMeshing<T>::_markFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Fraction marking (admissible) for refinement...\n";
    T cummulErrMarked = T(0);
    T errorMarkSum = m_refParam * m_totalError;

    // Accumulation operator for neighborhoods (boxContainer)
    auto accumulate_error = [](const T & val, const typename boxContainer::value_type & b)
    { return val + b.error(); };

    auto loop_action = [this,&elements,&predicates,&cummulErrMarked,&errorMarkSum,&elMarked,&accumulate_error]
                    (const index_t & index)
    {
        // Get the box
        HBox_ptr box = elements.at(index);
        if (_boxPtr(*box)->marked())
        {
            return false;
        }

        // Get the neighborhoods
        typename HBox::Container neighborhood = HBoxUtils::toContainer(HBoxUtils::markAdmissible(*box,m_m));
        _setContainerProperties(neighborhood);
        T neighborhoodError = std::accumulate(neighborhood.begin(),neighborhood.end(),0.0,accumulate_error);

        // if the errorMarkSum is exceeded with the current contribution, return true
        // NOTE: this can mean that this element suddenly has a large contribution! i.e., larger then the next element in lin
        // Maybe remove?
        // if (cummulErrMarked + neighborhoodError > errorMarkSum)
        //     return true;

        // Check all elements in the neighborhood if they satisfy the predicates
        if (_checkBoxes(neighborhood,predicates))
        {
            cummulErrMarked += neighborhoodError;
            _addAndMark(neighborhood,elMarked);
        }
        // return false;
        return (cummulErrMarked > errorMarkSum);
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_refPermutation.cbegin(),m_refPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark fraction] Marked "<<elMarked.totalSize()<<" elements with marked error = "<<cummulErrMarked<<" and threshold = "<<errorMarkSum<<((it==m_refPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}

template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if<!_coarsen && !_admissible, void>::type
gsAdaptiveMeshing<T>::_markFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Fraction marking (not admissible) for refinement...\n";
    T cummulErrMarked = T(0);
    T errorMarkSum = m_refParam * m_totalError;

    auto loop_action = [this,&elements,&predicates,&cummulErrMarked,&errorMarkSum,&elMarked]
                    (const index_t & index)
    {
        HBox_ptr box = elements.at(index);

        // if the errorMarkSum is exceeded with the current contribution, return true
        // if (cummulErrMarked + box->error() > errorMarkSum)
        //     return true;

        if (_checkBox(*box,predicates))
        {
            cummulErrMarked += box->error();
            _addAndMark(*box,elMarked);
        }
        return (cummulErrMarked > errorMarkSum);
        // return false;
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_refPermutation.cbegin(),m_refPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark fraction] Marked "<<elMarked.totalSize()<<" elements with marked error = "<<cummulErrMarked<<" and threshold = "<<errorMarkSum<<((it==m_refPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}

template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if< _coarsen &&  _admissible, void>::type
gsAdaptiveMeshing<T>::_markProjectedFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Projected fraction (admissible) marking for coarsening...\n";
    T targetError = m_crsParamExtra;
    if (m_totalError > targetError)
        return;
    if (m_uniformCrsError < targetError) // then all elements should be refined
    {
        this->_markFraction_impl<_coarsen,_admissible>(elements,predicates,elMarked);
        return;
    }
    T projectedError = m_totalError;

    // Accumulation operator for neighborhoods (boxContainer)
    auto accumulate_improvement = [](const T & val, const typename boxContainer::value_type & b)
    { return val + b.projectedErrorCrs() - b.error(); };

    auto loop_action = [this,&elements,&predicates,&projectedError,&targetError,&elMarked,&accumulate_improvement]
                    (const index_t & index)
    {
        // Get the box
        HBox_ptr box = elements.at(index);
        // Continue if the box is already marked or if it does not satsfy all checks
        if (_boxPtr(*box)->marked() || !_checkBox(*box,predicates))
            return false;

        // Get the neighborhoods
        typename HBox::Container sibs = box->getSiblings();
        _setContainerProperties(sibs);
        if (!gsHBoxUtils<2,T>::allActive(sibs))
            return false;

        HBoxContainer siblings(sibs);
        T siblingSetBack = std::accumulate(sibs.begin(),sibs.end(),0.0,accumulate_improvement);

        // Check all siblings if they satisfy the predicates
        if (_checkBoxes(sibs,predicates))
        {
            projectedError += siblingSetBack;
            _addAndMark(sibs,elMarked);
            projectedError += box->projectedErrorCrs() - box->error();
            _addAndMark(*box,elMarked);
        }

        // return false;
        return (projectedError > targetError);
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_crsPermutation.cbegin(),m_crsPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark projected fraction] Marked "<<elMarked.totalSize()<<" elements with projected error = "<<projectedError<<" and target error = "<<targetError<<((it==m_crsPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}


template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if< _coarsen && !_admissible, void>::type
gsAdaptiveMeshing<T>::_markProjectedFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    GISMO_NO_IMPLEMENTATION;
    // gsDebug<<"Projected fraction marking for coarsening...\n";
    // T projectedError = m_totalError;
    // T targetError = m_crsParam;
    // if (projectedError > targetError)
    //     return;

    // auto loop_action = [this,&elements,&predicates,&projectedError,&targetError,&elMarked]
    //                 (const index_t & index)
    // {
    //     // Get the box
    //     HBox_ptr box = elements.at(index);

    //     // Check box if it satisfy the predicates
    //     if (_checkBox(*box,predicates))
    //     {
    //         projectedError += box->projectedImprovement();
    //         _addAndMark(*box,elMarked);
    //     }
    //     // return false;
    //     return (projectedError > targetError);
    // };

    // std::find_if(m_crsPermutation.cbegin(),m_crsPermutation.cend(),loop_action);
    // elMarked = HBoxUtils::Unique(elMarked);
    // gsDebug<<"[Mark fraction] Marked "<<elMarked.totalSize()<<" elements with projected error = "<<projectedError<<" and target error = "<<targetError<<"\n";
}

template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if<!_coarsen &&  _admissible, void>::type
gsAdaptiveMeshing<T>::_markProjectedFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Projected fraction (admissible) marking for refinement...\n";
    T targetError = m_refParamExtra;
    if (m_totalError < targetError)
        return;
    if (m_uniformRefError > targetError) // then all elements should be refined
    {
        this->_markFraction_impl<_coarsen,_admissible>(elements,predicates,elMarked);
        return;
    }
    T projectedError = m_totalError;

    // Accumulation operator for neighborhoods (boxContainer)
    auto accumulate_improvement = [](const T & val, const typename boxContainer::value_type & b)
    {
        return val + b.error() - b.projectedErrorRef();
    };

    auto loop_action = [this,&elements,&predicates,&projectedError,&targetError,&elMarked,&accumulate_improvement]
                    (const index_t & index)
    {
        // Get the box
        HBox_ptr box = elements.at(index);
        if (_boxPtr(*box)->marked())
        {
            return false;
        }

        // Get the neighborhoods
        typename HBox::Container neighborhood = HBoxUtils::toContainer(HBoxUtils::markAdmissible(*box,m_m));
        HBoxContainer neighborhoodtmp = HBoxContainer(neighborhood);
        _setContainerProperties(neighborhood);
        T neighborhoodImprovement = std::accumulate(neighborhood.begin(),neighborhood.end(),0.0,accumulate_improvement);

        // Check all elements in the neighborhood if they satisfy the predicates
        if (_checkBoxes(neighborhood,predicates))
        {
            projectedError -= neighborhoodImprovement;
            _addAndMark(neighborhood,elMarked);
        }
        // return false;
        return (projectedError < targetError);
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_refPermutation.cbegin(),m_refPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark projected fraction] Marked "<<elMarked.totalSize()<<" elements with projected error = "<<projectedError<<" and target error = "<<targetError<<((it==m_refPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}

template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if<!_coarsen && !_admissible, void>::type
gsAdaptiveMeshing<T>::_markProjectedFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Projected fraction (not admissible) marking for refinement...\n";
    T targetError = m_refParamExtra;
    if (m_totalError < targetError)
        return;
    if (m_uniformRefError > targetError) // then all elements should be refined
    {
        this->_markFraction_impl<_coarsen,_admissible>(elements,predicates,elMarked);
        return;
    }
    T projectedError = m_totalError;

    auto loop_action = [this,&elements,&predicates,&projectedError,&targetError,&elMarked]
                    (const index_t & index)
    {
        HBox_ptr box = elements.at(index);

        // if the errorMarkSum is exceeded with the current contribution, return true
        // if (cummulErrMarked + box->error() > errorMarkSum)
        //     return true;

        if (_checkBox(*box,predicates))
        {
            projectedError -= box->projectedImprovement();
            _addAndMark(*box,elMarked);
        }
        return (projectedError < targetError);
        // return false;
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_refPermutation.cbegin(),m_refPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark projected fraction] Marked "<<elMarked.totalSize()<<" elements with projected error = "<<projectedError<<" and target error = "<<targetError<<((it==m_refPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}

template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if< _coarsen &&  _admissible, void>::type
gsAdaptiveMeshing<T>::_markPercentage_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Percentage marking for coarsening...\n";
    // Total number of elements:
    size_t NE = elements.size();
    // Compute the index from which the refinement should start,
    // once the vector is sorted.
    index_t NR = cast<T,index_t>( math::floor( m_crsParam * T(NE) ) );

    index_t nmarked = 0;
    auto loop_action = [this,&elements,&predicates,&NR,&nmarked,&elMarked]
                    (const index_t & index)
    {
        HBox * box = elements.at(index);

        // Continue if the box is already marked or if it does not satsfy all checks
        if (_boxPtr(*box)->marked() || !_checkBox(*box,predicates))
            return false;

        // Get the neighborhoods
        typename HBox::Container sibs = box->getSiblings();
        _setContainerProperties(sibs);
        if (!gsHBoxUtils<2,T>::allActive(sibs))
            return false;

        // Check all children if they satisfy the predicates
        if (_checkBoxes(sibs,predicates))
        {
            nmarked += sibs.size();
            _addAndMark(sibs,elMarked);
            nmarked += 1;
            _addAndMark(*box,elMarked);
        }
        return (nmarked > NR);
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_crsPermutation.cbegin(),m_crsPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark percentage] Marked "<<elMarked.totalSize()<<", ("<<nmarked<<") elements ("<<(T)nmarked/NE*100<<"%"<<" of NE "<<NE<<") and threshold = "<<NR<<" ("<<m_crsParam*100<<"%)"<<((it==m_crsPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}

template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if< _coarsen && !_admissible, void>::type
gsAdaptiveMeshing<T>::_markPercentage_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Percentage marking for coarsening...\n";
    // Total number of elements:
    size_t NE = elements.size();
    // Compute the index from which the refinement should start,
    // once the vector is sorted.
    index_t NR = cast<T,index_t>( math::floor( m_crsParam * T(NE) ) );

    index_t nmarked = 0;
    auto loop_action = [this,&elements,&predicates,&NR,&nmarked,&elMarked]
                    (const index_t & index)
    {
        HBox * box = elements.at(index);

        if (_checkBox(*box,predicates))
        {
            nmarked += 1;
            _addAndMark(*box,elMarked);
        }
        return (nmarked > NR);
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_crsPermutation.cbegin(),m_crsPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark percentage] Marked "<<elMarked.totalSize()<<", ("<<nmarked<<") elements ("<<(T)nmarked/NE*100<<"%"<<" of NE "<<NE<<") and threshold = "<<NR<<" ("<<m_crsParam*100<<"%)"<<((it==m_crsPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}

template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if<!_coarsen &&  _admissible, void>::type
gsAdaptiveMeshing<T>::_markPercentage_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Percentage marking (admissible) for refinement...\n";
    // Total number of elements:
    size_t NE = elements.size();
    // Compute the index from which the refinement should start,
    // once the vector is sorted.
    index_t NR = cast<T,index_t>( math::floor( m_refParam * T(NE) ) );

    index_t nmarked = 0;
    auto loop_action = [this,&elements,&predicates,&nmarked,&NR,&elMarked]
                    (const index_t & index)
    {
        // Get the box
                HBox_ptr box = elements.at(index);

        if (_boxPtr(*box)->marked())
            return false;

        typename HBox::Container neighborhood = HBoxUtils::toContainer(HBoxUtils::markAdmissible(*box,m_m));
        _setContainerProperties(neighborhood);
        // Check all elements in the neighborhood if they satisfy the predicates
        if (_checkBoxes(neighborhood,predicates))
        {
            nmarked += neighborhood.size();
            _addAndMark(neighborhood,elMarked);
        }
        return (nmarked > NR);
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_refPermutation.cbegin(),m_refPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark percentage] Marked "<<elMarked.totalSize()<<", ("<<nmarked<<") elements ("<<(T)nmarked/NE*100<<"%"<<" of NE "<<NE<<") and threshold = "<<NR<<" ("<<m_refParam*100<<"%)"<<((it==m_refPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}

template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if<!_coarsen && !_admissible, void>::type
gsAdaptiveMeshing<T>::_markPercentage_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Percentage marking (not admissible) for refinement...\n";
    // Total number of elements:
    size_t NE = elements.size();
    // Compute the index from which the refinement should start,
    // once the vector is sorted.
    index_t NR = cast<T,index_t>( math::floor( m_refParam * T(NE) ) );

    index_t nmarked = 0;
    auto loop_action = [&elements,&predicates,&NR,&nmarked,&elMarked]
                    (const index_t & index)
    {
        HBox * box = elements.at(index);
        bool check = true;
        for (typename std::vector<gsHBoxCheck<2,T>*>::const_iterator errIt = predicates.begin(); errIt!=predicates.end(); errIt++)
            check &= (*errIt)->check(*box);

        if (check)
        {
            elMarked.add(*box);
            nmarked++;
        }

        return (nmarked > NR);
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_refPermutation.cbegin(),m_refPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark percentage] Marked "<<elMarked.totalSize()<<", ("<<nmarked<<") elements ("<<(T)nmarked/NE*100<<"%"<<" of NE "<<NE<<") and threshold = "<<NR<<" ("<<m_refParam*100<<"%)"<<((it==m_refPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}

template <class T>
T gsAdaptiveMeshing<T>::_maxError(const boxMapType & elements)
{
    auto larger_than = [](const typename boxMapType::value_type & a, const typename boxMapType::value_type & b)
    {
        return (a.second->error() < b.second->error());
    };

    // First, conduct a brutal search for the maximum local error
    T maxErr = (std::max_element(elements.begin(), elements.end(), larger_than ))->second->error();

    return maxErr;
}

template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if< _coarsen &&  _admissible, void>::type
gsAdaptiveMeshing<T>::_markThreshold_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Threshold marking for coarsening...\n";
    GISMO_ASSERT(m_crsParam<=1 && m_crsParam>=0,"Coarsening parameter must be a percentage!");

    T Thr = m_crsParam * m_maxError;
    T current = 0;
    gsHBoxCheck<2,T> * thres_predicate = new gsLargerErrCompare<2,T>(Thr);

    auto loop_action = [this,&elements,&thres_predicate,&predicates,&elMarked,&current]
                    (const index_t & index)
    {
        HBox_ptr box = elements.at(index);

        if (thres_predicate->check(*box))
        {
            current = std::max(current,box->error());
            return true;
        }

        // Continue if the box is already marked or if it does not satsfy all checks
        if (_boxPtr(*box)->marked() || !_checkBox(*box,predicates))
            return false;


        // Get the neighborhoods
        typename HBox::Container sibs = box->getSiblings();
        _setContainerProperties(sibs);
        if (!gsHBoxUtils<2,T>::allActive(sibs))
            return false;

        // Check all siblings if they satisfy the predicates
        if (_checkBoxes(sibs,predicates))
        {
            _addAndMark(sibs,elMarked);
            _addAndMark(*box,elMarked);
        }
        return false;
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_crsPermutation.cbegin(),m_crsPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark threshold] Marked "<<elMarked.totalSize()<<" elements with largest error "<<current<<" and treshold = "<<Thr<<((it==m_crsPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}

template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if< _coarsen && !_admissible, void>::type
gsAdaptiveMeshing<T>::_markThreshold_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Threshold marking for coarsening...\n";
    GISMO_ASSERT(m_crsParam<=1 && m_crsParam>=0,"Coarsening parameter must be a percentage!");

    T Thr = m_crsParam * m_maxError;
    T current = 0;
    gsHBoxCheck<2,T> * thres_predicate = new gsLargerErrCompare<2,T>(Thr);

    auto loop_action = [this,&elements,&thres_predicate,&predicates,&elMarked,&current]
                    (const index_t & index)
    {
        HBox_ptr box = elements.at(index);

        if (thres_predicate->check(*box))
        {
            current = std::max(current,box->error());
            return true;
        }

        // Check all children if they satisfy the predicates
        if (_checkBox(*box,predicates))
            _addAndMark(*box,elMarked);
        return false;
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_crsPermutation.cbegin(),m_crsPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark threshold] Marked "<<elMarked.totalSize()<<" elements with largest error "<<current<<" and treshold = "<<Thr<<((it==m_crsPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}

template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if<!_coarsen &&  _admissible, void>::type
gsAdaptiveMeshing<T>::_markThreshold_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Threshold marking (admissible) for refinement...\n";
    GISMO_ASSERT(m_refParam<=1 && m_refParam>=0,"Refinement parameter must be a percentage!");

    T Thr = m_refParam * m_maxError;
    T current = 0;
    gsHBoxCheck<2,T> * thres_predicate = new gsSmallerErrCompare<2,T>(Thr);

    auto loop_action = [this,&elements,&thres_predicate,&predicates,&elMarked,&current]
                    (const index_t & index)
    {
        HBox_ptr box = elements.at(index);

        if (thres_predicate->check(*box))
        {
            current = std::max(current,box->error());
            return true;
        }

        if (_boxPtr(*box)->marked())
        {
            return false;
        }

        typename HBox::Container neighborhood = HBoxUtils::toContainer(HBoxUtils::markAdmissible(*box,m_m));
        _setContainerProperties(neighborhood);

        // Check all elements in the neighborhood if they satisfy the predicates
        if (_checkBoxes(neighborhood,predicates))
            _addAndMark(neighborhood,elMarked);
        return false;
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_refPermutation.cbegin(),m_refPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark threshold] Marked "<<elMarked.totalSize()<<" elements with largest error "<<current<<" and treshold = "<<Thr<<((it==m_refPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}

template <class T>
template<bool _coarsen,bool _admissible>
typename std::enable_if<!_coarsen && !_admissible, void>::type
gsAdaptiveMeshing<T>::_markThreshold_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<2,T> *> predicates, typename gsAdaptiveMeshing<T>::HBoxContainer & elMarked) const
{
    gsDebug<<"Threshold marking (not admissible) for refinement...\n";
    GISMO_ASSERT(m_refParam<=1 && m_refParam>=0,"Refinement parameter must be a percentage!");

    T Thr = m_refParam * m_maxError;
    T current = 0;
    gsHBoxCheck<2,T> * thres_predicate = new gsSmallerErrCompare<2,T>(Thr);

    auto loop_action = [this,&elements,&thres_predicate,&predicates,&elMarked,&current]
                    (const index_t & index)
    {
        HBox_ptr box = elements.at(index);

        if (thres_predicate->check(*box))
        {
            current = std::max(current,box->error());
            return true;
        }

        if (_checkBox(*box,predicates))
            _addAndMark(*box,elMarked);

        return false;
    };

    std::vector<index_t>::const_iterator it = std::find_if(m_refPermutation.cbegin(),m_refPermutation.cend(),loop_action);
    elMarked = HBoxUtils::Unique(elMarked);
    gsDebug<<"[Mark threshold] Marked "<<elMarked.totalSize()<<" elements with largest error "<<current<<" and treshold = "<<Thr<<((it==m_refPermutation.end()) ? " (maximum number marked)" : "")<<"\n";
}

template<class T>
void gsAdaptiveMeshing<T>::defaultOptions()
{
    m_options.addInt("CoarsenRule","Rule used for coarsening: 1=GARU, 2=PUCA, 3=BULK.",1);
    m_options.addInt("RefineRule","Rule used for refinement: 1=GARU, 2=PUCA, 3=BULK.",1);
    m_options.addReal("CoarsenParam","Parameter used for coarsening",0.1);
    m_options.addReal("RefineParam","Parameter used for refinement",0.1);
    m_options.addReal("CoarsenParamExtra","Extra parameter used for coarsening",0.1);
    m_options.addReal("RefineParamExtra","Extra parameter used for refinement",0.1);

    m_options.addInt("Convergence_alpha","Estimated convergence parameter of he error, for alpha*p+beta convergence",-1);
    m_options.addInt("Convergence_beta","Estimated convergence parameter of he error, for alpha*p+beta convergence",-1);

    m_options.addInt("CoarsenExtension","Extension coarsening",0);
    m_options.addInt("RefineExtension","Extension refinement",0);

    m_options.addInt("MaxLevel","Maximum refinement level",3);

    m_options.addInt("Admissibility","Admissibility region, 0=T-admissibility (default), 1=H-admissibility",0);
    m_options.addSwitch("Admissible","Mark the admissible region",true);
    m_options.addInt("Jump","Jump parameter m",2);
}

template<class T>
void gsAdaptiveMeshing<T>::getOptions()
{
    switch (m_options.askInt("CoarsenRule",3))
    {
        case 1:
            m_crsRule = GARU;
            break;
        case 2:
            m_crsRule = PUCA;
            break;
        case 3:
            m_crsRule = BULK;
            break;
        case 4:
            m_crsRule = PBULK;
            break;
        default:
            GISMO_ERROR("Coarsening marking strategy unknown");
            break;
    }

    switch (m_options.askInt("RefineRule",3))
    {
        case 1:
            m_refRule = GARU;
            break;
        case 2:
            m_refRule = PUCA;
            break;
        case 3:
            m_refRule = BULK;
            break;
        case 4:
            m_refRule = PBULK;
            break;
        default:
            GISMO_ERROR("Refinement marking strategy unknown");
            break;
    }

    m_crsParam = m_options.getReal("CoarsenParam");
    m_crsParamExtra = m_options.getReal("CoarsenParamExtra");
    m_refParam = m_options.getReal("RefineParam");
    m_refParamExtra = m_options.getReal("RefineParamExtra");

    m_crsExt = m_options.getInt("CoarsenExtension");
    m_refExt = m_options.getInt("RefineExtension");

    m_maxLvl = m_options.getInt("MaxLevel");

    m_admissible = m_options.getSwitch("Admissible");
    m_m = m_options.getInt("Jump");

    m_alpha=m_options.askInt("Convergence_alpha",-1);
    m_beta=m_options.askInt("Convergence_beta",-1);
}

template<class T>
void gsAdaptiveMeshing<T>::container_into(const std::vector<T> & elError, HBoxContainer & result)
{
    result.clear();
    this->_assignErrors(m_boxes,elError);
    for (typename boxMapType::iterator it = m_boxes.begin(); it!=m_boxes.end(); it++)
        result.add(*it->second);
}


template<class T>
void gsAdaptiveMeshing<T>::markRef_into(const std::vector<T> & elError, HBoxContainer & elMarked)
{
    elMarked.clear();
    this->_assignErrors(m_boxes,elError);
    if (m_refRule!=PBULK)
        m_refPermutation = this->_sortPermutation(m_boxes); // Index of the lowest error is first
    else
        m_refPermutation = this->_sortPermutationProjectedRef(m_boxes); // Index of the lowest error is first

    // To do:
    // - sort the errors for coarsened and refined in separate functions
    // try again
    // - computeError(primalL,dualL,dualH) class in gsThinSHellAssemblerDWR which derives from gsFunctionSet


    std::reverse(m_refPermutation.begin(),m_refPermutation.end()); // Index of the highest error is first

    std::vector<gsHBoxCheck<2,T> *> predicates;
    _refPredicates_into(predicates);

    if (m_admissible)
        _markElements<false,true>( elError, m_refRule, predicates, elMarked);//,flag [coarse]);
    else
        _markElements<false,false>( elError, m_refRule, predicates, elMarked);//,flag [coarse]);

    for (typename std::vector<gsHBoxCheck<2,T>*>::iterator pred=predicates.begin(); pred!=predicates.end(); pred++)
        delete *pred;
}

template<class T>
void gsAdaptiveMeshing<T>::markCrs_into(const std::vector<T> & elError, const HBoxContainer & markedRef, HBoxContainer & elMarked)
{
    elMarked.clear();
    this->_assignErrors(m_boxes,elError);

    if (m_crsRule!=PBULK)
        m_crsPermutation = this->_sortPermutation(m_boxes); // Index of the lowest error is first
    else
        m_crsPermutation = this->_sortPermutationProjectedCrs(m_boxes); // Index of the lowest error is first
    gsDebugVar(m_crsPermutation.size());

    std::vector<gsHBoxCheck<2,T> *> predicates;
    if (markedRef.totalSize()==0 || !m_admissible)
        _crsPredicates_into(predicates);
    else
        _crsPredicates_into(markedRef,predicates);

    if (m_admissible)
        _markElements<true,true>( elError, m_crsRule, predicates, elMarked);//,flag [coarse]);
    else
        _markElements<true,false>( elError, m_crsRule, predicates, elMarked);//,flag [coarse]);

    gsDebugVar(m_crsPermutation.size());

    for (typename std::vector<gsHBoxCheck<2,T>*>::iterator pred=predicates.begin(); pred!=predicates.end(); pred++)
        delete *pred;
}

template<class T>
void gsAdaptiveMeshing<T>::markCrs_into(const std::vector<T> & elError, HBoxContainer & elMarked)
{
    HBoxContainer container;
    this->markCrs_into(elError,container,elMarked);
}

template<class T>
void gsAdaptiveMeshing<T>::markRef(const std::vector<T> & errors)
{
    markRef_into( errors, m_markedRef);//,flag [coarse]);
}

template<class T>
void gsAdaptiveMeshing<T>::markCrs(const std::vector<T> & errors)
{
    markCrs_into( errors, m_markedRef, m_markedCrs);//,flag [coarse]);
}

template<class T>
bool gsAdaptiveMeshing<T>::refine(const HBoxContainer & markedRef)
{
    bool refine;
    if ((refine = markedRef.totalSize()>0))
        _refineMarkedElements(markedRef,m_refExt);

    return refine;
}

template<class T>
bool gsAdaptiveMeshing<T>::unrefine(const HBoxContainer & markedCrs)
{
    bool coarsen;
    if ((coarsen = markedCrs.totalSize()>0))
        _unrefineMarkedElements(markedCrs,m_crsExt);
    return coarsen;
}

template<class T>
bool gsAdaptiveMeshing<T>::refineAll()
{
    HBoxContainer ref;
    for (typename boxMapType::iterator it = m_boxes.begin(); it!=m_boxes.end(); it++)
        ref.add(*it->second);

    this->refine(ref);

    return true;
}

template<class T>
bool gsAdaptiveMeshing<T>::unrefineAll()
{
    HBoxContainer crs;
    for (typename boxMapType::iterator it = m_boxes.begin(); it!=m_boxes.end(); it++)
        crs.add(*it->second);

    this->unrefine(crs);

    return true;
}

// template<class T>
// void gsAdaptiveMeshing<T>::flatten(const index_t level)
// {
//     _flattenElementsToLevel(level);
// }

// template<class T>
// void gsAdaptiveMeshing<T>::unrefineThreshold(const index_t level)
// {
//     _unrefineElementsThreshold(level);
// }

template<class T>
void gsAdaptiveMeshing<T>::_refineMarkedElements(   const HBoxContainer & markedRef,
                                                    index_t refExtension)
{
    gsBasis<T> * basis = nullptr;

    gsMultiPatch<T> * mp;
    gsMultiBasis<T> * mb;

    for (index_t pn=0; pn < m_input->nPieces(); ++pn )// for all patches
    {
        if ( (mp = dynamic_cast<gsMultiPatch<T>*>(m_input)) != nullptr ) basis = &(mp->basis(pn));
        if ( (mb = dynamic_cast<gsMultiBasis<T>*>(m_input)) != nullptr ) basis = &(mb->basis(pn));
        GISMO_ENSURE(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");

        // if (m_options.getSwitch("Admissible"))
        // {
        //     if (m_options.getInt("Admissibility")==0)
        //     {
        //         if      ( gsTHBSplineBasis<2,T>  * g = dynamic_cast<gsTHBSplineBasis<2,T> *>( basis ) )
        //             marked.markTadmissible(m_m);
        //         else if (  gsHBSplineBasis<2,T>  * g = dynamic_cast< gsHBSplineBasis<2,T> *>( basis ) )
        //             marked.markHadmissible(m_m);
        //         else // if basis type unknown
        //             marked.markHadmissible(m_m);
        //     }
        //     else if (m_options.getInt("Admissibility")==1)
        //             marked.markTadmissible(m_m);
        //     else if (m_options.getInt("Admissibility")==2)
        //             marked.markHadmissible(m_m);
        //     else // if basis type unknown
        //         GISMO_ERROR("Admissibility type unknown or basis type not recognized");

        //     HBoxContainer container(marked.toUnitBoxes());
        //     gsDebugVar(container);
        //     std::vector<index_t> boxes = container.toRefBoxes();
        //     if ((mp = dynamic_cast<gsMultiPatch<T>*>(m_input)))
        //         mp->patch(pn).refineElements(container.toRefBoxes());
        //     else if ((mb = dynamic_cast<gsMultiBasis<T>*>(m_input)))
        //         mb->basis(pn).refineElements(container.toRefBoxes());
        //     else
        //         GISMO_ERROR("No gsMultiPatch or gsMultiBasis found; don't know what to refine");
        // }
        // else
        // {
            gsHBoxContainer<2,T> container = markedRef.patch(pn);
            container.toUnitBoxes();
            if (refExtension==0)
                if ((mp = dynamic_cast<gsMultiPatch<T>*>(m_input)))
                    mp->patch(pn).refineElements( container.toRefBoxes(pn) );
                else if ((mb = dynamic_cast<gsMultiBasis<T>*>(m_input)))
                    mb->basis( pn).refineElements( container.toRefBoxes(pn) );
                else
                    GISMO_ERROR("No gsMultiPatch or gsMultiBasis found");
            else
                if ((mp = dynamic_cast<gsMultiPatch<T>*>(m_input)))
                {
                    // Refine all of the found refBoxes in this patch
                    std::vector<index_t> elements = mp->patch(pn).basis().asElements(container.toCoords(pn), refExtension);
                    mp->patch(pn).refineElements( elements );
                }
                else if ((mb = dynamic_cast<gsMultiBasis<T>*>(m_input)))
                {
                    // Refine all of the found refBoxes in this patch
                    mb->basis( pn).refine(container.toCoords(pn), refExtension );
                }
                else
                    GISMO_ERROR("No gsMultiPatch or gsMultiBasis found");
        // }

    }
}

template<class T>
void gsAdaptiveMeshing<T>::_unrefineMarkedElements(     const HBoxContainer & markedCrs,
                                                        index_t crsExtension)
{
    gsBasis<T> * basis = nullptr;

    gsMultiPatch<T> * mp;
    gsMultiBasis<T> * mb;

    for (index_t pn=0; pn < m_input->nPieces(); ++pn )// for all patches
    {
        if ( (mp = dynamic_cast<gsMultiPatch<T>*>(m_input)) ) basis = &(mp->basis(pn));
        if ( (mb = dynamic_cast<gsMultiBasis<T>*>(m_input)) ) basis = &(mb->basis(pn));
        GISMO_ENSURE(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");

        // if (m_options.getSwitch("Admissible"))
        // {
        //     if (m_options.getInt("Admissibility")==0)
        //     {
        //         if      ( gsTHBSplineBasis<2,T>  * g = dynamic_cast<gsTHBSplineBasis<2,T> *>( basis ) )
        //             marked.markTadmissible(m_m);
        //         else if (  gsHBSplineBasis<2,T>  * g = dynamic_cast< gsHBSplineBasis<2,T> *>( basis ) )
        //             marked.markHadmissible(m_m);
        //         else // if basis type unknown
        //             marked.markHadmissible(m_m);
        //     }
        //     else if (m_options.getInt("Admissibility")==1)
        //             marked.markTadmissible(m_m);
        //     else if (m_options.getInt("Admissibility")==2)
        //             marked.markHadmissible(m_m);
        //     else // if basis type unknown
        //         GISMO_ERROR("Admissibility type unknown or basis type not recognized");

        //     HBoxContainer container(marked.toUnitBoxes());
        //     if ((mp = dynamic_cast<gsMultiPatch<T>*>(m_input)))
        //         mp->patch(pn).unrefineElements(container.toCrsBoxes());
        //     else if ((mb = dynamic_cast<gsMultiBasis<T>*>(m_input)))
        //         mb->basis(pn).unrefineElements(container.toCrsBoxes());
        //     else
        //         GISMO_ERROR("No gsMultiPatch or gsMultiBasis found; don't know what to refine");
        // }
        // else
        // {
        gsHBoxContainer<2,T> container = markedCrs.patch(pn);
        container.toUnitBoxes();
        if (crsExtension==0)
            if ((mp = dynamic_cast<gsMultiPatch<T>*>(m_input)))
                mp->patch(pn).unrefineElements( container.toCrsBoxes(pn) );
            else if ((mb = dynamic_cast<gsMultiBasis<T>*>(m_input)))
                mb->basis( pn).unrefineElements( container.toCrsBoxes(pn) );
            else
                GISMO_ERROR("No gsMultiPatch or gsMultiBasis found");
        else
            if ((mp = dynamic_cast<gsMultiPatch<T>*>(m_input))!= nullptr)
            {
                // Unrefine all of the found refBoxes in this patch
                std::vector<index_t> elements = mp->patch(pn).basis().asElementsUnrefine(markedCrs.toCoords(pn), crsExtension);
                mp->patch(pn).unrefineElements( elements );
            }
            else if ((mb = dynamic_cast<gsMultiBasis<T>*>(m_input))!= nullptr)
            {
                // Unrefine all of the found refBoxes in this patch
                mb->unrefine( pn, markedCrs.toCoords(pn), crsExtension );
            }
            else
                GISMO_ERROR("No gsMultiPatch or gsMultiBasis found");
        // }

    }
}

template<class T>
typename gsAdaptiveMeshing<T>::HBoxContainer gsAdaptiveMeshing<T>::_toContainer( const std::vector<bool> & bools) const
{
    HBoxContainer container;

    #pragma omp parallel
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
        const int nt  = omp_get_num_threads();
        index_t patch_cnt = 0;
#endif

        index_t c = 0;
        gsBasis<T> * basis = nullptr;

        gsMultiPatch<T> * mp;
        gsMultiBasis<T> * mb;
        typename gsBasis<T>::domainIter domIt;
        gsHDomainIterator<T,2> * domHIt = nullptr;
        for (index_t patchInd=0; patchInd < m_input->nPieces(); ++patchInd)
        {
            // Initialize domain element iterator
            if ( (mp = dynamic_cast<gsMultiPatch<T>*>(m_input))!= nullptr ) basis = &(mp->basis(patchInd));
            if ( (mb = dynamic_cast<gsMultiBasis<T>*>(m_input))!= nullptr ) basis = &(mb->basis(patchInd));
            GISMO_ASSERT(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");
            // for all elements in patch pn
            domIt  = basis->makeDomainIterator();
            domHIt = dynamic_cast<gsHDomainIterator<T,2> *>(domIt.get());
            GISMO_ENSURE(domHIt!=nullptr,"Domain should be hierarchical");

#ifdef _OPENMP
            c = patch_cnt + tid;
            patch_cnt += domHIt->numElements();// a bit costly
            for ( domHIt->next(tid); domHIt->good(); domHIt->next(nt) )
#else
            for (; domHIt->good(); domHIt->next() )
#endif
            {
                if (bools[c])
#pragma omp critical (gsAdaptiveMeshingtoContainer)
                    container.add(HBox(domHIt,patchInd));

#               ifdef _OPENMP
                c += nt;
#               else
                c++;
#               endif
            }
        }
    }

    return container;
}

template<class T>
index_t gsAdaptiveMeshing<T>::numBlocked() const
{
    gsMaxLvlCompare<2,T> comp(m_maxLvl);
    index_t numBlocked = 0;
    for (typename boxMapType::const_iterator it=m_boxes.cbegin(); it!=m_boxes.cend(); it++)
        numBlocked += comp.check(*it->second);

    return numBlocked;
}

template<class T>
index_t gsAdaptiveMeshing<T>::numElements() const
{
    return m_indices.size();
}

template<class T>
void gsAdaptiveMeshing<T>::assignErrors(const std::vector<T> & elError)
{
    this->_assignErrors(m_boxes,elError);
}

template<class T>
T gsAdaptiveMeshing<T>::blockedError() const
{
    gsMaxLvlCompare<2,T> comp(m_maxLvl);
    T error = 0;
    for (typename boxMapType::const_iterator it=m_boxes.cbegin(); it!=m_boxes.cend(); it++)
    {
        if (!(comp.check(*it->second)))
            error += it->second->error();
    }

    return error;
}

template<class T>
T gsAdaptiveMeshing<T>::nonBlockedError() const
{
    gsMaxLvlCompare<2,T> comp(m_maxLvl);
    T error = 0;
    for (typename boxMapType::const_iterator it=m_boxes.cbegin(); it!=m_boxes.cend(); it++)
    {
        if (comp.check(*it->second))
            error += it->second->error();
    }

    return error;
}

// The functions below are deprecated

/**
 * @brief      Flattens everything to \a level
 *
 * @param      m_input  The m_input structure (gsMultiBasis or gsMultiPatch)
 * @param[in]  level  The target level
 *
 */
// template<class T>
// void gsAdaptiveMeshing<T>::_flattenElementsToLevel(const index_t level)
// {
//     // Get all the elements of which the level exceeds level
//     std::vector<bool> elMarked;
//     _markLevelThreshold(level,elMarked);

//     const index_t dim = m_input->domainDim();

//     // numMarked: Number of marked cells on current patch, also currently marked cell
//     // poffset  : offset index for the first element on a patch
//     // globalCount: counter for the current global element index
//     index_t numMarked, poffset = 0, globalCount = 0;

//     // crsBoxes: contains marked boxes on a given patch
//     gsMatrix<T> crsBoxes;

//     gsBasis<T> * basis = nullptr;

//     gsMultiPatch<T> * mp;
//     gsMultiBasis<T> * mb;
//     for (index_t pn=0; pn < m_input->nPieces(); ++pn )// for all patches
//     {
//         if ( (mp = dynamic_cast<gsMultiPatch<T>*>(m_input)) ) basis = &(mp->basis(pn));
//         if ( (mb = dynamic_cast<gsMultiBasis<T>*>(m_input)) ) basis = &(mb->basis(pn));
//         GISMO_ASSERT(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");
//         // Get number of elements to be refined on this patch
//         const index_t numEl = basis->numElements();
//         numMarked = std::count_if(elMarked.begin() + poffset,
//                                   elMarked.begin() + poffset + numEl,
//                                   GS_BIND2ND(std::equal_to<bool>(), true) );

//         poffset += numEl;
//         crsBoxes.resize(dim, 2*numMarked);
//         numMarked = 0;// counting current patch element to be refined

//         // for all elements in patch pn
//         typename gsBasis<T>::domainIter domIt = basis->makeDomainIterator();
//         for (; domIt->good(); domIt->next())
//         {
//             if( elMarked[ globalCount++ ] ) // refine this element ?
//             {
//                 // Construct degenerate box by setting both
//                 // corners equal to the center
//                 crsBoxes.col(2*numMarked  ) =
//                         crsBoxes.col(2*numMarked+1) = domIt->centerPoint();

//                 // Advance marked cells counter
//                 numMarked++;
//             }
//         }

//         std::vector<index_t> elements;
//         elements = basis->asElementsUnrefine(crsBoxes,0);
//         const index_t offset = 2*dim+1;
//         index_t diff = 0;
//         GISMO_ASSERT(elements.size()%offset==0,"Element boxes have wrong size. Boxes should have size "<<offset<<" per box");
//         for (size_t k = 0; k!=elements.size()/offset; k++)
//         {
//             // gsDebug<<"Old Box:\n";
//             // for (index_t kk = 0; kk!=2*dim+1; kk++)
//             // {
//             //     gsDebug<<elements[offset*k+kk]<<",";
//             // }
//             // gsDebug<<"\n";

//             if ((diff = elements[offset*k] - level) > 0)
//             {
//                 elements[offset*k] = level;
//                 for (index_t kk = 0; kk!=2*dim; kk++)
//                 {

//                     elements[offset*k+kk+1] = elements[offset*k+kk+1] >> diff;
//                     // Remove comment below
//                     // // yes, exactly: if   i is index in old level a and needs to become a level b box, b<a, then the new index j is (in c++ computation):
//                     // j = i >> a-b;
//                     // // this is division by 2^{a-b}  using bit-operations
//                 }
//             }
//             // gsDebug<<"New Box:\n";
//             // for (index_t kk = 0; kk!=2*dim+1; kk++)
//             // {
//             //     gsDebug<<elements[offset*k+kk]<<",";
//             // }
//             // gsDebug<<"\n";
//         }

//         if ((mp = dynamic_cast<gsMultiPatch<T>*>(m_input)))
//         {
//             mp->patch(pn).unrefineElements( elements );
//         }
//         else if ((mb = dynamic_cast<gsMultiBasis<T>*>(m_input)))
//         {
//             // Refine all of the found refBoxes in this patch
//             mb->unrefineElements(pn, elements );
//         }
//         else
//             GISMO_ERROR("No gsMultiPatch or gsMultiBasis found");
//     }
// }

/**
 * @brief      Coarsens everything to that is lower than \a level
 *
 * note: DOES NOT TAKE INTO ACCOUNT SUITABLE GRADING
 *
 * @param      m_input  The m_input structure (gsMultiBasis or gsMultiPatch)
 * @param[in]  level  The target level
 *
 */
// template<class T>
// void gsAdaptiveMeshing<T>::_unrefineElementsThreshold(
//                                                 const index_t level)
// {
//     // Get all the elements of which the level exceeds level
//     std::vector<bool> elMarked;
//     _markLevelThreshold(level,elMarked);

//     const index_t dim = m_input->domainDim();

//     // numMarked: Number of marked cells on current patch, also currently marked cell
//     // poffset  : offset index for the first element on a patch
//     // globalCount: counter for the current global element index
//     index_t numMarked, poffset = 0, globalCount = 0;

//     // crsBoxes: contains marked boxes on a given patch
//     gsMatrix<T> crsBoxes;

//     gsBasis<T> * basis = nullptr;

//     gsMultiPatch<T> * mp;
//     gsMultiBasis<T> * mb;
//     for (index_t pn=0; pn < m_input->nPieces(); ++pn )// for all patches
//     {
//         if ( (mp = dynamic_cast<gsMultiPatch<T>*>(m_input)) ) basis = &(mp->basis(pn));
//         if ( (mb = dynamic_cast<gsMultiBasis<T>*>(m_input)) ) basis = &(mb->basis(pn));
//         GISMO_ASSERT(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");
//         // Get number of elements to be refined on this patch
//         const index_t numEl = basis->numElements();
//         numMarked = std::count_if(elMarked.begin() + poffset,
//                                   elMarked.begin() + poffset + numEl,
//                                   GS_BIND2ND(std::equal_to<bool>(), true) );

//         poffset += numEl;
//         crsBoxes.resize(dim, 2*numMarked);
//         numMarked = 0;// counting current patch element to be refined

//         // for all elements in patch pn
//         typename gsBasis<T>::domainIter domIt = basis->makeDomainIterator();
//         for (; domIt->good(); domIt->next())
//         {
//             if( elMarked[ globalCount++ ] ) // refine this element ?
//             {
//                 // Construct degenerate box by setting both
//                 // corners equal to the center
//                 crsBoxes.col(2*numMarked  ) =
//                         crsBoxes.col(2*numMarked+1) = domIt->centerPoint();

//                 // Advance marked cells counter
//                 numMarked++;
//             }
//         }

//         std::vector<index_t> elements;
//         elements = basis->asElementsUnrefine(crsBoxes,0);
//         const index_t offset = 2*dim+1;
//         index_t diff = 0;
//         GISMO_ASSERT(elements.size()%offset==0,"Element boxes have wrong size. Boxes should have size "<<offset<<" per box");
//         for (size_t k = 0; k!=elements.size()/offset; k++)
//         {
//             // gsDebug<<"Old Box:\n";
//             // for (index_t kk = 0; kk!=2*dim+1; kk++)
//             // {
//             //     gsDebug<<elements[offset*k+kk]<<",";
//             // }
//             // gsDebug<<"\n";

//             if ((diff = elements[offset*k] - level) > 0)
//             {
//                 elements[offset*k] = level;
//                 for (index_t kk = 0; kk!=2*dim; kk++)
//                 {

//                     elements[offset*k+kk+1] = elements[offset*k+kk+1] >> 1;

//                     // // yes, exactly: if   i is index in old level a and needs to become a level b box, b<a, then the new index j is (in c++ computation):
//                     // j = i >> a-b;
//                     // // this is division by 2^{a-b}  using bit-operations
//                 }
//             }
//             // gsDebug<<"New Box:\n";
//             // for (index_t kk = 0; kk!=2*dim+1; kk++)
//             // {
//             //     gsDebug<<elements[offset*k+kk]<<",";
//             // }
//             // gsDebug<<"\n";
//         }

//         if ((mp = dynamic_cast<gsMultiPatch<T>*>(m_input)))
//         {
//             mp->patch(pn).unrefineElements( elements );
//         }
//         else if ((mb = dynamic_cast<gsMultiBasis<T>*>(m_input)))
//         {
//             // Refine all of the found refBoxes in this patch
//             mb->unrefineElements(pn, elements );
//         }
//         else
//             GISMO_ERROR("No gsMultiPatch or gsMultiBasis found");
//     }
// }

// template <class T>
// void gsAdaptiveMeshing<T>::_getElLevels(  std::vector<index_t> & elLevels)
// {
//     // Now just check for each element, whether the level
//     // is above the target level or not, and mark accordingly.
//     //
//     index_t Nelements = 0;
//     if ( gsMultiPatch<T> * mp = dynamic_cast<gsMultiPatch<T>*>(m_input) ) Nelements = gsMultiBasis<T>(*mp).totalElements();
//     if ( gsMultiBasis<T> * mb = dynamic_cast<gsMultiBasis<T>*>(m_input) ) Nelements = mb->totalElements();
//     GISMO_ASSERT(Nelements!=0,"Number of elements is zero? It might be that the m_input is not a gsMultiBasis or gsMultiPatch");

//     elLevels.resize(Nelements);

//     gsBasis<T> * basis = nullptr;

//     #pragma omp parallel
//     {
// #ifdef _OPENMP
//         const int tid = omp_get_thread_num();
//         const int nt  = omp_get_num_threads();
//         index_t patch_cnt = 0;
// #endif

//         index_t c = 0;
//         for (index_t patchInd=0; patchInd < m_input->nPieces(); ++patchInd)
//         {
//             // Initialize domain element iterator
//             if ( gsMultiPatch<T> * mp = dynamic_cast<gsMultiPatch<T>*>(m_input) ) basis = &(mp->basis(patchInd));
//             if ( gsMultiBasis<T> * mb = dynamic_cast<gsMultiBasis<T>*>(m_input) ) basis = &(mb->basis(patchInd));
//             GISMO_ASSERT(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");
//             // for all elements in patch pn
//             typename gsBasis<T>::domainIter domIt = basis->makeDomainIterator();
//             gsHDomainIterator<T,2> * domHIt = nullptr;
//             domHIt = dynamic_cast<gsHDomainIterator<T,2> *>(domIt.get());
//             GISMO_ENSURE(domHIt!=nullptr,"Domain should be 2 dimensional for flattening");

// #ifdef _OPENMP
//             c = patch_cnt + tid;
//             patch_cnt += domHIt->numElements();// a bit costy
//             for ( domHIt->next(tid); domHIt->good(); domHIt->next(nt) )
// #else
//             for (; domHIt->good(); domHIt->next() )
// #endif
//             {
// #               ifdef _OPENMP
//                 elLevels[c] = domHIt->getLevel();
//                 c += nt;
// #               else
//                 elLevels[c++] = domHIt->getLevel();
// #               endif
//             }
//         }
//     }
// }

// template <class T>
// void gsAdaptiveMeshing<T>::_markLevelThreshold(  index_t level, HBoxContainer & elMarked)
// {
//     // Now just check for each element, whether the level
//     // is above the target level or not, and mark accordingly.
//     //
//     index_t Nelements = 0;
//     if ( gsMultiPatch<T> * mp = dynamic_cast<gsMultiPatch<T>*>(m_input) ) Nelements = gsMultiBasis<T>(*mp).totalElements();
//     if ( gsMultiBasis<T> * mb = dynamic_cast<gsMultiBasis<T>*>(m_input) ) Nelements = mb->totalElements();
//     GISMO_ASSERT(Nelements!=0,"Number of elements is zero? It might be that the m_input is not a gsMultiBasis or gsMultiPatch");

//     elMarked.resize(Nelements);

//     gsBasis<T> * basis = nullptr;

//     #pragma omp parallel
//     {
// #ifdef _OPENMP
//         const int tid = omp_get_thread_num();
//         const int nt  = omp_get_num_threads();
//         index_t patch_cnt = 0;
// #endif

//         index_t c = 0;
//         for (index_t patchInd=0; patchInd < m_input->nPieces(); ++patchInd)
//         {
//             // Initialize domain element iterator
//             if ( gsMultiPatch<T> * mp = dynamic_cast<gsMultiPatch<T>*>(m_input) ) basis = &(mp->basis(patchInd));
//             if ( gsMultiBasis<T> * mb = dynamic_cast<gsMultiBasis<T>*>(m_input) ) basis = &(mb->basis(patchInd));
//             GISMO_ASSERT(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");
//             // for all elements in patch pn
//             typename gsBasis<T>::domainIter domIt = basis->makeDomainIterator();
//             gsHDomainIterator<T,2> * domHIt = nullptr;
//             domHIt = dynamic_cast<gsHDomainIterator<T,2> *>(domIt.get());
//             GISMO_ENSURE(domHIt!=nullptr,"Domain should be 2 dimensional for flattening");

// #ifdef _OPENMP
//             c = patch_cnt + tid;
//             patch_cnt += domHIt->numElements();// a bit costy
//             for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
// #else
//             for (; domIt->good(); domIt->next() )
// #endif
//             {
// #               ifdef _OPENMP
//                 elMarked[c] = (domHIt->getLevel() > level);
//                 c += nt;
// #               else
//                 elMarked[c++] = (domHIt->getLevel() > level);
// #               endif
//             }
//         }
//     }
// }


} // namespace gismo
