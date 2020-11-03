/** @file gsAdaptiveRefUtils.h

    @brief Provides generic routines for adaptive refinement.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

#pragma once


#include <iostream>

#include <gsIO/gsIOUtils.h>

namespace gismo
{

enum MarkingStrategy
{
    GARU=1,
    PUCA=2,
    BULK=3
};

template <class T>
void gsMarkThreshold( const std::vector<T> & elError, T refParameter, std::vector<bool> & elMarked)
{
    // First, conduct a brutal search for the maximum local error
    const T maxErr = *std::max_element(elError.begin(), elError.end() );

    // Compute the threshold:
    const T Thr = refParameter * maxErr;

    elMarked.resize( elError.size() );
    // Now just check for each element, whether the local error
    // is above the computed threshold or not, and mark accordingly.

    typename std::vector<T>::const_iterator err = elError.begin();
    for(std::vector<bool>::iterator i = elMarked.begin(); i!=  elMarked.end(); ++i, ++err)
        *i = ( *err >= Thr );
}

template <class T>
void gsMarkPercentage( const std::vector<T> & elError, T refParameter, std::vector<bool> & elMarked)
{
    T Thr = T(0);

    // Total number of elements:
    size_t NE = elError.size();
    // The vector of local errors will need to be sorted,
    // which will be done on a copy:
    std::vector<T> elErrCopy = elError;

    // Compute the index from which the refinement should start,
    // once the vector is sorted.
    size_t idxRefineStart = cast<T,size_t>( math::floor( refParameter * T(NE) ) );
    // ...and just to be sure we are in range:
    if( idxRefineStart == elErrCopy.size() )
    {
        GISMO_ASSERT(idxRefineStart >= 1, "idxRefineStart can't get negative");
        idxRefineStart -= 1;
    }

    // Sort the list using bubblesort.
    // After each loop, the largest elements are at the end
    // of the list. Since we are only interested in the largest elements,
    // it is enough to run the sorting until enough "largest" elements
    // have been found, i.e., until we have reached indexRefineStart
    size_t lastSwapDone = elErrCopy.size() - 1;
    size_t lastCheckIdx = lastSwapDone;

    bool didSwap;
    T tmp;
    do{
        didSwap = false;
        lastCheckIdx = lastSwapDone;
        for( size_t i=0; i < lastCheckIdx; i++)
            if( elErrCopy[i] > elErrCopy[i+1] )
            {
                tmp = elErrCopy[i];
                elErrCopy[i] = elErrCopy[i+1];
                elErrCopy[i+1] = tmp;

                didSwap = true;
                lastSwapDone = i;
            }
    }while( didSwap && (lastSwapDone+1 >= idxRefineStart ) );

    // Compute the threshold:
    Thr = elErrCopy[ idxRefineStart ];
    elMarked.resize( elError.size() );
    // Now just check for each element, whether the local error
    // is above the computed threshold or not, and mark accordingly.
    for( size_t i=0; i < elError.size(); i++)
        ( elError[i] >= Thr ? elMarked[i] = true : elMarked[i] = false );
}


template <class T>
void gsMarkFraction( const std::vector<T> & elError, T refParameter, std::vector<bool> & elMarked)
{
    T Thr = T(0);

    // The vector of local errors will need to be sorted,
    // which will be done on a copy:
    std::vector<T> elErrCopy = elError;

    // Compute the sum, i.e., the global/total error
    T totalError = T(0);
    for( size_t i = 0; i < elErrCopy.size(); ++i)
        totalError += elErrCopy[i];

    // We want to mark just enough cells such that their
    // cummulated errors add up to a certain fraction
    // of the total error.
    T errorMarkSum = (1-refParameter) * totalError;
    T cummulErrMarked = 0;

    T tmp;
    GISMO_ASSERT(elErrCopy.size() >= 1, "elErrCopy needs at least 1 element");
    size_t lastSwapDone = elErrCopy.size() - 1;
    do{
        for( size_t i=0; i < lastSwapDone; i++)
            if( elErrCopy[i] > elErrCopy[i+1] )
            {
                tmp = elErrCopy[i];
                elErrCopy[i] = elErrCopy[i+1];
                elErrCopy[i+1] = tmp;
            }

        cummulErrMarked += elErrCopy[ lastSwapDone  ];
        lastSwapDone -= 1;

    }while( cummulErrMarked < errorMarkSum && lastSwapDone > 0 );

    // Compute the threshold:
    Thr = elErrCopy[ lastSwapDone + 1 ];
    elMarked.resize( elError.size() );
    // Now just check for each element, whether the local error
    // is above the computed threshold or not, and mark accordingly.
    for( size_t i=0; i < elError.size(); i++)
        ( elError[i] >= Thr ? elMarked[i] = true : elMarked[i] = false );
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
 * Let \f$\rho\f$ denote the input parameter \em refParameter.
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
 * \param refParameter parameter \f$ \rho \f$ for refinement criterion (see above).\n
 * \f$\rho = 0\f$ corresponds to global refinement,\n
 * \f$ \rho=1\f$ corresponds to (almost) no refinement.
 * \param[out] elMarked std::vector of Booleans indicating whether the corresponding element is marked or not.
 *
 * \ingroup Assembler
 */
template <class T>
void gsMarkElementsForRef( const std::vector<T> & elError, int refCriterion, T refParameter, std::vector<bool> & elMarked)
{
    switch (refCriterion)
    {
    case GARU:
        gsMarkThreshold(elError,refParameter,elMarked);
        break;
    case PUCA:
        gsMarkPercentage(elError,refParameter,elMarked);
        break;
    case BULK:
        gsMarkFraction(elError,refParameter,elMarked);
        break;
    default:
        GISMO_ERROR("unknown marking strategy");
    }

}




/** \brief Refine a gsMultiBasis, based on a vector of element-markings.
 *
 * Given the vector of element-markings (see gsRefineMarkedElements()),
 * the corresponding element
 * in the mesh underlying \em basis is refined.
 *
 * It is possible to extend the refinement to the
 * neighbouring elements of the marked area by setting
 * the parameter \a refExtension.\n
 * This parameter is given as number of cells at the level
 * of the marked element \em before refinement.
 *
 * \remarks
 * The ordering/numbering of the elements is implicitly defined by
 * the numbering of the patches in gsMultiBasis, and
 * by the gsDomainIterator of the respective patch-wise basis!
 *
 * \param basis gsMultiBasis to be refined adaptively.
 * \param elMarked std::vector of Booleans indicating
 * for each element of the mesh underlying \em basis, whether it should be refined or not.
 * \param refExtension Specifies how large the refinement extension
 * should be. Given as number of cells at the level \em before refinement.
 *
 * \ingroup Assembler
 *
 * \todo Make gsRefineMarkedElements a member of gsMultiBasis and propagate to gsBasis
 */
template <class T>
void gsRefineMarkedElements(gsMultiBasis<T> & basis,
                            const std::vector<bool> & elMarked,
                            index_t refExtension = 0)
{
    const short_t dim = basis.dim();

    // numMarked: Number of marked cells on current patch, also currently marked cell
    // poffset  : offset index for the first element on a patch
    // globalCount: counter for the current global element index
    index_t numMarked, poffset = 0, globalCount = 0;

    // refBoxes: contains marked boxes on a given patch
    gsMatrix<T> refBoxes;

    for (size_t pn=0; pn < basis.nBases(); ++pn )// for all patches
    {
        // Get number of elements to be refined on this patch
        const size_t numEl = basis[pn].numElements();
        numMarked = std::count_if(elMarked.begin() + poffset,
                                  elMarked.begin() + poffset + numEl,
                                  GS_BIND2ND(std::equal_to<bool>(), true) );
        poffset += numEl;
        refBoxes.resize(dim, 2*numMarked);
        //gsDebugVar(numMarked);
        numMarked = 0;// counting current patch element to be refined

        // for all elements in patch pn
        typename gsBasis<T>::domainIter domIt = basis.basis(pn).makeDomainIterator();
        for (; domIt->good(); domIt->next())
        {
            if( elMarked[ globalCount++ ] ) // refine this element ?
            {
                // Construct degenerate box by setting both
                // corners equal to the center
                refBoxes.col(2*numMarked  ) =
                        refBoxes.col(2*numMarked+1) = domIt->centerPoint();

                // Advance marked cells counter
                numMarked++;
            }
        }
        // Refine all of the found refBoxes in this patch
        basis.refine( pn, refBoxes, refExtension );
    }
}



} // namespace gismo
