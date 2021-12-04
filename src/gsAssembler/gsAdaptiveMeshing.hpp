/** @file gsAdaptiveMeshing.hpp

    @brief Provides class for adaptive refinement.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-)
*/

#pragma once

namespace gismo
{


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
void gsAdaptiveMeshing<T>::_markElements( const std::vector<T> & elError, int refCriterion, T refParameter, std::vector<bool> & elMarked)
{
    switch (refCriterion)
    {
    case GARU:
        _markThreshold(elError,refParameter,elMarked);
        break;
    case PUCA:
        _markPercentage(elError,refParameter,elMarked);
        break;
    case BULK:
        _markFraction(elError,refParameter,elMarked);
        break;
    default:
        GISMO_ERROR("unknown marking strategy");
    }

}

template <class T>
void gsAdaptiveMeshing<T>::_markFraction( const std::vector<T> & elError, T refParameter, std::vector<bool> & elMarked)
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


template <class T>
void gsAdaptiveMeshing<T>::_markPercentage( const std::vector<T> & elError, T refParameter, std::vector<bool> & elMarked)
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
void gsAdaptiveMeshing<T>::_markThreshold( const std::vector<T> & elError, T refParameter, std::vector<bool> & elMarked)
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
void gsAdaptiveMeshing<T>::_markLevelThreshold( gsFunctionSet<T> * input, index_t level, std::vector<bool> & elMarked)
{
    // Now just check for each element, whether the level
    // is above the target level or not, and mark accordingly.
    //
    index_t Nelements = 0;
    if ( gsMultiPatch<T> * mp = dynamic_cast<gsMultiPatch<T>*>(input) ) Nelements = gsMultiBasis<T>(*mp).totalElements();
    if ( gsMultiBasis<T> * mb = dynamic_cast<gsMultiBasis<T>*>(input) ) Nelements = mb->totalElements();
    GISMO_ASSERT(Nelements!=0,"Number of elements is zero? It might be that the input is not a gsMultiBasis or gsMultiPatch");

    elMarked.resize(Nelements);

    gsBasis<T> * basis = nullptr;

    #pragma omp parallel
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
        const int nt  = omp_get_num_threads();
        index_t patch_cnt = 0;
#endif

        index_t c = 0;
        for (index_t patchInd=0; patchInd < input->nPieces(); ++patchInd)
        {
            // Initialize domain element iterator
            if ( gsMultiPatch<T> * mp = dynamic_cast<gsMultiPatch<T>*>(input) ) basis = &(mp->basis(patchInd));
            if ( gsMultiBasis<T> * mb = dynamic_cast<gsMultiBasis<T>*>(input) ) basis = &(mb->basis(patchInd));
            GISMO_ASSERT(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");
            // for all elements in patch pn
            typename gsBasis<T>::domainIter domIt = basis->makeDomainIterator();
            gsHDomainIterator<T,2> * domHIt = nullptr;
            domHIt = dynamic_cast<gsHDomainIterator<T,2> *>(domIt.get());
            GISMO_ENSURE(domHIt!=nullptr,"Domain should be 2 dimensional for flattening");

#ifdef _OPENMP
            c = patch_cnt + tid;
            patch_cnt += domHIt->numElements();// a bit costy
            for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
            for (; domIt->good(); domIt->next() )
#endif
            {
#               ifdef _OPENMP
                elMarked[c] = (domHIt->getLevel() > level);
                c += nt;
#               else
                elMarked[c++] = (domHIt->getLevel() > level);
#               endif
            }
        }
    }
}


template<class T>
void gsAdaptiveMeshing<T>::defaultOptions()
{
    m_options.addInt("CoarsenRule","Rule used for coarsening: 1=GARU, 2=PUCA, 3=BULK.",1);
    m_options.addInt("RefineRule","Rule used for refinement: 1=GARU, 2=PUCA, 3=BULK.",1);
    m_options.addReal("CoarsenParam","Parameter used for coarsening",0.5);
    m_options.addReal("RefineParam","Parameter used for refinement",-1);

    m_options.addInt("CoarsenExtension","Extension coarsening",0);
    m_options.addInt("RefineExtension","Extension refinement",0);

    m_options.addInt("MaxLevel","Maximum refinement level",10);
}

template<class T>
void gsAdaptiveMeshing<T>::getOptions()
{
    if (m_options.getInt("CoarsenRule")==1)
        m_crsRule = GARU;
    else if (m_options.getInt("CoarsenRule")==2)
        m_crsRule = PUCA;
    else if (m_options.getInt("CoarsenRule")==3)
        m_crsRule = BULK;
    else
        GISMO_ERROR("Marking strategy unknown");

    if (m_options.getInt("RefineRule")==1)
        m_refRule = GARU;
    else if (m_options.getInt("RefineRule")==2)
        m_refRule = PUCA;
    else if (m_options.getInt("RefineRule")==3)
        m_refRule = BULK;
    else
        GISMO_ERROR("Marking strategy unknown");

    m_crsParam = m_options.getReal("CoarsenParam");
    m_refParam = m_options.getReal("RefineParam");

    m_crsExt = m_options.getInt("CoarsenExtension");
    m_refExt = m_options.getInt("RefineExtension");

    m_maxLvl = m_options.getInt("MaxLevel");
}

template<class T>
void gsAdaptiveMeshing<T>::mark(const std::vector<T> & errors)
{
    std::vector<T> refErrors(errors.size());
    std::vector<T> crsErrors;

    if (m_crsParam!=-1)
        crsErrors.resize(errors.size());

    if (m_crsParam==-1)
    {
        #pragma omp parallel
        {
        #   ifdef _OPENMP
            const int tid = omp_get_thread_num();
            const int nt  = omp_get_num_threads();
            for (size_t k = tid; k<errors.size(); k+=nt)
        #   else
            for (size_t k = 0; k<errors.size(); k++)
        #   endif
            {
                crsErrors[k] = -math::abs(errors[k]);
                refErrors[k] =  math::abs(errors[k]);
            }
        }
    }
    else
    {
        #pragma omp parallel
        {
        #   ifdef _OPENMP
            const int tid = omp_get_thread_num();
            const int nt  = omp_get_num_threads();
            for (size_t k = tid; k<errors.size(); k+=nt)
        #   else
            for (size_t k = 0; k<errors.size(); k++)
        #   endif
            {
                crsErrors[k] = -math::abs(errors[k]); // DO NOT TAKE ABS HERE BUT IN DWR
                refErrors[k] =  math::abs(errors[k]);
            }
        }
    }

    _markElements( refErrors, m_refRule, m_refParam, m_markedRef);//,flag [coarse]);

    if (m_crsParam!=-1)
        _markElements( crsErrors, m_refRule, m_refParam, m_markedCrs);//,flag [coarse]);

}

template<class T>
void gsAdaptiveMeshing<T>::refine(const std::vector<bool> & markedRef)
{
    GISMO_ASSERT(markedRef.size()!=0,"Mark vector is empty!");

    _refineMarkedElements(m_input,markedRef,m_refExt);
}

template<class T>
void gsAdaptiveMeshing<T>::unrefine(const std::vector<bool> & markedCrs)
{
    GISMO_ASSERT(markedCrs.size()!=0,"Mark vector is empty!");

    _unrefineMarkedElements(m_input,markedCrs,m_crsExt);
}

template<class T>
void gsAdaptiveMeshing<T>::adapt(const std::vector<bool> & markedRef,const std::vector<bool> & markedCrs)
{
    GISMO_ASSERT(markedRef.size()!=0,"Mark vector is empty!");
    GISMO_ASSERT(markedCrs.size()!=0,"Mark vector is empty!");

    _processMarkedElements(m_input,markedRef,markedCrs,m_refExt,m_crsExt);
}

template<class T>
void gsAdaptiveMeshing<T>::flatten(const index_t level)
{
    _flattenElements(m_input,level);
}

template<class T>
void gsAdaptiveMeshing<T>::_refineMarkedElements( gsFunctionSet<T> * input,
                                                    const std::vector<bool> & elMarked,
                                                    index_t refExtension)
{
    // gsMultiPatch<T> * mp;
    // GISMO_ASSERT(mp = dynamic_cast<gsMultiPatch<T>*>(&bases),"No gsMultiBasis!");
    const int dim = input->domainDim();

    // numMarked: Number of marked cells on current patch, also currently marked cell
    // poffset  : offset index for the first element on a patch
    // globalCount: counter for the current global element index
    int numMarked, poffset = 0, globalCount = 0;

    // refBoxes: contains marked boxes on a given patch
    gsMatrix<T> refBoxes;

    gsBasis<T> * basis = nullptr;

    for (index_t pn=0; pn < input->nPieces(); ++pn )// for all patches
    {
        if ( gsMultiPatch<T> * mp = dynamic_cast<gsMultiPatch<T>*>(input) ) basis = &(mp->basis(pn));
        if ( gsMultiBasis<T> * mb = dynamic_cast<gsMultiBasis<T>*>(input) ) basis = &(mb->basis(pn));
        GISMO_ASSERT(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");
        // Get number of elements to be refined on this patch
        const int numEl = basis->numElements();
        numMarked = std::count_if(elMarked.begin() + poffset,
                                  elMarked.begin() + poffset + numEl,
                                  GS_BIND2ND(std::equal_to<bool>(), true) );

        poffset += numEl;
        refBoxes.resize(dim, 2*numMarked);
        numMarked = 0;// counting current patch element to be refined

        // for all elements in patch pn
        typename gsBasis<T>::domainIter domIt = basis->makeDomainIterator();
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

        gsMultiPatch<T> * mp;
        gsMultiBasis<T> * mb;
        if ((mp = dynamic_cast<gsMultiPatch<T>*>(input)))
        {
            std::vector<index_t> elements = mp->patch(pn).basis().asElements(refBoxes, refExtension);
            mp->patch(pn).refineElements( elements );
        }
        else if ((mb = dynamic_cast<gsMultiBasis<T>*>(input)))
        {
            mb->refine( pn, refBoxes, refExtension );
        }
        else
            GISMO_ERROR("No gsMultiPatch or gsMultiBasis found");

        // GISMO_ASSERT(mp = dynamic_cast<gsMultiPatch<T>*>(&bases),"No gsMultiBasis!");
        // // Refine all of the found refBoxes in this patch
        // std::vector<index_t> elements = mp->patch(pn).basis().asElements(refBoxes, refExtension);
        // mp->patch(pn).refineElements( elements );
    }
}

template<class T>
void gsAdaptiveMeshing<T>::_unrefineMarkedElements(   gsFunctionSet<T> * input,
                                                        const std::vector<bool> & elMarked,
                                                        index_t extension)
{
    const short_t dim = input->domainDim();

    // numMarked: Number of marked cells on current patch, also currently marked cell
    // poffset  : offset index for the first element on a patch
    // globalCount: counter for the current global element index
    index_t numMarked, poffset = 0, globalCount = 0;

    // refBoxes: contains marked boxes on a given patch
    gsMatrix<T> refBoxes;

    gsBasis<T> * basis = nullptr;
    for (index_t pn=0; pn < input->nPieces(); ++pn )// for all patches
    {
        if ( gsMultiPatch<T> * mp = dynamic_cast<gsMultiPatch<T>*>(input) ) basis = &(mp->basis(pn));
        if ( gsMultiBasis<T> * mb = dynamic_cast<gsMultiBasis<T>*>(input) ) basis = &(mb->basis(pn));
        GISMO_ASSERT(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");
        // Get number of elements to be refined on this patch
        const size_t numEl = basis->numElements();
        numMarked = std::count_if(elMarked.begin() + poffset,
                                  elMarked.begin() + poffset + numEl,
                                  GS_BIND2ND(std::equal_to<bool>(), true) );
        poffset += numEl;
        refBoxes.resize(dim, 2*numMarked);
        numMarked = 0;// counting current patch element to be refined

        // for all elements in patch pn
        typename gsBasis<T>::domainIter domIt = basis->makeDomainIterator();
        for (; domIt->good(); domIt->next())
        {
            if( elMarked[ globalCount++ ] ) // refine this element ?
            {
                // Construct degenerate box by setting both
                // corners equal to the center
                refBoxes.col(2*numMarked  ) = domIt->lowerCorner();
                refBoxes.col(2*numMarked+1) = domIt->upperCorner();

                // Advance marked cells counter
                numMarked++;
            }
        }

        gsMultiPatch<T> * mp;
        gsMultiBasis<T> * mb;
        if ((mp = dynamic_cast<gsMultiPatch<T>*>(input)))
        {
            // Refine all of the found refBoxes in this patch
            std::vector<index_t> elements = mp->patch(pn).basis().asElementsUnrefine(refBoxes, extension);
            mp->patch(pn).unrefineElements( elements );
        }
        else if ((mb = dynamic_cast<gsMultiBasis<T>*>(input)))
        {
            // Refine all of the found refBoxes in this patch
            mb->unrefine( pn, refBoxes, extension );
        }
        else
            GISMO_ERROR("No gsMultiPatch or gsMultiBasis found");
    }
}

template<class T>
void gsAdaptiveMeshing<T>::_processMarkedElements(gsFunctionSet<T> * input,
                                                    const std::vector<bool> & elRefined,
                                                    const std::vector<bool> & elCoarsened,
                                                    index_t refExtension,
                                                    index_t crsExtension)
{
    const short_t dim = input->domainDim();

    // numMarked: Number of marked cells on current patch, also currently marked cell
    // poffset  : offset index for the first element on a patch
    // globalCount: counter for the current global element index
    index_t numRefined, numCoarsened, poffset = 0, globalCount = 0;

    // refBoxes: contains marked boxes on a given patch
    gsMatrix<T> refBoxes, crsBoxes;

    gsBasis<T> * basis = nullptr;
    for (index_t pn=0; pn < input->nPieces(); ++pn )// for all patches
    {
        if ( gsMultiPatch<T> * mp = dynamic_cast<gsMultiPatch<T>*>(input) ) basis = &(mp->basis(pn));
        if ( gsMultiBasis<T> * mb = dynamic_cast<gsMultiBasis<T>*>(input) ) basis = &(mb->basis(pn));
        GISMO_ASSERT(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");
        // Get number of elements to be refined on this patch
        const size_t numEl = basis->numElements();
        numRefined = std::count_if(elRefined.begin() + poffset,
                                  elRefined.begin() + poffset + numEl,
                                  GS_BIND2ND(std::equal_to<bool>(), true) );

        numCoarsened = std::count_if(elCoarsened.begin() + poffset,
                                  elCoarsened.begin() + poffset + numEl,
                                  GS_BIND2ND(std::equal_to<bool>(), true) );


        poffset += numEl;

        refBoxes.resize(dim, 2*numRefined);
        crsBoxes.resize(dim, 2*numCoarsened);
        numRefined = 0;// counting current patch element to be refined
        numCoarsened = 0;// counting current patch element to be refined

        // for all elements in patch pn
        typename gsBasis<T>::domainIter domIt = basis->makeDomainIterator();
        for (; domIt->good(); domIt->next())
        {
            if( elRefined[ globalCount ] ) // refine this element ?
            {
                // Construct degenerate box by setting both
                // corners equal to the center
                refBoxes.col(2*numRefined  ) = domIt->lowerCorner();
                refBoxes.col(2*numRefined+1) = domIt->upperCorner();

                // Advance marked cells counter
                numRefined++;
            }
            if( elCoarsened[ globalCount ] ) // refine this element ?
            {
                // Construct degenerate box by setting both
                // corners equal to the center
                crsBoxes.col(2*numCoarsened  ) = domIt->lowerCorner();
                crsBoxes.col(2*numCoarsened+1) = domIt->upperCorner();

                // Advance marked cells counter
                numCoarsened++;
            }

            globalCount++;
        }

        gsMultiPatch<T> * mp;
        gsMultiBasis<T> * mb;
        if ((mp = dynamic_cast<gsMultiPatch<T>*>(input)))
        {
            std::vector<index_t> elements;
            // Unrefine all of the found refBoxes in this patch
            elements = mp->patch(pn).basis().asElementsUnrefine(crsBoxes, crsExtension);
            mp->patch(pn).unrefineElements( elements );

            // Refine all of the found refBoxes in this patch
            elements = mp->patch(pn).basis().asElements(refBoxes, refExtension);
            mp->patch(pn).refineElements( elements );
        }
        else if ((mb = dynamic_cast<gsMultiBasis<T>*>(input)))
        {
            // Refine all of the found refBoxes in this patch
            mb->unrefine( pn, crsBoxes, crsExtension);
            mb->refine( pn, refBoxes, refExtension );
        }
        else
            GISMO_ERROR("No gsMultiPatch or gsMultiBasis found");
    }
}

template<class T>
void gsAdaptiveMeshing<T>::_flattenElements(  gsFunctionSet<T> * input,
                                                const index_t level)
{
    // Get all the elements of which the level exceeds level
    std::vector<bool> elMarked;
    _markLevelThreshold(input,level,elMarked);

    const int dim = input->domainDim();

    // numMarked: Number of marked cells on current patch, also currently marked cell
    // poffset  : offset index for the first element on a patch
    // globalCount: counter for the current global element index
    int numMarked, poffset = 0, globalCount = 0;

    // crsBoxes: contains marked boxes on a given patch
    gsMatrix<T> crsBoxes;

    gsBasis<T> * basis = nullptr;

    for (index_t pn=0; pn < input->nPieces(); ++pn )// for all patches
    {
        if ( gsMultiPatch<T> * mp = dynamic_cast<gsMultiPatch<T>*>(input) ) basis = &(mp->basis(pn));
        if ( gsMultiBasis<T> * mb = dynamic_cast<gsMultiBasis<T>*>(input) ) basis = &(mb->basis(pn));
        GISMO_ASSERT(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");
        // Get number of elements to be refined on this patch
        const int numEl = basis->numElements();
        numMarked = std::count_if(elMarked.begin() + poffset,
                                  elMarked.begin() + poffset + numEl,
                                  GS_BIND2ND(std::equal_to<bool>(), true) );

        poffset += numEl;
        crsBoxes.resize(dim, 2*numMarked);
        numMarked = 0;// counting current patch element to be refined

        // for all elements in patch pn
        typename gsBasis<T>::domainIter domIt = basis->makeDomainIterator();
        for (; domIt->good(); domIt->next())
        {
            if( elMarked[ globalCount++ ] ) // refine this element ?
            {
                // Construct degenerate box by setting both
                // corners equal to the center
                crsBoxes.col(2*numMarked  ) =
                        crsBoxes.col(2*numMarked+1) = domIt->centerPoint();

                // Advance marked cells counter
                numMarked++;
            }
        }

        std::vector<index_t> elements;
        elements = basis->asElementsUnrefine(crsBoxes,0);
        const int offset = 2*dim+1;
        index_t diff = 0;
        GISMO_ASSERT(elements.size()%offset==0,"Element boxes have wrong size. Boxes should have size "<<offset<<" per box");
        for (size_t k = 0; k!=elements.size()/offset; k++)
        {
            // gsDebug<<"Old Box:\n";
            // for (index_t kk = 0; kk!=2*dim+1; kk++)
            // {
            //     gsDebug<<elements[offset*k+kk]<<",";
            // }
            // gsDebug<<"\n";

            if ((diff = elements[offset*k] - level) > 0)
            {
                elements[offset*k] = level;
                for (index_t kk = 0; kk!=2*dim; kk++)
                {

                    elements[offset*k+kk+1] = elements[offset*k+kk+1] >> diff;

                    // // yes, exactly: if   i is index in old level a and needs to become a level b box, b<a, then the new index j is (in c++ computation):
                    // j = i >> a-b;
                    // // this is division by 2^{a-b}  using bit-operations
                }
            }
            // gsDebug<<"New Box:\n";
            // for (index_t kk = 0; kk!=2*dim+1; kk++)
            // {
            //     gsDebug<<elements[offset*k+kk]<<",";
            // }
            // gsDebug<<"\n";
        }
        gsMultiPatch<T> * mp;
        gsMultiBasis<T> * mb;
        if ((mp = dynamic_cast<gsMultiPatch<T>*>(input)))
        {
            mp->patch(pn).unrefineElements( elements );
        }
        else if ((mb = dynamic_cast<gsMultiBasis<T>*>(input)))
        {
            // Refine all of the found refBoxes in this patch
            mb->unrefineElements(pn, elements );
        }
        else
            GISMO_ERROR("No gsMultiPatch or gsMultiBasis found");
    }
}

} // namespace gismo
