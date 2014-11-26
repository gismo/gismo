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
#include <gismo.h>

using namespace std;
using namespace gismo;


// DOCUMENTATION WILL FOLLOW SHORTLY

/** @file Contains some auxiliary functions for cell marking
 */

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
 * <b>refCriterion = 1</b>:\n
 * Threshold computed based on the largest of all appearing local errors:
 * \f[ \Theta = \rho \cdot \max_K \{ \eta_K \} \f]
 * The actual number of marked elements can vary in each refinement step,
 * depending on the distribution of the error.
 *
 * <b>refCriterion = 2</b>:\n
 * In each step, a certain percentage of all elements are marked.
 * \f[ \Theta = (1-\rho)\cdot 100\ \textrm{-percentile of}\ \{ \eta_K \}_K \f]
 * For example, if \f$\rho = 0.8\f$, those 20% of all elements which have the
 * largest local errors are marked for refinement.
 *
 * <b>refCriterion = 3</b>:\n
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
 */
template <class T>
void gsMarkElementsForRef( const std::vector<T> & elError, int refCriterion, T refParameter, std::vector<bool> & elMarked)
{
    T Thr = T(0);

    if( refCriterion == 1 )
    {
        // First, conduct a brutal search for the maximum local error
        T maxErr = 0;
        for( unsigned i = 0; i < elError.size(); ++i)
            if( maxErr < elError[i] )
                 maxErr = elError[i];

        // Compute the threshold:
        Thr = refParameter * maxErr;
    }
    else if ( refCriterion == 2 )
    {
        // Total number of elements:
        unsigned NE = elError.size();
        // The vector of local errors will need to be sorted,
        // which will be done on a copy:
        std::vector<T> elErrCopy = elError;

        // Compute the index from which the refinement should start,
        // once the vector is sorted.
        unsigned idxRefineStart = static_cast<unsigned>( floor( refParameter * T(NE) ) );
        // ...and just to be sure we are in range:
        if( idxRefineStart == elErrCopy.size() )
            idxRefineStart -= 1;

       // Sort the list using bubblesort.
       // After each loop, the largest elements are at the end
       // of the list. Since we are only interested in the largest elements,
       // it is enough to run the sorting until enough "largest" elements
       // have been found, i.e., until we have reached indexRefineStart
       unsigned lastSwapDone = elErrCopy.size() - 1;
       unsigned lastCheckIdx = lastSwapDone;

       bool didSwap;
       T tmp;
       do{
           didSwap = false;
           lastCheckIdx = lastSwapDone;
           for( unsigned i=0; i < lastCheckIdx; i++)
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
    }
    else if( refCriterion == 3 )
    {
        // The vector of local errors will need to be sorted,
        // which will be done on a copy:
        std::vector<T> elErrCopy = elError;

        // Compute the sum, i.e., the global/total error
        T totalError = T(0);
        for( unsigned i = 0; i < elErrCopy.size(); ++i)
            totalError += elErrCopy[i];

        // We want to mark just enough cells such that their
        // cummulated errors add up to a certain fraction
        // of the total error.
        T errorMarkSum = (1-refParameter) * totalError;
        T cummulErrMarked = 0;

        T tmp;
        unsigned lastSwapDone = elErrCopy.size() - 1;
        do{
            for( unsigned i=0; i < lastSwapDone; i++)
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
    }

    elMarked.resize( elError.size() );

    // Now just check for each element, whether the local error
    // is above the computed threshold or not, and mark accordingly.
    for( unsigned i=0; i < elError.size(); i++)
        ( elError[i] >= Thr ? elMarked[i] = true : elMarked[i] = false );

} // gsMarkCells


/** \brief Refine a gsMultiBasis, based on a vector of element-markings.
 *
 * Given the vector of element-markings, the corresponding element
 * in the mesh underlying \em basis is refined.
 *
 * \remarks
 * The order/numbering of the elements is implicitly defined by
 * the numbering of the patches in gsMultiBasis, and
 * by the gsDomainIterator of the respective patch-wise basis!
 *
 * \param basis gsMultiBasis to be refined adaptively.
 * \param elMarked std::vector of Booleans indicating
 * for each element of the mesh underlying \em basis, whether it should be refined or not.
 *
 *
 */
template <class T>
void gsRefineMarkedElements( gsMultiBasis<T> & basis, std::vector<bool> & elMarked)
{
    int globalCount = 0;

    // refBoxes will contain a gsMatrix for each of the marked elements.
    std::vector< gsMatrix<T> > refBoxes;

    // Collect the coordinates of all elements of the gsMultiBasis which are marked.
    for (unsigned pn=0; pn < basis.nBases(); ++pn )// for all patches
    {
        typename gsBasis<T>::domainIter domIt = basis.basis(pn).makeDomainIterator();
        for (; domIt->good(); domIt->next())
        {
            if( elMarked[ globalCount ] )
            {
                //gsVector<T> ctr = domIt->centerPoint();
                gsVector<T> low = domIt->lowerCorner();
                gsVector<T> upp = domIt->upperCorner();

                // The refBox is not given by the actual corners of the
                // marked cells, but in the following form due to
                // implementational reasons.
                // In the refinement using the function gsBasis::refine(gsMatrix),
                // the function uniqueFindSpan is used. This causes problems,
                // if some of the corners concide with knot lines.
                gsMatrix<T> refBox( low.size(), 2 );
                for( unsigned i=0; i < low.size(); ++i )
                {
                refBox(i,0) = 0.75 * low[i] + 0.25 * upp[i];
                refBox(i,1) = 0.25 * low[i] + 0.75 * upp[i];
                }
                refBoxes.push_back( refBox );
            }
            globalCount += 1;
        }

        //std::cout << "Refining " << refBoxes.size() << " elements" << std::endl;

        // Refine all of the found refBoxes.
        for( unsigned i = 0; i < refBoxes.size(); i++ )
            basis.refine( pn, refBoxes[i] );
    }
}

}
