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

#include <gismo.h>

using namespace std;
using namespace gismo;


// DOCUMENTATION WILL FOLLOW SHORTLY

namespace gismo
{


template <class T>
void gsMarkCells( const std::vector<T> & elError, int refCriterion, T refParameter, std::vector<bool> & elMarked)
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
        unsigned NE = elError.size();
        std::vector<T> elErrCopy = elError;

        unsigned idxRefineStart = static_cast<unsigned>( floor( refParameter * T(NE) ) );
        if( idxRefineStart == elErrCopy.size() )
            idxRefineStart -= 1;

       // Sort the list using bubblesort.
       // After each loop, the largest elements are at the end
       // of the list. Since we are only interested in the largest elements,
       // it is enough to run the sorting until enough "largest" elements
       // have been found.
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

       Thr = elErrCopy[ idxRefineStart ];
    }
    else if( refCriterion == 3 )
    {
        std::vector<T> elErrCopy = elError;

        T totalError = T(0);
        for( unsigned i = 0; i < elErrCopy.size(); ++i)
            totalError += elErrCopy[i];

        T errorReduce = (1-refParameter) * totalError;
        T cummulErrRed = 0;

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

            cummulErrRed += elErrCopy[ lastSwapDone  ];
            lastSwapDone -= 1;

        }while( cummulErrRed < errorReduce && lastSwapDone > 0 );

        Thr = elErrCopy[ lastSwapDone + 1 ];
    }

    elMarked.resize( elError.size() );

    for( unsigned i=0; i < elError.size(); i++)
        ( elError[i] >= Thr ? elMarked[i] = true : elMarked[i] = false );

} // gsMarkCells

template <class T>
void gsRefineMarkedCells( gsMultiBasis<T> & basis, std::vector<bool> & elMarked)
{


std::cout << "ARU 130, nBases = " << basis.nBases() << std::endl;
std::cout << "basis(0) = " << basis[0] << std::endl;

    int globalCount = 0;

    std::vector< gsMatrix<T> > refBoxes;

    for (unsigned pn=0; pn < basis.nBases(); ++pn )// for all patches
    {
        typename gsBasis<T>::domainIter domIt = basis.basis(pn).makeDomainIterator();
        for (; domIt->good(); domIt->next())
        {
            if( elMarked[ globalCount ] )
            {
                gsVector<T> ctr = domIt->centerPoint();
                gsVector<T> low = domIt->lowerCorner();
                gsVector<T> upp = domIt->upperCorner();
                gsMatrix<T> refBox( ctr.size(), 2 );
                for( unsigned i=0; i < ctr.size(); ++i )
                {
                refBox(i,0) = 0.5 * low[i] + 0.5 * ctr[i];
                refBox(i,1) = 0.5 * ctr[i] + 0.5 * upp[i];
                }
                refBoxes.push_back( refBox );
            }
            globalCount += 1;
        }


        std::cout << " number of boxes: " << refBoxes.size() << std::endl;

        for( unsigned i = 0; i < refBoxes.size(); i++ )
        {
            //std::cout  << "-------- calling refine " << std::endl << refBoxes[i] << std::endl;
            basis.refine( pn, refBoxes[i] );

        }


    }


    std::cout << "\nARU 156, nBases = " << basis.nBases() << std::endl;
    std::cout << "basis(0) = " << basis[0] << std::endl;

}












}
