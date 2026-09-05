/** @file gsHSplinesParaview.hpp

    @brief Paraview export of gsHSplines types (moved from gsIO/gsWriteParaview,
    modularization stream S3 step A8: type-specific visualization lives
    with the type's module, the base IO module stays type-blind).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsHSplines/gsHSplinesParaview.h>
#include <gsIO/gsParaviewCollection.h>

#define PLOT_PRECISION 12 // mirrors gsWriteParaview.hpp

// NOTE: deliberately NOT including gsIO/gsWriteParaview.hpp - helper
// templates are declared in gsWriteParaview.h and their real_t
// instantiations are exported from the gsIO module (a second inclusion
// of that .hpp would poison the exported symbol visibility).

namespace gismo
{

/// Writes a single \ref gsHBox \a box to a file with name \a fn
template<short_t d, class T>
void writeSingleHBox(const gsHBox<d,T> & box, std::string const & fn)
{
    box.computeCoordinates();
    gsVector<index_t,d> np;
    np.setConstant(2);
    gsGridIterator<T,CUBE,d> grid(box.getCoordinates(),np);
    gsMatrix<T> points = grid.toMatrix();
    gsMatrix<T> values(3,points.cols());
    values.row(0).setConstant(box.level());
    values.row(1).setConstant(box.error());
    values.row(2).setConstant(box.projectedErrorRef());
    gsWriteParaviewTPgrid(points,values,np,fn);
}

/// Writes a single \ref gsHBox \a box to a file with name \a fn
template<short_t d, class T>
void gsWriteParaview(const gsHBox<d,T> & box, std::string const & fn, short_t mode)
{
    box.computeCoordinates();
    switch (mode)
    {
        case 1:
            gsWriteParaview(box.getCoordinates(), fn, box.error());
            break;
        case 2:
            gsWriteParaview(box.getCoordinates(), fn, box.projectedErrorRef());
            break;
        default:
            gsWriteParaview(box.getCoordinates(), fn, (T)box.level());
            break;
    }
}

template<short_t d, class T>
void gsWriteParaview(const gsHBoxContainer<d,T> & boxes, std::string const & fn, short_t mode)
{
    gsMatrix<T> boxCoords(d,boxes.totalSize()*2);
    gsVector<T> boxValues(boxes.totalSize());
    boxCoords.setZero();
    index_t i=0;
    std::string fileName;
    for (typename gsHBoxContainer<d,T>::cHIterator Hit = boxes.cbegin(); Hit!=boxes.cend(); Hit++)
        for (typename gsHBoxContainer<d,T>::cIterator Cit = Hit->cbegin(); Cit!=Hit->cend(); Cit++, i++)
        {
            Cit->computeCoordinates();
            boxCoords.middleCols(i*2,2) = Cit->getCoordinates();
            switch (mode)
            {
                case 1:
                    boxValues(i) = Cit->error();
                    break;
                case 2:
                    boxValues(i) = Cit->projectedErrorRef();
                    break;
                default:
                    boxValues(i) = Cit->level();
                    break;
            }
        }
    gsWriteParaview(boxCoords, fn, boxValues);
}

} // namespace gismo

#undef PLOT_PRECISION
