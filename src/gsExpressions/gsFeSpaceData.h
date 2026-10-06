/** @file gsFeSpaceData.h

    @brief Defines a data structure for the gsFeSpace

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
               H.M. Verhelst
*/

#pragma once

namespace gismo
{
namespace expr
{

/**
 * @brief Struct containing information for matrix assembly
 * @ingroup Expressions
 * @tparam T The expression type
 */
template<class T>
struct gsFeSpaceData
{
    gsFeSpaceData(const gsFunctionSet<T> & _fs, index_t _dim, index_t _id):
    fs(&_fs), dim(give(_dim)), id(give(_id)) { }

    const gsFunctionSet<T> * fs;
    index_t dim, id;
    gsDofMapper mapper;
    gsMatrix<T> fixedDofs;
    index_t cont; //int. coupling

    bool valid() const
    {
        GISMO_ASSERT(nullptr!=fs, "Invalid pointer.");
        return static_cast<size_t>(fs->size()*dim)==mapper.mapSize();
    }

    void init()
    {
        GISMO_ASSERT(nullptr!=fs, "Invalid pointer.");
        if (const gsMultiBasis<T> * mb =
            dynamic_cast<const gsMultiBasis<T>*>(fs) )
            mapper = gsDofMapper(*mb, dim );
        else if (const gsBasis<T> * b =
                 dynamic_cast<const gsBasis<T>*>(fs) )
            mapper = gsDofMapper(*b, dim );
        mapper.finalize();
        fixedDofs.clear();
        cont = -1;
    }
};

}// namespace expr
}// namespace gismo