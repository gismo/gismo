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

#include <gsAssembler/gsDofMapperCreator.h>

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
    fs(&_fs), dim(give(_dim)), id(give(_id)), mapperInstalled(false),
    installedFs(nullptr), installedDim(0) { }

    const gsFunctionSet<T> * fs;
    index_t dim, id;
    gsDofMapper mapper;
    gsMatrix<T> fixedDofs;
    index_t cont; //int. coupling

    /// True when \a mapper was installed by gsFeSpace::setupMapper rather than
    /// built by init(). An installed mapper may be ragged (different dof counts
    /// per component), so its size cannot be checked against fs->size()*dim.
    bool mapperInstalled;
    /// Source and component count \a mapper was installed against; the flag is
    /// only honoured while both still hold, since gsExprAssembler::getSpace may
    /// rebind fs/dim of an existing space in place.
    const gsFunctionSet<T> * installedFs;
    index_t installedDim;

    /// True when every component of \a map holds the same number of dofs,
    /// i.e. the layout of every non-ragged mapper.
    static bool uniformComponents(const gsDofMapper & map)
    {
        for (index_t c = 1; c < map.numComponents(); ++c)
            if (map.totalSize(c) != map.totalSize(0)) return false;
        return true;
    }

    bool valid() const
    {
        GISMO_ASSERT(nullptr!=fs, "Invalid pointer.");
        if (mapperInstalled && fs==installedFs && dim==installedDim)
            return true;
        return static_cast<size_t>(fs->size()*dim)==mapper.mapSize();
    }

    void init()
    {
        GISMO_ASSERT(nullptr!=fs, "Invalid pointer.");
        if (mapperInstalled && fs==installedFs && dim==installedDim)
            return;
        if (const gsMultiBasis<T> * mb =
            dynamic_cast<const gsMultiBasis<T>*>(fs) )
            mapper = createMapper(*mb, dim, /*conforming=*/false);
        else if (const gsBasis<T> * b =
                 dynamic_cast<const gsBasis<T>*>(fs) )
            mapper = createMapper(*b, dim, /*conforming=*/false);
        mapper.finalize();
        fixedDofs.clear();
        cont = -1;
        mapperInstalled = false;
        installedFs = nullptr;
    }
};

}// namespace expr
}// namespace gismo