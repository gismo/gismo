/** @file FeSpaceData.h

    @brief Defines a data structure for a finite element space in the new expression module.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsCore/gsDofMapper.h>
#include <gsCore/gsFunction.h>
#include <gsCore/gsMultiBasis.h>

namespace gismo
{
namespace Expr
{

/**
 * @brief Struct containing information for matrix assembly for the new expression module.
 * @ingroup gsNewExpressions
 */
template<class T>
struct FeSpaceData
{
    FeSpaceData(const gsFunctionSet<T>& basis, index_t dimension, index_t space_id)
        : fs(&basis), dim(dimension), id(space_id), cont(-1), initialized(false) {}

    void init()
    {
        GISMO_ASSERT(fs != nullptr, "FunctionSet pointer is null.");
        if (const gsMultiBasis<T>* mb = dynamic_cast<const gsMultiBasis<T>*>(fs))
        {
            mapper = gsDofMapper(*mb, dim);
        }
        else if (const gsBasis<T>* b = dynamic_cast<const gsBasis<T>*>(fs))
        {
            mapper = gsDofMapper(*b, dim);
        }
        else
        {
            GISMO_ERROR("FeSpaceData can only be initialized from gsBasis or gsMultiBasis.");
        }
        mapper.finalize();
        fixedDofs.setZero(mapper.boundarySize(), 1); // Initialize fixedDofs vector
        initialized = true;
    }

    bool isInitialized() const
    {
        return initialized && mapper.isFinalized();
    }

    const gsFunctionSet<T>* fs;
    index_t dim;
    index_t id;
    index_t cont;           // Interface continuity (-1: default, 0: C0 conforming)
    bool initialized;       // Whether space has been set up
    gsDofMapper mapper;
    gsMatrix<T> fixedDofs;
};

} // namespace Expr
} // namespace gismo
