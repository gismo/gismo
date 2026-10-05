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
    fs(&_fs), dim(give(_dim)), id(give(_id)), cont(-1),
    mapperStorage(gsDofMapper::storage::dense) { }

    const gsFunctionSet<T> * fs;
    index_t dim, id;
    gsDofMapper mapper;
    gsMatrix<T> fixedDofs;
    index_t cont; //int. coupling
    gsDofMapper::storage mapperStorage; ///< storage mode of the mappers built by init()

    bool valid() const
    {
        GISMO_ASSERT(nullptr!=fs, "Invalid pointer.");
        return static_cast<size_t>(fs->size()*dim)==mapper.mapSize();
    }

    /// Rejects, in every build type, a mapper that the expression
    /// evaluator cannot consume.  The evaluator holds one source basis
    /// per space and indexes every component with that basis' actives,
    /// so a mapper whose components were built from distinct bases, or
    /// whose components differ in size on some patch, would be indexed
    /// with the wrong local basis functions.
    ///
    /// The decision is gsDofMapper::hasUniformComponents(), which
    /// relies on the declared hasDistinctComponentSpaces() rather than
    /// on sizes: a Raviart-Thomas pair on an isotropic mesh has equal
    /// component sizes on every patch.  Mappers from the single-basis
    /// creators and default-constructed ones always pass.
    static void ensureUsableByUniformEvaluator(const gsDofMapper & mapper)
    {
        if (mapper.hasUniformComponents())
            return;

        GISMO_ENSURE(!mapper.hasDistinctComponentSpaces(),
                     "The dof-mapper was built from distinct per-component bases "
                     "(hasDistinctComponentSpaces()), which the expression "
                     "assembler cannot evaluate: it uses one basis for all "
                     "components of a space.");

        const index_t nComp = mapper.numComponents();
        if (gsDofMapper::GlobalIdentity == mapper.layout())
        {
            for (index_t c = 1; c != nComp; ++c)
                GISMO_ENSURE(mapper.totalSize(c) == mapper.totalSize(0),
                             "The dof-mapper has components of unequal size, which "
                             "the expression assembler cannot evaluate: component "
                             <<c<<" has "<<mapper.totalSize(c)<<" dofs, component 0 has "
                             <<mapper.totalSize(0)<<".");
        }
        else
        {
            for (index_t p = 0; p != static_cast<index_t>(mapper.numPatches()); ++p)
                for (index_t c = 1; c != nComp; ++c)
                    GISMO_ENSURE(mapper.patchSize(p,c) == mapper.patchSize(p,0),
                                 "The dof-mapper has components of unequal size, which "
                                 "the expression assembler cannot evaluate: component "
                                 <<c<<" has "<<mapper.patchSize(p,c)<<" dofs on patch "<<p
                                 <<", component 0 has "<<mapper.patchSize(p,0)<<".");
        }
        GISMO_ERROR("The dof-mapper cannot be evaluated by the expression assembler.");
    }

    void init()
    {
        GISMO_ASSERT(nullptr!=fs, "Invalid pointer.");
        if (const gsMultiBasis<T> * mb =
            dynamic_cast<const gsMultiBasis<T>*>(fs) )
            mapper = createMapper(*mb, dim, /*conforming=*/false,
                                  /*finalize=*/false, mapperStorage);
        else if (const gsBasis<T> * b =
                 dynamic_cast<const gsBasis<T>*>(fs) )
            mapper = createMapper(*b, dim, /*conforming=*/false,
                                  /*finalize=*/false, mapperStorage);
        mapper.finalize();
        fixedDofs.clear();
        cont = -1;
    }
};

}// namespace expr
}// namespace gismo