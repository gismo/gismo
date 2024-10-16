/** @file gsMapperUtils.h

    @brief utils for manipolating gsWeightMapper objects

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#pragma once

#include <gsMSplines/gsWeightMapper.h>

namespace gismo {

/**
 * @brief reorderMapperTarget
 *        permutes the target indices of mapper according to permutation.
 *        If permutation does not contain all indices then the order of
 *        the remaining indices is preserved and the specified indices
 *        are appended to the end.
 *
 *        If an index is repeated in permutation the behaviour is unspecified
 *        and can change in the future because of implementation changes.

 * @param mapper
 * @param permutation
 * @param permMatrix : if specified, the pointed to matrix is the permutation matrix
 * @return the position of the first specified index
 */

typedef gsEigen::PermutationMatrix<Dynamic,Dynamic,index_t> gsPermutationMatrix;

template <class T>
index_t reorderMapperTarget (gsWeightMapper<T> &mapper, const std::vector<index_t>& permutation, gsPermutationMatrix* permMatrix=NULL);



/**
 * @brief  combineMappers
 *         Given a set of mappers it creates a new mapper that combines all.
 *         It can be used to construct a global mapper for the cartesian product of spaces
 *         or to map both unknown and test space.
 * @param  mappers
 * @param  shifts, is a vector containing the shift of source indices for each space applied in the use of this mapper.
 * @param  shift, if true the target indices of each spaces are shifted (this is the default).
 *         This is the case of representing the cartesian product of spaces.
 *         If false the target indices are preserved; this is needed to share the same mapper
 *         for both unknown and test spaces.
 */
template <class T>
gsWeightMapper<T>* combineMappers (const std::vector<gsWeightMapper<T>*> &mappers, std::vector<index_t> &shifts, bool needShifting=true);

template <class T>
void combineMappers (const std::vector<gsWeightMapper<T>*> &mappers, bool needShifting, std::vector<index_t> &shifts, gsWeightMapper<T>& result);


}

