/** @file gsPatchwiseFunction.h

    @brief Provides declaration of a gsPatchwiseFunction class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFunction.h>

namespace gismo
{

/** @brief A function depending on an index \a i, typically referring
    to a patch/sub-domain.

    This class is abstract, derived classes implement it in different
    ways.

    \tparam T argument and value type
    
    \ingroup Core
*/
template <class T>
class gsPatchwiseFunction : public gsFunctionSet<T>
{
public:

    virtual ~gsPatchwiseFunction() { }    

    GISMO_UPTR_FUNCTION_PURE(gsPatchwiseFunction, clone)

    operator const gsFunction<T> & ()
    { 
        return piece(0); 
    }

public:

    virtual const gsFunction<T> & piece(const index_t i) const = 0;

    virtual index_t size() const = 0;

public:

    inline const gsFunction<T> & operator [] (const index_t i) const 
    {
        return piece(i);
    }

    std::ostream &print(std::ostream &os) const
    {
        os << "Patch-wise function.\n";
        return os; 
    }

    friend std::ostream & operator<<(std::ostream &os, const gsPatchwiseFunction &pwf)
    {
        return pwf.print(os);
    }
};

} // namespace gismo

