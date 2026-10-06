/** @file gsComposition.h

    @brief Defines an expression for a composed variable

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
 * @brief Expression for the composition of a variable and a geometry map
 * @ingroup Expressions
 * @tparam T The expression type
 */
template<class T>
class gsComposition : public symbol_expr< gsComposition<T> >
{ //comp(f,G)
    friend class gismo::gsExprHelper<T>;
    typedef symbol_expr< gsComposition<T> > Base;
    typename gsGeometryMap<T>::Nested_t _G;
protected:
    explicit gsComposition(const gsGeometryMap<T> & G, index_t _d = 1)
    : Base(_d), _G(G) { }
public:
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    typename gsMatrix<T>::constColumn
    eval(const index_t k) const { return this->m_fd->values[0].col(k); }

    const gsGeometryMap<T> & inner() const { return _G;};

    void parse(gsExprHelper<T> & evList) const
    {
        //evList.add(_G); //done in gsExprHelper
        evList.add(*this);
        this->data().flags |= NEED_VALUE|NEED_ACTIVE;
        //_G.data().flags  |= NEED_VALUE; //done in gsExprHelper
    }
};


}// namespace expr
}// namespace gismo