/** @file BaseObject.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gsCore/gsFunctionSet.h>
#include <gsNewExpressions/SpaceObject.h>

#pragma once

namespace gismo
{
namespace Expr
{

template <class T, enum SpaceType _Space, size_t _Order>
struct ExpressionTraits<SolutionObject<T, _Space, _Order>>
{
    typedef T Scalar;
    static constexpr size_t Order = _Order;
    static constexpr SpaceType Space = SpaceType::None; // the solution has SpaceType::None, its underlying space has _Space
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = false;
};


template <class T, enum SpaceType _Space, size_t _Order>
class SolutionObject : public BaseObject<SolutionObject<T, _Space, _Order>>
{
    using Base = BaseObject<SolutionObject<T, _Space, _Order>>;

public:
    // Expose the static traits publicly
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    typedef typename Base::Scalar Scalar;

private:
    const SpaceObject<T,_Space,_Order> * m_Space;
    gsMatrix<Scalar> * m_solVector;
public:
    SolutionObject(const SpaceObject<T,_Space,_Order> & space,
                gsMatrix<Scalar> & solVector,
                std::string label="u")
    :
    Base(space.domainDim(),space.sizes(), label),
    m_Space(&space),
    m_solVector(&solVector)
    {
    }

    ExpressionResult<Scalar> eval(const index_t k) const
    {
        // TODO: SolutionObject needs proper implementation
        // For now, return placeholder
        ExpressionResult<Scalar> result(1, 1);
        result(0, 0) = gsMatrix<Scalar>();
        return result;
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
    }

    const SpaceObject<T,_Space,_Order> & getSpace() const
    {
        // gsInfo<<"m_Space.Space = "<<m_Space->space<<"\n";
        // gsInfo<<"m_Space.Order = "<<m_Space->order<<"\n";
        return *m_Space;
    }

    void print(std::ostream & os) const
    {
        os<<Base::label_;
    }
    // void print(std::ostream & os) const
    // {
    //     os<<"Solution: "<<Base::label_
    //       <<", solVector size "<<m_solVector->rows()<<"x"<<m_solVector->cols()<<")"
    //       <<", space (id="<<m_Space->id()<<"): "<<getSpace();
    // }

};

template <class T, enum SpaceType Sol_Space, size_t Sol_Order, typename SpaceObj>
auto variation(const SolutionObject<T,Sol_Space,Sol_Order> & sol,const SpaceObj & space)
-> VariationObject<SolutionObject<T,Sol_Space,Sol_Order>,SpaceObj>
{
    return VariationObject<SolutionObject<T,Sol_Space,Sol_Order>,SpaceObj>(sol,space, sol.getSpace().id() == space.id());
}

}//namespace Expr
}//namespace gismo