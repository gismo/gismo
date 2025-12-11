/** @file VariationObject.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gsCore/gsFunctionSet.h>

#pragma once

namespace gismo
{
namespace Expr
{

template <class Expr, class SpaceObj>
struct ExpressionTraits<VariationObject<Expr, SpaceObj>>
{
    // Passthrough traits from the wrapped SpaceObj
    typedef typename SpaceObj::Scalar Scalar;
    static constexpr size_t Order = SpaceObj::Order;
    static constexpr size_t Space = SpaceObj::Space; 
    
    // It is not constant because it depends on the basis functions
    static constexpr bool IsConstant = false;
    
    // We start with 0 derivatives (gradient implementation would increment this)
    static constexpr size_t Deriv = 0;
};

template <class Expr, class SpaceObj>
class VariationObject : public BaseObject<VariationObject<Expr, SpaceObj>>
{
    using Base = BaseObject<VariationObject<Expr, SpaceObj>>;
    using Scalar = typename SpaceObj::Scalar;

public:
    // Expose traits
    using Base::Order;
    using Base::Space;

private:
    const SpaceObj & m_space;
    const bool m_matches;

public:
    VariationObject(const Expr & expr, const SpaceObj & space, bool matches)
    : 
    Base(space.domainDim(), space.sizes(), "δ_{" + space.label() + "}" + expr.label()),
    m_space(space),
    m_matches(matches)
    {
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        if (m_matches) {
            m_space.parse(helper);
        }
    }

    template <size_t S = Space>
    typename std::enable_if<S == SpaceType::Test || S == SpaceType::Trial || S == SpaceType::None, ExpressionValue<Scalar>>::type
    eval(const index_t k) const
    {
        if (m_matches) {
            return m_space.eval(k);
        }
        
        // STRUCTURAL ZERO
        // If IDs don't match, we return a Zero.
        // Assuming ExpressionValue has a constructor for scalar zero or default constructor:
        return ExpressionValue<Scalar>(0); 
    }

/*     // 3. STRUCTURAL INFO
    // Exposed so your assembler can skip integration entirely if true
    bool isZero() const { return !m_matches; } */
    
    size_t id() const { return m_space.id(); }

    void print(std::ostream & os) const
    {
        if (m_matches) os << Base::label();
        else os << "0";
    }
};

template <typename Expr, typename MySpaceObj, typename SpaceObj>
VariationObject<VariationObject<Expr, MySpaceObj>, SpaceObj> variation(const VariationObject<Expr, MySpaceObj> & sol,const SpaceObj & space)
{
    GISMO_ERROR("Not implemented.");
}

}//namespace Expr
}//namespace gismo

