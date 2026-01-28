/** @file BaseObject.h

@brief

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): H.M.Verhelst
*/

#pragma once

namespace gismo
{
namespace Expr
{

template <class T, size_t _Order, bool _IsConstant>
struct ExpressionTraits<VariableObject<T, _Order, _IsConstant>>
{
    typedef T Scalar;
    static constexpr size_t Order = _Order;
    static constexpr SpaceType Space = SpaceType::None;
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = _IsConstant;
};

template <class T, size_t _Order, bool _IsConstant = false>
class VariableObject : public BaseObject<VariableObject<T, _Order, _IsConstant>>
{
    using Base = BaseObject<VariableObject<T, _Order, _IsConstant>>;
    using Base::Deriv_;

public:
    // Expose the static traits publicly
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    using Base::sizes_;
    typedef typename Base::Scalar Scalar;

private:
    const gsFunctionSet<Scalar> * m_fs;
    const gsFuncData<Scalar>    * m_fd;
public:

    VariableObject(size_t domainDim, const std::array<size_t, Order> & input_sizes, std::string label=(Order==0)?"f":"F")
    :
    Base(domainDim, input_sizes, label),
    m_fs(NULL), m_fd(NULL)
    {
    }

    ExpressionResult<Scalar> eval(const index_t k) const
    {
        ExpressionResult<Scalar> result(1, 1);
        result(0, 0) = m_fd->values[0].col(k)/* .blockDiag(Order+1) */;
        return result;
    }

    const gsFunctionSet<T> & source() const {return *m_fs;}
    const gsFuncData<T>    & data()   const
    {
        GISMO_ASSERT(NULL!=m_fd, "VariableObject: invalid data "<< m_fs <<","<<m_fd);
        return *m_fd;
    }

    void setSource(const gsFunctionSet<Scalar> & fs) { m_fs = &fs;}
    void setData(const gsFuncData<Scalar> & val) { m_fd = &val;}

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        helper.add(*this);
        m_fd->flags |= NEED_VALUE;
        if (Deriv_ > 0)
            m_fd->flags |= NEED_DERIV;
        if (Deriv_ > 1)
            m_fd->flags |= NEED_DERIV2;
    }

    void print(std::ostream & os) const
    {
        os<<Base::label_;
        _print_arguments(os);
    }

protected:

    void _print_arguments(std::ostream & os) const
    {
        GISMO_UNUSED(os);
        // os<<"(";
        // for (size_t d=0; d!=this->domainDim(); d++)
        // {
        //     os<<"x"<<d;
        //     if (d!=this->domainDim()-1)
        //         os<<",";
        // }
        // os<<")";
    }
};
}//namespace Expr
}//namespace gismo