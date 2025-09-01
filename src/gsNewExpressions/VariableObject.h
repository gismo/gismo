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

template <class T, size_t _order, bool _isConstant>
struct ExpressionTraits<VariableObject<T, _order, _isConstant>>
{
    typedef T Scalar;
    static constexpr size_t order = _order;
    static constexpr size_t space = Space::None;
    static constexpr size_t deriv = 0;
    static constexpr bool isConstant = _isConstant;
};

template <class T, size_t _order, bool _isConstant = false>
class VariableObject : public BaseObject<T, _order, _isConstant, Space::None>
{
    using Base = BaseObject<T, _order, _isConstant, 0>;
    using Base::deriv_;
public:

    typedef typename ExpressionTraits<ConstantObject<T, _order>>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<ConstantObject<T, _order>>::order;
    static constexpr size_t space = ExpressionTraits<ConstantObject<T, _order>>::space;
    static constexpr size_t deriv = ExpressionTraits<ConstantObject<T, _order>>::deriv;
    static constexpr bool isConstant = ExpressionTraits<ConstantObject<T, _order>>::isConstant;

private:
    const gsFunctionSet<Scalar> * m_fs;
    const gsFuncData<Scalar>    * m_fd;
public:
    // Default constructor for Eigen compatibility
    VariableObject() : BaseObject<Scalar, _order, _isConstant, 0>(0, std::array<size_t, order>{}), m_fs(NULL), m_fd(NULL) {}

    // Constructor for zero initialization (used by BlockDiag)
    VariableObject(int zero_value) : BaseObject<Scalar, _order, _isConstant, 0>(0, std::array<size_t, order>{}), m_fs(NULL), m_fd(NULL)
    {
        static_cast<void>(zero_value); // Suppress unused parameter warning
    }

    VariableObject(size_t domainDim, const std::array<size_t, order> & input_sizes)
    :
    BaseObject<Scalar, _order, _isConstant, 0>(domainDim, input_sizes),
    m_fs(NULL), m_fd(NULL)
    {
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        return m_fd->values[0].col(k).blockDiag(order+1);
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
        if (deriv_ > 0)
            m_fd->flags |= NEED_DERIV;
        if (deriv_ > 1)
            m_fd->flags |= NEED_DERIV2;
    }

    void print(std::ostream & os) const override
    {
        _print_impl<order>(os);
    }

protected:

    void _print_arguments(std::ostream & os) const
    {
        os<<"(";
        for (size_t d=0; d!=this->domainDim(); d++)
        {
            os<<"x"<<d;
            if (d!=this->domainDim()-1)
                os<<",";
        }
        os<<")";
    }

    template <size_t _ORDER>
    typename std::enable_if<_ORDER==0,void>::type
    _print_impl(std::ostream & os) const
    {
        os<<"f";
        _print_arguments(os);
    }

    template <size_t _ORDER>
    typename std::enable_if<_ORDER!=0,void>::type
    _print_impl(std::ostream & os) const
    {
        os<<"F";
        _print_arguments(os);
    }
};
}//namespace Expr
}//namespace gismo