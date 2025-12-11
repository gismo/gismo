/** @file BaseObject.h

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

template <class T, size_t _Space, size_t _Order>
struct ExpressionTraits<SpaceObject<T, _Space, _Order>>
{
    typedef T Scalar;
    static constexpr size_t Order = _Order;
    static constexpr size_t Space = _Space;
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = false;
};

template <class T, size_t _Space, size_t _Order>
class SpaceObject : public BaseObject<SpaceObject<T, _Space, _Order>>
{
    static_assert(_Space==SpaceType::None || _Space==SpaceType::Test || _Space==SpaceType::Trial, "SpaceObject can only be None, Test or Trial space.");

    using Base = BaseObject<SpaceObject<T, _Space, _Order>>;
    using Base::Deriv_;

public:
    // Expose the static traits publicly
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    typedef typename Base::Scalar Scalar;

private:
    const gsFunctionSet<Scalar> * m_fs;
    const gsFuncData<Scalar>    * m_fd;
    size_t m_id;

public:
    SpaceObject(size_t domainDim, const std::array<size_t, Order> & input_sizes, size_t id,
                std::string label=(Space==SpaceType::Test)?"φ":(Space==SpaceType::Trial)?"ψ":"UNKNOWN_SPACE")
    :
    Base(domainDim, input_sizes, label),
    m_fs(NULL), m_fd(NULL), m_id(id)
    {
    }

    // Specialization for Test space using enable_if
    template <size_t S = _Space>
    typename std::enable_if<S == SpaceType::Test, ExpressionValue<Scalar>>::type
    eval(const index_t k) const
    {
        // Get the number of active basis functions
        index_t numActive = m_fd->values[0].rows();
        
        // Test space: cardinality is (numActive, 1) - no runtime branching
        ExpressionValue<Scalar> result(numActive, 1);
        
        // Fill in the values for each basis function
        for (index_t i = 0; i < numActive; ++i)
        {
            result(i, 0) = m_fd->values[0].col(k).segment(i, 1).blockDiag(Order ? this->sizes()[0] : 1);
        }
        
        return result;
    }

    // Specialization for Trial space using enable_if
    template <size_t S = _Space>
    typename std::enable_if<S == SpaceType::Trial, ExpressionValue<Scalar>>::type
    eval(const index_t k) const
    {
        // Get the number of active basis functions
        index_t numActive = m_fd->values[0].rows();
        
        // Trial space: cardinality is (1, numActive) - no runtime branching
        ExpressionValue<Scalar> result(1, numActive);
        
        // Fill in the values for each basis function
        for (index_t i = 0; i < numActive; ++i)
        {
            result(0, i) = m_fd->values[0].col(k).segment(i, 1).blockDiag(Order ? this->sizes()[0] : 1);
        }
        
        return result;
    }

    // Specialization for None space using enable_if
    template <size_t S = _Space>
    typename std::enable_if<S == SpaceType::None, ExpressionValue<Scalar>>::type
    eval(const index_t k) const
    {
        // Get the number of active basis functions
        index_t numActive = m_fd->values[0].rows();
        
        // None space: cardinality is (1, 1) - no runtime branching
        ExpressionValue<Scalar> result(1, 1);
        
        // Fill in the values (single entry)
        for (index_t i = 0; i < numActive; ++i)
        {
            result(0, 0) = m_fd->values[0].col(k).segment(i, 1).blockDiag(Order ? this->sizes()[0] : 1);
        }
        
        return result;
    }

    const gsFunctionSet<T> & source() const {return *m_fs;}
    const gsFuncData<T>    & data()   const
    {
        GISMO_ASSERT(NULL!=m_fd, "SpaceObject: invalid data "<< m_fs <<","<<m_fd);
        return *m_fd;
    }

    void setSource(const gsFunctionSet<Scalar> & fs) { m_fs = &fs;}
    void setData(const gsFuncData<Scalar> & val) { m_fd = &val;}

    size_t id() const { return m_id; }

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
        os<<"(";
        for (size_t d=0; d!=this->domainDim(); d++)
        {
            os<<"x"<<d;
            if (d!=this->domainDim()-1)
                os<<",";
        }
        os<<")";
    }
};

}//namespace Expr
}//namespace gismo

