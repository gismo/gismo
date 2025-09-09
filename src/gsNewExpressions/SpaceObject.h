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

template <class T, size_t _space, size_t _order>
class SpaceObject : public BaseObject<SpaceObject<T, _space, _order>>
{
    using Base = BaseObject<SpaceObject<T, _space, _order>>;
    typedef typename Base::Scalar Scalar;
    using Base::order;
    using Base::space;
    using Base::deriv_;

private:
    const gsFunctionSet<Scalar> * m_fs;
    const gsFuncData<Scalar>    * m_fd;
    size_t m_id;
public:
    SpaceObject(size_t domainDim, const std::array<size_t, order> & input_sizes, size_t id,
                std::string label=(space==Space::Test)?"φ":(space==Space::Trial)?"ψ":"UNKNOWN_SPACE")
    :
    Base(domainDim, input_sizes, label),
    m_fs(NULL), m_fd(NULL), m_id(id)
    {
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        return m_fd->values[0].col(k).blockDiag(order);
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
        if (deriv_ > 0)
            m_fd->flags |= NEED_DERIV;
        if (deriv_ > 1)
            m_fd->flags |= NEED_DERIV2;
    }

    const SpaceObject<Scalar,_space,_order> & rowVar() const {return (_space==Space::Test  || _space==Space::Both) ? *this : NullObject<T,order>::get();}
    const SpaceObject<Scalar,_space,_order> & colVar() const {return (_space==Space::Trial || _space==Space::Both) ? *this : NullObject<T,order>::get();}

    void print(std::ostream & os) const
    {
        os<<Base::label_;
        _print_arguments(os);
        // _print_impl<space>(os);
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