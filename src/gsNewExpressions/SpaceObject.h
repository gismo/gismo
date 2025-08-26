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
    class SpaceObject : public BaseObject<T, _order, false, _space>
    {
        using Base = BaseObject<T, _order, false, _space>;
        using Base::deriv_;
    public:
        typedef T Scalar; // Define Scalar type for this expression
        static constexpr size_t order = Base::order;
        static constexpr bool isConstant = Base::isConstant;
        static constexpr size_t space = Base::space;

    private:
        const gsFunctionSet<Scalar> * m_fs;
        const gsFuncData<Scalar>    * m_fd;
    public:
        SpaceObject(size_t domainDim, const std::array<size_t, order> & input_sizes)
        :
        BaseObject<Scalar, _order, false, _space>(domainDim, input_sizes),
        m_fs(NULL), m_fd(NULL)
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
            _print_impl<space>(os);
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

        template <size_t _SPACE>
        typename std::enable_if<_SPACE==Space::Test,void>::type
        _print_impl(std::ostream & os) const
        {
            os<<"φ";
            _print_arguments(os);
        }

        template <size_t _SPACE>
        typename std::enable_if<_SPACE==Space::Trial,void>::type
        _print_impl(std::ostream & os) const
        {
            os<<"ψ";
            _print_arguments(os);
        }

        template <size_t _SPACE>
        typename std::enable_if<_SPACE!=Space::Test && _SPACE!=Space::Trial,void>::type
        _print_impl(std::ostream & os) const
        {
            os<<"[Space object unknown]";
        }

    };
}//namespace Expr
}//namespace gismo