/** @file col_expr.h

    @brief Defines the column expression

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

/// Side tag of a symbol expression on a face or interface.
/// left/right evaluate one-sidedly; jump/avg are two-sided symbols whose
/// gsFuncData is the row-wise stack of the two sides (see gsExprHelper).
struct symbolSide
{
    enum mode { left = 0, right = 1, jump = 2, avg = 3 };
};

/**
 * @brief base expression
 * @todo Documentation
 * @ingroup Expressions
 * @tparam E The type of the expression
 */
template<class E>
class symbol_expr : public _expr<E>
{
public:
    typedef typename expr_traits<E>::Scalar Scalar;

    friend class gismo::gsExprHelper<Scalar>;
protected:
    const gsFunctionSet<Scalar> * m_fs; ///< Evaluation source for this FE variable
    const gsFuncData<Scalar>    * m_fd; ///< Temporary variable storing flags and evaluation data
    index_t m_d;                   ///< Dimension of this (scalar or vector) variable
    symbolSide::mode m_side; ///< which side(s) of a face this symbol is evaluated on

public:

    /// Returns whether this expression is evaluated across an interface
    bool isAcross() const { return symbolSide::right == m_side; }

    /// Returns whether this symbol carries data from BOTH sides of a face
    /// (jump/avg), stacked by gsExprHelper::precompute(const boundaryInterface&)
    bool isTwoSided() const
    { return symbolSide::jump == m_side || symbolSide::avg == m_side; }

    /// Returns the side tag (left/right/jump/avg) of this symbol
    symbolSide::mode sideMode() const { return m_side; }

    E right() const
    {
        E ac(this->derived());
        ac.m_fs = m_fs;//needed?
        ac.m_side = symbolSide::right;
        return ac;
    }

    E left() const
    {
        E ac(this->derived());
        ac.m_fs = m_fs;
        ac.m_side = symbolSide::left;
        return ac;
    }

    /// The jump [[u]] = u_left - u_right of this symbol across a face,
    /// evaluated by stacking the two sides' gsFuncData row-wise
    /// (gsExprHelper::precompute(const boundaryInterface&)). Only valid
    /// inside a face/interface loop (assembleSkeleton/assembleGhost/assembleIfc).
    E jump() const
    {
        E ac(this->derived());
        ac.m_fs = m_fs;
        ac.m_side = symbolSide::jump;
        return ac;
    }

    /// The average {u} = (u_left + u_right)/2 of this symbol across a face,
    /// evaluated by stacking the two sides' gsFuncData row-wise
    /// (gsExprHelper::precompute(const boundaryInterface&)). Only valid
    /// inside a face/interface loop (assembleSkeleton/assembleGhost/assembleIfc).
    E avg() const
    {
        E ac(this->derived());
        ac.m_fs = m_fs;
        ac.m_side = symbolSide::avg;
        return ac;
    }

    /// Returns the function source
    const gsFunctionSet<Scalar> & source() const {return *m_fs;}

    /// Returns the function data
    const gsFuncData<Scalar> & data() const
    {
        GISMO_ASSERT(NULL!=m_fd, "FuncData member not registered "<<this<<"/"<< m_fs);
        return *m_fd;
    }

    // used by FeSpace, FeVariable, ..
    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(*this);
        this->m_fd->flags |= NEED_VALUE | NEED_ACTIVE;
    }

    index_t cardinality_impl() const
    {
        GISMO_ASSERT(this->data().actives.rows()!=0,"Cardinality depends on the NEED_ACTIVE flag");
        return m_d * this->data().actives.rows();
    }

    //public for now due to Bc
    void setSource(const gsFunctionSet<Scalar> & fs) { m_fs = &fs;}

private:
    void setData(const gsFuncData<Scalar> & val) { m_fd = &val;}
    void setDim(index_t _d) { m_d = _d; }
    void clear() { m_fs = NULL; }

protected:
    explicit symbol_expr(index_t _d)
    : m_fs(NULL), m_fd(NULL), m_d(_d), m_side(symbolSide::left) { }

public:
    bool isValid() const { return NULL!=m_fd && NULL!=m_fs; }

    // component
    // expr comp(const index_t i) const { return comp_expr<Scalar>(*this,i); }
    // eval(k).col(i)

    // The evaluation return rows for (basis) functions and columns
    // for (coordinate) components
    MatExprType eval(const index_t k) const
    { return m_fd->values[0].col(k).blockDiag(m_d); } //!!
    //{ return m_fd->values[0].col(k); }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    index_t rows() const
    {
        GISMO_ASSERT(NULL!=m_fs, "FeVariable: Function source not registered");
        return m_fs->targetDim();
    }

    index_t cols() const { return m_d; }

    void print(std::ostream &os) const { os << "u"; }

public:

    /// Returns the vector dimension of the FE variable
    index_t dim() const { return m_d;}

    /// Returns the target dimension of the FE variable
    /// before vector-replication
    index_t targetDim() const { return m_fs->targetDim(); }

    /// Returns the parameter domain dimension the FE variable
    index_t parDim() const { return m_fs->domainDim(); }

    index_t cSize()  const
    {
        GISMO_ASSERT(0!=m_fd->values[0].size(),"Probable error.");
        return m_fd->values[0].rows();
    } // coordinate size
};

}// namespace expr
}// namespace gismo