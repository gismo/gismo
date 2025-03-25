/** @file precomputed_expr.h

    @brief Defines expression precomputed_expr

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


namespace gismo
{

namespace expr
{

//template<class E>
template<class T>
class precomputed_expr : public _expr<precomputed_expr<T> >
{
public:
    typedef T Scalar;

    enum{Space = 1, ScalarValued=0, ColBlocks=0};

    //friend class gismo::gsExprHelper<Scalar>;
protected:

    const gsFeSpace<Scalar> m_rowvar, m_colvar;
    std::vector<gsMatrix<T> > m_data; // per element
    size_t m_curId;
    
public:

    // used by FeSpace, FeVariable, ..
    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(*this);
    }

    index_t cardinality_impl() const
    {
        return 0;
    }

private:
    void clear() { m_data.clear(); }

protected:

    template <typename E>
    explicit precomputed_expr(_expr<E> const& u) :
    m_rowvar(u.rowVar()), m_colvar(u.colVar())
    {  }

public:
    //bool isValid() const { return NULL!=m_fd && NULL!=m_fs; }

    gsMatrix<T> elementValues(size_t id)
    {
        return m_data[id];
    }
    
    MatExprType eval(const index_t k) const
    { return m_data[m_curId]; }

    const gsFeSpace<Scalar> & rowVar() const {return m_rowvar;}
    const gsFeSpace<Scalar> & colVar() const {return m_colvar;}

    index_t rows() const
    {
        return 0;
    }

    index_t cols() const { return 0; }

    void print(std::ostream &os) const { os << "u"; }
};


} //namespace expr
} //namespace gismo

