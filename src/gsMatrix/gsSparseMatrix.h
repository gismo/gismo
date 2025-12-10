/** @file gsSparseMatrix.h

    @brief Provides declaration of the gsSparseMatrix class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

# pragma once

// Assumes that Eigen library has been already included

namespace gismo
{

/**
   @brief Class that provides a container for triplets (i,j,value) to be
   filled in a sparse matrix.

   Constructing a sparse matrix from triplets is much faster than
   inserting directly.  Use gsSparseMatrix().setFrom(gsSparseEntries)
   to pass the triplets to the matrix.

   \tparam T coefficient type
   \ingroup Matrix
*/
template<typename T>
class gsSparseEntries : public std::vector<gsEigen::Triplet<T,index_t> >
{
public:
    typedef gsEigen::Triplet<T,index_t> Triplet;
    typedef std::vector<gsEigen::Triplet<T,index_t> > Base;

    typedef typename Base::iterator       iterator;
public:

    inline void add( int i, int j, T value )
    { this->push_back( Triplet(i,j,value) ); }

    inline void addSorted(int i, int j, T value = 0)
    {
        Triplet t(i,j,value);
        iterator pos = std::lower_bound(Base::begin(), Base::end(), t, compTriplet );
        if ( pos == Base::end() || (pos->row() != t.row() && pos->col() !=t.col()) )// If not found
            Base::insert(pos, t);
    }

protected:
    // comparison operator for triplet struct, used to introduce columnwise lex-ordering.
    struct _compTriplet
    {
        bool operator() (const Triplet & left, const Triplet & right)
        {
            return (left.col() < right.col()) ||
                    (left.col() == right.col() && left.row() < right.row()) ;
        }
    } compTriplet ;

};

/**
    @brief Iterator over the non-zero entries of a sparse matrix

    This class is similar to gsEigen::SparseMatrix::InnerIteretor but in
    addition it is default-constructible and assignable.
*/
template<typename T, int _Options, typename _Index>
class gsSparseMatrixIter
{
    typedef gsEigen::SparseMatrix<T,_Options,_Index> SparseMatrix;
    static const int IsRowMajor = SparseMatrix::IsRowMajor;

public:
    gsSparseMatrixIter()
    : m_values(NULL), m_indices(NULL),
      m_end(NULL)   , m_outer(0)
    { }

    gsSparseMatrixIter(const SparseMatrix& mat, const _Index outer)
    : m_outer(outer)
    {
        const _Index oind = mat.outerIndexPtr()[outer];
        m_values  = mat.valuePtr()      + oind;
        m_indices = mat.innerIndexPtr() + oind;
        m_end = mat.isCompressed() ?
            mat.innerIndexPtr() + mat.outerIndexPtr()[outer+1] :
            mat.innerIndexPtr() + oind + mat.innerNonZeroPtr()[outer];
    }

    static gsSparseMatrixIter end(const SparseMatrix& mat, const _Index outer)
    {
        gsSparseMatrixIter o(mat,outer);
        o.m_values += o.m_end - o.m_indices;
        o.m_indices = o.m_end;
        return o;
    }

    inline gsSparseMatrixIter& operator++()
    { ++m_values; ++m_indices; return *this; }

    inline gsSparseMatrixIter& operator+=(size_t a)
    { m_values+=a; m_indices+=a; return *this; }

    inline gsSparseMatrixIter& operator+(size_t a)
    {
        gsSparseMatrixIter tmp(*this);
        return tmp+=a;
    }

    inline gsSparseMatrixIter& operator--()
    { --m_values; --m_indices; return *this; }

    inline gsSparseMatrixIter& operator-=(index_t a)
    { m_values-=a; m_indices-=a; return *this; }

    inline gsSparseMatrixIter& operator-(size_t a)
    {
        gsSparseMatrixIter tmp(*this);
        return tmp-=a;
    }

    inline const T& operator[](size_t i) const
    { return *(m_values+i); }

    inline T& operator[](size_t i)
    { return const_cast<T&>(*(m_values+i)); }

    inline const T& value() const { return *m_values; }
    inline T& valueRef() { return const_cast<T&>(*m_values); }

    inline _Index index() const { return *m_indices; }
    inline _Index outer() const { return m_outer; }
    inline _Index row() const { return IsRowMajor ? m_outer : index(); }
    inline _Index col() const { return IsRowMajor ? index() : m_outer; }

    inline operator bool() const { return (m_indices < m_end); }

protected:
    const T      * m_values;
    const _Index * m_indices;
    const _Index * m_end;
    _Index m_outer;
};

/** @brief Sparse matrix class, based on gsEigen::SparseMatrix.
 *
 * See http://eigen.tuxfamily.org/dox/group__SparseQuickRefPage.html
 * for Eigen's sparse matrix manipulations and
 * http://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html for
 * documentation of the gsEigen::SparseMatrix class.
 *
 * Remarks:
 *
 * An entry of the gsSparseMatrix can be accessed by <b>coeff( index
 * row, index col )</b> or just with the operator <b>( index row, index col )</b>.\n
 * An entry can be changed with either
 * <b>coeffRef( index row, index col)</b> or operator <b>( index row, index col )</b>.\n
 *
 *   \tparam T coefficient type
 *   \tparam _Option zero is ColMajor order.
 *   \tparam _Index index type
 *  \ingroup Matrix
 */

// Export the result to a file: saveAsBitmap(...);

template<typename T, int _Options, typename _Index>
class gsSparseMatrix : public gsEigen::SparseMatrix<T,_Options,_Index>
{
public:
    typedef gsEigen::SparseMatrix<T,_Options,_Index> Base;

    typedef gsSparseMatrixIter<T,_Options,_Index> iterator;

    // Type pointing to a block of the sparse matrix
    typedef typename gsEigen::Block<Base> Block;

    // Type pointing to a const block of the sparse matrix
    typedef typename gsEigen::Block<const Base> constBlock;

    // Type pointing to a block view of the sparse matrix
    typedef gsMatrixBlockView<Base> BlockView;

    // Type pointing to a block view of the sparse matrix
    typedef gsMatrixBlockView<const Base> constBlockView;

    /// Shared pointer for gsSparseMatrix
    typedef memory::shared_ptr<gsSparseMatrix> Ptr;

    /// Unique pointer for gsSparseMatrix
    typedef memory::unique_ptr<gsSparseMatrix> uPtr;

    /// Type of the full view of the matrix, for the case when only
    /// the lower diagonal part is stored
    typedef typename gsEigen::SparseSelfAdjointView<Base, Lower> fullView;

    /// Type of the full view of the matrix, for the case when only
    /// the lower diagonal part is stored
    typedef typename gsEigen::SparseSelfAdjointView<const Base, Lower> constFullView;

public:
    gsSparseMatrix() ;

    gsSparseMatrix(_Index rows, _Index cols) ;

    /// This constructor allows constructing a gsSparseMatrix from Eigen expressions
    template<typename OtherDerived>
    gsSparseMatrix(const gsEigen::EigenBase<OtherDerived>& other)  : Base(other) { }

    /// This constructor allows constructing a gsSparseMatrix from a selfadjoint view
    template<typename OtherDerived, unsigned int UpLo>
    gsSparseMatrix(const gsEigen::SparseSelfAdjointView<OtherDerived, UpLo>& other)
    : Base(other) { }

    /// This constructor allows constructing a gsSparseMatrix from Eigen expressions
    template<typename OtherDerived>
    gsSparseMatrix(const gsEigen::MatrixBase<OtherDerived>& other)  : Base(other) { }

    /// This constructor allows constructing a gsSparseMatrix from another sparse expression
    template<typename OtherDerived>
    gsSparseMatrix(const gsEigen::SparseMatrixBase<OtherDerived>& other)  : Base(other) { }

    /// This constructor allows constructing a gsSparseMatrix from Eigen expressions
    template<typename OtherDerived>
    gsSparseMatrix(const gsEigen::ReturnByValue<OtherDerived>& other)  : Base(other) { }

#if !EIGEN_HAS_RVALUE_REFERENCES
    // swap assignment operator
    gsSparseMatrix & operator=(gsSparseMatrix other)
    {
        this->swap(other);
        return *this;
    }

    template<typename OtherDerived, int a>
    gsSparseMatrix & operator=(const gsEigen::SparseSymmetricPermutationProduct<OtherDerived, a>& other)
    {
        this->Base::operator=(other);
        return *this;
    }
#else
#  ifdef _MSC_VER
    template <class EigenExpr>
    gsSparseMatrix& operator= (const EigenExpr & other)
    {
        this->Base::operator=(other);
        return *this;
    }
#  else
    using Base::operator=;
#  endif

    // Avoid default keyword for MSVC<2013
    // https://msdn.microsoft.com/en-us/library/hh567368.aspx
    gsSparseMatrix(const gsSparseMatrix& other) : Base(other)
    { Base::operator=(other); }
    gsSparseMatrix& operator= (const gsSparseMatrix & other)
    { Base::operator=(other); return *this; }

    gsSparseMatrix(gsSparseMatrix&& other)
    { gsSparseMatrix::operator=(std::forward<gsSparseMatrix>(other)); }

    gsSparseMatrix & operator=(gsSparseMatrix&& other)
    {
        this->swap(other);
        other.clear();
        return *this;
    }

#endif

    /**
       \brief This function returns a smart pointer to the
       matrix. After calling it, the matrix object becomes empty, ie
       the size of the matrix is 0
    */
    uPtr moveToPtr()
    {
        uPtr m(new gsSparseMatrix);
        m->swap(*this);
        return m;
    }

    /// \brief Returns an iterator to the first non-zero elemenent of
    /// column \ a outer (or row \a outer if the matrix is RowMajor)
    inline iterator begin(const index_t outer) const { return iterator(*this,outer);}

    /// \brief Returns an iterator to the  end position of
    /// column \ a outer (or row \a outer if the matrix is RowMajor)
    inline iterator end(const index_t outer) const { return iterator::end(*this,outer);}

    void clear()
    {
        this->resize(0,0);
        this->data().squeeze();
    }

    std::pair<index_t,index_t> dim() const
    { return std::make_pair(this->rows(), this->cols() ); }

    void reservePerColumn(const index_t nz)
    {
        Base::reserve(gsVector<index_t>::Constant(this->outerSize(), nz));
    }

    void setFrom( gsSparseEntries<T> const & entries) ;

    inline T   at (_Index i, _Index j ) const { return this->coeff(i,j); }
    inline T & at (_Index i, _Index j ) { return this->coeffRef(i,j); }

    inline T    operator () (_Index i, _Index j ) const { return this->coeff(i,j); }
    inline T  & operator () (_Index i, _Index j ) { return this->coeffRef(i,j); }

    /// Add to entry (\a i, \a j) the value \a val, but not an
    /// explicit zero
    inline void addTo(_Index i, _Index j, const T val)
    {
        if (0!=val) this->coeffRef(i,j) += val;
    }

    /// Insert to entry (\a i, \a j) the value \a val, but not an
    /// explicit zero
    inline void insertTo(_Index i, _Index j, const T val)
    {
        if (0!=val)
            this->insert(i,j) = val;
    }

    inline bool isExplicitZero(_Index row, _Index col) const
    {
        const _Index outer = gsSparseMatrix::IsRowMajor ? row : col;
        const _Index inner = gsSparseMatrix::IsRowMajor ? col : row;
        const _Index end = this->m_innerNonZeros ?
            this->m_outerIndex[outer] + this->m_innerNonZeros[outer] : this->m_outerIndex[outer+1];
        const _Index id = this->m_data.searchLowerIndex(this->m_outerIndex[outer], end-1, inner);
        return !((id<end) && ( this->m_data.index(id)==inner));
    }

    /// Adds an explicit zero, only if (row,col) is not in the matrix
    inline void insertExplicitZero(_Index row, _Index col)
    {
        GISMO_ASSERT(row>=0 && row<this->rows() && col>=0 && col<this->cols(), "Invalid row/col index.");
        const _Index outer = gsSparseMatrix::IsRowMajor ? row : col;
        const _Index inner = gsSparseMatrix::IsRowMajor ? col : row;
        const _Index start = this->m_outerIndex[outer];
        const _Index end   = this->m_innerNonZeros ?
            this->m_outerIndex[outer] + this->m_innerNonZeros[outer] : this->m_outerIndex[outer+1];
        if(end<=start)
        { this->insert(row,col); return; }
        const _Index p = this->m_data.searchLowerIndex(start,end-1,inner);
        if(! ((p<end) && (this->m_data.index(p)==inner)) )
            this->insert(row,col);
    }

    inline T& coeffUpdate(_Index row, _Index col)
    {
        GISMO_ASSERT(row>=0 && row<this->rows() && col>=0 && col<this->cols(), "Invalid row/col index.");
        const _Index outer = gsSparseMatrix::IsRowMajor ? row : col;
        _Index start = this->m_outerIndex[outer];
        _Index end = this->m_innerNonZeros ? this->m_outerIndex[outer] + this->m_innerNonZeros[outer]
            : this->m_outerIndex[outer+1];
        const _Index inner = gsSparseMatrix::IsRowMajor ? col : row;
        const _Index p = this->m_data.searchLowerIndex(start,end-1,inner);
        if((p<end) && (this->m_data.index(p)==inner))
            return this->m_data.value(p);
        GISMO_ERROR("(row,col) = ("<< row <<","<<col<<") is not in the matrix (sparsity pattern).");
    }

    /// Return a block view of the matrix with \a rowSizes and \a colSizes
    BlockView blockView(const gsVector<index_t> & rowSizes,
                        const gsVector<index_t> & colSizes)
    {
        return BlockView(*this, rowSizes, colSizes);
    }

    /// Return a const block view of the matrix with \a rowSizes and \a colSizes
    constBlockView blockView(const gsVector<index_t> & rowSizes,
                        const gsVector<index_t> & colSizes) const
    {
        return constBlockView(*this, rowSizes, colSizes);
    }

    /// Returns a pointer wrapped as a gsAsConstVector, which contains
    /// the number of non-zero entries per column. Note that the
    /// matrix must be uncompressed format for this to work
    gsAsConstVector<_Index> nonZerosPerCol()
    {
        if ( this->isCompressed() )
        {
            gsWarn<<"nonZerosPerCol(): Uncompressing the gsSparseMatrix.\n";
            this->uncompress();
        }
        return gsAsConstVector<_Index>(this->innerNonZeroPtr(),this->innerSize());
    }

    std::string printSparsity() const
    {
        std::ostringstream os;
        os <<"Sparsity: "<< std::fixed << std::setprecision(2)
           <<(double)100*this->nonZeros()/this->size() <<'%'<<", nnz: "<<this->nonZeros() <<"\n";
        for (index_t i = 0; i!=this->rows(); ++i)
        {
            for (index_t j = 0; j!=this->cols(); ++j)
                os << ( isExplicitZero(i,j) ? "\u00B7" :
                        ( 0 == this->coeff(i,j) ? "o" : "x") );
            os<<"\n";
        }
        return os.str();
    }

    /// Allocates the sparsity pattern
    //template<typename InputIterator>
    //void setPattern(const InputIterator& begin, const InputIterator& end)

    void rrefInPlace();

    /// Returns a set of (inner) indices, which consists of all
    /// rows/columns of the matrix with non-zero coefficients at
    /// columns/rows \a outer
    template<class container>
    container innerOf(const container & outer) const
    {
        std::vector<bool> v(this->innerSize(), false);
        gsSparseMatrix<>::iterator mIt;
        for (typename container::const_iterator k =
                 outer.begin(); k!=outer.end(); ++k )
            for( mIt = this->begin(*k); mIt; ++mIt)
                v[mIt.index()] = (0!=mIt.value());

        container inner;
        inner.reserve( std::count(v.begin(), v.end(), true) );
        std::vector<bool>::iterator it = std::find(v.begin(),v.end(), true);
        while (v.end() != it)
        {
            inner.push_back(std::distance(v.begin(),it));
            it = std::find(++it,v.end(), true);
        }
        return inner;
    }

    /// Returns the result of multiplication of \a this and \a other,
    /// where other has rows/columns indexed by \a inner and \a outer
    template<class container>
    gsMatrix<T> multiplyBy(const container & inner,
                           const container & outer,
                           const gsMatrix<T> & other)
    {
        // todo: better implementation
        return submatrix(inner,outer) * other;
    }

    /// Returns the submatrix consisting of the rows and columns
    /// indexed by the vector containers \a inner and \a outer
    template<class container>
    gsMatrix<T> submatrix(const container & inner,
                          const container & outer) const
    {
        const index_t nr = inner.size();
        const index_t nc = outer.size();
        gsMatrix<T> result(nr, nc);
        result.setZero();
        for (index_t c=0; c!=nc; ++c)
        {
            iterator it = begin(outer[c]); //sorted iteration
            for (index_t r=0; r!=nr && it;)
            {

                while(it && inner[r] > it.index()) ++it;
                if (it) while(r!=nr && inner[r] < it.index()) ++r;
                if (r!=nr)
                {
                    result(r,c) = it.value();
                    ++r; ++it;
                }
            }
        }
        return result;
    }

    gsVector<index_t> nonZerosPerInner(index_t upto = std::numeric_limits<index_t>::max())  const
    {
        upto = math::min(upto, this->outerSize());
        gsVector<index_t> nz(upto);
        index_t * v = nz.data();
        for (index_t i = 0; i != upto; ++i, ++v)
            *v = this->innerVector(i).nonZeros();
        return nz;
    }

    /// Returns the Kronecker product of \a this with \a other
    gsSparseMatrix kron(const gsSparseMatrix& other)  const
    {
        const index_t r  = this->rows(), c = this->cols();
        const index_t ro = other.rows(), co = other.cols();
        gsSparseMatrix result(r*ro, c*co);
        if (0 == result.size()) return result;
        result.reserve(this->nonZerosPerInner()
                       .kron(other.nonZerosPerInner()));

        iterator it1, it2;
        for (index_t k1=0; k1 != (gsSparseMatrix::IsRowMajor?r:c); ++k1)
            for (it1 = this->begin(k1); it1; ++it1)
                for (index_t k2=0; k2 != (gsSparseMatrix::IsRowMajor?ro:co); ++k2)
                    for (it2 = other.begin(k2); it2; ++it2)
                    {
                        const index_t i = it1.row() * ro + it2.row();
                        const index_t j = it1.col() * co + it2.col();
                        result.insert(i,j) = it1.value() * it2.value();
                    }
        result.makeCompressed();
        return result;
    }

private:

    /*
      The inherited setZero() destroys column-nonzeros structure and
      can make further computations very slow.  Almost equivalent to
      M.setZero() is:
      M = gsSparseMatrix(M.rows(), M.cols());

      void setZero();
    */

}; // class gsSparseMatrix


//template<class T>
//inline void gsSparseEntries<T>::add( int const& i, int const& j, const T & value)
//        { this->push_back( Triplet(i,j,value) ); }


template<typename T, int _Options, typename _Index> inline
gsSparseMatrix<T, _Options, _Index>::gsSparseMatrix() : Base() { }

template<typename T, int _Options, typename _Index> inline
gsSparseMatrix<T, _Options, _Index>::gsSparseMatrix(_Index rows, _Index cols) : Base(rows,cols) { }

template<typename T, int _Options, typename _Index> inline
void gsSparseMatrix<T, _Options, _Index>::setFrom( gsSparseEntries<T> const & entries)
{ this->setFromTriplets(entries.begin(),entries.end() ); }

template<typename T, int _Options, typename _Index> void
gsSparseMatrix<T, _Options, _Index>::rrefInPlace()
{
    gsMatrix<T,1,Dynamic> R;
    index_t c_i, c_j;
    const index_t nr = this->rows();
    const index_t nc = this->cols();

    gsVector<index_t> piv(nr);

    for (index_t i = 0; i!=nr; ++i)
    {
        // pull row(i)
        R = this->row(i);
        for (index_t j = 0; j!=i; ++j)
        {
            c_j = piv[j];
            if (c_j != nc )
                R.tail(nc-c_j) -= R.at(c_j) * this->row(j).tail(nc-c_j);
        }

            // store row pivot
            c_i = 0; while (c_i!= nc && 0 == R.at(c_i)) ++c_i;
            piv[i] = c_i;
            if (c_i == nc ) continue;
            R.tail(nc-c_i-1).array() /= R.at(c_i);
            R.at(c_i) = 1;
            // push row(i)
            this->row(i) = R.sparseView();
        }

        // sort rows wrt pivots
        bool didSwap;
        c_i = nr - 1; // last swap done
        do
        {
            didSwap = false;
            c_j = c_i;

            for( index_t i=0; i != c_j; i++)
                if( piv[i] > piv[i+1] )
                {
                    std::swap(piv[i], piv[i+1]);
                    //this->row(i).swap(this->row(i+1));
                    didSwap = true;
                    c_i = i;
                }
        } while(didSwap);

        // todo: precompute nz (piv[nz])

        // lower part
        for (index_t i = 0; i!=nr; ++i)
        {
            // pull row(i)
            R = this->row(i);
            c_i = piv[i];
            if (c_i == nc ) break;

            for (index_t j = i+1; j!=nr; ++j) // nr - zr
            {
                c_j = piv[j];
                if ( c_j == nc ) break;
                R.tail(nc-c_j) -= R.at(c_j) * this->row(j).tail(nc-c_j);
            }
            // push row(i)
            this->row(i) = R.sparseView();
        }
}

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsSparseMatrix
   */
  template<typename T>
  void pybind11_init_gsSparseMatrix(pybind11::module &m, const std::string & typestr)
  {
    using Class = gsSparseMatrix<T>;
    std::string pyclass_name = std::string("gsSparseMatrix") + typestr;
    pybind11::class_<Class>(m, pyclass_name.c_str(), pybind11::buffer_protocol(), pybind11::dynamic_attr())
    // Constructors
    .def(pybind11::init<>())
    .def(pybind11::init<index_t, index_t>())
    // Member functions
    .def("size",      &Class::size)
    .def("rows",      &Class::rows)
    .def("cols",      &Class::cols)
    // .def("transpose", &Class::transpose)
    // .def("addTo",     &Class::addTo)
    // .def("insertTo",  &Class::insertTo)
    .def("toDense",   &Class::toDense)
    ;
  }

#endif // GISMO_WITH_PYBIND11

} // namespace gismo


namespace gsEigen { namespace internal {
template<typename T, int _Options, typename _Index>
struct traits<gismo::gsSparseMatrix<T,_Options,_Index> >:
gsEigen::internal::traits<gsEigen::SparseMatrix<T,_Options,_Index> > { };
} }

/* *****************************************************************
#ifdef GISMO_BUILD_LIB
#ifdef gsMatrix_EXPORT
#undef  EXTERN_CLASS_TEMPLATE
#define EXTERN_CLASS_TEMPLATE CLASS_TEMPLATE_INST
#endif
namespace gismo
{
EXTERN_CLASS_TEMPLATE gsSparseMatrix<real_t>;
}
#endif
// *****************************************************************
*/
