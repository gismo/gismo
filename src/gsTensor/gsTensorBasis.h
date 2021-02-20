/** @file gsTensorBasis.h

    @brief Provides declaration of TensorBasis class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsCore/gsBoundary.h>

namespace gismo
{

/** 
 *  @brief Abstract base class for tensor product bases.
 *
 *   \param d dimension of the parameter domain
 *   \param Basis_t type of the coordinate-wise bases
 *
 *   \ingroup basis
 *   \ingroup Tensor
 */

template<short_t d, class T>
class gsTensorBasis : public gsBasis<T>
{
public: 
    typedef memory::shared_ptr< gsTensorBasis > Ptr;
    typedef memory::unique_ptr< gsTensorBasis > uPtr;

    typedef gsTensorBasis<d,T> Self_t;

    typedef gsBasis<T> Basis_t;

    typedef Basis_t CoordinateBasis;

    /// Coefficient type
    typedef T Scalar_t;

    /// Dimension of the parameter domain
    static const short_t Dim = d;

    /// Iterators on coordinate bases
    typedef Basis_t** iterator;
    typedef Basis_t* const* const_iterator;

public:
    
    gsTensorBasis() : m_bases() { }

    ~gsTensorBasis() { freeAll(m_bases, m_bases+d); }

    gsTensorBasis( const gsTensorBasis & o);
    gsTensorBasis& operator=( const gsTensorBasis & o);
    
#if EIGEN_HAS_RVALUE_REFERENCES
    gsTensorBasis(gsTensorBasis&& other)
    { 
        util::copy(other.m_bases, other.m_bases+d, m_bases);
        std::fill (other.m_bases, other.m_bases+d, nullptr);
    }
    gsTensorBasis & operator=(gsTensorBasis&&other)
    {
        freeAll(m_bases, m_bases+d);
        util::copy(other.m_bases, other.m_bases+d, m_bases);
        std::fill (other.m_bases, other.m_bases+d, nullptr);
        return *this;
    }
#endif
    bool isValid() const { return std::find(m_bases,m_bases+d,
                                            static_cast<Basis_t*>(0)) == m_bases+d; }
    
    /// Constructor 2D (takes ownership of the passed bases)
    gsTensorBasis( Basis_t* x,  Basis_t* y);
    // template<class U> gsTensorBasis(   
    // typename enable_if<d==2 && is_same<U,Basis_t*>::value,U>::type x, U y) 
    // { m_bases[0] = x; m_bases[1] = y; }

    /// Constructor 3D (takes ownership of the passed bases)
    gsTensorBasis( Basis_t* x,  Basis_t* y, Basis_t* z ) ;
    
    /// Constructor 4D (takes ownership of the passed bases)
    gsTensorBasis( Basis_t* x,  Basis_t* y, Basis_t* z, Basis_t* w ) ;
    
    /// Constructor nD (takes ownership of the passed bases)
    explicit gsTensorBasis(iterator it) 
    {
        for (short_t i = 0; i < d; ++i)
            m_bases[i] = *(it++);
    }
    
public:

    // Returns the dimension of the basis
    short_t domainDim() const { return Dim; }

    /// Returns the number of elements in the basis
    index_t size() const 
    {
        index_t r=1;
        for (short_t i = 0; i < d; ++i)
            r *= m_bases[i]->size();
        return r; 
    }

    // Look at gsBasis class for a description
    size_t numElements() const
    {
        size_t nElem = m_bases[0]->numElements();
        for (short_t dim = 1; dim < d; ++dim)
            nElem *= m_bases[dim]->numElements();
        return nElem;
    }

    // Look at gsBasis class for a description
    size_t numElements(boxSide const & s) const
    {
        const short_t dir =  s.direction();
        size_t nElem = 1;
        for (short_t dim = 0; dim < d; ++dim)
        {
            if(dim == dir)
                continue;
            nElem *= m_bases[dim]->numElements();
        }
        return nElem;
    }

    // Look at gsBasis class for a description
    size_t elementIndex(const gsVector<T> & u ) const
    {
        GISMO_ASSERT( u.rows() == d, "Wrong vector dimension");

        size_t ElIndex = m_bases[d-1]->elementIndex( u.col(d-1) );
        for ( short_t i=d-2; i>=0; --i )
            ElIndex = ElIndex * m_bases[i]->numElements() 
                    + m_bases[i]->elementIndex( u.col(i) );

        return ElIndex;        
    }

    // Look at gsBasis class for a description
    gsMatrix<T> elementInSupportOf(index_t j) const;
    
    /// Returns the number of elements (component wise)
    void numElements_cwise(gsVector<unsigned>& result) const
    {
        result.resize(d);
        for (short_t dim = 0; dim < d; ++dim)
            result(dim) = static_cast<unsigned>(m_bases[dim]->numElements());
    }

    /// Returns the anchors (graville absissae) that represent the members of the basis
    void anchors_into(gsMatrix<T>& result) const;

    /// Returns the anchors (graville absissae) that represent the members of the basis
    void anchor_into(index_t i, gsMatrix<T>& result) const;

    // TODO: Why is this documentation not in gsBasis?
    /**
     * \brief Returns the indices of active basis functions at points
     * <em>u</em>, as a list of indices, in <em>result</em>. A
     * function is said to be <em>active</em> in a point if this point
     * lies in the closure of the function's support.
     *
     * \par Tensor indexing in result
     * Assume that the parameter domain is three dimensional.
     * Let <em>n1</em>, <em>n2</em>, and <em>n3</em> denote the number of \em univariate basis
     * functions in the first, second and third coordinate direction, respectively.\n
     * Let the <em>trivariate</em> tensor product basis function <em>B_I</em> be defined by\n
     * <em>B_I(x,y,z) = B_i(x) * B_j(y) * B_k(z)</em>.\n Then, the index \em I, which is
     * returned in \em result, is computed as \n
     * <em>I = i + j * n1 + k * n1*n2</em>.\n
     * Examples:\n
     * I <-> (i,j,k)\n
     * 0 <-> (0,0,0)\n
     * 1 <-> (1,0,0)\n
     * 2 <-> (2,0,0)\n
     * ...\n
     * (n1-1) <-> (n1-1,0,0)\n
     * n1 <-> (0,1,0)\n
     * n1+1 <-> (1,1,0)\n
     * n1+2 <-> (2,1,0)\n
     * ...\n
     * n1*n2-1 <-> (n1,n2,0) \n
     * n1*n2 <-> (0,0,1) \n
     * n1*n2+1 <-> (1,0,1) \n
     * ...\n
     * n1*n2*n3-1 <-> (n1,n2,n3) \n
     *
     * \param[in] u  gsMatrix containing evaluation points. Each column represents one evaluation point.
     * \param[out]  result For every column \a i of \a u, a column containing the
     *   active basis functions at evaluation point <em>u</em>.col(<em>i</em>)
     *
     */
    virtual void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const;

    // Look at gsBasis class for documentation 
    bool isActive(const index_t i, const gsVector<T>& u) const;

    /// Returns a box with the coordinate-wise active functions
    /// \param u evaluation points
    /// \param low lower left corner of the box
    /// \param upp upper right corner of the box
    void active_cwise(const gsMatrix<T> & u, gsVector<index_t,d>& low,
                      gsVector<index_t,d>& upp ) const;

    // Look at gsBasis class for documentation 
    virtual void connectivity(const gsMatrix<T> & nodes, gsMesh<T> & mesh) const;

    /// Returns the indices of the basis functions that touch the domain
    /// boundary
    gsMatrix<index_t> allBoundary( ) const ;

    /// Returns the indices of the basis functions that touch the domain
    /// boundary
    gsMatrix<index_t> boundaryOffset(boxSide const & s, index_t offset) const;

    index_t functionAtCorner(boxCorner const & c) const;

    /// Returns the components for a basis on the face \a s 
    void getComponentsForSide(boxSide const & s, std::vector<Basis_t*> & rr) const;

    // see gsBasis for doxygen documentation
    // Returns a bounding box for the basis' domain
    gsMatrix<T> support() const ;

    // see gsBasis for doxygen documentation
    // Returns a bounding box for the support of the ith basis function
    gsMatrix<T> support(const index_t & i ) const ;

    // see gsBasis for doxygen documentation
    // Evaluates the non-zero basis functions (and optionally their
    // first k derivatives) at value u into result
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const ;

    // see gsBasis for doxygen documentation
    // Evaluate the i-th basis function at all columns of the matrix
    // (or vector) u
    void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const ;

    /// Evaluate an element of the space given by coefs at points u
    virtual void eval_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const;

    // see gsBasis for doxygen documentation
    // Evaluate the nonzero basis functions and their derivatives up to
    // order n at all columns of u
    virtual void evalAllDers_into(const gsMatrix<T> & u, int n,
                                  std::vector<gsMatrix<T> >& result) const;

    // see gsBasis for doxygen documentation
    // Evaluates the gradient the non-zero basis functions at value u.
    virtual void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    // Evaluates the second derivatives of the non-zero basis functions at value u.
    virtual void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

private:
    // Internal function
    //
    // values: array of std::vectors of gsMatrix<T>
    // values[i], i = 0,...,d-1, contains the result of
    // the univariate evalAllDers_into, corresponding to coordinate direction i.
    //
    // values[i] is a std::vector< gsMatrix<T> >
    // values[i][j] contains the j-th derivatives in coordinate direction i
    //
    // size: gsVector of length d, size[i] contains the number of basis functions
    // in coordinate direction i.
    static void deriv2_tp(const std::vector< gsMatrix<T> > values[],
                   const gsVector<unsigned, d> & size,
                   gsMatrix<T>& result);

public:
    // see gsBasis for doxygen documentation
    // Evaluate the i-th basis function derivative at all columns of
    virtual void derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const ;

    virtual void deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const ;

    // Evaluates the (partial) derivatives of an element given by coefs at (the columns of) u.
    //void deriv_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for documentation 
    typename gsBasis<T>::domainIter makeDomainIterator() const
    {
        return typename gsBasis<T>::domainIter(new gsTensorDomainIterator<T,d>(*this) );
    }

    // Look at gsBasis class for documentation 
    typename gsBasis<T>::domainIter makeDomainIterator(const boxSide & s) const
    {
        return ( s == boundary::none ? 
                 typename gsBasis<T>::domainIter(new gsTensorDomainIterator<T,d>(*this)) :
                 typename gsBasis<T>::domainIter(new gsTensorDomainBoundaryIterator<T,d>(*this, s))
                );
    }

    // Look at gsBasis class for documentation 
    virtual typename gsGeometry<T>::uPtr interpolateAtAnchors(gsMatrix<T> const& vals) const;

    /// Interpolates values on a tensor-grid of points, given in
    /// tensor form (d coordinate-wise vectors). Samples \a vals
    /// should be ordered as the tensor-basis coefficients
    typename gsGeometry<T>::uPtr interpolateGrid(gsMatrix<T> const& vals,
                                    std::vector<gsMatrix<T> >const& grid) const;

    /// Prints the object as a string, pure virtual function of gsTensorBasis.
    virtual std::ostream &print(std::ostream &os) const = 0;

    // Look at gsBasis class for documentation 
    virtual void uniformRefine(int numKnots = 1, int mul=1)
    {
        for (short_t j = 0; j < d; ++j)
            m_bases[j]->uniformRefine(numKnots,mul);
    }

    /** \brief Refine elements defined by their tensor-index.
     *
     * In a tensor mesh, each element has a unique index computed
     * as follows:
     *
     * Let \f$n_i\f$ denote the number of basis functions in the
     * <em>i</em>-th component, and let \f$k_i\f$ denote the index
     * of an element in the (1-dimensional) mesh of the <em>i</em>-th
     * component. The global index of element \f$(a,b,c)\f$ is
     * given by \f$a + b \cdot n_1 + c\cdot n_1 \cdot n_2\f$.
     *
     * \param[in] elements vector of unsigned containing the
     * indices of the elements that should be refined (see above).
     */
    void refineElements(std::vector<index_t> const & elements);

    /// Refine the basis uniformly and perform knot refinement for the
    /// given coefficient vector
    void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots=1, int mul=1);

    /// Refine the basis uniformly and produce a sparse matrix which
    /// maps coarse coefficient vectors to refined ones
    void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots=1, int mul=1);
    
    // Look at gsBasis class for documentation 
    virtual void uniformCoarsen(int numKnots = 1)
    {
        for (short_t j = 0; j < d; ++j)
            m_bases[j]->uniformCoarsen(numKnots);
    }

    // Look at gsBasis class for documentation
    void uniformCoarsen_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots=1);

    // Look at gsBasis class for documentation 
    virtual void degreeElevate(short_t const & i = 1, short_t const dir = -1)
    { 
        if (dir == -1)
        {
            for (short_t j = 0; j < d; ++j)
                m_bases[j]->degreeElevate(i);
        }
        else 
        {
            GISMO_ASSERT( dir < this->dim(),
                          "Invalid basis component requested" );
            m_bases[dir]->degreeElevate(i);
        }
    }

    // Look at gsBasis class for documentation
    virtual void degreeIncrease(short_t const & i = 1, short_t const dir = -1)
    {
        if (dir == -1)
        {
            for (short_t j = 0; j < d; ++j)
                m_bases[j]->degreeIncrease(i);
        }
        else
        {
            GISMO_ASSERT( dir < this->dim(),
                          "Invalid basis component requested" );
            m_bases[dir]->degreeIncrease(i);
        }
    }

    // Look at gsBasis class for documentation
    virtual void degreeDecrease(short_t const & i = 1, short_t const dir = -1)
    {
        if (dir == -1)
        {
            for (short_t j = 0; j < d; ++j)
                m_bases[j]->degreeDecrease(i);
        }
        else
        {
            GISMO_ASSERT( dir < this->dim(),
                          "Invalid basis component requested" );
            m_bases[dir]->degreeDecrease(i);
        }
    }

    // Look at gsBasis class for documentation 
    virtual void degreeReduce(short_t const & i = 1, short_t const dir = -1)
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT( dir < this->dim(),
                      "Invalid basis component requested" );
        for (short_t j = 0; j < d; ++j)
            m_bases[j]->degreeReduce(i);
    }

    /// Get a const-iterator to the beginning of the bases vector
    /// \return an iterator to the beginning of the bases vector
    const_iterator begin() const
    { return &m_bases[0]; }
  
    /// Get a const-iterator to the end of the  bases vector
    /// \return an iterator to the end of the  bases vector
    const_iterator end() const
    { return &m_bases[d]; }
  
    /// Get an iterator to the beginning of the  bases vector
    /// \return an iterator to the beginning of the  bases vector
    iterator begin()
    { return &m_bases[0]; }

    /// Get an iterator to the end of the  bases vector
    /// \return an iterator to the end of the  bases vector
    iterator end()
    { return &m_bases[d]; }

    /// The number of basis functions in the direction of the k-th parameter component
    index_t size(short_t k) const { return m_bases[k]->size(); }

    /// The number of basis functions in the direction of the k-th parameter component
    template<int s>
    void size_cwise(gsVector<index_t,s> & result) const
    {
        result.resize(d);
        for ( short_t k = 0; k!=d; ++k )
            result[k] = m_bases[k]->size(); 
    }
    
    /**\brief 
       Returns all the basis functions with tensor-numbering \param k in direction \param dir
    
    ## Detailed explanation: ##
    Tensor-numbering in N-variate tensor-product basis means that each basis function
    is assigned an identifier (i_0, i_1, ..., i_{N-1}).
    This function returns indices of basis functions with i_dir = k
    and the returned indices are numbering of the basis functions in the basis
    (i.e., 0,1, ..., basis.size() ).
    ## Example: ##
    Bivariate tensor-product basis functions have tensor numbering (a,b).
    Calling dir=0, k=1 gives all functions with tensor-numbering (1,b).
    Calling dir=1, k=3 gives all functions with tensor-numbering (a,3).
    */
    gsMatrix<index_t> coefSlice(short_t dir, index_t k) const;

    /// Returns the degree of the basis wrt variable \a i 
    short_t degree(short_t i) const
    { 
        return m_bases[i]->degree(0); 
    }

    short_t maxDegree() const
    { 
        short_t td = m_bases[0]->degree(0);
        // take maximum of coordinate bases degrees
        for (short_t k=1; k!=d; ++k)
            td = math::max(td, m_bases[k]->degree(0));
        return td;
    }
    
    short_t minDegree() const
    { 
        short_t td = m_bases[0]->degree(0);
        // take minimum of coordinate bases degrees
        for (short_t k=1; k!=d; ++k)
            td = math::min(td, m_bases[k]->degree(0));
        return td;
    }
    
    short_t totalDegree() const
    { 
        short_t td = 0;
        for (short_t k=0; k!=d; ++k)
            td = + m_bases[k]->degree(0);
        return td;
    }

    gsVector<int> cwiseDegree() const
    {
        gsVector<int> deg(d);
        for (short_t k=0; k!=d; ++k)
            deg[k] = m_bases[k]->degree(0);
        return deg;
    }

    /// Returns the global index of the basis function created by
    /// components of indices i,j,k (for 2d or 3d only)
    inline unsigned index(unsigned i, unsigned j, unsigned k=0) const;

    /// Returns the stride for dimension dir
    inline unsigned stride(short_t dir) const;

    /// Returns the strides for all dimensions
    void stride_cwise(gsVector<index_t,d> & result) const 
    { 
        //result.resize(d);
        result[0] = 1;
        for ( short_t i=1; i != d; ++i )
            result[i] = result[i-1] * m_bases[i-1]->size();
    }

    /// Returns the global index of the basis function created by
    /// components of indices given in the vector v
    inline index_t index(gsVector<index_t,d> const & v) const;
    //  inline unsigned index(gsVector<unsigned>         & v) const;

    /// \brief Returns the tensor index of the basis function with
    /// global index \a m.
    inline gsVector<index_t, d> tensorIndex(const index_t& m) const
    {
        gsVector<index_t, d> ind;
        int mm = m;
        for (short_t i = 0; i<d; ++i )
        {
            ind(i)= mm % size(i);
            mm -= ind(i);
            mm /= size(i);
        }
        return ind;
    }
    
    void swapDirections(const unsigned i, const unsigned j)
    {
        GISMO_ASSERT( static_cast<int>(i) < Dim && static_cast<int>(j) < Dim,
                      "Invalid basis components "<<i<<" and "<<j<<" requested" );
        std::swap(m_bases[i],m_bases[j]);
    }

    /// \brief Returns true iff the basis function with multi-index
    /// \em ind is on the boundary
    inline bool indexOnBoundary(const gsVector<index_t, d> & ind) const 
    {
        for ( short_t i = 0; i < d; ++i )
            if ( ind[i] == size(i)-1 )
                return true;
        return ((0 == ind.array()).any() );
    }

    /// \brief Returns true iff the basis function indexed \a m is on
    /// the boundary
    inline bool indexOnBoundary(const index_t m) const 
    {
        return ( indexOnBoundary( tensorIndex(m) ) );
    }

    // see gsBasis for documentation
    void matchWith(const boundaryInterface & bi, const gsBasis<T> & other,
                   gsMatrix<index_t> & bndThis, gsMatrix<index_t> & bndOther) const;

    /// Get the minimum mesh size, as expected for inverse inequalities
    virtual T getMinCellLength() const;
    
    /// Get the maximum mesh size, as expected for approximation error estimates
    virtual T getMaxCellLength() const;
    
    Basis_t& x() const 
    { 
        return *m_bases[0]; 
    }

    Basis_t& y() const { 
        if (d > 1) return *m_bases[1]; 
        else 
            GISMO_ERROR("gsTensorBasis has no y component"); 
    }
  
    Basis_t& z() const { 
        if (d > 2) 
            return *m_bases[2]; 
        else
            GISMO_ERROR("gsTensorBasis has no z component"); 
    }

    Basis_t& component(short_t dir)
    { 
        GISMO_ASSERT( dir < Dim,
                      "Invalid basis component requested" );
        return *m_bases[dir];
    }
    
    const Basis_t & component(short_t dir) const
    { 
        GISMO_ASSERT( dir < Dim,
                      "Invalid basis component requested" );
        return *m_bases[dir];
    }

    //inline int trueSize(int k) const { return m_bases[k]->trueSize(); }

// Data members
protected:

    void swap(gsTensorBasis& o)
    { std::swap_ranges(m_bases, m_bases+d, o.m_bases); }
    
    Basis_t* m_bases[d];

}; // class gsTensorBasis

// Next line disallows instantization of gsTensorBasis<0,T>
template<typename T> class gsTensorBasis<0,T>
{using T::GISMO_ERROR_gsTensorBasis_cannot_have_dimension_zero;};

/* 
 *  @brief 
 *  Class for a Tensor product spline space of dimension 1.
 *  This specialization is mainly for compatibility.
 *
 *   \tparam T coefficient type
 *
 *  \ingroup Tensor
 */
template<class T>
class gsTensorBasis<1,T> : public gsBasis<T>
{
public: 
    typedef memory::shared_ptr< gsTensorBasis<1,T> > Ptr;
    typedef memory::unique_ptr< gsTensorBasis<1,T> > uPtr;

    static const short_t Dim = 1;

    typedef gsBasis<T> Base;

    typedef gsTensorBasis<1,T> Self_t;

    typedef gsBasis<T> Basis_t;

    /// Coefficient type
    typedef T Scalar_t;

    typedef gsBasis<T> CoordinateBasis;
       
    /// Iterators on coordinate bases
    typedef Basis_t** iterator;
    typedef Basis_t* const* const_iterator;

    typedef gsBasis<T> ** base_iterator;
public:

    /// \brief Default empty constructor
    gsTensorBasis() : Base()
    { m_address = this;}

    /// \brief Copy Constructor
    gsTensorBasis( const gsTensorBasis & o) 
    : Basis_t(o)
    { 
        m_address = this;
    }
    
    /// Assignment opearator
    gsTensorBasis& operator=( const gsTensorBasis & o)
    { 
        this->Base::operator=(o);
        m_address = this;
        return *this;
    }
    
    // Destructor
    ~gsTensorBasis() 
    { 
        m_address = NULL;
    }

    #if EIGEN_HAS_RVALUE_REFERENCES
    gsTensorBasis(gsTensorBasis&& other)
    { gsTensorBasis::operator=(std::forward<gsTensorBasis>(other)); }
    gsTensorBasis & operator=(gsTensorBasis&&other)
    {
        other.m_address = m_address;
        other.m_address = nullptr;
        return *this;
    }
#endif
    bool isValid() const { return static_cast<Basis_t*>(0) != m_address; }
    
    /// \brief Constructor by basis pointers (takes ownership of the
    /// passed bases)
    explicit gsTensorBasis(Base * x) 
    : Base(*x)
    {
        m_address = this; 
        delete x;
        x = NULL;
    }
    
    /// \brief Constructor by basis pointers (takes ownership of the
    /// passed bases)
    explicit gsTensorBasis(base_iterator it) 
    //: Basis_t(*static_cast<Basis_t*>(*it))
    : Basis_t(**it)
    {
        m_address = this; 
        delete *it;
        *it = NULL;
    }
    
public:

    short_t dim() const { return 1;}

    /// Returns a box with the coordinate-wise active functions
    /// \param u evaluation points
    /// \param low lower left corner of the box
    /// \param upp upper right corner of the box   
    void active_cwise(const gsMatrix<T> & u, 
                      gsVector<index_t,1>& low,
                      gsVector<index_t,1>& upp ) const
    { 
        gsMatrix<index_t> act;
        this->active_into(u, act);
        low[0]= act(0,0);
        upp[0]= act(act.size()-1, 0 );
    }

    /// Returns the strides for all dimensions
    void stride_cwise(gsVector<index_t,1> & result) const 
    { 
        result[0] = 1;
    }
    
    /// Get a const-iterator to the beginning of the bases vector
    /// \return an iterator to the beginning of the bases vector
    const_iterator begin() const
    { return &m_address; }
    
    /// Get a const-iterator to the end of the  bases vector
    /// \return an iterator to the end of the  bases vector
    const_iterator end() const
    { return (&m_address)+1; }
    
    /// Get an iterator to the beginning of the  bases vector
    /// \return an iterator to the beginning of the  bases vector
    iterator begin()
    { return &m_address; }
    
    /// Get an iterator to the end of the  bases vector
    /// \return an iterator to the end of the  bases vector
    iterator end()
    { return (&m_address)+1; }
    
    // Unhide/forward gsBasis<T>::size(), since the following
    // overload with size(k) automatically hides it in this class
    // Note that MSVC 2010 produces compilation error if we just
    // do a "using gsBasis<T>::size"
    index_t size() const = 0;

    /// \brief The number of basis functions in the direction of the k-th
    /// parameter component
    index_t size(short_t k) const
    {
        GISMO_UNUSED(k);
        GISMO_ASSERT(k==0, "Invalid direction");
        return this->size();
    }

    /// \brief The number of basis functions in the direction of the k-th
    /// parameter component
    void size_cwise(gsVector<index_t,1> & result) const 
    { result[0] = this->size(); }

    gsVector<int> cwiseDegree() const
    {
        gsVector<int> deg(1);
        deg[0] = this->degree(0);
        return deg;
    }
    
    void swapDirections(const unsigned i, const unsigned j)
    {
        GISMO_UNUSED(i);
        GISMO_UNUSED(j);
        GISMO_ASSERT( static_cast<int>(i) == 0 && static_cast<int>(j) == 0,
                      "Invalid basis components "<<i<<" and "<<j<<" requested" );
    }

    /// Returns all the basis functions with tensor-numbering \a k in direction \a dir 
    gsMatrix<index_t> coefSlice(short_t dir, index_t k) const
    {
        GISMO_UNUSED(dir);
        GISMO_UNUSED(k);
        GISMO_ASSERT(dir == 0, "Invalid direction");
        GISMO_ASSERT(k < this->size(), "Invalid index");
        // return 0 or size()-1
        GISMO_NO_IMPLEMENTATION
     }
    
    inline unsigned index(unsigned i) const
    { return i; }

    /// \todo remove
    inline unsigned index(unsigned , unsigned ) const
    { GISMO_ERROR("The basis is 1D"); }

    /// Returns the stride for dimension dir
    inline unsigned stride(short_t dir) const
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT(dir==0,"Invalid direction");
        return 1; 
    }

    /// Returns the components for a basis on the face \a s 
    void getComponentsForSide(boxSide const &, std::vector<gsBasis<T>*> & rr) const
    { rr.clear(); }

    /// Returns the global index of the basis function created by
    /// components of indices given in the vector v
    inline index_t index(gsVector<index_t,1> const & v) const
    { return v[0]; }
    
    /// Returns the tensor index of the basis function with global index
    /// \a m
    inline gsVector<index_t,1> tensorIndex(const index_t& m) const
    {
        return gsVector<index_t,1>::Constant(1,m);
    }

    const Basis_t& x() const 
    { 
        return *this; 
    }

    Basis_t & component(short_t i)
    {
        GISMO_UNUSED(i);
        GISMO_ASSERT(i==0,"Invalid component requested");
        return *this; 
    }

    const Basis_t & component(short_t i) const
    {
        GISMO_UNUSED(i);
        GISMO_ASSERT(i==0,"Invalid component requested");
        return *this; 
    }
    
private:
    
    /// Keeps the address of the object (for iterator compatibility with d>1)
    Basis_t * m_address;
    
}; // class gsTensorBasis<1,T>

/* ******************************************** */
/* ******************************************** */

template<short_t d, class Basis_t >
inline index_t gsTensorBasis<d,Basis_t>::index(gsVector<index_t,d> const & v) const
{
    index_t ind;

    ind = v(d-1) ;//compute global index in the tensor product
    for ( int i=d-2; i>=0; --i )
        ind = ind * size(i) + v(i) ;
    return ind;
}

template<short_t d, class Basis_t >
inline unsigned gsTensorBasis<d,Basis_t>::index(unsigned i, unsigned j, unsigned k ) const
{
    return size(0) * (size(1) * k + j) + i;
}


template<short_t d, class Basis_t >
inline unsigned gsTensorBasis<d,Basis_t>::stride(short_t dir) const
{
    GISMO_ASSERT( dir>=0 &&  dir< this->dim(), 
                  "Something went wrong with requested direction." );
    unsigned s(1);
    for ( short_t i=0; i<dir; ++i )
        s *= size(i);
    return s;
}


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorBasis.hpp)
#endif


