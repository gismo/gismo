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
#include <gsUtils/gsSortedVector.h>

namespace gismo
{

// Forward declaration
template<unsigned d, class Basis_t > class gsTensorBasis;

/** 
 *  Abstract base class for tensor product bases.
 *
 *   \param d dimension of the parameter domain
 *   \param Basis_t type of the coordinate-wise bases
 *
 *   \ingroup basis
 */

template<unsigned d, class Basis_t >
class gsTensorBasis : public gsBasis<typename Basis_t::Scalar_t >  
{
public: 
    /// Coefficient type
    typedef typename Basis_t::Scalar_t Scalar_t;
    typedef typename Basis_t::Scalar_t T;

    /// Coordinate basis type
    typedef Basis_t CoordinateBasis;

    /// Associated geometry type
    typedef typename gsTraits<Basis_t,d>::TensorGeometryType GeometryType;

    /// Associated Boundary basis type
    typedef typename gsTraits<Basis_t,d>::TensorBoundaryType BoundaryBasisType;

    /// Dimension of the parameter domain
    static const int Dim = d;

    /// Iterators on coordinate bases
    typedef Basis_t** iterator;
    typedef Basis_t* const* const_iterator;

public:
    
    /// Default empty constructor
    gsTensorBasis()
    {
        for (unsigned i = 0; i < d; ++i)
            m_bases[i] = new Basis_t;
    }
    
    explicit gsTensorBasis( Basis_t* x)
    { gsWarn<< "Invalid constructor: Tried to make tensor basis of parDim="<<d<<" using only one component.\n"; }

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
    gsTensorBasis( std::vector<Basis_t* > const & bb ) ;
    
    /// Copy Constructor
    gsTensorBasis( const gsTensorBasis & o);
    
    /// Assignment opearator
    gsTensorBasis& operator=( const gsTensorBasis & o);
    
    // Destructor
    ~gsTensorBasis() 
    { 
        for (unsigned i = 0; i < d; ++i)
            delete m_bases[i];
    }
    
public:

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

    // Returns the dimension of the basis
    int dim() const { 
        return Dim; 
    }

    /// Returns the number of elements in the basis
    int size() const {
        unsigned r=1;
        for (unsigned i = 0; i < d; ++i)
            r *= m_bases[i]->Basis_t::size();
        return r; 
    }

    // Look at gsBasis class for a description
    int numElements() const 
    {
        int nElem = m_bases[0]->numElements();
        for (unsigned dim = 1; dim < d; ++dim)
            nElem *= m_bases[dim]->numElements();
        return nElem;
    }

    // Look at gsBasis class for a description
    int elementIndex(const gsVector<T> & u ) const
    {
        GISMO_ASSERT( u.rows() == d, "Wrond vector dimension");

        int ElIndex = m_bases[d-1]->elementIndex(u[d-1]);
        for ( int i=d-2; i>=0; --i )
            ElIndex = ElIndex * m_bases[i]->numElements() 
                    + m_bases[i]->elementIndex(u[i]);

        return ElIndex;        
    }

    /// Returns the number of elements (component wise)
    void numElements_cwise(gsVector<unsigned>& result) const
    {
        result.resize(d);
        for (unsigned dim = 0; dim < d; ++dim)
            result(dim) = static_cast<unsigned>(m_bases[dim]->numElements());
    }

    /// Returns the anchors (graville absissae) that represent the members of the basis
    void anchors_into(gsMatrix<T>& result) const;

    /**
     * \brief Returns the indices of active (non-zero) basis functions
     * at points <em>u</em>, as a list of indices, in <em>result</em>.
     *
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
    void active_into(const gsMatrix<T> & u, gsMatrix<unsigned>& result) const;

    /// Returns a box with the coordinate-wise active functions
    /// \param u evaluation points
    /// \param low lower left corner of the box
    /// \param upp upper right corner of the box
    void active_cwise(const gsMatrix<T> & u, gsVector<unsigned,d>& low, 
                      gsVector<unsigned,d>& upp ) const;

    // Look at gsBasis class for documentation 
    virtual void connectivity(const gsMatrix<T> & nodes, gsMesh<T> & mesh) const;

    /// Returns the indices of the basis functions that touch the domain
    /// boundary
    gsMatrix<unsigned> * boundary( ) const ;

    /// Returns the indices of the basis functions that touch the domain
    /// boundary
    gsMatrix<unsigned> * boundary(boxSide const & s ) const ;

    /// Returns a basis for the face s gsTensorBasis<T,d-1,Basis_t>
    BoundaryBasisType * boundaryBasis(boxSide const & s ) const ;

    /// Returns a bounding box for the basis' domain
    gsMatrix<T> support() const ;

    /// Returns a bounding box for the support of the ith basis function
    gsMatrix<T> support( const unsigned & i ) const ;

    /// Evaluates the non-zero basis functions (and optionally their
    /// first k derivatives) at value u into result
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const ;

    /// Evaluate the i-th basis function at all columns of the matrix
    /// (or vector) u
    void evalSingle_into(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result) const ;

    /// Evaluate an element of the space given by coefs at points u
    virtual void eval_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const;

    /// Evaluate the nonzero basis functions and their derivatives up to
    /// order n at all columns of u
    virtual void evalAllDers_into(const gsMatrix<T> & u, int n, gsMatrix<T>& result) const;

    /// Evaluates the gradient the non-zero basis functions at value u.
    virtual void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    /// Evaluates the second derivatives of the non-zero basis functions at value u.
    virtual void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    /// Evaluate the i-th basis function derivative at all columns of
    virtual void derivSingle_into(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result) const ;

    virtual void deriv2Single_into(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result) const ;

    /// Evaluates the (partial) derivatives of an element given by coefs at (the columns of) u.
    void deriv_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const ;

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
    virtual gsGeometry<T> * interpolate(gsMatrix<T> const& vals) const;

    /// Interpolates values on a tensor-grid of points, given in
    /// tensor form (d coordinate-wise vectors). Samples \a vals
    /// should be ordered as the tensor-basis coefficients
    gsGeometry<T> * interpolateGrid(gsMatrix<T> const& vals,
                                    std::vector<gsMatrix<T> >const& grid) const;

    /// Prints the object as a string, pure virtual function of gsTensorBasis.
    virtual std::ostream &print(std::ostream &os) const = 0;

    // Look at gsBasis class for documentation 
    virtual void uniformRefine(int numKnots = 1)
    {
        for (unsigned j = 0; j < d; ++j)
            m_bases[j]->uniformRefine(numKnots);
    }

    // Look at gsBasis class for documentation 
    void refineElements(std::vector<unsigned> const & elements)
    {
        gsSortedVector<unsigned> elIndices[d];
        unsigned tmp, mm;

        // Get coordinate wise element indices
        for ( typename  std::vector<unsigned>::const_iterator
              it = elements.begin(); it != elements.end(); ++it )
        {
            mm = *it;
            for (unsigned i = 0; i<d; ++i )
            {
                const unsigned nEl_i = m_bases[i]->numElements();
                tmp = mm % nEl_i;
                mm = (mm - tmp) / nEl_i;
                elIndices[i].push_sorted_unique(tmp);
            }
        }

        // Refine in each coordinate
        // Element refinement propagates along knot-lines
        for (unsigned i = 0; i<d; ++i )
        {
            m_bases[i]->refineElements(elIndices[i]);
        }
    }

    /// Refine the basis uniformly and perform knot refinement for the
    /// given coefficient vector
    void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots=1);

    /// Refine the basis uniformly and produce a sparse matrix which
    /// maps coarse coefficient vectors to refined ones
    void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots=1);

    // Look at gsBasis class for documentation 
    virtual void degreeElevate(int const & i = 1)
    { 
        for (unsigned j = 0; j < d; ++j)
            m_bases[j]->degreeElevate(i);
    }
    
    // Look at gsBasis class for documentation 
    virtual void degreeElevateComponent(unsigned dir,int const & i = 1)
    {
        GISMO_ASSERT( static_cast<int>(dir) < this->dim(),
                      "Invalid basis component requested" );
        m_bases[dir]->degreeElevate(i);
    }

    // Look at gsBasis class for documentation 
    virtual void degreeReduce(int const & i = 1)
    { 
        for (unsigned j = 0; j < d; ++j)
            m_bases[j]->degreeReduce(i);
    }

    // Look at gsBasis class for documentation 
    virtual void setDegree(int const & i)
    { 
        for (unsigned j = 0; j < d; ++j)
            m_bases[j]->setDegree(i);
    }


//////////////////////////////////////////////////
// Additional members for Tensor Basis
//////////////////////////////////////////////////

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
    int size(int k) const { return m_bases[k]->Basis_t::size(); }

    /// The number of basis functions in the direction of the k-th parameter component
    void size_cwise(gsVector<unsigned,d> & result) const 
    { 
        for ( unsigned k = 0; k!=d; ++k )
            result[k] = m_bases[k]->Basis_t::size(); 
    }

    /// Returns all the basis functions with tensor-numbering \param k in direction \param dir
    ///
    /// ## Detailed explanation: ##
    /// Tensor-numbering in N-variate tensor-product basis means that each basis function
    /// is assigned an identifier (i_0, i_1, ..., i_{N-1}).
    /// This function returns indices of basis functions with i_dir = k
    /// and the returned indices are numbering of the basis functions in the basis
    /// (i.e., 0,1, ..., basis.size() ).
    /// ## Example: ##
    /// Bivariate tensor-product basis functions have tensor numbering (a,b).
    /// Calling dir=0, k=1 gives all functions with tensor-numbering (1,b).
    /// Calling dir=1, k=3 gives all functions with tensor-numbering (a,3).
    typename gsMatrix<unsigned>::uPtr slice(int dir, int k) const;

    /// Returns the degree of the basis wrt variable \a i 
    int degree(int i) const { 
        return m_bases[i]->degree(); 
    }

    int degree() const 
    { 
        int td = 0;
        // take maximum of coordinate bases degrees
        for (unsigned k=0; k!=d; ++k)
            td = math::max(td, m_bases[k]->degree());
        return td;
    }
    
    int maxDegree() const 
    { 
        int td = m_bases[0]->degree();
        // take maximum of coordinate bases degrees
        for (unsigned k=1; k!=d; ++k)
            td = math::max(td, m_bases[k]->degree());
        return td;
    }
    
    int minDegree() const 
    { 
        int td = m_bases[0]->degree();
        // take minimum of coordinate bases degrees
        for (unsigned k=1; k!=d; ++k)
            td = math::min(td, m_bases[k]->degree());
        return td;
    }
    
    int totalDegree() const 
    { 
        int td = 0;
        for (unsigned k=0; k!=d; ++k)
            td = + m_bases[k]->degree();
        return td;
    }


    /// Returns the global index of the basis function created by
    /// components of indices i,j,k (for 2d or 3d only)
    inline unsigned index(unsigned i, unsigned j, unsigned k=0) const;

    /// Returns the stride for dimension dir
    inline unsigned stride(int dir) const;

    /// Returns the strides for all dimensions
    void stride_cwise(gsVector<int> & result) const 
    { 
        result.resize(d);
        result[0] = 1;
        for ( unsigned i=1; i != d; ++i )
            result[i] = result[i-1] * m_bases[i]->Basis_t::size();
    }

    /// Returns the global index of the basis function created by
    /// components of indices given in the vector v
    inline unsigned index(gsVector<unsigned,d> const & v) const;
    //  inline unsigned index(gsVector<unsigned>         & v) const;

    /// Returns the tensor index of the basis function with global index \a m.
    inline gsVector<unsigned, d> tensorIndex(const unsigned& m) const {
        gsVector<unsigned, d> ind;
        int mm = m;
        for (unsigned i = 0; i<d; ++i )
        {
            ind(i)= mm % size(i);
            mm -= ind(i);
            mm /= size(i);
        }
        return ind;
    }

    /// Returns true iff the basis function with multi-index \ind is on
    /// the boundary
    inline bool indexOnBoundary(const gsVector<unsigned, d> & ind) const 
    {
        for ( unsigned i = 0; i < d; ++i )
            if ( ind[i] == static_cast<unsigned>(size(i)-1) )
                return false;
        return ( (ind.array() > 0).all() );
    }

    /// Returns true iff the basis function indexed \a m is on the
    /// boundary
    inline bool indexOnBoundary(const unsigned m) const 
    {
        return ( indexOnBoundary( tensorIndex(m) ) );
    }

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

    Basis_t& component(unsigned dir) const 
    { 
        GISMO_ASSERT( static_cast<int>(dir) < this->dim(), 
                      "Invalid basis component requested" );
        return *m_bases[dir];
    }

    inline int trueSize(int k) const { return m_bases[k]->Basis_t::trueSize(); }

// Data members
private:

    Basis_t* m_bases[d];

}; // class gsTensorBasis


/** 
 *  Class for a Tensor product spline space of dimension 1.
 *  This specialization is mainly for compatibility.
 *
 *   \param T coefficient type
 *   \param Basis_t type of the (single) coordinate-wise bases
 */
    
template<class Basis_t>
class gsTensorBasis<1,Basis_t> : public Basis_t
{
public: 
    typedef Basis_t Base;

    /// Coefficient type
    typedef typename Basis_t::Scalar_t Scalar_t;
    typedef typename Basis_t::Scalar_t T;
    
    /// Coordinate basis type
    typedef Basis_t CoordinateBasis;

    /// Associated Boundary basis type
    typedef typename Basis_t::BoundaryBasisType BoundaryBasisType;
    
    /// Iterators on coordinate bases
    typedef Basis_t** iterator;
    typedef Basis_t* const* const_iterator;
    
public:

    /// Default empty constructor
    gsTensorBasis() : Basis_t()
    { m_address = this;}

    explicit gsTensorBasis(Basis_t * x) : Basis_t(*x)
    { 
        m_address = this;
        delete x;
    }
    
    gsTensorBasis( Basis_t* x,  Basis_t*  y)
    { gsWarn<<"Invalid constructor.\n"; }
    
    gsTensorBasis( Basis_t* x,  Basis_t* y, Basis_t* z )
    { gsWarn<<"Invalid constructor.\n"; }
    
    gsTensorBasis( Basis_t* x,  Basis_t* y, Basis_t* z, Basis_t* w )
    { gsWarn<<"Invalid constructor.\n"; }
    
    // Constructor by basis pointers (takes ownership of the passed bases)
    gsTensorBasis( std::vector<Basis_t* > const & bb ) 
    : Basis_t(*bb[0])
    { 
        GISMO_ASSERT(bb.size() == 1, "Invalid number of bases in 1D constructor");
        m_address = this; 
        delete bb.front();
    }
    
    
    /// Copy Constructor
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
    
public:
    
    /// Returns a box with the coordinate-wise active functions
    /// \param u evaluation points
    /// \param low lower left corner of the box
    /// \param upp upper right corner of the box   
    void active_cwise(const gsMatrix<T> & u, gsVector<unsigned,1>& low, 
                      gsVector<unsigned,1>& upp ) const
    { 
        gsMatrix<unsigned> act;
        this->active_into(u, act);
        low[0]= act(0,0);
        upp[0]= act(act.size()-1, 0 );
    }
    
    /// Get a const-iterator to the beginning of the bases vector
    /// \return an iterator to the beginning of the bases vector
    const_iterator begin() const
    { return &m_address; }
    
    /// Get a const-iterator to the end of the  bases vector
    /// \return an iterator to the end of the  bases vector
    const_iterator end() const
    { return &(m_address)+1; }
    
    /// Get an iterator to the beginning of the  bases vector
    /// \return an iterator to the beginning of the  bases vector
    iterator begin()
    { return &m_address; }
    
    /// Get an iterator to the end of the  bases vector
    /// \return an iterator to the end of the  bases vector
    iterator end()
    { return &(m_address)+1; }
    
    /// The number of basis functions in the direction of the k-th parameter component
    int size(int k) const { return Basis_t::size(); }
    
    int size() const {return Basis_t::size(); }
    
    /// The number of basis functions in the direction of the k-th parameter component
    void size_cwise(gsVector<unsigned,1> & result) const 
    { result[0] = Basis_t::size(); }
    
    /// Returns all the basis functions with tensor-numbering \a k in direction \a dir 
    typename gsMatrix<unsigned>::uPtr slice(int dir, int k) const
    {
        GISMO_ASSERT(dir == 0, "Invalid direction");
        GISMO_ASSERT(k < size(), "Invalid index");
        // return 0 or size()-1
        GISMO_NO_IMPLEMENTATION
     }
    
    inline unsigned index(unsigned i) const
    { return i; }

    /// \todo remove
    inline unsigned index(unsigned i, unsigned j, unsigned k=0) const
    { return i; }
    
    /// Returns the stride for dimension dir
    inline unsigned stride(int dir) const 
    { 
        GISMO_ASSERT(dir==0,"Invalid direction");
        return 1; 
    }

    /// Returns the global index of the basis function created by
    /// components of indices given in the vector v
    inline unsigned index(gsVector<unsigned,1> const & v) const
    { return v[0]; }
    
    /// Returns the tensor index of the basis function with global index
    /// \a m
    inline gsVector<unsigned, 1> tensorIndex(const unsigned& m) const 
    {
        return gsVector<unsigned, 1>::Constant(1,m);
    }

    Basis_t& x() const 
    { 
        return *m_address; 
    }
    
    Basis_t& component(unsigned i) const 
    {
        GISMO_ASSERT(i==0,"Invalid component requested");
        return *m_address; 
    }
    
private:
    
    /// Keeps the address of the object (for compatibility with d>1)
    Basis_t * m_address;
    
}; // class gsTensorBasis1D
    

//////////////////////////////////////////////////
//////////////////////////////////////////////////

template<unsigned d, class Basis_t >
inline unsigned gsTensorBasis<d,Basis_t>::index(gsVector<unsigned,d> const & v) const
{
    unsigned ind;

    ind = v(d-1) ;//compute global index in the tensor product
    for ( int i=d-2; i>=0; --i )
        ind = ind * trueSize(i) + v(i) ;
    return ind;
}

template<unsigned d, class Basis_t > 
inline unsigned gsTensorBasis<d,Basis_t>::index(unsigned i, unsigned j, unsigned k ) const
{
    return size(0) * (size(1) * k + j) + i;
}


template<unsigned d, class Basis_t >
inline unsigned gsTensorBasis<d,Basis_t>::stride(int dir) const
{
    GISMO_ASSERT( dir>=0 &&  dir< this->dim(), "Something went with requested direction." );
    unsigned s(1);
    for ( int i=0; i<dir; ++i )
        s *= size(i);
    return s;
}


}; // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorBasis.hpp)
#endif


