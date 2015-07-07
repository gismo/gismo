/** @file gsBasis.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

#define GISMO_MAKE_GEOMETRY_NEW                                         \
    virtual gsGeometry<T> * makeGeometry( const gsMatrix<T> & coefs ) const \
    { return new GeometryType(*this, coefs); }                          \
    virtual gsGeometry<T> * makeGeometry( gsMovable< gsMatrix<T> > coefs ) const \
    { return new GeometryType(*this, coefs); }


namespace gismo
{

/// Traits class to define basis-related types.
template <class Basis_t, unsigned d = 1>
struct gsTraits
{
// To do: add here generic gsBasis traits..
};

/** 
    \brief A basis represents a family of \em scalar basis functions
    defined over a common parameter domain.

    Among the most important classes in \gismo are those representing
    bases, for instance, a B-spline basis, a tensor product NURBS
    basis, and so on.

    A basis can be viewed as a collection of scalar basis functions

    \f[ B :\ \mathbb R^d \to \mathbb R \f]

    defined over a common parameter domain of dimension *d*. For the
    bases that G+Smo deals with, parameter domains are usually of
    tensor product type and therefore just *d*-dimensional boxes.
        
    Ths gsBasis inteface provides many member functions to evaluate
    their basis functions as well as their derivatives at sets of
    points in the parameter domain.

    The parameter domain has dimension \em d, which can be queried
    at runtime using gsBasis::dim().
    Implementations of gsBasis typically have a dimension which is
    statically known at compile time, but we do not expose the
    dimension as a template parameter of gsBasis in order to make
    dimension-independent code possible to write without excessive templating.

    The main job of a gsBasis is to allow evaluation of its basis functions
    and their derivatives at arbitrary points of the parameter domain.
    All evaluation functions in gsBasis take a
    matrix \em u as an argument which specifies where the basis functions
    should be evaluated. This matrix should have \em d rows, and every
    column specifies one point of the parameter domain at which the
    basis should be evaluated.

    The number of basis functions in a basis is called its \em size
    and can be queried using gsBasis::size().
    
    Most bases have basis functions with local support, i.e., they
    are nonzero only in a small region of the parameter domain.
    A basis function is called \em active at a point \em u if
    \em u is contained in its support. The indices of the active
    basis functions at one or more points can be found by calling
    gsBasis::active().

    \tparam T coefficient type

    \ingroup basis
    \ingroup Core
*/

template<class T>
class gsBasis
{

public: 
    /// Shared pointer for gsBasis
    typedef memory::shared_ptr< gsBasis > Ptr;
    //typedef memory::unique_ptr< gsBasis > LocalPtr;

    typedef T Scalar_t;

    static const bool IsRational = false;

    typedef typename gsMatrix<T>::uPtr        uMatrixPtr;

    typedef memory::auto_ptr< gsDomainIterator<T> > domainIter;

protected:
    
    /// Constructor for derived classes
    gsBasis() { }

    /// Constructor for derived classes
    gsBasis(const gsBasis& other) { }

    /// Assignment for derived classes
    gsBasis & operator=(const gsBasis& other){return *this;}

public:

    /// Destructor
    virtual ~gsBasis() { }

public:

    /*
      Member functions with non-virtual implementations
      (override the _into versions in derived classes).
    */

    /// Returns the \a i-th basis function as a gsFunction.
    //
    /// Note that the gsBasisFun object only holds a reference to the current
    /// basis, so it is invalidated when the basis is destroyed.
    gsBasisFun<T> function(unsigned i) const;

    /// @name Evaluation functions
    /// @{

    /// Evaluate the nonzero basis functions at points \a u.
    uMatrixPtr eval(const gsMatrix<T> & u) const
    // gsMatrix<T> eval(const gsMatrix<T> & u) const
    {
        gsMatrix<T> * result = new gsMatrix<T>;
        this->eval_into(u, *result);
        //return give(result);
        return uMatrixPtr(result);
    }

    /// \brief Evaluate the derivatives of the nonzero basis functions at points \a u.
    ///
    /// See deriv_into() for detailed documentation.
    /// \todo Rename to grad
    uMatrixPtr deriv(const gsMatrix<T> & u) const
    {
        gsMatrix<T> * result = new gsMatrix<T>;
        this->deriv_into(u, *result);
        return uMatrixPtr(result);
    }

    /** \brief Evaluates the second derivatives of active (i.e., non-zero) basis at points \a u.
     *
     * See documentation for deriv2_into()
     * (the one without input parameter \em coefs) for details.
     *
     * \param[in] u     Evaluation points in columns.
     * \return  For every column of \a u, a column containing the second derivatives.
     * See documentation for deriv2_into() (the one without input parameter \em coefs) for details.
     */
    uMatrixPtr deriv2(const gsMatrix<T> & u ) const 
    {
        gsMatrix<T> * result = new gsMatrix<T>;
        this->deriv2_into(u, *result);
        return uMatrixPtr(result);
    }

    /// @brief Evaluate the nonzero basis functions and their derivatives up to
    /// order \a n at points \a u.
    uMatrixPtr evalAllDers(const gsMatrix<T> & u, int n ) const
    {
        gsMatrix<T> *result = new gsMatrix<T>;
        this->evalAllDers_into(u, n, *result);
        return uMatrixPtr(result);
    }

    /// Evaluate a single basis function \a i at points \a u.
    uMatrixPtr evalSingle(unsigned i, const gsMatrix<T> & u) const
    {
        gsMatrix<T> * result = new gsMatrix<T>;
        this->evalSingle_into(i, u, *result);
        return uMatrixPtr(result);
    }

    /// Evaluate a single basis function \a i derivative at points \a u.
    uMatrixPtr derivSingle(unsigned i, const gsMatrix<T> & u) const
    {
        gsMatrix<T> *result = new gsMatrix<T>;
        this->derivSingle_into(i, u, *result);
        return uMatrixPtr(result);
    }

    /// Evaluate the second derivative of a single basis function \a i at points \a u.
    uMatrixPtr deriv2Single(unsigned i, const gsMatrix<T> & u) const
    {
        gsMatrix<T> * result = new gsMatrix<T>;
        this->deriv2Single_into(i, u, *result);
        return uMatrixPtr(result);
    }

    /// @brief Returns the indices of active (nonzero) basis functions at
    /// points \a u, as a list of indices.
    typename gsMatrix<unsigned>::uPtr active(const gsMatrix<T> & u) const
    {
        typename gsMatrix<unsigned>::uPtr result ( new gsMatrix<unsigned> );
        this->active_into(u, *result);
        return result;
    }

    /** \brief Number of active basis functions at an arbitrary parameter value.
     *
     *  Usually, this is used for getting the active functions on one
     *  element, assuming that this number doesn't change for different
     *  parameters inside the element.
     */
    typename gsVector<unsigned>::uPtr numActive(const gsMatrix<T> & u) const
    {
        typename gsVector<unsigned>::uPtr result ( new gsVector<unsigned> );
        this->numActive_into(u, *result);
        return result;
    }

    /// @}



    /** @name Geometry evaluation functions
     *
     * These functions evaluate not the individual basis functions of the basis,
     * but a geometry object which is represented by a coefficient matrix w.r.t.\ this
     * basis object. For the format of the coefficient matrix, see gsGeometry.
     *
     * These functions have default implementations which simply compute the basis
     * function values and perform linear combination, but they may be overridden in
     * derived classes if a higher-performance implementation is possible.
     */
    /// @{

    /** \brief Evaluate the function described by \a coefs at points \a u.
     *
     * This function has a default implementation that may be overridden
     * in derived classes for higher performance.
     * 
     * \param u     evaluation points as \em m column vectors
     * \param coefs coefficient matrix describing the geometry in this basis, \em n columns
     * \return      a matrix of size <em>n x m</em> with one function value as a column vector
     *              per evaluation point
     */
    uMatrixPtr eval(const gsMatrix<T> & u, const gsMatrix<T> & coefs) const
    {
        gsMatrix<T> * result = new gsMatrix<T>;
        this->evalFunc_into(u, coefs, *result);
        return uMatrixPtr(result);
    }

    /** \brief Evaluate the function described by \a coefs at points \a u,
     * i.e., evaluates a linear combination of coefs x BasisFunctions, into \a result.
     *
     * This function has a default implementation that may be overridden
     * in derived classes for higher performance.
     * 
     * \param u     evaluation points as \em m column vectors
     * \param coefs coefficient matrix describing the geometry in this basis, \em n columns
     * \param[out] result  a matrix of size <em>n x m</em> with one function value as a column vector
     *              per evaluation point
     */
    virtual void evalFunc_into(const gsMatrix<T> & u, 
                               const gsMatrix<T> & coefs, 
                               gsMatrix<T>& result) const;


    /** @brief Evaluate the derivatives of the function described by \a coefs at points \a u.
     * 
     * \param u     evaluation points as \em m column vectors
     * \param coefs coefficient matrix describing the geometry in this basis, \em n columns
     * \return   For every column of \a u, the result matrix will contain
     *   one Jacobian matrix of size <em>n x d</em>, such that the total size of
     *   the result is <em>n x (d * m)</em>
     */
    uMatrixPtr deriv(const gsMatrix<T> & u, const gsMatrix<T> & coefs) const
    {
        gsMatrix<T> * result = new gsMatrix<T>;
        this->derivFunc_into(u, coefs, *result);
        return uMatrixPtr(result);
    }

    /** \brief Evaluate the derivatives of the function described by \a coefs at points \a u.
     *
     * Evaluates a linear combination of <em>coefs</em>*<em>BasisFunctionDerivatives</em>, into \a result.
     *
     * This function has a default implementation that may be overridden
     * in derived classes for higher performance.
     * 
     * Let the function \f$f: \mathbb R^d \to \mathbb R^m \f$ be described by the coefficients \em coefs,
     * i.e.,\n
     * each evaluation point is in \f$\mathbb R^d\f$, and\n
     * each coefficient is a point in \f$\mathbb R^m\f$.
     *
     * The <em>n</em> <b>evaluation points u</b> are given in a gsMatrix of size <em>d</em> x <em>n</em>.
     * Each \em column of \em u represents one evaluation point.
     *
     * The <em>k</em> <b>coefficients coefs</b> are given as a gsMatrix of size <em>k</em> x <em>m</em>.
     * Each \em row of \em coefs represents one coefficient in \f$\mathbb R^m\f$.
     *
     * The gsMatrix <b>result</b> contains the following data:\n
     * For every column of \em u, the matrix \em result contains
     * one Jacobian matrix of size <em>m</em> x <em>d</em> "next to each other", such that the total size of
     * \a result is <em>m</em> x <em>(d * n)</em>.
     *
     * <b>Example 1:</b>\n
     * Let \f$f(s,t)\f$ be a bivariate scalar function, \f$f:\mathbb R^2 \to \mathbb R\f$ (i.e., d=2, m=1),
     * and let the evaluation point \f$ u_i\f$ be represented by the <em>i</em>-th column of \em u.\n
     * Then, \em result has the form
     * \f[
     \left( \begin{array}{ccccccc}
     \partial_s f(u_1) & \partial_t f(u_1) & \partial_s f(u_2) & \partial_t f(u_2) & \partial_s f(u_3) & \ldots &  \partial_t f(u_n)
     \end{array}
     \right)
     \f]
     * <b>Example 2:</b>\n
     * Let \f$f(s,t) = ( f_1(s,t), f_2(s,t), f_3(s,t) )\f$ represent a surface in space, \f$f:\mathbb R^2 \to \mathbb R^3\f$ (i.e., d=2, m=3),
     * and let the evaluation point \f$ u_i\f$ be represented by the <em>i</em>-th column of \em u.\n
     * Then, \em result has the form
     * \f[
     \left( \begin{array}{ccccccc}
     \partial_s f_1(u_1) & \partial_t f_1(u_1) & \partial_s f_1(u_2) & \partial_t f_1(u_2) & \partial_s f_1(u_3) & \ldots &  \partial_t f_1(u_n) \\
     \partial_s f_2(u_1) & \partial_t f_2(u_1) & \partial_s f_2(u_2) & \partial_t f_2(u_2) & \partial_s f_2(u_3) & \ldots &  \partial_t f_2(u_n) \\
     \partial_s f_3(u_1) & \partial_t f_3(u_1) & \partial_s f_3(u_2) & \partial_t f_3(u_2) & \partial_s f_3(u_3) & \ldots &  \partial_t f_3(u_n)
     \end{array}
     \right)
     \f]
     *
     * \param[in] u     Evaluation points as \em d x <em>n</em>-matrix.
     * \param[in] coefs Coefficient matrix describing the geometry in this basis as <em>k</em> x <em>m</em>-matrix.
     * \param[in,out] result   For every column of \em u, the matrix \em result will contain
     *   one Jacobian matrix of size <em>m</em> x <em>d</em> "next to each other", such that the total size of
     *   \a result is <em>m</em> x <em>(d * n)</em>
     *
     * where\n
     *\em d is the dimension of the parameter domain\n
     *\em m is the dimension of the physical domain\n
     *\em n is the number of evaluation points\n
     *\em k is the number of coefficients
     */
    virtual void derivFunc_into(const gsMatrix<T> & u, 
                                const gsMatrix<T> & coefs, 
                                gsMatrix<T>& result ) const;


    /** \brief Evaluates the second derivatives of the
     * function described by \a coefs at points \a u.
     *
     * See documentation for deriv2_into() (the one with input parameter \em coefs) for details.
     *
     * \param[in] u     Evaluation points in columns.
     * \param[in] coefs Coefficient matrix describing the geometry in this basis.
     * \return  For every column of \a u, a column containing the second derivatives.
     * See documentation for deriv2_into() (the one with input parameter \em coefs) for details.
     */
    //  second derivatives as follows (eg. for a surface (x,y) --> (f_1, f_2, f_3):
    //     *   ( dxxf_1 dyyf_1 dxyf_1 dxxf_2 dyyf_2 dxy2f_2 dxxf_3 dyyf_3 dxyf_3 )^T
    uMatrixPtr deriv2(const gsMatrix<T> & u, const gsMatrix<T> & coefs) const
    {
        gsMatrix<T> *result = new gsMatrix<T>;
        this->deriv2Func_into(u, coefs, *result);
        return uMatrixPtr(result);
    }


    /** \brief Evaluates the second derivatives of the function described by \a coefs at points \a u.
     *
     * ...i.e., evaluates a linear combination of
     * \em coefs * <em>(2nd derivatives of basis functions)</em>, into \a result.
     *
     * <b>Evaluation points \em u</b> are given as gsMatrix of size \em ParDim x \em n, where\n
     * \em ParDim is the dimension of the parameter domain and\n
     * \em n is the number of evaluation points.\n
     * Each column of \em u corresponds to the coordinates of one evaluation point.\n
     * \n
     * The <b>coefficients \em coefs</b> are given as gsMatrix of size \em N x \em PhysDim, where\n
     * \em N is the number of points = number of basis functions and\n
     * \em PhysDim is the dimension of the physical domain.\n
     * Each row of \em coefs corresponds to the coordinates of one control point.
     *
     * Let the function \f$ f: \mathbb R^3 \to \mathbb R^3\f$ be given by
     * \f[ f = ( f_1, f_2, f_3)^T = \sum_{i=1}^N c_i B_i(x,y,z), \f]
     * where \f$ B_i(x,y,z)\f$ are scalar basis functions and \f$c_i\f$ are the
     * corresponding (<em>PhysDim</em>-dimensional) control points.
     * Then, for each column in \em u, the corresponding column in <b>\em result</b> represents
     * \f[ (
     * \partial_{xx}\ f_1, \partial_{yy}\ f_1, \partial_{zz}\ f_1,
     * \partial_{xy}\ f_1, \partial_{xz}\ f_1, \partial_{yz}\ f_1,
     * \partial_{xx}\ f_2, \partial_{yy}\ f_2, \ldots ,
     * \partial_{xz}\ f_3, \partial_{yz}\ f_3)^T.
     * \f]
     * at the respective evaluation point.
     *
     * \param[in] u     Evaluation points in columns (see above for format).
     * \param[in] coefs Coefficient matrix describing the geometry in this basis.
     * \param[in,out]  result For every column of \a u, a column containing the
     *   second derivatives at the respective point in the format described above.
     *
     * This function has a default implementation that may be overridden
     * in derived classes for higher performance.\n
     * See also deriv2() (the one with input parameter \em coefs).
     *
     */
    virtual void deriv2Func_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const;

    /// @}


    /*
      Member functions that need to be redefined in the derived class.
    */

    /// Returns the dimension \em d of the parameter space.
    virtual int dim() const = 0;

    /**
     * @brief
     * Returns the anchor points that represent the members of the basis.
     * There is exactly one anchor point per basis function.
     *
     * The exact definition of the anchor points depends on the particular basis.
     * For instance, for a Bspline basis these are the Greville abscissae.
     * In general, evaluating a function at the anchor points should provide
     * enough information to interpolate that function using this basis.
     */
    uMatrixPtr anchors() const
    {
        gsMatrix<T> * result = new gsMatrix<T>;
        this->anchors_into(*result);
        return uMatrixPtr(result);
    }

    /**
     * @brief
     * Returns the anchor points that represent the members of the basis in \a result.
     * There is exactly one anchor point per basis function.
     *
     * The exact definition of the anchor points depends on the particular basis.
     * For instance, for a Bspline basis these are the Greville abscissae.
     * In general, evaluating a function at the anchor points should provide
     * enough information to interpolate that function using this basis.
     */
    virtual void anchors_into(gsMatrix<T>& result) const;

    /// Returns the anchor point for member \a i of the basis.
    virtual void anchor_into(unsigned i, gsMatrix<T>& result) const;

    /// Returns the connectivity structure of the basis
    /// The returned mesh has the anchor points as vertices
    virtual void connectivity(gsMesh<T> & mesh) const;

    /// Returns the connectivity structure of the basis
    /// The returned mesh has vertices the rows of matrix \a nodes
    virtual void connectivity(const gsMatrix<T> & nodes, gsMesh<T> & mesh) const;

    /// @name Evaluation functions
    /// @{

    /** \brief Returns the indices of active (non-zero) basis functions
     * at points <em>u</em>, as a list of indices, in <em>result</em>.
     *
     * \param[in] u  gsMatrix containing evaluation points. Each column represents one evaluation point.
     * \param[out]  result For every column \a i of \a u, a column containing the indices of the
     *   active basis functions at evaluation point <em>u</em>.col(<em>i</em>).
     */
    virtual void active_into(const gsMatrix<T> & u, gsMatrix<unsigned>& result) const;

    /// Returns the number of active (nonzero) basis functions at points \a u in \a result.
    virtual void numActive_into(const gsMatrix<T> & u, gsVector<unsigned>& result) const;

    /** \brief Returns the matrix <em>result</em> of active
     * coefficients at points <em>u</em>, each row being one
     * coefficient. The order of the rows is the same as active_into
     * and eval_into functions.
     *
     * \param[in]   u      gsVector containing an evaluation point.
     * \param[in]   coefs  gsMatrix is a coefficient matrix with as many rows as the size of the basis
     * \param[out]  result For every column \a i of \a u, a column containing the indices of the
     *   active basis functions at evaluation point <em>u</em>.col(<em>i</em>).
     */
    virtual void activeCoefs_into(const gsVector<T> & u, const gsMatrix<T> & coefs, 
                                  gsMatrix<T>& result) const;
    
    /// @}

    /// Returns the indices of the basis functions that are nonzero at the domain boundary.
    virtual gsMatrix<unsigned> * allBoundary( ) const;

    /// Returns the indices of the basis functions that are nonzero at the domain boundary.
    /// If an offset is provided (the default is zero), it will return the indizes of the basis
    /// functions having this offset to the provided boxSide. Note that the offset cannot be
    /// bigger than the size of the basis in the direction orthogonal to boxSide.
    virtual gsMatrix<unsigned> * boundaryOffset(boxSide const & s, unsigned offset) const;

    /// Returns the indices of the basis functions that are nonzero at the domain boundary as single-column-matrix.
    gsMatrix<unsigned> * boundary(boxSide const & s) const
    { return this->boundaryOffset(s,0); }

    virtual unsigned functionAtCorner(boxCorner const & c) const;

    /// Returns the boundary basis for side s.
    virtual gsBasis<T> * boundaryBasis(boxSide const & s) const;

    /// @brief Returns (a bounding box for) the domain of the whole basis.
    ///
    /// Returns a dx2 matrix, containing the two diagonally extreme
    /// corners of a hypercube.
    virtual gsMatrix<T> support() const ;

    /// @brief Returns (a bounding box for) the support of the i-th basis function.
    ///
    /// Returns a dx2 matrix, containing the two diagonally extreme
    /// corners of a hypercube.
    virtual gsMatrix<T> support(const unsigned & i) const;

    /// \brief Returns an interval that contains the parameter values in
    /// direction \a dir.
    ///
    /// Returns a 1x2 matrix, containing the two endpoints of the interval.
    gsMatrix<T> supportInterval(unsigned dir) const;

    /// @name Evaluation functions
    /// @{

    /**
    \brief Evaluates nonzero basis functions at point \a u into \a result.
    
    Let...\n
    \a d denote the dimension of the parameter domain.\n
    \a k denote the number of active (i.e., non-zero) basis functions (see active_into()).
    \a n denote the number of evaluation points.\n
    
    The \a n <b>evaluation points \a u</b> are given in a gsMatrix of size <em>d</em> x <em>n</em>.
    Each column of \a u represents one evaluation point.\n
    \n
    The gsMatrix <b>\a result</b> contains the computed function values in the following form:\n
    Column \a j of \a result corresponds to one evaluation point (specified by the <em>j</em>-th column of \a u).
    The column contains the values of all active functions "above" each other.\n
    
    For example, for scalar basis functions \a Bi : (x,y,z)-> R, a column represents\n
    (B1, B2, ... , Bn)^T,\n
    where the order the basis functions \a Bi is as returned by active() and active_into().
    
    \param[in] u Evaluation points given as gsMatrix of size <em>d</em> x <em>n</em>.
    See above for details.
    \param[in,out] result gsMatrix of size <em>k</em> x <em>n</em>.
    See above for details.
    */    
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// Evaluate the \a i-th basis function at points \a u into \a result.
    virtual void evalSingle_into(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /**
    \brief Evaluates the first partial derivatives of the nonzero basis function.
    
    Let...\n
    \a d denote the dimension of the parameter domain.\n
    \a k denote the number of active (i.e., non-zero) basis functions (see active_into()).
    \a n denote the number of evaluation points.\n
    
    The \a n <b>evaluation points \a u</b> are given in a gsMatrix of size <em>d</em> x <em>n</em>.
    Each column of \a u represents one evaluation point.\n
    \n
    The gsMatrix <b>\a result</b> contains the computed derivatives in the following form:\n
    Column \a j of \a result corresponds to one evaluation point (specified by the <em>j</em>-th column of \a u).
    The column contains the gradients of all active functions "above" each other.\n
    
    For example, for scalar basis functions \a Bi : (x,y,z)-> R, a column represents\n
    (dx B1, dy B1, dz B1, dx B2, dy B2, dz B2, ... , dx Bn, dy Bn, dz Bn)^T,\n
    where the order the basis functions \a Bi is as returned by active() and active_into().
    
    \param[in] u Evaluation points given as gsMatrix of size <em>d</em> x <em>n</em>.
    See above for details.
    \param[in,out] result gsMatrix of size <em>(k*d)</em> x <em>n</em>.
    See above for details.
    
    \todo Rename to _ grad_into
    */
    virtual void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    /// \brief Evaluates the (partial) derivatives of the i-th basis function
    /// at points \a u into \a result.
    ///
    /// See deriv_into() for detailed documentation.
    ///
    /// \todo rename grad_into
    virtual void derivSingle_into(unsigned i, 
                                  const gsMatrix<T> & u, 
                                  gsMatrix<T>& result ) const;

    /** \brief Evaluate the second derivatives of all active basis function at points \a u.
     *
     * Input parameter \em u is a gsMatrix of size \em d x \em n, where\n
     * \em d is the dimension of the parameter domain and\n
     * \em n is the number of evaluation points.\n
     * Each column of \em u corresponds to the coordinates of one evaluation point.\n
     * \n
     * \em result is a gsMatrix of size <em>(N * d)</em> x \em n, where\n
     * \em N is the number of active basis functions at the evaluation point.\n
     * Each column of \em result corresponds to a column of \em u. It contains the
     * "pure" and the mixed derivatives for each active basis function, "above" each other.\n
     * \n
     ** <b>Example (bivariate):</b> Let \f$B_i(x,y)$\f be bivariate basis functions,
     * and let the functions with indices <em>3,4,7, and 8</em> be active at an evaluation
     * point \em u. Then, the corresponding column of \em result represents:\n
     * \f$ (
     * \partial_{xx}\, B_3(u), \partial_{yy}\, B_3(u), \partial_{xy}\, B_3(u),
     * \partial_{xx}\, B_4(u), \partial_{yy}\, B_4(u), \partial_{xy}\, B_4(u),
     * \partial_{xx}\, B_7(u), ... , \partial_{xy}\, B_8(u) )^T \f$\n
     * \n
     * <b>Example (trivariate):</b> Let \f$B_i(x,y,z)\f$ be trivariate basis functions,
     * and let the functions with indices <em>3,4,7, and 8</em> be active at an evaluation
     * point \em u. Then, the corresponding column of \em result represents:\n
     * \f$(
     * \partial_{xx}\, B_3(u), \partial_{yy}\, B_3(u), \partial_{zz}\, B_3(u),
     * \partial_{xy}\, B_3(u), \partial_{xz}\, B_3(u), \partial_{yz}\, B_3(u),
     * \partial_{xx}\, B_4(u), ... , \partial_{yz}\, B_8(u) )^T \f$
     *
     * \param[in] u   Evaluation points in columns (see above for format).
     * \param[in,out] result For every column of \a u, a column containing the
     *   second derivatives as described above.
     *
     * See also deriv2() (the one without input parameter \em coefs).
     */
    virtual void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    /// @brief Evaluate the (partial) derivatives of the \a i-th basis function
    /// at points \a u into \a result.
    // hessianSingle_into
    virtual void deriv2Single_into(unsigned i, 
                                   const gsMatrix<T> & u, 
                                   gsMatrix<T>& result ) const;

    /// @brief Evaluate the nonzero basis functions and their derivatives up
    /// to order \a n at points \a u into \a result.
    virtual void evalAllDers_into(const gsMatrix<T> & u, int n, gsMatrix<T>& result) const;

    /// @brief Evaluate the basis function \a i and its derivatives up
    /// to order \a n at points \a u into \a result.
    virtual void evalAllDersSingle_into(unsigned i, const gsMatrix<T> & u, 
                                        int n, gsMatrix<T>& result) const;

    /// @brief Evaluate the (partial) derivative(s) of order \a n the
    /// \a i-th basis function at points \a u into \a result.
    virtual void evalDerSingle_into(unsigned i, const 
                                    gsMatrix<T> & u, int n, 
                                    gsMatrix<T>& result) const;

    /// Compute the Laplacian of all nonzero basis functions at points \a u.
    virtual gsMatrix<T> * laplacian(const gsMatrix<T> & u ) const;

    /// @}

    /// @brief Clone this basis, making a deep copy.
    virtual gsBasis * clone() const = 0;

    /// @brief Create an empty basis of the derived type and return a
    /// pointer to it
    virtual gsBasis * create() const;

    /// Return a tensor basis of \a this and \a other
    virtual gsBasis * tensorize(const gsBasis & other) const;

    /// Clone the source of this basis in case of rational basis, same
    /// as clone() otherwise
    virtual gsBasis<T> * makeNonRational() const { return clone(); }

    /// @brief Create a gsGeometry of proper type for this basis with the
    /// given coefficient matrix.
    virtual gsGeometry<T> * makeGeometry( const gsMatrix<T> & coefs ) const = 0;

    /// @brief Create a gsGeometry of proper type for this basis,
    /// taking ownership of the coefficient matrix.
    virtual gsGeometry<T> * makeGeometry( gsMovable< gsMatrix<T> > coefs ) const = 0;

    /// @brief Create a domain iterator for the computational mesh of
    /// this basis, that points to the first element of the domain
    virtual domainIter makeDomainIterator() const;

    /// @brief Create a boundary domain iterator for the computational
    /// mesh this basis, that points to the first element on the
    /// boundary of the domain
    virtual domainIter makeDomainIterator(const boxSide & s) const;

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const = 0;

    /// Prints the object as a string with extended details.
    virtual std::string detail() const 
    {
        // By default just uses print(..)
        std::ostringstream os;
        print(os);
        return os.str();
    }

    /*
      Member functions that may be implemented or not in the derived class
    */

    /// @brief The number of basis functions in this basis.
    virtual int size() const = 0;

    /// @brief The number of elements.
    virtual int numElements() const;

    /// @brief The number of elements on side \a s.
    virtual int numElements(boxSide const & s) const;

    /// Returns an index for the element which contains point \a u
    virtual int elementIndex(const gsVector<T> & u ) const;

    /// For a tensor product basis, return the 1-d basis for the \a i-th parameter component.
    virtual gsBasis<T>& component(unsigned i) const;

    /** @brief Refine the basis on the area defined by the matrix \a boxes.
     *
     * \a boxes is a <em>d</em> x <em>n</em>-matrix (\a n even), where
     * \a d is the dimension of the parameter domain.\n
     * \a n must be even, and every 2 successive columns in the matrix
     * define a box in the parameter domain (the first column represents the coordinates
     * of the lower corner, the second column the coordinates of the upper corner).
     *
     * <b>Example:</b> The input of the matrix
     * \f[ \left[\begin{array}{cccc} 0 & 0.2 & 0.8 & 1 \\
     0.4 & 0.6 & 0.2 & 0.4
     \end{array} \right] \f]
     * results in refinement of the two boxes
     * \f$[0,0.2]\times[0.4,0.6]\f$ and \f$[0.8,1]\times[0.2,0.4]\f$.
     * \param[in] boxes gsMatrix of size <em>d</em> x <em>n</em>, see above
     * for description of size and meaning.
     */
    virtual void refine(gsMatrix<T> const & boxes);

    /** @brief Refine the basis to levels and in the areas defined by
     * \a boxes with an extension.
     *
     * As of now (03.Oct.2014), only used for hierarchical
     * tensor basis. See gsHTensorBasis for detailed documentation.
     */
    virtual void refine(gsMatrix<T> const & boxes, int refExt);

    /** @brief Refine the basis to levels and in the areas defined by \a boxes.
     *
     * As of now (03.Oct.2014), only used for hierarchical
     * tensor basis. See gsHTensorBasis for detailed documentation.
     */
    virtual void refineElements(std::vector<unsigned> const & boxes);

    /// @brief Refine the basis uniformly by inserting \a numKnots new
    /// knots with multiplicity \a mul on each knot span
    virtual void uniformRefine(int numKnots = 1, int mul=1);

    /// @brief Refine the basis uniformly and adjust the given matrix
    /// of coefficients accordingly
    virtual void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots = 1, int mul=1);

    /// @brief Refine the basis uniformly and produce a sparse matrix which
    /// maps coarse coefficient vectors to refined ones
    virtual void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, 
                                            int numKnots = 1, int mul=1);

    /// @brief Elevate the degree of the basis by the given amount, preserve smoothness.
    virtual void degreeElevate(int const & i = 1, int const dir = -1);

    /// @brief Elevate the degree of the basis by the given amount, preserve knots multiplicity.
    virtual void degreeIncrease(int const & i = 1, int const dir = -1);

    /// @brief Reduce the degree of the basis by the given amount.
    virtual void degreeReduce(int const & i = 1);

    /// @brief Set the degree of the basis.
    virtual void setDegree(int const& i);

    /// @brief Reduces the continuity of the basis along element boundaries
    virtual void reduceContinuity(int const & i = 1);

    /// Return the gsDomain which represents the parameter domain of
    /// this basis. Currently unused.
    virtual gsDomain<T> * domain() const;

    /// @brief If the basis is of polynomial or piecewise polynomial
    /// type, then this function returns the maximum polynomial degree.
    virtual int maxDegree() const;

    /// @brief If the basis is of polynomial or piecewise polynomial
    /// type, then this function returns the minimum polynomial degree.
    virtual int minDegree() const;

    /// @brief If the basis is of polynomial or piecewise polynomial
    /// type, then this function returns the total polynomial degree.
    virtual int totalDegree() const;

    /// @brief Degree with respect to the i-th variable.
    /// If the basis is a tensor product of (piecewise)
    /// polynomial bases, then this function returns the polynomial
    /// degree of the \a i-th component.
    virtual int degree(int i) const;

    /// @brief Applies interpolation given the parameter values \a pts
    /// and values \a vals.
    gsGeometry<T> * interpolateData(gsMatrix<T> const& vals,
                                    gsMatrix<T> const& pts ) const;

    /// @brief Applies interpolation of values \a pts using the
    /// anchors as parameter points.  May be reimplemented in derived
    /// classes with more efficient algorithms. (by default uses
    /// interpolateData(pts,vals)
    virtual gsGeometry<T> * interpolateAtAnchors(gsMatrix<T> const& vals) const;
    
    //gsGeometry<T> * projectL2(gsFunction<T> const & func) const;

    /// @brief Computes the collocation matrix w.r.t. points \a u.
    ///
    /// The collocation matrix is a sparse matrix with \em u.cols rows
    /// and \em size() columns. The entry \em (i,j) is the value of
    /// basis function \em j at evaluation point \em i.
    void collocationMatrix(gsMatrix<T> const& u, gsSparseMatrix<T> & result) const;

    /// Reverse the basis
    virtual void reverse();

protected:

    // inline void getLinearCombination(
    // const gsMatrix<T>         & scalars,
    // const gsMatrix<T> * const & coefs, 
    // const gsMatrix<unsigned>  & indices, 
    // gsMatrix<T>&                result );

}; // class gsBasis


//////////////////////////////////////////////////
//////////////////////////////////////////////////


/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsBasis<T>& b)
{return b.print(os); }


}; // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBasis.hpp)
#endif

