
/** @file gsBasis.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFunctionSet.h>

#define GISMO_MAKE_GEOMETRY_NEW    \
virtual memory::unique_ptr<gsGeometry<T> > makeGeometry( gsMatrix<T>coefs ) const      \
    { return memory::unique_ptr<gsGeometry<T> >(new GeometryType(*this, give(coefs))); }

namespace gismo
{

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
class gsBasis : public gsFunctionSet<T>
{

public:

    typedef gsFunctionSet<T> Base;

    /// Shared pointer for gsBasis
    typedef memory::shared_ptr< gsBasis > Ptr;

    /// Unique pointer for gsBasis
    typedef memory::unique_ptr< gsBasis > uPtr;

    typedef T Scalar_t;

    static const bool IsRational = false;

    typedef memory::unique_ptr< gsDomainIterator<T> > domainIter;

public:

    gsBasis();

    gsBasis(const gsBasis& other);

    virtual ~gsBasis();

public:

    /*
      Member functions with non-virtual implementations
      (override the _into versions in derived classes).
    */

    const gsBasis<T> & piece(const index_t k) const
    {
        GISMO_ENSURE(0==k, "Single basis is defined on single subdomain, received: "<<k );
        return *this;
    }

    /// Returns the \a i-th basis function as a gsFunction.
    //
    /// Note that the gsBasisFun object only holds a reference to the current
    /// basis, so it is invalidated when the basis is destroyed.
    gsBasisFun<T> function(index_t i) const;

    /// @name Evaluation functions
    /// @{

    /// Evaluate a single basis function \a i at points \a u.
    gsMatrix<T> evalSingle(index_t i, const gsMatrix<T> & u) const
    {
        gsMatrix<T> result;
        this->evalSingle_into(i, u, result);
        return result;
    }

    /// Evaluate a single basis function \a i derivative at points \a u.
    gsMatrix<T> derivSingle(index_t i, const gsMatrix<T> & u) const
    {
        gsMatrix<T> result;
        this->derivSingle_into(i, u, result);
        return result;
    }

    /// Evaluate the second derivative of a single basis function \a i at points \a u.
    gsMatrix<T> deriv2Single(index_t i, const gsMatrix<T> & u) const
    {
        gsMatrix<T> result;
        this->deriv2Single_into(i, u, result);
        return result;
    }

    /** \brief Number of active basis functions at an arbitrary parameter value.
     *
     *  Usually, this is used for getting the active functions on one
     *  element, assuming that this number doesn't change for different
     *  parameters inside the element.
     */
    gsVector<index_t> numActive(const gsMatrix<T> & u) const
    {
        gsVector<index_t> result;
        this->numActive_into(u, result);
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
    gsMatrix<T> evalFunc(const gsMatrix<T> & u, const gsMatrix<T> & coefs) const
    {
        gsMatrix<T> result;
        this->evalFunc_into(u, coefs, result);
        return result;
    }

    /** \brief Evaluate the function described by \a coefs at points \a u,
     * i.e., evaluates a linear combination of coefs x BasisFunctions, into \a result.
     *
     * This function has a default implementation that may be overridden
     * in derived classes for higher performance.
     *
     * \param u     evaluation points as \em N column vectors
     * \param coefs coefficient matrix describing the geometry in this basis, \em n columns
     * \param[out] result  a matrix of size <em>n x N</em> with one function value as a column vector
     *              per evaluation point
     */
    virtual void evalFunc_into(const gsMatrix<T> & u,
                               const gsMatrix<T> & coefs,
                               gsMatrix<T>& result) const;


    /** @brief Evaluate the derivatives of the function described by \a coefs at points \a u.
     *
     * \param u     evaluation points as \em N column vectors
     * \param coefs coefficient matrix describing the geometry in this basis, \em n columns
     * \return   For every column of \a u, the result matrix will contain
     *   one Jacobian matrix of size <em>d * n</em>, such that the total size of
     *   the result is <em>n x (d * n) x N</em>
     */
    gsMatrix<T> derivFunc(const gsMatrix<T> & u, const gsMatrix<T> & coefs) const
    {
        gsMatrix<T> result;
        this->derivFunc_into(u, coefs, result);
        return result;
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
     * The <em>N</em> <b>evaluation points u</b> are given in a gsMatrix of size <em>d</em> x <em>N</em>.
     * Each \em column of \em u represents one evaluation point.
     *
     * The <em>K</em> <b>coefficients coefs</b> are given as a gsMatrix of size <em>K</em> x <em>m</em>.
     * Each \em row of \em coefs represents one coefficient in \f$\mathbb R^m\f$.
     *
     * The gsMatrix <b>result</b> contains the following data:\n
     * For every column of \em u, the corresponding column in the matrix \em result contains
     * the gradients of the \em m components of the function above each other.
     * Hence, the size of \em result is <em>(d*m)</em> x <em>N</em>.
     *
     * <b>Example 1:</b>\n
     * Let \f$f(s,t)\f$ be a bivariate scalar function, \f$f:\mathbb R^2 \to \mathbb R\f$ (i.e., d=2, m=1),
     * and let the evaluation point \f$ u_i\f$ be represented by the <em>i</em>-th column of \em u.\n
     * Then, \em result has the form
     * \f[
     \left( \begin{array}{cccc}
     \partial_s f(u_1) & \partial_s f(u_2) & \ldots &  \partial_t f(u_{N}) \\
     \partial_t f(u_1) & \partial_t f(u_2) & \ldots &  \partial_t f(u_{N})
     \end{array}
     \right)
     \f]
     * <b>Example 2:</b>\n
     * Let \f$f(s,t) = ( f_1(s,t), f_2(s,t), f_3(s,t) )\f$ represent a surface in space, \f$f:\mathbb R^2 \to \mathbb R^3\f$ (i.e., d=2, m=3),
     * and let the evaluation point \f$ u_i\f$ be represented by the <em>i</em>-th column of \em u.\n
     * Then, \em result has the form
     * \f[
     \left( \begin{array}{ccccccc}
     \partial_s f_1(u_1) & \partial_s f_1(u_2) & \ldots &  \partial_s f_1(u_N) \\
     \partial_t f_1(u_1) & \partial_t f_1(u_2) & \ldots &  \partial_t f_1(u_N) \\
     \partial_s f_2(u_1) & \partial_s f_2(u_2) & \ldots &  \partial_s f_2(u_N) \\
     \partial_t f_2(u_1) & \partial_t f_2(u_2) & \ldots &  \partial_t f_2(u_N) \\
     \partial_s f_3(u_1) & \partial_s f_3(u_2) & \ldots &  \partial_s f_3(u_N) \\
     \partial_t f_3(u_1) & \partial_t f_3(u_2) & \ldots &  \partial_t f_3(u_N) \\
     \end{array}
     \right)
     \f]
     *
     * \param[in] u     Evaluation points as \em d x <em>N</em>-matrix.
     * \param[in] coefs Coefficient matrix describing the geometry in this basis as <em>K</em> x <em>m</em>-matrix.\n
     * \em K should equal the size() of the basis, i.e., the number basis functions.
     * \param[in,out] result gsMatrix of size <em>d*m</em> x <em>N</em>, see above for format.
     *
     * where\n
     *\em d is the dimension of the parameter domain\n
     *\em m is the dimension of the physical domain\n
     *\em N is the number of evaluation points\n
     *\em K is the number of coefficients
     */
    virtual void derivFunc_into(const gsMatrix<T> & u,
                                const gsMatrix<T> & coefs,
                                gsMatrix<T>& result ) const;

    /** \brief Evaluate the Jacobian of the function described by \a coefs at points \a u.
        Jacobian matrices are stacked in blocks
     */
    virtual void jacobianFunc_into(const gsMatrix<T> & u,
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
    gsMatrix<T> deriv2Func(const gsMatrix<T> & u, const gsMatrix<T> & coefs) const
    {
        gsMatrix<T> result;
        this->deriv2Func_into(u, coefs, result);
        return result;
    }

    /** \brief Evaluates the second derivatives of the function described by \a coefs at points \a u.
     *
     * ...i.e., evaluates a linear combination of
     * \em coefs * <em>(2nd derivatives of basis functions)</em>, into \a result.
     *
     * <b>Evaluation points \em u</b> are given as gsMatrix of size \em d x \em N, where\n
     * \em d is the dimension of the parameter domain and\n
     * \em N is the number of evaluation points.\n
     * Each column of \em u corresponds to the coordinates of one evaluation point.\n
     * \n
     * The <b>coefficients \em coefs</b> are given as gsMatrix of size \em N x \em n, where\n
     * \em N is the number of points = number of basis functions and\n
     * \em n is the dimension of the physical domain.\n
     * Each row of \em coefs corresponds to the coordinates of one control point.
     *
     * Let the function \f$ f: \mathbb R^3 \to \mathbb R^3\f$ be given by
     * \f[ f = ( f_1, f_2, f_3)^T = \sum_{i=1}^N c_i B_i(x,y,z), \f]
     * where \f$ B_i(x,y,z)\f$ are scalar basis functions and \f$c_i\f$ are the
     * corresponding (<em>m</em>-dimensional) control points.
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
     * \param[out]  result For every column of \a u, a column containing the
     *   second derivatives at the respective point in the format described above.
     *
     * This function has a default implementation that may be overridden
     * in derived classes for higher performance.\n
     * See also deriv2() (the one with input parameter \em coefs).
     *
     */
    virtual void deriv2Func_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const;

    /**
     * @brief Evaluates all derivatives up to order \em n of the function described by \a coefs at points \a u.
     *
     * <b>Evaluation points \em u</b> are given as gsMatrix of size \em d x \em N, where\n
     * \em d is the dimension of the parameter domain and\n
     * \em N is the number of evaluation points.\n
     * Each column of \em u corresponds to the coordinates of one evaluation point.\n
     * \n
     * The <b>coefficients \em coefs</b> are given as gsMatrix of size \em K x \em n, where\n
     * \em K is the number of (active) basis functions (=size()) and\n
     * \em n is the dimension of the physical domain.\n
     * Each row of \em coefs corresponds to the coordinates of one control point.
     *
     * \em result is a std::vector, where the entry <em>result[i]</em> contains the
     * gsMatrix corresponding to the <em>i</em>-th derivatives. The format of the
     * respective entry is as in\n
     * evalFunc_into()\n
     * derivFunc_into()\n
     * deriv2Func_into()\n
     *
     * \todo finish documentation
     *
     * @param[in] u
     * @param[in] coefs
     * @param[in] n
     * @param[out] result
     */
    virtual void evalAllDersFunc_into(const gsMatrix<T> & u,
                                      const gsMatrix<T> & coefs,
                                      const unsigned n,
                                      std::vector< gsMatrix<T> >& result ) const;

    /**
     * @brief Computes the linear combination \em coefs * <em> values( actives ) </em>
     *
     * \todo documentation
     *
     * @param[in] coefs gsMatrix of size \em K x \em m, where \em K should equal size() of the basis (i.e., the number of basis functions).
     * @param[in] actives gsMatrix of size \em numAct x \em numPts
     * @param[in] values gsMatrix of size <em>stride*numAct</em> x \em numPts
     * @param[out] result gsMatrix of size \em stride x \em numPts
     */
    static void linearCombination_into(const gsMatrix<T> & coefs,
                                       const gsMatrix<index_t> & actives,
                                       const gsMatrix<T> & values,
                                       gsMatrix<T> & result);

    /// @}

    inline short_t dim() const {return this->domainDim();}

    /*
      Member functions that need to be redefined in the derived class.
    */

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
    gsMatrix<T> anchors() const
    {
        gsMatrix<T> result;
        this->anchors_into(result);
        return result;
    }

    gsMatrix<T> anchor(index_t i) const
    {
        gsMatrix<T> result;
        this->anchor_into(i, result);
        return result;
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
    virtual void anchor_into(index_t i, gsMatrix<T>& result) const;

    /// Returns the connectivity structure of the basis
    /// The returned mesh has the anchor points as vertices
    virtual void connectivityAtAnchors(gsMesh<T> & mesh) const;

    /// Returns the connectivity structure of the basis
    /// The returned mesh has vertices the rows of matrix \a nodes
    virtual void connectivity(const gsMatrix<T> & nodes, gsMesh<T> & mesh) const;

    /// @name Evaluation functions
    /// @{

    /** \brief Returns the indices of active basis functions at points
     * <em>u</em>, as a list of indices, in <em>result</em>. A
     * function is said to be <em>active</em> in a point if this point
     * lies in the closure of the function's support.
     *
     * \param[in] u  gsMatrix containing evaluation points. Each column represents one evaluation point.
     * \param[out]  result For every column \a i of \a u, a column containing the indices of the
     *   active basis functions at evaluation point <em>u</em>.col(<em>i</em>).
     */
    virtual void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const;

    /// \brief Returns the number of active (nonzero) basis functions at points \a u in \a result.
    virtual void numActive_into(const gsMatrix<T> & u, gsVector<index_t>& result) const;

    /// \brief Returns true if there the point \a u with non-zero
    /// value or derivatives when evaluated at the basis function \a i
    virtual bool isActive(const index_t i, const gsVector<T> & u) const;

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
    virtual gsMatrix<index_t> allBoundary( ) const;

    /// Returns the indices of the basis functions that are nonzero at the domain boundary.
    /// If an offset is provided (the default is zero), it will return the indizes of the basis
    /// functions having this offset to the provided boxSide. Note that the offset cannot be
    /// bigger than the size of the basis in the direction orthogonal to boxSide.
    virtual gsMatrix<index_t> boundaryOffset(boxSide const & s, index_t offset) const;

    /// Returns the indices of the basis functions that are nonzero at the domain boundary as single-column-matrix.
    gsMatrix<index_t> boundary(boxSide const & s) const
    { return this->boundaryOffset(s,0); }

    virtual index_t functionAtCorner(boxCorner const & c) const;

#ifdef __DOXYGEN__
    /// @brief Returns the boundary basis for side s.
    gsBasis<T>::uPtr boundaryBasis(boxSide const & s);
#endif
    GISMO_UPTR_FUNCTION_DEC(gsBasis<T>, boundaryBasis, boxSide const &)

    /// @brief Returns the basis that corresponds to the component
    virtual uPtr componentBasis(boxComponent b) const;

    /// @brief Returns the basis that corresponds to the component
    ///
    /// @param b           The component
    /// @param indices     The row vector where the indices are stored to
    /// @param noBoundary  If true, the transfer matrix does not include parts belonging to lower-order
    ///                    components (i.e., edges without corners or faces without corners and edges)
    virtual uPtr componentBasis_withIndices(boxComponent b, gsMatrix<index_t>& indices, bool noBoundary = true) const;

    /// @brief Returns (a bounding box for) the domain of the whole basis.
    ///
    /// Returns a dx2 matrix, containing the two diagonally extreme
    /// corners of a hypercube.
    virtual gsMatrix<T> support() const;

    /// @brief Returns (a bounding box for) the support of the i-th basis function.
    ///
    /// Returns a dx2 matrix, containing the two diagonally extreme
    /// corners of a hypercube.
    virtual gsMatrix<T> support(const index_t & i) const;

    /// \brief Returns an interval that contains the parameter values in
    /// direction \a dir.
    ///
    /// Returns a 1x2 matrix, containing the two endpoints of the interval.
    gsMatrix<T> supportInterval(index_t dir) const;

    /// @name Evaluation functions
    /// @{

    /**
    \brief Evaluates nonzero basis functions at point \a u into \a result.

    Let...\n
    \a d denote the dimension of the parameter domain.\n
    \a K denote the number of active (i.e., non-zero) basis functions (see active_into()).
    \a N denote the number of evaluation points.\n

    The \a n <b>evaluation points \a u</b> are given in a gsMatrix of size <em>d</em> x <em>N</em>.
    Each column of \a u represents one evaluation point.\n
    \n
    The gsMatrix <b>\a result</b> contains the computed function values in the following form:\n
    Column \a j of \a result corresponds to one evaluation point (specified by the <em>j</em>-th column of \a u).
    The column contains the values of all active functions "above" each other.\n

    For example, for scalar basis functions \a Bi : (x,y,z)-> R, a column represents\n
    (B1, B2, ... , BN)^T,\n
    where the order the basis functions \a Bi is as returned by active() and active_into().

    \param[in] u Evaluation points given as gsMatrix of size <em>d</em> x <em>N</em>.
    See above for details.
    \param[in,out] result gsMatrix of size <em>K</em> x <em>N</em>.
    See above for details.
    */
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// Evaluate the \a i-th basis function at points \a u into \a result.
    virtual void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /**
    \brief Evaluates the first partial derivatives of the nonzero basis function.

    Let \n
    \a d denote the dimension of the parameter domain.\n
    \a K denote the number of active (i.e., non-zero) basis functions (see active_into()).\n
    \a N denote the number of evaluation points.\n

    The \a N <b>evaluation points \a u</b> are given in a gsMatrix of size <em>d</em> x <em>N</em>.
    Each column of \a u represents one evaluation point.\n
    \n
    The gsMatrix <b>\a result</b> contains the computed derivatives in the following form:\n
    Column \a j of \a result corresponds to one evaluation point (specified by the <em>j</em>-th column of \a u).
    The column contains the gradients of all active functions "above" each other.\n

    For example, for scalar basis functions \a \f$B_i : (x,y,z)-> R\f$, a column represents\n
    \f$(dx B_1, dy B_1, dz B_1, dx B_2, dy B_2, dz B_2, ... , dx B_n, dy B_N, dz B_N)^T\f$,\n
    where the order the basis functions \a \f$B_i\f$ is as returned by active() and active_into().

    \param[in] u Evaluation points given as gsMatrix of size <em>d</em> x <em>N</em>.
    See above for details.
    \param[in,out] result gsMatrix of size <em>(K*d)</em> x <em>N</em>.
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
    virtual void derivSingle_into(index_t i,
                                  const gsMatrix<T> & u,
                                  gsMatrix<T>& result ) const;

    /** \brief Evaluate the second derivatives of all active basis function at points \a u.
     *
     * Input parameter \em u is a gsMatrix of size \em d x \em N, where\n
     * \em d is the dimension of the parameter domain and\n
     * \em N is the number of evaluation points.\n
     * Each column of \em u corresponds to the coordinates of one evaluation point.\n
     * \n
     * \em result is a gsMatrix of size <em>(K * d)</em> x \em N, where\n
     * \em K is the number of active basis functions at the evaluation point.\n
     * Each column of \em result corresponds to a column of \em u. It contains the
     * "pure" and the mixed derivatives for each active basis function, "above" each other.\n
     * \n
     ** <b>Example (bivariate):</b> Let \f$B_i(x,y)\f$, <em>d = 2</em> be bivariate basis functions,
     * and let the functions with indices <em>3,4,7, and 8</em> (K = 4) be active at an evaluation
     * point \em u. Then, the corresponding column of \em result represents:\n
     * \f$ (
     * \partial_{xx}\, B_3(u), \partial_{yy}\, B_3(u), \partial_{xy}\, B_3(u),
     * \partial_{xx}\, B_4(u), \partial_{yy}\, B_4(u), \partial_{xy}\, B_4(u),
     * \partial_{xx}\, B_7(u), ... , \partial_{xy}\, B_8(u) )^T \f$\n
     * \n
     * <b>Example (trivariate):</b> Let \f$B_i(x,y,z)\f$, <em>d = 3</em> be trivariate basis functions,
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
    virtual void deriv2Single_into(index_t i,
                                   const gsMatrix<T> & u,
                                   gsMatrix<T>& result ) const;

    /** @brief Evaluate the nonzero basis functions and their derivatives up
     to order \a n at points \a u into \a result.

     The derivatives (the 0-th derivative is the function value) are stored
     in a \em result. \em result is a std::vector, where <em>result[i]</em> is a gsMatrix
     which contains the <em>i</em>-th derivatives.

     The entries in <em>result[0]</em>, <em>result[1]</em>, and <em>result[2]</em> are ordered as in
     eval_into(), deriv_into(), and deriv2_into(), respectively. For <em>i > 2</em>, the
     derivatives are stored in lexicographical order, e.g. for order <em>i = 3</em> and dimension <em>2</em>
     the derivatives are stored as follows:
     \f$ \partial_{xxx}, \, \partial_{xxy}, \, \partial_{xyy}, \, \partial_{yyy}.\, \f$\n

     \param[in] u Evaluation points, each column corresponds to one evaluation point.
     \param[in] n All derivatives up to order \em n are computed and stored
     in \b result.
     \param[in,out] result See above for format.
    */
    virtual void evalAllDers_into(const gsMatrix<T> & u, int n,
                                  std::vector<gsMatrix<T> >& result) const;

    /// @brief Evaluate the basis function \a i and its derivatives up
    /// to order \a n at points \a u into \a result.
    virtual void evalAllDersSingle_into(index_t i, const gsMatrix<T> & u,
                                        int n, gsMatrix<T>& result) const;

    /// @brief Evaluate the (partial) derivative(s) of order \a n the
    /// \a i-th basis function at points \a u into \a result.
    virtual void evalDerSingle_into(index_t i, const
                                    gsMatrix<T> & u, int n,
                                    gsMatrix<T>& result) const;

    /// Compute the Laplacian of all nonzero basis functions at points \a u.
    virtual gsMatrix<T> laplacian(const gsMatrix<T> & u ) const;

    /// @}

    GISMO_UPTR_FUNCTION_PURE(gsBasis, clone)

    /// @brief Create a gsGeometry of proper type for this basis with the
    /// given coefficient matrix.
    virtual memory::unique_ptr<gsGeometry<T> > makeGeometry(gsMatrix<T> coefs) const = 0;

    /// @brief Create an empty basis of the derived type and return a
    /// pointer to it
    virtual typename gsBasis::uPtr create() const;

    /// Return a tensor basis of \a this and \a other
    virtual typename gsBasis::uPtr tensorize(const gsBasis & other) const;

    /// Applicable for rational bases: returns the underlying "source"
    /// (non-rational) basis
    virtual const gsBasis & source () const { return *this; }

    /// Applicable for rational bases: returns the underlying "source"
    /// (non-rational) basis
    virtual gsBasis & source () { return *this; }

    /// Clone the source of this basis in case of rational basis, same
    /// as clone() otherwise
    virtual memory::unique_ptr<gsBasis<T> > makeNonRational() const { return clone(); }

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

    /// @brief The number of elements.
    virtual size_t numElements() const;

    /// @brief The number of elements on side \a s.
    // fixme: default arg = none
    virtual size_t numElements(boxSide const & s) const;

    /// @brief Returns an index for the element which contains point \a u
    virtual size_t elementIndex(const gsVector<T> & u ) const;

    /// @brief For a tensor product basis, return the (const) 1-d
    /// basis for the \a i-th parameter component.
    virtual const gsBasis<T> & component(short_t i) const;

    /// @brief For a tensor product basis, return the 1-d basis for
    /// the \a i-th parameter component.
    virtual gsBasis<T> & component(short_t i);

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
     *
     * \param[in] boxes gsMatrix of size <em>d</em> x <em>n</em>, see above
     * for description of size and meaning.
     *
     * \param[in] refExt Extension to be applied to the refinement boxes
     */
    virtual void refine(gsMatrix<T> const & boxes, int refExt = 0);

    virtual std::vector<index_t> asElements(gsMatrix<T> const & boxes, int refExt = 0) const;

    /** @brief Refinement function, with different sytax for different basis.
     *
     * See documentation of\n
     * gsTensorBasis::refineElements()\n
     * gsHTensorBasis::refineElements()
     *
     */
    virtual void refineElements(std::vector<index_t> const & boxes);

    /** @brief Refine basis and geometry coefficients to levels.
     *
     * Refines the basis as well as the coefficients. The refinement and the format of the
     * input depend on the implementation of refineElements().
     */
    virtual void refineElements_withCoefs(gsMatrix<T> & coefs,std::vector<index_t> const & boxes);

    /// @brief Refine the basis uniformly by inserting \a numKnots new
    /// knots with multiplicity \a mul on each knot span
    virtual void uniformRefine(int numKnots = 1, int mul=1);

    /// @brief Refine the basis uniformly
    ///
    /// The function simultainously updates the vector \a coefs, representing a function
    /// in the bases, such that its new version represents the same function.
    ///
    /// This function is equivalent to
    /// \code
    /// gsSparseMatrix<T,RowMajor> transfer;
    /// basis->uniformRefine_withTransfer(transfer, numKnots, mul);
    /// coefs = transfer * coefs;
    /// \endcode
    ///
    /// \sa gsBasis::uniformRefine
    virtual void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots = 1, int mul = 1);

    /// @brief Refine the basis uniformly
    ///
    /// The function writes a sparse matrix into the variable \a transfer that indicates
    /// how the functions on the coarse grid are represented as linear combinations as fine
    /// grid functions
    ///
    /// \sa gsBasis::uniformRefine
    virtual void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor>& transfer,
                                            int numKnots = 1, int mul = 1);

    /// @brief Coarsen the basis uniformly by removing groups of \a
    /// numKnots consecutive knots, each knot removed \a mul times
    ///
    /// This function is the oposite of gsBasis::uniformRefine
    ///
    /// The execution of
    /// \code
    /// basis->uniformRefine (nKnots, mul)
    /// basis->uniformCoarsen(nKnots);
    /// \endcode
    /// results in no overall change in "basis". However,
    /// \code
    /// basis->uniformCoarsen(nKnots);
    /// basis->uniformRefine (nKnots, mul)
    /// \endcode
    /// is not guaranteed to keep "basis" unchanged.
    ///
    /// \sa gsBasis::uniformRefine
    virtual void uniformCoarsen(int numKnots = 1);

    /// @brief Coarsen the basis uniformly and produce a sparse matrix which
    /// maps coarse coefficient vectors to refined ones
    ///
    /// The function writes a sparse matrix into the variable \a transfer that indicates
    /// how the functions on the coarse grid are represented as linear combinations as fine
    /// grid functions
    ///
    /// \sa gsBasis::uniformCoarsen
    virtual void uniformCoarsen_withTransfer(gsSparseMatrix<T,RowMajor>& transfer,
                                             int numKnots = 1);

    /// @brief Elevate the degree of the basis by the given amount, preserve smoothness.
    virtual void degreeElevate(short_t const & i = 1, short_t const dir = -1);

    /// @brief Reduce the degree of the basis by the given amount, preserve smoothness.
    virtual void degreeReduce(short_t const & i = 1, short_t const dir = -1);

    /// @brief Elevate the degree of the basis by the given amount, preserve knots multiplicity.
    virtual void degreeIncrease(short_t const & i = 1, short_t const dir = -1);

    /// @brief Lower the degree of the basis by the given amount, preserving knots multiplicity.
    virtual void degreeDecrease(short_t const & i = 1, short_t const dir = -1);

    /// @brief Set the degree of the basis (either elevate or
    /// reduce) in order to have degree equal to \a i wrt to each variable
    void setDegree(short_t const& i);

    /// @brief Set the degree of the basis (either increase or
    /// decrecee) in order to have degree equal to \a i
    void setDegreePreservingMultiplicity(short_t const& i);

    /// @brief Elevates the continuity of the basis along element boundaries
    virtual void elevateContinuity(int const & i = 1);

    /// @brief Reduces the continuity of the basis along element boundaries
    virtual void reduceContinuity(int const & i = 1);

    /// Return the gsDomain which represents the parameter domain of
    /// this basis. Currently unused.
    virtual gsDomain<T> * domain() const;

    /// @brief If the basis is of polynomial or piecewise polynomial
    /// type, then this function returns the maximum polynomial degree.
    virtual short_t maxDegree() const;

    /// @brief If the basis is of polynomial or piecewise polynomial
    /// type, then this function returns the minimum polynomial degree.
    virtual short_t minDegree() const;

    /// @brief If the basis is of polynomial or piecewise polynomial
    /// type, then this function returns the total polynomial degree.
    virtual short_t totalDegree() const;

    /// @brief Degree with respect to the i-th variable.
    /// If the basis is a tensor product of (piecewise)
    /// polynomial bases, then this function returns the polynomial
    /// degree of the \a i-th component.
    virtual short_t degree(short_t i) const;

    /// @brief Applies interpolation given the parameter values \a pts
    /// and values \a vals.
    memory::unique_ptr<gsGeometry<T> > interpolateData(gsMatrix<T> const& vals,
                                    gsMatrix<T> const& pts ) const;

    /// @brief Applies interpolation of values \a pts using the
    /// anchors as parameter points.  May be reimplemented in derived
    /// classes with more efficient algorithms. (by default uses
    /// interpolateData(pts,vals)
    virtual memory::unique_ptr<gsGeometry<T> > interpolateAtAnchors(gsMatrix<T> const& vals) const;

    //gsGeometry<T> * projectL2(gsFunction<T> const & func) const;

    /// @brief Computes the collocation matrix w.r.t. points \a u.
    ///
    /// The collocation matrix is a sparse matrix with \em u.cols rows
    /// and \em size() columns. The entry \em (i,j) is the value of
    /// basis function \em j at evaluation point \em i.
    void collocationMatrix(gsMatrix<T> const& u, gsSparseMatrix<T> & result) const;

    /// Reverse the basis
    virtual void reverse();

    /// \brief Computes the indices of DoFs that match on the
    /// interface \a bi. The interface is assumed to be a common face
    /// between this patch and \a other.
    /// The output is two lists of indices \a bndThis and \a bndOther,
    /// with indices that match one-to-one on the boundary \a bi.
    virtual void matchWith(const boundaryInterface & bi, const gsBasis<T> & other,
                           gsMatrix<index_t> & bndThis, gsMatrix<index_t> & bndOther) const;


    /// Get the minimum mesh size, as expected for inverse inequalities
    virtual T getMinCellLength() const;

    /// Get the maximum mesh size, as expected for approximation error estimates
    virtual T getMaxCellLength() const;

protected:

    // inline void getLinearCombination(
    // const gsMatrix<T>         & scalars,
    // const gsMatrix<T> * const & coefs,
    // const gsMatrix<index_t>  & indices,
    // gsMatrix<T>&                result );

}; // class gsBasis


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBasis.hpp)
#endif
