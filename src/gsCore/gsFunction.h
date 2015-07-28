/** @file gsFunction.h

    @brief Provides declaration of Function abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/
 
#pragma once

#include <gsCore/gsForwardDeclarations.h>

namespace gismo
{

/** \brief A function from a <em>d</em>-dimensional domain to an
    <em>m</em>-dimensional image.

    Implementations of gsFunction must at the very least implement the evaluation
    function gsFunction::eval_into(). It is also recommended to specify the
    source and target dimensions by overriding gsFunction::domainDim() and
    gsFunction::targetDim().

    The functions for the derivatives may either be overridden or left
    as the default implementations, which use finite differences.

    \section func_eval_members Evaluation members
    
    All evaluation functions take a matrix \em u as an argument which
    specifies where the function should be evaluated. This matrix
    should have \em d rows, and every column specifies one point of
    the domain at which the function should be evaluated.
    
    Here is an overview over the different evaluation procedures available:
    
    Name of procedure            | Evaluate what
    -----------------------------|-------------------------------
    \c eval(u)                   | value
    \c deriv(u)                  | first derivative(s)
    \c deriv2(u)                 | second derivative(s)
    
    All evaluation functions also provide a version suffixed with \c _into
    which takes a matrix reference as an additional output parameter into which
    the result will be stored.

    \tparam T arithmetic type

    \ingroup function
    \ingroup Core
*/

template<class T>
class gsFunction
{
public:

    typedef typename memory::auto_ptr<gsMatrix<T> >        uMatrixPtr;

public:

    /// Clones the function object, making a deep copy.
    virtual gsFunction * clone() const;

    /// Destructor
    virtual ~gsFunction() { }

public:

    /// The dimension of the function domain, i.e., the source space.
    virtual int domainDim() const;

    /// The dimension of the target space.
    virtual int targetDim() const;

    /// @brief Returns (a bounding box for) the support of the function.
    ///
    /// Returns a dx2 matrix, containing the two diagonally extreme
    /// corners of a hypercube.
    virtual gsMatrix<T> support() const;


    /**
        @name Evaluation functions
        @anchor Evaluation_functions

        These functions allow to evaluate the function as well as its derivatives
        at one or more points in the parameter space. See also \ref func_eval_members.

        @{
    */

    /*
      Member functions with non-virtual implementations
      (override the _into versions in derived classes).
    */

    /// Evaluate the function, see eval_into() for details.
    uMatrixPtr eval(const gsMatrix<T>& u) const;

    /// Evaluate the derivatives, see deriv_into() for details.
    //Returns a matrix of size targetDim() x (domainDim() * u.cols()), gradients are stored as row vectors
    uMatrixPtr deriv(const gsMatrix<T>& u) const;

    // Evaluate the gradient
    //uMatrixPtr grad(const gsMatrix<T>& u) const        { return deriv(u); /* should return gradients as column vectors for 1D functions */ }
  
    /// Evaluate the second derivatives
    uMatrixPtr deriv2(const gsMatrix<T>& u) const;

    /** \brief Evaluate the function at points \a u into \a result.
     *
     * Let \em d be the dimension of the source space ( d = domainDim() ).\n
     * Let \em D be the dimension of the image/target space ( D = targetDim() ).\n
     * Let \em n denote the number of evaluation points.
     *
     * \param[in] u gsMatrix of size <em>d</em> x <em>n</em>, where each
     * column of \em u represents one evaluation point.
     * \param[out] result gsMatrix of size <em>D</em> x <em>n</em>, where each
     * column of \em u represents the result of the function at the
     * respective valuation point.
     */
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const = 0;

    /// Evaluate the function for component \a comp in the target dimension at points \a u into \a result.
    virtual void eval_component_into(const gsMatrix<T>& u, 
                                     const index_t comp, 
                                     gsMatrix<T>& result) const;

    /** \brief Evaluate derivatives of the function at points \a u into \a result.
     *
     * Let \em d be the dimension of the source space ( d = domainDim() ).\n
     * Let \em D be the dimension of the image/target space ( D = targetDim() ).\n
     * Let \em n denote the number of evaluation points.
     *
     * Let \f$ f:\mathbb R^2 \rightarrow \mathbb R^3 \f$, i.e.,
     * \f$ f(x,y) = ( f_1(x,y), f_2(x,y), f_3(x,y) )^T\f$,\n
     * and let
     * \f$ u = ( u_1, \ldots, u_n) = ( (x_1,y_1)^T, \ldots, (x_n,y_n)^T )\f$.\n
     * Then, \em result is of the form
     * \f[
     \left[
     \begin{array}{ccccc}
        \partial_x f_1(u_1) & \partial_y f_1(u_1) & \partial_x f_1(u_2)
           & \ldots & \partial_y f_1(u_n) \\
        \partial_x f_2(u_1) & \partial_y f_2(u_1) & \partial_x f_2(u_2)
           & \ldots & \partial_y f_2(u_n) \\
        \partial_x f_3(u_1) & \partial_y f_3(u_1) & \partial_x f_3(u_2)
           & \ldots & \partial_y f_3(u_n)
     \end{array}
     \right]
     \f]
     *
     * \param[in] u gsMatrix of size <em>d</em> x <em>n</em>, where each
     * column of \em u represents one evaluation point.
     * \param[out] result gsMatrix of size <em>D</em> x <em>(d*n)</em>.
     * Each row of \em result corresponds to one component in the target
     * space and contains the gradients for each evaluation point,
     * as row vectors, one after the other (see above for details on the format).
     *
     * \warning By default, gsFunction uses central finite differences
     * with h=0.00001! One must override this function in derived
     * classes to get proper results.
     */
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    // temporary: towards sanitizing the format of deriv_into ( point u.col(i) --> result.col(i) )
    void jacobian_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    uMatrixPtr jacobian(const gsMatrix<T>& u) const;

    void grad_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

 
    /** @brief Evaluate second derivatives of the function at points \a u into \a result.
     *
     * Let \em d be the dimension of the source space ( d = domainDim() ).\n
     * Let \em D be the dimension of the image/target space ( D = targetDim() ).\n
     * Let \em n denote the number of evaluation points.
     *
     * \param[in] u gsMatrix of size <em>d</em> x <em>n</em>, where each
     * column of \em u represents one evaluation point.
     * \param[out] result gsMatrix of size <em>(S*D)</em> x <em>n</em>,
     * where <em>S=d*(d+1)/2</em>.\n
     * Each column in \em result corresponds to one point (i.e., one column in \em u)\n
     * and contains the following values (for <em>d=3</em>, <em>D=3</em>):\n
     * \f$ (\partial_{xx} f_1, \partial_{yy} f_1, \partial_{zz} f_1, \partial_{xy} f_1,
       \partial_{xz} f_1, \partial_{yz} f_1, \partial_{xx} f_2,\ldots,\partial_{yz} f_3 )^T\f$
     * \warning By default uses central finite differences with h=0.00001!
     * One must override this function in derived
     * classes to get proper results.
     */
    virtual void deriv2_into( const gsMatrix<T>& u, gsMatrix<T>& result ) const;
  
    /// Evaluates the Hessian (matrix of second partial derivatives) of
    /// coordinate \a coord at points \a u.
    virtual uMatrixPtr hess(const gsMatrix<T>& u, unsigned coord = 0) const;

    /// @brief Evaluate the Laplacian at points \a u.
    ///
    /// By default uses central finite differences with h=0.00001
    virtual gsMatrix<T>* laplacian( const gsMatrix<T>& u ) const;
  
    /// @}

    /// @brief Computes the L2-distance between this function and the
    /// field and a function \a func
    virtual T distanceL2(gsFunction<T> const & func) const
    { GISMO_NO_IMPLEMENTATION }

    /// Newton-Raphson method to find a solution of the equation f(\a
    /// arg) = \a value with starting vector \a arg.
    /// If the point cannot be inverted the corresponding parameter
    /// values will be undefined
    int newtonRaphson(const gsVector<T> & value,
                      gsVector<T> & arg,
                      bool withSupport = true, 
                      const T accuracy = 1e-6,
                      int max_loop = 100 ) const;
    
    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const
    {
        os << "gsFunction.\n"; return os; 
    }

}; // class gsFunction


//////////////////////////////////////////////////
//////////////////////////////////////////////////

/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsFunction<T>& b)
{return b.print(os); };


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFunction.hpp)
#endif
