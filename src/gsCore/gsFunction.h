/** @file gsFunction.h

    @brief Provides declaration of Function abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/
 
#pragma once

#include <gsCore/gsFunctionSet.h>

namespace gismo
{

/** \brief A function \f$f:\mathbb{R}^n\rightarrow\mathbb{R}^m\f$
 * from a <em>n</em>-dimensional domain to an
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
class gsFunction : public gsFunctionSet<T>
{
public:
    typedef gsFunctionSet<T> Base;

    /// Shared pointer for gsFunction
    typedef memory::shared_ptr< gsFunction > Ptr;

    /// Unique pointer for gsFunction
    typedef memory::unique_ptr< gsFunction > uPtr;
    
    using Base::support;
    using Base::domainDim;
    using Base::targetDim;

    GISMO_UPTR_FUNCTION_NO_IMPLEMENTATION(gsFunction, clone)

    virtual const gsFunction & piece(const index_t k) const
    {
        GISMO_ENSURE(0==k, "Single function of type "<< typeid(*this).name() <<" is defined on single subdomain, received: "<<k<<". Is piece(k) implemented?" );
        return *this; 
    }

    /// Returns the scalar function giving the i-th coordinate of this function
    gsFuncCoordinate<T> coord(const index_t c) const;

    void active_into (const gsMatrix<T>  & u, gsMatrix<index_t> &result) const
    { result.setConstant(1,u.cols(),0); }
    
    /**
        @name Evaluation functions
        @anchor Evaluation_functions

        These functions allow one to evaluate the function as well as its derivatives
        at one or more points in the parameter space. See also \ref func_eval_members.

        @{
    */

    /*
      Member functions with non-virtual implementations
      (override the _into versions in derived classes).
    */

    /** \brief Evaluate the function at points \a u into \a result.
     *
     * Let \em n be the dimension of the source space ( n = domainDim() ).\n
     * Let \em m be the dimension of the image/target space ( m = targetDim() ).\n
     * Let \em N denote the number of evaluation points.
     *
     * \param[in] u gsMatrix of size <em>n</em> x <em>N</em>, where each
     * column of \em u represents one evaluation point.
     * \param[out] result gsMatrix of size <em>m</em> x <em>N</em>, where each
     * column of \em u represents the result of the function at the
     * respective valuation point.
     */
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const = 0;

    /// Evaluate the function for component \a comp in the target dimension at points \a u into \a result.
    virtual void eval_component_into(const gsMatrix<T>& u, 
                                     const index_t comp, 
                                     gsMatrix<T>& result) const;

    /** \brief Evaluate derivatives of the function
     * \f$f:\mathbb{R}^n\rightarrow\mathbb{R}^m\f$
     * at points \a u into \a result.
     *
     * Let \em n be the dimension of the source space ( n = domainDim() ).\n
     * Let \em m be the dimension of the image/target space ( m = targetDim() ).\n
     * Let \em N denote the number of evaluation points.
     *
     * Let \f$ f:\mathbb R^2 \rightarrow \mathbb R^3 \f$, i.e.,
     * \f$ f(x,y) = ( f^{(1)}(x,y), f^{(2)}(x,y), f^{(3)}(x,y) )^T\f$,\n
     * and let
     * \f$ u = ( u_1, \ldots, u_N) = ( (x_1,y_1)^T, \ldots, (x_N, y_N)^T )\f$.\n
     * Then, \em result is of the form
     * \f[
     \left[
     \begin{array}{cccc}
        \partial_x f^{(1)}(u_1) & \partial_x f^{(1)}(u_2)
           & \ldots & \partial_x f^{(1)}(u_N) \\
        \partial_y f^{(1)}(u_1) & \partial_y f^{(1)}(u_2)
           & \ldots & \partial_y f^{(1)}(u_N) \\
        \partial_x f^{(2)}(u_1) & \partial_x f^{(2)}(u_2)
           & \ldots & \partial_x f^{(2)}(u_N) \\
        \partial_y f^{(2)}(u_1) & \partial_y f^{(2)}(u_2)
           & \ldots & \partial_x f^{(2)}(u_N) \\
        \partial_x f^{(3)}(u_1) & \partial_x f^{(3)}(u_2)
           & \ldots & \partial_x f^{(3)}(u_N)\\
        \partial_y f^{(3)}(u_1) & \partial_y f^{(3)}(u_2)
           & \ldots & \partial_y f^{(3)}(u_N)
     \end{array}
     \right]
     \f]
     *
     * \param[in] u gsMatrix of size <em>n</em> x <em>N</em>, where each
     * column of \em u represents one evaluation point.
     * \param[out] result gsMatrix of size <em>(n * m)</em> x <em>N</em>.
     * Each row of \em result corresponds to one component in the target
     * space and contains the gradients for each evaluation point,
     * as row vectors, one after the other (see above for details on the format).
     *
     * \warning By default, gsFunction uses central finite differences
     * with h=0.00001! One must override this function in derived
     * classes to get proper results.
     */
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /** @brief Computes for each point \a u a block of \a result
     * containing the Jacobian matrix
     */
    virtual void jacobian_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /** @brief Computes for each point \a u a block of \a result
     * containing the divergence matrix
     */
    void div_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    gsMatrix<T> jacobian(const gsMatrix<T>& u) const;

    /** @brief Evaluate second derivatives of the function at points \a u into \a result.
     *
     * Let \em n be the dimension of the source space ( n = domainDim() ).\n
     * Let \em m be the dimension of the image/target space ( m = targetDim() ).\n
     * Let \em N denote the number of evaluation points.
     *
     * \param[in] u gsMatrix of size <em>n</em> x <em>N</em>, where each
     * column of \em u represents one evaluation point.
     * \param[out] result gsMatrix of size <em>(S*m)</em> x <em>N</em>,
     * where <em>S=n*(n+1)/2</em>.\n
     * Each column in \em result corresponds to one point (i.e., one column in \em u)\n
     * and contains the following values (for <em>n=3</em>, <em>m=3</em>):\n
     * \f$ (\partial_{xx} f^{(1)}, \partial_{yy} f^{(1)}, \partial_{zz} f^{(1)}, \partial_{xy} f^{(1)},
       \partial_{xz} f^{(1)}, \partial_{yz} f^{(1)}, \partial_{xx} f^{(2)},\ldots,\partial_{yz} f^{(3)} )^T\f$
     * \warning By default uses central finite differences with h=0.00001!
     * One must override this function in derived
     * classes to get proper results.
     */
    virtual void deriv2_into( const gsMatrix<T>& u, gsMatrix<T>& result ) const;

    virtual void hessian_into(const gsMatrix<T>& u, gsMatrix<T>& result,
                              index_t coord = 0) const;

    /// Evaluates the Hessian (matrix of second partial derivatives) of
    /// coordinate \a coord at points \a u.
    virtual gsMatrix<T> hessian(const gsMatrix<T>& u, index_t coord = 0) const
    {
        gsMatrix<T> res;
        hessian_into(u,res,coord);
        return res;
    }

    /// @brief Evaluate the Laplacian at points \a u.
    ///
    /// By default uses central finite differences with h=0.00001
    virtual gsMatrix<T> laplacian( const gsMatrix<T>& u ) const;
  
    /// @}

    /// @brief Computes the L2-distance between this function and the
    /// field and a function \a func
    virtual T distanceL2(gsFunction<T> const &) const
    { GISMO_NO_IMPLEMENTATION }

    /// Newton-Raphson method to find a solution of the equation f(\a
    /// arg) = \a value with starting vector \a arg.
    /// If the point cannot be inverted the corresponding parameter
    /// values will be undefined
    int newtonRaphson(const gsVector<T> & value,
                      gsVector<T> & arg,
                      bool withSupport = true, 
                      const T accuracy = 1e-6,
                      int max_loop = 100,
                      double damping_factor = 1) const;

    gsMatrix<T> argMin(const T accuracy = 1e-6,//index_t coord = 0
                       int max_loop = 100,
                       double damping_factor = 1) const;
    
    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const
    {
        os << "gsFunction.\n"; return os; 
    }

    /**
       @brief Computes map function data

       This function evaluates the functions and their derivatives at
       the points \a InOut.points and writes them in the corresponding
       fields of \a InOut. Which field to write (and what to compute)
       is controlled by the \a InOut.flags (see also gsMapData).
       This is intended for parametrizations only and it works on
       functions sets of cardinality 1 only.

       @param[in,out] InOut
     */
    virtual void computeMap(gsMapData<T> & InOut) const;

    index_t size() const { return 1;}

private:

    template<int mode>
    int newtonRaphson_impl(
        const gsVector<T> & value,
        gsVector<T> & arg, bool withSupport = true,
        const T accuracy = 1e-6, int max_loop = 100,
        double damping_factor = 1, T scale = 1.0) const;

}; // class gsFunction


/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsFunction<T>& b)
{return b.print(os); }


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFunction.hpp)
#endif
