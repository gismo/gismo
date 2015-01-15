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

    /// Evaluate the function
    uMatrixPtr eval(const gsMatrix<T>& u) const;

    /// Evaluate the derivatives. Returns a matrix of size targetDim() x (domainDim() * u.cols()), gradients are stored as row vectors
    uMatrixPtr deriv(const gsMatrix<T>& u) const;

    // Evaluate the gradient
    //uMatrixPtr grad(const gsMatrix<T>& u) const        { return deriv(u); /* should return gradients as column vectors for 1D functions */ }
  
    /// Evaluate the second derivatives
    uMatrixPtr deriv2(const gsMatrix<T>& u) const;

    /// Evaluate the function at points \a u into \a result.
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const = 0;

    /// Evaluate the function for component \a comp in the target dimension at points \a u into \a result.
    virtual void eval_component_into(const gsMatrix<T>& u, 
                                     const index_t comp, 
                                     gsMatrix<T>& result) const;

    /// @brief Evaluate derivatives of the function at points \a u into \a result.
    ///
    /// By default uses central finite differences with h=0.00001
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
  
    /// @brief Evaluate second derivatives of the function at points \a u into \a result.
    ///
    /// By default uses central finite differences with h=0.00001
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
    {
        GISMO_NO_IMPLEMENTATION
            }

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
