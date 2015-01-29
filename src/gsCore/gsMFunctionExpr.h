/** @file gsMFunctionExpr.h

    @brief Provides declaration of FunctionExpr class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): Mantzaflaris, J. Sogn
*/

#pragma once

# include <gsCore/gsFunction.h>


namespace gismo
{

template<typename T> class gsMFunctionExprPrivate;

/** 
    \brief Class defining a multivariate function given by string
    mathematical expressions.

    \ingroup function
    
    \ingroup Core
*/
    
template<typename T>
class gsMFunctionExpr : public gsFunction<T>
{
    
public:
    
    /// Default empty constructor
    gsMFunctionExpr(); 
    
    /// Constructor by expression strings
    gsMFunctionExpr(const std::vector<std::string> & expression_strings, int ddim);

    /// Constructor by one expression strings
    gsMFunctionExpr(const std::string & expression_string1,
                    int ddim = 2);

    /// Constructor by two expression strings
    gsMFunctionExpr(const std::string & expression_string1, 
                    const std::string & expression_string2,
                    int ddim = 2);

    /// Constructor by three expression strings
    gsMFunctionExpr(const std::string & expression_string1, 
                    const std::string & expression_string2,
                    const std::string & expression_string3,
                    int ddim = 3);

    /// Constructor by four expression strings
    gsMFunctionExpr(const std::string & expression_string1,
                    const std::string & expression_string2,
                    const std::string & expression_string3,
                    const std::string & expression_string4,
                    int ddim);

    /// Constructor by six expression strings
    gsMFunctionExpr(const std::string & expression_string1,
                    const std::string & expression_string2,
                    const std::string & expression_string3,
                    const std::string & expression_string4,
                    const std::string & expression_string5,
                    const std::string & expression_string6,
                    int ddim);

    /// Constructor by nine expression strings
    gsMFunctionExpr(const std::string & expression_string1,
                    const std::string & expression_string2,
                    const std::string & expression_string3,
                    const std::string & expression_string4,
                    const std::string & expression_string5,
                    const std::string & expression_string6,
                    const std::string & expression_string7,
                    const std::string & expression_string8,
                    const std::string & expression_string9,
                    int ddim);

    ~gsMFunctionExpr();
    
    gsMFunctionExpr(const gsMFunctionExpr& other);
    
    gsMFunctionExpr& operator=(const gsMFunctionExpr& other);
    
    gsMFunctionExpr * clone() const
    { return new gsMFunctionExpr(*this); }
    
private:

    // initialize this instance from its m_string member
    void init();
    
public:
    
    virtual int domainDim() const;
    
    virtual int targetDim() const;
    
    /// Evaluate the expression (overrided from gsFunction)
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    
    /// Evaluate the expression for component \a comp in the target dimension (overrided from gsFunction)
    virtual void eval_component_into(const gsMatrix<T>& u, const index_t comp, gsMatrix<T>& result) const;

    /** \brief Evaluator the gradient
     *
     * \param[in] u gsMatrix of evaluation points.
     * Each column of \em u corresponds to one evaluation point.
     * \param[out] result gsMatrix of size
     * \em targetDim x <em>domainDim * u.cols()</em>, where\n
     * one row of \em results corresponds to one component of the function.
     * Within this row, the transposed gradients for the respective points
     * "lie next to each other".
     *
     * Note:\n
     * \em targetDim is the dimension of the result of the function,\n
     * \em domainDim is the dimension of the parameter space.
     *
     * \warning The gradients are computed
     * NUMERICALLY with a default stepsize!!!
     *
     */
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
  
    /// returns the last value computed
    T value(int k) const;
    
    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const;
    
// Data members
private:

  gsMFunctionExprPrivate<T> * my;

}; // class gsMFunctionExpr


//////////////////////////////////////////////////
//////////////////////////////////////////////////


}; // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMFunctionExpr.hpp)
#endif
