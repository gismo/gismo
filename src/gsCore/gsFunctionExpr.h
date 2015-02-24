/** @file gsFunctionExpr.h

    @brief Provides declaration of FunctionExpr class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/
 
#pragma once

# include <gsCore/gsFunction.h>

namespace gismo
{


template<typename T> class gsFunctionExprPrivate;

/** 
    @brief Class defining a real function given by a string
    mathematical expression.

    \ingroup function
    \ingroup Core
*/

template<typename T>
class gsFunctionExpr : public gsFunction<T>
{
    
public:
    
    /// Default empty constructor
    gsFunctionExpr(); 
  
    /// Constructor by an expression string
    gsFunctionExpr(std::string expression_string);

    ~gsFunctionExpr();
  
    gsFunctionExpr(const gsFunctionExpr& other);
  
    gsFunctionExpr& operator=(const gsFunctionExpr& other);
  
    gsFunctionExpr * clone() const
    { return new gsFunctionExpr(*this); }

private:

    // initialize this instance from its m_string member
    void init();

    static void stringReplace(std::string& str, const std::string& oldStr, const std::string& newStr)
    {
        size_t pos = 0;
        while((pos = str.find(oldStr, pos)) != std::string::npos)
        {
            str.replace(pos, oldStr.length(), newStr);
            pos += newStr.length();
        }
    }

    
public:

    // Documented in gsFunction class
    virtual int domainDim() const;

    // Documented in gsFunction class
    virtual int targetDim() const             { return 1; }
  
    //set "constants"
    void set_x (T const & v) const;
    void set_y (T const & v) const;
    void set_z (T const & v) const;
    void set_w (T const & v) const;
    void set_u (T const & v) const;
    void set_v (T const & v) const;
  
    void set_x_der (std::string expression_string);
  
    void set_y_der (std::string expression_string);
  
    /// Evaluate the expression (overrided from gsFunction)
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// Evaluate the expression for component \a comp in the target
    /// dimension (overrided from gsFunction)
    virtual void eval_component_into(const gsMatrix<T>& u, 
                                     const index_t comp, 
                                     gsMatrix<T>& result) const;

    /// Evaluate the gradient
    virtual void deriv_into(const gsMatrix<T>& u, 
                            gsMatrix<T>& result) const;
  
    /// Evaluate the HESSIAN matrix
    typename gsFunction<T>::uMatrixPtr hess(const gsMatrix<T>& u, unsigned coord = 0) const;
  
    /// Evaluate the LAPLACIAN
    /// By default uses central finite differences with h=0.00001
    gsMatrix<T> * laplacian(const gsMatrix<T>& u) const;
  
    ///Mixed derivative wrt variables k and j
    gsMatrix<T> * mderiv(const gsMatrix<T>& u, const index_t &k, const index_t &j ) const;
  
  
    /// returns the last value computed
    T value() const;
  
    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const;
  
    /// returns the last value computed
    //T varList() const

  
// Data members
private:

    gsFunctionExprPrivate<T> * my;

}; // class gsFunctionExpr


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFunctionExpr.hpp)
#endif
