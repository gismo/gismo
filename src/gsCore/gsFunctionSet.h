/** @file gsFunctionSet.h

    @brief This is the interface of all objects that computes functions
    on points like gsBasis, gsGeometry and gsFunctions.

    Related to this is the gsFuncData object that is a cache of computed
    values.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, A. Mantzaflaris, J. Vogl
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
//#include <gsCore/gsForwardDeclarations.h> // included by gsLinearAlgebra.h

// The following macros are a design pattern to solve the problem with
// unique::memory pointers as return value of virtual functions in base/derived
// classes. It is expected that a class where this macros are used is derived
// from gsFunctionSet or its derivatives. It assumes that a concrete
// implementation has the suffix "_impl". From outside that class, some can
// call that function by its name and get back a pointer inside a
// memory::unique_ptr (aka. uPtr) of the correct type. If casts are needed
// afterward, use memory::convert_ptr<toType>(from).

// NOTE: One has to add for the DOXYGEN documentation the signature by hand
// and "comment out" the _impl method in the hpp file. This can be done with
// the preprocessor definition __DOXYGEN__. Commenting out is done with @cond
// till @endcond as a doxygen comment.

// Example:

// #ifdef __DOXYGEN__
// /// @briefSome Method
// /// @param a
// /// @param b
// /// @return gsBasis<T>::uPtr
// gsBasis<T>::uPtr someMethod(int a, const int b)
// #endif
//
// GISMO_UPTR_FUNCTION_PURE(gsBasis<T>, someMethod, int, const int)

// /// @cond
// gsBasis<T>* someClass<T>::someMethod_impl(int a, const int b) { .... }
// /// @endcond

// Helper macros for counting arguments, works till highest number in PP_RSEQ_N
// Call with PP_NARG(__VA_ARGS__)
// For upgrade, add values to PP_ARG_N, PP_RSEQ_N and PP_COMMASEQ_N
#define PP_EXPAND(x) x

#define PP_ARG_N(_1,_2,_3,N,...) N
#define PP_RSEQ_N() 3,2,1,0
#define PP_NARG_(...) PP_EXPAND(PP_ARG_N(__VA_ARGS__))
#define PP_COMMASEQ_N()  1,1,0,0
#define PP_COMMA(...) ,
#define PP_HASCOMMA(...) PP_NARG_(__VA_ARGS__,PP_COMMASEQ_N())
#define PP_NARG(...) PP_NARG_HELPER1(PP_HASCOMMA(__VA_ARGS__),PP_HASCOMMA(PP_COMMA __VA_ARGS__ ()),PP_NARG_(__VA_ARGS__, PP_RSEQ_N()))
#define PP_NARG_HELPER1(a,b,N) PP_NARG_HELPER2(a, b, N)
#define PP_NARG_HELPER2(a,b,N) PP_NARG_HELPER3_ ## a ## b(N)
#define PP_NARG_HELPER3_01(N) 0
#define PP_NARG_HELPER3_00(N) 1
#define PP_NARG_HELPER3_11(N) N

// Declaration prototypes. Followed by: { ";" , "= 0;" , "{ ... }" }
#define __DECn(n, type, name, ...)  __DEC ## n(type, name, __VA_ARGS__)
#define __DEC0(type, name, void)    private: virtual type * name##_impl() const
#define __DEC1(type, name, t1)      private: virtual type * name##_impl(t1 n1) const
#define __DEC2(type, name, t1, t2)  private: virtual type * name##_impl(t1 n1, t2 n2) const

// Definition prototypes
#define __DEFn(n, type, name, ...)  __DEF ## n(type, name, __VA_ARGS__)
#define __DEF0(type, name, void)    public:  inline memory::unique_ptr< type > name() const { return memory::unique_ptr< type >(name##_impl()); }
#define __DEF1(type, name, t1)      public:  inline memory::unique_ptr< type > name(t1 n1) const { return memory::unique_ptr< type >(name##_impl(n1)); }
#define __DEF2(type, name, t1, t2)  public:  inline memory::unique_ptr< type > name(t1 n1, t2 n2) const { return memory::unique_ptr< type >(name##_impl(n1, n2)); }

// Declaration of virtual function
// 1st: return type
// 2nd: function name
// 3rd: types of parameter arguments
// Definition of real the implementation in the form "type * name_impl(...)"
// should be done in the corresponding .hpp file.
#define GISMO_UPTR_FUNCTION_DEC(type, name, ...) \
        GISMO_UPTR_FUNCTION_DEC_(PP_NARG(__VA_ARGS__), type, name, __VA_ARGS__)
#define GISMO_UPTR_FUNCTION_DEC_(n, type, name, ...) \
        __DECn(n, type, name, __VA_ARGS__); \
        __DEFn(n, type, name, __VA_ARGS__)

// Declaration and start of definition of virtual function
// 1st: return type
// 2nd: function name
// 3rd: types of parameter arguments
// must be finished with a block of { return type * your implementation }
// don't forget that your in "private:" afterward
#define GISMO_UPTR_FUNCTION_DEF(type, name, ...) \
        GISMO_UPTR_FUNCTION_DEF_(PP_NARG(__VA_ARGS__), type, name, __VA_ARGS__)
#define GISMO_UPTR_FUNCTION_DEF_(n, type, name, ...) \
        __DEFn(n, type, name, __VA_ARGS__) \
        __DECn(n, type, name, __VA_ARGS__)

// Declaration of pure virtual function
// 1st: return type
// 2nd: function name
// 3rd: types of parameter arguments
#define GISMO_UPTR_FUNCTION_PURE(type, name, ...) \
        GISMO_UPTR_FUNCTION_PURE_(PP_NARG(__VA_ARGS__), type, name, __VA_ARGS__)
#define GISMO_UPTR_FUNCTION_PURE_(n, type, name, ...) \
        __DECn(n, type, name, __VA_ARGS__) = 0; \
        __DEFn(n, type, name, __VA_ARGS__)

// Declaration and definition with GISMO_NO_IMPLEMENTATION
// 1st: return type
// 2nd: function name
// 3rd: types of parameter arguments
#define GISMO_UPTR_FUNCTION_NO_IMPLEMENTATION(type, name, ...) \
        GISMO_UPTR_FUNCTION_NO_IMPLEMENTATION_(PP_NARG(__VA_ARGS__), type, name, __VA_ARGS__)
#define GISMO_UPTR_FUNCTION_NO_IMPLEMENTATION_(n, type, name, ...) \
        __DECn(n, type, name, __VA_ARGS__) { GISMO_NO_IMPLEMENTATION } \
        __DEFn(n, type, name, __VA_ARGS__)

// Declaration, definition and implementation of clone function
// 1st: return type
#define GISMO_CLONE_FUNCTION(type) \
        __DEC0(type, clone, void) { return new type(*this); } \
        __DEF0(type, clone, void)

namespace gismo {

/**
   \brief Interface for the set of functions defined on a domain
   (the total number of functions in the set equals to \f$S\f$ )

   All G+Smo objects that evaluate function[s] derive from this object.
   Examples are gsMatrix, gsFunction and gsGeometry.

   The available evaluation procedures are:

   Name of procedure            | Evaluate what
   -----------------------------|-------------------------------
   \c eval   (points)           | value
   \c deriv  (points)           | first derivative(s)
   \c deriv2 (points)           | second derivative(s)
   <!--
   \c div    (points)           | divergence
   \c curl   (points)           | curl / rotor
   \c laplacian (points)        | laplacian 
    -->

   There are two classes of functions for evaluation that differ in the arguments and in
   the returned values. The argument can either be

   1) a matrix;\n
   2) a gsTransform object representing the parametrization of the domain.

   In the first case the the matrix specifies a set of points of a size $N$ (one point per column) and it
   is the caller responsibility to distinguish between functions defined on the physical and
   on the parametric domain. In the second case the evaluator will fetch the appropriate points
   depending on whether it is defined on the parametric or in the physical space.
   In complex situations: in particular when the physical domain of a function is the parametric
   domain of another the automatic mechanism is not guaranteed to work.

   The result is a matrix in which the i-th column contains the requested values
   or derivatives at the i-th point, i.e. the point whose coordinate are the i-th column
   of the input matrix or of the matrix supplied by the gsTransform object.
   On each column the data is grouped in blocks corresponding to different functions,
   so that that if the requested evaluation contains s values
   \f$ v_1, \ldots, v_s \f$ for each pair \f$ (f_i, p_j) \f$ , \f$ i = 1, \ldots, S \f$, \f$ j = 1, \ldots, N \f$ ,
   of (function, point) the output matrix looks like

   \f[
   \left[
   \begin{array}{ccccc}
   v_1(f_1,p_1) & v_1(f_1,p_2) & \ldots & v_1(f_1,p_N)\\
   v_2(f_1,p_1) & v_2(f_1,p_2) & \ldots & v_2(f_1,p_N)\\
   \vdots       & \vdots       &        & \vdots\\
   v_s(f_1,p_1) & v_s(f_1,p_2) & \ldots & v_s(f_1,p_N)\\
   v_1(f_2,p_1) & v_1(f_2,p_2) & \ldots & v_1(f_2,p_N)\\
   \vdots       & \vdots       &        & \vdots\\
   v_s(f_S,p_1) & v_s(f_S,p_2) & \ldots & v_s(f_S,p_N)
   \end{array}
   \right]
   \f]

   Any implementation of the gsFunctionSet interface must:

   1) overload the virtual function int size() that returns the
   number of functions in the set;

   2) overload the needed evaluation functions with gsMatrix
   argument: eval_into, deriv_into, deriv2_into; the one that
   are not implemented will fail at runtime printing, which
   function need implementing to the console

   and possibly:

   3) write optimized versions of div, curl, laplacian, and ect.
   evaluation. By default they are computed by evaluating
   the derivatives and applying the definition.

   \tparam T type for real numbers

   \ingroup  function
*/
template <typename T>
class gsFunctionSet
{
public:
    /// Shared pointer for gsFunctionSet
    typedef memory::shared_ptr< gsFunctionSet > Ptr;

    /// Unique pointer for gsFunctionSet
    typedef memory::unique_ptr< gsFunctionSet > uPtr;

public:

    gsFunctionSet();

    gsFunctionSet(const gsFunctionSet &);

    virtual ~gsFunctionSet();

#ifdef __DOXYGEN__
    /// @brief Clone methode. Produceds a deep copy inside a uPtr.
    uPtr clone();
#endif
    GISMO_UPTR_FUNCTION_NO_IMPLEMENTATION(gsFunctionSet, clone)

    /// @brief Returns the piece(s) of the function(s) at subdomain \a k
    virtual const gsFunctionSet & piece(const index_t) const {return *this;}

    /// @brief Helper which casts and returns the k-th piece of this
    /// function set as a gsFunction
    const gsFunction<T> & function(const index_t k) const;

    /// @brief Helper which casts and returns the k-th piece of this
    /// function set as a gsBasis
    const gsBasis<T> & basis(const index_t k) const;

public:

    /*
     @brief Returns (a bounding box for) the support of the function(s).
    
     Returns either a zero-sized matrix or a dx2 matrix, containing
     the two diagonally extreme corners of a hypercube.

     If the returned matrix is empty, it is understood that the domain
     of definition is the whole of \f$R^{domainDim}\f$
    */
    virtual gsMatrix<T> support() const;

    /**
      @brief Indices of active (non-zero) function(s) for each point.
     
      The columns are sorted in increasing order, if on a point there
      are less active then the number of rows in the result matrix
      (some other point has more actives) then the rest of the column
      is filled with 0s.

      @param u
      @param result
     */
    virtual void active_into (const gsMatrix<T>  & u, gsMatrix<index_t> &result) const;
    
public:
    /**
       @brief Evaluates the function(s).

       For scalar valued functions \f$f_1, \ldots, f_S\f$ from \f$\mathbb{R}^n\rightarrow\mathbb{R}\f$ format is:
       \f[
       \left[
       \begin{array}{ccccc}
       f_1(p_1) & f_1(p_2) & \ldots & f_1(p_N)\\
       f_2(p_1) & f_2(p_2) & \ldots & f_2(p_N)\\
       \vdots       & \vdots       &        & \vdots\\
       f_S(p_1) & f_S(p_2) & \ldots & f_S(p_N)
       \end{array}
       \right]
       \f]
       For vector valued functions function \f$f_1, \ldots, f_S\f$ from \f$\mathbb{R}^n\rightarrow\mathbb{R}^m\f$ the format is:
       \f[
       \left[
       \begin{array}{ccccc}
       f_1^1(p_1) & f_1^{(1)}(p_2) & \ldots & f_1^{(1)}(p_N)\\
       f_1^2(p_1) & f_1^{(2)}(p_2) & \ldots & f_1^{(2)}(p_N)\\
       \vdots     & \vdots     &        & \vdots\\
       f_1^{(m)}(p_1) & f_1^{(m)}(p_2) & \ldots & f_1^{(m)}(p_N)\\
       f_2^{(1)}(p_1) & f_2^{(1)}(p_2) & \ldots & f_2^{(1)}(p_N)\\
       \vdots     & \vdots     &        & \vdots\\
       f_S^{(m)}(p_1) & f_S^{(m)}(p_2) & \ldots & f_S^{(m)}(p_N)
       \end{array}
       \right]
       \f]
       where \f$f^{(i)}_j\f$ is the \f$i\f$-th component of function \f$f_j\f$ of the set.
       @param u
       @param result
    */
    virtual void eval_into      (const gsMatrix<T>  & u, gsMatrix<T> &result) const;

    /**
       @brief First derivatives.

       For scalar valued functions \f$f_1, \ldots, f_S\f$ from \f$\mathbb{R}^n\rightarrow\mathbb{R}\f$ format is:
       \f[
       \left[
       \begin{array}{ccccc}
       \partial_{1}f_1(p_1) & \partial_{1}f_1(p_2) & \ldots & \partial_{1}f_1(p_N)\\
       \partial_{2}f_1(p_1) & \partial_{2}f_1(p_2) & \ldots & \partial_{2}f_1(p_N)\\
       \vdots       & \vdots       &        & \vdots\\
       \partial_{k}f_1(p_1) & \partial_{k}f_1(p_2) & \ldots & \partial_{k}f_1(p_N)\\
       \partial_{1}f_2(p_1) & \partial_{1}f_2(p_2) & \ldots & \partial_{1}f_2(p_N)\\
       \vdots       & \vdots       &        & \vdots\\
       \partial_{k}f_S(p_1) & \partial_{k}f_S(p_2) & \ldots & \partial_{k}f_S(p_N)\\
       \end{array}
       \right]
       \f]
       For vector valued functions function \f$f_1, \ldots, f_S\f$ from \f$\mathbb{R}^n\rightarrow\mathbb{R}^{m}\f$ the format is:
       \f[
       \left[
       \begin{array}{ccccc}
       \partial_{1}f_1^1(p_1) & \partial_{1}f_1^1(p_2) & \ldots & \partial_{1}f_1^1(p_N)\\
       \partial_{2}f_1^1(p_1) & \partial_{1}f_2^1(p_2) & \ldots & \partial_{1}f_2^1(p_N)\\
       \vdots     & \vdots     &        & \vdots\\
       \partial_{k}f_1^1(p_1) & \partial_{k}f_1^1(p_2) & \ldots & \partial_{k}f_1^1(p_N)\\
       \partial_{1}f_1^2(p_1) & \partial_{1}f_1^2(p_2) & \ldots & \partial_{1}f_1^2(p_N)\\
       \vdots     & \vdots     &        & \vdots\\
       \partial_{k}f_1^2(p_1) & \partial_{k}f_1^2(p_2) & \ldots & \partial_{k}f_1^2(p_N)\\
       \partial_{1}f_2^1(p_1) & \partial_{1}f_2^1(p_2) & \ldots & \partial_{1}f_2^1(p_N)\\
       \vdots     & \vdots     &        & \vdots\\
       \partial_{k}f_S^{(m)}(p_1) & \partial_{k}f_S^{(m)}(p_2) & \ldots & \partial_{k}f_S^{(m)}(p_N)
       \end{array}
       \right]
       \f]
       where \f$f^{(i)}_j\f$ is the \f$i\f$-th component of function \f$f_j\f$ of the set.
       @param u
       @param result
    */
    virtual void deriv_into     (const gsMatrix<T>  &  u, gsMatrix<T> &result) const;

    /**
     * @brief Second derivatives.
     *
     For scalar valued functions \f$f_1, \ldots, f_S\f$ from \f$\mathbb{R}^n\rightarrow\mathbb{R}\f$ format is:
     \f[
     \left[
     \begin{array}{ccccc}
     \partial_{1}\partial_{1}f_1(p_1) & \partial_{1}\partial_{1}f_1(p_2) & \ldots & \partial_{1}\partial_{1}f_1(p_N)\\
     \partial_{2}\partial_{2}f_1(p_1) & \partial_{2}\partial_{2}f_1(p_2) & \ldots & \partial_{2}\partial_{2}f_1(p_N)\\
     \vdots       & \vdots       &        & \vdots\\
     \partial_{k}\partial_{k}f_1(p_1) & \partial_{k}\partial_{k}f_1(p_2) & \ldots & \partial_{k}\partial_{k}f_1(p_N)\\
     \partial_{1}\partial_{2}f_1(p_1) & \partial_{1}\partial_{2}f_1(p_2) & \ldots & \partial_{1}\partial_{2}f_1(p_N)\\
     \partial_{1}\partial_{3}f_1(p_1) & \partial_{1}\partial_{3}f_1(p_2) & \ldots & \partial_{1}\partial_{3}f_1(p_N)\\
     \vdots       & \vdots       &        & \vdots\\
     \partial_{1}\partial_{k}f_1(p_1) & \partial_{1}\partial_{k}f_1(p_2) & \ldots & \partial_{1}\partial_{k}f_1(p_N)\\
     \partial_{2}\partial_{3}f_1(p_1) & \partial_{2}\partial_{3}f_1(p_2) & \ldots & \partial_{2}\partial_{3}f_1(p_N)\\
     \vdots       & \vdots       &        & \vdots\\
     \partial_{2}\partial_{k}f_1(p_1) & \partial_{2}\partial_{k}f_1(p_2) & \ldots & \partial_{2}\partial_{k}f_1(p_N)\\
     \partial_{3}\partial_{4}f_1(p_1) & \partial_{3}\partial_{4}f_1(p_2) & \ldots & \partial_{3}\partial_{4}f_1(p_N)\\
     \vdots       & \vdots       &        & \vdots\\
     \partial_{k-1}\partial_{k}f_1(p_1) & \partial_{k-1}\partial_{k}f_1(p_2) & \ldots & \partial_{k-1}\partial_{k}f_1(p_N)\\
     \partial_{1}\partial_{1}f_2(p_1) & \partial_{1}\partial_{1}f_2(p_2) & \ldots & \partial_{1}\partial_{1}f_2(p_N)\\
     \vdots       & \vdots       &        & \vdots\\
     \partial_{k-1}\partial_{k}f_S(p_1) & \partial_{k-1}\partial_{k}f_S(p_2) & \ldots & \partial_{k-1}\partial_{k}f_S(p_N)\\
     \end{array}
     \right]
     \f]
     For vector valued functions function \f$f_1, \ldots, f_S\f$ from \f$\mathbb{R}^n\rightarrow\mathbb{R}^{m}\f$ the format is:
     \f[
     \left[
     \begin{array}{ccccc}
     \partial_{1}\partial_{1}f_1^{(1)}(p_1) & \partial_{1}\partial_{1}f_1^{(1)}(p_2) & \ldots & \partial_{1}\partial_{1}f_1^{(1)}(p_N)\\
     \partial_{2}\partial_{2}f_1^{(1)}(p_1) & \partial_{2}\partial_{2}f_1^{(1)}(p_2) & \ldots & \partial_{2}\partial_{2}f_1^{(1)}(p_N)\\
     \vdots       & \vdots       &        & \vdots\\
     \partial_{k}\partial_{k}f_1^{(1)}(p_1) & \partial_{k}\partial_{k}f_1^{(1)}(p_2) & \ldots & \partial_{k}\partial_{k}f_1^{(1)}(p_N)\\
     \partial_{1}\partial_{2}f_1^{(1)}(p_1) & \partial_{1}\partial_{2}f_1^{(1)}(p_2) & \ldots & \partial_{1}\partial_{2}f_1^{(1)}(p_N)\\
     \partial_{1}\partial_{3}f_1^{(1)}(p_1) & \partial_{1}\partial_{3}f_1^{(1)}(p_2) & \ldots & \partial_{1}\partial_{3}f_1^{(1)}(p_N)\\
     \vdots       & \vdots       &        & \vdots\\
     \partial_{k-1}\partial_{k}f_1{(1)}(p_1) & \partial_{k-1}\partial_{k}f_1^{(1)}(p_2) & \ldots & \partial_{k-1}\partial_{k}f_1^{(1)}(p_N)\\
     \partial_{1}\partial_{1}f_1^{(2)}(p_1) & \partial_{1}\partial_{1}f_1^{(2)}(p_2) & \ldots & \partial_{1}\partial_{1}f_1^{(2)}(p_N)\\
     \vdots       & \vdots       &        & \vdots\\
     \partial_{k-1}\partial_{k}f_1^{(m)}(p_1) & \partial_{k-1}\partial_{k}f_1^{(m)}(p_2) & \ldots & \partial_{k-1}\partial_{k}f_1^{(m)}(p_N)\\
     \vdots       & \vdots       &        & \vdots\\
     \partial_{k-1}\partial_{k}f_S^{(m)}(p_1) & \partial_{k-1}\partial_{k}f_S^{(m)}(p_2) & \ldots & \partial_{k-1}\partial_{k}f_S^{(m)}(p_N)\\
     \end{array}
     \right]
     \f]
     where \f$ f^{(i)}_j\f$ is the \f$i\f$-th component of function \f$ f_j\f$ of the set.

     @param u
     @param result
    */
    virtual void deriv2_into    (const gsMatrix<T>  &  u, gsMatrix<T> &result) const;

    /// @brief Evaluate the nonzero functions and their derivatives up
    /// to order \a n. If n is -1 then no computation is performed.
    virtual void evalAllDers_into(const gsMatrix<T> & u, int n,
                                  std::vector<gsMatrix<T> > & result) const;

    /// Evaluate the function, \see eval_into()
    gsMatrix<T> eval(const gsMatrix<T>& u) const;

    /// Evaluate the derivatives, \see deriv_into()
    gsMatrix<T> deriv(const gsMatrix<T>& u) const;

    /** \brief Evaluates the second derivatives of active (i.e., non-zero) basis at points \a u.
     
        See documentation for deriv2_into()
        (the one without input parameter \em coefs) for details.
        
        \see deriv2_into()
        
        \param[in] u     Evaluation points in columns.
        \return  For every column of \a u, a column containing the second derivatives.
        See documentation for deriv2_into() (the one without input parameter \em coefs) for details.
    */
    gsMatrix<T> deriv2(const gsMatrix<T>& u) const;


    /// @brief Returns the indices of active (nonzero) functions at
    /// points \a u, as a list of indices.
    /// \sa active_into()
    gsMatrix<index_t> active(const gsMatrix<T> & u) const
    {
        gsMatrix<index_t> rvo;
        this->active_into(u, rvo);
        return rvo;
    }

    /*
       @brief Divergence (of vector-valued functions).

       This makes sense only for vector valued functions whose target dimension is a multiple of the domain dimension
       \f[\left[
       \begin{array}{ccccc}
       \text{div} f_1(p_1) & \text{div} f_1(p_2) & \ldots & \text{div} f_1(p_N)\\
       \text{div} f_2(p_1) & \text{div} f_2(p_2) & \ldots & \text{div} f_2(p_N)\\
       \vdots       & \vdots       &        & \vdots\\
       \text{div} f_S(p_1) & \text{div} f_S(p_2) & \ldots & \text{div} f_S(p_N)
       \end{array}
       \right]
       \f]
       If the target dimension is an integer multiple of the domain dimension, \f$\text{div} f(p)\f$ contains
       <em>(target dim)/(domain dim)</em> rows, each row is the divergence of a <em>(domain dim)</em>-tuple of entries of \f$f\f$.
       Equivalently, \f$f\f$ is interpreted as a matrix valued function in column maior standard and
       the divergence corresponding vector valued divergence is returned.

       If the target dimension is not an integer multiple of the domain dimension calling this function results in an error.

       @param u
       @param result
    virtual void div_into       (const gsMatrix<T>  & u, gsMatrix<T> &result) const;
    */


    /*
       @brief Curl (of vector-valued functions).

       This makes sense only for vector-valued functions whose target dimension is a multiple of the domain dimension.
       \f[\left[
       \begin{array}{ccccc}
       \text{curl} f_1(p_1) & \text{curl} f_1(p_2) & \ldots & \text{curl} f_1(p_N)\\
       \text{curl} f_2(p_1) & \text{curl} f_2(p_2) & \ldots & \text{curl} f_2(p_N)\\
       \vdots       & \vdots       &        & \vdots\\
       \text{curl} f_S(p_1) & \text{curl} f_S(p_2) & \ldots & \text{curl} f_S(p_N)
       \end{array}
       \right]
       \f]
       If the target dimension is an integer multiple of the domain dimension, \f$\text{curl} f(p)\f$ contains
       <em>(target dim)/(domain dim)</em> rows, each row is the rotor of a <em>(domain dim)</em>-tuple of entries of \f$f\f$.
       Equivalently, \f$ f\f$ is interpreted as a matrix valued function in column maior standard and
       the rotor corresponding vector valued rotor is returned.

       If the target dimension is not an integer multiple of the domain dimension calling this function results in an error.
       @param u
       @param result
    virtual void curl_into      (const gsMatrix<T> &  u, gsMatrix<T> &result) const;
    */


    /*
       @brief Laplacian (of vector-valued functions).

       \f[\left[
       \begin{array}{ccccc}
       \Delta f_1(p_1) & \Delta f_1(p_2) & \ldots & \Delta f_1(p_N)\\
       \Delta f_2(p_1) & \Delta f_2(p_2) & \ldots & \Delta f_2(p_N)\\
       \vdots       & \vdots       &        & \vdots\\
       \Delta f_S(p_1) & \Delta f_S(p_2) & \ldots & \Delta f_S(p_N)
       \end{array}
       \right]
       \f]
       If the target dimension is different from 1, then each \f$\Delta f_i(p_j)\f$ is a block of target
       dimension rows containing the component wise Laplacian.
       @param u
       @param result
    virtual void laplacian_into (const gsMatrix<T>  &  u, gsMatrix<T> &result) const;
    */



    /**
       @brief Computes function data

       This function evaluates the functions and their derivatives at
       the points \a in and writes them in the corresponding fields of \a out.
       Which field to write (and what to compute) is controlled
       by the \a out.flags (see also gsFuncData).
     
       The input points \a in are expected to be compatible with the
       implementation/representation of the function, i.e. they should
       be points inside the domain of definitition of the function
       
       @param[in] in
       @param[out] out
     */
    virtual void compute(const gsMatrix<T> & in, gsFuncData<T> & out) const;

public:
    /**
       @brief Dimension of the (source) domain.
       @return For \f$f:\mathbb{R}^n\rightarrow\mathbb{R}^m\f$, returns \f$n\f$.
     */
    virtual short_t domainDim () const = 0;

    /**
       @brief Dimension of the target space.
       @return For \f$f:\mathbb{R}^n\rightarrow\mathbb{R}^m\f$, returns \f$m\f$.
     */
    virtual short_t targetDim () const {return 1;}

    /*
      @brief Dimension of domain and target
      @return the pair of integers domainDim() and targetDim()
    */
    std::pair<short_t, short_t> dimensions() const {return std::make_pair(domainDim(),targetDim());}
    
    /**
      @brief size
     
      \warning gsFunction and gsGeometry have size() == 1. This should
      not be confused with the size eg. of gsGeometry::basis(), which
      is the number of basis functions in the basis
     
      @return the size of the function set: the total number of functions
     
     */
    virtual index_t size() const //= 0;
    {GISMO_NO_IMPLEMENTATION}

    /**
       @brief Number of pieces in the domain of definition
     */
    virtual index_t nPieces() const {return 1;}

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const// = 0;
    {
        os << "gsFunctionSet\n";
        return os; 
    }

};

/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsFunctionSet<T>& b)
{return b.print(os); }

} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFunctionSet.hpp)
#endif
