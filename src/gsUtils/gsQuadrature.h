//
// gsQuadrature.h
//
// Functions for creating iterated, tensor, and
// tensor iterated quadrature rules from basic
// 1-D quadrature rules.
//
// Clemens Hofreither
//

#pragma once

#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{

///////////////////////////////////////////////////////////////////
////// Generic functions for mapping reference quadrature rules
///////////////////////////////////////////////////////////////////


/** \brief
    Create an iterated (composite) 1-D quadrature rule for the given
    \a intervals based on \a Rule.

    \tparam T   coefficient type
    \tparam Rule the basic quadrature rule to apply in each interval
    \param[out] nodes   the resulting quadrature nodes
    \param[out] weights the resulting quadrature weights
    \param n    number of quadrature points to use per interval
    \param intervals    \em m+1 scalars defining \em m intervals in which to create the rule
*/
template<class T, void (*Rule)(gsVector<T>&, gsVector<T>&, int, T, T)>
void iteratedQuadratureRule(gsVector<T>& nodes, gsVector<T>& weights, int n,
				const std::vector<T>& intervals);

/// Create an iterated (composite) 1-D quadrature rule for the given
/// \a intervals based on \a Rule (outputs nodes as row vector)
template<class T, void (*Rule)(gsVector<T>&, gsVector<T>&, int, T, T)>
void iteratedQuadratureRule(gsMatrix<T>& nodes, gsVector<T>& weights, int n,
				const std::vector<T>& intervals);

/// Create an n-D tensor rule from the given basic Rule
template<class T, void (*Rule)(gsVector<T>&, gsVector<T>&, int, T, T)>
void tensorQuadratureRule(gsMatrix<T>& ngrid, gsVector<T>& wgrid,
			  const gsVector<int>& numNodes,
			  const gsVector<T>& lower, const gsVector<T>& upper);

/** \brief Create a tensor iterated quadrature rule for the d-cube [\a
 *  lower, \a upper], with \a numNodes[i] nodes per interval in
 *  direction \a i and iterated over \a intervals[i], based on \a Rule
 */
template<class T, void (*Rule)(gsVector<T>&, gsVector<T>&, int, T, T)>
void tensorIteratedQuadratureRule(gsMatrix<T>& ngrid, gsVector<T>& wgrid,
				  const gsVector<int>& numNodes,
				  const std::vector< std::vector<T> >& intervals);


///////////////////////////////////////////////////////////////////
////// Mappers for the Gauss quadrature rule
///////////////////////////////////////////////////////////////////

/** \brief
    Create an iterated 1-D Gauss quadrature rule for the given \a intervals.

    \tparam T   coefficient type
    \param[out] nodes   the resulting quadrature nodes
    \param[out] weights the resulting quadrature weights
    \param n    number of quadrature points to use per interval
    \param intervals    \em m+1 scalars defining \em m intervals in which to create the rule
*/
template<class T>
void iteratedGaussRule(gsMatrix<T>& nodes, gsVector<T>& weights, int n,
			   const std::vector<T>& intervals);

/// Create a tensor Gauss quadrature rule for the d-cube [\a lower, \a
/// upper], with \a numNodes[i] nodes in direction \a i
template<class T>
void tensorGaussRule(gsMatrix<T>& ngrid, gsVector<T>& wgrid, const gsVector<int>& numNodes,
			 const gsVector<T>& lower, const gsVector<T>& upper);

/** \brief Create a tensor iterated Gauss quadrature rule for the
 *  d-cube [\a lower, \a upper], with \a numNodes[i] nodes per
 *  interval in direction \a i and iterated over \a intervals[i]
 */
template<class T>
void tensorIteratedGaussRule(gsMatrix<T>& ngrid, gsVector<T>& wgrid,
				 const gsVector<int>& numNodes,
				 const std::vector< std::vector<T> >& intervals);


/// Create a quasi-uniform tensor iterated Gauss rule with
/// approximately numEvals total evaluations
template<class T>
void uniformGaussRule(gsMatrix<T>& ngrid, gsVector<T>& wgrid, int numEvals,
			  int degree, const gsVector<T>& lower, const gsVector<T>& upper);


///////////////////////////////////////////////////////////////////
////// Reference 1D Quadrature rules
///////////////////////////////////////////////////////////////////


///    Computes the gauss nodes (x) and weights (w)
///    N must be between 1 and 33 or 63/64/65, 127/128/129,
///    255/256/257.
template<class T>
void gsRefGaussRule (gsVector<T>& x, gsVector<T>& w, int n, T a = T(-1), T b = T(1));

/// Computes the isogeometric half-point macroelement rule of ..Hughes
template<class T> inline
void gsHalfPointRule (unsigned const& n, gsVector<T>& x, gsVector<T>& w,
			  T const & a= T(-1), T const & b= T(1) )
{ GISMO_ERROR("Not implemented."); }

///    Computes the quadrature rule for periodic uniform B-splines
template<class T> inline
void gsPeriodicRule (unsigned const& n, gsVector<T>& x, gsVector<T>& w,
			 T const & a= T(-1), T const & b= T(1) )
{ GISMO_ERROR("Not implemented."); }


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsQuadrature.hpp)
#endif
