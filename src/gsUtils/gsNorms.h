//
// gsNorms.h
//
// Functions for computing various norms
// and distances between functions in norms.
//
// Clemens Hofreither
//

#pragma once

#include <gsCore/gsForwardDeclarations.h>

// check which ones are really needed:

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFunction.h>
#include <gsCore/gsConstantFunction.h>
#include <gsCore/gsField.h>
#include <gsCore/gsGeometry.h>
#include <gsUtils/gsPointGrid.h>
#include <gsCore/gsMultiBasis.h>

#include <gsCore/gsDomainIterator.h>
#include <gsTensor/gsTensorDomainIterator.h>
#include <gsTensor/gsTensorDomainBoundaryIterator.h>

#include <gsCore/gsGeometryEvaluator.h>

/** @file gsNorms.h
 * Global functions for norm computations.
 * \ingroup Utils
 */

namespace gismo
{

/*
 * In all these norm functions,
 *  parametrized=false means the function is defined on the physical geometry,
 *  parametrized=true  means the function is defined on the parametric domain.
 */

 struct boundaryInterface;
/*
 * L2 norms
 */

template <typename T>
T computeL2Norm(const gsGeometry<T>& geo, const gsFunction<T>& u, bool isParametrized_u, int numEvals=1000);

template <typename T>
T computeL2Norm(const gsField<T>& u, int numSamples=1000);

template <typename T>
T computeL2Distance(const gsGeometry<T>& geo, const gsFunction<T>& u, bool isParametrized_u, const gsFunction<T>& v, bool isParametrized_v, int numEvals=1000);

template <typename T>
T computeL2Distance(const gsField<T>& u, const gsFunction<T>& v, bool isParametrized_v, int numEvals=1000);

template <typename T>
T computeL2Distance(const gsField<T>& u, const gsField<T>& v, int numEvals=1000);


/// \brief L2-distance between \em func and \em v.
///
/// ...i.e., the L2-norm of the difference ( func - v ).
///
/// \param[in] patch Geometry information on a single patch.
/// \param[in] func Function (parameterized, i.e., an "isogeometric" function).
/// \param[in] v Function defined on the parameter domain or physical domain.
/// Specified by \em v_isParam.
/// \param[in] v_isParam Specifies whether \em v is defined on the
/// parameter domain (<em>v_isParam = true</em>) or on the
/// physical domain (<em>v_isParam = false</em>).
///
/// \returns Scalar, L2-norm of the difference.
///
/// \ingroup Utils
///
template <typename T>
T igaL2Distance(const gsGeometry<T>& patch, 
                const gsGeometry<T>& func, 
                const gsFunction<T>& v, 
                bool v_isParam = false);
template <typename T>
T igaL2Distance(const gsGeometry<T>& patch,
                const gsFunction<T>& func,
                const gsFunction<T>& v,
                const gsBasis<T>& B,
                bool v_isParam);


/// \brief Element-wise L2-distance between \em func and \em v
/// for isogeometric parameterizations.
///
/// Returns a matrix \em Errs of the <b>squared</b> element-wise L2-norms of the difference ( func - v ).
///
/// The entries of the matrix \em Errs are as follows:\n
/// <em>Errs(i,0)</em> = squared L2-norm of (func-v) on element \em i.\n
/// <em>Errs(i,1,...,d)</em> = coordinates of lower corner of the element.\n
/// <em>Errs(i,d+1,...,2*d)</em> = coordinates of upper corner of the element.\n
/// ...where \em d denotes the dimension of the parameter domain.
///
/// \param[in] patch Geometry information on a single patch.
/// \param[in] func Function (parameterized, i.e.,
/// an "isogeometric" function)
/// \param[in] v Function defined on the parameter domain
/// or physical domain.
/// Specified by \em v_isParam.
/// \param[in] v_isParam Specifies whether \em v is defined on the
/// parameter domain (<em>v_isParam = true</em>) or on the
/// physical domain (<em>v_isParam = false</em>).
///
/// \returns Errs gsMatrix of size <em>NE</em> x <em>( 2*d+1 )</em>, where\n
/// \em NE is the number of elements.\n
/// \em d is the dimension of the parameter domain.\n
/// See above for format of entries.
///
/// \ingroup Utils
///
template <typename T>
gsMatrix<T> igaL2DistanceEltWiseSq(const gsGeometry<T>& patch,
                const gsGeometry<T>& func,
                const gsFunction<T>& v,
                bool v_isParam = false);


/// \brief Element-wise L2-distance between \em u and \em v.
///
/// Returns a vector \em Errs.
/// The entry <em>Errs[k]</em> is a gsMatrix containing the
/// <b>squared</b> element-wise L2-norms of the difference ( u - v ) on patch \em k.
///
/// Let <em>Errs_k = Errs[k]</em>.
/// The entries of the matrix \em Errs_k are as follows:\n
/// <em>Errs_k(i,0)</em> = squared L2-norm of (func-v) on element \em i.\n
/// <em>Errs_k(i,1,...,d)</em> = coordinates of lower corner of the element.\n
/// <em>Errs_k(i,d+1,...,2*d)</em> = coordinates of upper corner of the element.\n
/// ...where \em d denotes the dimension of the parameter domain.
///
/// \param[in] u gsField with geometry and function information, possibly multi-patch.
/// \param[in] v Function defined on the parameter domain
/// or physical domain (which is specified by \em v_isParam).
/// \param[in] v_isParam Specifies whether \em v is defined on the
/// parameter domain (<em>v_isParam = true</em>) or on the
/// physical domain (<em>v_isParam = false</em>).
///
/// \returns Errs gsVector of length \em N, where \em N is the number of patches (a.k.a. subdomains) of \em u.\n
/// Each entry of the vector is a gsMatrix of size <em>NE</em> x <em>( 2*d+1 )</em>, where\n
/// \em NE is the number of elements.\n
/// \em d is the dimension of the parameter domain.\n
/// See above for format of entries.
///
/// \ingroup Utils
///
template <typename T>
gsVector< gsMatrix<T> > igaFieldL2DistanceEltWiseSq(const gsField<T>& u, const gsFunction<T>& v, bool v_isParam);

/// \brief L2-distance between \em u and \em v.
///
/// ...i.e., the L2-norm of the difference ( func - v ).
///
/// \param[in] u gsField with geometry and function information, possibly multi-patch.
/// \param[in] v Function defined on the parameter domain or physical domain
/// (which is specified by \em v_isParam).
/// \param[in] v_isParam Specifies whether \em v is defined on the
/// parameter domain (<em>v_isParam = true</em>) or on the
/// physical domain (<em>v_isParam = false</em>).
///
/// \returns L2-norm of the difference.
///
/// \ingroup Utils
///
template <typename T>
T igaFieldL2Distance(const gsField<T>& u, 
                     const gsFunction<T>& v, 
                     bool v_isParam = false);
template <typename T>
T igaFieldL2Distance(const gsField<T>& u,
                     const gsFunction<T>& v,
                     const gsMultiBasis<T>& B,
                     bool v_isParam = false);
template <typename T>
T igaH1DistanceOnElt( const typename gsGeometryEvaluator<T>::uPtr & geoEval ,
                      const gsFunction<T> & func,
                      const gsFunction<T>& v,
                      const bool & v_isParam,
                      const typename gsBasis<T>::domainIter & domIt,
					  const gsQuadRule<T> & quRule);


// Auxiliary functions for igaL2Distance() and igaL2DistanceEltWiseSq().
template <typename T>
T igaL2DistanceOnElt( const typename gsGeometryEvaluator<T>::uPtr & geoEval ,
                      const typename gsGeometryEvaluator<T>::uPtr & funcEval,
                      const gsFunction<T>& v,
                      const bool & v_isParam,
                      const typename gsBasis<T>::domainIter & domIt,
					  const gsQuadRule<T> & quRule);
template <typename T>
T igaL2DistanceOnElt( const typename gsGeometryEvaluator<T>::uPtr & geoEval ,
                      const gsFunction<T> & func,
                      const gsFunction<T>& v,
                      const bool & v_isParam,
                      const typename gsBasis<T>::domainIter & domIt,
					  const gsQuadRule<T> & quRule);

 /*
  * Maximum norms
  */

template <typename T>
T computeMaximumNorm(const gsFunction<T>& f, const gsVector<T>& lower, const gsVector<T>& upper, int numSamples=1000);

template <typename T>
T computeMaximumNorm(const gsGeometry<T>& geo, const gsFunction<T>& u, bool isParametrized_u, int numSamples=1000);

template <typename T>
T computeMaximumNorm(const gsField<T>& u, int numSamples=1000);

template <typename T>
T computeMaximumDistance(const gsFunction<T>& f1, const gsFunction<T>& f2, const gsVector<T>& lower, const gsVector<T>& upper, int numSamples=1000);

template <typename T>
T computeMaximumDistance(const gsGeometry<T>& geo, const gsFunction<T>& u, bool isParametrized_u, const gsFunction<T>& v, bool isParametrized_v, int numSamples=1000);

template <typename T>
T computeMaximumDistance(const gsField<T>& u, const gsFunction<T>& v, bool isParametrized_v, int numSamples=1000);


/*
 * H^1 semi-norm
 */

template <typename T>
T igaH1Distance(const gsGeometry<T>& patch, 
                const gsGeometry<T>& func, 
                const gsFunction<T>& v, 
                bool v_isParam = false);
template <typename T>
T igaH1Distance(const gsGeometry<T>& patch,
                const gsFunction<T>& func,
                const gsFunction<T>& v,
                const gsBasis<T>& B,
                bool v_isParam);

template <typename T>
T igaFieldH1Distance(const gsField<T>& u,
                     const gsFunction<T>& v,
                     const gsMultiBasis<T>& B,
                     bool v_isParam = false);

template <typename T>
gsMatrix<T> igaH1DistanceEltWiseSq(const gsGeometry<T>& patch,
                const gsGeometry<T>& func,
                const gsFunction<T>& v,
                bool v_isParam = false);

/*
 * Auxiliary function for igaH1Distance() and igaH1DistanceEltWiseSq().
 */
template <typename T>
T igaH1DistanceOnElt( const typename gsGeometryEvaluator<T>::uPtr & geoEval ,
                      const typename gsGeometryEvaluator<T>::uPtr & funcEval,
                      const gsFunction<T>& v,
                      const bool & v_isParam,
                      const typename gsBasis<T>::domainIter & domIt,
                      const int d,
					  const gsQuadRule<T> & quRule);

template <typename T>
T igaFieldH1Distance(const gsField<T>& u, 
                     const gsFunction<T>& v, 
                     bool v_isParam = false);

template <typename T>
gsVector< gsMatrix<T> > igaFieldH1DistanceEltWiseSq(const gsField<T>& u, const gsFunction<T>& v, bool v_isParam);


/*
 * DG norm
 */

template <typename T>
T igaDGDistanceJump(const gsGeometry<T>& patch1, const gsGeometry<T>& patch2,
		    const gsGeometry<T>& func1,  const gsGeometry<T>& func2, // approximati solution
                    const gsFunction<T>& v1, const gsFunction<T>& v2,	// exact solution
		    const boundaryInterface & bi, // interface
		    const T mu,
                    bool v_isParam);

template <typename T>
T igaFieldDGDistance(const gsField<T>& u, const gsFunction<T>& v, bool v_isParam= false);

} // namespace gismo



#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsNorms.hpp)
#endif

