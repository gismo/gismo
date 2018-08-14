#include <gsCore/gsTemplateTools.h>

#include <gsUtils/gsNorms.h>
#include <gsUtils/gsNorms.hpp>

#define T real_t
#define uZ unsigned
#define Z int

namespace gismo
{

/*
 * Norm & distance computation
 */

TEMPLATE_INST
T computeL2Norm<T>(const gsGeometry<T>& geo, const gsFunction<T>& u, bool isParametrized_u, int numEvals);

TEMPLATE_INST
T computeL2Norm<T>(const gsField<T>& u, int numSamples);

TEMPLATE_INST
T computeL2Distance<T>(const gsGeometry<T>& geo, const gsFunction<T>& u, bool isParametrized_u, const gsFunction<T>& v, bool isParametrized_v, int numEvals);

TEMPLATE_INST
T computeL2Distance<T>(const gsField<T>& u, const gsFunction<T>& v, bool isParametrized_v, int numEvals);

TEMPLATE_INST
T computeL2Distance<T>(const gsField<T>& u, const gsField<T>& v, int numEvals);

TEMPLATE_INST
T igaFieldL2Distance(const gsField<T>& u,
                     const gsFunction<T>& v,
                     bool v_isParam);

TEMPLATE_INST
T igaFieldL2Distance(const gsField<T>& u,
                     const gsFunction<T>& v,
		     const gsMultiBasis<T>& B,
                     bool v_isParam);

TEMPLATE_INST
T igaL2Distance(const gsGeometry<T>& patch, 
                const gsGeometry<T>& func, 
                const gsFunction<T>& v, 
                bool v_isParam);

TEMPLATE_INST
T igaL2Distance(const gsGeometry<T>& patch, 
                const gsFunction<T>& func, 
                const gsFunction<T>& v, 
		const gsBasis<T>& B,
                bool v_isParam);

TEMPLATE_INST
gsVector< gsMatrix<T> > igaFieldL2DistanceEltWiseSq(const gsField<T>& u, const gsFunction<T>& v, bool v_isParam);

TEMPLATE_INST
gsMatrix<T> igaL2DistanceEltWiseSq(const gsGeometry<T>& patch,
                const gsGeometry<T>& func,
                const gsFunction<T>& v,
                bool v_isParam);

TEMPLATE_INST
T igaL2DistanceOnElt( const gsGeometryEvaluator<T>::uPtr & geoEval ,
                      const gsGeometryEvaluator<T>::uPtr & funcEval,
                      const gsFunction<T>& v,
                      const bool & v_isParam,
                      const gsBasis<T>::domainIter & domIt,
					  const gsQuadRule<T> & quRule);

TEMPLATE_INST
T igaH1Distance(const gsGeometry<T>& patch, 
                const gsGeometry<T>& func, 
                const gsFunction<T>& v, 
                bool v_isParam);

TEMPLATE_INST
T igaFieldH1Distance(const gsField<T>& u,
                     const gsFunction<T>& v,
                     const gsMultiBasis<T>& B,
                     bool v_isParam);
TEMPLATE_INST
T igaH1Distance(const gsGeometry<T>& patch, 
                const gsFunction<T>& func, 
                const gsFunction<T>& v, 
                const gsBasis<T>& B,
                bool v_isParam);
TEMPLATE_INST
T igaH1DistanceOnElt( const gsGeometryEvaluator<T>::uPtr & geoEval ,
                      const gsFunction<T> & func,
                      const gsFunction<T>& v,
                      const bool & v_isParam,
                      const gsBasis<T>::domainIter & domIt,
					  const gsQuadRule<T> & quRule);


TEMPLATE_INST
gsMatrix<T> igaH1DistanceEltWiseSq(const gsGeometry<T>& patch,
                const gsGeometry<T>& func,
                const gsFunction<T>& v,
                bool v_isParam );

TEMPLATE_INST
T igaFieldH1Distance(const gsField<T>& u, 
                     const gsFunction<T>& v, 
                     bool v_isParam);


TEMPLATE_INST
T igaH1DistanceOnElt( const gsGeometryEvaluator<T>::uPtr & geoEval ,
                      const gsGeometryEvaluator<T>::uPtr & funcEval,
                      const gsFunction<T>& v,
                      const bool & v_isParam,
                      const gsBasis<T>::domainIter & domIt,
                      const int d,
					  const gsQuadRule<T> & quRule);

TEMPLATE_INST
gsVector< gsMatrix<T> > igaFieldH1DistanceEltWiseSq(const gsField<T>& u, const gsFunction<T>& v, bool v_isParam);






TEMPLATE_INST
T computeMaximumNorm<T>(const gsFunction<T>& f, const gsVector<T>& lower, const gsVector<T>& upper, int numSamples);

TEMPLATE_INST
T computeMaximumNorm<T>(const gsGeometry<T>& geo, const gsFunction<T>& u, bool isParametrized_u, int numSamples);

TEMPLATE_INST
T computeMaximumNorm<T>(const gsField<T>& u, int numSamples);

TEMPLATE_INST
T computeMaximumDistance<T>(const gsFunction<T>& f1, const gsFunction<T>& f2, const gsVector<T>& lower, const gsVector<T>& upper, int numSamples);

TEMPLATE_INST
T computeMaximumDistance<T>(const gsGeometry<T>& geo, const gsFunction<T>& u, bool isParametrized_u, const gsFunction<T>& v, bool isParametrized_v, int numSamples);

TEMPLATE_INST
T computeMaximumDistance<T>(const gsField<T>& u, const gsFunction<T>& v, bool isParametrized_v, int numSamples);

TEMPLATE_INST
T igaDGDistanceJump(const gsGeometry<T>& patch1, const gsGeometry<T>& patch2, const gsGeometry<T>& func1,  const gsGeometry<T>& func2, 
                    const gsFunction<T>& v1, const gsFunction<T>& v2, const boundaryInterface & bi, const T mu, bool v_isParam);

TEMPLATE_INST
T igaFieldDGDistance(const gsField<T>& u, const gsFunction<T>& v, bool v_isParam);

} // namespace gismo

#undef T
#undef uZ
#undef Z
