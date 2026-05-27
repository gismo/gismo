#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsDimMacro.h>

#include <gsNurbs/gsBoehm.h>
#include <gsNurbs/gsBoehm.hpp>
#include <gsNurbs/gsKnotVector.h>

#include <gsMatrix/gsFiberMatrix.h>

#define T real_t

namespace gismo
{

// gsBoehm

TEMPLATE_INST
void gsBoehm<T, gsKnotVector<T>, gsMatrix<T> >(
    gsKnotVector<T> & knots,
    gsMatrix<T> & coefs,
    T val,
    int r,
    bool update_knots
    );

TEMPLATE_INST
void gsBoehmSingle<T, gsKnotVector<T>, gsMatrix<T> >(
    gsKnotVector<T> & knots,
    gsMatrix<T> & coefs,
    T val,
    bool update_knots
    );

// gsBoehmSingle (v2)

TEMPLATE_INST
void gsBoehmSingle<T, gsKnotVector<T>::iterator, gsMatrix<T> >(
    gsKnotVector<T>::iterator knot,
    gsMatrix<T> & coefs,
    int p,
    T val
    );

// gsBoehmRefine gsKnotVector + gsMatrix + iterator / const_iterator
//               gsKnotVector + gsFiberMatrix + const_iterator

TEMPLATE_INST
void gsBoehmRefine<gsKnotVector<T>,
                   gsMatrix<T>,
                   std::vector<T>::const_iterator>(
    gsKnotVector<T> & knots,
    gsMatrix<T> & coefs,
    int p,
    std::vector<T>::const_iterator valBegin,
    std::vector<T>::const_iterator valEnd,
    bool update_knots
    );

TEMPLATE_INST
void gsBoehmRefine<gsKnotVector<T>,
                   gsMatrix<T>,
                   std::vector<T>::iterator>(
    gsKnotVector<T> & knots,
    gsMatrix<T> & coefs,
    int p,
    std::vector<T>::iterator valBegin,
    std::vector<T>::iterator valEnd,
    bool update_knots
    );

TEMPLATE_INST
void gsBoehmRefine<gsKnotVector<T>,
                   gsFiberMatrix<T,RowMajor>,
                   std::vector<T>::const_iterator>(
    gsKnotVector<T> & knots,
    gsFiberMatrix<T,RowMajor> & coefs,
    int p,
    std::vector<T>::const_iterator valBegin,
    std::vector<T>::const_iterator valEnd,
    bool update_knots
    );

// gsTensorBoehm

TEMPLATE_INST
void gsTensorBoehm<T, gsKnotVector<T>, gsMatrix<T> >(
        gsKnotVector<T>& knots,
        gsMatrix<T>& coefs,
        T val,
        int direction,
        gsVector<unsigned> str,
        int r,
        bool update_knots);

// gsTensorBoehmRefine

TEMPLATE_INST
void gsTensorBoehmRefine<gsKnotVector<T>,
                         gsMatrix<T>,
                         std::vector<T>::const_iterator>(
        gsKnotVector<T>& knots,
        gsMatrix<T>& coefs,
        const int direction,
        gsVector<unsigned> str,
        std::vector<T>::const_iterator valBegin,
        std::vector<T>::const_iterator valEnd,
        bool update_knots);

// gsTensorBoehmRefineLocal

#define INST_BOEHM_REFINE_LOCAL(D) \
TEMPLATE_INST \
void gsTensorBoehmRefineLocal<D, \
                              gsKnotVector<T>, \
                              gsMatrix<T>, \
                              std::vector<T>::const_iterator>( \
        gsKnotVector<T>& knots, \
        const unsigned index, \
        gsMatrix<T>& coefs, \
        gsVector<index_t, D>& nmb_of_coefs, \
        const gsVector<index_t, D>& act_size_of_coeffs, \
        const gsVector<index_t, D>& size_of_coefs, \
        const unsigned direction, \
        std::vector<T>::const_iterator valBegin, \
        std::vector<T>::const_iterator valEnd, \
        const bool update_knots);

GISMO_DIM_FOREACH(INST_BOEHM_REFINE_LOCAL)
#undef INST_BOEHM_REFINE_LOCAL

// =============================================================================
// gsTensorInsertKnotDegreeTimes
// =============================================================================

#define INST_INSERT_KNOT(D) \
TEMPLATE_INST \
void gsTensorInsertKnotDegreeTimes<D, \
                        T, \
                        gsKnotVector<T>, \
                        gsMatrix<T> >( \
        const gsKnotVector<T>& knots, \
        gsMatrix<T>& coefs, \
        const gsVector<index_t, D>& size_of_coefs, \
        T val, \
        const unsigned direction, \
        gsVector<index_t, D>& start, \
        gsVector<index_t, D>& end);

GISMO_DIM_FOREACH_FROM2(INST_INSERT_KNOT)
#undef INST_INSERT_KNOT


} // end namespace gismo

#undef T
