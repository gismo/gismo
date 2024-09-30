#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsBoehm.h>
#include <gsNurbs/gsBoehm.hpp>
#include <gsNurbs/gsKnotVector.h>

#include <gsMatrix/gsSparseRows.hpp>

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
//               gsKnotVector + gsSparseRows + const_iterator

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
                   gsSparseRows<T>,
                   std::vector<T>::const_iterator>(
    gsKnotVector<T> & knots,
    gsSparseRows<T> & coefs,
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

TEMPLATE_INST
void gsTensorBoehmRefineLocal<1,
                              gsKnotVector<T>,
                              gsMatrix<T>,
                              std::vector<T>::const_iterator>(
        gsKnotVector<T>& knots,
        const unsigned index,
        gsMatrix<T>& coefs,
        gsVector<index_t, 1>& nmb_of_coefs,
        const gsVector<index_t, 1>& act_size_of_coeffs,
        const gsVector<index_t, 1>& size_of_coefs,
        const unsigned direction,
        std::vector<T>::const_iterator valBegin,
        std::vector<T>::const_iterator valEnd,
        //const unsigned number_of_iterations,
        const bool update_knots);


TEMPLATE_INST
void gsTensorBoehmRefineLocal<2,
                              gsKnotVector<T>,
                              gsMatrix<T>,
                              std::vector<T>::const_iterator>(
        gsKnotVector<T>& knots,
        const unsigned index,
        gsMatrix<T>& coefs,
        gsVector<index_t, 2>& nmb_of_coefs,
        const gsVector<index_t, 2>& act_size_of_coeffs,
        const gsVector<index_t, 2>& size_of_coefs,
        const unsigned direction,
        std::vector<T>::const_iterator valBegin,
        std::vector<T>::const_iterator valEnd,
        //const unsigned number_of_iterations,
        const bool update_knots);

TEMPLATE_INST
void gsTensorBoehmRefineLocal<3,
                              gsKnotVector<T>,
                              gsMatrix<T>,
                              std::vector<T>::const_iterator>(
        gsKnotVector<T>& knots,
        const unsigned index,
        gsMatrix<T>& coefs,
        gsVector<index_t, 3>& nmb_of_coefs,
        const gsVector<index_t, 3>& act_size_of_coeffs,
        const gsVector<index_t, 3>& size_of_coefs,
        const unsigned direction,
        std::vector<T>::const_iterator valBegin,
        std::vector<T>::const_iterator valEnd,
        //const unsigned number_of_iterations,
        const bool update_knots);

TEMPLATE_INST
void gsTensorBoehmRefineLocal<4,
                              gsKnotVector<T>,
                              gsMatrix<T>,
                              std::vector<T>::const_iterator>(
        gsKnotVector<T>& knots,
        const unsigned index,
        gsMatrix<T>& coefs,
        gsVector<index_t, 4>& nmb_of_coefs,
        const gsVector<index_t, 4>& act_size_of_coeffs,
        const gsVector<index_t, 4>& size_of_coefs,
        const unsigned direction,
        std::vector<T>::const_iterator valBegin,
        std::vector<T>::const_iterator valEnd,
        //const unsigned number_of_iterations,
        const bool update_knots);

// =============================================================================
// gsTensorInsertKnotDegreeTimes
// =============================================================================

TEMPLATE_INST
void gsTensorInsertKnotDegreeTimes<3,
                        T,
                        gsKnotVector<T>,
                        gsMatrix<T> >(
        const gsKnotVector<T>& knots,
        gsMatrix<T>& coefs,
        const gsVector<index_t, 3>& size_of_coefs,
        T val,
        const unsigned direction,
        gsVector<index_t, 3>& start,
        gsVector<index_t, 3>& end);

TEMPLATE_INST
void gsTensorInsertKnotDegreeTimes<2,
                        T,
                        gsKnotVector<T>,
                        gsMatrix<T> >(
        const gsKnotVector<T>& knots,
        gsMatrix<T>& coefs,
        const gsVector<index_t, 2>& size_of_coefs,
        T val,
        const unsigned direction,
        gsVector<index_t, 2>& start,
        gsVector<index_t, 2>& end);


} // end namespace gismo

#undef T
