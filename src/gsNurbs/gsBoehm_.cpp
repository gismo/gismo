#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsBoehm.h>
#include <gsNurbs/gsBoehm.hpp>
#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsCompactKnotVector.h>

#include <gsMatrix/gsSparseRows.hpp>

#define T real_t

namespace gismo
{

// gsBoehm: gsKnotVector + gsCompactKnotVector

TEMPLATE_INST
void gsBoehm<T, gsKnotVector<T>, gsMatrix<T> >(
    gsKnotVector<T> & knots,
    gsMatrix<T> & coefs,
    T val,
    int r,
    bool update_knots
    );

/*
TEMPLATE_INST
void gsBoehm<T, gsCompactKnotVector<T>, gsMatrix<T> >(
    gsCompactKnotVector<T> & knots,
    gsMatrix<T> & coefs,
    T val,
    int r,
    bool update_knots
    );

*/

// gsBoehmSingle gsKnotVector + gsCompactKnotVector

TEMPLATE_INST
void gsBoehmSingle<T, gsKnotVector<T>, gsMatrix<T> >(
    gsKnotVector<T> & knots,
    gsMatrix<T> & coefs,
    T val,
    bool update_knots
    );

/*
TEMPLATE_INST
void gsBoehmSingle<T, gsCompactKnotVector<T>, gsMatrix<T> >(
    gsCompactKnotVector<T> & knots,
    gsMatrix<T> & coefs,
    T val,
    bool update_knots
    );
*/


// gsBoehmSingle (v2)

TEMPLATE_INST
void gsBoehmSingle<T, gsKnotVector<T>::iterator, gsMatrix<T> >(
    gsKnotVector<T>::iterator knot,
    gsMatrix<T> & coefs,
    int p,
    T val
    );



// gsBoehmRefine gsKnotVector + gsMatrix + iterator / const_iterator
//               gsCompactKnotVector + gsMatrix + iterator / const_iterator
//               gsCompactKnotVector + gsSparseRows + const_iterator
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

/*
TEMPLATE_INST
void gsBoehmRefine<gsCompactKnotVector<T>,
                   gsMatrix<T>,
                   std::vector<T>::iterator>(
    gsCompactKnotVector<T> & knots,
    gsMatrix<T> & coefs,
    int p,
    std::vector<T>::iterator valBegin,
    std::vector<T>::iterator valEnd,
    bool update_knots
    );

TEMPLATE_INST
void gsBoehmRefine<gsCompactKnotVector<T>,
                   gsMatrix<T>,
                   std::vector<T>::const_iterator>
(
    gsCompactKnotVector<T> & knots,
    gsMatrix<T> & coefs,
    int p,
    std::vector<T>::const_iterator valBegin,
    std::vector<T>::const_iterator valEnd,
    bool update_knots
    );

TEMPLATE_INST
void gsBoehmRefine<gsCompactKnotVector<T>,
                   gsSparseRows<T>,
                   std::vector<T>::const_iterator>(
    gsCompactKnotVector<T> & knots,
    gsSparseRows<T> & coefs,
    int p,
    std::vector<T>::const_iterator valBegin,
    std::vector<T>::const_iterator valEnd,
    bool update_knots
    );
*/

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

/*
TEMPLATE_INST
void gsTensorBoehm<T, gsCompactKnotVector<T>, gsMatrix<T> >(
        gsCompactKnotVector<T>& knots,
        gsMatrix<T>& coefs,
        T val,
        int direction,
        gsVector<unsigned> str,
        int r,
        bool update_knots);
*/

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


/*
TEMPLATE_INST
void gsTensorBoehmRefine<gsCompactKnotVector<T>,
                         gsMatrix<T>,
                         std::vector<T>::const_iterator>(
        gsCompactKnotVector<T>& knots,
        gsMatrix<T>& coefs,
        const int direction,
        gsVector<unsigned> str,
        std::vector<T>::const_iterator valBegin,
        std::vector<T>::const_iterator valEnd,
        bool update_knots);
*/


// gsTensorBoehmRefineLocal

TEMPLATE_INST
void gsTensorBoehmRefineLocal<1,
                              gsCompactKnotVector<T>,
                              gsMatrix<T>,
                              std::vector<T>::const_iterator>(
        gsCompactKnotVector<T>& knots,
        const unsigned index,
        gsMatrix<T>& coefs,
        gsVector<unsigned, 1>& nmb_of_coefs,
        const gsVector<unsigned, 1>& act_size_of_coeffs,
        const gsVector<unsigned, 1>& size_of_coefs,
        const unsigned direction,
        std::vector<T>::const_iterator valBegin,
        std::vector<T>::const_iterator valEnd,
        //const unsigned number_of_iterations,
        const bool update_knots);


TEMPLATE_INST
void gsTensorBoehmRefineLocal<2,
                              gsCompactKnotVector<T>,
                              gsMatrix<T>,
                              std::vector<T>::const_iterator>(
        gsCompactKnotVector<T>& knots,
        const unsigned index,
        gsMatrix<T>& coefs,
        gsVector<unsigned, 2>& nmb_of_coefs,
        const gsVector<unsigned, 2>& act_size_of_coeffs,
        const gsVector<unsigned, 2>& size_of_coefs,
        const unsigned direction,
        std::vector<T>::const_iterator valBegin,
        std::vector<T>::const_iterator valEnd,
        //const unsigned number_of_iterations,
        const bool update_knots);

TEMPLATE_INST
void gsTensorBoehmRefineLocal<3,
                              gsCompactKnotVector<T>,
                              gsMatrix<T>,
                              std::vector<T>::const_iterator>(
        gsCompactKnotVector<T>& knots,
        const unsigned index,
        gsMatrix<T>& coefs,
        gsVector<unsigned, 3>& nmb_of_coefs,
        const gsVector<unsigned, 3>& act_size_of_coeffs,
        const gsVector<unsigned, 3>& size_of_coefs,
        const unsigned direction,
        std::vector<T>::const_iterator valBegin,
        std::vector<T>::const_iterator valEnd,
        //const unsigned number_of_iterations,
        const bool update_knots);

TEMPLATE_INST
void gsTensorBoehmRefineLocal<4,
                              gsCompactKnotVector<T>,
                              gsMatrix<T>,
                              std::vector<T>::const_iterator>(
        gsCompactKnotVector<T>& knots,
        const unsigned index,
        gsMatrix<T>& coefs,
        gsVector<unsigned, 4>& nmb_of_coefs,
        const gsVector<unsigned, 4>& act_size_of_coeffs,
        const gsVector<unsigned, 4>& size_of_coefs,
        const unsigned direction,
        std::vector<T>::const_iterator valBegin,
        std::vector<T>::const_iterator valEnd,
        //const unsigned number_of_iterations,
        const bool update_knots);

// =============================================================================
// gsTensorInsertKnotDegreeTimes
// =============================================================================

/*
TEMPLATE_INST
void gsTensorInsertKnotDegreeTimes<3,
                        T,
                        gsCompactKnotVector<T>,
                        gsMatrix<T> >(
        const gsCompactKnotVector<T>& knots,
        gsMatrix<T>& coefs,
        const gsVector<unsigned, 3>& size_of_coefs,
        T val,
        const unsigned direction,
        gsVector<unsigned, 3>& start,
        gsVector<unsigned, 3>& end);
*/

TEMPLATE_INST
void gsTensorInsertKnotDegreeTimes<3,
                        T,
                        gsKnotVector<T>,
                        gsMatrix<T> >(
        const gsKnotVector<T>& knots,
        gsMatrix<T>& coefs,
        const gsVector<unsigned, 3>& size_of_coefs,
        T val,
        const unsigned direction,
        gsVector<unsigned, 3>& start,
        gsVector<unsigned, 3>& end);

/*
TEMPLATE_INST
void gsTensorInsertKnotDegreeTimes<2,
                        T,
                        gsCompactKnotVector<T>,
                        gsMatrix<T> >(
        const gsCompactKnotVector<T>& knots,
        gsMatrix<T>& coefs,
        const gsVector<unsigned, 2>& size_of_coefs,
        T val,
        const unsigned direction,
        gsVector<unsigned, 2>& start,
        gsVector<unsigned, 2>& end);
*/

TEMPLATE_INST
void gsTensorInsertKnotDegreeTimes<2,
                        T,
                        gsKnotVector<T>,
                        gsMatrix<T> >(
        const gsKnotVector<T>& knots,
        gsMatrix<T>& coefs,
        const gsVector<unsigned, 2>& size_of_coefs,
        T val,
        const unsigned direction,
        gsVector<unsigned, 2>& start,
        gsVector<unsigned, 2>& end);


} // end namespace gismo

#undef T
