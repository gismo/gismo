/** @file gsRemapInterface.h

    @brief Provides a mapping between the corresponding sides of two patches sharing an interface

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Seiler, R. Schneckenleitner, S. Takacs
    Created on: 2018-06-12
*/


#pragma once

#include <gsCore/gsAffineFunction.h>
#include <gsCore/gsBoundary.h>

namespace gismo {

/** @brief Provides a mapping between the corresponding sides of two patches sharing an interface
  *
  * Have two patches \f$ \Omega_1 = G_1( \widehat \Omega_1 )\f$ and \f$ \Omega_2 = G_2( \widehat \Omega_2 )\f$
  * with geometry functions \f$ G_1 \f$ and \f$ G_2 \f$ and parameter domains \f$ \widehat \Omega_1 \f$ and
  * \f$ \widehat \Omega_2 \f$.
  *
  * This class represents an interface \f[ \Gamma := G_1( \widehat S_1 ) \cap G_2( \widehat S_2 ), \f] where
  * \f$ \widehat S_1 \f$ and \f$ \widehat S_2 \f$ are sides of the corresponding parameter domains (i.e., whole
  * edges or whole faces). The indices of the patches, as well as the corresponding sides (as \a boxSide)
  * are prescribed by the incoming \a boundaryInterface.
  *
  * The pre-images of the interface \f$ \Gamma \f$ are \f$ \widehat \Gamma_1 := G_1^{-1}( \Gamma ) \f$
  * and \f$ \widehat \Gamma_2 := G_2^{-1}( \Gamma ) \f$.
  *
  * We say the interface is \em matching if \f[ \Gamma = G_1( \widehat S_1 ) = G_2( \widehat S_2 ). \f]
  * In this case, the pre-images of the interface \f$ \Gamma \f$ are the whole sides, i.e.,
  * \f$ \widehat \Gamma_1 = \widehat S_1 \f$ and \f$ \widehat \Gamma_2 = \widehat S_2 \f$. Otherwise,
  * the pre-images are only subsets of the sides.
  *
  * We say the interface is \em affine if \f$ G_1^{-1} \circ G_2 \f$ or, euqivalently \f$ G_2^{-1} \circ G_1 \f$
  * is an affine-linear function. In this case, this class uses that affine-linear function.
  *
  * If the mapping is found out not to be affine-linear, this class constructs an approximation of the interface
  * map. This is available only if \f$ \Omega_1 \f$ and \f$ \Omega_2 \f$ are two dimensional, i.e., the interface
  * itself is one dimensional. We refer to this as the 2D case.
  *
  * @ingroup Assembler
 **/
template <class T>
class gsRemapInterface : public gsFunction<T>
{
public:

    /// Shared pointer for gsRemapInterface
    typedef memory::shared_ptr< gsRemapInterface > Ptr;

    /// Unique pointer for gsRemapInterface
    typedef memory::unique_ptr< gsRemapInterface > uPtr;

public:

    enum {
        neverAffine = -1,  ///< Instructs constructor not to use an affine mapping
        alwaysAffine = 0   ///< Instructs constructor always to use an affine mapping
    };

    /// @brief Constructor
    ///
    /// @param  mp                 The multi-patch object
    /// @param  mb                 The multi-basis object
    /// @param  bi                 The boundary interface (specifying the patches and their \a patchSide s)
    /// @param  checkAffine        The number of interior points (per direction) in point grid used to
    ///                            check if the mapping is affine. If set to i, the grid consits of i interor
    ///                            point and the two boundary points per direction, so \f$ (i+2)^d \f$ points.
    ///                            For \a alwaysAffine or \a notAffine, no checks are performed.
    /// @param  equalityTolerance  Points are considered equal iff difference is smaller.
    /// @param  newtonTolerance    Tolerance for Newton solvers (should be significantly smaller than
    ///                            equalityTolerance)
    gsRemapInterface(const gsMultiPatch<T> & mp,
                     const gsMultiBasis<T> & mb,
                     const boundaryInterface & bi,
                     index_t checkAffine = 1,
                     T equalityTolerance = 1e-5,
                     T newtonTolerance = 1e-8);

public:

    /// @brief Returns the interface map
    ///
    /// Interfaces map \f$ \widehat \Gamma_1 \rightarrow \widehat \Gamma_2 \f$ that represents
    /// \f$ G_2^{-1} \circ G_1 \f$
    const typename gsFunction<T>::Ptr & interfaceMap() const { return m_intfMap; }

    /// @brief Evaluates the interface map
    ///
    /// Interfaces map \f$ \widehat \Gamma_1 \rightarrow \widehat \Gamma_2 \f$ that represents
    /// \f$ G_2^{-1} \circ G_1 \f$
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// @brief Returns parameter dimension of the domains
    virtual short_t domainDim() const { return m_g1->geoDim(); }

    /// @brief Returns a domain iterator
    ///
    /// The domain iterator lives on \f$ \widehat \Gamma_1 \f$. Its break points are the union of
    /// the brakpoints of the basis on \f$ \widehat \Omega_1 \f$ and the breakpoints of the basis
    /// on \f$ \widehat \Omega_2 \f$, mapped to \f$ \widehat \Omega_1 \f$.
    typename gsDomainIterator<T>::uPtr makeDomainIterator() const;

    /// @brief Returns true iff the interface is matching
    bool isMatching() const { return m_isMatching; }

    /// @brief Returns true iff the interface is affine
    bool isAffine() const { return m_isAffine; }

    /// @brief Returns the break points used in \ref makeDomainIterator
    const std::vector< std::vector<T> > & breakPoints() const { return m_breakpoints; }

    /// @brief Prints the state of the object
    virtual std::ostream & print(std::ostream& os) const override;

private:

    /// Computes the box which represents the intersection of sides of incoming patches
    void computeBoundingBox();

    /// Checks if affine mapping between the incoming patches is correct
    bool checkIfAffine(index_t steps);

    /// Constructs the breakpoints \a m_breakpoints
    void constructBreaks();

    /// Helper to compute the closest point to lti on the other patch via Newton's method
    gsMatrix<T> closestPoint(const gsMatrix<T> b_null, const gsGeometry<T> & R, const gsMatrix<T> & lti);

    /// Member to enrich a matrix of 1D points to a matrix of m_domain.geoDim() points
    void enrichToVector(boxSide boundarySide, const gsGeometry<T> & geo, const gsMatrix<T> & intervals, gsMatrix<T> & pts);

    /// Check if the incoming evaluation points are out of bounds because of rounding errors
    gsMatrix<T> checkIfInBound(const gsMatrix<T> & u) const;

    /// Constructs the reparametrization \a m_intfMap in the non-affine case
    void constructReparam();

private:
    const gsGeometry<T> * m_g1;                       ///< Geometry of first patch
    const gsGeometry<T> * m_g2;                       ///< Geometry of second patch

    const gsBasis<T> * m_b1;                          ///< Basis on first patch
    const gsBasis<T> * m_b2;                          ///< Basis on second patch

    boundaryInterface m_bi;                           ///< Corresponding boundary interface

    bool m_isMatching;                                ///< True iff the interface is matching
    bool m_isAffine;                                  ///< True iff the interface is affine

    std::vector< std::vector<T> > m_breakpoints;      ///< Union of breakpoints of both bases

    // TODO: This is only the interface map in the affine case, otherwise it is the fitting curve.
    typename gsFunction<T>::Ptr m_intfMap;            ///< The interface map itself

    /// @brief The bounds of the box that represents \f$ \widehat \Gamma_1 \f$
    ///
    /// If \f$ \widehat \Gamma_1 = (a,b)\times (c,d) \times (e,f) \f$, the matrix
    /// has the structure
    /// \f[ \begin{pmatrix} a & c & e \\ b & d & f \end{pmatrix} \f]
    gsMatrix<T> m_parameterBounds1;

    /// @brief The bounds of the box that represents \f$ \widehat \Gamma_2 \f$
    ///
    /// See \ref m_parameterBounds1
    gsMatrix<T> m_parameterBounds2;

    T m_equalityTolerance;                            ///< Tolerance for considering points to be equal
    T m_newtonTolerance;                              ///< Tolerance for Newton solver

}; // End gsRemapInterface

/// @brief   Prints the state of the object
/// @relates gsRemapInterface
template <class T>
inline std::ostream & operator<<(std::ostream& os, const gsRemapInterface<T>& remapIf)
{ return remapIf.print(os); }

} // End namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRemapInterface.hpp)
#endif
