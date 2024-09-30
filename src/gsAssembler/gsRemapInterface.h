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
#include <gsIO/gsOptionList.h>

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
    /// @param  opt                Options; see below:
    ///
    ///  | Option                  | Meaning
    ///  |------------------------ | --------------------------------------------------------------------------------------
    ///  | CheckAffine             | The number of interior points (per direction) in point grid used to
    ///  |                         | check if the mapping is affine. If set to i, the grid consits of i interior
    ///  |                         | points and the two boundary points per direction, so \f$ (i+2)^d \f$ points.
    ///  |                         | For \a alwaysAffine (=0) or \a neverAffine (-1), no checks are performed; default: 1
    ///  | NumSamplePoints         | Number of sample points for computing fitting curve (if not affine); default: 11
    ///  | IntervalsOfFittingCurve | Number of knot space in the fitting curve (if not affine); default: 5
    ///  | DegreeOfFittingCurve    | Spline degree of fitting curve (if not affine); default: 3
    ///  | EqualityTolerance       | Points are considered equal iff difference is smaller; default 1e-5
    ///  | NewtonTolerance         | Tolerance for Newton solvers (should be significantly smaller than EqualityTolerance);
    ///  |                         | default: 1e-10
    gsRemapInterface(const gsMultiPatch<T>   & mp,
                     const gsMultiBasis<T>   & mb,
                     const boundaryInterface & bi,
                     const gsOptionList      & opt = defaultOptions() );

    /// @Returns default options
    ///
    /// See constructor for their meaning
    static gsOptionList defaultOptions()
    {
        gsOptionList result;
        result.addInt ( "CheckAffine",             "The number of interior points (per direction) in point grid used to "
                                                   "check if the mapping is affine",                                      1     );
        result.addInt ( "NumSamplePoints",         "Number of sample points for computing fitting curve (if not affine)", 11    );
        result.addInt ( "IntervalsOfFittingCurve", "Number of knot space in the fitting curve (if not affine)",           5     );
        result.addInt ( "DegreeOfFittingCurve",    "Spline degree of fitting curve (if not affine)",                      3     );
        result.addReal( "EqualityTolerance",       "Points are considered equal iff difference is smaller",               1e-5  );
        result.addReal( "NewtonTolerance",         "Tolerance for Newton solvers",                                        1e-10 );
        return result;
    }

public:

    /// @brief Returns the interface map iff affine, otherwise the fitting curve
    ///
    /// Interfaces map \f$ \widehat \Gamma_1 \rightarrow \widehat \Gamma_2 \f$ that represents
    /// \f$ G_2^{-1} \circ G_1 \f$
    const typename gsFunction<T>::Ptr & interfaceMap() const { return m_intfMap; }

    /// @brief Evaluates the interface map
    ///
    /// Interfaces map \f$ \widehat \Gamma_1 \rightarrow \widehat \Gamma_2 \f$ that represents
    /// \f$ G_2^{-1} \circ G_1 \f$
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const;

    /// Returns parameter dimension of the domains
    virtual short_t domainDim() const { return m_g1->domainDim(); }

    /// @brief Returns a domain iterator
    ///
    /// The domain iterator lives on \f$ \widehat \Gamma_1 \f$. Its break points are the union of
    /// the brakpoints of the basis on \f$ \widehat \Omega_1 \f$ and the breakpoints of the basis
    /// on \f$ \widehat \Omega_2 \f$, mapped to \f$ \widehat \Omega_1 \f$.
    typename gsDomainIterator<T>::uPtr makeDomainIterator() const;

    /// Returns true iff the interface is matching
    bool isMatching() const { return m_isMatching; }

    /// Returns true iff the interface is affine
    bool isAffine() const { return m_isAffine; }

    /// Returns the break points used in \ref makeDomainIterator
    const std::vector< std::vector<T> > & breakPoints() const { return m_breakpoints; }

    /// Prints the state of the object
    virtual std::ostream & print(std::ostream& os) const;

private:

    /// Computes the box which represents the intersection of sides of incoming patches
    void constructInterfaceBox();

    /// Estimates error between geometry1 and geometry2 after repearametrization on physical domain
    T estimateReparamError(index_t steps) const;

    /// Constructs the fitting curve \a m_intfMap in the non-affine case
    void constructFittingCurve(index_t numSamplePoints, index_t intervalsOfFittingCurve,
            index_t degreeOfFittingCurve);

    /// Constructs the breakpoints \a m_breakpoints
    void constructBreaks();

private:
    const gsGeometry<T> * m_g1;                       ///< Geometry of first patch
    const gsGeometry<T> * m_g2;                       ///< Geometry of second patch

    const gsBasis<T> * m_b1;                          ///< Basis on first patch
    const gsBasis<T> * m_b2;                          ///< Basis on second patch

    boundaryInterface m_bi;                           ///< Corresponding boundary interface

    bool m_isMatching;                                ///< True iff the interface is matching
    bool m_isAffine;                                  ///< True iff the interface is affine

    std::vector< std::vector<T> > m_breakpoints;      ///< Union of breakpoints of both bases

    typename gsFunction<T>::Ptr m_intfMap;            ///< Iff affine, interface map, otherwise the fitting curve

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
inline std::ostream & operator<<(std::ostream & os, const gsRemapInterface<T> & remapIf)
{ return remapIf.print(os); }

} // End namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRemapInterface.hpp)
#endif
