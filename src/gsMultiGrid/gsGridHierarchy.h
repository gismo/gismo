/** @file gsGridHierarchy.h

    @brief Coarsening algorithms for knot vectors and bases.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#pragma once

#include <gsCore/gsDofMapper.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

/** @brief
    Grid Hierarchy

    This class allows to construct a grid hierarchy and stors a grid hierarchy (vector of
    bases, local transfer matrices and transfer matrices).

    \ingroup Solver
*/
template< typename T >
class gsGridHierarchy
{

public:

    /// @brief This function sets up a multigrid hierarchy by uniform refinement
    ///
    /// @param mBasis                    The gsMultiBasis to be refined (initial basis)
    /// @param boundaryConditions        The boundary conditions
    /// @param assemblerOptions          A gsOptionList defining a "DirichletStrategy" and a "InterfaceStrategy"
    /// @param levels                    The number of levels
    /// @param numberOfKnotsToBeInserted The number of knots to be inserted, defaulted to 1
    /// @param multiplicityOfKnotsToBeInserted The multiplicity of the knots to be inserted, defaulted to 1
    ///
    /// \ingroup Solver
    static gsGridHierarchy buildByRefinement(
        gsMultiBasis<T> mBasis,
        const gsBoundaryConditions<T>& boundaryConditions,
        const gsOptionList& assemblerOptions,
        index_t levels,
        index_t numberOfKnotsToBeInserted = 1,
        index_t multiplicityOfKnotsToBeInserted = 1
        );

    /// @brief This function sets up a multigrid hierarchy by uniform refinement
    ///
    /// @param mBasis                    The gsMultiBasis to be refined (initial basis)
    /// @param boundaryConditions        The boundary conditions
    /// @param options                   A gsOptionList defining the necessary infomation
    ///
    /// \ingroup Solver
    static gsGridHierarchy buildByRefinement(
        gsMultiBasis<T> mBasis,
        const gsBoundaryConditions<T>& boundaryConditions,
        const gsOptionList& options
        )
    {
        return gsGridHierarchy::buildByRefinement(
            give(mBasis),
            boundaryConditions,
            options,
            options.askInt( "Levels", 3 ),
            options.askInt( "NumberOfKnotsToBeInserted", 1 ),
            options.askInt( "MultiplicityOfKnotsToBeInserted", 1 )
        );
    }


    /// @brief This function sets up a grid hierarchy by coarsening
    ///
    /// @param mBasis                    The gsMultiBasis to be coarsened (initial basis)
    /// @param boundaryConditions        The boundary conditions
    /// @param assemblerOptions          A gsOptionList defining a "DirichletStrategy" and a "InterfaceStrategy"
    /// @param levels                    The maximum number of levels
    /// @param degreesOfFreedom          Number of dofs in the coarsest grid in the grid hierarchy
    ///
    /// The algorithm terminates if either the number of levels is reached or the number of degrees of freedom
    /// is below the given threshold.
    ///
    static gsGridHierarchy buildByCoarsening(
        gsMultiBasis<T> mBasis,
        const gsBoundaryConditions<T>& boundaryConditions,
        const gsOptionList& assemblerOptions,
        index_t levels,
        index_t degreesOfFreedom = 0
        );

    /// @brief This function sets up a grid hierarchy by coarsening
    ///
    /// @param mBasis                    The gsMultiBasis to be coarsened (initial basis)
    /// @param boundaryConditions        The boundary conditions
    /// @param options                   A gsOptionList defining the necessary infomation
    ///
    static gsGridHierarchy buildByCoarsening(
        gsMultiBasis<T> mBasis,
        const gsBoundaryConditions<T>& boundaryConditions,
        const gsOptionList& options
    )
    {
        return gsGridHierarchy::buildByCoarsening(
            give(mBasis),
            boundaryConditions,
            options,
            options.askInt( "Levels", 3 ),
            options.askInt( "DegreesOfFreedom", 0 )
        );
    }

    /// Get the default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addInt( "DirichletStrategy", "Method for enforcement of Dirichlet BCs [11..14]", 11 );
        opt.addInt( "InterfaceStrategy", "Method of treatment of patch interfaces [0..3]", 1  );
        opt.addInt( "Levels", "Number of levels to be constructed in the grid hierarchy", 3 );
        opt.addInt( "DegreesOfFreedom",   "Number of dofs in the coarsest grid in the grid hierarchy (only buildByCoarsening)", 0 );
        opt.addInt( "NumberOfKnotsToBeInserted", "The number of knots to be inserted (only buildByRefinement)", 1 );
        opt.addInt( "MultiplicityOfKnotsToBeInserted",   "The multiplicity of the knots to be inserted (only buildByRefinement)", 1 );
        return opt;
    }

    /// Reset the object (to save memory)
    void clear()
    {
        //m_boundaryConditions.clear();
        //m_assemblerOptions.clear();
        m_mBases.clear();
        m_transferMatrices.clear();
        m_localTransferMatrices.clear();
    }

    /// Get the vector of multi bases (by reference)
    const std::vector< gsMultiBasis<T> >& getMultiBases() const
    { return m_mBases; }
    /// Get the vector of multi bases
    gsGridHierarchy& moveMultiBasesTo( std::vector< gsMultiBasis<T> >& o )
    { o = give(m_mBases); return *this; }

    /// Get the vector of transfer matrices (by reference)
    const std::vector< gsSparseMatrix<T, RowMajor> >& getTransferMatrices() const
    { return m_transferMatrices; }
    /// Get the vector of transfer matrices
    gsGridHierarchy& moveTransferMatricesTo( std::vector< gsSparseMatrix<T, RowMajor> >& o )
    { o = give(m_transferMatrices); return *this; }

    /// Get the vector of local transfer matrices (by reference)
    const std::vector< std::vector< gsSparseMatrix<T, RowMajor> > >& getLocalTransferMatrices() const
    { return m_localTransferMatrices; }
    /// Get the vector of local transfer matrices
    gsGridHierarchy& moveLocalTransferMatricesTo( std::vector< std::vector< gsSparseMatrix<T, RowMajor> > >& o )
    { o = give(m_localTransferMatrices); return *this; }

private:
    gsBoundaryConditions<T> m_boundaryConditions;
    gsOptionList m_assemblerOptions;

    std::vector< gsMultiBasis<T> > m_mBases;
    std::vector< gsSparseMatrix<T, RowMajor> > m_transferMatrices;
    std::vector< std::vector< gsSparseMatrix<T, RowMajor> > > m_localTransferMatrices;
};

/// @brief This function refines a gsMultiBasis uniformly and provides the transfer matrix.
///
/// @param mBasis[in]                          The gsMultiBasis to be refined
/// @param boundaryConditions[in]              The boundary conditions
/// @param assemblerOptions[in]                A gsOptionList defining a "DirichletStrategy" and a "InterfaceStrategy"
/// @param numberOfKnotsToBeInserted[in]       The number of knots to be inserted (typically 1), cf. the corresponding parameter in gsBasis
/// @param multiplicityOfKnotsToBeInserted[in] The multiplicity of the knots to be inserted (typically 1), cf. the corresponding parameter in gsBasis
/// @param refinedMBasis[out]                  The refined gsMultiBasis
/// @param transferMatrix[out]                 The transfer matrix for the free dofs
/// @param localTransferMatrices[out]          A vector of the local transfer matrices
///
/// \ingroup Solver
template <typename T>
void uniformRefine_withTransfer(
    const gsMultiBasis<T>& mBasis,
    const gsBoundaryConditions<T>& boundaryConditions,
    const gsOptionList& assemblerOptions,
    index_t numberOfKnotsToBeInserted,
    index_t multiplicityOfKnotsToBeInserted,
    gsMultiBasis<T> & refinedMBasis,
    gsSparseMatrix<T, RowMajor>& transferMatrix,
    std::vector< gsSparseMatrix<T, RowMajor> >& localTransferMatrices
    );

/// @brief This function refines a gsMultiBasis uniformly and provides the transfer matrix.
///
/// @param mBasis[in]                          The gsMultiBasis to be refined
/// @param boundaryConditions[in]              The boundary conditions
/// @param assemblerOptions[in]                A gsOptionList defining a "DirichletStrategy" and a "InterfaceStrategy"
/// @param numberOfKnotsToBeInserted[in]       The number of knots to be inserted (typically 1), cf. the corresponding parameter in gsBasis
/// @param multiplicityOfKnotsToBeInserted[in] The multiplicity of the knots to be inserted (typically 1), cf. the corresponding parameter in gsBasis
/// @param refinedMBasis[out]                  The refined gsMultiBasis
/// @param transferMatrix[out]                 The transfer matrix for the free dofs
///
/// \ingroup Solver
template <typename T>
inline void uniformRefine_withTransfer(
    const gsMultiBasis<T>& mBasis,
    const gsBoundaryConditions<T>& boundaryConditions,
    const gsOptionList& assemblerOptions,
    index_t numberOfKnotsToBeInserted,
    index_t multiplicityOfKnotsToBeInserted,
    gsMultiBasis<T> & refinedMBasis,
    gsSparseMatrix<T, RowMajor>& transferMatrix
    )
{
    std::vector< gsSparseMatrix<T, RowMajor> > localTransferMatrices;
    uniformRefine_withTransfer(
        mBasis, boundaryConditions, assemblerOptions, numberOfKnotsToBeInserted, multiplicityOfKnotsToBeInserted,
        refinedMBasis, transferMatrix, localTransferMatrices
    );
}

/// @brief This function refines a gsMultiBasis uniformly and provides the transfer matrix.
/// with \a refinedKnots = 1 and \a mult = 1
///
/// @param mBasis[in]                    The gsMultiBasis to be refined
/// @param boundaryConditions[in]              The boundary conditions
/// @param assemblerOptions[in]                A gsOptionList defining a "DirichletStrategy" and a "InterfaceStrategy"
/// @param refinedMBasis[out]                  The refined gsMultiBasis
/// @param transferMatrix[out]                 The transfer matrix for the free dofs
///
/// \ingroup Solver
template <typename T>
inline void uniformRefine_withTransfer(
    const gsMultiBasis<T>& mBasis,
    const gsBoundaryConditions<T>& boundaryConditions,
    const gsOptionList& assemblerOptions,
    gsMultiBasis<T> & refinedMBasis,
    gsSparseMatrix<T, RowMajor>& transferMatrix
    )
{
    uniformRefine_withTransfer( mBasis, boundaryConditions, assemblerOptions, 1, 1, refinedMBasis, transferMatrix );
}

/// @brief This function coarsens a knot vector
///
/// @param kv[in]                        The initial knot vector
/// @param removedKnots[out]             The remoed knots
/// @return                              The coarsened knot vector
///
/// \ingroup Solver
template <typename T>
gsKnotVector<T> coarsenKnotVector(const gsKnotVector<T>& kv, std::vector<T>& removedKnots);

/// @brief This function coarsens a tensor B-spline basis
///
/// @param b[in]                         The initial basis
/// @param removedKnots[out]             The remoed knots
/// @return                              The coarsened basis
///
/// \ingroup Solver
template <unsigned d, typename T>
typename gsTensorBSplineBasis<d, T>::uPtr coarsenTensorBasis(const gsTensorBSplineBasis<d, T>& b, std::vector< std::vector<T> > & removedKnots);

/// @brief This function coarsens a tensor B-spline basis
///
/// @param b[in]                         The initial basis
/// @param removedKnots[out]             The remoed knots
/// @return                              The coarsened basis
///
/// \ingroup Solver
template <typename T>
typename gsBasis<T>::uPtr coarsenBasis(const gsBasis<T> & b, std::vector< std::vector<T> >& removedKnots);

/// @brief This function coarsens a tensor B-spline basis
///
/// @param b[in]                         The initial basis
/// @param transferMatrix[out]           The transfer matrix
/// @return                              The coarsened basis
///
/// \ingroup Solver
template <typename T>
typename gsBasis<T>::uPtr coarsenBasis_withTransfer(const gsBasis<T>& b, gsSparseMatrix<T, RowMajor>& transferMatrix);

/// @brief This function coarsens the bases of a gsMultiBasis and provides the transfer matrix.
///
/// @param mBasis[in]                          The gsMultiBasis to be refined
/// @param boundaryConditions[in]              The boundary conditions
/// @param assemblerOptions[in]                A gsOptionList defining a "DirichletStrategy" and a "InterfaceStrategy"
/// @param coarsenedMBasis[out]                A pointer to the coarsened gsMultiBasis
/// @param transferMatrix[out]                 The transfer matrix for the free dofs
/// @param localTransferMatrices[out]          A vector of the local transfer matrices
///
/// \ingroup Solver
template <typename T>
void coarsenMultiBasis_withTransfer(
    const gsMultiBasis<T>& mBasis,
    const gsBoundaryConditions<T>& boundaryConditions,
    const gsOptionList& assemblerOptions,
    gsMultiBasis<T> & coarsenedMBasis,
    gsSparseMatrix<T, RowMajor>& transferMatrix,
    std::vector< gsSparseMatrix<T, RowMajor> >& localTransferMatrices
    );

/// @brief This function coarsens the bases of a gsMultiBasis and provides the transfer matrix.
///
/// @param mBasis[in]                          The gsMultiBasis to be refined
/// @param boundaryConditions[in]              The boundary conditions
/// @param assemblerOptions[in]                A gsOptionList defining a "DirichletStrategy" and a "InterfaceStrategy"
/// @param coarsenedMBasis[out]                A pointer to the coarsened gsMultiBasis
/// @param transferMatrix[out]                 The transfer matrix for the free dofs
///
/// \ingroup Solver
template <typename T>
inline void coarsenMultiBasis_withTransfer(
    const gsMultiBasis<T>& mBasis,
    const gsBoundaryConditions<T>& boundaryConditions,
    const gsOptionList& assemblerOptions,
    gsMultiBasis<T> & coarsenedMBasis,
    gsSparseMatrix<T, RowMajor>& transferMatrix
    )
{
    std::vector< gsSparseMatrix<T, RowMajor> > localTransferMatrices;
    coarsenMultiBasis_withTransfer( mBasis, boundaryConditions, assemblerOptions, coarsenedMBasis, transferMatrix, localTransferMatrices );
}

/// @brief This function takes local transfer matrices (per patch) and combines them
/// using given DofMappers to a global transfer matrix. Simultanously, this
/// function restricts the matrices to the free dofs, e.g., Dirichlet dofs might be
/// eliminated
///
/// @param localTransferMatrices[in]     The local and full (also non-free dofs) transfer matrices per patch
/// @param coarseMapper[in]              The DofMapper on the coarse grid
/// @param fineMapper[in]                The DofMapper on the fine grid
/// @param transferMatrix[out]           The combined transfer matrix restricted to the free dofs
///
/// \ingroup Solver
template <typename T>
void combineTransferMatrices(
    const std::vector< gsSparseMatrix<T, RowMajor> >& localTransferMatrices,
    const gsDofMapper& coarseMapper,
    const gsDofMapper& fineMapper,
    gsSparseMatrix<T, RowMajor>& transferMatrix
    );


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsGridHierarchy.hpp)
#endif
