/** @file gsIetiMapper.h

    @brief Algorithms that help with assembling the matrices required for IETI-solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsCore/gsMultiBasis.h>
#include <gsAssembler/gsExprAssembler.h>

namespace gismo
{

/** @brief Ieti Mapper
 *
 *  Algorithms that help with assembling the matrices required for IETI-Solvers
 *
 *  This class is written to work with the expression assembler and with assemblers
 *  derived from \a gsAssembler. If applied to a system, it is expected that
 *  individual instances of this class are used for each of the variables.
 *
 *  The objects of this class are initialized using a global dof mapper and a
 *  vector that contains function values for the eliminated variables (usually
 *  for the Dirichlet boundary).
 *
 *  This class then allows to obtain the jump matrices (\ref jumpMatrix), the
 *  patch-local dof mappers (\ref dofMapperLocal) and the patch-local function
 *  values for the eliminated dofs (\ref fixedPart). The member function
 *  \ref initFeSpace allows to pass dof mapper and the function values for the
 *  eliminated dofs to a variable object of the \a gsExprAssembler.
 *
 *  Moreover, this class allows to construct the primal degrees of freedom,
 *  which are then handled by the class \a gsPrimalSystem.
 *
 *  Finally, the member \ref constructGlobalSolutionFromLocalSolutions allows
 *  the combination of the patch-local solutions to a global one.
 *
 *  @ingroup Solver
*/
template< typename T >
class gsIetiMapper
{
    typedef gsMatrix<T>                Matrix;        ///< Matrix type
    typedef gsSparseVector<T>          SparseVector;  ///< Sparse matrix type
    typedef gsSparseMatrix<T,RowMajor> JumpMatrix;    ///< Sparse matrix type for jumps

public:
    /// @brief Default constructor
    gsIetiMapper() : m_multiBasis(NULL), m_nPrimalDofs(0), m_status(0) {}

    /// @brief Create the ieti mapper
    ///
    /// @param multiBasis     The corresponding \a gsMultiBasis
    /// @param gsDofMapper    The dof mapper for the global problem
    /// @param fixedPart      The function values for the eliminated dofs
    ///
    /// The multibasis object has to outlive the ieti mapper.
    gsIetiMapper(
        const gsMultiBasis<T>& multiBasis,
        gsDofMapper dofMapperGlobal,
        const Matrix& fixedPart
    )
    { init(multiBasis, give(dofMapperGlobal), fixedPart); }

    /// @brief Init the ieti mapper after default construction
    ///
    /// @param multiBasis     The corresponding \a gsMultiBasis
    /// @param gsDofMapper    The dof mapper for the global problem
    /// @param fixedPart      The function values for the eliminated dofs
    ///
    /// Instances of the class can be reused using this member function.
    void init(
        const gsMultiBasis<T>& multiBasis,
        gsDofMapper dofMapperGlobal,
        const Matrix& fixedPart
    );

    /// @brief Apply the required changes to a space object of the expression
    /// assembler
    ///
    /// It is assumed that the space object is fully functioning, i.e., the
    /// setup member has been called.
    ///
    /// This function exposes the \ref dofMapperLocal and the \ref fixedPart to
    /// the space.
    void initFeSpace(typename gsExprAssembler<T>::space u, index_t k)
    {
        GISMO_ASSERT( m_status&1,
            "gsIetiMapper: The class has not been initialized." );
        GISMO_ASSERT( u.mapper().size() ==  m_dofMapperLocal[k].size(),
            "gsIetiMapper::initFeSpace: The sizes do not agree." );
        const_cast<expr::gsFeSpace<T>&>(u).mapper() = m_dofMapperLocal[k];
        const_cast<expr::gsFeSpace<T>&>(u).fixedPart() = m_fixedPart[k];
    }

    /// @brief This function computes the jump matrices
    ///
    /// @param fullyRedundant  Compute the jump matrices in a fullyRedundant way;
    ///                        if false, then no redundancy
    /// @param excludeCorners  Ignore corners for jump matrices. This makes sense
    ///                        if the corners are chosen as primal dofs
    void computeJumpMatrices(bool fullyRedundant, bool excludeCorners);

    /// @brief Set up the corners as primal dofs
    void cornersAsPrimals();

    /// @brief Set up interface averages as primal dofs
    ///
    /// @param geo             \a gsMultiPatch object describing the geometry
    /// @param d               The dimension of the interfaces to be considered:
    ///                        d=1 yields edge averages, d=2 yields face averages
    ///
    /// Adds contributions only iff the interface (edge, face,...) belongs to at
    /// least two patches. This function also works if d=geo.dim(), i.e., for the
    /// patch averages. In this case, the requirement that it contributes to
    /// at least two patches does not apply. For corners, use \ref cornersAsPrimals
    /// instead.
    void interfaceAveragesAsPrimals(const gsMultiPatch<T>& geo, short_t d);

    /// @brief With this function, the caller can register more primal constraints
    ///
    /// All primal constraints in the vector are considered to be one single primal
    /// dof. Each component of the vector contains a pair consisting of the patch
    /// index and the vector representing the primal constraint.
    void customPrimalConstraints( std::vector< std::pair<index_t,SparseVector> > data );

    /// @brief Construct the global solution from a vector of patch-local ones
    Matrix constructGlobalSolutionFromLocalSolutions( const std::vector<Matrix>& localContribs );

    /// @brief Returns a list of dofs that are (on the coarse level) coupled
    ///
    /// @param patch   Number of the patch
    std::vector<index_t> skeletonDofs( index_t patch ) const;

public:
    /// @brief Returns the number of Lagrange multipliers.
    index_t nLagrangeMultipliers()
    {
        GISMO_ASSERT(! m_jumpMatrices.empty(),
            "gsIetiMapper: Number of Lagrange multipliers not yet known.");
        return m_jumpMatrices[0].rows();
    }

    /// @brief Returns the size of the primal problem (number of primal dofs)
    index_t nPrimalDofs() const                                            { return m_nPrimalDofs;                }

    /// @brief Returns the primalConstraints (as vectors) for the given patch
    ///
    /// These vectors form the matrix \f$ C_k \f$ in the local saddle point
    /// system, see  \a gsPrimalSystem
    const std::vector<SparseVector> & primalConstraints(index_t k) const   { return m_primalConstraints[k];       }

    /// @brief Returns the indices of the primal dofs that are associated to
    /// the primal constraints for the given patch
    const std::vector<index_t> & primalDofIndices(index_t k) const         { return m_primalDofIndices[k];        }

    /// @brief Returns the jump matrix \f$ B_k \f$ for the given patch
    const JumpMatrix& jumpMatrix(index_t k) const                          { return m_jumpMatrices[k];            }

    /// @brief The global dof mapper
    const gsDofMapper& dofMapperGlobal() const                             { return m_dofMapperGlobal;            }

    /// @brief The dof mapper for the given patch
    const gsDofMapper& dofMapperLocal(index_t k) const                     { return m_dofMapperLocal[k];          }

    /// @brief The function values for the eliminated dofs on the given patch
    const Matrix& fixedPart(index_t k) const                               { return m_fixedPart[k];               }

private:

    static gsSparseVector<T> assembleAverage( const gsGeometry<T>& geo,
        const gsBasis<T>& basis, const gsDofMapper& dm,
        boxComponent bc );   ///< Assembles for \ref interfaceAveragesAsPrimals

private:
    const gsMultiBasis<T>*                        m_multiBasis;          ///< Pointer to the respective multibasis
    gsDofMapper                                   m_dofMapperGlobal;     ///< The global dof mapper
    std::vector<gsDofMapper>                      m_dofMapperLocal;      ///< A vector of the patch-local dof mappers
    std::vector<Matrix>                           m_fixedPart;           ///< The values for the elminated (Dirichlet) dofs
    std::vector<JumpMatrix>                       m_jumpMatrices;        ///< The jump matrices
    index_t                                       m_nPrimalDofs;         ///< The number of primal dofs already created
    std::vector< std::vector<SparseVector> >      m_primalConstraints;   ///< The primal constraints
    std::vector< std::vector<index_t> >           m_primalDofIndices;    ///< The primal dof indices for each of the primal constraints
    unsigned                                      m_status;              ///< A status flag that is checked by assertions
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsIetiMapper.hpp)
#endif
