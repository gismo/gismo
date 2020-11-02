/** @file gsVisitorDg.h

    @brief A DG interface visitor for the Poisson problem .

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer, A. Mantzflaris, S. Moore, S. Takacs
*/

#pragma once

namespace gismo
{
/** @brief
    Implementation of a interface condition for the
    discontinuous Galerkin Assembler.

    It sets up an assembler and assembles the system patch wise and
    combines the patch-local stiffness matrices into a global system.
    Dirichlet boundary can also be imposed weakly (i.e. Nitsche).
*/

template <class T>
class gsVisitorMixed
{
public:

   /** \brief Visitor for adding the interface conditions for the
     * interior penalty method of the Poisson problem.
     *
     * This visitor adds the following term to the left-hand side (bilinear form).
     * \f[ - \{\nabla u\} \cdot \mathbf{n} [ v ]
     *     - \{\nabla v\} \cdot \mathbf{n} [ u ]
     *     + \alpha [ u ][  v ] \f].
     * Where \f[v\f] is the test function and \f[ u \f] is trial function.
     */

   gsVisitorMixed()
    {}

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T>    & rule)
    {
        gsVector<index_t> numQuadNodes( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i) // to do: improve
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// NB!

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;

    }

    // Evaluate on element.
    inline void evaluate(const gsBasis<T>       & B1, // to do: more unknowns
                         const gsGeometry<T>    & geo1,
                         const gsBasis<T>       & B2, // to do: more unknowns
                         const gsMatrix<T>      & quNodes1)
    {
        md.points = quNodes1;

        // Compute the active basis functions
        B1.active_into(md.points.col(0), actives1);
        B2.active_into(md.points.col(0), actives2);

        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();

        // Evaluate basis functions and their first derivatives
        B1.evalAllDers_into( md.points, 2, basisData1);
        B2.evalAllDers_into( md.points, 2, basisData2);

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geo1.computeMap(md);

        // Initialize local matrices
        localMat12.setZero(numActive1, numActive2);
        localMat21.setZero(numActive2, numActive1);
    }

    // assemble on element
    inline void assemble(gsDomainIterator<T>    & element1,
                         gsDomainIterator<T>    & element2,
                         gsVector<T>            & quWeights)
    {

        //gsMatrix<T> & basisVals1  = basisData1[0];
        gsMatrix<T> & basisGrads1 = basisData1[1];
        gsMatrix<T> & basis2ndDerivs1 = basisData1[2];

        //gsMatrix<T> & basisVals2  = basisData2[0];
        gsMatrix<T> & basisGrads2 = basisData2[1];
        gsMatrix<T> & basis2ndDerivs2 = basisData2[2];

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            const T weight = quWeights[k] * md.measure(k);

            // Compute physical laplacian at k as a 1 x numActive matrix
            if(md.dim.first == md.dim.second)
            {
                transformLaplaceHgrad(md, k, basisGrads1, basis2ndDerivs1, physBasisLaplace1);
                transformLaplaceHgrad(md, k, basisGrads2, basis2ndDerivs2, physBasisLaplace2);

                localMat12.noalias() += weight * (physBasisLaplace1.transpose() * physBasisLaplace2);
                localMat21.noalias() += weight * (physBasisLaplace2.transpose() * physBasisLaplace1);

            }
        }
    }

    inline void localToGlobal(const int                         patch1,
                              const int                         patch2,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>               & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives1, patch1, actives1);
        system.mapColIndices(actives2, patch2, actives2);

        //gsInfo << localMat12.dim() << " : " << localMat21.dim() << " : " << actives1.transpose() << " : " << actives2.transpose() << "\n";

        system.pushToMatrix(localMat12, actives1, actives2, eliminatedDofs[0], 0, 0);
        system.pushToMatrix(localMat21, actives2, actives1, eliminatedDofs[0], 0, 0);

    }

protected:
    // Right hand side
    const gsFunction<T> * rhs_ptr;

protected:
    // Basis values

    std::vector<gsMatrix<T> > basisData1, basisData2;
    gsMatrix<T>        physBasisLaplace1, physBasisLaplace2;

    gsMatrix<unsigned> actives1, actives2;


protected:
    // Local values of the right hand side
    gsMatrix<T> rhsVals;

protected:
    // Local matrices
    gsMatrix<T> localMat12, localMat21;

    gsMapData<T> md;
};


} // namespace gismo
