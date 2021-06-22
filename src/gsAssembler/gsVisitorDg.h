/** @file gsVisitorDg.h

    @brief Visitor for adding the interface conditions for the interior
    penalty methods of the Poisson problem.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer, A. Mantzflaris, S. Moore, S. Takacs
*/

#pragma once

namespace gismo
{
   /** \brief Visitor for adding the interface conditions for the interior
     * penalty methods of the Poisson problem.
     *
     * This visitor adds the following term to the left-hand side (bilinear form).
     * \f[
     *     s(u,v) :=
     *        - \gamma \{\nabla u\} \cdot \mathbf{n} [ v ]
     *        - \beta  \{\nabla v\} \cdot \mathbf{n} [ u ]
     *        + \alpha (h_k^{-1} + h_\ell^{-1}) [ u ][ v ],
     * \f]
     * where \f$ v \f$ is the test function and \f$ u \f$ is trial function,
     * \f$ [u] = u^{(k)} - u^{(\ell)} \f$ denotes the jump accross the interface
     * and \f$ \{ u \} = (u^{(k)} + u^{(\ell)})/2\f$ denotes the average between
     * the two patches.
     *
     * The default values are \f$ \gamma=\beta= 1 \f$ and \f$ \alpha = p^2 \f$,
     * which corresponds to the Symmetric Interior Penalty discontinuous Galerkin
     * (SIPG) method.
     *
     * An analogous visitor for handling the Dirichlet boundary conditions weakly,
     * is the \a gsVisitorNitsche.
     *
     * \ingroup Assembler
     */

template <class T>
class gsVisitorDg
{
public:

    /// Constructor
    gsVisitorDg()
    {}

    /// Constructor. The given Pde is ignored.
    gsVisitorDg(const gsPde<T> &)
    {}

    /// Default options
    static gsOptionList defaultOptions()
    {
        gsOptionList options;
        options.addReal  ("DG.Alpha",        "",    -1);
        options.addReal  ("DG.Beta",         "",     1);
        options.addReal  ("DG.Gamma",        "",     1);
        options.addSwitch("DG.NonSymmetric", "", false);
        return options;
    }

    /// Initialize
    void initialize(const gsBasis<T> & basis1,
                    const gsBasis<T> & basis2,
                    const boundaryInterface & bi,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {

        side1 = bi.first().side();

        // Setup Quadrature
        rule = gsQuadrature::get(basis1, options, side1.direction());

        // Compute penalty parameter
        const int deg = math::max( basis1.maxDegree(), basis2.maxDegree() );

        // If not given, run automatically
        m_alpha = options.askReal("DG.Alpha",-1);
        if (m_alpha<0)
            m_alpha = (deg + basis1.dim()) * (deg + 1) * T(2);

        m_beta  = options.askReal("DG.Beta" , 1);
        m_gamma = options.askReal("DG.Gamma", 1);
        m_nonSymm = options.askSwitch("DG.NonSymmetric", false);

        // Set Geometry evaluation flags
        md1.flags = md2.flags = NEED_VALUE|NEED_JACOBIAN|NEED_GRAD_TRANSFORM;
    }

    /// Evaluate on element
    inline void evaluate(const gsBasis<T>       & B1, // to do: more unknowns
                         const gsGeometry<T>    & geo1,
                         const gsBasis<T>       & B2, // to do: more unknowns
                         const gsGeometry<T>    & geo2,
                         const gsMatrix<T>      & quNodes1,
                         const gsMatrix<T>      & quNodes2)
    {
        md1.points = quNodes1;
        md2.points = quNodes2;
        // Compute the active basis functions
        B1.active_into(md1.points.col(0), actives1);
        B2.active_into(md2.points.col(0), actives2);
        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();

        // Evaluate basis functions and their first derivatives
        B1.evalAllDers_into( md1.points, 1, basisData1);
        B2.evalAllDers_into( md2.points, 1, basisData2);

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geo1.computeMap(md1);
        geo2.computeMap(md2);

        // Initialize local matrices
        B11.setZero(numActive1, numActive1); B12.setZero(numActive1, numActive2);
        E11.setZero(numActive1, numActive1); E12.setZero(numActive1, numActive2);
        B22.setZero(numActive2, numActive2); B21.setZero(numActive2, numActive1);
        E22.setZero(numActive2, numActive2); E21.setZero(numActive2, numActive1);
    }

    /// Assemble on element
    inline void assemble(gsDomainIterator<T>    & element1,
                         gsDomainIterator<T>    & element2,
                         gsVector<T>            & quWeights)
    {
        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute the outer normal vector from patch1
            outerNormal(md1, k, side1, unormal);

            // Integral transformation and quadrature weight (patch1)
            // assumed the same on both sides
            const T weight = quWeights[k] * unormal.norm();

            // Compute the outer unit normal vector from patch1 in place
            unormal.normalize();

            // Take blocks of values and derivatives of basis functions
            //const typename gsMatrix<T>::Column val1 = basisData1.col(k);//bug
            const typename gsMatrix<T>::Block val1 = basisData1[0].block(0,k,numActive1,1);
            gsMatrix<T> & grads1 = basisData1[1];// all grads
            //const typename gsMatrix<T>::Column val2 = basisData2.col(k);//bug
            const typename gsMatrix<T>::Block val2 = basisData2[0].block(0,k,numActive2,1);
            gsMatrix<T> & grads2 = basisData2[1];// all grads

            // Transform the basis gradients
            transformGradients(md1, k, grads1, phGrad1);
            transformGradients(md2, k, grads2, phGrad2);

            // Compute element matrices
            const T c1     = weight * T(0.5);
            N1.noalias()   = unormal.transpose() * phGrad1;
            N2.noalias()   = unormal.transpose() * phGrad2;
            // TODO: m_beta and m_gamma
            // TODO: ieti part
            B11.noalias() += c1 * ( val1 * N1 );
            B12.noalias() += c1 * ( val1 * N2 );
            B22.noalias() -= c1 * ( val2 * N2 );
            B21.noalias() -= c1 * ( val2 * N1 );

            const T h1     = element1.getCellSize();
            const T h2     = element2.getCellSize();
            // Maybe, the h should be scaled with the patch diameter, since its the h from the parameter domain.
            const T c2     = weight * m_alpha * 2*(1./h1 + 1./h2);

            E11.noalias() += c2 * ( val1 * val1.transpose() );
            E12.noalias() += c2 * ( val1 * val2.transpose() );
            E22.noalias() += c2 * ( val2 * val2.transpose() );
            E21.noalias() += c2 * ( val2 * val1.transpose() );
        }
    }

    /// Adds the contirbutions to the sparse system
    inline void localToGlobal(const index_t                     patch1,
                              const index_t                     patch2,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>               & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives1, patch1, actives1);
        system.mapColIndices(actives2, patch2, actives2);

        m_localRhs1.setZero(actives1.rows(),system.rhsCols());
        m_localRhs2.setZero(actives2.rows(),system.rhsCols());

        if (m_nonSymm) // TODO: should be done differently
        {
            system.push((-B11 - B11.transpose() + E11)/T(2), m_localRhs1,actives1,actives1,eliminatedDofs.front(),0,0);
            system.push( -B21                   - E21 /T(2), m_localRhs2,actives2,actives1,eliminatedDofs.front(),0,0);
            system.push(      - B21.transpose() - E12 /T(2), m_localRhs1,actives1,actives2,eliminatedDofs.front(),0,0);
            system.push(                          E22 /T(2), m_localRhs2,actives2,actives2,eliminatedDofs.front(),0,0);
        }
        else
        {
            system.push(-B11 - B11.transpose() + E11, m_localRhs1,actives1,actives1,eliminatedDofs.front(),0,0);
            system.push(-B21 - B12.transpose() - E21, m_localRhs2,actives2,actives1,eliminatedDofs.front(),0,0);
            system.push(-B12 - B21.transpose() - E12, m_localRhs1,actives1,actives2,eliminatedDofs.front(),0,0);
            system.push(-B22 - B22.transpose() + E22, m_localRhs2,actives2,actives2,eliminatedDofs.front(),0,0);
        }
    }
    
private:

    // Penalty constant
    T m_alpha, m_beta, m_gamma;

    bool m_nonSymm;

    // Side
    boxSide side1;

private:

    // Basis values etc
    std::vector<gsMatrix<T> > basisData1, basisData2;
    gsMatrix<T>       phGrad1   , phGrad2;
    gsMatrix<index_t> actives1  , actives2;

    // Outer normal
    gsVector<T> unormal;

    // Auxiliary element matrices
    gsMatrix<T> B11, B12, E11, E12, N1,
                B22, B21, E22, E21, N2;

    gsMatrix<T> m_localRhs1, m_localRhs2;

    gsMapData<T> md1, md2;
};


} // namespace gismo
