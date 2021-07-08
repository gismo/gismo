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
     * This visitor adds the following term to the left-hand side (bilinear form):
     * \f[
     *     s_{k,\ell}(u,v) :=
     *        - \alpha \int_{\Gamma^{(k,\ell)}}  \{\nabla u\} \cdot \mathbf{n} [ v ]
     *        - \beta  \int_{\Gamma^{(k,\ell)}}  \{\nabla v\} \cdot \mathbf{n} [ u ]
     *        + \delta ( h_k^{-1} + h_\ell^{-1}) \int_{\Gamma^{(k,\ell)}}  [ u ][ v ],
     * \f]
     * where \f$ v \f$ is the test function and \f$ u \f$ is trial function,
     * \f$ [u] = u^{(k)} - u^{(\ell)} \f$ denotes the jump accross the interface
     * and \f$ \{ u \} = (u^{(k)} + u^{(\ell)})/2\f$ denotes the average between
     * the two patches.
     *
     * The default values for DG.Alpha and DG.Beta are \f$\alpha=\beta=1\f$,
     * which corresponds to the Symmetric Interior Penalty discontinuous Galerkin
     * (SIPG) method. The default value for DG.Penalty is -1, which yields
     * \f$\delta=4(p+d)(p+1)\f$. If a positive value for DG.Penalty is chosen,
     * that value is taken for \f$\delta\f$.
     *
     * The (normal) grid sizes are estimated based on the geometry function. Use
     * the option DG.ParameterGridSize to just use the grid size on the parameter
     * domain.
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
    : m_pde(NULL)
    {}

    /// @brief Constructor
    ///
    /// @param pde     Reference to \a gsPde object
    gsVisitorDg(const gsPde<T> & pde)
    : m_pde(&pde)
    {}

    /// Default options
    static gsOptionList defaultOptions()
    {
        gsOptionList options;
        options.addReal  ("DG.Alpha",
                          "Parameter alpha for dG scheme; use 1 for SIPG and NIPG.",                              1);
        options.addReal  ("DG.Beta",
                          "Parameter beta for dG scheme; use 1 for SIPG and -1 for NIPG",                         1);
        options.addReal  ("DG.Penalty",
                          "Penalty parameter penalty for dG scheme; if negative, default 4(p+d)(p+1) is used.",  -1);
        options.addSwitch("DG.ParameterGridSize",
                          "Use grid size on parameter domain for the penalty term.",                          false);
        return options;
    }

    /// Initialize
    void initialize(const gsBasis<T> & basis1,
                    const gsBasis<T> & basis2,
                    const boundaryInterface & bi,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        m_side1 = bi.first();
        m_side2 = bi.second();

        // Setup Quadrature
        rule = gsQuadrature::get(basis1, options, m_side1.direction());

        m_penalty     = options.askReal("DG.Penalty",-1);
        // If not given, use default
        if (m_penalty<0)
        {
            const index_t deg = math::max( basis1.maxDegree(), basis2.maxDegree() );
            m_penalty = T(4) * (deg + basis1.dim()) * (deg + 1);
        }

        m_alpha     = options.askReal("DG.Alpha", 1);
        m_beta      = options.askReal("DG.Beta" , 1);

        if (options.getSwitch("DG.ParameterGridSize"))
        {
            m_h1    = basis1.getMinCellLength();
            m_h2    = basis2.getMinCellLength();
        }
        else
        {
            GISMO_ENSURE (m_pde, "gsVisitorDg::initialize: No PDE given.");
            m_h1    = estimateSmallestPerpendicularCellSize(
                          basis1,
                          m_pde->patches()[m_side1.patch],
                          m_side1
                      );
            m_h2    = estimateSmallestPerpendicularCellSize(
                          basis2,
                          m_pde->patches()[m_side2.patch],
                          m_side2
                      );
        }

        // Set Geometry evaluation flags
        md1.flags = md2.flags = NEED_VALUE|NEED_JACOBIAN|NEED_GRAD_TRANSFORM;
    }

    /// Evaluate on element
    inline void evaluate(const gsBasis<T>       & B1,
                         const gsGeometry<T>    & geo1,
                         const gsBasis<T>       & B2,
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
        B1.evalAllDers_into(md1.points, 1, basisData1);
        B2.evalAllDers_into(md2.points, 1, basisData2);

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
            outerNormal(md1, k, m_side1, unormal);

            // Multiply quadrature weight by the geometry measure
            // Integral transformation and quadrature weight (patch1)
            // assumed the same on both sides
            const T weight = quWeights[k] * unormal.norm();

            // Compute the outer unit normal vector from patch1 in place
            unormal.normalize();

            // Take blocks of values and derivatives of basis functions
            const typename gsMatrix<T>::Block val1 = basisData1[0].block(0,k,numActive1,1);
            gsMatrix<T> & grads1 = basisData1[1];// all grads
            const typename gsMatrix<T>::Block val2 = basisData2[0].block(0,k,numActive2,1);
            gsMatrix<T> & grads2 = basisData2[1];// all grads

            // Transform the basis gradients
            transformGradients(md1, k, grads1, phGrad1);
            transformGradients(md2, k, grads2, phGrad2);

            // Compute element matrices
            const T c1     = weight / T(2);
            N1.noalias()   = unormal.transpose() * phGrad1;
            N2.noalias()   = unormal.transpose() * phGrad2;

            B11.noalias() += c1 * ( val1 * N1 );
            B12.noalias() += c1 * ( val1 * N2 );
            B21.noalias() -= c1 * ( val2 * N1 );
            B22.noalias() -= c1 * ( val2 * N2 );

            const T c2     = weight * m_penalty * (1./m_h1 + 1./m_h2);

            E11.noalias() += c2 * ( val1 * val1.transpose() );
            E12.noalias() += c2 * ( val1 * val2.transpose() );
            E22.noalias() += c2 * ( val2 * val2.transpose() );
            E21.noalias() += c2 * ( val2 * val1.transpose() );

        }
    }

    /// Adds the contributions to the sparse system
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

        system.push(-m_alpha*B11 - m_beta*B11.transpose() + E11, m_localRhs1,actives1,actives1,eliminatedDofs.front(),0,0);
        system.push(-m_alpha*B21 - m_beta*B12.transpose() - E21, m_localRhs2,actives2,actives1,eliminatedDofs.front(),0,0);
        system.push(-m_alpha*B12 - m_beta*B21.transpose() - E12, m_localRhs1,actives1,actives2,eliminatedDofs.front(),0,0);
        system.push(-m_alpha*B22 - m_beta*B22.transpose() + E22, m_localRhs2,actives2,actives2,eliminatedDofs.front(),0,0);
    }

    /// Estimates the gird size perpendicular to the given side on the physical domain, as required for SIPG and Nitsche
    static T estimateSmallestPerpendicularCellSize(const gsBasis<T> & basis,
                                                   const gsGeometry<T> & geo,
                                                   patchSide side)
    {
        T result = 0;
        bool first = true;
        typename gsDomainIterator<T>::uPtr domIt = basis.makeDomainIterator(side);

        for (gsDomainIterator<T>& element = *domIt; element.good(); element.next() )
        {
            gsVector<T> parameterCenter = element.centerPoint();
            const gsMatrix<T> center1 = geo.eval(parameterCenter);

            parameterCenter(side.direction()) += ( side.parameter() == 0 ? 1 : -1 )
                                                 * element.getPerpendicularCellSize();
            const gsMatrix<T> center2 = geo.eval(parameterCenter);

            const real_t diff = (center1 - center2).norm();

            if (first || result > diff)
                result = diff;

            first = false;
        }
        return result;
    }



private:

    /// The underlying PDE
    const gsPde<T> * m_pde;

    /// Parameter \f$\alpha\f$ for the bilinear form
    T m_alpha;

    /// Parameter \f$\beta\f$ for the bilinear form
    T m_beta;

    /// Parameter \f$\delta\f$ for the bilinear form
    T m_penalty;

    /// Side on first patch that corresponds to interface
    patchSide m_side1;

    /// Side on second patch that corresponds to interface
    patchSide m_side2;

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

    // Grid sizes
    T m_h1, m_h2;
};


} // namespace gismo
