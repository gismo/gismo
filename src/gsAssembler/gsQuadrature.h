/** @file gsQuadrature.h

    @brief Creates a variety of quadrature rules

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsIO/gsOptionList.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsAssembler/gsLobattoRule.h>
#include <gsAssembler/gsNewtonCotesRule.h>
#include <gsAssembler/gsPatchRule.h>
#include <gsAssembler/gsOverIntegrateRule.h>
#include <gsAssembler/gsGaussRule.h>

#include <gsCore/gsFunction.h>
#include <gsDomain/gsDomainIterator.h>
#include <gsDomain/gsImplicitTrimmedDomain.h>

#ifdef gsAlgoim_ENABLED
#include <gsAlgoim/gsAlgoimRule.h>
#endif

namespace gismo
{

template<class T>
class gsCutCellRule;

template<class T>
class gsOctreeCutCellRule;

template<class T>
class gsCutCellSurfaceRule;

template<class T>
class gsOctreeCutCellSurfaceRule;

/**
 * @file gsQuadrature.h
 * @brief This file contains the definition and implementation of various quadrature rules and methods for constructing quadrature rules based on input options.
 *
 * The quadrature rules supported include:
 * - Gauss-Legendre quadrature
 * - Gauss-Lobatto quadrature
 * - Patch-wise quadrature rule (Johannessen 2017)
 *
 * The file provides functions to construct quadrature rules, retrieve quadrature nodes, and evaluate them using geometries.
 *
 * Reference:
 * Johannessen, K. A. (2017). Optimal quadrature for univariate and tensor product splines.
 * Computer Methods in Applied Mechanics and Engineering, 316, 84–99.
 * https://doi.org/10.1016/j.cma.2016.04.030
 *
 * @typedef Real
 * Type alias for GISMO_COEFF_TYPE.
 *
 * @enum rule
 * Enumeration of quadrature rule types.
 *
 * @function get
 * Constructs a quadrature rule based on input options.
 *
 * @function getPtr
 * Constructs a quadrature rule based on input options and returns a unique pointer to the quadrature rule.
 *
 * @function getUnivariate
 * Constructs a univariate quadrature rule based on input options.
 *
 * @function numNodes
 * Computes the number of nodes for the quadrature rule based on the domain and input options.
 *
 * @function getAllNodes
 * Retrieves all quadrature nodes for the given basis or domain.
 *
 * @function getAllNodes
 * Retrieves all quadrature nodes for a specified side of a given basis.
 *
 * @function getAllNodes
 * Retrieves all quadrature nodes for a specified side of a basis and evaluates them using a geometry.
 *
 * @function getAllNodes
 * Retrieves all quadrature nodes for multiple sides of a given basis.
 *
 * @function getAllNodes
 * Collects and evaluates all quadrature nodes for multiple sides of a given basis.
 *
 * @function getAllNodes
 * Collects all quadrature nodes for a multi-basis.
 *
 * @function getAllNodes
 * Collects all quadrature nodes for a multi-basis and evaluates them using a multi-patch geometry.
 */
struct gsQuadrature
{
    typedef GISMO_COEFF_TYPE Real;

    /// Quadrature rule types
    enum rule
    {
        GaussLegendre = 1, ///< Gauss-Legendre quadrature
        GaussLobatto  = 2, ///< Gauss-Lobatto quadrature
        PatchRule     = 3, ///< Patch-wise quadrature rule  (Johannessen 2017)
        CutCellRule   = 11,
        AlgoimRule    = 12,
        OctreeRule    = 13
    };
    /*
    Reference:
        Johannessen, K. A. (2017). Optimal quadrature for univariate and tensor product splines.
        Computer Methods in Applied Mechanics and Engineering, 316, 84–99.
        https://doi.org/10.1016/j.cma.2016.04.030
    */

    /// Constructs a quadrature rule based on input \a options
    template<class T>
    static gsQuadRule<T> get(const gsBasis<T> & basis,
                             const gsOptionList & options, short_t fixDir = -1)
    {
        return get<T>(*basis.domain(),options,fixDir);
    }

    template<class T>
    static gsQuadRule<T> get(const gsDomain<T> & domain,
                             const gsOptionList & options, short_t fixDir = -1)
    {
        const index_t qu  = options.askInt("quRule", GaussLegendre);
        const Real    quA = options.getReal("quA");
        const index_t quB = options.getInt ("quB");
        const gsVector<index_t> nnodes = numNodes(domain,quA,quB,fixDir);
        return get<T>(qu, nnodes);
    }

    /// Constructs a quadrature rule based on input \a options
    template<class T>
    static typename gsQuadRule<T>::uPtr
                      getPtr(const gsBasis<T> & basis,
                             const gsOptionList & options, short_t fixDir = -1)
    {
        return getPtr<T>(*basis.domain(),options,fixDir);
    }

    /// Constructs a quadrature rule based on input \a options
    template<class T>
    static typename gsQuadRule<T>::uPtr
                      getPtr(const gsDomain<T> & domain,
                             const gsOptionList & options, short_t fixDir = -1)
    {
        const index_t qu   = options.askInt("quRule", GaussLegendre);
        const Real    quA  = options.getReal("quA");
        const index_t quB  = options.getInt ("quB");
        const bool    over = options.askSwitch ("overInt", false);  // use overintegration?

        if ( (qu==GaussLegendre || qu==GaussLobatto) )
        {
            if (!over)
            {
                switch (qu)
                {
                    case GaussLegendre :
                        return gsGaussRule<T>::make(numNodes(domain,quA,quB,fixDir));
                    case GaussLobatto :
                        return gsLobattoRule<T>::make(numNodes(domain,quA,quB,fixDir));
                    default:
                        GISMO_ERROR("Invalid Quadrature rule request ("<<qu<<")");
                };
            }
            else
            {
                /*
                    Uses quadrature rule with quA and quB for the interior
                    elements and one with quAb and quBb for the boundary elements
                */
                const Real    quAb  = options.askReal("quAb",quA+1);
                const index_t quBb  = options.askInt ("quBb",quB);

                const gsVector<index_t> nnodesI = numNodes(domain,quA,quB,fixDir);
                const gsVector<index_t> nnodesB = numNodes(domain,quAb,quBb,fixDir);

                std::vector<gsQuadRule<T> > quInterior(nnodesI.size());
                std::vector<gsQuadRule<T> > quBoundary(nnodesB.size());

                for (index_t d = 0; d != nnodesI.size(); d++)
                {
                    quInterior[d] = getUnivariate<T>(qu,nnodesI[d]);
                    quBoundary[d] = getUnivariate<T>(qu,nnodesB[d]);
                }

                return gsOverIntegrateRule<T>::make(domain,quInterior,quBoundary);
            }
        }
        else if (qu==PatchRule)
        {
            // quA: Order of the target space
            // quB: Regularity of the target space
            return gsPatchRule<T>::make(domain,cast<T,index_t>(quA),quB,over,fixDir);
        }
        else if (qu==CutCellRule)
        {
            // quDim: -1 = volumetric (mask phi>0); >=0 = surface (integrate phi==0)
            const index_t quDim = options.askInt("quDim", -1);
            return makeCutCellPtr<T>(domain, quA, quB, fixDir, quDim);
        }
        else if (qu==OctreeRule)
        {
            return makeOctreePtr<T>(domain, options, quA, quB, fixDir);
        }
#ifdef gsAlgoim_ENABLED
        else if (qu==AlgoimRule)
        {
            // Forward the option list so the immersed boundary assembly can
            // request surface quadrature (dim == D) via the "quDim" option.
            return gsAlgoimGenericRule<T>::make(domain, options);
        }
#endif
        else
        {
            GISMO_ERROR("Quadrature with index "<<qu<<" unknown.");
        }
    }

    /// Constructs a quadrature rule based on input \a options
    template<class T>
    static inline gsQuadRule<T> get(index_t qu, gsVector<index_t> const & numNodes, unsigned digits = 0)
    {
        switch (qu)
        {
        case GaussLegendre :
            return gsGaussRule<T>(numNodes, digits);
        case GaussLobatto :
            return gsLobattoRule<T>(numNodes, digits);
        default:
            GISMO_ERROR("Invalid Quadrature rule request ("<<qu<<")");
        };
    }

    /// Constructs a quadrature rule based on input \a options
    template<class T>
    static inline gsQuadRule<T> getUnivariate(index_t qu, index_t numNodes, unsigned digits = 0)
    {
        switch (qu)
        {
        case GaussLegendre :
            return gsGaussRule<T>(numNodes, digits);
        case GaussLobatto :
            return gsLobattoRule<T>(numNodes, digits);
        default:
            GISMO_ERROR("Invalid Quadrature rule request ("<<qu<<")");
        };
    }

    /// Computes and integer quA*deg_i + quB where deg_i is the degree
    /// of \a domain
    template<class T>
    static gsVector<index_t> numNodes(const gsBasis<T> & basis,
                               const Real quA, const index_t quB, short_t fixDir = -1)
    {
        return numNodes(*basis.domain(),quA,quB,fixDir);
    }

    /// Computes and integer quA*deg_i + quB where deg_i is the degree
    /// of \a domain
    template<class T>
    static gsVector<index_t> numNodes(const gsDomain<T> & domain,
                               const Real quA, const index_t quB, short_t fixDir = -1)
    {
        const short_t d  = domain.dim();
        GISMO_ASSERT( fixDir < d && fixDir>-2, "Invalid input fixDir = "<<fixDir);
        gsVector<index_t> nnodes(d);

        if (-1==fixDir)
            fixDir = d;
        else
            nnodes[fixDir] = 1;

        short_t i;
        for(i=0; i!=fixDir; ++i )
            //note: +0.5 for rounding
            nnodes[i] = cast<Real,index_t>(quA * domain.degree(i) + quB + 0.5);
        for(++i; i<d; ++i )
            nnodes[i] = cast<Real,index_t>(quA * domain.degree(i) + quB + 0.5);
        return nnodes;
    }

    // template<class T>
    // static std::pair<gsMatrix<T>,gsVector<T>> getAllNodesAndWeights(const gsBasis<T> & basis,
    //                          const gsOptionList & options)
    // {

    // }

    /**
     * @brief Retrieves all quadrature nodes for the given basis.
     *
     * This function computes and returns all quadrature nodes for a given
     * \p domain, using the provided \p options to determine the quadrature rules.
     *
     * @tparam T          Real type.
     * @param[in] domain  The domain for which quadrature nodes are computed.
     * @param[in] options Options specifying the quadrature rule.
     * @return            A matrix where each column represents a quadrature node in the parametric domain.
     */
    template<class T>
    static gsMatrix<T> getAllNodes(const gsDomain<T> & domain,
                                   const gsOptionList & options)
    {
        typename gsBasis<T>::domainIter domIt    = domain.beginAll();
        typename gsBasis<    T>::    domainIter domItEnd = domain.endAll();

        index_t     quadSize = 0;
        typename gsQuadRule<T>::uPtr QuRule;
        QuRule = getPtr(domain, options);

        for (; domIt<domItEnd; ++domIt )
        {
            QuRule = gsQuadrature::getPtr(domain, options);
            quadSize+=QuRule->numNodes();
        }

        gsMatrix<T> result(domain.dim(),quadSize);

        index_t offset = 0;
        gsMatrix<T> nodes;
        gsVector<T> weights;

        domIt = domain.beginAll();
        for (; domIt<domItEnd; ++domIt )
        {
            QuRule = gsQuadrature::getPtr(domain, options);
            // Map the Quadrature rule to the element
            QuRule->mapTo( domIt.lowerCorner(), domIt.upperCorner(),
                           nodes, weights);
            result.block(0,offset,domain.dim(),QuRule->numNodes()) = nodes;
            offset += QuRule->numNodes();
        }
        return result;
    }

    /**
     * @brief Retrieves all quadrature nodes for the given basis and evaluates them using a geometry.
     *
     * @tparam T        Real type.
     * @param[in] basis The basis for which quadrature nodes are computed.
     * @param[in] geom  The geometry used to evaluate the quadrature nodes.
     * @param[in] options Options specifying the quadrature rule.
     * @return result   A matrix of quadrature nodes, where each column corresponds to a quadrature node.
     */
    template<class T>
    static gsMatrix<T> getAllNodes(const gsBasis<T> & basis,
                                   const gsOptionList & options)
    {
        return getAllNodes(*basis.domain(),options);
    }

    /**
     * @brief Get all quadrature nodes for a specified side of a given basis.
     *
     * @tparam T        Real type.
     * @param[in] basis     The basis for which the quadrature nodes are to be collected.
     * @param[in] options   Quadrature rule.
     * @param[in] side      The side of the basis.
     * @return result   A matrix of quadrature nodes, where each column corresponds to a quadrature node.
     */
    template<class T>
    static gsMatrix<T> getAllNodes( const gsBasis<T> & basis,
                                    const gsOptionList & options,
                                    const patchSide side)
    {
        return getAllNodes(*basis.domain(),options,side);
    }

    /**
     * @brief Get all quadrature nodes for a specified side of a given domain.
     *
     * @tparam T        Real type.
     * @param[in] domain    The domain for which the quadrature nodes are to be collected.
     * @param[in] options   Quadrature rule.
     * @param[in] side      The side of the domain.
     * @return result   A matrix of quadrature nodes, where each column corresponds to a quadrature node.
     */
    template<class T>
    static gsMatrix<T> getAllNodes( const gsDomain<T> & domain,
                                    const gsOptionList & options,
                                    const patchSide side)
    {
        typename gsBasis<T>::domainIter domIt    = domain.beginBdr(side.side());
        typename gsBasis<T>::domainIter domItEnd = domain.endBdr(side.side());

        index_t quadSize = 0;
        typename gsQuadRule<T>::uPtr QuRule;
        QuRule = getPtr(domain, options, side.side().direction());

        // First pass: count boundary elements
        for (; domIt<domItEnd; ++domIt )
            quadSize += QuRule->numNodes();

        gsMatrix<T> result(domain.dim(), quadSize);

        // Second pass: collect boundary quadrature points
        index_t offset = 0;
        gsMatrix<T> nodes;
        gsVector<T> weights;
        domIt = domain.beginBdr(side.side());
        for (; domIt<domItEnd; ++domIt )
        {
            QuRule = gsQuadrature::getPtr(domain, options, side.side().direction());
            // Map the Quadrature rule to the element
            QuRule->mapTo( domIt.lowerCorner(), domIt.upperCorner(),
                            nodes, weights);
            result.block(0,offset,domain.dim(),QuRule->numNodes()) = nodes;
            offset += QuRule->numNodes();
        }
        return result;
    }

    /**
     * @brief Retrieves all quadrature nodes for multiple sides of a given domain.
     */
    template<class T>
    static gsMatrix<T> getAllNodes(const gsBasis<T> & basis,
                                   const gsOptionList & options,
                                   const std::vector<patchSide> & sides)
    {
        return getAllNodes(*basis.domain(),options,sides);
    }

    /**
     * @brief Collects all quadrature nodes for a multi-basis.
     */
    template<class T>
    static gsMatrix<T> getAllNodes( const gsDomain<T> & domain,
                                    const gsOptionList & options,
                                    const std::vector<patchSide> & sides)
    {
        std::vector<gsMatrix<T>> nodes(sides.size());
        index_t cols = 0;
        for (size_t s = 0; s != sides.size(); s++)
        {
            nodes[s] = getAllNodes(domain,options,sides[s]);
            cols += nodes[s].cols();
        }
        gsMatrix<T> result(domain.dim(),cols);
        cols = 0;

        for (size_t s = 0; s != sides.size(); s++)
        {
            result.block(0,cols,nodes[s].rows(),nodes[s].cols()) = nodes[s];
            cols += nodes[s].cols();
        }

        return result;
    }

private:
    template<class T>
    static typename gsQuadRule<T>::uPtr
    makeCutCellPtr(const gsDomain<T> & domain,
                   const Real quA,
                   const index_t quB,
                   short_t fixDir,
                   index_t quDim = -1)
    {
        const gsVector<index_t> nnodes = numNodes(domain,quA,quB,fixDir);
        gsGaussRule<T> baseRule(nnodes);

        const gsFunction<T> * levelSet = nullptr;
        if (const auto * d1 = dynamic_cast<const gsImplicitTrimmedDomain<1,T>*>(&domain))
            levelSet = &d1->implicitFunction();
        else if (const auto * d2 = dynamic_cast<const gsImplicitTrimmedDomain<2,T>*>(&domain))
            levelSet = &d2->implicitFunction();
        else if (const auto * d3 = dynamic_cast<const gsImplicitTrimmedDomain<3,T>*>(&domain))
            levelSet = &d3->implicitFunction();

        if (!levelSet)
        {
            static bool warned = false;
            if (!warned)
            {
                gsWarn << "CutCellRule requested on a non-implicit domain; "
                       << "falling back to GaussRule for this integration context.\n";
                warned = true;
            }
            return typename gsQuadRule<T>::uPtr(new gsGaussRule<T>(baseRule));
        }

        if (quDim >= 0) // surface quadrature on the interface phi == 0
            return typename gsQuadRule<T>::uPtr(
                new gsCutCellSurfaceRule<T>(*levelSet, nnodes.maxCoeff()));

        return typename gsQuadRule<T>::uPtr(new gsCutCellRule<T>(baseRule, *levelSet));
    }

    template<class T>
    static typename gsQuadRule<T>::uPtr
    makeOctreePtr(const gsDomain<T> & domain,
                  const gsOptionList & options,
                  const Real quA,
                  const index_t quB,
                  short_t fixDir)
    {
        gsGaussRule<T> baseRule(numNodes(domain, quA, quB, fixDir));

        const gsFunction<T> * levelSet = nullptr;
        if (const auto * d1 = dynamic_cast<const gsImplicitTrimmedDomain<1,T>*>(&domain))
            levelSet = &d1->implicitFunction();
        else if (const auto * d2 = dynamic_cast<const gsImplicitTrimmedDomain<2,T>*>(&domain))
            levelSet = &d2->implicitFunction();
        else if (const auto * d3 = dynamic_cast<const gsImplicitTrimmedDomain<3,T>*>(&domain))
            levelSet = &d3->implicitFunction();

        if (!levelSet)
        {
            static bool warned = false;
            if (!warned)
            {
                gsWarn << "OctreeRule requested on a non-implicit domain; "
                       << "falling back to GaussRule for this integration context.\n";
                warned = true;
            }
            return typename gsQuadRule<T>::uPtr(new gsGaussRule<T>(baseRule));
        }

        const index_t levels = options.askInt("octLevels", 0);
        const index_t quDim  = options.askInt("quDim", -1);
        if (quDim >= 0) // surface quadrature on the interface phi == 0
            return typename gsQuadRule<T>::uPtr(
                new gsOctreeCutCellSurfaceRule<T>(*levelSet,
                    numNodes(domain,quA,quB,fixDir).maxCoeff(), levels));

        return typename gsQuadRule<T>::uPtr(new gsOctreeCutCellRule<T>(baseRule, *levelSet, levels));
    }

};

/**
    \brief Quadrature wrapper that suppresses contributions where a level-set is negative.

    The wrapped quadrature rule is first mapped to the current element and then
    evaluated with the provided level-set. Nodes with negative level-set value
    receive zero weight.

    \ingroup Assembler
*/
template<class T>
class gsCutCellRule GISMO_FINAL : public gsQuadRule<T>
{
public:
    typedef memory::unique_ptr<gsCutCellRule> uPtr;

    gsCutCellRule(const gsQuadRule<T> & rule, const gsFunction<T> & levelSet)
        : m_rule(rule), m_levelSet(levelSet)
    {
    }

    static uPtr make(const gsQuadRule<T> & rule, const gsFunction<T> & levelSet)
    {
        return uPtr(new gsCutCellRule(rule, levelSet));
    }

    using gsQuadRule<T>::mapTo;
    void mapTo(const gsVector<T>& lower, const gsVector<T>& upper,
               gsMatrix<T>& nodes, gsVector<T>& weights) const override
    {
        gsMatrix<T> vals;
        m_rule.mapTo(lower, upper, nodes, weights);
        m_levelSet.eval_into(nodes, vals);
        for (index_t i = 0; i < vals.cols(); ++i)
            if (vals(0, i) > 0)
                weights[i] = 0;
    }

private:
    gsQuadRule<T> m_rule;
    const gsFunction<T> & m_levelSet;
};

template<class T>
class gsOctreeCutCellRule GISMO_FINAL : public gsQuadRule<T>
{
public:
    typedef memory::unique_ptr<gsOctreeCutCellRule> uPtr;

    gsOctreeCutCellRule(const gsQuadRule<T> & rule,
                        const gsFunction<T> & levelSet,
                        index_t levels)
        : m_rule(rule), m_levelSet(levelSet), m_levels(levels < 0 ? 0 : levels)
    {
    }

    using gsQuadRule<T>::mapTo;

    void mapTo(const gsVector<T>& lower, const gsVector<T>& upper,
               gsMatrix<T>& nodes, gsVector<T>& weights) const override
    {
        nodes.resize(lower.size(), 0);
        weights.resize(0);
        _mapBox(lower, upper, m_levels, nodes, weights);
    }

private:
    void _append(gsMatrix<T>& nodes, gsVector<T>& weights,
                 const gsMatrix<T>& addNodes, const gsVector<T>& addWeights) const
    {
        if (addNodes.cols() == 0)
            return;

        const index_t oldCols = nodes.cols();
        nodes.conservativeResize(nodes.rows(), oldCols + addNodes.cols());
        nodes.block(0, oldCols, nodes.rows(), addNodes.cols()) = addNodes;

        const index_t oldSize = weights.size();
        weights.conservativeResize(oldSize + addWeights.size());
        weights.segment(oldSize, addWeights.size()) = addWeights;
    }

    void _mapBox(const gsVector<T>& lower, const gsVector<T>& upper,
                 index_t level, gsMatrix<T>& nodes, gsVector<T>& weights) const
    {
        const index_t dim = lower.size();
        const index_t numCorners = (index_t(1) << dim);

        gsMatrix<T> samplePoints(dim, numCorners + 1);
        for (index_t mask = 0; mask < numCorners; ++mask)
        {
            for (index_t j = 0; j < dim; ++j)
                samplePoints(j, mask) = (mask & (index_t(1) << j)) ? upper[j] : lower[j];
        }
        samplePoints.col(numCorners) = (lower + upper) / T(2);

        gsMatrix<T> values;
        m_levelSet.eval_into(samplePoints, values);

        bool allInside = true;
        bool allOutside = true;
        for (index_t i = 0; i < values.cols(); ++i)
        {
            const T v = values(0, i);
            if (v > T(0))
                allInside = false;
            if (v < T(0))
                allOutside = false;
        }

        if (allOutside)
            return;

        if (allInside || level == 0)
        {
            gsMatrix<T> leafNodes;
            gsVector<T> leafWeights;
            m_rule.mapTo(lower, upper, leafNodes, leafWeights);

            if (!allInside && leafNodes.cols() != 0)
            {
                gsMatrix<T> leafValues;
                m_levelSet.eval_into(leafNodes, leafValues);
                for (index_t i = 0; i < leafValues.cols(); ++i)
                    if (leafValues(0, i) > T(0))
                        leafWeights[i] = T(0);
            }

            _append(nodes, weights, leafNodes, leafWeights);
            return;
        }

        const gsVector<T> mid = (lower + upper) / T(2);
        const index_t childCount = (index_t(1) << dim);
        for (index_t mask = 0; mask < childCount; ++mask)
        {
            gsVector<T> childLower = lower;
            gsVector<T> childUpper = upper;
            for (index_t j = 0; j < dim; ++j)
            {
                if (mask & (index_t(1) << j))
                    childLower[j] = mid[j];
                else
                    childUpper[j] = mid[j];
            }
            _mapBox(childLower, childUpper, level - 1, nodes, weights);
        }
    }

    gsQuadRule<T> m_rule;
    const gsFunction<T> & m_levelSet;
    index_t m_levels;
};


// --- Helpers for contour-based surface quadrature (2D, marching squares) -----

/// \brief Extracts piecewise-linear segments of the interface {phi==0} inside
/// the axis-aligned 2D box [lower,upper], appending endpoint pairs to \a segs.
/// Linear approximation per cell (marching squares); the saddle case is
/// resolved with the centre sample. \ingroup Assembler
template<class T>
void gsExtractInterfaceSegments(const gsVector<T>           & lower,
                                const gsVector<T>           & upper,
                                const gsFunction<T>         & levelSet,
                                std::vector<std::pair<gsVector<T>,gsVector<T> > > & segs)
{
    GISMO_ASSERT(lower.size()==2, "Contour surface quadrature is 2D-only.");

    // Corners in CCW order: c0=(x0,y0) c1=(x1,y0) c2=(x1,y1) c3=(x0,y1)
    gsMatrix<T> C(2,4);
    C.col(0) << lower[0], lower[1];
    C.col(1) << upper[0], lower[1];
    C.col(2) << upper[0], upper[1];
    C.col(3) << lower[0], upper[1];

    gsMatrix<T> vals;
    levelSet.eval_into(C, vals);

    // Edges (corner index pairs) in CCW order
    static const int E[4][2] = {{0,1},{1,2},{2,3},{3,0}};
    std::vector<gsVector<T> > x;
    for (int e=0; e!=4; ++e)
    {
        const T pa = vals(0,E[e][0]);
        const T pb = vals(0,E[e][1]);
        if ( (pa<T(0)) != (pb<T(0)) )           // sign change along the edge
        {
            const T denom = pa - pb;
            const T t = (denom!=T(0)) ? pa/denom : T(0.5);
            x.push_back(C.col(E[e][0]) + t*(C.col(E[e][1])-C.col(E[e][0])));
        }
    }

    if (x.size()==2)
        segs.push_back(std::make_pair(x[0], x[1]));
    else if (x.size()==4)                        // saddle: disambiguate by centre
    {
        gsMatrix<T> ctr(2,1);
        ctr << (lower[0]+upper[0])/T(2), (lower[1]+upper[1])/T(2);
        gsMatrix<T> cval;
        levelSet.eval_into(ctr, cval);
        if (cval(0,0) < T(0)) // centre inside
        {
            segs.push_back(std::make_pair(x[0], x[3]));
            segs.push_back(std::make_pair(x[1], x[2]));
        }
        else
        {
            segs.push_back(std::make_pair(x[0], x[1]));
            segs.push_back(std::make_pair(x[2], x[3]));
        }
    }
    // 0 crossings: interface does not cross this cell -> nothing appended
}

/// \brief Places \a n1d-point 1D Gauss nodes on each segment with genuine
/// arc-length weights (parametric). \ingroup Assembler
template<class T>
void gsPlaceGaussOnSegments(const std::vector<std::pair<gsVector<T>,gsVector<T> > > & segs,
                            index_t n1d,
                            gsMatrix<T> & nodes,
                            gsVector<T> & weights)
{
    gsGaussRule<T> g1d(gsVector<index_t>::Constant(1, math::max(index_t(1),n1d)));
    gsMatrix<T> gn;
    gsVector<T> gw;
    g1d.mapTo(gsVector<T>::Zero(1), gsVector<T>::Constant(1,T(1)), gn, gw); // on [0,1]

    const index_t k = gw.size();
    nodes.resize(2, k*static_cast<index_t>(segs.size()));
    weights.resize(k*static_cast<index_t>(segs.size()));

    index_t c = 0;
    for (size_t s=0; s!=segs.size(); ++s)
    {
        const gsVector<T> d = segs[s].second - segs[s].first;
        const T len = d.norm();
        for (index_t j=0; j!=k; ++j, ++c)
        {
            nodes.col(c) = segs[s].first + gn(0,j)*d;
            weights[c]   = gw[j]*len;          // arc-length weight
        }
    }
}

/**
    \brief Cut-cell surface quadrature: linear contour of {phi==0} per cell with
    1D Gauss points and arc-length weights. Counterpart to gsCutCellRule (volume).
    \ingroup Assembler
*/
template<class T>
class gsCutCellSurfaceRule GISMO_FINAL : public gsQuadRule<T>
{
public:
    gsCutCellSurfaceRule(const gsFunction<T> & levelSet, index_t n1d)
    : m_levelSet(levelSet), m_n1d(n1d) { }

    using gsQuadRule<T>::mapTo;
    void mapTo(const gsVector<T>& lower, const gsVector<T>& upper,
               gsMatrix<T>& nodes, gsVector<T>& weights) const override
    {
        nodes.resize(lower.size(), 0);
        weights.resize(0);
        if (lower.size()!=2) return;           // contour surface rule is 2D-only

        std::vector<std::pair<gsVector<T>,gsVector<T> > > segs;
        gsExtractInterfaceSegments(lower, upper, m_levelSet, segs);
        if (segs.empty()) return;
        gsPlaceGaussOnSegments(segs, m_n1d, nodes, weights);
    }

private:
    const gsFunction<T> & m_levelSet;
    index_t m_n1d;
};

/**
    \brief Octree-refined cut-cell surface quadrature: subdivides straddling
    cells \a levels times, then a linear contour + 1D Gauss per finest leaf, so
    the piecewise-linear interface better resolves curved boundaries.
    \ingroup Assembler
*/
template<class T>
class gsOctreeCutCellSurfaceRule GISMO_FINAL : public gsQuadRule<T>
{
public:
    gsOctreeCutCellSurfaceRule(const gsFunction<T> & levelSet, index_t n1d, index_t levels)
    : m_levelSet(levelSet), m_n1d(n1d), m_levels(levels<0?0:levels) { }

    using gsQuadRule<T>::mapTo;
    void mapTo(const gsVector<T>& lower, const gsVector<T>& upper,
               gsMatrix<T>& nodes, gsVector<T>& weights) const override
    {
        nodes.resize(lower.size(), 0);
        weights.resize(0);
        if (lower.size()!=2) return;

        std::vector<std::pair<gsVector<T>,gsVector<T> > > segs;
        _collect(lower, upper, m_levels, segs);
        if (segs.empty()) return;
        gsPlaceGaussOnSegments(segs, m_n1d, nodes, weights);
    }

private:
    void _collect(const gsVector<T>& lower, const gsVector<T>& upper, index_t level,
                  std::vector<std::pair<gsVector<T>,gsVector<T> > > & segs) const
    {
        // Sample 4 corners + centre to test whether the interface crosses
        gsMatrix<T> S(2,5);
        S.col(0) << lower[0], lower[1];
        S.col(1) << upper[0], lower[1];
        S.col(2) << upper[0], upper[1];
        S.col(3) << lower[0], upper[1];
        S.col(4) << (lower[0]+upper[0])/T(2), (lower[1]+upper[1])/T(2);

        gsMatrix<T> V;
        m_levelSet.eval_into(S, V);
        bool allIn = true, allOut = true;
        for (index_t i=0; i!=V.cols(); ++i)
        {
            if (V(0,i) > T(0)) allIn  = false;
            if (V(0,i) < T(0)) allOut = false;
        }
        if (allIn || allOut) return;           // interface not in this box

        if (level==0)
        {
            gsExtractInterfaceSegments(lower, upper, m_levelSet, segs);
            return;
        }

        const gsVector<T> mid = (lower + upper) / T(2);
        for (index_t mask=0; mask<4; ++mask)
        {
            gsVector<T> cl = lower, cu = upper;
            for (index_t j=0; j!=2; ++j)
            {
                if (mask & (index_t(1) << j)) cl[j] = mid[j];
                else                          cu[j] = mid[j];
            }
            _collect(cl, cu, level-1, segs);
        }
    }

    const gsFunction<T> & m_levelSet;
    index_t m_n1d;
    index_t m_levels;
};


}// namespace gismo
