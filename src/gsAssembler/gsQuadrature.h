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
#include <gsAssembler/gsMomentRule.h>
#include <gsAssembler/gsGaussRule.h>

#include <gsCore/gsFunction.h>
#include <gsDomain/gsDomainIterator.h>
#include <gsDomain/gsImplicitTrimmedDomain.h>

#include <cmath>   // std::lround, used by makeMomentFittingPtr below

#ifdef gsAlgoim_ENABLED
#include <gsAlgoim/gsAlgoimRule.h>
// Retained: the DEFAULT rule wrapped by MomentFittingRule is the adaptive one.
// gsAlgoimRule.h only offers gsAlgoimGenericRule, a different rule with
// different cut weights, so dropping this include would move the numbers.
#include <gsAlgoim/gsAlgoimAdaptiveRule.h>
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
        OctreeRule    = 13,
        MomentFittingRule = 14 ///< Solve-free nodal moment fitting (core, see
                               ///< gsMomentRule) around the cut-cell rule
                               ///< selected by the option "quMomentUnderlying"
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
        const typename gsDomain<T>::Ptr dom = basis.domain();
        return get<T>(*dom,options,fixDir,_basisDegrees(basis,*dom));
    }

    /// Constructs a quadrature rule based on input \a options
    /// \param[in] degrees Per-direction polynomial degrees of the space being
    ///        integrated; empty means fall back to \c gsDomain::degree(i),
    ///        which is a guess and a deprecated path -- supply the degrees
    ///        from the basis whenever you have one (see the KNOWN DESIGN
    ///        DEFECT note at gsDomain.h:232-275).
    /// \pre \a degrees is empty or has exactly domain.dim() non-negative entries.
    ///      Checked by GISMO_ASSERT only: in release builds a malformed vector
    ///      silently yields wrong quadrature, so callers must size it correctly.
    template<class T>
    static gsQuadRule<T> get(const gsDomain<T> & domain,
                             const gsOptionList & options, short_t fixDir = -1,
                             const gsVector<short_t> & degrees = gsVector<short_t>())
    {
        const index_t qu  = options.askInt("quRule", GaussLegendre);
        const Real    quA = options.getReal("quA");
        const index_t quB = options.getInt ("quB");
        const gsVector<index_t> nnodes = numNodes(domain,quA,quB,fixDir,degrees);
        return get<T>(qu, nnodes);
    }

    /// Constructs a quadrature rule based on input \a options
    template<class T>
    static typename gsQuadRule<T>::uPtr
                      getPtr(const gsBasis<T> & basis,
                             const gsOptionList & options, short_t fixDir = -1)
    {
        const typename gsDomain<T>::Ptr dom = basis.domain();
        return getPtr<T>(*dom,options,fixDir,_basisDegrees(basis,*dom));
    }

    /// Constructs a quadrature rule based on input \a options
    /// \param[in] degrees Per-direction polynomial degrees of the space being
    ///        integrated; empty means fall back to \c gsDomain::degree(i),
    ///        which is a guess and a deprecated path -- supply the degrees
    ///        from the basis whenever you have one (see the KNOWN DESIGN
    ///        DEFECT note at gsDomain.h:232-275).
    /// \pre \a degrees is empty or has exactly domain.dim() non-negative entries.
    ///      Checked by GISMO_ASSERT only: in release builds a malformed vector
    ///      silently yields wrong quadrature, so callers must size it correctly.
    template<class T>
    static typename gsQuadRule<T>::uPtr
                      getPtr(const gsDomain<T> & domain,
                             const gsOptionList & options, short_t fixDir = -1,
                             const gsVector<short_t> & degrees = gsVector<short_t>())
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
                        return gsGaussRule<T>::make(numNodes(domain,quA,quB,fixDir,degrees));
                    case GaussLobatto :
                        return gsLobattoRule<T>::make(numNodes(domain,quA,quB,fixDir,degrees));
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

                const gsVector<index_t> nnodesI = numNodes(domain,quA,quB,fixDir,degrees);
                const gsVector<index_t> nnodesB = numNodes(domain,quAb,quBb,fixDir,degrees);

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
            // degrees is deliberately NOT used here: gsPatchRule::make takes
            // the target-space order directly as quA and never consults
            // domain.degree() (or the degrees vector) for it.
            return gsPatchRule<T>::make(domain,cast<T,index_t>(quA),quB,over,fixDir);
        }
        else if (qu==CutCellRule)
        {
            // quDim: -1 = volumetric (mask phi>0); >=0 = surface (integrate phi==0)
            const index_t quDim = options.askInt("quDim", -1);
            return makeCutCellPtr<T>(domain, quA, quB, fixDir, quDim, degrees);
        }
        else if (qu==OctreeRule)
        {
            return makeOctreePtr<T>(domain, options, quA, quB, fixDir, degrees);
        }
        else if (qu==MomentFittingRule)
        {
            // Core rule (gsMomentRule); the wrapped rule is selected by the
            // option "quMomentUnderlying", see makeMomentFittingPtr.
            return makeMomentFittingPtr<T>(domain, options, quA, quB, fixDir, degrees);
        }
#ifdef gsAlgoim_ENABLED
        else if (qu==AlgoimRule)
        {
            // Forward the option list so the immersed boundary assembly can
            // request surface quadrature (dim == D) via the "quDim" option.
            return gsAlgoimGenericRule<T>::make(domain, options, degrees);
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
    /// of \a basis
    template<class T>
    static gsVector<index_t> numNodes(const gsBasis<T> & basis,
                               const Real quA, const index_t quB, short_t fixDir = -1)
    {
        const typename gsDomain<T>::Ptr dom = basis.domain();
        return numNodes(*dom,quA,quB,fixDir,_basisDegrees(basis,*dom));
    }

    /// Computes and integer quA*deg_i + quB where deg_i is either the
    /// per-direction \a degrees, or -- when \a degrees is empty (the
    /// deprecated fallback) -- the guessed degree of \a domain.
    /// \param[in] degrees Per-direction polynomial degrees of the space being
    ///        integrated; empty means fall back to \c gsDomain::degree(i),
    ///        which is a guess and a deprecated path -- supply the degrees
    ///        from the basis whenever you have one (see the KNOWN DESIGN
    ///        DEFECT note at gsDomain.h:232-275).
    /// \pre \a degrees is empty or has exactly domain.dim() non-negative entries.
    ///      Checked by GISMO_ASSERT only: in release builds a malformed vector
    ///      silently yields wrong quadrature, so callers must size it correctly.
    template<class T>
    static gsVector<index_t> numNodes(const gsDomain<T> & domain,
                               const Real quA, const index_t quB, short_t fixDir = -1,
                               const gsVector<short_t> & degrees = gsVector<short_t>())
    {
        const short_t d  = domain.dim();
        GISMO_ASSERT( fixDir < d && fixDir>-2, "Invalid input fixDir = "<<fixDir);
        GISMO_ASSERT(0==degrees.size() || degrees.size()==d,
                     "gsQuadrature: degrees has size "<<degrees.size()
                     <<", expected 0 (domain fallback) or "<<d);
        GISMO_ASSERT(0==degrees.size() || degrees.minCoeff()>=0,
                     "gsQuadrature: negative entry in degrees.");
        gsVector<index_t> nnodes(d);

        if (-1==fixDir)
            fixDir = d;
        else
            nnodes[fixDir] = 1;

        short_t i;
        for(i=0; i!=fixDir; ++i )
            //note: +0.5 for rounding
            nnodes[i] = cast<Real,index_t>(quA * (0==degrees.size() ? domain.degree(i) : degrees[i]) + quB + 0.5);
        for(++i; i<d; ++i )
            nnodes[i] = cast<Real,index_t>(quA * (0==degrees.size() ? domain.degree(i) : degrees[i]) + quB + 0.5);
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
     * @param[in] degrees Per-direction polynomial degrees of the space being
     *        integrated; empty means fall back to \c gsDomain::degree(i),
     *        which is a guess and a deprecated path -- supply the degrees
     *        from the basis whenever you have one (see the KNOWN DESIGN
     *        DEFECT note at gsDomain.h:232-275).
     * \pre \a degrees is empty or has exactly domain.dim() non-negative entries.
     *      Checked by GISMO_ASSERT only: in release builds a malformed vector
     *      silently yields wrong quadrature, so callers must size it correctly.
     * @return            A matrix where each column represents a quadrature node in the parametric domain.
     */
    template<class T>
    static gsMatrix<T> getAllNodes(const gsDomain<T> & domain,
                                   const gsOptionList & options,
                                   const gsVector<short_t> & degrees = gsVector<short_t>())
    {
        typename gsBasis<T>::domainIter domIt    = domain.beginAll();
        typename gsBasis<    T>::    domainIter domItEnd = domain.endAll();

        index_t     quadSize = 0;
        typename gsQuadRule<T>::uPtr QuRule;
        QuRule = getPtr(domain, options, -1, degrees);

        for (; domIt<domItEnd; ++domIt )
        {
            QuRule = gsQuadrature::getPtr(domain, options, -1, degrees);
            quadSize+=QuRule->numNodes();
        }

        gsMatrix<T> result(domain.dim(),quadSize);

        index_t offset = 0;
        gsMatrix<T> nodes;
        gsVector<T> weights;

        domIt = domain.beginAll();
        for (; domIt<domItEnd; ++domIt )
        {
            QuRule = gsQuadrature::getPtr(domain, options, -1, degrees);
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
        const typename gsDomain<T>::Ptr dom = basis.domain();
        return getAllNodes(*dom,options,_basisDegrees(basis,*dom));
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
        const typename gsDomain<T>::Ptr dom = basis.domain();
        return getAllNodes(*dom,options,side,_basisDegrees(basis,*dom));
    }

    /**
     * @brief Get all quadrature nodes for a specified side of a given domain.
     *
     * @tparam T        Real type.
     * @param[in] domain    The domain for which the quadrature nodes are to be collected.
     * @param[in] options   Quadrature rule.
     * @param[in] side      The side of the domain.
     * @param[in] degrees Per-direction polynomial degrees of the space being
     *        integrated; empty means fall back to \c gsDomain::degree(i),
     *        which is a guess and a deprecated path -- supply the degrees
     *        from the basis whenever you have one (see the KNOWN DESIGN
     *        DEFECT note at gsDomain.h:232-275).
     * \pre \a degrees is empty or has exactly domain.dim() non-negative entries.
     *      Checked by GISMO_ASSERT only: in release builds a malformed vector
     *      silently yields wrong quadrature, so callers must size it correctly.
     * @return result   A matrix of quadrature nodes, where each column corresponds to a quadrature node.
     */
    template<class T>
    static gsMatrix<T> getAllNodes( const gsDomain<T> & domain,
                                    const gsOptionList & options,
                                    const patchSide side,
                                    const gsVector<short_t> & degrees = gsVector<short_t>())
    {
        typename gsBasis<T>::domainIter domIt    = domain.beginBdr(side.side());
        typename gsBasis<T>::domainIter domItEnd = domain.endBdr(side.side());

        index_t quadSize = 0;
        typename gsQuadRule<T>::uPtr QuRule;
        QuRule = getPtr(domain, options, side.side().direction(), degrees);

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
            QuRule = gsQuadrature::getPtr(domain, options, side.side().direction(), degrees);
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
        const typename gsDomain<T>::Ptr dom = basis.domain();
        return getAllNodes(*dom,options,sides,_basisDegrees(basis,*dom));
    }

    /**
     * @brief Collects all quadrature nodes for a multi-basis.
     *
     * @param[in] degrees Per-direction polynomial degrees of the space being
     *        integrated; empty means fall back to \c gsDomain::degree(i),
     *        which is a guess and a deprecated path -- supply the degrees
     *        from the basis whenever you have one (see the KNOWN DESIGN
     *        DEFECT note at gsDomain.h:232-275).
     * \pre \a degrees is empty or has exactly domain.dim() non-negative entries.
     *      Checked by GISMO_ASSERT only: in release builds a malformed vector
     *      silently yields wrong quadrature, so callers must size it correctly.
     */
    template<class T>
    static gsMatrix<T> getAllNodes( const gsDomain<T> & domain,
                                    const gsOptionList & options,
                                    const std::vector<patchSide> & sides,
                                    const gsVector<short_t> & degrees = gsVector<short_t>())
    {
        std::vector<gsMatrix<T>> nodes(sides.size());
        index_t cols = 0;
        for (size_t s = 0; s != sides.size(); s++)
        {
            nodes[s] = getAllNodes(domain,options,sides[s],degrees);
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

    /// \brief One-shot flag for the moment-fitting exactness warning.
    /// The warning fires once per binary image (this function is
    /// implicitly inline, so a shared library and an executable linking
    /// against it each get their own copy of the static) so a long assembly
    /// is not drowned in repetitions, and the flag is shared by every call
    /// to this function within that image; reset it to false (test hook) to
    /// re-arm it deterministically.
    static bool & momentExactnessWarned()
    { static bool warned = false; return warned; }

private:
    /// Per-direction degrees of \a basis for integration over \a domain, as
    /// required by the quadrature factories.  This is the ACCURATE source of
    /// the degree; the gsDomain overloads' fallback to gsDomain::degree() is
    /// a guess.
    ///
    /// The result is sized from \a domain -- the same object the consumers
    /// (numNodes, _maxDegree) check it against -- so it can never be a
    /// wrong-length vector.  If the basis's parametric dimension disagrees
    /// with domain.dim(), an EMPTY vector is returned and the caller keeps
    /// the (deprecated, guessing) gsDomain::degree() fallback: a short or
    /// long vector would be read out of bounds in Release builds, where the
    /// GISMO_ASSERT at numNodes is compiled out. For every basis family in
    /// the tree the two dimensions agree, so this branch is defensive only
    /// (see also gsExprHelper<T>::quadratureDegrees, which takes the same
    /// decision for the same reason).
    template<class T>
    static gsVector<short_t> _basisDegrees(const gsBasis<T> & basis,
                                           const gsDomain<T> & domain)
    {
        const short_t d = basis.domainDim();
        if (domain.dim() != d)
            return gsVector<short_t>();
        gsVector<short_t> degs(d);
        for (short_t i = 0; i != d; ++i)
            degs[i] = basis.degree(i);
        return degs;
    }

    /// Largest degree relevant for \a domain: max over \a degrees, or the
    /// domain's guessed degree when \a degrees is empty (deprecated path).
    /// The max is the conservative choice for exactness: with n_i = quA*deg_i+quB,
    /// n_i-(2 deg_i+1) = (quA-2) deg_i + quB - 1, so for quA < 2 the largest
    /// degree is the binding direction.
    template<class T>
    static short_t _maxDegree(const gsDomain<T> & domain,
                              const gsVector<short_t> & degrees)
    {
        GISMO_ASSERT(0==degrees.size() || degrees.size()==domain.dim(),
                     "gsQuadrature: degrees has size "<<degrees.size()
                     <<", expected 0 (domain fallback) or "<<domain.dim());
        GISMO_ASSERT(0==degrees.size() || degrees.minCoeff()>=0,
                     "gsQuadrature: negative entry in degrees.");
        return 0==degrees.size() ? domain.degree() : degrees.maxCoeff();
    }

    template<class T>
    static typename gsQuadRule<T>::uPtr
    makeCutCellPtr(const gsDomain<T> & domain,
                   const Real quA,
                   const index_t quB,
                   short_t fixDir,
                   index_t quDim = -1,
                   const gsVector<short_t> & degrees = gsVector<short_t>())
    {
        const gsVector<index_t> nnodes = numNodes(domain,quA,quB,fixDir,degrees);

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
            return gsGaussRule<T>::make(nnodes);
        }

        if (quDim >= 0) // surface quadrature on the interface phi == 0
            return typename gsQuadRule<T>::uPtr(
                new gsCutCellSurfaceRule<T>(*levelSet, nnodes.maxCoeff()));

        return typename gsQuadRule<T>::uPtr(
            new gsCutCellRule<T>(gsGaussRule<T>::make(nnodes), *levelSet));
    }

    template<class T>
    static typename gsQuadRule<T>::uPtr
    makeOctreePtr(const gsDomain<T> & domain,
                  const gsOptionList & options,
                  const Real quA,
                  const index_t quB,
                  short_t fixDir,
                  const gsVector<short_t> & degrees = gsVector<short_t>())
    {
        const gsVector<index_t> nnodes = numNodes(domain, quA, quB, fixDir, degrees);

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
            return gsGaussRule<T>::make(nnodes);
        }

        const index_t levels = options.askInt("octLevels", 0);
        const index_t quDim  = options.askInt("quDim", -1);
        if (quDim >= 0) // surface quadrature on the interface phi == 0
            return typename gsQuadRule<T>::uPtr(
                new gsOctreeCutCellSurfaceRule<T>(*levelSet,
                    nnodes.maxCoeff(), levels));

        return typename gsQuadRule<T>::uPtr(
            new gsOctreeCutCellRule<T>(gsGaussRule<T>::make(nnodes), *levelSet, levels));
    }

    /// \brief Solve-free nodal moment fitting (gsMomentRule) around a cut-cell
    /// rule built on the implicit function of \a domain.
    ///
    /// The wrapped rule is selected by the assembler-level option
    /// "quMomentUnderlying" (a gsQuadrature::rule value; 0 = use the default,
    /// which is the adaptive Algoim rule where gsAlgoim is available and
    /// CutCellRule otherwise). quA and quB size the OUTPUT tensor grid exactly
    /// as before; maxDepth, indicator, indicatorTol and LipschitzConstant are
    /// forwarded to the adaptive Algoim rule when that one is wrapped.
    ///
    /// Surface quadrature (quDim >= 0) and PatchRule are refused, see below.
    template<class T>
    static typename gsQuadRule<T>::uPtr
    makeMomentFittingPtr(const gsDomain<T> & domain,
                         const gsOptionList & options,
                         const Real quA,
                         const index_t quB,
                         short_t fixDir,
                         const gsVector<short_t> & degrees = gsVector<short_t>())
    {
        // quDim: -1 = volumetric (phi<0); >=0 = surface (phi==0), see gsAlgoimRule
        const index_t quDim = options.askInt("quDim", -1);
        GISMO_ENSURE(quDim < 0,
            "MomentFittingRule cannot wrap a surface rule (quDim = " << quDim << "). "
            "Surface points lie on the zero level set, a manifold of lower dimension "
            "than the element; compressing them onto a volume tensor grid is "
            "meaningless. Use quRule = CutCellRule/OctreeRule/AlgoimRule with "
            "quDim >= 0 directly for interface integrals.");

        // Which rule is compressed? 0 means "use the default".
        index_t inner = options.askInt("quMomentUnderlying", 0);
        if (0 == inner)
#ifdef gsAlgoim_ENABLED
            inner = AlgoimRule;      // -> gsAlgoimAdaptiveRule, see below
#else
            inner = CutCellRule;
#endif

        GISMO_ENSURE(inner != PatchRule,
            "MomentFittingRule cannot wrap PatchRule: gsPatchRule is patch-global "
            "(gsPatchRule.h:258 m_maps, :249-250 its own m_nodes/m_weights), so an "
            "element's output is a slice of a global rule, not an independent local "
            "rule. It already uses FEWER points than Gauss-per-element, so "
            "compressing it onto a (2p+1)^d grid would increase the point count and "
            "discard the global optimality that is the rule's whole purpose.");

        // The legacy "alpha" key is the FCM fictitious-domain blend
        // w <- (1-a)w + a*w^G, NOT the multiplicative integrand field alpha of
        // gsMomentRule (see its file header). It is not supported on this path;
        // dropping it silently would change results without notice.
        if (options.askReal("alpha", 0.0) > 0)
        {
            static bool warned = false;
            if (!warned)
            {
                gsWarn << "MomentFittingRule: the option alpha = "
                       << options.askReal("alpha", 0.0)
                       << " is the FCM fictitious-domain blend, which this path "
                       << "does not support (gsMomentRule's alpha is a different "
                       << "quantity: a multiplicative field inside the integrand). "
                       << "It is IGNORED here.\n";
                warned = true;
            }
        }

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
                gsWarn << "MomentFittingRule requested on a non-implicit domain; "
                       << "falling back to GaussRule for this integration context.\n";
                warned = true;
            }
            // Unwrapped on purpose: compressing a Gauss rule onto a Gauss grid
            // of the same element is a no-op at best.
            return gsGaussRule<T>::make(numNodes(domain, quA, quB, fixDir, degrees));
        }

        // Forwarded quadrature counts; getPtr() always supplies them, the
        // arguments are the fallback for direct callers.
        const Real    quA_used = options.askReal("quA", quA);
        const index_t quB_used = options.askInt ("quB", quB);

        // --- the wrapped (compressed) rule ---------------------------------
        typename gsQuadRule<T>::uPtr innerRule;
        if (CutCellRule == inner)
            innerRule = makeCutCellPtr<T>(domain, quA, quB, fixDir, -1, degrees);
        else if (OctreeRule == inner)
            innerRule = makeOctreePtr<T>(domain, options, quA, quB, fixDir, degrees);
        else if (AlgoimRule == inner)
        {
#ifdef gsAlgoim_ENABLED
            gsOptionList o = gsAlgoimAdaptiveRule<T>::defaultOptions();
            o.setInt   ("dim", -1);        // volumetric; surface refused above
            o.setReal  ("quA", quA_used);
            o.setInt   ("quB", quB_used);
            // Adaptive-rule controls (kept at their defaults if unset)
            o.setInt   ("maxDepth",          options.askInt   ("maxDepth",o.getInt   ("maxDepth")));
            o.setString("indicator",         options.askString("indicator",o.getString("indicator")));
            o.setReal  ("indicatorTol",      options.askReal  ("indicatorTol",
                                                               o.getReal("indicatorTol")));
            o.setReal  ("LipschitzConstant", options.askReal  ("LipschitzConstant",
                                                               o.getReal("LipschitzConstant")));
            innerRule = typename gsQuadRule<T>::uPtr(
                new gsAlgoimAdaptiveRule<T>(*levelSet, _maxDegree(domain, degrees), o));
#else
            GISMO_ERROR("MomentFittingRule: quMomentUnderlying = AlgoimRule ("
                        << AlgoimRule << ") requires the optional gsAlgoim module, "
                        "which is not enabled in this build (gsAlgoim_ENABLED "
                        "undefined). Use CutCellRule or OctreeRule instead.");
#endif
        }
        else
            GISMO_ERROR("MomentFittingRule cannot wrap quadrature rule "<<inner<<".");

        // --- output grid ----------------------------------------------------
        // Exactness warning on the assembler-facing path. gsMomentRule itself
        // does not warn -- it is a geometry-free compressor and has no notion of
        // the discretization degree -- so this is the only place where the
        // caller-controlled order can be checked against the operator
        // requirement n >= 2*degree+1 (see the gsMomentRule file header). The
        // one-shot flag lives behind momentExactnessWarned() and is shared by
        // every call to this function, which is exactly why the accessor
        // exists -- it keeps a long assembly from drowning in repetitions,
        // and can be reset by tests to re-arm it deterministically.
        const short_t deg = _maxDegree(domain, degrees);
        const long nRaw = std::lround(static_cast<double>(quA_used)
                                      * static_cast<double>(deg))
                        + static_cast<long>(quB_used);
        const index_t n = nRaw > 0 ? static_cast<index_t>(nRaw) : 1;
        if (n < 2 * static_cast<index_t>(deg) + 1)
        {
            bool & warned = momentExactnessWarned();
            if (!warned)
            {
                gsWarn << "MomentFittingRule: output order n = " << n
                       << " points/direction, but degree " << deg
                       << " operators need n >= 2*degree+1 = "
                       << (2 * static_cast<index_t>(deg) + 1)
                       << "; moment fitting is exact only up to degree n-1. "
                       << "Fine for measure/volume integrands, under-integrating "
                       << "for operator assembly -- remedy: quA = 2.0.\n";
                warned = true;
            }
        }

        // alpha is deliberately nullptr here: see the FCM note above.
        return gsMomentRule<T>::make(give(innerRule),
                                     gsVector<index_t>::Constant(domain.dim(), n));
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

    /// Takes OWNERSHIP of \a rule: the gsQuadRule hierarchy has no clone(), so a
    /// polymorphic rule stored by value would be sliced down to the base affine map.
    gsCutCellRule(typename gsQuadRule<T>::uPtr rule, const gsFunction<T> & levelSet)
        : m_rule(give(rule)), m_levelSet(levelSet)
    {
        GISMO_ASSERT(m_rule, "gsCutCellRule: null wrapped rule.");
    }

    /// Takes ownership of \a rule, see the constructor.
    static uPtr make(typename gsQuadRule<T>::uPtr rule, const gsFunction<T> & levelSet)
    {
        return uPtr(new gsCutCellRule(give(rule), levelSet));
    }

    using gsQuadRule<T>::mapTo;
    void mapTo(const gsVector<T>& lower, const gsVector<T>& upper,
               gsMatrix<T>& nodes, gsVector<T>& weights) const override
    {
        gsMatrix<T> vals;
        m_rule->mapTo(lower, upper, nodes, weights);
        m_levelSet.eval_into(nodes, vals);
        for (index_t i = 0; i < vals.cols(); ++i)
            if (vals(0, i) > 0)
                weights[i] = 0;
    }

private:
    typename gsQuadRule<T>::uPtr m_rule;
    const gsFunction<T> & m_levelSet;
};

template<class T>
class gsOctreeCutCellRule GISMO_FINAL : public gsQuadRule<T>
{
public:
    typedef memory::unique_ptr<gsOctreeCutCellRule> uPtr;

    /// Takes OWNERSHIP of \a rule: the gsQuadRule hierarchy has no clone(), so a
    /// polymorphic rule stored by value would be sliced down to the base affine map.
    gsOctreeCutCellRule(typename gsQuadRule<T>::uPtr rule,
                        const gsFunction<T> & levelSet,
                        index_t levels)
        : m_rule(give(rule)), m_levelSet(levelSet), m_levels(levels < 0 ? 0 : levels)
    {
        GISMO_ASSERT(m_rule, "gsOctreeCutCellRule: null wrapped rule.");
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
            m_rule->mapTo(lower, upper, leafNodes, leafWeights);

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

    typename gsQuadRule<T>::uPtr m_rule;
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
