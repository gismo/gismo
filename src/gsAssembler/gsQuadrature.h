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

namespace gismo
{

/// Helper class for obtaining a quadrature rule
struct gsQuadrature
{
    typedef GISMO_COEFF_TYPE Real;
  
    /// Quadrature rule types
    enum rule
    {
        GaussLegendre = 1, ///< Gauss-Legendre quadrature
        GaussLobatto  = 2, ///< Gauss-Lobatto quadrature
        PatchRule     = 3  ///< Patch-wise quadrature rule  (Johannessen 2017)

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
        const index_t qu  = options.askInt("quRule", GaussLegendre);
        const Real    quA = options.getReal("quA");
        const index_t quB = options.getInt ("quB");
        const gsVector<index_t> nnodes = numNodes(basis,quA,quB,fixDir);
        return get<T>(qu, nnodes);
    }

    /// Constructs a quadrature rule based on input \a options
    template<class T>
    static typename gsQuadRule<T>::uPtr
                      getPtr(const gsBasis<T> & basis,
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
                        return gsGaussRule<T>::make(numNodes(basis,quA,quB,fixDir));
                    case GaussLobatto :
                        return gsLobattoRule<T>::make(numNodes(basis,quA,quB,fixDir));
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

                const gsVector<index_t> nnodesI = numNodes(basis,quA,quB,fixDir);
                const gsVector<index_t> nnodesB = numNodes(basis,quAb,quBb,fixDir);

                std::vector<gsQuadRule<T> > quInterior(nnodesI.size());
                std::vector<gsQuadRule<T> > quBoundary(nnodesB.size());

                for (index_t d = 0; d != nnodesI.size(); d++)
                {
                    quInterior[d] = getUnivariate<T>(qu,nnodesI[d]);
                    quBoundary[d] = getUnivariate<T>(qu,nnodesB[d]);
                }

                return gsOverIntegrateRule<T>::make(basis,quInterior,quBoundary);
            }
        }
        else if (qu==PatchRule)
        {
            // quA: Order of the target space
            // quB: Regularity of the target space
            return gsPatchRule<T>::make(basis,cast<T,index_t>(quA),quB,over,fixDir);
        }
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
        const short_t d  = basis.dim();
        GISMO_ASSERT( fixDir < d && fixDir>-2, "Invalid input fixDir = "<<fixDir);
        gsVector<index_t> nnodes(d);

        if (-1==fixDir)
            fixDir = d;
        else
            nnodes[fixDir] = 1;

        short_t i;
        for(i=0; i!=fixDir; ++i )
            //note: +0.5 for rounding
            nnodes[i] = cast<Real,index_t>(quA * basis.degree(i) + quB + 0.5);
        for(++i; i<d; ++i )
            nnodes[i] = cast<Real,index_t>(quA * basis.degree(i) + quB + 0.5);
        return nnodes;
    }

    // template<class T>
    // static std::pair<gsMatrix<T>,gsVector<T>> getAllNodesAndWeights(const gsBasis<T> & basis,
    //                          const gsOptionList & options)
    // {

    // }

    // HV DOCUMENT THIS
    template<class T>
    static gsMatrix<T> getAllNodes(const gsBasis<T> & basis,
                             const gsOptionList & options)
    {
        typename gsBasis<real_t>::domainIter domIt = basis.makeDomainIterator();

        index_t quadSize = 0;
        typename gsQuadRule<T>::uPtr QuRule;
        QuRule = getPtr(basis, options);

        for (; domIt->good(); domIt->next())
        {
            QuRule = gsQuadrature::getPtr(basis, options);
            quadSize+=QuRule->numNodes();
        }

        gsMatrix<T> result(basis.domainDim(),quadSize);

        index_t offset = 0;
        gsMatrix<T> nodes;
        gsVector<T> weights;
        for (domIt->reset(); domIt->good(); domIt->next() )
        {
            QuRule = gsQuadrature::getPtr(basis, options);
            // Map the Quadrature rule to the element
            QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                           nodes, weights);
            result.block(0,offset,basis.domainDim(),QuRule->numNodes()) = nodes;
            offset += QuRule->numNodes();
        }
        return result;
    }

    template<class T>
    static gsMatrix<T> getAllNodes(const gsBasis<T> & basis,
                             const gsOptionList & options, const patchSide side)
    {
        typename gsBasis<real_t>::domainIter domIt = basis.makeDomainIterator(side);

        index_t quadSize = 0;
        typename gsQuadRule<T>::uPtr QuRule;
        QuRule = getPtr(basis, options,side.side().direction());

        for (; domIt->good(); domIt->next())
        {
            QuRule = gsQuadrature::getPtr(basis, options, side.side().direction());
            quadSize+=QuRule->numNodes();
        }

        gsMatrix<T> result(basis.domainDim(),quadSize);

        index_t offset = 0;
        gsMatrix<T> nodes;
        gsVector<T> weights;
        for (domIt->reset(); domIt->good(); domIt->next() )
        {
            QuRule = gsQuadrature::getPtr(basis, options, side.side().direction());
            // Map the Quadrature rule to the element
            QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                           nodes, weights);
            result.block(0,offset,basis.domainDim(),QuRule->numNodes()) = nodes;
            offset += QuRule->numNodes();
        }
        return result;
    }

    template<class T>
    static gsMatrix<T> getAllNodes(const gsBasis<T> & basis, const gsGeometry<T> & geom,
                             const gsOptionList & options, const patchSide side)
    {
        gsMatrix<T> nodes = getAllNodes(basis,options,side);
        return geom.eval(nodes);
    }

    template<class T>
    static gsMatrix<T> getAllNodes(const gsBasis<T> & basis,
                             const gsOptionList & options, const std::vector<patchSide> sides)
    {
        std::vector<gsMatrix<T>> nodes(sides.size());
        index_t cols = 0;
        for (index_t s = 0; s != sides.size(); s++)
        {
            nodes[s] = getAllNodes(basis,options,sides[s]);
            cols += nodes[s].cols();
        }
        gsMatrix<T> result(basis.domainDim(),cols);
        cols = 0;

        for (index_t s = 0; s != sides.size(); s++)
        {
            result.block(0,cols,nodes[s].rows(),nodes[s].cols()) = nodes[s];
            cols += nodes[s].cols();
        }

        return result;
    }

    template<class T>
    static gsMatrix<T> getAllNodes(const gsBasis<T> & basis, const gsGeometry<T> & geom,
                             const gsOptionList & options, const std::vector<patchSide> sides)
    {
        gsMatrix<T> nodes = getAllNodes(basis,options,sides);
        return geom.eval(nodes);
    }

    template<class T>
    static gsMatrix<T> getAllNodes(const gsMultiBasis<T> & bases,
                             const gsOptionList & options)
    {
        return getAllNodes(bases,options);
    }

    template<class T>
    static std::vector<gsMatrix<T>> getAllNodes(const gsMultiBasis<T> & bases,
                             const gsOptionList & options, const std::vector<patchSide> sides)
    {
        GISMO_ASSERT(bases.nBases()==sides.size(),"Number of bases must be equal to the number of fixed directions");
        std::vector<gsMatrix<T>> nodes(bases.nBases());
        for (size_t p = 0; p != bases.nBases(); p++)
            nodes[p] = getAllNodes(bases.basis(p),options,sides[p]);

        return nodes;
    }

    template<class T>
    static gsMatrix<T> getAllNodes(const gsMultiBasis<T> & bases, const gsMultiPatch<T> & mp,
                             const gsOptionList & options, const std::vector<patchSide> sides)
    {
        std::vector<gsMatrix<T>> nodes = getAllNodes(bases,options,sides);
        index_t cols = 0;
        for (size_t p = 0; p != nodes.size(); p++)
            cols += nodes[p].cols();

        gsMatrix<T> result(mp.targetDim(),cols);
        cols = 0;
        for (size_t p = 0; p != nodes.size(); p++)
        {
            result.block(0,cols,mp.targetDim(),nodes[p].cols()) = mp.patch(p).eval(nodes[p]);
            cols += nodes[p].cols();
        }

        return result;
    }
};


}// namespace gismo
