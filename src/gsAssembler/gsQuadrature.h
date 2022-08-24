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
#include <gsAssembler/gsPatchRule.h>
#include <gsAssembler/gsOverIntegrateRule.h>

namespace gismo
{

/// Helper class for obtaining a quadrature rule
struct gsQuadrature
{
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
        const T       quA = options.getReal("quA");
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
        const T       quA  = options.getReal("quA");
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
                const T       quAb  = options.askReal("quAb",quA+1);
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
                               const T quA, const index_t quB, short_t fixDir = -1)
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
            nnodes[i] = cast<T,index_t>(quA * basis.degree(i) + quB + 0.5);
        for(++i; i<d; ++i )
            nnodes[i] = cast<T,index_t>(quA * basis.degree(i) + quB + 0.5);
        return nnodes;
    }
};

}// namespace gismo
