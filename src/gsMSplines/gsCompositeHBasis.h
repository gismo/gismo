/** @file gsCompositeHBasis.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsMSplines/gsCompositeIncrSmoothnessBasis.h>
#include <gsMSplines/gsCompositeMapFactoryHB2D.h>
#include <gsHSplines/gsHTensorBasis.h>

namespace gismo
{

/** @brief
      A univariate Lagrange basis.

      \tparam T coefficient type

      \ingroup basis
  */

template<short_t d, class T>
class gsCompositeHBasis : public gsCompositeIncrSmoothnessBasis<d,T>
{
public:
    /// Shared pointer for gsCompositeHBasis
    typedef memory::shared_ptr< gsCompositeHBasis > Ptr;

    /// Unique pointer for gsCompositeHBasis
    typedef memory::unique_ptr< gsCompositeHBasis > uPtr;

    /// Dimension of the parameter domain
    static const int Dim = d;

    typedef index_t indexType;
    typedef gsCompositeIncrSmoothnessBasis<d,T> Base;
    typedef gsHTensorBasis<d,T> BasisType;
    typedef typename std::vector<gsBasis<T> *>::const_iterator ConstBasisIter;

protected:
    using Base::m_bases;
    using Base::m_topol;
    using Base::m_incrSmoothnessDegree;
    using Base::m_mapper;
    using Base::m_vertices;
    using Base::m_distances;
    using Base::m_minDist;
public:
    using Base::nPatches;
    using Base::updateTopol;
    using Base::degree;
    using Base::repairPatches;
    using Base::maxDegree;
protected:
    using Base::_getPatch;
    using Base::_getPatchIndex;
    using Base::_setDistanceOfAllVertices;
    using Base::_setMapping;
    using Base::_initVertices;

//TODO: gsCompactKnotVector to gsKnotVector,
//      shift knots by constant
//      read in of thb multipatch
//      compute topology of thb multipatch

public:
    /// Default empty constructor
    gsCompositeHBasis()
    { }

    gsCompositeHBasis (std::vector<BasisType *> const & bases, gsBoxTopology const & topol,
                       int increaseSmoothnessLevel=-1, int minEVDistance=-1);

    gsCompositeHBasis (std::vector<BasisType *> const & bases, gsBoxTopology const & topol,
                       std::vector<gsMatrix<T> * > & coefs,int increaseSmoothnessLevel=-1, int minEVDistance=-1);

    gsCompositeHBasis (BasisType const & base, gsBoxTopology const & topol);

    gsCompositeHBasis( gsMultiPatch<T> const & mp, int increaseSmoothnessLevel = -1,
                       int minEVDistance=-1);

    gsCompositeHBasis( gsMultiBasis<T> const & mb,gsBoxTopology const & topol, int increaseSmoothnessLevel = -1,
                       int minEVDistance=-1);

    gsCompositeHBasis( const gsCompositeHBasis& other );

    gsCompositeHBasis<d,T>& operator=( const gsCompositeHBasis& other );

    ~gsCompositeHBasis()
    { 

    } //destructor

private:
    //////////////////////////////////////////////////
    // functions for initializing and updating
    //////////////////////////////////////////////////

    gsCompositeMapFactoryHB2D<d,T,gsCompositeMapperMatrix<T> > * _getMapFactory()
    {
        return new gsCompositeMapFactoryHB2D<d,T,gsCompositeMapperMatrix<T> >(m_incrSmoothnessDegree,& m_topol,this);
    }

    bool _checkTopologyWithBases() const
    { return true; }

public:
    //////////////////////////////////////////////////
    // general functions for interacting with this class
    //////////////////////////////////////////////////

    /// Clone function. Used to make a copy of a derived basis
    GISMO_CLONE_FUNCTION(gsCompositeHBasis)

    gsHTensorBasis<d,T> & basis(size_t i)
    { return static_cast<gsHTensorBasis<d,T>&>(*m_bases[i]); }

    const gsHTensorBasis<d,T> & basis(size_t i) const 
    { return static_cast<const gsHTensorBasis<d,T>&>(*m_bases[i]); }

    unsigned basisFunctionsOnSide(const patchSide& ps) const;

    virtual bool isLocallyConnected(indexType i,indexType j) const;

public:
    //////////////////////////////////////////////////
    // functions for evaluating and derivatives
    //////////////////////////////////////////////////

    /// Returns the number of active (nonzero) basis functions at points \a u in \a result.
    void numActive_into(const unsigned patch,const gsMatrix<T> & u, gsVector<unsigned>& result) const
    {
        GISMO_UNUSED(patch); GISMO_UNUSED(u); GISMO_UNUSED(result);
        GISMO_NO_IMPLEMENTATION
    }

public:
    //////////////////////////////////////////////////
    // functions for refinement
    //////////////////////////////////////////////////

    void refine(const index_t patch, gsMatrix<T> const & boxes, bool updateBasis = true);

    void refineElements(const index_t patch, std::vector<index_t> const & boxes, bool updateBasis = true);

    void refine_withCoefs(gsMatrix<T>& localCoef, const index_t patch, gsMatrix<T> const & boxes, bool updateBasis = true);

    void refineElements_withCoefs(gsMatrix<T>& localCoef, const index_t patch, std::vector<index_t> const & boxes, bool updateBasis = true);

    /// @brief Refine the are defined by \em boxes
    /// on patch \em k with extension \em refExt.
    ///
    /// See gsHTensorBasis::refine() for further documentation.
    virtual void refineWithExtension(const index_t patch,gsMatrix<T> const & boxes, int refExt = 0,bool updateBasis = true);

//private:
    //////////////////////////////////////////////////
    // private helper functions for the refinement
    //////////////////////////////////////////////////

    void repairPatches(std::vector<gsMatrix<T> *> & coefs, index_t startFromPatch = -1);

    bool _innerBoxesAreSuitable(const index_t patch,
                                std::vector<index_t>& boxes);

    bool _boxesMatchNeighbours(const index_t patch,
                               std::vector<index_t>& boxes, std::vector<index_t>& checkPatches);

    void _addBoundaryBox(const index_t patch,const boxSide s,const int start, const int end,const unsigned level, std::vector<index_t> & boxes, std::vector<bool> & sideToCheck);

    void _addFunBox(const index_t patch,const unsigned uMin,const unsigned vMin,const unsigned uMax,const unsigned vMax,const unsigned level, std::vector<index_t> & boxes);

    void _addBox(const index_t patch,const unsigned uMin,const unsigned vMin,const unsigned uMax,const unsigned vMax,const unsigned level, std::vector<index_t> & boxes);

protected:
    //////////////////////////////////////////////////
    // functions going back and forth between absolute val and parametric value for C^0 parts
    //////////////////////////////////////////////////

    void _endpointsOfActiveBoundaryFunctions(patchSide const & ps,bool orient,std::vector<T>& endpoints) const;

    T findParameter(patchSide const & ps,patchCorner const & pc,unsigned nrBasisFuncs) const;

}; // class gsCompositeHBasis

}
