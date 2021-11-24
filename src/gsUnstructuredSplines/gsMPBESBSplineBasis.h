/** @file gsMPBESBSplineBasis.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsMSplines/gsMPBESBasis.h>
#include <gsMSplines/gsMPBESMapB2D.h>
#include <gsNurbs/gsKnotVector.h>

namespace gismo
{

/** @brief
      A univariate Lagrange basis.

      \tparam T coefficient type

      \ingroup basis
  */

template<short_t d, class T>
class gsMPBESBSplineBasis : public gsMPBESBasis<d,T>
{
public:
    /// Shared pointer for gsMPBESBSplineBasis
    typedef memory::shared_ptr< gsMPBESBSplineBasis > Ptr;

    /// Unique pointer for gsMPBESBSplineBasis
    typedef memory::unique_ptr< gsMPBESBSplineBasis > uPtr;

    /// Dimension of the parameter domain
    static const short_t Dim = d;

    typedef index_t                                             indexType;
    typedef gsMPBESBasis<d,T>                                   Base;
    typedef gsTensorBSplineBasis<d,T>                           BasisType;
    typedef typename std::vector<BasisType *>                   BasisContainer;
    typedef typename std::vector<gsBasis<T>* >::const_iterator  ConstBasisIter;
    typedef typename std::vector<gsBasis<T>* >::iterator        BasisIter;
    typedef typename std::vector<gsMatrix<T> *>::const_iterator ConstMatrixPtrIter;

protected:
    using Base::m_bases;
    using Base::m_topol;
    using Base::m_incrSmoothnessDegree;
    using Base::m_mapper;
    using Base::m_vertices;
    using Base::m_distances;
    using Base::m_minDist;
public:
    using Base::updateTopol;
    using Base::degree;
    using Base::global_coef_to_local_coef;
    using Base::local_coef_to_global_coef;
    using Base::exportToPatches;
    using Base::nPatches;
    using Base::isSpecialVertex;
    using Base::repairPatches;
    using Base::maxDegree;
protected:
    using Base::_getPatch;
    using Base::_getPatchIndex;
    using Base::_initVertices;
    using Base::_setDistanceOfAllVertices;

public:
    /// Default empty constructor
    gsMPBESBSplineBasis()
    { }

    gsMPBESBSplineBasis (const BasisContainer & bases, gsBoxTopology const & topol,
                             int increaseSmoothnessLevel=-1, int minEVDistance=-1);

    gsMPBESBSplineBasis (const BasisContainer & bases, gsBoxTopology const & topol,std::vector<gsMatrix<T> * > coefs,
                             int increaseSmoothnessLevel=-1, int minEVDistance=-1);

    gsMPBESBSplineBasis( gsMultiPatch<T> const & mp, int increaseSmoothnessLevel = -1,
                             int minEVDistance=-1);

    gsMPBESBSplineBasis( const gsMPBESBSplineBasis& other );

    gsMPBESBSplineBasis<d,T> & operator=( const gsMPBESBSplineBasis& other );

    ~gsMPBESBSplineBasis()
    { } //destructor

private:
    //////////////////////////////////////////////////
    // functions for initializing and updating
    //////////////////////////////////////////////////

    bool _checkTopologyWithBases() const
    {
//        T eps=1e-6;
        bool consistent = true;
//        int nrOfBases = m_bases.size();
//        consistent = consistent && nrOfBases>0;
//        GISMO_ASSERT(nrOfBases>0,"empty list of bases");
//        consistent = consistent && nrOfBases==m_topol.size();
//        GISMO_ASSERT(nrOfBases==m_topol.size(),"bases and topol not of same size.");
//        unsigned i = 0;
//        for(ConstBasisIter it=m_bases.begin();it!=m_bases.end();++it)
//        {
//            int deg_u = basis(i).degree(0);
//            int deg_v = basis(i).degree(1);
//            unsigned u_dim = basis(i).size(0), v_dim = basis(i).size(1);
//            consistent = consistent && u_dim >= static_cast<unsigned>(deg_u*2);
//            consistent = consistent && v_dim >= static_cast<unsigned>(deg_v*2);
//            if(u_dim<static_cast<unsigned>(deg_u*2))
//                GISMO_ERROR("patch is too small: u_dim too small");
//            if(v_dim<static_cast<unsigned>(deg_v*2))
//                GISMO_ERROR("patch is too small: v_dim too small");
//            i++;
//        }
//        std::vector<boundaryInterface> interfaces = m_topol.interfaces();
//        for(unsigned i=0;i<interfaces.size();i++)
//        {
//            patchSide first() = interfaces[i].first();
//            patchSide second() = interfaces[i].second();
//            int nrOfBasisFuncs1=basis(first().patch).size(!direction(first().side));
//            int nrOfBasisFuncs2=basis(second().patch).size(!direction(second().side));
//            consistent = consistent && nrOfBasisFuncs1==nrOfBasisFuncs2;
//            if(!(nrOfBasisFuncs1==nrOfBasisFuncs2))
//                GISMO_ERROR("number of basisfunctions on interface are different");
//            int deg_1 = m_bases[first().patch]->degree(!direction(first().side));
//            int deg_2 = m_bases[first().patch]->degree(!direction(second().side));
//            consistent = consistent && deg_1==deg_2;
//            if(!(deg_1==deg_2))
//                GISMO_ERROR("degrees of patches on interface are different");
//            gsKnotVector<T> kv1 = basis(first().patch).knots(!direction(first().side));
//            gsKnotVector<T> kv2 = basis(second().patch).knots(!direction(second().side));
//            if(kv1.size()!=kv2.size())
//                GISMO_ERROR("different knot vector lengths");
//            if(!interfaces[i].orient[0])
//                kv2.reverse();
//            for(int j = 0;j<kv1.size();j++)
//            {
//                if(std::abs(kv1[j]-kv2[j])>eps)
//                {
//                    GISMO_ERROR("different knot elements.");
//                }
//            }
//        }
        return consistent;
    }

protected:
    void _setMapping();

public:
    //////////////////////////////////////////////////
    // general functions for interacting with this class
    //////////////////////////////////////////////////

    /// Clone function. Used to make a copy of a derived basis
    GISMO_CLONE_FUNCTION(gsMPBESBSplineBasis)

    gsTensorBSplineBasis<d,T> & basis(size_t i)
    { return static_cast<gsTensorBSplineBasis<d,T>&>(*m_bases[i]); }

    const gsTensorBSplineBasis<d,T> & basis(size_t i) const
    { return static_cast<const gsTensorBSplineBasis<d,T>&>(*m_bases[i]); }

    unsigned basisFunctionsOnSide(const patchSide& ps) const;

    virtual bool isLocallyConnected(indexType i,indexType j) const;

public:
    //////////////////////////////////////////////////
    // functions for evaluating and derivatives
    //////////////////////////////////////////////////

    /// Returns the number of active (nonzero) basis functions at points \a u in \a result.
    void numActive_into(const index_t patch,const gsMatrix<T> & u, gsVector<index_t>& result) const;

public:
    //////////////////////////////////////////////////
    // functions for refinement
    //////////////////////////////////////////////////

    void refine(const index_t patch,const std::vector<T>& knots_u, const std::vector<T>& knots_v,bool updateBasis = true);

    void refine(const index_t patch, gsMatrix<T> const & boxes, bool updateBasis = true);

    void refineElements(const index_t patch, std::vector<index_t> const & boxes, bool updateBasis = true);

    void refine_withCoefs(gsMatrix<T>& localCoef, const index_t patch,const std::vector<T>& knots_u, const std::vector<T>& knots_v,
                          bool updateBasis = true);

    void refine_withCoefs(gsMatrix<T>& coefs, const index_t patch, gsMatrix<T> const & boxes,
                          bool updateBasis = true);

    void refineElements_withCoefs(gsMatrix<T>& coefs, const index_t patch, std::vector<index_t> const & boxes,
                          bool updateBasis = true);

private:
    //////////////////////////////////////////////////
    // private helper functions for the refinement
    //////////////////////////////////////////////////

    void _boxesVectorToMatrix(const std::vector<index_t> & boxes, gsMatrix<T> & mat_boxes);

    void _boxToRefineKnots(const index_t patch,gsMatrix<T> const & boxes,std::vector<std::vector<T> > & refineKnots);

    // takes a knotvector, a standard vector of new knots which will be inserted, an indicator par which tells
    // if the start or the end of the knotvector has to be checked for new knots, and a distance indicating the distance
    // from start/end where new knots are looked for. Returns the number of new knots inserted in the knotspan, going
    // distance knots from the beginning or end of the given knotvector.
    unsigned getNrOfSpecialKnots(const gsKnotVector<T> kv,const std::vector<T>& new_knots,bool par,int distance);

    void repairPatches(std::vector<gsMatrix<T> *> & coefs,
                        index_t startFromPatch = -1);

    bool _knotsMatchNeighbours(index_t patch,std::vector<std::vector<T> >& knotsToInsert,
                               std::vector<index_t>& checkPatches,
                               T eps = 1e-10);

protected:
    //////////////////////////////////////////////////
    // functions going back and forth between absolute val and parametric value for C^0 parts
    //////////////////////////////////////////////////

    T findParameter(patchSide const & ps,patchCorner const & pc,unsigned nrBasisFuncs) const;

}; // class gsMPBESBSplineBasis

}

