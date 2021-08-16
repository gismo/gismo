/** @file gsCompositeIncrSmoothnessBasis.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsUtils/gsMesh/gsMesh.h>
#include <gsMSplines/gsMappedBasis.h>
#include <gsMSplines/gsMapFactory.h>

namespace gismo
{

template<class T> class gsCompositeMapperMatrix
{
    // (!) this class does not exist (!)
};

/** \brief
    Purely abstract class gsMappedBasis, which gives means of combining basis functions to new, global ones.

    This is done by having a set of basis functions (as a collection of basis patches) plus a
    topology, which specifies the relation of these patches. A mapping is then created, which
    combines basis functions with certain weights, to achieve a (predefined) increased smoothness degree
    over the interfaces.

    The function getMappedSingleBasis() then provides a class implementing gsBasis, which can be used also
    as a basis for a geometry.

    To split the gsMappedBasis into Patches again, the function exportToPatches() can be used.

    The notion of a localIndex is used for the accumulated indizes of the basis functions of the
    patches, whereas the globalIndex is used for the new, combined 'global' basis functions.
    PatchIndex is the index of a local basis function in its patch,
    so: localIndex = patchIndex + _firstLocalIndex(patch).

    \tparam d dimension of the parameter domains
    \tparam T coefficient type

    \ingroup basis
  */

template<short_t d, class T>
class gsCompositeIncrSmoothnessBasis : public gsMappedBasis<d,T>
{
private:
    /// Dimension of the parameter domains
    static const short_t Dim = d;

    typedef T weightType;
    typedef index_t indexType; //indizes of gsMatrix
    typedef std::vector<std::pair<int,int> >::iterator step_iter;
    typedef gsBasis<T> BasisType;
    typedef typename std::vector<BasisType *>::const_iterator ConstBasisIter;
    typedef typename std::vector<BasisType *>::iterator BasisIter;
    typedef typename std::vector<gsMatrix<T> *>::const_iterator ConstMatrixPtrIter;
    typedef typename gsSparseMatrix<weightType,0,indexType>::InnerIterator InIterMat;
    typedef typename std::vector<indexType> IndexContainer;
    typedef typename std::vector<indexType>::const_iterator ConstIndexIter;
    typedef typename std::vector<weightType> WeightContainer;
    typedef typename std::vector<weightType>::const_iterator ConstWeightIter;

    typedef memory::unique_ptr< gsDomainIterator<T> > domainIter;

protected:
    using gsMappedBasis<d,T>::m_mapper;
    using gsMappedBasis<d,T>::m_topol;
    using gsMappedBasis<d,T>::m_bases;
public:
    /// Shared pointer for gsCompositeIncrSmoothnessBasis
    typedef memory::shared_ptr< gsCompositeIncrSmoothnessBasis > Ptr;

    /// Unique pointer for gsCompositeIncrSmoothnessBasis
    typedef memory::unique_ptr< gsCompositeIncrSmoothnessBasis > uPtr;

    using gsMappedBasis<d,T>::size;
    using gsMappedBasis<d,T>::degree;
    using gsMappedBasis<d,T>::nPatches;

public:
    /// Default empty constructor
    gsCompositeIncrSmoothnessBasis() : gsMappedBasis<d,T>()
    {
//        m_mapper = NULL;
    }

    /// constructor, deletes all the bases and the mapper
    virtual ~gsCompositeIncrSmoothnessBasis()
    {

    } //destructor

public:
    //////////////////////////////////////////////////
    // functions for initializing and updating
    //////////////////////////////////////////////////

    /// updates the mapping of this basis (f.e. after a knot insertion)
    void updateTopol()
    {
        //_setDistanceOfAllVertices(0);
        _setMapping();
    }

protected:

    /// getter for m_bases[i]
    BasisType* getBasePointer(int i)
    { return m_bases[i]; }

    /// create a new mapping of the local basisfunctions
    void _setMapping();

    /// initializes the m_vertices field
    void _initVertices();

    /// initializes the m_distances field
    void _setDistanceOfAllVertices();

    /// gives back a mapfactory, which will create the map
    virtual gsMapFactory * _getMapFactory() = 0;

    /// Checks the gsMappedBasis for consistency
    bool _check() const;

    /// Checks the gsMappedBasis for consistency
    virtual bool _checkTopologyWithBases() const = 0;

public:
    //////////////////////////////////////////////////
    // general functions for interacting with this class
    //////////////////////////////////////////////////

    /// Returns the dimension \em d of the parameter space.
    short_t dim() const { return Dim; }

    /// Clone function. Used to make a copy of a derived basis
    GISMO_UPTR_FUNCTION_PURE(gsCompositeIncrSmoothnessBasis, clone)

    /// Prints the object as a string with extended details.
    std::string detail() const
    {
        // By default just uses print(..)
        std::ostringstream os;
        print(os);
        return os.str();
    }

    /// Prints the object to the stream
    std::ostream & print(std::ostream & os) const
    { os << m_mapper->asMatrix().toDense() <<"\n"; return os; }

    /// Returns the amount of basis functions on a given side of a given patch
    virtual unsigned basisFunctionsOnSide(const patchSide& ps) const = 0;

    // Look at gsBasis class for a description
    void connectivity(const gsMatrix<T> & nodes,gsMesh<T> & mesh) const;

    bool isConnected(indexType i,indexType j) const;

    virtual bool isLocallyConnected(indexType i,indexType j) const = 0;

public:
    //////////////////////////////////////////////////
    // getters for the private fields
    //////////////////////////////////////////////////

    /// getter for m_incrSmoothnessDegree
    int getIncrSmoothnessDegree() const
    { return m_incrSmoothnessDegree; }

    unsigned getMinDist() const
    { return m_minDist; }

public:
    //////////////////////////////////////////////////
    // functions for refinement
    //////////////////////////////////////////////////

    /** makes a uniform refinement of all the bases, with a defined number of new knots
     *  inserted in every knot span.
     *
     * \param[in] numKnots : the amount of knots which are inserted in every knot span
     * \param[in] updateBasis : if true, a new mapping will be constructed for the new basis functions,
     *                          only set to false if one wants to insert more boxes
     */
    void uniformRefine(int numKnots = 1, int mul=1,bool updateBasis = true);

    /** makes a uniform refinement of all the bases, with a defined number of new knots
     *  inserted in every knot span.
     *
     * \param[in,out] localCoefs : the coefficients to the local basis functions
     * \param[in] numKnots : the amount of knots which are inserted in every knot span
     * \param[in] updateBasis : if true, a new mapping will be constructed for the new basis functions,
     *                          only set to false if one wants to insert more boxes
     */
    void uniformRefine_withCoefs(gsMatrix<T>& localCoefs, int numKnots = 1, int mul=1,
                                 bool updateBasis = true);

    /** inserts the given boxes in the specified patch. The concrete implementation
     *  is left to the implementing classes of this abstract parent class.
     *
     * \param[in] patch : the patch number we look at
     * \param[in] boxes : the boxes to be inserted
     * \param[in] updateBasis : if true, a new mapping will be constructed for the new basis functions,
     *                          only set to false if one wants to insert more boxes
     */
    virtual void refine(const index_t patch, gsMatrix<T> const & boxes, bool updateBasis = true) = 0;

    /** inserts the given boxes in the specified patch. The concrete implementation
     *  is left to the implementing classes of this abstract parent class.
     *
     * \param[in] patch : the patch number we look at
     * \param[in] boxes : the boxes to be inserted
     * \param[in] updateBasis : if true, a new mapping will be constructed for the new basis functions,
     *                          only set to false if one wants to insert more boxes
     */
    virtual void refineElements(const index_t patch, std::vector<index_t> const & boxes, bool updateBasis = true) = 0;

    /** inserts the given boxes in the specified patch, while also updating the given coefficients. The concrete implementation
     *  is left to the implementing classes of this abstract parent class.
     *
     * \param[in,out] localCoefs : the coefficients to the local basis functions
     * \param[in] patch : the patch number we look at
     * \param[in] boxes : the boxes to be inserted
     * \param[in] updateBasis : if true, a new mapping will be constructed for the new basis functions,
     *                          only set to false if one wants to insert more boxes
     */
    virtual void refine_withCoefs(gsMatrix<T>& localCoefs, const index_t patch, gsMatrix<T> const & boxes,
                                  bool updateBasis = true) = 0;

    /** inserts the given boxes in the specified patch, while also updating the given coefficients. The concrete implementation
     *  is left to the implementing classes of this abstract parent class.
     *
     * \param[in,out] localCoefs : the coefficients to the local basis functions
     * \param[in] patch : the patch number we look at
     * \param[in] boxes : the boxes to be inserted
     * \param[in] updateBasis : if true, a new mapping will be constructed for the new basis functions,
     *                          only set to false if one wants to insert more boxes
     */
    virtual void refineElements_withCoefs(gsMatrix<T>& localCoefs, const index_t patch, std::vector<index_t> const & boxes,
                                  bool updateBasis = true) = 0;

    /** degree elevates all the bases preserving smoothness
     *
     * \param[in] amount : the amount of elevation the bases will go through
     * \param[in] updateBasis : if true, a new mapping will be constructed for the new basis functions,
     *                          only set to false if one wants to apply more changes to the bases before
     *                          calculating a new mapping
     */
    void degreeElevate( index_t amount = 1 , bool updateBasis = true);

    /** degree elevates all the bases preserving multiplicity
     *
     * \param[in] amount : the amount of elevation the bases will go through
     * \param[in] updateBasis : if true, a new mapping will be constructed for the new basis functions,
     *                          only set to false if one wants to apply more changes to the bases before
     *                          calculating a new mapping
     */
    void degreeIncrease( index_t amount = 1, index_t dir=-1, bool updateBasis = true);

    /// Adjusts the patches to match each other at the interfaces by knot insertion
    /// (will be called automatically by refinement functions)
    virtual void repairPatches(std::vector<gsMatrix<T> *> & coefs,
                               index_t startFromPatch = -1) = 0;

    void repairPatches(gsMatrix<T> & localCoef, index_t startFromPatch=-1);

    /// Adjusts the patches to match each other at the interfaces by knot insertion
    /// (will be called automatically by refinement functions)
    void repairPatches(index_t startFromPatch = -1);

public:
    //////////////////////////////////////////////////
    // functions for smoothing the basis
    //////////////////////////////////////////////////

    /** smooths a given edge connected to the corner back to the the minimal amount of C^0 continuity,
     *  specified by m_minDist.
     *
     * @param[in] pc : patchCorner, specifying the corner which should be smoothed
     * @param[in] ps : patchSide, specifying the edge, connected to the corner, which should be smoothed along
     * @param[in] updateBasis : if true, a new mapping will be constructed for the new basis functions,
     *                          only set to false if one wants to insert more boxes
     */
    void smoothCornerEdge(const patchCorner&pc,const patchSide& ps,bool updateBasis = true);

    /** smooths all edges back to the the minimal amount of C^0 continuity,
     *  specified by m_minDist.
     */
    void smoothEverything();

public:
    //////////////////////////////////////////////////
    // functions for working with weights on interfaces
    //////////////////////////////////////////////////

    /** gets the weight at the given patchSide.
     *
     * \param[in] ps : patchSide
     */
    T getWeight(const patchSide & ps) const;

    /** sets the weight at the given patchSide to value weight
     *
     * \param[in] ps : patchSide
     * \param[in] weight : T
     */
    bool setWeight(const patchSide & ps, const T weight);

public:
    //////////////////////////////////////////////////
    // functions for working with special vertices and it's C^0 distance
    //////////////////////////////////////////////////

    /// sets a patchCorner to be a special vertex, which means edges connected
    /// to this corner will have a c^0 continuity
    void setC0(patchCorner pc);

    /// gives back true, if the given patchCorner is a special vertex
    bool isSpecialVertex(const patchCorner & pc) const;

    /// gives back the parametric c^0 distance of the edge \a ps starting from corner \a pc
    T getParametricDistanceOfVertex(const patchCorner& pc,const patchSide& ps) const;

protected:
    /// finds the parametric c^0 distance of the edge \a ps starting from corner \a pc
    /// if the number of c^0 basis functions is given in \a nrBasisFuncs
    virtual T findParameter(patchSide const & ps,patchCorner const & pc,unsigned nrBasisFuncs) const = 0;

private:
    //////////////////////////////////////////////////
    // struct for distances from special vertices
    //////////////////////////////////////////////////

    /** \brief
        Private stract that has the purpose of storing distance information of c^0 parts around special vertices.

        One struct ressembles an edge between to vertices, which is able to provide the distance at both its ends.
      */

    struct distances
    {
        boundaryInterface interface;
        patchCorner corner1;
        T parametricDistance1;
        patchCorner corner2;
        T parametricDistance2;

        /// Construct distances by a given interface, two patchCorners and the associated
        /// gsMappedBasis
        distances(const boundaryInterface&  iface, const patchCorner& pc1,
                  const patchCorner& pc2,const gsCompositeIncrSmoothnessBasis<d,T>& basis);

        /// checks if this distances struct ressembles the interface given
        bool isDistancesOfInterface(const boundaryInterface& bi) const
        {
            return (interface==bi) || (interface.getInverse()==bi);
        }

        /// sets the parametric distance of the edge starting from \a pc
        /// to get the c^0 distance from \a absoluteVal basisfunctions.
        void setParamDist(unsigned absoluteVal,const patchCorner& pc,const gsCompositeIncrSmoothnessBasis<d,T>& basis);

        /// gets the parametric distance from the corner \a pc
        T getParamDist(const patchCorner& pc,const gsCompositeIncrSmoothnessBasis<d,T>& basis) const;

        /// determines the right values for the two distances, only used in the constructer
        void _determineValues(patchSide side,patchSide ls,patchSide rs,int dist,unsigned degree,unsigned max,
                              unsigned& left,unsigned& right,const gsCompositeIncrSmoothnessBasis<d,T>& basis) const;
    };

    // Data members
protected:
    /// smoothness degree that is tried to achive over patch interfaces
    short_t m_incrSmoothnessDegree;
    /// minimal C^0 distance from special (extraordinary) vertices, specified in basisfunctions
    unsigned m_minDist;
    /// vector storing the weights for interfaces
    std::vector<std::pair<patchSide,T> > m_patchSideWeights;
    /// vector storing all the inner vertices, with a flag if it is an Extraordinary vertex or an ordinary vertex.
    std::vector<std::pair<patchCorner,bool> > m_vertices; //true = EV, false = OV
    /// vector of distances objects, that store C^0 distances from special vertices on the edges
    std::vector<distances> m_distances;
}; // class gsMappedBasis

}


//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCompositeIncrSmoothnessBasis.hpp)
#endif
