/** @file gsTHBSplineBasis.h

    @brief Provides declaration of THBSplineBasis class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.


    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#pragma once

#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsTHBSpline.h>


namespace gismo
{  

/**
 * \brief
 * Truncated hierarchical B-spline basis.
 *
 * \param d the dimension of the parameter domain
 * \param T the coefficient type
 *
 * \ingroup basis
 * \ingroup HSplines
 */
template<short_t d, class T>
class gsTHBSplineBasis : public gsHTensorBasis<d,T>
{
public:
    /// @brief Associated geometry type.
    typedef gsTHBSpline<d,T> GeometryType;
    
    typedef typename gsHTensorBasis<d,T>::CMatrix CMatrix;

    typedef typename gsHTensorBasis<d,T>::cmatIterator cmatIterator;
    
    typedef typename gsHTensorBasis<d,T>::tensorBasis tensorBasis;

    /// @brief Shared pointer for gsTHBSplineBasis.
    typedef memory::shared_ptr< gsTHBSplineBasis > Ptr;

    /// @brief Unique pointer for gsTHBSplineBasis.
    typedef memory::unique_ptr< gsTHBSplineBasis > uPtr;

    /// @brief Associated Boundary basis type.
    typedef typename
    util::conditional<d==1, gsConstantBasis<T>, gsTHBSplineBasis<static_cast<short_t>(d-1),T> >::type BoundaryBasisType;

    using gsHTensorBasis<d, T>::flatTensorIndexOf;

    // polygon lines in parameter domain
    // the stucture is [levels [ line [ segments [ x y z w ] ] ] ],
    // where x y z w describes segment from (x, y) to (z, w)
    typedef typename std::vector< std::vector< std::vector< std::vector<T> > > > Polylines;

    // axis aligned bounding boxes in parameter domain
    // the structure is [levels [ boxes [ low_x low_y upp_x upp_y] ] ]
    // where the box is define by lower left corner (low_x, low_y) and upper right 
    // corner (upp_x, upp_y)
    typedef typename std::vector< std::vector< std::vector<index_t> > > AxisAlignedBoundingBox;

    //typedef typename std::vector< gsAabb<T> > AxisAlignedBoundingBox;

    // trimming curves in parameter domain
    // the stucture is [level [connected componenet [ line [ segments [ x y z w ] ] ] ] ],
    // where x y z w describes segment from (x, y) to (z, w), and first line 
    // indicates outer loop and next line describes holes
    typedef typename std::vector< std::vector< std::vector< std::vector< std::vector<T> > > > > TrimmingCurves;

public:

    gsTHBSplineBasis()
    { 
        representBasis(); 
    }

    /// @brief Constructor out of a Tensor BSpline Basis
    gsTHBSplineBasis(gsTensorBSplineBasis<d,T> const&  tbasis) 
    : gsHTensorBasis<d,T>(tbasis) 
    { representBasis(); }

    gsTHBSplineBasis(gsTensorBSplineBasis<d,T> const&  tbasis, 
                     const std::vector<index_t> & boxes) 
    : gsHTensorBasis<d,T>(tbasis, boxes)
    { representBasis(); }

    gsTHBSplineBasis(gsTensorBSplineBasis<d,T> const&  tbasis, 
                     gsMatrix<T> const & boxes)
    : gsHTensorBasis<d,T>(tbasis, boxes) 
    {  representBasis(); }

    gsTHBSplineBasis( gsTensorBSplineBasis<d,T> const&  tbasis, 
                      gsMatrix<T> const & boxes, 
                      const std::vector<index_t> & levels)
    : gsHTensorBasis<d,T>(tbasis, boxes, levels)
    {  representBasis(); }

    /// @brief Constructor out of a tensor BSpline Basis
    gsTHBSplineBasis(gsBasis<T> const&  tbasis)
        : gsHTensorBasis<d,T>(tbasis)
    {  representBasis(); }

#ifdef __DOXYGEN__
    /// @brief Gives back the boundary basis at boxSide s
    typename BoundaryBasisType::uPtr boundaryBasis(boxSide const & s);
#endif
    GISMO_UPTR_FUNCTION_DEF(BoundaryBasisType, boundaryBasis, boxSide const &)
    {
        return basisSlice(n1.direction(),n1.parameter());
    }

public:
    /// @brief Gives back the basis at a slice in \a dir_fixed at \a par
    BoundaryBasisType * basisSlice(index_t dir_fixed,T par ) const;

    // Look at gsBasis class for documentation
    void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result)const;

    // Look at gsBasis class for documentatation
    void deriv2Single_into(index_t i,
                           const gsMatrix<T>& u,
                           gsMatrix<T>& result) const;

    // Look at gsBasis class for documentation
    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
 
    // Look at gsBasis class for documentation
    void derivSingle_into(index_t i,
                             const gsMatrix<T> & u,
                          gsMatrix<T>& result) const;

    // look at eval_into
    void fastEval_into(const gsMatrix<T>& u,
                       gsMatrix<T>& result) const
    {
        gsMatrix<index_t> indices;
        this->active_into(u, indices);

        result.setZero(indices.rows(), u.cols());

        const unsigned maxLvl = this->m_tree.getMaxInsLevel() + 1;
        std::vector< gsMatrix<T> > tmpResults(maxLvl, gsMatrix<T>());
        std::vector< gsMatrix<index_t> > tmpActive(maxLvl, gsMatrix<index_t>());
        gsVector<int> processed(maxLvl);

        for (int pt = 0; pt != u.cols(); pt++)
        {
            processed.setZero();

            for (int ind = 0; ind != indices.rows(); ind++)
            {
                unsigned index = indices(ind, pt);
                if (ind != 0 && index == 0)
                    break;

                unsigned lvl = getPresLevelOfBasisFun(index);

                if (processed(lvl) == 0)
                {
                    this->m_bases[lvl]->eval_into(u.col(pt), tmpResults[lvl]);
                    this->m_bases[lvl]->active_into(u.col(pt), tmpActive[lvl]);
                    processed(lvl) = 1;
                }

                if (m_is_truncated[index] == -1)
                {
                    index_t flatTenIndx = this->flatTensorIndexOf(index, lvl);
                    int localIndex = -1;
                    for (int row = 0; row != tmpActive[lvl].rows(); ++row)
                    {
                        if (tmpActive[lvl](row, 0) == flatTenIndx)
                        {
                            localIndex = row;
                            break;
                        }
                    }
                    result(ind, pt) = tmpResults[lvl](localIndex, 0);
                }
                else // basis function is truncated
                {
                    const gsMatrix<T>& basis = tmpResults[lvl];
                    const gsMatrix<index_t>& active = tmpActive[lvl];


                    const gsSparseVector<T>& coefs = getCoefs(index);
                    T tmp = coefs(active(0, 0)) * basis(0, 0);
                    for (int i = 1; i < active.rows(); i++)
                    {
                        tmp += coefs(active(i, 0)) * basis(i, 0);
                    }
                    result(ind, pt) = tmp;
                }
            }
        }
    }

    // look at deriv_into
    void fastDeriv_into(const gsMatrix<T>& u,
                        gsMatrix<T>& result) const
    {
        gsMatrix<index_t> indices;
        this->active_into(u, indices);

        result.setZero(indices.rows() * d, u.cols());

        const unsigned maxLvl = this->m_tree.getMaxInsLevel() + 1;
        std::vector< gsMatrix<T> > tmpDeriv( maxLvl, gsMatrix<T>());
        std::vector< gsMatrix<index_t> > tmpActive(maxLvl, gsMatrix<index_t>());
        gsVector<int> processed(maxLvl);

        for (index_t pt = 0; pt != u.cols(); pt++)
        {
            processed.setZero();

            for (int ind = 0; ind != indices.rows(); ind++)
            {
                index_t index = indices(ind, pt);
                if (ind != 0 && index == 0)
                    break;

                unsigned lvl = getPresLevelOfBasisFun(index);

                if (processed(lvl) == 0)
                {
                    this->m_bases[lvl]->deriv_into(u.col(pt), tmpDeriv[lvl]);
                    this->m_bases[lvl]->active_into(u.col(pt), tmpActive[lvl]);
                    processed(lvl) = 1;
                }

                if (m_is_truncated[index] == -1)
                {
                    index_t flatTenIndx = this->flatTensorIndexOf(index, lvl);
                    int localIndex = -1;
                    for (int row = 0; row != tmpActive[lvl].rows(); ++row)
                    {
                        if (tmpActive[lvl](row, 0) == flatTenIndx)
                        {
                            localIndex = row;
                            break;
                        }
                    }

                    result.block(ind * d, pt, d, 1) =
                            tmpDeriv[lvl].block(localIndex * d, 0, d, 1);
                }
                else // basis function is truncated
                {
                    const gsMatrix<T>& basis = tmpDeriv[lvl];
                    const gsMatrix<index_t>& active = tmpActive[lvl];
                    const gsSparseVector<T>& coefs = getCoefs(index);

                    for (unsigned dim = 0; dim != d; dim++) // for all deric
                    {
                        for (index_t i = 0; i != active.rows(); ++i)
                        {
                            result(ind * d + dim, pt) +=
                                    coefs(active(i, 0)) * basis(i * d + dim, 0);
                        }
                    }
                }
            }
        }
    }

    // look at deriv2_into
    void fastDeriv2_into(const gsMatrix<T>& u,
                         gsMatrix<T>& result) const
    {
        gsMatrix<index_t> indices;
        this->active_into(u, indices);
        const unsigned numDers = (d * (d + 1)) / 2;

        result.setZero(indices.rows() * numDers, u.cols());

        const unsigned maxLvl = this->m_tree.getMaxInsLevel() + 1;
        std::vector< gsMatrix<T> > tmpDeriv2(maxLvl, gsMatrix<T>());
        std::vector< gsMatrix<index_t> > tmpActive(maxLvl, gsMatrix<index_t>());
        gsVector<int> processed(maxLvl);

        for (int pt = 0; pt != u.cols(); pt++)
        {
            processed.setZero();

            for (int ind = 0; ind != indices.rows(); ind++)
            {

                unsigned index = indices(ind, pt);
                if (ind != 0 && index == 0)
                    break;

                unsigned lvl = getPresLevelOfBasisFun(index);

                if (processed(lvl) == 0)
                {
                    this->m_bases[lvl]->deriv2_into(u.col(pt), tmpDeriv2[lvl]);
                    this->m_bases[lvl]->active_into(u.col(pt), tmpActive[lvl]);
                    processed(lvl) = 1;
                }

                if (m_is_truncated[index] == -1)
                {
                    index_t flatTenIndx = this->flatTensorIndexOf(index, lvl);
                    int localIndx = -1;
                    for (int row = 0; row != tmpActive[lvl].rows(); ++row)
                    {
                        if (tmpActive[lvl](row, 0) == flatTenIndx)
                        {
                            localIndx = row;
                            break;
                        }
                    }

                    result.block(ind * numDers, pt, numDers, 1) =
                            tmpDeriv2[lvl].block(localIndx * numDers, 0, numDers, 1);
                }
                else // basis function is truncated
                {
                    const gsMatrix<T>& basis = tmpDeriv2[lvl];
                    const gsMatrix<index_t>& active = tmpActive[lvl];
                    const gsSparseVector<T>& coefs = getCoefs(index);

                    for (unsigned der = 0; der != numDers; der++)
                    {
                        for (index_t i = 0; i != active.rows(); ++i)
                        {
                            result(ind * numDers + der, pt) +=
                                    coefs(active(i, 0)) *
                                    basis(i * numDers + der, 0);
                        }
                    }
                }
            }
        }
    }

    // Look at gsBasis class for documentation
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    // Because of overriding one of the "eval_into" functions, all
    // functions in the base class with this name are hidden from the
    // derived class: Compiler does not search the base class as soon
    // as the function name is found in the derived class.
    // Therefore we need to add a "using" declaration to "reveal" the
    // function again. THis brings into scope all definitions of
    // eval_into from the base class.
    using gsBasis<T>::eval_into;

    /// @brief Returns the number of truncated basis functions
    unsigned numTruncated() const
    { return m_presentation.size(); }

    bool isTruncated(unsigned i) const
    {
        return (this->m_is_truncated[i] != -1);
    }

    /// @brief Returns an iterator to the representation of the first truncated basis function
    typename std::map<index_t, gsSparseVector<T> >::const_iterator truncatedBegin() const
    { return m_presentation.begin(); }

    /// @brief Returns an iterator past the last truncated basis function
    typename std::map<index_t, gsSparseVector<T> >::const_iterator truncatedEnd() const
    { return m_presentation.end(); }

    /// @brief Returns sparse representation of the i-th basis function.
    const gsSparseVector<T>& getCoefs(unsigned i) const
    {
        if (this->m_is_truncated[i] == -1)
        {
            GISMO_ERROR("This basis function has no sparse representation. "
                        "It is not truncated.");
        }
        else
        {
            return this->m_presentation.find(i)->second;
        }
    }


    void evalSingle_into(index_t i,
                         const gsMatrix<T>& u,
                         gsMatrix<T>& result) const;

private:

    unsigned getPresLevelOfBasisFun(const unsigned index) const
    {
        if (m_is_truncated[index] == -1)
        {
            return this->levelOf(index);
        }
        else
        {
            return m_is_truncated[index];
        }
    }

    /// @brief Computes and saves representation of all basis functions.
    void representBasis(); // rename: precompute coeffs


    /// @brief Computes representation of j-th basis function on pres_level and
    /// saves it.
    ///
    /// @param j     index of basis function
    /// @param pres_level levet at which we want to present j-th basis function
    /// @param finest_low "low index" of support of j-th basis function (finest grid)
    /// @param finest_high "high index" of support of j-th basis function (finest grid)
    void _representBasisFunction(const unsigned j,
                                const unsigned pres_level,
                                const gsVector<index_t, d>& finest_low,
                                const gsVector<index_t, d>& finest_high);



    /// @brief Saves a presentation of the j-th basis function. Presentation is given
    /// by the coefficients coefs. The coefficients corresponds to the
    /// apropriate BSplines at presentation level.
    ///
    /// @param coefs coefficients
    /// @param act_size_of_coefs size of the coefficients
    /// @param j global tensor index of a basis function
    /// @param pres_level presentation level
    /// @param finest_low "low index" of support of j-th basis function
    ///        (at the finest grid)
    void _saveNewBasisFunPresentation(const gsMatrix<T>& coefs,
                                      const gsVector<index_t, d>& act_size_of_coefs,
                                      const unsigned j,
                                      const unsigned pres_level,
                                      const gsVector<index_t, d>& finest_low);



    /// @brief Computes tensor index of a basis function on a finer level (new_level)
    /// which is presented with tensor index (index) at a coarse level (level).
    ///
    /// @param index tensor index of the basis function on a coarse level
    /// @param level coarse level
    /// @param fin_low "low index" of the support of the basis function
    /// @param new_level finer level
    ///
    /// @return global tensor index of a basis function on the finer level
    unsigned _basisFunIndexOnLevel(const gsVector<index_t, d>& index,
                                   const unsigned level,
                                   const gsVector<index_t, d>& fin_low,
                                   const unsigned new_level);



    /// @brief Performs truncation.
    ///
    /// @param coefs coefficients refined from level - 1 to level
    /// @param act_size_of_coefs actuall size of the coefficients
    /// @param size_of_coefs size of the non-zero coeffs
    /// @param level at which we truncate
    /// @param bspl_vec_ti B-Spline tensor index vector of the original function
    ///                    (we want to trancate this basis function)
    /// @param bspl_vec_ti_level level of the basis function given with vector
    ///                          bspl_vec_ti
    /// @param finest_low "low index" of the support of the basis function
    void _truncate(gsMatrix<T>& coefs,
                   const gsVector<index_t, d>& act_size_of_coefs,
                   const gsVector<index_t, d>& size_of_coefs,
                   const unsigned level,
                   const gsVector<index_t, d>& bspl_vec_ti,
                   const unsigned bspl_vec_ti_level,
                   const gsVector<index_t, d>& finest_low);


    /// @brief We get current size of the coefficients. Function updates this sizes
    /// accordingly to the refinement from coarser to finer grid. Assumption is
    /// that we are doing refinement from the index finest_low till the index
    /// highest_low.
    ///
    /// @param clevel coarse level
    /// @param flevel finer level
    /// @param finest_low "low index" of the support of a basis function (finest grid)
    /// @param finest_high "high index" of the support of a basis function (finest
    ///        grid)
    /// @param size_of_coefs size of the coefficients
    ///
    /// @return number of all new coefficients
    unsigned _updateSizeOfCoefs(const unsigned clevel,
                                const unsigned flevel,
                                const gsVector<index_t, d>& finest_low,
                                const gsVector<index_t, d>& finest_high,
                                gsVector<index_t, d>& size_of_coefs);

public:

  /// @brief Returns the dimension of the parameter space
  short_t domainDim() const { return d; }

    GISMO_CLONE_FUNCTION(gsTHBSplineBasis)

  /// @brief Prints the object as a string.
  std::ostream &print(std::ostream &os) const
  {
      os << "Truncated ";
      gsHTensorBasis<d,T>::printBasic(os);
      //this->printCharMatrix(os);
      return os;
  }

    GISMO_MAKE_GEOMETRY_NEW


  /**
   * @brief Returns the B-spline representation of a THB-spline subpatch.
   * @param b1  bottom left corner of the box (vector of indices with respect to the gsCompactKnotVector of the highest possible level)
   * @param b2  top right corner of the box (vector of indices with respect to the gsCompactKnotVector of the highest possible level)
   * @param level   level of the box
   * @param geom_coef control points of the THB-spline geometry
   * @param[out] cp control points of the B-spline patch
   * @param[out] k1 knot vector of the B-spline patch (first dimension)
   * @param[out] k2 knot vector of the B-spline patch (second dimension)
  */
  void getBsplinePatchGlobal(gsVector<index_t> b1, gsVector<index_t> b2, unsigned level, const gsMatrix<T>& geom_coef, gsMatrix<T>& cp, gsKnotVector<T>& k1, gsKnotVector<T>& k2) const;

  /**
   * @brief Return the list of B-spline patches to represent a THB-spline geometry.
   * @param geom_coef control points of the THB-spline geometry
   * @param[out] cp  control points of all B-spline patches stacked on top of each other
   * @param[out] b1  bottom left corners of the box (vector of indices with respect to the gsKnotVector of the highest possible level)
   * @param[out] b2  top right corners of the box (vector of indices with respect to the gsKnotVector of the highest possible level)
   * @param[out] level levels of the boxes (level[i]: level of the i-th box,)
   * @param[out] nvertices number of control points (nvertices[i,j]: number of control points in j-direction for the i-th box)
  */
  void getBsplinePatches(const gsMatrix<T>& geom_coef, gsMatrix<T>& cp, gsMatrix<index_t>& b1, gsMatrix<index_t>& b2, gsVector<index_t>& level, gsMatrix<index_t>& nvertices) const;

  /**
   * @brief Return a multipatch structure of B-splines
   * @param geom_coef control points of the THB-spline geometry
  */
  gsMultiPatch<T> getBsplinePatchesToMultiPatch(const gsMatrix<T>& geom_coef) const;

  /**
     * @brief Return the list of B-spline patches to represent a THB-spline geometry.
     * @param geom_coef control points of the THB-spline geometry
     * @param[out] cp  control points of all B-spline patches stacked on top of each other
     * @param[out] b1  bottom left corners of the box (vector of indices with respect to the gsKnotVector of the highest possible level)
     * @param[out] b2  top right corners of the box (vector of indices with respect to the gsKnotVector of the highest possible level)
     * @param[out] level levels of the boxes (level[i]: level of the i-th box,)
     * @param[out] nvertices number of control points (nvertices[i,j]: number of control points in j-direction for the i-th box)
     * @param[out] trim_curves the trimming curves for parasolid vector<connected_component<polylines<segments<T> > > > where the first polyline is the outer curve and the rest are holes
    */
    void getBsplinePatches_trimming(const gsMatrix<T>& geom_coef, gsMatrix<T>& cp, gsMatrix<index_t>& b1, gsMatrix<index_t>& b2, gsVector<index_t>& level, gsMatrix<index_t>& nvertices,
                           std::vector<std::vector<std::vector< std::vector<T> > > >& trim_curves) const;

    /**
       * @brief Return a multipatch structure of B-splines
       * @param geom_coef control points of the THB-spline geometry
       * @param[out] trim_curves the trimming curves for parasolid vector<connected_component<polylines<segments<T> > > > where the first polyline is the outer curve and the rest are holes
      */
    gsMultiPatch<T> getBsplinePatchesToMultiPatch_trimming(
            const gsMatrix<T>& geom_coef,
            std::vector<std::vector<std::vector< std::vector<T> > > >& trim_curves) const;

    /*
       * @brief Return the connected components of domain levels in knot vector indices (the boundary of a CC and the holes in the corresponding component)
       * @param[out] level levels of the boxes (level[i]: level of the i-th box,)
       * @param[out] connectedComponents the connected components in format vector<connected_component<polylines<segments<index_t> > > > where the first polyline is the outer curve and the rest are holes
    */
    void getConnectedComponents(std::vector<std::vector<std::vector< std::vector<index_t> > > >& connectedComponents, gsVector<index_t>& level) const;


   /// @brief returns transfer matrices betweend the levels of the given hierarchical spline
   void transferbyLvl (std::vector<gsSparseMatrix<T> >& result);

    /// @brief Decomposes domain of the THB-Spline-Basis into partitions.
    ///
    /// Each partiotion describes an area of the same level, this area is contained in 
    /// a bounding box (boundaryAABB) and has its own trimming curves (trimmCurves).
    /// Look above for the definition of types (AxisAlignedBoundingBox, TrimmingCurves).
    ///
    /// @param[out] boundaryAABB axis aligned bounding boxes -- "each" trimCurve has it own box
    /// @param[out] trimCurves trimming curves, each trimming curve describes a region of the same level
    void decomposeDomain(typename gsTHBSplineBasis::AxisAlignedBoundingBox& boundaryAABB,
                         typename gsTHBSplineBasis::TrimmingCurves& trimCurves) const;

    /// @brief Returns a tensor B-Spline patch defined by boundingBox.
    ///
    /// The B-Spline patch knots are the same as the THB-Spline-Basis knots from the input
    /// level. Geometry of the patch is defined via input coefficients.
    gsTensorBSpline<d, T>
    getBSplinePatch(const std::vector<index_t>& boundingBox,
                    const unsigned level,
                    const gsMatrix<T>& geomCoefs) const;

private:
    /**
     * @brief Initialize the characteristic and coefficient
     * matrices and the internal bspline representations.
    **/
    void update_structure() 
    {
        gsHTensorBasis<d,T>::update_structure(); 
        representBasis();
    }

    /**
      @brief Returns a representation of \a thbCoefs as tensor-product
      B-spline coefficientes \a lvlCoefs at level \a level.

      The representation is valid only at areas of the parameter
      domain with level less or equal to \a level

      @param[in] thbCoefs The input coefficients corresponding to basis function in this THB
      @param[in] level the level of tensor-product basis to be computed
      @param[out] lvlCoefs coefficients in tensor-product basis of level \a level
     */
    // todo: rename as: representAtLevel
    void globalRefinement(const gsMatrix<T> & thbCoefs, int level, 
                          gsMatrix<T> & lvlCoefs) const;

    gsSparseMatrix<T> coarsening(const std::vector<gsSortedVector<index_t> >& old,
                           const std::vector<gsSortedVector<index_t> >& n,
                           const gsSparseMatrix<T,RowMajor> & transfer) const;

    gsSparseMatrix<T> coarsening_direct( const std::vector<gsSortedVector<index_t> >& old,
                                   const std::vector<gsSortedVector<index_t> >& n, 
                                   const std::vector<gsSparseMatrix<T,RowMajor> >& transfer) const;

    gsSparseMatrix<T> coarsening_direct2( const std::vector<gsSortedVector<index_t> >& old,
                                   const std::vector<gsSortedVector<index_t> >& n,
                                   const std::vector<gsSparseMatrix<T,RowMajor> >& transfer) const;
    

    // ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    // Utility functions for decomposeDomain
    // ......................................................................
    
    /// @brief Checks if the first box is completely inside second box.
    /// 
    /// Utility function for the decomposeDomain member function.
    bool isFirstBoxCompletelyInsideSecond(const std::vector<index_t>& firstBox,
					  const std::vector<index_t>& secondBox) const
    {
	return secondBox[0] < firstBox[0] && secondBox[1] < firstBox[1] &&
	       firstBox[2] < secondBox[2] && firstBox[3] < secondBox[3];   
    }
    

    /// @brief Checks if the boxes are the same
    ///
    /// Utility function for the decomposeDomain member function.
    bool areBoxesTheSame(const std::vector<index_t>& firstBox,
			 const std::vector<index_t>& secondBox) const
    {
	return firstBox[0] == secondBox[0] && firstBox[1] == secondBox[1] &&
  	       firstBox[2] == secondBox[2] && firstBox[3] == secondBox[3];
    }
    
    /// @brief Breaks the cycles of polylines and returns updated polylines.
    ///
    /// domainBoundaryParams can return polylines with cycles, this function 
    /// decomposes the cycles and returns polylines without cycles. Look 
    /// above for the definition of types (AxisAlignedBoundingBox, TrimmingCurves).
    ///
    /// @param aabb axis aligned bounding boxes each polyline has it own box
    /// @param polylines polyline describing area with the same level
    void breakCycles(typename gsTHBSplineBasis::AxisAlignedBoundingBox& aabb,
		     typename gsTHBSplineBasis::Polylines& polylines) const;


    /// @brief Identify if the polyline can be split into two cycles.
    ///
    /// @param[in] polyline description of polyline
    /// @param[out] pt the point where two cycles meet (if there are two cycles)
    /// @return the index of the segment with the point 
    index_t identifyCycle(const std::vector< std::vector< T> >& polyline,
			  std::pair<T, T>& pt) const;

    
    /// @brief Breaks polyline into two parts.
    ///
    /// Function splits input polyline (with al least two cycles) at given point 
    /// into two cycles (part1 and part2). 
    ///
    /// @param[in] polyline input polyline with at least two cycles
    /// @param[in] segment the index of the segment of the polyline where the point is
    /// @param[in] pt the point where two cyles in polyline meet
    /// @param[out] part1 first cycle of polyline
    /// @param[out] part2 second cycle of polyline
    void breakPolylineIntoTwoParts(const std::vector< std::vector<T> >& polyline,
				   const index_t segment, 
				   const std::pair<T, T>& pt,
				   std::vector< std::vector<T> >& part1,
				   std::vector< std::vector<T> >& part2) const;


    /// @brief Finds new axis aligned bounding box for given polyline.
    void findNewAABB(const std::vector< std::vector<T> >& polyline,
		     std::vector<index_t>& aabb) const;

	

private:

    // m_is_truncated(j)
    // if -1   : j-th basis function is not truncated,
    // if  \em k>0 : j-th basis function is truncated and its representation level is \em k
    gsVector<int> m_is_truncated;


    // m_presentation[j] is presentation of the j-th basis function in terms of
    // B-Splines at level m_is_truncated[j]
    //
    // if m_is_truncated[j] is equal to -1, then there is no entry
    // m_presentation[j]
    std::map<index_t, gsSparseVector<T> > m_presentation;

    using gsHTensorBasis<d,T>::m_bases;
    using gsHTensorBasis<d,T>::m_xmatrix;
    using gsHTensorBasis<d,T>::m_xmatrix_offset;
    using gsHTensorBasis<d,T>::m_deg;
};
/**
 * End of class gsTHBSplineBasis definition
 */


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTHBSplineBasis.hpp)
#endif

