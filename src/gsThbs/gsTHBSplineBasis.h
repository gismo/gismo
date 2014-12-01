/** @file gsTHBSplineBasis.h

    @brief Provides declaration of THBSplineBasis class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#pragma once

#include <gsThbs/gsHTensorBasis.h>
#include <gsThbs/gsTHBSpline.h>


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
 */
template<unsigned d, class T>
class gsTHBSplineBasis : public gsHTensorBasis<d,T>
{
public:
    /// Associated geometry type
    typedef gsTHBSpline<d,T> GeometryType;
    
    typedef typename gsHTensorBasis<d,T>::CMatrix CMatrix;
    
    typedef memory::shared_ptr< gsTHBSplineBasis > Ptr;

    using gsHTensorBasis<d, T>::flatTensorIndexOf;

public:

    /// Constructor out of a 2D Tensor BSpline Basis
    gsTHBSplineBasis(gsTensorBSplineBasis<d,T> const&  tbasis, int nlevels = 10) : gsHTensorBasis<d,T>(tbasis, nlevels) {
        initialize();
    }

    gsTHBSplineBasis(gsTensorBSplineBasis<d,T> const&  tbasis, std::vector<unsigned> boxes, int nlevels = 10) : gsHTensorBasis<d,T>(tbasis, nlevels, boxes) {
        initialize();
    }

    gsTHBSplineBasis(gsTensorBSplineBasis<d,T> const&  tbasis, gsMatrix<T> const & boxes, int nlevels = 10) : gsHTensorBasis<d,T>(tbasis, nlevels, boxes) {
        initialize();
    }

    gsTHBSplineBasis( gsTensorBSplineBasis<d,T> const&  tbasis, gsMatrix<T> const & boxes, std::vector<unsigned int> levels, int nlevels = 10): gsHTensorBasis<d,T>(tbasis, nlevels, boxes, levels){
        initialize();
    }
    /// Constructor out of a tensor BSpline Basis
    gsTHBSplineBasis(gsBasis<T> const&  tbasis, int nlevels = 10)
        : gsHTensorBasis<d,T>(tbasis, nlevels)
    {
        // initialize(); // is done in the base constructor
    }

    ~gsTHBSplineBasis()
    {

    }

public:

    // Look at gsBasis class for documentation
    void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result)const;

    // Look at gsBasis class for documentatation
    void deriv2Single_into(unsigned i,
                           const gsMatrix<T>& u,
                           gsMatrix<T>& result) const;



    // Look at gsBasis class for documentation
    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
 
    // Look at gsBasis class for documentation
    void derivSingle_into(unsigned i,
                             const gsMatrix<T> & u,
                          gsMatrix<T>& result) const;

    // look at eval_into
    void fastEval_into(const gsMatrix<T>& u,
                       gsMatrix<T>& result) const
    {
        gsMatrix<unsigned> indices;
        this->active_into(u, indices);

        result.setZero(indices.rows(), u.cols());

        const unsigned maxLvl = this->m_tree.getMaxInsLevel() + 1;
        std::vector< gsMatrix<T> > tmpResults(maxLvl, gsMatrix<T>());
        std::vector< gsMatrix<unsigned> > tmpActive(maxLvl, gsMatrix<unsigned>());
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
                    unsigned flatTenIndx = this->flatTensorIndexOf(index, lvl);
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
                    const gsMatrix<unsigned>& active = tmpActive[lvl];


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
        gsMatrix<unsigned> indices;
        this->active_into(u, indices);

        result.setZero(indices.rows() * d, u.cols());

        const unsigned maxLvl = this->m_tree.getMaxInsLevel() + 1;
        std::vector< gsMatrix<T> > tmpDeriv( maxLvl, gsMatrix<T>());
        std::vector< gsMatrix<unsigned> > tmpActive(maxLvl, gsMatrix<unsigned>());
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
                    this->m_bases[lvl]->deriv_into(u.col(pt), tmpDeriv[lvl]);
                    this->m_bases[lvl]->active_into(u.col(pt), tmpActive[lvl]);
                    processed(lvl) = 1;
                }

                if (m_is_truncated[index] == -1)
                {
                    unsigned flatTenIndx = this->flatTensorIndexOf(index, lvl);
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
                    const gsMatrix<unsigned>& active = tmpActive[lvl];
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
        gsMatrix<unsigned> indices;
        this->active_into(u, indices);
        const unsigned numDers = (d * (d + 1)) / 2;

        result.setZero(indices.rows() * numDers, u.cols());

        const unsigned maxLvl = this->m_tree.getMaxInsLevel() + 1;
        std::vector< gsMatrix<T> > tmpDeriv2(maxLvl, gsMatrix<T>());
        std::vector< gsMatrix<unsigned> > tmpActive(maxLvl, gsMatrix<unsigned>());
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
                    unsigned flatTenIndx = this->flatTensorIndexOf(index, lvl);
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
                    const gsMatrix<unsigned>& active = tmpActive[lvl];
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


    /// Returns sparse representation of the i-th basis function.
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


    void evalSingle_into(unsigned i,
                         const gsMatrix<T>& u,
                         gsMatrix<T>& result) const;

    ///look into gsHTensorBasis
    void increaseMultiplicity(index_t lvl, int dir, T knotValue, int mult = 1);

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
    void representBasis();


    /// Computes representation of j-th basis function on pres_level and
    /// saves it.
    ///
    /// @param j     index of basis function
    /// @param pres_level levet at which we want to present j-th basis function
    /// @param low "low index" of support of j-th basis function (finest grid)
    /// @param high "high index" of support of j-th basis function (finest grid)
    void _representBasisFunction(const unsigned j,
                                const unsigned pres_level,
                                const gsVector<unsigned, d>& finest_low,
                                const gsVector<unsigned, d>& finest_high);



    /// Saves a presentation of the j-th basis function. Presentation is given
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
                                      const gsVector<unsigned, d>& act_size_of_coefs,
                                      const unsigned j,
                                      const unsigned pres_level,
                                      const gsVector<unsigned, d>& finest_low);



    /// Computes tensor index of a basis function on a finer level (new_level)
    /// which is presented with tensor index (index) at a coarse level (level).
    ///
    /// @param index tensor index of the basis function on a coarse level
    /// @param level coarse level
    /// @param fin_low "low index" of the support of the basis function
    /// @param new_level finer level
    ///
    /// @return global tensor index of a basis function on the finer level
    unsigned _basisFunIndexOnLevel(const gsVector<unsigned, d>& index,
                                   const unsigned level,
                                   const gsVector<unsigned, d>& fin_low,
                                   const unsigned new_level);



    /// Performs truncation.
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
                   const gsVector<unsigned, d>& act_size_of_coefs,
                   const gsVector<unsigned, d>& size_of_coefs,
                   const unsigned level,
                   const gsVector<unsigned, d>& bspl_vec_ti,
                   const unsigned bspl_vec_ti_level,
                   const gsVector<unsigned, d>& finest_low);


    /// We get current size of the coefficients. Function updates this sizes
    /// accordingly to the refinement from coarser to finer grid. Assumption is
    /// that we are doing refinement from the index finest_low till the index
    /// highest_low.
    ///
    /// @param clevel coarse level
    /// @param flevel finer level
    /// @param flow "low index" of the support of a basis function (finest grid)
    /// @param chigh "high index" of the support of a basis function (finest
    ///        grid)
    /// @param size_of_coefs size of the coefficients
    ///
    /// @return number of all new coefficients
    unsigned _updateSizeOfCoefs(const unsigned clevel,
                                const unsigned flevel,
                                const gsVector<unsigned, d>& finest_low,
                                const gsVector<unsigned, d>& finest_high,
                                gsVector<unsigned, d>& size_of_coefs);

public:

  /// Returns the dimension of the parameter space
  int dim() const 
  { return d; }



  virtual gsTHBSplineBasis* clone() const
    { return new gsTHBSplineBasis(*this); }

  /// Prints the object as a string.
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
  void getBsplinePatchGlobal(gsVector<unsigned> b1, gsVector<unsigned> b2, unsigned level, const gsMatrix<T>& geom_coef, gsMatrix<T>& cp, gsCompactKnotVector<T>& k1, gsCompactKnotVector<T>& k2) const;  

  /**
   * @brief Return the list of B-spline patches to represent a THB-spline geometry.
   * @param geom_coef control points of the THB-spline geometry
   * @param[out] cp  control points of all B-spline patches stacked on top of each other
   * @param[out] b1  bottom left corners of the box (vector of indices with respect to the gsCompactKnotVector of the highest possible level)
   * @param[out] b2  top right corners of the box (vector of indices with respect to the gsCompactKnotVector of the highest possible level)
   * @param[out] level levels of the boxes (level[i]: level of the i-th box,)
   * @param[out] nvertices number of control points (nvertices[i,j]: number of control points in j-direction for the i-th box)
  */
  void getBsplinePatches(const gsMatrix<T>& geom_coef, gsMatrix<T>& cp, gsMatrix<unsigned>& b1, gsMatrix<unsigned>& b2, gsVector<unsigned>& level, gsMatrix<unsigned>& nvertices) const;

  /**
   * @brief Return a multipatch structure of B-splines
   * @param geom_coef control points of the THB-spline geometry
  */
  gsMultiPatch<T> getBsplinePatchesToMultiPatch(const gsMatrix<T>& geom_coef) const;

  /**
     * @brief Return the list of B-spline patches to represent a THB-spline geometry.
     * @param geom_coef control points of the THB-spline geometry
     * @param[out] cp  control points of all B-spline patches stacked on top of each other
     * @param[out] b1  bottom left corners of the box (vector of indices with respect to the gsCompactKnotVector of the highest possible level)
     * @param[out] b2  top right corners of the box (vector of indices with respect to the gsCompactKnotVector of the highest possible level)
     * @param[out] level levels of the boxes (level[i]: level of the i-th box,)
     * @param[out] nvertices number of control points (nvertices[i,j]: number of control points in j-direction for the i-th box)
     * @param[out] trim_curves the trimming curves for parasolid vector<connected_component<polylines<segments<T> > > > where the first polyline is the outer curve and the rest are holes
    */
    void getBsplinePatches_trimming(const gsMatrix<T>& geom_coef, gsMatrix<T>& cp, gsMatrix<unsigned>& b1, gsMatrix<unsigned>& b2, gsVector<unsigned>& level, gsMatrix<unsigned>& nvertices,
                           std::vector<std::vector<std::vector< std::vector<T> > > >& trim_curves) const;

    /**
       * @brief Return a multipatch structure of B-splines
       * @param geom_coef control points of the THB-spline geometry
       * @param[out] trim_curves the trimming curves for parasolid vector<connected_component<polylines<segments<T> > > > where the first polyline is the outer curve and the rest are holes
      */
    gsMultiPatch<T> getBsplinePatchesToMultiPatch_trimming(
            const gsMatrix<T>& geom_coef,
            std::vector<std::vector<std::vector< std::vector<T> > > >& trim_curves) const;

    /**
       * @brief Return the connected components of domain levels in knot vector indices (the boundary of a CC and the holes in the corresponding component)
       * @param[out] level levels of the boxes (level[i]: level of the i-th box,)
       * @param[out] connectedComponents the connected components in format vector<connected_component<polylines<segments<unsigned int> > > > where the first polyline is the outer curve and the rest are holes
      */
    void getConnectedComponents(std::vector<std::vector<std::vector< std::vector<unsigned int> > > >& connectedComponents, gsVector<unsigned>& level) const;



  // see gsBasis.h for documentation
  /*
   * @brief Insert the given boxes into the quadtree.
   * @param boxes matrix of size (dim x (num(boxes) * 2)) - subsequent columns in the matrix define a box in the parameter space
   */
  void refine(gsMatrix<T> const & boxes);

  // see gsHTensorBasis.h for documentation
  void refineElements(std::vector<unsigned> const & boxes);

  // see gsBasis.h for documentation
  void uniformRefine(int numKnots = 1);

  /**
   * @brief Initializes the \a m_cmatrix with 0 for evaluation of basis functions.
   */
  void update_cmatrix(std::vector< std::map<unsigned,T> > & m_cmatrix) const;

  /**
   * @brief Initialize the m_cmatrix up to \a c_level
   * with the \a geom_coeff of the geometry in the
   * direction specified by \a col.
   *
   * @param geom_coeff control points of the geometry
   * @param col direction (0, 1, 2 = x, y, z)
   * @param c_level the maximum level of interest
   */
   void update_cmatrix(const gsMatrix<T>& geom_coeff, int col, int c_level, 
                       std::vector< std::map<unsigned,T> > & m_cmatrix) const;

   ///returns transfer matrices betweend the levels of the given hierarchical spline
   void transferbyLvl (std::vector<gsMatrix<T> >& result);


private:
    /**
     * @brief Initialize the characteristic and coefficient
     * matrices and the internal bspline representations.
    **/
    void initialize();

    /**
     * @brief Returns the coefficients computed by Boehm algorithm (called by \ref getBsplinePatchGlobal).
     * @param level maximum refinement level
     * @param[out] coeffs coefficients obtained by knot insertion
     */
    void globalRefinement(int level, gsMatrix<T>& coeffs,
                          std::vector< std::map<unsigned,T> > & m_cmatrix ) const;

    /**
     * @brief return_cp_1D converts the coefficient matrix mat in the given direction to a column of the control points matrix
     * @param mat coefficient matrix
     * @param direction direction (0, 1, 2 = x, y, z)
     * @param cp control points
     */
    void return_cp_1D(const gsMatrix<T> & mat, int direction, gsMatrix<T>& cp)const;

    gsMatrix<T> coarsening(const std::vector<gsSortedVector<unsigned> >& old,
                           const std::vector<gsSortedVector<unsigned> >& n,
                           const gsSparseMatrix<T,RowMajor> & transfer);

    gsMatrix<T> coarsening_direct( const std::vector<gsSortedVector<unsigned> >& old,
                                   const std::vector<gsSortedVector<unsigned> >& n, 
                                   const std::vector<gsSparseMatrix<T,RowMajor> >& transfer);

private:

    // m_is_truncated(j)
    // if -1: j-th basis function is not truncated,
    // if  1: j-th basis function is truncated and its representation is in 1 level
    // if  2: j-th basis function is truncated and its representation is in 2 level
    // ...
    gsVector<int> m_is_truncated;


    // m_presentation[j] is presentation of the j-th basis function in terms of
    // B-Splines at level m_is_truncated[j]
    //
    // if m_is_truncated[j] is equal to -1, then there is no entry
    // m_presentation[j]
    std::map<unsigned, gsSparseVector<T> > m_presentation;
};
/**
 * End of class gsTHBSplineBasis definition
 */


} // namespace gismo

#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsTHBSplineBasis.hpp)
#endif

