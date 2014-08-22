
#pragma once

#include <vector>
//#include <set>
#include <map>
//#include <unordered_map> come on, lets start using C++11

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
    void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result)const
    {
        gsMatrix<unsigned> indices;
        this->active_into(u, indices);

        const unsigned numDers = (d * (d + 1)) / 2;
        gsMatrix<T> res(numDers, 1); // result of deriv2Single_into

        result.resize(indices.rows() * numDers, u.cols());
        result.fill(0);

        for (int i = 0; i < indices.cols(); i++)
        {
            for (int j = 0; j < indices.rows(); j++)
            {
                unsigned index = indices(j, i);
                if (j != 0 && index == 0)
                    break;

                this->deriv2Single_into(index, u.col(i), res);

                result.block(j * numDers, i, numDers, 1) = res;
            }
        }
    }

    // Look at gsBasis class for documentatation
    void deriv2Single_into(unsigned i,
                           const gsMatrix<T>& u,
                           gsMatrix<T>& result) const;



    // Look at gsBasis class for documentation
    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {

        gsMatrix<unsigned> indices;
        this->active_into(u, indices);
        gsMatrix<T> res(d, 1);

        result.resize(indices.rows() * d, u.cols());
        result.fill(0);


        for (int i = 0; i < indices.cols(); i++)
        {
            for (int j = 0; j < indices.rows(); j++)
            {

                unsigned index = indices(j, i);
                if (j != 0 && index == 0)
                    break;

                this->derivSingle_into(index, u.col(i), res);

                //                std::cout << "u.col(i): \n" << u.col(i) << std::endl;

                //                std::cout << "res: " << res.rows() << " x "
                //                          << res.cols() << std::endl;

                result.block(j * d, i, d, 1) = res;

//                for (unsigned dim = 0; dim < d; dim++)
//                {
//                    result(dim + j * d, i) = res(dim);
//                }
            }
        }
    }

 
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

        const unsigned maxLvl = this->m_tree.max_ins_level + 1;
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

        const unsigned maxLvl = this->m_tree.max_ins_level + 1;
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

        const unsigned maxLvl = this->m_tree.max_ins_level + 1;
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
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
    {

        gsMatrix<unsigned> indices;
        gsMatrix<T> res(1, 1);
        this->active_into(u, indices);

        result.resize(indices.rows(), u.cols());
        result.fill(0);

        //std::cout << "before for loop" << std::endl;

        for (int i = 0; i < indices.cols(); i++)
        {
            //std::cout << i << " / " << indices.cols() << std::endl;

            for (int j = 0; j < indices.rows(); j++)
            {
                //std::cout << j << " / " << indices.rows() << std::endl;

                unsigned index = indices(j, i);
                if (j != 0 && index == 0)
                    break;

                //std::cout << "index: " << index << std::endl;

                evalSingle_into(index, u.col(i), res);

                //std::cout << "end evalSingle_into_v2" << std::endl;
                result(j, i) = res(0, 0);
            }
        }
    }


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





   ///// ///////////////////////////////////////////////////////////////////////
   ////  DEVELOPMENT TILL HERE
   /// ////////////////////////////////////////////////////////////////////////



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

  /**
   * @brief Insert the given boxes into the quadtree.
   * @param boxes matrix of size (dim x (num(boxes) * 2)) - subsequent columns in the matrix define a box in the parameter space
   */
  void refine(gsMatrix<T> const & boxes);

  /**
   * @brief Insert the given boxes into the quadtree.
   * @param boxes vector of size (num(boxes)*5) - each box is defined by 5
   * elements in boxes (level of the box followed by 4 coordinates in
   * index space which define the left bottom and right top corner of
   * the box)
   */
  void refine(std::vector<unsigned> const & boxes);

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
   std::vector<gsMatrix<T> > transferbyLvl ();


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

    gsMatrix<T> coarsening(std::vector<gsSortedVector<unsigned> > old,  std::vector<gsSortedVector<unsigned> > n, const gsSparseMatrix<T,RowMajor> & transfer);
    gsMatrix<T> coarsening_direct( std::vector<gsSortedVector<unsigned> > old,std::vector<gsSortedVector<unsigned> > n, const std::vector<gsSparseMatrix<T,RowMajor> >& transfer);

private:

    // =========================================================================
    // some experimantal variables - under development

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

    // end under development
    // =========================================================================


};
/**
 * End of class gsTHBSplineBasis definition
 */


} // namespace gismo

#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsTHBSplineBasis.hpp)
#endif

////representation of the basis function in the terms of functions in the highest level-initialization step
///**
// * @brief Returns the coefficients computed by Boehm algorithm.
// * @param u parameter values
// * @return coefficients obtained by knot insertion
// */
//std::vector<gsMatrix<T>  * > rep_on_finest_level(const gsMatrix<T> & u)const;

///**
// * @brief Returns the coefficients computed by Boehm algorithm (recursive part).
// * @param u parameter values
// * @param R coefficient matrix
// * @param level current level
// * @param dir_x shift in x direction
// * @param dir_y shift in y direction
// * @return coefficients obtained by knot insertion
// */
//gsMatrix<T> * rep_on_finest_level_recur(const gsMatrix<T> & u, gsMatrix<T> & R, int level, int dir_x, int dir_y)const;


//    /**
//     * @brief Evaluation of a B-spline basis function 1D.
//     * @param i number of the function
//     * @param parameter parameter value
//     * @param deg degree
//     * @param m_knots knot vector
//     * @return B-spline value
//     */
//    T eval_one_direction(unsigned &i, T &parameter,const int &deg,const gsCompactKnotVector<T> &m_knots) const;
//    /**
//     * @brief Derivative of a B-spline basis function on the highest level.
//     * @param function number of the THB-spline basis function in hid
//     * @param parameter parameter values
//     * @param coef coefficients obtained by knot insertion
//     * @param[out] derivs  B-spline derivative multiplied by \a coef
//     * @return B-spline derivative with respect to any direction
//     */
//    std::vector<T> hier_derivSingle_into(int function, gsVector<T> parameter, T const &coef,gsMatrix<T> & derivs)const;

//    /**
//     * @brief Derivative of a THB-spline basis function.
//     * @param function number of the THB-spline basis function in hid
//     * @param parameter parameter values
//     * @param[out] derivs THB-spline derivative
//     */
//    void hier_derivSingle_into(int function, gsVector<T> parameter, gsMatrix<T> & derivs)const;
//    /**
//     * @brief Second derivative of a B-spline basis function on the highest level.
//     * @param function number of the THB-spline basis function in hid
//     * @param parameter parameter values
//     * @param coef coefficients obtained by knot insertion
//     * @param[out] derivs B-spline second derivative multiplied by \a coef
//     * @return B-spline second derivative with respect to any direction
//     */
//    std::vector<T> hier_deriv2Single_into(int function, gsVector<T> parameter, T const &coef,gsMatrix<T> & derivs)const;

//    /**
//     * @brief Second derivative of a THB-spline basis function.
//     * @param function number of the THB-spline basis function in hid
//     * @param parameter parameter values
//     * @param[out] derivs THB-spline second derivative
//     */
//    void hier_deriv2Single_into(int function, gsVector<T> parameter, gsMatrix<T> & derivs)const;

//  /**
//   * @brief Evaluates a THB-spline basis function.
//   * @param i number of the function in cid numbering
//   * @param u parameter values
//   * @param[out] result THB-spline values
//   */
//  void evalSingle_into_v1(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

//  /**
//   * @brief Evaluates the derivative of a THB-spline basis function.
//   * @param i number of the function in cid numbering
//   * @param u parameter values
//   * @param[out] result derivative values
//  */
//  void derivSingle_into_v1(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result ) const;

//  /**
//   * @brief Evaluates the second derivative of a THB-spline basis function.
//   * @param i number of the function in cid numbering
//   * @param u parameter values
//   * @param[out] result second derivative values
//  */
//  void deriv2Single_into_v1(unsigned & i, const gsMatrix<T> & u, gsMatrix<T>& result ) const;


///// Evaluate the nonzero basis functions of the highest level at all columns of the matrix (or vector) u
//// TO DO: Better evaluate not the basis functions, but their linera combination, ie result is 1 x u.cols()
//void eval_surf_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
//{
//    gsMatrix<T> temp( 2,1) ;
//    std::vector<std::vector<T> > temp_output;//collects the outputs//TODO this is always d+1 times d+1, it is possible to set result directly
//    temp_output.resize(u.cols());//each column a point
//    result.resize((this->m_deg[0]+1)*(this->m_deg[1]+1),u.cols());
//    std::vector<int>spans;
//    spans.resize(this->m_bases[0]->dim() );
//    for(int i = 0; i < u.cols();i++){//for each point compute the values of the functions from all levels
//        int a = this->get_max_inserted_level();
//        //find the lower left function which act on the point u(i)
//        for(unsigned int j = 0; j < spans.size(); j++){
//            spans[j] = this->m_bases[a]->component(j).knots().Uniquefindspan(u(j,i));
//        }
//        for(unsigned int j = 0; j < spans.size(); j++){
//            spans[j] = this->m_bases[a]->component(j).knots().lastKnotIndex(spans[j]) - (this->m_deg[j]);
//        }
//        temp = u.col(i);
//        if(hl.size()<af){
//            std::vector<gsMatrix<T>  *> temp1  = rep_on_finest_level(u);
//            hl.push_back(temp1[0]);
//            while (temp1.size()-1 )
//            {
//                delete temp1.back();
//                temp1.pop_back();
//            }
//            for(int j = spans[0]; j <= spans[0]+this->m_deg[0]; j++){
//                for(int k = spans[1]; k <= spans[1]+this->m_deg[1]; k++){
//                    //call the evaluation of a single basis function
//                    //on thi highest level for all basis functions
//                    //on the point u(i)
//                    //repalce with tensor basis eval / coef??
//                    //result( (j-spans[0]) * (this->m_deg[1]+1) + (k-spans[1]),i) = hier_evalSingle_into(this->max_size[a] + k * (this->m_bases[a]->size(0)) + j , temp, (*hl[af-1])(j-spans[0],k-spans[1]) );
//                    gsMatrix<T> ev;
//                    this->m_bases[a]->evalSingle_into(k * (this->m_bases[a]->size(0)) + j, temp, ev);
//                    result( (j-spans[0]) * (this->m_deg[1]+1) + (k-spans[1]),i) = (*hl[af-1])(j-spans[0],k-spans[1]) * ev(0,0);//hier_evalSingle_into(this->max_size[a] + k * (this->m_bases[a]->size(0)) + j , temp, (*hl[af-1])(j-spans[0],k-spans[1]) );
//                }
//            }

//        }else{
//            for(int j = spans[0]; j <= spans[0]+this->m_deg[0]; j++){
//                for(int k = spans[1]; k <= spans[1]+this->m_deg[1]; k++){
//                    //call the evaluation of a single basis function
//                    //on thi highest level for all basis functions
//                    //on the point u(i)
//                    //result( (j-spans[0]) * (this->m_deg[1]+1) + (k-spans[1]),i) = hier_evalSingle_into(this->max_size[a] + k * (this->m_bases[a]->size(0)) + j , temp, (*hl[af-1])(j-spans[0],k-spans[1]) );
//                    gsMatrix<T> ev;
//                    this->m_bases[a]->evalSingle_into(k * (this->m_bases[a]->size(0)) + j, temp, ev);
//                    result( (j-spans[0]) * (this->m_deg[1]+1) + (k-spans[1]),i) = (*hl[af-1])(j-spans[0],k-spans[1]) * ev(0,0);
//                }
//            }

//        }
//    }
//}

///*
// * Computes all the active functions on a given point and for each
// * of them calls the evalSingle_into();
// */
//void eval_into_v1(const gsMatrix<T> & u, gsMatrix<T>& result) const {
//    gsVector<T,2> temp;
//    std::vector<std::vector<T> > temp_output;//collects the outputs
//    temp_output.resize(u.cols());//each column a point
//    gsMatrix<unsigned> act;
//    this->hier_active_into(u,act);
//    // TODO order the input so all points which are on the same cell
//    // in the highest level are evaluated after each other

//    std::vector<int> spans;
//    spans.resize(this->m_bases[0]->dim() );
//    for(int i = 0; i < u.cols();i++){
//    //for each point compute the values of the functions from all levels
//    //for all functions in act colums call eval for one function
//        temp[0] = u(0,i);
//        temp[1] = u(1,i);
//        int a = this->get_max_inserted_level();
//        //check if the prewious and the current evaluated points are in the same cell in the highest level
//        for(unsigned int j = 0; j < spans.size(); j++){
//            spans[j] = this->m_bases[a]->component(j).knots().Uniquefindspan(u(j,i));
//        }
//        for(unsigned int j = 0; j < spans.size(); j++){
//            spans[j] = this->m_bases[a]->component(j).knots().lastKnotIndex(spans[j]) - (this->m_deg[j]);
//        }
//        if( (ks[0] == spans[0])&&(ks[1] == spans[1]) ){ //if the cells are the same we do not have to recompute the coeffcients on the highest level
//            af = 0;
//        }else{
//          af =0;
//          //deleteing the highest level coeficients- the current and
//          //the previous point are nont in the same cell
//          freeAll( hl );

//          ks[0] = spans[0];
//          ks[1] = spans[1];
//        }
//        for(index_t j = 0; j < act.rows();j++){
//            if( act(j,i) != this->undefined_value ){
//                if(this->get_max_inserted_level()==0){
//                    // if max_inserted_level is 0, then call the B-spline evaluation
//                    //int level  = this->get_level(function);
//                    //unsigned k = function - this->max_size[level]; // go back to the tensor-product id
//                    // Evaluate the tensor-basis
//                    gsMatrix<T> ev;
//                    this->m_bases[0]->evalSingle_into(act(j,i), temp, ev);
//                    temp_output[i].push_back(ev(0,0));
//                    //temp_output[i].push_back(hier_evalSingle_into(act(j,i),temp,1));
//                }else{
//                    // otherwise, standard THB-representation on the finest level
//                    af++;
//                    //temp_output[i].push_back(hier_evalSingle_into(act(j,i),temp));
//                    gsMatrix<T> ev(1,1);
//                    evalSingle_into(this->hidTocid(act(j,i)),temp,ev);
//                    temp_output[i].push_back(ev(0,0));
//                }
//            }
//        }
//    }
//    this->copy_to_matrix(temp_output,result);
//    temp_output.clear();
//    freeAll( hl );
//    af = 0;
//    ks[0] =-1;
//    ks[1] =-1;
//}

///**
// * @brief Computes the derivatives of basis functions of the highest level for all points in the matrix \a u.
// * @param u
// * @param[out] result
// */
//void deriv_surf_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const{
//    //vector<vector<T> > temp_output;//collects the outputs//TODO this is always d+1 times d+1, it is possible to set result directly
//    //temp_output.resize(u.cols());//each column a point
//    int directions = 2;
//    gsMatrix<T> temp_output;  //temporary matrix of coefficients
//    result.resize((this->m_deg[0]+1)*(this->m_deg[1]+1)*directions,u.cols());

//    //vector<gsMatrix<T>  *> temp1  = rep_on_finest_level(u);
//    std::vector<int> spans;
//    spans.resize(this->m_bases[0]->dim() );
//    for(int i = 0; i < u.cols();i++){//for each point compute the values of the functions from all levels
//        int a = this->get_max_inserted_level();
//        //finthe lower left basis function which act on the point u(i)
//        for(unsigned int j = 0; j < spans.size(); j++){
//            spans[j] = this->m_bases[a]->component(j).knots().Uniquefindspan(u(j,i));
//        }
//        for(unsigned int j = 0; j < spans.size(); j++){
//            spans[j] = this->m_bases[a]->component(j).knots().lastKnotIndex(spans[j]) - (this->m_deg[j]);
//        }
//        gsMatrix<T> temp = u.col(i);
//        int index = 0;
//        if(hl.size()<af){
//            std::vector<gsMatrix<T>  *> temp1  = rep_on_finest_level(u);
//            hl.push_back(temp1[0]);
//            for(int j = spans[0]; j <= spans[0]+this->m_deg[0]; j++){
//                for(int k = spans[1]; k <= spans[1]+this->m_deg[1]; k++){
//                    //call the derivative of a single basis function on tei highest level for all basis functions on the point u(i)

//                    hier_derivSingle_into(this->max_size[a] + k * (this->m_bases[a]->size(0)) + j , temp, (*hl[af-1])(j-spans[0],k-spans[1]), temp_output);
//                    result(index,i) = temp_output(0,0);
//                    result(index+1,i) = temp_output(0,1);
//                    index +=2;
//                }
//            }

//        }else{
//            for(int j = spans[0]; j <= spans[0]+this->m_deg[0]; j++){
//                for(int k = spans[1]; k <= spans[1]+this->m_deg[1]; k++){
//                    //call the derivative of a single basis function on tei highest level for all basis functions on the point u(i)
//                    hier_derivSingle_into(this->max_size[a] + k * (this->m_bases[a]->size(0)) + j , temp, (*hl[af-1])(j-spans[0],k-spans[1]), temp_output);
//                    result(index,i) = temp_output(0,0);
//                    result(index+1,i) = temp_output(0,1);
//                    index +=2;
//                }
//            }
//        }
//    }
//}

///**
// * @brief Computes the derivatives of all basis functions acting on poinst in \a u.
// * @param u
// * @param[out] result
// */
//void deriv_into_v1(const gsMatrix<T> & u, gsMatrix<T>& result) const {
//    gsVector<T,2> temp;
//    std::vector<std::vector<T> > temp_output;
//    gsMatrix<T> T_deriv;
//    temp_output.resize(u.cols());//each column a point
//    gsMatrix<unsigned> act;
//    this->hier_active_into(u,act);
//    std::vector<int> spans;
//    spans.resize(this->m_bases[0]->dim() );
//    for(int i = 0; i < u.cols();i++){//for each poin compute the values of the functions from all levels
//    //for all functions in act colums call eval for one function
//        temp[0] = u(0,i);
//        temp[1] = u(1,i);

//        int a = this->get_max_inserted_level();
//        for(unsigned int j = 0; j < spans.size(); j++){
//            spans[j] = this->m_bases[a]->component(j).knots().Uniquefindspan(u(j,i));
//        }
//        for(unsigned int j = 0; j < spans.size(); j++){
//            spans[j] = this->m_bases[a]->component(j).knots().lastKnotIndex(spans[j]) - (this->m_deg[j]);
//        }

//        if( (ks[0] == spans[0])&&(ks[1] == spans[1]) ){//here we test if we can use the highes level coefficient matrices from the previous derivative
//            af = 0;
//        }else{
//          af =0;
//          freeAll( hl );

//          ks[0] = spans[0];
//          ks[1] = spans[1];
//        }

//        for(index_t j = 0; j < act.rows();j++){
//            if( act(j,i) !=this->undefined_value ){
//                //void hier_derivSingle_into(int function, gsVector<T> parameter, gsMatrix<T> & derivs)const;
//                if(this->get_max_inserted_level()==0){
//                    hier_derivSingle_into(act(j,i),temp,1, T_deriv);
//                    temp_output[i].push_back( T_deriv(0,0) );
//                    temp_output[i].push_back( T_deriv(0,1) );
//                }else{
//                    af++;
//                    hier_derivSingle_into(act(j,i),temp, T_deriv);
//                    temp_output[i].push_back( T_deriv(0,0) );
//                    temp_output[i].push_back( T_deriv(0,1) );
//                }
//            }
//        }
//    }
//    this->copy_to_matrix(temp_output,result);
//    temp_output.clear();

//    freeAll( hl );

//    af = 0;
//    ks[0] =-1;
//    ks[1] =-1;
//}

//void deriv2_into_v1(const gsMatrix<T> & u, gsMatrix<T>& result) const {
//    std::vector<std::vector<T> > temp_output;
//    std::vector<std::vector<T> > temp_output_xy;
//    std::vector<unsigned int> t;
//    gsMatrix<T> T_deriv;
//    temp_output.resize(u.cols());//each column a point
//    temp_output_xy.resize(u.cols());//each column a point
//    int level, function;

//    gsMatrix<unsigned> act;
//    this->hier_active_into(u,act);
//    for(int i = 0; i < u.cols();i++){//for each poin compute the values of the functions from all levels
//    //for all functions in act colums call deriv2 for one function
//        for(index_t j = 0; j < act.rows();j++){
//            if( act(j,i) !=this->undefined_value ){
//                function = act(j,i);
//                level = this->get_level(function);
//                t = this->get_x_y(function, level);
//                //this->m_cmatrix[level][ this->fromPair(t[0], t[1], level) ] = 1;
//                this->m_cmatrix[level][ this->m_bases[level]->index(t[0],t[1]) ] = 1;
//                deriv_surf_into2(u.col(i), T_deriv);
//                for(int k = 3; k < T_deriv.rows();k++ ){
//                    if ((k%3)==0){
//                        T_deriv(0,0) +=T_deriv(k,0);
//                    }else{
//                        if((k%3)==1){
//                          T_deriv(1,0) +=T_deriv(k,0);
//                        }else{
//                          T_deriv(2,0) +=T_deriv(k,0);
//                        }

//                    }
//                }

//                temp_output[i].push_back( T_deriv(0,0) );
//                temp_output[i].push_back( T_deriv(1,0) );
//                temp_output[i].push_back( T_deriv(2,0) );
//                T_deriv.resize(1,1);
//                //this->m_cmatrix[level][ this->fromPair(t[0], t[1], level) ] = 0;
//                this->m_cmatrix[level][ this->m_bases[level]->index(t[0],t[1]) ] = 0;
//            }
//        }
//    }
//    this->copy_to_matrix(temp_output, result);
//}


//void deriv_surf_into2(const gsMatrix<T> & u, gsMatrix<T>& result ) const{
//    int directions = 3;
//    gsMatrix<T> temp_output;  //temporary matrix of coefficients
//    result.resize((this->m_deg[0]+1)*(this->m_deg[1]+1)*directions,u.cols());
//    std::vector<gsMatrix<T>  *> temp1;
//    if(this->get_max_inserted_level() > 0){
//       temp1  = rep_on_finest_level(u);
//    }else{
//        temp1.push_back(new gsMatrix<T>());
//        temp1[0]->setOnes(this->m_deg[0]+1, this->m_deg[1]+1);
//    }

//    std::vector<unsigned> spans;
//    spans.resize(this->m_bases[0]->dim() );
//    for(int i = 0; i < u.cols();i++){//for each point compute the values of the functions from all levels
//        int a = this->get_max_inserted_level();
//        for(unsigned int j = 0; j < spans.size(); j++){
//            spans[j] = this->m_bases[a]->component(j).knots().Uniquefindspan(u(j,i));
//        }
//        for(unsigned int j = 0; j < spans.size(); j++){
//            spans[j] = this->m_bases[a]->component(j).knots().lastKnotIndex(spans[j]) - (this->m_deg[j]);
//        }
//        gsMatrix<T> temp = u.col(i);
//        int index = 0;

//        gsMatrix<T>  temp_x = this->m_bases[a]->component(0).deriv2(temp.row(0));
//        gsMatrix<T>  temp_y = this->m_bases[a]->component(1).deriv2( temp.row(1));
//        // gsMatrix<T> * temp_x= this->T_bsp.deriv2(temp.row(0));
//        // gsMatrix<T> * temp_y = this->T_bsp1.deriv2( temp.row(1));

//        for(unsigned j = spans[0]; j <= spans[0]+this->m_deg[0]; j++){
//            for(unsigned k = spans[1]; k <= spans[1]+this->m_deg[1]; k++){
//                //hier_deriv2Single_into(this->max_size[a] + j * (this->m_cvs[1][a].size()-(this->m_deg[1]+1)) + k , temp, (*temp1[i])(j-span_x,k-span_y), temp_output);
//                //result(index,i) = temp_output(0,0);
//                //result(index+1,i) = temp_output(0,1);
//                //result(index+2,i) = temp_output(0,2);
//                std::vector<T> dxy = hier_derivSingle_into(this->max_size[a] + k * (this->m_bases[a]->size(0)) + j , temp, (*temp1[i])(j-spans[0],k-spans[1]), temp_output);
//                result(index,i) = temp_x(j-spans[0],0) * eval_one_direction(k,temp(1,0),this->m_deg[1],this->m_bases[a]->component(1).knots())*(*temp1[i])(j-spans[0],k-spans[1]);//xx deriv;

//                result(index+1,i) = temp_y(k-spans[1],0) * eval_one_direction(j,temp(0,0),this->m_deg[0],this->m_bases[a]->component(0).knots())*(*temp1[i])(j-spans[0],k-spans[1]);//yy deriv

//                result(index+2,i) = dxy[0]*dxy[1]*(*temp1[i])(j-spans[0],k-spans[1]);
//                index +=directions;
//            }
//        }
//    }
//    freeAll( temp1 );
//}
