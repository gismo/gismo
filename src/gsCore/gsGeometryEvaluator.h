/** @file gsGeometryEvaluator.h

    @brief Provides declaration of GeometryEvaluator interface.
    
    This file is part of the G+Smo library-
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, C. Hofreither, J. Sogn, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsBoundary.h>

namespace gismo
{


/**
    \brief interface to compute a geometry on sets of points

    the idea of greometry evaluators is that they act as a local
    cache for computed values. In this way it is possible to pass
    a unique argument to functions containing all the required values.

    All accessors to the stored values come in two versions: plural and
    singular. The plural returns a matrix with the computed values on all
    points, the singular takes as additional parameter that represent the
    point index and returns the associated subblock.

    \tparam T the coefficient type
**/
template <typename T>
class gsGeometryEvaluator
{
public:

    /**
       \brief Constructor using flags that define what should be
       evaluated, the parametric dimension, and geometric dimension.
       
       See gismo::gsNeedEnum for available flags.
    */
    gsGeometryEvaluator(unsigned flags, unsigned parDim, unsigned geoDim)
    : m_flags(flags), m_parDim(parDim), m_geoDim(geoDim)
    {}
    
    virtual ~gsGeometryEvaluator() {}

public:
    
    /**
       \brief change the values to compute

        This interface is provided in such a way that it is possible
        split the logic of the inizialization of the evaluator.
        (?)
    **/
    void addFlags (unsigned newFlags) {this->setFlags(m_flags|newFlags);}

    /**
        \brief set the values to compute

        This interface is provided in such a way that it is possible
        split the logic ot the inizialization of the evaluator.
    **/
    virtual void setFlags (unsigned newFlags) {m_flags = newFlags;}

    /**
        \brief returns which values are computed
    **/
    unsigned getFlags () {return m_flags;}

public:

    /**
       \brief Evaluate the geometry's quantities (specified by flags) at points \em u.
       
        The flags which are passed with the constructor specify which quantities
        should be computed (e.g., Jacobian, determinant of Jacobian, ...).\n
        By calling <em>evaluateAt(u)</em>, these quantities are computed at points \em u and stored. They
        can then be accessed by the corresponding accessors.
       
        \param[in] u gsMatrix of size <em>d</em> x \em N, where\n
        \em d is the dimension of the parameter domain, and\n
        \em N is the number of evaluation points.\n
        Each column of \em u corresponds to one evaluation point.
    **/
    virtual void evaluateAt(const gsMatrix<T>& u) = 0;

    /// Get the physical coordinates (i.e., the image of the
    /// evaluation points in the physical space).
    const gsMatrix<T>& values() const
    { GISMO_ASSERT(m_flags & NEED_VALUE, "Geometry values not computed"); return m_values; }

    /// Get the physical coordinates (i.e., the image of the
    /// evaluation points in the physical space).
    const typename gsMatrix<T>::constColumn  value(index_t k) const
    { GISMO_ASSERT(m_flags & NEED_VALUE, "Geometry values not computed"); return m_values.col(k); }

    /**
      \brief Get the geometry Jacobians.
    
      This function returns a matrix which contains the Jacobians for each evaluation point
      "next to each other".\n
    
      Let
      \f[ f:\mathbb R^2 \to \mathbb R^3, \quad f = ( f_1(s,t), f_2(s,t), f_3(s,t)).\f]
      Let the evaluation be at \f$N\f$ points \f$u_1,u_2,\ldots,u_N\f$.
      Then, the returned Jacobian \em J is a <em>3</em> x <em>2N</em> matrix of the following entries:
      \f[
      J = \left(
      \begin{array}{ccccc}
      \partial_s f_1(u_1) & \partial_t f_1(u_1) & \partial_s f_1(u_2) & \ldots & \partial_t f_1(u_N)\\
      \partial_s f_2(u_1) & \partial_t f_2(u_1) & \partial_s f_2(u_2) & \ldots & \partial_t f_1(u_N)\\
      \partial_s f_3(u_1) & \partial_t f_3(u_1) & \partial_s f_3(u_2) & \ldots & \partial_t f_1(u_N)
      \end{array}
      \right)
      \f]
    
      \param[out] J gsMatrix of size <em>physDim</em> x <em>parDim*N</em>, where\n
      \em physDim is the dimension of the physical domain, and\n
      \em parDim is the dimension of the parameter domain, and\n
      \em N is the number of evaluation points.
    */
    const gsMatrix<T>& jacobians() const
    { GISMO_ASSERT(m_flags & NEED_JACOBIAN, "Jacobians not computed"); return m_jacobians; }

    /// Returns the Jacobian at point \a k
    const typename gsMatrix<T>::constColumns  jacobian(index_t k) const
    {
        GISMO_ASSERT(m_flags & NEED_JACOBIAN, "Jacobians not computed");
        return m_jacobians.middleCols(k*m_parDim, m_parDim);
    }

    /**
        \brief get the pull back of the destination measure

        this is the weight in the integral on the parametric domain
        that is caused by the change of variables.

        if the dimension ot the parametric domain equals that of the target
        space this is the determinant of the jacobian f$\det J^f$.
        If the geometry is a curve it is the modulus of the normal.
        If it is a surface in R^3 it is the length of the normal.
        In general it can always be computed as f$\det(J^t J)^{1/2}f$
    **/
    const gsVector<T>&  measures() const
    { GISMO_ASSERT(m_flags & NEED_MEASURE, "det(J) not computed"); return m_measures; }

    /// Returns the measure at point \a k
    inline T  measure(index_t k) const
    { GISMO_ASSERT(m_flags & NEED_MEASURE, "det(J) not computed"); return m_measures(k); }

    /// Returns the (signed) Jacobian determinant at point \a k
    inline T  jacDet(index_t k) const
    { 
        GISMO_ASSERT(m_flags & NEED_MEASURE, "det(J) not computed"); 
        return m_orientation * m_measures(k); 
    }

    /**
        \brief Returns the orientation of the geometry.

        This is always 1 if the codimension is not 0.
        This returns garbage if the geometry is not oriented.

        The orientation is equal to the sign of the Jacobian
        determinant. It is assumed that this sign is the same for the
        whole of the domain. It is computed on a generic point inside
        the parameter domaim.

        For instance, the Jacobian determinant at point \a k can be
        retrieved as
        
        this->orientation() * this->measure(k);

        since this->measure(k) returns the absolute value of the
        Jacobian determinant in the co-dimension zero case.

    **/
    int orientation() const
    { 
        return m_orientation;
    }

    /** \brief Returns the second derivaties.
       
       Example for 2D case: Column k contains all the second
       derivatives at point \ k, stored as one big column:

       f1_uu
       f1_vv
       f1_uv
       f2_uu
       f2_vv
       f2_uv
       
       Note that the segment of the vector that refers to f1 has length (ParDim + (ParDim*(ParDim-1))/2).

    */
    const gsMatrix<T>& derivs2() const
    { 
        GISMO_ASSERT(m_flags & NEED_2ND_DER, "2nd derivatives not computed");
        return m_2ndDers; 
    }

    /// Returns the second derivatives (1D only)
    const typename gsMatrix<T>::constColumn  deriv2(index_t k) const
    {
        GISMO_ASSERT(m_flags & NEED_2ND_DER, "2nd derivatives not computed");
        return m_2ndDers.col(k);
    }

    /** \brief Returns the transformation of a gradient through this geometry mapping.
    *
    * In the case that the parameter domain and the physical domain have the same dimension,
    * this function returns a matrix which contains the transposed inverse of the Jacobians
    * for each evaluation point
    * "next to each other".\n
    *
    * Let \f$ f:\mathbb R^d \to \mathbb R^d\f$,
    * let the evaluation be at \f$N\f$ points \f$u_1,u_2,\ldots,u_N\f$, and let the
    * Jacobian w.r.t. point \f$ u_i \f$ be denoted by \f$ J_i \f$.
    * Then, the returned matrix \em invJ is a <em>d</em> x <em>d*N</em> matrix of the form
    * \f[
    * invJ = \left(
    * \begin{array}{cccc}
    * J_1^{-T} & J_2^{-T} & \ldots & J_N^{-T}
    * \end{array}
    * \right)
    * \f]
    *
    * If the parameter and the physical domains are not of equal
    * dimension, the result is a non-square matrix taking parameter
    * gradients to physical gradient.
    *
    * \return gsMatrix of size <em>d</em> x <em>d*N</em>, where\n    * 
    * \em d is the dimension of the parameter domain, and\n
    * \em N is the number of evaluation points.\n
    *
    * In the co-dimension zero case, the <em>k</em>-th <em>d</em> x
    * <em>d</em>-block corresponds to the transposed inverse of the
    * Jacobian at point \f$u_k\f$.
    *
    */
    const gsMatrix<T>& gradTransforms() const
    { 
        GISMO_ASSERT(m_flags & NEED_GRAD_TRANSFORM, "J^-1 not computed");
        return m_jacInvs; 
    }

    /**
    * \brief Returns the transform matrix (from parameter to physical
    * domain) of the gradient at point \a k
    *
    * This is equal to the transposed inverse of the Jacobian matrix
    * for the co-dimension zero case.
    */
    const typename gsMatrix<T>::constColumns gradTransform(index_t k) const
    { 
        GISMO_ASSERT(m_flags & NEED_GRAD_TRANSFORM, "J^-1 not computed");
        return m_jacInvs.middleCols( k*m_parDim, m_parDim);
    }

public:


    /**
      \brief Transforms parametric gradients to a gradients on the physical domain.
    
      The gradient information on the parameter domain at a certain points
      is given in \em <b>allGrads</b> in the following
      format:\n
      Each column of \em allGrads corresponds to one evaluation point. In this column, the
      gradients of all active (i.e., non-zero) basis functions are stored "one above the other".\n
      Example: Let \f$B_i(s,t), i = 1,...,9\f$ be a set of bivariate basis functions.
      Then, a column of \em allGrads reads\n
      \f[
      ( \partial_s B_1, \partial_t B_1, \partial_s B_2, \partial_t B_2, \ldots, \partial_t B_9 )^T.
      \f]
      Hence, \em <b>allGrads</b> is a gsMatrix of size <em>(n * ParDim)</em> x <em>K</em>,
      where \n
      ...\em n denotes the number of active functions, \n
      ...<em>ParDim</em> denotes the dimension of the parameter domain, and\n
      ...\em K denotes the number of columns of \em allGrads.\n
      \n
      These gradients in the <em>k</em>-th column of \em allGrads
      are transformed to the physical domain and stored in \em <b>trfGradsK</b>.
      \em trfGradsK is of size <em>ParDim</em> x <em>n</em>. The <em>i</em>-th column
      of \em trfGradsK corresponds to the transformed gradient of the <em>i</em>-th function.\n
      \n
      Example: Let \f$x\f$ and \f$y\f$ be the coordinates in the phyiscal domain. Then,
      in the above example, \em <b>trfGradsK</b> is given by\n
      \f[
      \left(\begin{array}{cccc}
      \partial_x B_1 & \partial_x B_2 & \ldots & \partial_x B_9 \\
      \partial_y B_1 & \partial_y B_2 & \ldots & \partial_y B_9
      \end{array}\right)
      \f]
    
    
      \param k Indicates which column of \em allGrads should be transformed.
      \param[in] allGrads gsMatrix containing computed gradients in the format
      described above.
      \param[out] trfGradsK gsMatrix with the corresponding gradients on the
      physical domain in the format as described above.
    */
    // rename to transformGradsHgrad
    virtual void transformGradients(index_t k, const gsMatrix<T>& allGrads, gsMatrix<T>& trfGradsK)  const = 0;

    /**
      \brief Transforms parametric basis values at point \a k of a vector field, while perserving the divergence (Piola Transformation).
      NOT TESTED
    */
    // allValues: NumAct X TargetDim
    virtual void transformValuesHdiv(index_t k, const gsMatrix<T>& allValues, 
                                    gsMatrix<T>  & trfValues) const = 0;

    /**
      \brief Transforms parametric gradients to a gradients on the physical domain of a vector field while perserving the divergence.

      The transformation is sometimes called the contravariant Piola transform:
      \f[
        \mathcal{F}^{\text{div}}(\hat{\mathbf{u}}) =\frac{1}{\text{det}\, J} J\,\hat{\mathbf{u}} \circ F^{-1},
      \f]
      where \f$F\f$ is mapping from \f$\hat{\Omega}\f$ to \f$\Omega\f$ and \f$J\f$ 
      is the jacobian.

      The gradient information on the parameter domain at a certain points for a
      vector field is given in \em <b>allGrads_vec</b> in the following format:\n

      Each component of the vector \em allGrads_vec correspods to a component
      of the vector field.

      Each column of a matrix component in \em allGrads_vec corresponds to one
      evaluation point. In this column, the gradients of all active (i.e., non-zero)
      basis functions are stored "one above the other".\n


      Tha basis values information on the parameter domain at a certain point for a
      vector field is given in \em <b>allValues</b>:\n

      Each column of the matrix \em allValues corresponds to one vector component
      in the parametric domain. In this column, the value of all active (i.e., non-zero)
      basis functions are stored "one above the other".\n

      Example: Let \f$B_i(s,t), i = 1,...,9\f$ be a set of bivariate basis functions.
      Then, a column of a matrix component in \em allGrads_vec reads\n
      \f[
      ( \partial_s B_1, \partial_t B_1, \partial_s B_2, \partial_t B_2, \ldots, \partial_t B_9 )^T.
      \f]


      \param k Indicates which column of the matries in \em allGrads_vec should be transformed.
      \param[in] allGrads_vec std::vector of gsMatrix containing computed gradients in the format
      described above for each component of the vector field.
      \param[out] trfGradsK_vec std::vector of gsMatrix with the corresponding gradients on the
      physical domain in the format as described above for each component of the vector field.
    */
    virtual void transformGradsHdiv(index_t k, const gsMatrix<T>& allValues, // NumAct X TargetDim
                                    std::vector<gsMatrix<T> > const & allGrads_vec,
                                    std::vector<gsMatrix<T> > & trfGradsK_vec) const = 0;


    /// \brief Transforms parametric gradients to gradients on the physical domain.
    ///
    /// See the other transformGradients() for documentation.
    virtual void transformGradients(index_t k, const typename gsMatrix<T>::Block allGrads, gsMatrix<T>& trfGradsK) const = 0;

    /// Computes the outer normal on side \a s
    /// Assumes that points \a u are on the side 
    /// result is a vector containing the normal vector
    /// NB! Contrary to usual calculus, the normal vector does not have norm 1, but
    /// its norm represent the local surface area.
    virtual void outerNormal(index_t k, boundary::side s, gsVector<T> & result)  const = 0;

    /// Computes the normal vector of a co-dimension one geometry at evaluation point \a k
    virtual void normal(index_t k, gsVector<T> & result)  const = 0;

    /// Computes the divergence of the geometry at evaluation point \k
    virtual void divergence(gsVector<T> & result)  const = 0;

/*
    const          gsVector<T>              & divs        ()          const {return m_div;}
                   T                          div         (index_t k) const {return m_div(k);}

    const          gsMatrix<T>              & curls       ()          const {return m_curl;}
    const typename gsMatrix<T>::constColumn   curl        (index_t k) const {return m_curl.col(k);}

    const          gsMatrix<T>              & laplacians  ()          const {return m_lap;}
    const typename gsMatrix<T>::constColumn   laplacian   (index_t k) const {return m_lap.col(k);}

    const          gsMatrix<T>              & normals     ()          const {return m_normal;}
    const typename gsMatrix<T>::constColumn   normal      (index_t k) const {return m_normal.col(k);}
*/

protected:
    unsigned           m_flags;
    unsigned           m_parDim;
    unsigned           m_geoDim;// unused
    int                m_orientation;

    gsMatrix<T>        m_values;
    gsMatrix<T>        m_jacobians;
    gsMatrix<T>        m_jacInvs;
    gsVector<T>        m_measures;
    gsMatrix<T>        m_2ndDers;

/*
    gsVector<T>        m_div;
    gsMatrix<T>        m_curl;
    gsMatrix<T>        m_lap;
    gsMatrix<T>        m_normal;
*/
private:
    // disable copying
    gsGeometryEvaluator(const gsGeometryEvaluator& other);
    gsGeometryEvaluator& operator= (const gsGeometryEvaluator& other);
};


//////////////////////////////////////////////////
//////////////////////////////////////////////////

template <class T, int ParDim, int codim>
class gsGenericGeometryEvaluator : public gsGeometryEvaluator<T>
{
public:
    typedef T Scalar_t;
    static const int GeoDim = ParDim + codim;

public:
    gsGenericGeometryEvaluator(const gsGeometry<T> & geo, unsigned flags)
        : gsGeometryEvaluator<T>(flags,ParDim,GeoDim), m_geo(geo), m_maxDeriv(-1)
    {
        m_orientation = m_geo.orientation();
        setFlags(flags);
    }

    void setFlags (unsigned newFlags)
    {
        m_flags=newFlags;
        if (m_flags & (NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DIV 
                       | NEED_CURL | NEED_NORMAL | NEED_OUTER_NORMAL ))
            m_flags |= NEED_JACOBIAN;
        if (m_flags & NEED_LAPLACIAN)
            m_flags |= NEED_2ND_DER;
        if (m_flags & NEED_VALUE)
            m_maxDeriv=0;
        if (m_flags & NEED_JACOBIAN)
            m_maxDeriv=1;
        if (m_flags & NEED_2ND_DER)
            m_maxDeriv=2;
    }

    // Documentation at gsGeometryEvaluator::evaluateAt
    void evaluateAt(const gsMatrix<T>& u);

public:

    // Documentation at gsGeometryEvaluator::transformGradients
    void transformGradients(index_t k, const gsMatrix<T>& allGrads, gsMatrix<T>& trfGradsK) const
    {
        GISMO_ASSERT(this->m_flags & NEED_GRAD_TRANSFORM, "J^-1 not computed");
        GISMO_ASSERT(allGrads.rows() % ParDim == 0, "Invalid size of gradient matrix");

        const index_t numGrads = allGrads.rows() / ParDim;
        const gsAsConstMatrix<T,ParDim> grads_k(allGrads.col(k).data(), ParDim, numGrads);
        trfGradsK.noalias() = m_jacInvs.template block<GeoDim,ParDim>(0, k*ParDim) * grads_k;
    }

    // output: NumAct X TargetDim of basis values transform with the Piola transformation
    void transformValuesHdiv(index_t k, const gsMatrix<T>& allValues, // NumAct X TargetDim
                             gsMatrix<T>  & trfValues) const;

    // output: d x numAct matrix of divergence preserving transformed gradients at point k
    void transformGradsHdiv(index_t k, const gsMatrix<T>& allValues, // NumAct X TargetDim
                           std::vector<gsMatrix<T> > const & allGrads_vec,
                           std::vector<gsMatrix<T> > & trfGradsK_vec) const;

    // output: d x numAct matrix of transformed gradients at point k
    void transformGradients(index_t k, const typename gsMatrix<T>::Block allGrads, 
                            gsMatrix<T>& trfGradsK) const
    {
        GISMO_ASSERT(this->m_flags & NEED_GRAD_TRANSFORM, "J^-1 not computed");
        GISMO_ASSERT(allGrads.rows() % ParDim == 0, "Invalid size of gradient matrix");

        const index_t numGrads = allGrads.rows() / ParDim;
        const gsAsConstMatrix<T,ParDim> grads_k(allGrads.col(k).data(), ParDim, numGrads);
        trfGradsK.noalias() = m_jacInvs.template block<GeoDim,ParDim>(0, k*ParDim) * grads_k;

        /*
          const index_t numGrads = allGrads.rows() / ParDim;
          trfGradsK.resize( GeoDim, numGrads );
          for (index_t i = 0; i < numGrads; ++i)
          {
          trfGradsK.template block<GeoDim, 1>(0, i) =
          m_jacInvs.template block<GeoDim,ParDim>(0, k*ParDim).transpose() *
          allGrads.template block<ParDim, 1>(i * ParDim, k);
          }
          // */
    }

    /// \note This implementation silently assumes a tensor domain
    void outerNormal(index_t k, boundary::side s, gsVector<T> & result) const
    {
        GISMO_ASSERT(this->m_flags & NEED_JACOBIAN, "Jacobians not computed");

        const T   sgn = sideOrientation(s) * m_orientation;
        const int dir = direction(s);

        // assumes points u on boundary "s"
        //GISMO_ASSERT( ParDim == GeoDim, "Codim different than one");
        result.resize(GeoDim);
        const gsMatrix<T, ParDim, ParDim> Jk = 
            m_jacobians.template block<ParDim,ParDim>(0, k*ParDim);

        T alt_sgn = sgn;
        gsMatrix<T, ParDim-1, ParDim-1> minor;        
        for (int i = 0; i != ParDim; ++i) // for all components of the normal
        {
            Jk.firstMinor(i, dir, minor);
            result[i] = alt_sgn * minor.determinant();
            alt_sgn = -alt_sgn;
        }
    }

    /// \note This implementation silently assumes a tensor domain
    void normal(index_t k, gsVector<T> & result) const
    {
        GISMO_ASSERT(this->m_flags & NEED_JACOBIAN, "Jacobians not computed");
        GISMO_ASSERT( ParDim+1 == GeoDim, "Codimension should be equal to one");
        result.resize(ParDim+1);
        const gsMatrix<T, ParDim+1, ParDim> Jk = 
            m_jacobians.template block<ParDim+1,ParDim>(0, k*ParDim);

        T alt_sgn(1.0);
        gsMatrix<T, ParDim, ParDim> minor;
        for (int i = 0; i <= ParDim; ++i) // for all components of the normal vector
        {
            Jk.rowMinor(i, minor);
            result[i] = alt_sgn * minor.determinant();
            alt_sgn = -alt_sgn;
        }
    }

    /// \note This implementation silently assumes a tensor domain
    void divergence(gsVector<T> & result) const
    {
        GISMO_ASSERT(this->m_flags & NEED_JACOBIAN, "Jacobians not computed");
        GISMO_ASSERT( ParDim == GeoDim, "Codimension should be equal to zero");

        result.resize(m_numPts);

        for (index_t j=0; j < m_numPts; ++j)  // for all evaluation points
        {
            result[j] = ( m_jacobians.template block<ParDim,ParDim>(0, j*ParDim) ).trace();
        }
    }

    const gsGeometry<T> & geometry() const {return m_geo;}

private:
    void computeValues();
    void computeJacobians();
    //JS2 Not tested if it gives the correct values
    void compute2ndDerivs();

    void computeCurl();
    void computeLaplacian();
    void computeNormal();

    // disable copying
    gsGenericGeometryEvaluator(const gsGenericGeometryEvaluator& other);
    gsGenericGeometryEvaluator& operator=(const gsGenericGeometryEvaluator& other);

protected:

    using gsGeometryEvaluator<T>::m_orientation;
    using gsGeometryEvaluator<T>::m_values;
    using gsGeometryEvaluator<T>::m_jacobians;
    using gsGeometryEvaluator<T>::m_measures;
    using gsGeometryEvaluator<T>::m_jacInvs;
    using gsGeometryEvaluator<T>::m_2ndDers;
    using gsGeometryEvaluator<T>::m_flags;

/*
    using gsGeometryEvaluator<T>::m_lap;
    using gsGeometryEvaluator<T>::m_div;
    using gsGeometryEvaluator<T>::m_curl;
    using gsGeometryEvaluator<T>::m_normal;
*/

    gsMatrix<T>           m_basisVals;
    gsMatrix<unsigned>    m_active;
private:
    const gsGeometry<T> & m_geo;
    index_t m_numPts;
    int m_maxDeriv;
};


} // namespace gismo

#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsGeometryEvaluator.hpp)
#endif
