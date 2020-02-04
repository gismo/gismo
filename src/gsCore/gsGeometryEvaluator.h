/** @file gsGeometryEvaluator.h

    @brief Provides declaration of GeometryEvaluator interface.
    
    This file is part of the G+Smo library-
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, C. Hofreither, A. Mantzaflaris, J. Sogn
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsBoundary.h>


namespace gismo
{

template <class T=real_t>                class gsGeometryEvaluator;

/**
    \brief interface to compute a geometry on sets of points

    the idea of geometry evaluators is that they act as a local
    cache for computed values. In this way it is possible to pass
    a unique argument to functions containing all the required values.

    All accessors to the stored values come in two versions: plural and
    singular. The plural returns a matrix with the computed values on all
    points, the singular takes as additional parameter that represent the
    point index and returns the associated sub-block.

    \tparam T the coefficient type
    
    \ingroup Core
**/
template <typename T>
class gsGeometryEvaluator
{
public:
    /// Shared pointer for gsDomainIterator
    typedef memory::shared_ptr< gsGeometryEvaluator > Ptr;
    /// Unique pointer for gsDomainIterator
    typedef memory::unique_ptr< gsGeometryEvaluator > uPtr;
    
public:
    /**
       \brief Constructor using the geometry and the flags that define
       what should be evaluated.
       
       See gismo::gsNeedEnum for available flags.
    */
    gsGeometryEvaluator(const gsGeometry<T> & geo, unsigned flags)
    : m_geo(geo), m_flags(flags), m_parDim(geo.parDim())
    {}
    
    virtual ~gsGeometryEvaluator() {}

public:
    
    /**
       \brief change the values to compute

        This interface is provided in such a way that it is possible
        split the logic of the initialization of the evaluator.
        (?)
    **/
    void addFlags (unsigned newFlags) {this->setFlags(m_flags|newFlags);}

    /**
        \brief set the values to compute

        This interface is provided in such a way that it is possible
        split the logic of the initialization of the evaluator.
    **/
    virtual void setFlags (unsigned newFlags) {m_flags = newFlags;}

    /**
        \brief returns which values are computed
    **/
    unsigned getFlags () const {return m_flags;}

    /**
        \brief Returns the number of evaluation points that the
        evaluator currently holds.
    **/
    index_t numPoints () const {return m_numPts;}

    /**
        \brief Returns the number of evaluation points that the
        evaluator currently holds.
    **/
    index_t parDim () const {return m_geo.parDim();}

    const gsGeometry<T> & geometry() const {return m_geo;}

    inline size_t id() const {return m_geo.id();}

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
    
      \returns J gsMatrix of size <em>physDim</em> x <em>parDim*N</em>, where\n
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

        if the dimension of the parametric domain equals that of the target
        space this is the determinant of the jacobian \f$\det J\f$.
        If the geometry is a curve it is the modulus of the normal.
        If it is a surface in R^3 it is the length of the normal.
        In general it can always be computed as \f$\det(J^t J)^{1/2}\f$
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
        the parameter domain.

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

    /** \brief Returns the second derivatives.
       
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
    * \return gsMatrix of size <em>d</em> x <em>d*N</em>, where\n
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
      Example: Let \f$x\f$ and \f$y\f$ be the coordinates in the physical domain. Then,
      in the above example, \em <b>trfGradsK</b> is given by\n
      \f[
      \left(\begin{array}{cccc}
      \partial_x B_1 & \partial_x B_2 & \ldots & \partial_x B_9 \\
      \partial_y B_1 & \partial_y B_2 & \ldots & \partial_y B_9
      \end{array}\right)
      \f]
      \param[in] k Indicates which column of \em allGrads should be transformed.
      \param[in] allGrads gsMatrix containing computed gradients in the format
      described above.
      \param[out] trfGradsK gsMatrix with the corresponding gradients on the
      physical domain in the format as described above.
    */
    // rename to transformGradsHgrad
    virtual void transformGradients(index_t k, const gsMatrix<T>& allGrads, gsMatrix<T>& trfGradsK)  const = 0;

    /**
      \brief Transforms parametric basis values at point \a k of a vector field, while preserving the divergence (Piola Transformation).
      NOT TESTED
    */
    // allValues: NumAct X TargetDim
    virtual void transformValuesHdiv(index_t k, 
                                     const std::vector<gsMatrix<T> >& allValues, 
                                     gsMatrix<T>  & result) const = 0;

// S.K.: Commented declaration, because the function code in the hpp-file
// is completely commented out an the function is not used yet.
//
//    /*
//      \brief Transforms parametric gradients to a gradients on the physical domain of a vector field while preserving the divergence.

//      The transformation is sometimes called the contravariant Piola transform:
//      \f[
//        \mathcal{F}^{\text{div}}(\hat{\mathbf{u}}) =\frac{1}{\text{det}\, J} J\,\hat{\mathbf{u}} \circ F^{-1},
//      \f]
//      where \f$F\f$ is mapping from \f$\hat{\Omega}\f$ to \f$\Omega\f$ and \f$J\f$
//      is the jacobian.

//      The gradient information on the parameter domain at a certain points for a
//      vector field is given in \em <b>allGrads</b> in the following format:\n

//      Each component of the vector \em allGrads corresponds to a component
//      of the vector field.

//      Each column of a matrix component in \em allGrads corresponds to one
//      evaluation point. In this column, the gradients of all active (i.e., non-zero)
//      basis functions are stored "one above the other".\n

//      Example: Let \f$B_i(s,t), i = 1,...,9\f$ be a set of bivariate basis functions.
//      Then, a column of a matrix component in \em allGrads reads\n
//      \f[
//      ( \partial_s B_1, \partial_t B_1, \partial_s B_2, \partial_t B_2, \ldots, \partial_t B_9 )^T.
//      \f]

//      Tha basis values information on the parameter domain at a certain point for a
//      vector field is given in \em <b>allValues</b>:\n

//      Each column of the matrix \em allValues corresponds to one vector component
//      in the parametric domain. In this column, the value of all active (i.e., non-zero)
//      basis functions are stored "one above the other".\n


//      \param[in] k Indicates which column of the matrices in \em allegros_vec should be transformed.
//      \param[in] allValues std::vector of gsMatrix containing function values in the format
//      described above for each component of the vector field.
//      \param[in] allGrads std::vector of gsMatrix containing computed gradients in the format
//      described above for each component of the vector field.
//      \param[out] result std::vector of gsMatrix with the corresponding gradients on the
//      physical domain in the format as described above for each component of the vector field.
//    */
//    virtual void transformGradsHdiv(index_t k,
//                                    const std::vector<gsMatrix<T> >& allValues,
//                                    std::vector<gsMatrix<T> > const & allGrads,
//                                    gsMatrix<T> & result) const = 0;

    /// \brief Transforms parametric gradients to gradients on the physical domain.
    ///
    /// See the other transformGradients() for documentation.
    virtual void transformGradients(index_t k, const typename gsMatrix<T>::Block allGrads, gsMatrix<T>& trfGradsK) const = 0;

    /**
    \brief Transforms paramatric 1st and 2nd derivatives to Laplacians on the physical domain.

    The gradient information on the parameter domain at a certain points
    is given in \em <b>allGrads</b> in the following
    format:\n
    Each column of \em allGrads corresponds to one evaluation point. In this column, the
    gradients of all active (i.e., non-zero) basis functions are stored "one above the other".\n
    Example: Let \f$B_i(r,s,t), i = 1,...,9\f$ be a set of bivariate basis functions.
    Then, a column of \em allGrads reads\n
    \f[
    ( \partial_r B_1, \partial_s B_1, \partial_t B_1, \partial_r B_2, \partial_s B_2, \partial_t B_2, \ldots, \partial_t B_9 )^T.
    \f]
    Hence, \em <b>allGrads</b> is a gsMatrix of size <em>(n * ParDim)</em> x <em>K</em>,
    where \n
    ...\em n denotes the number of active functions, \n
    ...<em>ParDim</em> denotes the dimension of the parameter domain, and\n
    ...\em K denotes the number of columns of \em allGrads.\n
    \n
    The second derivative (the Hessian) on the parameter domain at a certain points
    is given in \em <b>allHessians</b> in the following
    format:\n
    Each column of \em allHessians corresponds to one evaluation point. In this column, the
    second derivatives of all active (i.e., non-zero) basis functions are stored "one above the other".\n
    Example: Let \f$B_i(r,s,t), i = 1,...,9\f$ be a set of bivariate basis functions.
    Then, a column of \em allHessians reads\n
    \f[
    ( \partial_{rr} B_1, \partial_{ss} B_1, \partial_{tt} B_1, \partial_{rs} B_1, \partial_{rt} B_1, \partial_{st} B_1, \partial_{rr} B_2, \partial_{ss} B_2, \partial_{tt} B_2, \partial_{rs} B_2, \partial_{rt} B_2, \partial_{st} B_2, \ldots, \partial_{st} B_9 )^T.
    \f]
    Hence, \em <b>allHessians</b> is a gsMatrix of size <em>(n * (ParDim + (ParDim*(ParDim-1))/2))</em> x <em>K</em>,
    where \n
    ...\em n denotes the number of active functions, \n
    ...<em>ParDim</em> denotes the dimension of the parameter domain, and\n
    ...\em K denotes the number of columns of \em allHessians.\n
    \n


    These gradients and Hessians in the <em>k</em>-th column of \em allGrads and \em allHessians
    are transformed to the physical domain and stored in \em <b>result</b>.
    \em result is of size <em>1</em> x <em>n</em>. The <em>i</em>-th column
    of \em result corresponds to the transformed Laplacian of the <em>i</em>-th function.\n
    \n
    Example: Let \f$x\f$, \f$y\f$ and \f$z\f$ be the coordinates in the physical domain. Then,
    in the above example, \em <b>result</b> is given by\n
    \f[
    \left(\begin{array}{cccc}
    \partial_{xx} B_1 + \partial_{yy} B_1 + \partial_{zz} B_1  & \partial_{xx} B_2 + \partial_{yy} B_2 + \partial_{zz} B_2  & \ldots & \partial_{xx} B_9 + \partial_{yy} B_9 + \partial_{zz} B_9
    \end{array}\right)
    \f]


    \param k Indicates which column of \em allGrads should be transformed.
    \param[in] allGrads gsMatrix containing computed gradients in the format
    described above.
    \param[in] allHessians gsMatrix containing computed second derivatives in the format
    described above.
    \param[out] result gsMatrix with the corresponding gradients on the
    physical domain in the format as described above.
    */
    virtual void transformLaplaceHgrad(index_t k,
                                       const gsMatrix<T> & allGrads,
                                       const gsMatrix<T> & allHessians,
                                       gsMatrix<T> & result) const = 0;

    /**
    \brief Transforms paramatric 1st and 2ed derivatives to 2nd derivatives on the physical domain.

    The gradient information on the parameter domain at a certain points
    is given in \em <b>allGrads</b> in the following
    format:\n
    Each column of \em allGrads corresponds to one evaluation point. In this column, the
    gradients of all active (i.e., non-zero) basis functions are stored "one above the other".\n
    Example: Let \f$B_i(r,s,t), i = 1,...,9\f$ be a set of bivariate basis functions.
    Then, a column of \em allGrads reads\n
    \f[
    ( \partial_r B_1, \partial_s B_1, \partial_t B_1, \partial_r B_2, \partial_s B_2, \partial_t B_2, \ldots, \partial_t B_9 )^T.
    \f]
    Hence, \em <b>allGrads</b> is a gsMatrix of size <em>(n * ParDim)</em> x <em>K</em>,
    where \n
    ...\em n denotes the number of active functions, \n
    ...<em>ParDim</em> denotes the dimension of the parameter domain, and\n
    ...\em K denotes the number of columns of \em allGrads.\n
    \n
    The second derivative (the Hessian) on the parameter domain at a certain points
    is given in \em <b>allHessians</b> in the following
    format:\n
    Each column of \em allHessians corresponds to one evaluation point. In this column, the
    second derivatives of all active (i.e., non-zero) basis functions are stored "one above the other".\n
    Example: Let \f$B_i(r,s,t), i = 1,...,9\f$ be a set of bivariate basis functions.
    Then, a column of \em allHessians reads\n
    \f[
    ( \partial_{rr} B_1, \partial_{ss} B_1, \partial_{tt} B_1, \partial_{rs} B_1, \partial_{rt} B_1, \partial_{st} B_1, \partial_{rr} B_2, \partial_{ss} B_2, \partial_{tt} B_2, \partial_{rs} B_2, \partial_{rt} B_2, \partial_{st} B_2, \ldots, \partial_{st} B_9 )^T.
    \f]
    Hence, \em <b>allHessians</b> is a gsMatrix of size <em>(n * (ParDim + (ParDim*(ParDim-1))/2))</em> x <em>K</em>,
    where \n
    ...\em n denotes the number of active functions, \n
    ...<em>ParDim</em> denotes the dimension of the parameter domain, and\n
    ...\em K denotes the number of columns of \em allHessians.\n
    \n


    These gradients and Hessians in the <em>k</em>-th column of \em allGrads and \em allHessians
    are transformed to the physical domain and stored in \em <b>result</b>.
    \em result is of size <em>number of Gradients</em> x <em>GeoDim*(GeoDim+1)/2</em>. Each row
    corresponds to one active function and each row gives the values of the second derivatives.\n
    \n

    \param k Indicates which column of \em allGrads should be transformed.
    \param[in] allGrads gsMatrix containing computed gradients in the format
    described above.
    \param[in] allHessians gsMatrix containing computed second derivatives in the format
    described above.
    \param[out] result gsMatrix with the corresponding second derivatives on the
    physical domain in the format as described above.
    */
    virtual void transformDeriv2Hgrad(  index_t k,
                            const gsMatrix<T> & allGrads,
                            const gsMatrix<T> & allHessians,
                            gsMatrix<T> & result) const = 0;

    /// Computes the outer normal on side \a s
    /// Assumes that points \a u are on the side 
    /// result is a vector containing the normal vector
    /// NB! Contrary to usual calculus, the normal vector does not have norm 1, but
    /// its norm represent the local surface area.
    virtual void outerNormal(index_t k, boxSide s, gsVector<T> & result)  const = 0;

    /// Computes the normal vector of a co-dimension one geometry at evaluation point \a k
    virtual void normal(index_t k, gsVector<T> & result)  const = 0;

    /// Computes the divergence of the geometry at evaluation point \a k
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

    const gsGeometry<T> & m_geo;
    unsigned           m_flags;

    index_t            m_numPts;
    unsigned           m_parDim;
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


template <class T, short_t ParDim, short_t codim>
class gsGenericGeometryEvaluator : public gsGeometryEvaluator<T>
{
public:
    typedef gsGeometryEvaluator<T> Base;
    typedef T Scalar_t;
    static const short_t GeoDim = ParDim + codim;

public:
    gsGenericGeometryEvaluator(const gsGeometry<T> & geo, unsigned flags)
    : gsGeometryEvaluator<T>(geo,flags), m_maxDeriv(-1)
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
    void transformValuesHdiv(index_t k, 
                             const std::vector<gsMatrix<T> >& allValues, // NumAct X TargetDim
                             gsMatrix<T>  & result) const;

    // output: d x numAct matrix of divergence preserving transformed gradients at point k
    void transformGradsHdiv(index_t k, 
                            const std::vector<gsMatrix<T> >& allValues,
                            std::vector<gsMatrix<T> > const & allGrads,
                            gsMatrix<T> & result) const;

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

    void transformLaplaceHgrad(index_t k,
                               const gsMatrix<T> & allgrads,
                               const gsMatrix<T> & allHessians,
                               gsMatrix<T> & result) const;

    void transformDeriv2Hgrad(  index_t k,
                            const gsMatrix<T> & allGrads,
                            const gsMatrix<T> & allHessians,
                            gsMatrix<T> & result) const;


    /// \note This implementation silently assumes a tensor domain
    void outerNormal(index_t k, boxSide s, gsVector<T> & result) const;

    /// \note This implementation silently assumes a tensor domain
    void normal(index_t k, gsVector<T> & result) const
    {
        GISMO_ASSERT(this->m_flags & NEED_JACOBIAN, "Jacobians not computed");
        GISMO_ASSERT( ParDim+1 == GeoDim, "Codimension should be equal to one");
        result.resize(ParDim+1);
        const gsMatrix<T, ParDim+1, ParDim> Jk = 
            m_jacobians.template block<ParDim+1,ParDim>(0, k*ParDim);

        T alt_sgn(1.0);
        typename gsMatrix<T,ParDim+1,ParDim>::RowMinorMatrixType minor;
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

    std::vector<gsMatrix<T> > m_basisVals;
    gsMatrix<unsigned>    m_active;
private:
    int m_maxDeriv;

    using Base::m_geo;
    using Base::m_numPts;
};


template <class T>
gsGeometryEvaluator<T> * evaluator1(unsigned flags, const gsGeometry<T> & g);
template <class T>
gsGeometryEvaluator<T> * evaluator2(unsigned flags, const gsGeometry<T> & g);
template <class T>
gsGeometryEvaluator<T> * evaluator3(unsigned flags, const gsGeometry<T> & g);
template <class T>
gsGeometryEvaluator<T> * evaluator4(unsigned flags, const gsGeometry<T> & g);

template <class T>
inline gsGeometryEvaluator<T> * getEvaluator(unsigned flags, const gsGeometry<T> & g)
{
    switch (g.parDim())
    {
        case 1:
            return evaluator1(flags, g);
        case 2:
            return evaluator2(flags, g);
        case 3:
            return evaluator3(flags, g);
        case 4:
            return evaluator4(flags, g);
        default:
            GISMO_ERROR("ParDim problem.");
    }
}

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsGeometryEvaluator.hpp)
#endif
