/** @file gsGeometry.h

    @brief Provides declaration of Geometry abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFunction.h>
#include <gsCore/gsBoundary.h>


#define GISMO_BASIS_ACCESSORS \
    Basis & basis() { return static_cast<Basis&>(*this->m_basis); } \
    const Basis & basis() const { return static_cast<const Basis&>(*this->m_basis); }
    // bool isProjective() const{ return Basis::IsRational; }

namespace gismo
{

/** 
    \brief
    Abstract base class representing a geometry map.

    The combination of a basis with a coefficient matrix of the proper
    size describes a function. Such objects are called geometries.


    Note that geometries are, essentially, functions mapping from the
    parameter domain to the physical space,

    \f[ G \ :\ \mathbb R^d \to \mathbb R^n \f]

    so that 

    \f[ G(\hat x) = x \qquad \hat x \in \mathbb R^d, \ x \in \mathbb R^n \f]

    and therefore they derive from gsFunction. A gsGeometry can thus be
    used wherever a gsFunction is expected.
    
    All geometry classes derive from the abstract base class
    gsGeometry, see the \ref geometry.  In general, there is a
    statically known one-to-one mapping between basis types and their
    associated geometry types. For instance, when you ask a
    gsBSplineBasis to create a geometry using some given coefficient
    matrix, you will get back an instance of gsBSpline. In turn, the
    basis that the gsBSpline object stores internally is of type
    gsBSplineBasis.

    The parameter domain of the basis is also the parameter domain of
    the geometry map, and has a certain source dimension \em d (see
    gsGeometry::parDim()). Depending on the parameter domain
    dimension, derived geometries inherit gsCurve, gsSurface, gsVolume
    or gsBulk, for dimension 1,2,3 or 4, respectively.

    The dimension \em n of the resulting geometry,
    i.e., of the image of the geometry map, is defined by the
    size of the given coefficients, and may be larger than \em d
    (see gsGeometry::geoDim()). $\,$

    For instance, for a B-spline basis (<em>d = 1</em>), we could
    have coefficients of dimension <em>n = 1</em>, which results
    in a parametrization of a 1-D interval, or we could have
    coefficients of size <em>n = 3</em>, which results in a B-spline
    curve in 3-D space.

    The coefficients are stored as an <em>s x n</em> matrix,
    where \em s is the size of the basis, i.e., the number of
    its basis functions.
    Every row of the coefficient matrix is an *n*-dimensional control
    point, i.e., the vector-valued coefficient for the associated
    basis function.

    Evaluation at parameter values is done with
    \ref func_eval_members "the Evaluation member functions" derived 
    from gsFunction.

    \tparam T coefficient type

    \ingroup function
    \ingroup geometry
    \ingroup Core
*/
template<class T>
class gsGeometry : public gsFunction<T>
{

public:
    /// Shared pointer for gsGeometry
    typedef memory::shared_ptr< gsGeometry > Ptr;

    /// Unique pointer for gsGeometry
    typedef memory::unique_ptr< gsGeometry > uPtr;

    typedef T Scalar_t;

public:

    /// @name Constructors
    /// @{
    
    /// @brief Default constructor.  Note: Derived constructors (except for
    /// the default) should assign \a m_basis to a valid pointer
    gsGeometry() :m_basis( NULL ), m_id(0)
    { }

    /// @brief Constructor by a basis and coefficient vector
    ///
    /// Coefficients are given by \em{give(coefs) and they are
    /// consumed, i.e. the \coefs variable will be empty after the call
    gsGeometry( const gsBasis<T> & basis, gsMatrix<Scalar_t> coefs) :
    m_basis(basis.clone().release()), m_id(0)
    {
        m_coefs.swap(coefs);
        GISMO_ASSERT( basis.size() == m_coefs.rows(), 
                      "The coefficient matrix of the geometry (rows="<<m_coefs.rows()
                      <<") does not match the number of basis functions in its basis("
                      << basis.size() <<").");
    }

    /// @brief Copy Constructor
    gsGeometry(const gsGeometry & o) 
    : m_coefs(o.m_coefs), m_basis(o.m_basis != NULL ? o.basis().clone().release() : NULL), m_id(o.m_id)
    { }

    /// @}

    gsGeometry& operator=( const gsGeometry & o)
    {
        if ( this != &o )
        {
            m_coefs = o.m_coefs;
            delete m_basis;
            m_basis = o.basis().clone().release() ;
            m_id = o.m_id;
        }
        return *this;
    }
    
    virtual ~gsGeometry() 
    {
        delete m_basis;
    }

#if EIGEN_HAS_RVALUE_REFERENCES
    gsGeometry(gsGeometry&& other) 
    : m_coefs(std::move(other.m_coefs)), m_basis(other.m_basis), 
      m_id(std::move(other.m_id))
    {
        other.m_basis = NULL;
    }
    gsGeometry & operator=(gsGeometry&& other)
    {
        m_coefs.swap(other.m_coefs); other.m_coefs.clear();
        delete m_basis;
        m_basis = other.m_basis; other.m_basis = NULL;
        m_id = std::move(other.m_id);
        return *this;
    }
#endif

public:

    /**
        @name Evaluation functions

        re-implemented from gsFunction, see also \ref func_eval_members
        @{
    */

    /** \brief Evaluate the function at points \a u into \a result.
     *
     * Let \em d be the dimension of the source space ( d = domainDim() ).\n
     * Let \em n be the dimension of the image/target space ( n = targetDim() ).\n
     * Let \em N denote the number of evaluation points.
     *
     * \param[in] u gsMatrix of size <em>d</em> x <em>N</em>, where each
     * column of \em u represents one evaluation point.
     * \param[out] result gsMatrix of size <em>n</em> x <em>N</em>, where each
     * column of \em u represents the result of the function at the
     * respective valuation point.
     */
    // Look at gsFunction class for documentation
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    { this->basis().evalFunc_into(u, m_coefs, result); }

    /** \brief Evaluate derivatives of the function
     * \f$f:\mathbb{R}^d\rightarrow\mathbb{R}^n\f$
     * at points \a u into \a result.
     *
     * Let \em d be the dimension of the source space ( d = domainDim() ).\n
     * Let \em n be the dimension of the image/target space ( n = targetDim() ).\n
     * Let \em N denote the number of evaluation points.
     *
     * Let \f$ f:\mathbb R^2 \rightarrow \mathbb R^3 \f$, i.e.,
     * \f$ f(x,y) = ( f^{(1)}(x,y), f^{(2)}(x,y), f^{(3)}(x,y) )^T\f$,\n
     * and let
     * \f$ u = ( u_1, \ldots, u_N) = ( (x_1,y_1)^T, \ldots, (x_N, y_N)^T )\f$.\n
     * Then, \em result is of the form
     * \f[
     \left[
     \begin{array}{cccc}
        \partial_x f^{(1)}(u_1) & \partial_x f^{(1)}(u_2)
           & \ldots & \partial_x f^{(1)}(u_N) \\
        \partial_y f^{(1)}(u_1) & \partial_y f^{(1)}(u_2)
           & \ldots & \partial_y f^{(1)}(u_N) \\
        \partial_x f^{(2)}(u_1) & \partial_x f^{(2)}(u_2)
           & \ldots & \partial_x f^{(2)}(u_N) \\
        \partial_y f^{(2)}(u_1) & \partial_y f^{(2)}(u_2)
           & \ldots & \partial_x f^{(2)}(u_N) \\
        \partial_x f^{(3)}(u_1) & \partial_x f^{(3)}(u_2)
           & \ldots & \partial_x f^{(3)}(u_N)\\
        \partial_y f^{(3)}(u_1) & \partial_y f^{(3)}(u_2)
           & \ldots & \partial_y f^{(3)}(u_N)
     \end{array}
     \right]
     \f]
     *
     * \param[in] u gsMatrix of size <em>d</em> x <em>N</em>, where each
     * column of \em u represents one evaluation point.
     * \param[out] result gsMatrix of size <em>(d * n)</em> x <em>N</em>.
     * Each row of \em result corresponds to one component in the target
     * space and contains the gradients for each evaluation point,
     * as row vectors, one after the other (see above for details on the format).
     *
     * \warning By default, gsFunction uses central finite differences
     * with h=0.00001! One must override this function in derived
     * classes to get proper results.
     */
    // Look at gsFunction class for documentation
    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    { this->basis().derivFunc_into(u, m_coefs, result); }


    /** @brief Evaluate second derivatives of the function at points \a u into \a result.
     *
     * Let \em d be the dimension of the source space ( d = domainDim() ).\n
     * Let \em n be the dimension of the image/target space ( n = targetDim() ).\n
     * Let \em N denote the number of evaluation points.
     *
     * \param[in] u gsMatrix of size <em>d</em> x <em>N</em>, where each
     * column of \em u represents one evaluation point.
     * \param[out] result gsMatrix of size <em>(S*n)</em> x <em>N</em>,
     * where <em>S=d*(d+1)/2</em>.\n
     * Each column in \em result corresponds to one point (i.e., one column in \em u)\n
     * and contains the following values (for <em>d=3</em>, <em>n=3</em>):\n
     * \f$ (\partial_{xx} f^{(1)}, \partial_{yy} f^{(1)}, \partial_{zz} f^{(1)}, \partial_{xy} f^{(1)},
       \partial_{xz} f^{(1)}, \partial_{yz} f^{(1)}, \partial_{xx} f^{(2)},\ldots,\partial_{yz} f^{(3)} )^T\f$
     * \warning By default uses central finite differences with h=0.00001!
     * One must override this function in derived
     * classes to get proper results.
     */
    // Look at gsFunctionSet class for documentation
    void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    { this->basis().deriv2Func_into(u, m_coefs, result); }

    void evalAllDers_into(const gsMatrix<T> & u, int n,
                          std::vector<gsMatrix<T> > & result) const
    { this->basis().evalAllDersFunc_into(u, m_coefs, n, result); }

    /// @}


    // Look at gsFunctionSet for documentation
    virtual void compute(const gsMatrix<T> & in, gsFuncData<T> & out) const;

    /// \brief Evaluates if the geometry orientation coincide with the
    /// ambient orientation.
    /// This is computed in the center of the parametrization and will
    /// fail to be meaningful if the geometry is singular.
    /// returns one if codimension is not zero
    int orientation() const
    {
        if ( parDim() == geoDim() )
        {
            const T val = gsFunction<T>::jacobian( parameterCenter() ).determinant();
            return (T(0) < val) - (val < T(0));
        }
        return 1;
    }

    /// @}

    /*************************************************************************/

    /// @name Accessors
    /// @{

    /// \brief Returns a const reference to the basis of the geometry.
    /// \note This function will return the derived concrete type of the basis.
    virtual const gsBasis<T> & basis() const = 0;

    /// \brief Returns a reference to the basis of the geometry.
    /// \note This function will return the derived concrete type of the basis.
    virtual       gsBasis<T> & basis()       = 0;

    /// Dimension of the ambient physical space (overriding gsFunction::targetDim())
    short_t targetDim() const { return this->coefDim(); }

    /// Dimension \em n of the coefficients (control points)
    short_t coefDim() const { return static_cast<short_t>(m_coefs.cols()); }

    /// Dimension \em n of the absent physical space
    short_t geoDim() const { return this->coefDim(); }

    /// Dimension \em d of the parameter domain (overriding gsFunction::domainDim()).
    virtual short_t domainDim() const { return this->basis().domainDim(); }

    /// Dimension \em d of the parameter domain (same as domainDim()).
    short_t parDim() const { return this->basis().domainDim(); }

    /// Co-dimension of the geometric object
    short_t coDim() const { return coefDim()-this->basis().domainDim(); }

    /// Returns the range of parameters (same as parameterRange())
    gsMatrix<T> support() const
    { return this->basis().support(); }

    /// Returns the range of parameters as a matrix with two columns, [lower upper]
    gsMatrix<T> parameterRange() const
    { return this->basis().support(); }

    /// Returns a "central" point inside inside the parameter domain
    virtual gsMatrix<T> parameterCenter() const
    { 
        // default impl. assumes convex support
        gsMatrix<T> S = this->basis().support();
        return ( S.col(0) + S.col(1) ) * T(0.5);
    }

    /// Get coordinates of the boxCorner \a bc in the parameter domain
    gsMatrix<T> parameterCenter( const boxCorner& bc );

    /// Get coordinates of the midpoint of the boxSide \a bs in the parameter domain
    gsMatrix<T> parameterCenter( const boxSide& bs );

    /// Get back the side of point \a u
    //boxSide sideOf(const gsVector<T> & u); //

    /// Check if points \a u also lie on the geometry and if required computes the if the points in \a u lie on one of the boundaries of the geometry
    /**
     *
     * @param u Matrix of points of the form geoDim() x #points
     * @param onG2 gsVector of booleans which indicate if the point is in the domain or outside
     * @param preIm Matrix of preimages of the points u
     * @param lookForBoundary if required the boundaries are computed the points in u lie on
     * @return A std::vector of boxSides containing the numbers of the sides or zero if the points are in the interior
     * If the flag lookForBoundary is not set then a vector containing anything will be returned
     */
    std::vector<boxSide> locateOn(const gsMatrix<T> & u, gsVector<bool> & onG2, gsMatrix<T> & preIm, bool lookForBoundary = false, real_t tol = 1.e-6) const; //

    // Whether the coefficients of this geometry are stored in projective or affine form
    //virtual bool isProjective() const = 0;

    /// @}

    /*************************************************************************/

    /** @name Coefficient access functions

        These functions allow direct access to the coefficient matrix of the geometry.
        @{
    */

    /// Returns the coefficient matrix of  the geometry
    /// Coefficient matrix of size coefsSize() x geoDim()
    // todo: coefsSize() x (geoDim() + 1) if projective
          gsMatrix<T> & coefs()       { return this->m_coefs; }

    /// Returns the coefficient matrix of  the geometry
    const gsMatrix<T> & coefs() const { return this->m_coefs; }

    /// Returns the i-th coefficient of the geometry as a row expression
    typename gsMatrix<T>::RowXpr       coef(index_t i)       { return m_coefs.row(i); }

    /// Returns the i-th coefficient of the geometry as a row expression
    typename gsMatrix<T>::ConstRowXpr  coef(index_t i) const { return m_coefs.row(i); }

    /// Returns the j-th coordinate of the i-th coefficient of the geometry
    T & coef(index_t i, index_t j)
    {
        GISMO_ASSERT( ((i < m_coefs.rows()) && (j < m_coefs.cols()) ),
                      "Coefficient or coordinate which is out of range requested.");
        return m_coefs(i,j);
    }

    /// Returns the j-th coordinate of the i-th coefficient of the geometry
    const T coef(index_t i, index_t j) const
    {
        GISMO_ASSERT( ((i < m_coefs.rows()) && (j < m_coefs.cols()) ),
                      "Coefficient or coordinate which is out of range requested.");
        return m_coefs(i,j);
    }

    /// Set the coefficient matrix of the geometry, taking ownership of the matrix
    void setCoefs(gsMatrix<T> cc) { this->m_coefs.swap(cc); }

    /// Return the number of coefficients (control points)
    unsigned coefsSize() const { return m_coefs.rows(); }
    // Warning: This can cause some clash while using periodic basis, since the ghost coefs are stored but ignored.



    /// @}

    /*************************************************************************/

    /** @name Transformation functions

        These functions apply various linear and affine transformations
        to the coefficients.
        @{
    */

    /// Apply the given square matrix to every control point.
    void linearTransform(const gsMatrix<T>& mat)
    {
        this->m_coefs = this->m_coefs * mat.transpose();
    }

    /// Apply 3D Rotation by \a angle radians around axis \a axis
    void rotate(T angle, const gsVector<T,3> & axis )
    {
        assert( geoDim() == 3 );
        Eigen::Transform<T,3,Eigen::Affine> 
            rot( Eigen::AngleAxis<T> (angle,axis.normalized()) );
        // To do: Simpler way to use transforms ?
        this->m_coefs = (this->m_coefs.rowwise().homogeneous() * 
                         rot.matrix().transpose() ).leftCols(3) ;
    }

    /// Apply 2D Rotation by \a angle radians
    void rotate(T angle)
    {
        GISMO_ASSERT( geoDim() == 2, "Only for 2D");
        Eigen::Rotation2D<T> rot(angle);
        this->m_coefs *= rot.matrix().transpose();
    }

    /// Apply Scaling by factor \a s
    void scale(T s, int coord = -1)
    {
        if ( coord == -1) // Uniform scaling
            this->m_coefs *= s;
        else if ( coord <geoDim() )//scale coordinate coord
            this->m_coefs.col(coord) *= s;
    }

    /// Apply Scaling coord-wise by a vector v
    void scale(gsVector<T> const & v)
    {
        GISMO_ASSERT( v.rows() == this->m_coefs.cols(), "Sizes do not agree." );
        this->m_coefs.array().rowwise() *= v.array().transpose();
    }

    /// Apply translation by vector v
    void translate(gsVector<T> const & v)
    {
        this->m_coefs.rowwise() += v.transpose();
    }

    /// \brief Returns the control point at corner \a c
    typename gsMatrix<T>::RowXpr
    coefAtCorner(boxCorner const & c);

    /// \brief Returns the control point at corner \a c
    typename gsMatrix<T>::ConstRowXpr
    coefAtCorner(boxCorner const & c) const;

    /// @}

    /*************************************************************************/

    /// @name Other miscellaneous functions
    /// @{

    /// Refine the geometry uniformly, inserting \a numKnots new knots into each knot span
    virtual void uniformRefine(int numKnots = 1, int mul=1) // todo: int dir = -1
    {
        this->basis().uniformRefine_withCoefs( m_coefs, numKnots, mul);
    }

    /** \brief Refines the basis and adjusts the coefficients to keep the geometry the same.
     *
     * The syntax of \em boxes depends on the implementation in the
     * underlying basis. See gsBasis::refineElements_withCoefs() for details.
     */
    void refineElements( std::vector<index_t> const & boxes )
    {
        this->basis().refineElements_withCoefs(this->m_coefs, boxes );
    }

    typename gsGeometry::uPtr coord(const index_t c) const {return this->basis().makeGeometry( this->coefs().col(c) ); }
    
    /// Embeds coefficients in 3D
    void embed3d()
    {
        embed(3);
    }

    /// Embeds coefficients in \a N dimension
    void embed(index_t N)
    { 
        GISMO_ASSERT( N > 0, "Embed dimension must be positive");

        const index_t nc = N - m_coefs.cols();

        if ( nc != 0 )
        {
            m_coefs.conservativeResize(Eigen::NoChange, N);
            if ( nc > 0 )
                m_coefs.rightCols(nc).setZero();
            else // nc < 0
            {
                gsWarn<<"Coefficients projected (deleted)..\n";
            }
        }
    }

    /// \brief Returns the degree wrt direction i
    short_t degree(const short_t & i) const
     //{ return this->basisComponent(i).degree(); };
     { return this->basis().degree(i); }

    /// \brief Elevate the degree by the given amount \a i for the
    /// direction \a dir. If \a dir is -1 then degree elevation is
    /// done for all directions
    virtual void degreeElevate(short_t const i = 1, short_t const dir = -1);

    /// \brief Reduces the degree by the given amount \a i for the
    /// direction \a dir. If \a dir is -1 then degree reduction is
    /// done for all directions
    virtual void degreeReduce(short_t const i = 1, short_t const dir = -1);
    
    /// Compute the Hessian matrix of the coordinate \a coord
    /// evaluated at points \a u
    virtual void hessian_into(const gsMatrix<T>& u, gsMatrix<T> & result,
                              index_t coord) const;
    
    /// Return the control net of the geometry
    void controlNet( gsMesh<T> & mesh) const
    { basis().connectivity(m_coefs, mesh); }

    /// @brief Computes the outer normals at parametric points \a u.
    ///
    /// Assumes that \a u is a list of points on the boundary of the geometry.
    void outerNormal_into(const gsMatrix<T>& u, gsMatrix<T> & result) const;

    /// Get boundary of this geometry as a vector of new gsGeometry instances
    std::vector<gsGeometry *> boundary() const;

    /// Get parametrization of boundary side \a s as a new gsGeometry uPtr.
    typename gsGeometry::uPtr boundary(boxSide const& s) const;

    GISMO_UPTR_FUNCTION_PURE(gsGeometry, clone)

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const
    {
        os << "Geometry "<< "R^"<< this->parDim() << 
            " --> R^"<< this->geoDim()<< ", #control pnts= "<< coefsSize() <<
            ": "<< coef(0) <<" ... "<< coef(this->coefsSize()-1); 
        os<<"\nBasis:\n" << this->basis() ;
        return os; 
    }

    /// Merge the given \a other geometry into this one.
    virtual void merge( gsGeometry * other );

    virtual void toMesh(gsMesh<T> & msh, int npoints) const;

    /// Updates the vertices of input \a mesh by evaluating the
    /// geometry at vertices.
    ///  Vertices of the new mesh are
    ///
    /// { geom(v) | v vertex of input mesh }
    ///
    void evaluateMesh(gsMesh<T>& mesh) const;


    /// Splits the geometry 2^d parts, where each direction is divided into two parts in
    /// in a uniform way, i.e., in the middle of the corresponding side. This method
    /// allocated new space for each new geometry, the original one stays unchanged.
    virtual std::vector<gsGeometry<T>* > uniformSplit(index_t dir =-1 ) const;

    /// @}

    /// Gives back an isoParametric slice of the geometry with fixed
    /// \a par in direction \a dim_fixed as an gsGeometrySlice object.
    gsGeometrySlice<T> getIsoParametricSlice(index_t dir_fixed, T par) const;

    /// Takes the physical \a points and computes the corresponding
    /// parameter values.  If the point cannot be inverted (eg. is not
    /// part of the geometry) the corresponding parameter values will be undefined
    virtual void invertPoints(const gsMatrix<T> & points, gsMatrix<T> & result,
                              const T accuracy = 1e-6,
                              const bool useInitialPoint = false) const;

    /// Returns the parameters of closest point to \a pt
    void closestPointTo(const gsVector<T> & pt,
                        gsVector<T> & result,
                        const T accuracy = 1e-6,
                        const bool useInitialPoint = false) const;

    /// Sets the patch index for this patch
    void setId(const size_t i) { m_id = i; }

    /// Returns the patch index for this patch
    size_t id() const { return m_id; }


protected:
    void swap(gsGeometry & other)
    {
        std::swap(m_basis, other.m_basis);
        m_coefs.swap(other.m_coefs);
        std::swap(m_id, other.m_id);
    }

protected:

    /// Coefficient matrix of size coefsSize() x geoDim()
    //todo: coefsSize() x (geoDim() + 1) if projective
    gsMatrix<T> m_coefs;

    /// Pointer to the basis of this geometry
    gsBasis<T> * m_basis;

    /// An auxiliary index for this geometry (eg. in case it is part
    /// of a multi-patch object)
    size_t m_id;

}; // class gsGeometry

/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsGeometry<T>& b)
{return b.print(os); }

} // namespace gismo

#include <gsCore/gsCurve.h>
#include <gsCore/gsSurface.h>
#include <gsCore/gsVolume.h>
#include <gsCore/gsBulk.h>

namespace gismo
{

// Generic traits for geometry with dimension known at compile time
template <short_t d, typename T>
struct gsGeoTraits
{
    typedef gsGeometry<T> GeometryBase;
};

// Traits for curve
template <typename T>
struct gsGeoTraits<1,T>
{
    typedef gsCurve<T> GeometryBase;
};

// Traits for surface
template <typename T>
struct gsGeoTraits<2,T>
{
    typedef gsSurface<T> GeometryBase;
};

// Traits for volume
template <typename T>
struct gsGeoTraits<3,T>
{
    typedef gsVolume<T> GeometryBase;
};

// Traits for bulk
template <typename T>
struct gsGeoTraits<4,T>
{
    typedef gsBulk<T> GeometryBase;
};

} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsGeometry.hpp)
#endif
