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

#include <gsCore/gsCurve.h>
#include <gsCore/gsSurface.h>
#include <gsCore/gsVolume.h>
#include <gsCore/gsBulk.h>

namespace gismo
{


template <unsigned d, typename T>
struct gsGeoTraits
{
    typedef gsGeometry<T> GeometryBase;
};

template <typename T>
struct gsGeoTraits<1,T>
{
    typedef gsCurve<T> GeometryBase;
};

template <typename T>
struct gsGeoTraits<2,T>
{
    typedef gsSurface<T> GeometryBase;
};

template <typename T>
struct gsGeoTraits<3,T>
{
    typedef gsVolume<T> GeometryBase;
};

template <typename T>
struct gsGeoTraits<4,T>
{
    typedef gsBulk<T> GeometryBase;
};


/** 
    \brief
    Abstract base class representing a geometry map, that is,
    a basis together with coefficients for that basis.

    The parameter domain of the basis is also the parameter domain
    of the geometry map, and has a certain source dimension \em d
    (see gsGeometry::parDim()).
    The dimension \em n of the resulting geometry,
    i.e., of the image of the geometry map, is defined by the
    size of the given coefficients, and may be larger than \em d
    (see gsGeometry::geoDim()).

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

    Here is an overview over the different evaluation procedures provided by gsGeometry:

    Name of procedure            | Which basis function   | Evaluate what
    -----------------------------|------------------------|-------------------------------
    \c eval(u)                   | all active at u        | value
    \c evalSingle(i, u)          | basis function i       | value
    \c deriv(u)                  | all active at u        | first derivative(s)
    \c derivSingle(i, u)         | basis function i       | first derivative(s)
    \c deriv2(u)                 | all active at u        | second derivative(s)
    \c deriv2Single(i, u)        | basis function i       | second derivative(s)
    \c evalAllDers(u, k)         | all active at u        | value and all derivatives up to order k
    \c evalAllDersSingle(i, u, k)| basis function i       | value and all derivatives up to order k
    \c evalDerSingle(i, u, k)    | basis function i       | k-th derivative (k=0, ... , p-1)

    ALl evaluation functions also provide a version suffixed with \c _into
    which takes a matrix reference as an additional output parameter into which
    the result will be stored.

    Note that gsGeometry derives from gsFunction and supports all its
    evaluation functions. A gsGeometry can thus be used wherever a
    gsFunction is expected.
    
    \f[ G(u) = x \qquad u \in \mathbb R^d, x \in \mathbb R^n \f]

    \tparam T coefficient type

    \ingroup geometry
*/

template<class T>
class gsGeometry : public gsFunction<T>
{

public: 
    /// Shared pointer for gsGeometry
    typedef memory::shared_ptr< gsGeometry > Ptr;
//  typedef memory::unique_ptr< gsGeometry > LocalPtr;

    typedef T Scalar_t;

    typedef memory::auto_ptr<gsGeometryEvaluator<T> > Evaluator;

public:

    /// @name Constructors
    /// @{

    /// Default empty constructor
    gsGeometry() { }

    /// Constructor which copies the given coefficient matrix \a coefs.
    gsGeometry( const gsMatrix<T> & coefs ) :
    m_coefs( coefs )
    { }

    /// Constructor which takes ownership of the given coefficient matrix \a coefs.
    gsGeometry( gsMovable< gsMatrix<T> > coefs ) :
    m_coefs( coefs )
    { }

    /// @}

public:

    /** @name Evaluation functions

        These functions allow to evaluate the geometry as well as its derivatives
        at one or multiple points of the parameter space.
        All evaluation functions of gsFunction, from which gsGeometry derives,
        are also supported.

        \note
        These functions generally do not have to be overridden in
        derived classes since the basis type will provide the proper implementation.

        @{
    */

    /// Evaluates the geometry into result
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    { this->basis().eval_into(u, m_coefs, result); }

    /// Evaluates the derivative matrix into result
    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    { this->basis().deriv_into(u, m_coefs, result); }

    /// Evaluates the second derivatives into result
    void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    { this->basis().deriv2_into(u, m_coefs, result); }

    /// \brief Evaluates the Jacobian (matrix of first partial derivatives).
    ///
    /// Per column of \a u, result is a <em>n x d</em> matrix;
    /// result is stored in a <em>n x (d*u.cols)</em> matrix
    ///
    /// Default implementation: sum(coefs x basis_derivs).
    typename gsMatrix<T>::uPtr jac(const gsMatrix<T>& u) const
    { return this->deriv(u); }

    /// \brief Evaluates if the geometry orientation coincide with the
    /// ambient orientation.
    /// This is computed in the center of the parametrization and will
    /// fail to be meaningful if the geometry is singular.
    /// returns one if codimension is not zero
    int orientation() const
    {
        if ( parDim() == geoDim() )
        {
            const T val = jac( parameterCenter() )->determinant();
            return (T(0) < val) - (val < T(0));
        }
        return 1;
    }

    /// @}

    /*************************************************************************/

    /// @name Accessors
    /// @{

    /// \brief Returns the basis of the geometry.
    ///
    /// In the derived class, substitute the return type by the actual basis type.
    virtual gsBasis<T> & basis() const = 0;

    /// Dimension \em n of the absent physical space
    virtual int geoDim() const = 0;

    /// Dimension \em n of the coefficients (control points)
    int coefDim() const { return m_coefs.cols(); }
    
    /// Dimension of the absent physical space (overriding gsFunction::targetDim())
    int targetDim() const { return this->geoDim(); }

    /// Dimension \em d of the parameter domain.
    virtual int parDim() const { return this->basis().dim(); }

    /// Co-dimension of the geometric object
    int coDim() const { return geoDim()-this->basis().dim(); }

    /// Dimension \em d of the parameter domain (overriding gsFunction::domainDim()).
    int domainDim() const { return this->basis().dim(); }

    /// Returns the range of parameters as a matrix with two colums, [lower upper]
    virtual gsMatrix<T> parameterRange() const
    { return this->basis().support(); }

    /// Returns the range of parameters (same as parameterRange())
    virtual gsMatrix<T> support() const
    { return this->basis().support(); }

    /// Returns a "central" point inside inside the parameter domain
    virtual gsMatrix<T> parameterCenter() const
    { 
        // default impl. assumes convex support
        gsMatrix<T> S = this->basis().support();
        return ( S.col(0) + S.col(1) ) * T(0.5);
    }

    /// Whether the coefficients of this geometry are stored in projective or affine form
    virtual bool isProjective() const = 0;

    /// @}

    /*************************************************************************/

    /** @name Coefficient access functions

        These functions allow direct access to the coefficient matrix of the geometry.
        @{
    */

    /// Returns the coefficient matrix of  the geometry
    /// Coefficient matrix of size coefsSize() x geoDim(), or coefsSize() x (geoDim() + 1) if projective
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

    /// Set the coefficient matrix of the geometry
    void setCoefs( const gsMatrix<T> & cc) { this->m_coefs = cc; }

    /// Set the coefficient matrix of the geometry, taking ownership of the matrix
    void setCoefs( gsMovable< gsMatrix<T> > cc) { this->m_coefs = cc; }

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
        assert( mat.rows() == geoDim() && mat.cols() == geoDim() );
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
        assert( geoDim() == 2 );
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
        this->m_coefs.col(0) *= v[0];
        this->m_coefs.col(1) *= v[1];
        this->m_coefs.col(2) *= v[2];
    }

    /// Apply translation by vector v
    void translate(gsVector<T> const & v)
    {
        this->m_coefs.rowwise() += v.transpose();
    }

//    /// Reverse the coefficients // Now in gsCurve
//    void reverse()

    /// @}

    /*************************************************************************/

    /// @name Other miscellaneous functions
    /// @{

    /// Refine the geometry uniformly, inserting \a numKnots new knots into each knot span
    virtual void uniformRefine(int numKnots = 1)
    {
        this->basis().uniformRefine_withCoefs( m_coefs, numKnots );
    }

    /// Embeds coefficients in 3D
    void embed3d()
    { 
        const int n = 3 - m_coefs.cols();

        if ( n != 0 )
        {
            m_coefs.conservativeResize(Eigen::NoChange, 3);
            if ( n > 0 )
                m_coefs.rightCols(n).setZero();
            else
            {
                gsWarn<<"Coefficients deleted..\n";
            }
        }
    }

    /// Elevate the degree by the given amount.
    virtual void degreeElevate(int const i = 1) = 0;
    
    /// Return the control net of the geometry
    void controlNet( gsMesh<T> & mesh) const
    { basis().connectivity(m_coefs, mesh); }

    /// @brief Computes the outer normals at parametric points \a u.
    ///
    /// Assumes that \a u is a list of points on the boundary of the geometry.
    void outerNormal_into(const gsMatrix<T>& u, gsMatrix<T> & result) const;

    /// Get boundary of this geometry as a vector of new gsGeometry instances
    std::vector<gsGeometry *> boundary() const;

    /// Get parametrization of boundary side \a s as a new gsGeometry
    gsGeometry * boundary(boundary::side const& s) const;

    /// Clone function. Makes a deep copy of the geometry object.
    virtual gsGeometry * clone() const = 0;

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const
    {
        os << "Geometry "<< "R^"<< this->parDim() << 
            " --> R^"<< this->geoDim()<< ", #control pnts= "<< coefsSize() <<
            ": "<< coef(0) <<" ... "<< coef(this->coefsSize()-1); 
        os<<"\nBasis:\n" << this->basis() ;
        return os; 
    }

    /// Returns an evaluator object with the given \a flags that
    /// provides geometry-related intrinsics
    virtual gsGeometryEvaluator<Scalar_t> * evaluator(unsigned flags) const;

    /// Merge the given \a other geometry into this one.
    virtual void merge( gsGeometry * other );

    virtual void toMesh(gsMesh<T> & msh, int npoints) const;

    /// Returns new mesh, with the same connections as the input \a mesh.
    ///  Vertices of the new mesh are
    ///
    /// { geom(v) | v vertex of input mesh }
    ///
    void evaluateMesh(gsMesh<T>& mesh) const;

    /// @}


protected:

    // Coefficient matrix of size coefsSize() x geoDim(), or coefsSize() x (geoDim() + 1) if projective
    gsMatrix<T> m_coefs;

}; // class gsGeometry

//////////////////////////////////////////////////
//////////////////////////////////////////////////


/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsGeometry<T>& b)
{return b.print(os); };


////////////////////////////////////////////////////////////////////////////////


/** \brief  Base class for geometries with statically known basis types.
 *
 * This is an internal implementation class which all concrete geometries
 * derive from. It allows the compile-time specification of the associated
 * basis type, as well as some optimizations which are possible due to
 * the statically known basis dimension.
 *
 * \tparam Basis_t  type of the basis for this geometry
 *
 * \ingroup geometry
 */
template<class Basis_t> // change to template<unsigned d, class T>  ?
class gsGenericGeometry : 
        public gsGeoTraits<Basis_t::Dim, 
                           typename Basis_t::Scalar_t>::GeometryBase
{

public: 
    // obtain the dimension of the parameter domain from the basis
    static const int ParDim = Basis_t::Dim;

    /// Shared pointer for gsGenericGeometry
    typedef memory::shared_ptr< gsGenericGeometry > Ptr;

    typedef typename Basis_t::Scalar_t Scalar_t;
    typedef typename Basis_t::Scalar_t T;

    typedef typename gsGeoTraits<ParDim,Scalar_t>::GeometryBase Base;

    typedef typename gsGeometry<Scalar_t>::Evaluator Evaluator;

    // for rational bases, we WILL store the coefficients in projective form later
    static const bool IsProjective = Basis_t::IsRational;

public:

    /// Default empty constructor
    gsGenericGeometry() : Base(), m_basis(NULL) { }

    gsGenericGeometry( const Basis_t & basis, const gsMatrix<Scalar_t> & coefs )
    : Base (coefs), m_basis( basis.clone() ) 
    { 
        GISMO_ASSERT( basis.size() == coefs.rows(), 
                      "The coefficient matrix of the geometry (rows="<<coefs.rows()<<") does not match the number of basis functions in its basis("<< basis.size() <<").");
    }
    
    gsGenericGeometry( const Basis_t & basis, gsMovable< gsMatrix<Scalar_t> > coefs )
    : Base (coefs), m_basis( basis.clone() ) 
    { 
        GISMO_ASSERT( basis.size() == this->m_coefs.rows(), 
                      "The coefficient matrix of the geometry (rows="<<this->m_coefs.rows()<<") does not match the number of basis functions in its basis("<< basis.size() <<").");
    }
    
    gsGenericGeometry( const gsMatrix<Scalar_t> & coefs )
    : Base (coefs), m_basis( 0 ) 
    { }

    gsGenericGeometry( gsMovable< gsMatrix<Scalar_t> > coefs )
    : Base (coefs), m_basis( 0 ) 
    { }

    /// Copy Constructor
    gsGenericGeometry( const gsGenericGeometry & o)
    : Base (o)
    {
        m_basis = o.basis().clone() ;
    }
    
    /// Assignment operator
    gsGenericGeometry& operator=( const gsGenericGeometry & o)
    {
        if ( this == &o )
            return *this;
        gsGeometry<Scalar_t>::operator= (o);
        if (m_basis) delete m_basis;
        m_basis = o.basis().clone() ;
        return *this;
    }
    
    /// Destructor
    virtual ~gsGenericGeometry() 
    {
        if (m_basis) delete m_basis; // NULL is the default initializer, need to check
    }
    
    
public:
    
    virtual bool isProjective() const
    { return IsProjective; }
    
    virtual int geoDim() const
    { 
        // TODO
        //return ( isProjective ? this->m_coefs.cols() - 1 : this->m_coefs.cols() )
        return this->m_coefs.cols();
    }
    
    virtual int parDim() const { return ParDim; }

    static void computeJacobians(const gsMatrix<Scalar_t>& coefs, const gsMatrix<unsigned>& active, typename gsMatrix<Scalar_t>::Block derivs, gsMatrix<Scalar_t>& J);

    virtual void deriv_into(const gsMatrix<Scalar_t>& u, gsMatrix<Scalar_t>& result) const;

    virtual typename gsMatrix<T>::uPtr hess(const gsMatrix<Scalar_t>& u, 
                                            unsigned coord = 0) const;

    /// Returns the basis, as a concrete type, of the geometry.
    virtual Basis_t & basis() const { return *this->m_basis; }

    /// Elevate the degree by the given amount.
    virtual void degreeElevate(int const i = 1);

// Data members
protected:
 
    // Basis of the geometry
    Basis_t  * m_basis ;

}; // class gsGenericGeometry

template<class Basis_t>
void gsGenericGeometry<Basis_t>::computeJacobians(const gsMatrix<Scalar_t>& coefs, const gsMatrix<unsigned>& active, typename gsMatrix<Scalar_t>::Block derivs, gsMatrix<Scalar_t>& J)
{
    unsigned n = coefs.cols();
    unsigned numPts = derivs.cols();       // at how many points to evaluate the gradients

    J.setZero(n, numPts * ParDim);

    for (unsigned j = 0; j < numPts; ++j)
        for (unsigned k=0; k<n; ++k ) // for all rows of the jacobian
            for ( index_t i=0; i< active.rows() ; i++ ) // for all non-zero basis functions)
            {
                J.template block<1,ParDim>(k,j*ParDim) +=  coefs(active(i,j),k ) * derivs.template block<ParDim,1>(i*ParDim, j).transpose(); 
            }
}


template<class Basis_t>
void gsGenericGeometry<Basis_t>::deriv_into(const gsMatrix<typename Basis_t::Scalar_t>& u, gsMatrix<Scalar_t>& result) const 
{  
    gsMatrix<Scalar_t> B;
    this->basis().deriv_into(u,B);     // col j = nonzero derivatives at column point u(..,j)
    gsMatrix<unsigned> ind;
    this->basis().active_into(u,ind);    // col j = indices of active functions at column point u(..,j)
  
    this->computeJacobians(this->m_coefs, ind, B.topRows( B.rows() ), result);
}

template<class Basis_t>
typename gsMatrix<typename Basis_t::Scalar_t> ::uPtr 
gsGenericGeometry<Basis_t>::hess(const gsMatrix<typename gsGenericGeometry<Basis_t>::Scalar_t>& u, unsigned coord) const 
{  
    static const unsigned d = ParDim;

    gsMatrix<Scalar_t> B, *DD = new gsMatrix<Scalar_t>(d,d);
    gsMatrix<Scalar_t,d,d> tmp;
    gsMatrix<unsigned> ind;

    // coefficient matrix row k = coef. of basis function k
    const gsMatrix<Scalar_t>& C = this->m_coefs; 
    // col j = nonzero second derivatives at column point u(..,j)
    this->basis().deriv2_into(u, B) ; 
    // col j = indices of active functions at column point u(..,j)
    this->basis().active_into(u, ind);  
  
    DD->setZero();
    unsigned j=0;// just one column
    //for ( unsigned j=0; j< u.cols(); j++ ) // for all points (columns of u)
    for ( index_t i=0; i< ind.rows() ; i++ ) // for all non-zero basis functions)
    {
        unsigned m=i*d;
        unsigned r= ind.rows()*d + i*d*(d-1)/2;
        //construct the Hessian of basis function ind(i,0)
        for (unsigned k=0; k<d; ++k ) // for all rows
        {
            tmp(k,k) = B(m+k,j);
            for (unsigned l=k+1; l<d; ++l ) // for all cols
                tmp(k,l) = tmp(l,k) = B(r++,0);
        }
        *DD += C(ind(i,j), coord) * tmp;
    }
  
    return typename gsMatrix<T>::uPtr(DD); 
}


template<class Basis_t>
void gsGenericGeometry<Basis_t>::degreeElevate(int const i) 
{
    gsBasis<T> * b = basis().clone();
    b->degreeElevate(i);
    
    gsMatrix<T> iVals, iPts = b->anchors();
    this->eval_into(iPts, iVals);
    gsGenericGeometry<Basis_t> * g = static_cast< gsGenericGeometry<Basis_t>*>(
        b->interpolate(iVals, iPts) );

    std::swap(m_basis, g->m_basis);
    g->coefs().swap(this->coefs());

    delete g;
    delete b;
}

} // namespace gismo


#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsGeometry.hpp)
#endif
