
#pragma once

#include <ostream>

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>
#include <gsThbs/gsTHBSplineBasis.h>
#include <gsCore/gsForwardDeclarations.h>

namespace gismo
{

/** \brief
    A truncated hierarchical B-Spline function, in \em d dimensions.

    This is the geometry type associated with gsTHBSplineBasis.

    R^d -> R

    \tparam T is the coefficient type

    \ingroup geometry
*/

template<unsigned d, class T>
class gsTHBSpline : public gsGenericGeometry<gsTHBSplineBasis<d,T> >
{
public:
  typedef gsTHBSplineBasis<d,T> Basis;
  typedef gsGenericGeometry< gsTHBSplineBasis<d,T> > Base;

  /// Shared pointer for gsHBSpline
  typedef memory::shared_ptr< gsTHBSpline<d,T> > Ptr;

public:

    /// Default empty constructor
    gsTHBSpline() { }

    /// Construct THB-Spline by basis functions and coefficient matrix
    gsTHBSpline( const Basis & basis, const gsMatrix<T> & coefs ) :
        Base( basis, coefs ) 
    {
        
    }

    /// Construct THB-Spline by basis functions and coefficient matrix
    gsTHBSpline( const Basis& basis, gsMatrix<T> & coefs ) :
        Base( basis, coefs ) 
    {

    }

    /// Construct B-Spline from a Tensor B-Spline
  gsTHBSpline( const gsTensorBSpline<d,T> & tbsp )
        {
          this->m_basis = new Basis(tbsp.basis(), 3);// 3 levels
          this->m_coefs = tbsp.coefs();
        }

    /// Copy constructor
    gsTHBSpline( const gsTHBSpline & other )
    {
        this->m_basis = other.basis().clone();
        this->m_coefs = other.coefs();
    }

    /// Clone the gsHBspline
    virtual gsTHBSpline * clone() const
        { return new gsTHBSpline(*this); }

    ~gsTHBSpline() { } //destructor

    //void deriv2_into(const gsMatrix<T> &u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const;

public:

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

    /// Returns the degree wrt direction i
    unsigned degree(const unsigned & i) const
     //{ return this->basisComponent(i).degree(); };
     { return this->basis().component(i).degree(); }

    /** Refines the basis and adjusts the coefficients to keep the geometry the same.
     * The indices of the boxes are the same as in gsHTensorBasis<>::refineElements_withCoefs.
     */
    void refineElements( std::vector<unsigned> const & boxes )
    {
        gsMatrix<> & coefs = this->m_coefs;
        this->basis().refineElements_withCoefs( coefs, boxes );
    }


//////////////////////////////////////////////////
/// Other member functions
//////////////////////////////////////////////////
public:
    ///get all the B-spline patches out of a THB-spline geometry
    //void getBsplinePatches(gsMatrix<unsigned>& b1, gsMatrix<unsigned>& b2, gsVector<unsigned>& level, std::vector< gsTensorBSpline<2> > & bpatches) const;
    void getBsplinePatches(gsMatrix<unsigned>& b1, gsMatrix<unsigned>& b2, gsVector<unsigned>& level) const;

    /// Refines the whole domain to the finest level present in the mesh. Returns the refined geometry as result.
    void convertToBSpline( gsTensorBSpline<d,T,gsCompactKnotVector<T> >& result );

private:
    ///get B-spline control points on a given box of a certain level by refining eveywhere
    void getBsplinePatchGlobal(gsVector<unsigned> b1, gsVector<unsigned> b2, unsigned l, gsTensorBSpline<2> geo) const;
    ///function for getBsplinePatchGlobal
    void globalRefinement(int level)const;
    ///initialization of cmatrix
    void initialize_cmatrix(int col, int c_level) const;
    ///convert the coefficient matrix mat in the given direction to a column of the control points matrix
    void return_cp_1D(const gsMatrix<T> & mat, int direction, gsMatrix<T>& cp)const;
}; // class gsTHBSpline


//////////////////////////////////////////////////
//////////////////////////////////////////////////




}; // namespace gismo


#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsTHBSpline.hpp)
#endif
