
#pragma once

#include <ostream>

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>
#include <gsNurbs/gsRatBezierBasis.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsUtils/gsBoehm.hpp>


namespace gismo
{

template <class T, class KnotVectorType> class gsRatBezier;
template <class T, class KnotVectorType> class gsRatBezierBasis;

/** \brief
    Class for a rational Bezier curve of one argument, with arbitrary target dimension.

    R -> R^n

    \todo Not sure if this is properly implemented, gsRatBezierBasis does not exist.

    \tparam T coefficient type

    \ingroup geometry
*/

template<class T, class KnotVectorType>
class gsRatBezier : public gsGenericGeometry< gsRatBezierBasis<T,KnotVectorType> >
{

public: 
  /// Coefficient type
  typedef T Scalar_t;

  typedef gsRationalBasis<gsBSplineBasis<T> > RBasis;
  typedef gsRatBezierBasis<T,KnotVectorType> Basis;
  typedef gsGenericGeometry<Basis> Base;

  /// Shared pointer for gsRatBezier
  typedef memory::shared_ptr< gsRatBezier<T> > Ptr;

public:

    /// Default empty constructor
    gsRatBezier() : Base() { };
    
    /// Construct NURBS by NURBS basis functions and coefficient matrix
    gsRatBezier( const Basis * basis, const gsMatrix<T> * coefs ) :
        Base( basis, coefs) { };

    /// Combatibilty constructor to avoid compiler errors, usually never called
    gsRatBezier( const RBasis * basis, const gsMatrix<T> * coefs ) 
    { assert(false); };

    /// Construct B-Spline by a knot vector, degree, weights and coefficient matrix
    gsRatBezier( const gsKnotVector<T>& KV, unsigned p, gsMatrix<T> * w, gsMatrix<T> & coefs )
      {
        this->m_basis = new Basis(KV, p, w);
        this->m_coefs.swap( coefs );
      }
    
    /// Construct B-Spline by a knot vector, degree and projective coefficient matrix
    gsRatBezier( const gsKnotVector<T>& KV, const unsigned & p, const gsMatrix<T> * pcoefs ) :
        Base( new Basis(KV,p, & pcoefs->rightCols(1) ), pcoefs) {
      // TO DO: divide pcoefs by the weights
    };
    
    /// Clone function. Used to make a copy of the (derived) geometry
    virtual gsRatBezier * clone() const
        { return new gsRatBezier(*this); };


public:

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

  /// Prints the object as a string.
  std::ostream &print(std::ostream &os) const
    { os << "NURBS curve "<< "of degree "<< this->basis().degree()
         << " over knots "<< this->basis().knots() <<",\n";
        os << "weights: ["<< this->weights().transpose()<< " ]\n ";
	os << "with control points "<< this->m_coefs << ".\n"; 
        return os; 
    };


//////////////////////////////////////////////////
// Additional members for univariate B-Splines
//////////////////////////////////////////////////

  /// Returns the starting value of the domain of the basis
  T domainStart() const { return this->basis().knots().first(); };

  /// Returns the starting value of the domain of the basis
  T domainEnd() const { return this->basis().knots().last(); };

  /// Access to i-th weight
  T & weight(int i) const { return this->basis().weight(i); }

  /// Returns the weights of the rational basis
  gsMatrix<T> & weights() const { return this->basis().weights(); }

    void isCompatible( gsGeometry<T> * other )
        {
            // compatible curves: same degree, same first/last p+1 knots
        };
    void makeCompatible( gsGeometry<T> * other )
        {
            // compatible curves: same degree, same first/last p+1 knots
        };

     void merge( gsGeometry<T> * other )
        {
            //check for BSpline
            gsRatBezier<T> *  bother = static_cast<gsRatBezier<T> *>( other );

            //check degree
            if ( this->basis().degree() != bother->basis().degree() )
            {std::cout<<"gsRatBezier: Cannot merge with different degree curve"<<"\n"; return;}
            
            //check geometric dimension

            //check that it touches *this curve

             // merge knot vector
            this->basis().knots().merge( bother->basis().knots() ) ;

            // merge coefficients
            int n= this->coefsSize();
            this->m_coefs.conservativeResize( n + bother->coefsSize() -1, Eigen::NoChange ) ;
            
            this->m_coefs.block( n,0,bother->coefsSize()-1,bother->geoDim() ) =
                bother->m_coefs.block( 1,0,bother->coefsSize()-1,bother->geoDim() ) ;

            // merge Weights
            this->weights().conservativeResize( n + bother->coefsSize() -1, Eigen::NoChange ) ;
            
            this->weights().block( n,0,bother->coefsSize()-1,1 ) =
                bother->weights().block( 1,0,bother->coefsSize()-1, 1 ) ;
        };
    
    /// Insert a knot without changing the curve, optionally, insert i times
    void insertKnot( T const& knot , int const& i = 1)
        {
            assert( i>0);
            // TO DO: There is also Oslo Algorithm and others

            std::cout <<"before \n"<< this->m_coefs.transpose() <<"\n";

            gsMatrix<T> tmp( this->m_coefs.rows(), this->m_coefs.cols()+1 );
            tmp.leftCols(this->m_coefs.cols()) =  
                this->weights().asDiagonal() * this->m_coefs;
            tmp.rightCols(1) = this->weights();

            gsBoehm( this->basis().knots(), tmp, knot, i); 

            int l = tmp.cols() -1;
            this->basis().set_weights( tmp.rightCols( 1 ) );
            for ( index_t k = 0; k< this->m_coefs.rows(); ++k)
                for ( index_t j = 0; j< this->m_coefs.cols(); ++j)
                  tmp(k,j) /= tmp(k,l);           
            //tmp.col(j).cwiseQuotient( *this->basis().m_weights );
            
            this->m_coefs =  tmp.leftCols( this->m_coefs.cols() );
        }

    
    void elevate( )
        {
            // uses knot insertion
            //gsMatrix<T> * c = this->coefs();
        };    

  

protected:
 // TODO Check function
 // check function: check the coefficient number, degree, knot vector ...


// Data members
private:

    bool projective;

}; // class gsRatBezier 


//////////////////////////////////////////////////
//////////////////////////////////////////////////




}; // namespace gismo
