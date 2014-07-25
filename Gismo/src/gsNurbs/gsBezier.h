// defined over a Bernstein basis
// list of coefficients with maching ends


#pragma once

#include <ostream>

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>
#include <gsNurbs/gsBernsteinBasis.h>


namespace gismo
{

/** \brief
    A sequence of Bezier curves of one argument, with arbitrary target dimension.

    The coefficients are a sequence of Bezier branches.
    For degree p, we get p*t+1 control points, where t is
    the number of branches.

    This is the geometry type associated with gsBernsteinBasis.
    
    \tparam T coefficient type

    \ingroup geometry
*/

template<class T>
class gsBezier : public gsGenericGeometry< gsBernsteinBasis<T> >
{

public:
    typedef gsBernsteinBasis<T> Basis;
    typedef gsGenericGeometry<Basis> Base;
    
    /// Default empty constructor
    gsBezier() : Base() { }
    
    /// Construct B-Spline by basis functions and coefficient matrix
    gsBezier( const Basis & basis, const gsMatrix<T> & coefs ) :
        Base( basis, coefs) { }
    
    /// Construct B-Spline by basis and coefficient matrix.
    gsBezier( const Basis & basis, gsMovable< gsMatrix<T> >  coefs ) :
        Base( basis, coefs )
        { }

    /// Construct Bezier curve by degree, domain and coefficient matrix
    gsBezier( const unsigned & p, const gsMatrix<T> & coefs, const T & u0=0, const T & u1= 1)
    {
        this->m_basis = new Basis(u0,u1,p);
        this->m_coefs = coefs;
    }
    
    /// Construct Bezier curve by degree, domain and coefficient matrix
    gsBezier( const unsigned & p, gsMovable< gsMatrix<T> >  coefs, const T & u0=0, const T & u1= 1)
    {
        this->m_basis = new Basis(u0,u1,p);
        this->m_coefs = coefs;
    }
    
    ~gsBezier() { }; //destructor
    

public:

//////////////////////////////////////////////////
// Virtual base members with a new implementation
//////////////////////////////////////////////////

  // Evaluates the geopetry  TO DO DeBoor algorithm
  // For now default implementation on gsGeometry
  //gsMatrix<T,n> * eval(gsMatrix<T,1>& u) const; 
  
  // Evaluates the jacobian 
  // Implementation on gsGeometry
  // defined in base class gsMatrix<T> * jac(gsMatrix<T,1>& ) const ;
  
  // Evaluates the hessian matrix
  // Implementation on gsGeometry TO DO
  //gsMatrix<T> * hess(const gsMatrix<T,d>& u) const { return NULL; } ;

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

  /// Clone function. Used to make a copy of the (derived) geometry
  virtual gsBezier * clone() const
    { return new gsBezier(*this); };

  /// Prints the object as a string.
  std::ostream &print(std::ostream &os) const
    { os << "Bezier curve "<< "of degree "<< 
            this->basis().degree()<< ", over [";
        for ( typename gsKnotVector<T>::const_iterator itr= 
                  this->basis().domain()->begin(); itr != this->basis().domain()->end(); ++itr )
            os << *itr << " ";
        os << "].";
        return os;
    };


//////////////////////////////////////////////////
// Additional members for univariate B-Splines
//////////////////////////////////////////////////


    T domainStart() const { return this->basis().m_breaks->first(); };
    
    /// Returns the starting value of the domain of the basis
    T domainEnd() const { return this->basis().m_breaks->clast(); };


    /// Elevate the degree of the curve by 1
    void elevate( )
        {

            gsMatrix<T> & c = this->coefs();
            int p = this->basis().degree() ;
            gsMatrix<T> t1(1,this->geoDim()), t2(1,this->geoDim());
            
            c.conservativeResize( p+2 , this->geoDim() ) ;
            c.row(p+1) = c.row(p) ;
            
            p++;
            t1 = t2 = c.row(0);
            for ( int i = 1; i< p; ++i )
            {
                t2 = c.row(i) ;
                c.row(i) = ( T(i)/T(p) )*t1 + ( T(p-i)/T(p) )* c.row(i) ;
                std::swap(t1, t2) ;
            }
            
            this->basis().setDegree(p) ;
            this->basis().knots().push_front( this->basis().knots().first() ) ;
            this->basis().knots().push_back( this->basis().knots().last() ) ;
        };    


protected:
 // TODO Check function
 // check function: check the coefficient number, degree, knot vector ...


// Data members
private:

}; // class gsBezier


//////////////////////////////////////////////////
//////////////////////////////////////////////////




}; // namespace gismo
