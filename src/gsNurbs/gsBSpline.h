/** @file gsBSpline.h

    @brief Represents a B-spline curve/function with one parameter

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>

#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsBoehm.h>
#include <gsNurbs/gsBSplineSolver.h>

#include <gsUtils/gsPointGrid.h>


namespace gismo
{
    
/** \brief
    A B-spline function of one argument, with arbitrary target dimension.

    This is the geometry type associated with gsBSplineBasis.
    
    \tparam T coefficient type

    \ingroup geometry
    \ingroup Nurbs
*/
    
    
template<class T, class KnotVectorType >
class gsBSpline : public gsGenericGeometry<gsBSplineBasis<T, KnotVectorType> >
{
    
public: 
    typedef gsBSplineBasis<T, KnotVectorType> Basis;
    typedef gsGenericGeometry< gsBSplineBasis<T, KnotVectorType> > Base;

    /// Shared pointer for gsBSpline
    typedef memory::shared_ptr< gsBSpline<T, KnotVectorType> > Ptr;
    
public:
    
    /// Default empty constructor.
    gsBSpline() { }
    
    /// Construct B-Spline by basis and coefficient matrix.
    gsBSpline( const Basis & basis, const gsMatrix<T> & coefs ) :
        Base( basis, coefs )
    {
        if( this->basis().isPeriodic() )
            this->basis().expandCoefs(m_coefs);
    }
        
    /// Construct B-Spline by basis and coefficient matrix.
    gsBSpline( const Basis & basis, gsMovable< gsMatrix<T> >  coefs ) :
        Base( basis, coefs )
    {
        if( this->basis().isPeriodic() )
            this->basis().expandCoefs(m_coefs);
    }
        
    /// Construct B-Spline by a knot vector and coefficient matrix.
    gsBSpline( const KnotVectorType & KV, const gsMatrix<T> & coefs, bool periodic = false )
    {
        this->m_basis = new Basis(KV);
        m_coefs = coefs;
            
        if( periodic )
        {
            const index_t sz = this->basis().size();

            this->basis().setPeriodic();

            if ( m_coefs.rows() == sz )
            {
                this->basis().trimCoefs(m_coefs);
            }
            else if  ( m_coefs.rows() == this->basis().size() )
            {
                this->basis().expandCoefs(m_coefs);
            }
            else
            {
                GISMO_ERROR("Wrong number of coefficients for periodic basis.");
            }
        }
        else // non-periodic
        {
            if( this->m_coefs.rows() + KV.degree() + 1 != static_cast<int>( KV.size() ) )
                gsWarn << "gsBSpline Warning: #Knots="<< KV.size()<< ", #coefs="<< this->m_coefs.rows() <<"\n";
        }
    }
        
    /// Construct B-Spline by a knot vector and coefficient matrix.
    gsBSpline( const KnotVectorType & KV, gsMovable< gsMatrix<T> > coefs, bool periodic = false )
    {
        this->m_basis = new Basis(KV);
        m_coefs = coefs;

        if( periodic )
        {
            const index_t sz = this->basis().size();

            this->basis().setPeriodic();

            if ( m_coefs.rows() == sz )
            {
                this->basis().trimCoefs(m_coefs);
            }
            else if  ( m_coefs.rows() == this->basis().size() )
            {
                this->basis().expandCoefs(m_coefs);
            }
            else
            {
                GISMO_ERROR("Wrong number of coefficients for periodic basis.");
            }
        }
        else // non-periodic
        {
            if( this->m_coefs.rows() + KV.degree() + 1 != int( KV.size() ) )
                gsWarn << "gsBSpline Warning: #Knots="<< KV.size()<< ", #coefs="<< this->m_coefs.rows() <<"\n";
        }
    }
    
    /// @brief Construct a B-spline from an interval and knot vector specification.
    /// \param u0 starting parameter
    /// \param u1 end parameter
    /// \param interior number of interior knots
    /// \param degree degree of the spline space
    /// \param coefs coefficients of the spline space
    /// \param mult_interior multiplicity at the interior knots
    ///
    /// \ingroup Nurbs
    gsBSpline(T u0, T u1, unsigned interior, int degree, gsMovable< gsMatrix<T> > coefs, unsigned mult_interior=1, bool periodic = false)
    {
        this->m_basis = new Basis(u0, u1, interior, degree, mult_interior, periodic );
        if( periodic )
        {
            GISMO_ASSERT( this->basis().numCrossingFunctions() == 1,
                          "Knot vector with clamped knots should lead to m_periodic == 1.");

            this->m_coefs = this->basis().perCoefs(coefs);
            GISMO_ASSERT( this->m_coefs.rows() >= this->basis().size(),
                          "You should give at least as many coefficients as the size of the basis after converting to periodic.");
            if( this->m_coefs.rows() > this->basis().size() )
                gsWarn << "Some of the last coefficients are going to be lost during conversion to periodic basis.\n";
        }
        else // non-periodic
        {
            this->m_coefs = coefs;
            GISMO_ASSERT( this->m_coefs.rows() == this->basis().size(),
                          "Number of coefficients does not match the size of the basis.");
        }
    }
    
    /// Clone function. Used to make a copy of the (derived) geometry
    virtual gsBSpline * clone() const
        { return new gsBSpline(*this); }
    

    gsBSpline( gsBSpline const & o ) : Base(o) { }
    
    ~gsBSpline() { } //destructor
    
    
public:
    
//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////
    
    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "BSpline curve "<< "of degree "<< this->basis().degree()<< ", "<<  this->basis().knots() <<".\n";
        os << "with control points "<< this->m_coefs.row(0)<<" ... "<<this->m_coefs.bottomRows(1) << ".\n";
        if( this->basis().isPeriodic() )
            os << "Periodic with overlay " << this->basis().numCrossingFunctions() << ".\n";
        return os;
    }
    
//////////////////////////////////////////////////
// Additional members for univariate B-Splines
//////////////////////////////////////////////////

    /// Returns the starting value of the domain of the basis
    T domainStart() const { return this->basis().knots().first(); }
    
    /// Returns the end value of the domain of the basis
    T domainEnd() const { return this->basis().knots().last(); }
    
    /// Returns a reference to the knot vector
    KnotVectorType & knots()             { return this->basis().knots(); }
    /// Returns a reference to the knot vector
    const KnotVectorType & knots() const { return this->basis().knots(); }

    /// Returns the degree of the B-spline
    int degree() const { return this->basis().degree(); }
    
    // Angelos, shall I implement the following two? It should be analogous to conversion to periodic.
    // compatible curves: same degree, same first/last p+1 knots
    void isCompatible( gsGeometry<T> * other )
      { GISMO_NO_IMPLEMENTATION }

    // compatible curves: same degree, same first/last p+1 knots
    void makeCompatible( gsGeometry<T> * other )
      { GISMO_NO_IMPLEMENTATION }

    /// Merge other B-spline into this one.
    void merge( gsGeometry<T> * other )
    {
        GISMO_ASSERT( this->basis().isPeriodic() == false, "Cannot merge a closed curve with anything." );
        //check for BSpline
        gsBSpline *  bother = static_cast<gsBSpline *>( other );

        //check degree
        bool elevateOther = false;
        int thisDegree = this->basis().degree();
        int otherDegree = bother->basis().degree();
        if (thisDegree > otherDegree)
        {
            // make a degree elevated copy of the other spline
            elevateOther = true;
            bother = bother->clone();
            bother->elevate(thisDegree - otherDegree);
        }
        else if (otherDegree > thisDegree)
        {
            // degree elevate this spline in place, since we are merging
            // in place.
            elevate(otherDegree - thisDegree);
        }

        //check geometric dimension

        //check that it touches *this curve

         // merge knot vector
        this->basis().knots().merge( bother->basis().knots() ) ;

        // merge coefficients
        int n= this->coefsSize();
        this->m_coefs.conservativeResize( n + bother->coefsSize() -1, Eigen::NoChange ) ;

        this->m_coefs.block( n,0,bother->coefsSize()-1,bother->geoDim() ) =
            bother->m_coefs.block( 1,0,bother->coefsSize()-1,bother->geoDim() ) ;

        if(elevateOther)
        {
            delete bother;
        }
    }

    /// Insert the given new knot (multiplicity \a i) without changing the curve.
    void insertKnot( T knot, int i = 1)
    {
        assert( i>0);
        //if ( i==1)
        //single knot insertion: Boehm's algorithm
        //else
        //knot with multiplicity:   Oslo algorithm
        if( this->basis().isPeriodic() )
        {
            int borderKnotMult = this->basis().borderKnotMult();
            KnotVectorType & knots = this->knots();
            unsigned deg = this->basis().degree();


            GISMO_ASSERT( knot != knots[deg] && knot != knots[knots.size() - deg - 1],
                          "You are trying to increase the multiplicity of the p+1st knot but the code is not ready for that.\n");

            // If we would be inserting to "passive" regions, we rather insert the knot into the mirrored part.
            // Adjustment of the mirrored knots is then desirable.
            if( knot < knots[deg - borderKnotMult + 1] )
            {
                knot += this->basis()._activeLength();
            }
            else if( knot > knots[knots.size() - deg + borderKnotMult - 2] )
            {
                knot -= this->basis()._activeLength();
            }
            // If necessary, we update the mirrored part of the knot vector.
            if((knot < knots[2*deg + 1 - borderKnotMult]) || (knot >= knots[knots.size() - 2*deg - 2 + borderKnotMult]))
                this->basis().enforceOuterKnotsPeriodic();

            // We copy some of the control points to pretend for a while that the basis is not periodic.
            //gsMatrix<T> trueCoefs = this->basis().perCoefs( this->coefs() );
            gsBoehm( this->basis().knots(), this->coefs(), knot, i );
            //this->coefs() = trueCoefs;
            //this->coefs().conservativeResize( this->basis().size(), this->coefs().cols() );
        }
        else // non-periodic
            gsBoehm( this->basis().knots(), this->coefs() , knot, i);
    }

    /// Insert the given new knots in the range \a [begin..end)
    /// without changing the curve.
    /// \param begin iterator pointing to the first knot to be inserted
    /// \param end iterator pointing to one position after the last
    /// knot to be inserted
    template <class It>
    void insertKnots( It begin, It end)
    {
        if( this->basis().isPeriodic() )
        {
            // We assume we got valid (i.e., non-NULL) iterators; I don't think we have a reasonable way to test it in GISMO_ASSERT.

            GISMO_ASSERT( (*begin > *(this->knots().begin()))
                          && (*(end-1) < *(this->knots().end()-1)),
                          "Please, ask me to insert knots inside the knot interval only." );

            // It can happen that user would ask us to insert knots outside the active range.
            // Should it happen, we shift the values and then sort the knot vector to remain non-decreasing.

            T activeLength = this->basis()._activeLength();
            T blue1 = *(this->basis().knots().begin()  + this->degree() );
            T blue2 = *(this->basis().knots().end()    - this->degree() );
            for( It it = begin; it != end; ++it )
            {
                if( *it < blue1 )
                    *it += activeLength;
                else if( *it > blue2 )
                    *it -= activeLength;
                else if( (*it == blue1) || (*it == blue2) )
                {
                    // I would like to erase it and go on but for that the iterator is not enough.
                    // ANGELOS, wouldn't it be better to put a GISMO_ASSERT here instead? Or throw an exception?
                    gsWarn << "Interrupting insertKnots(): Cannot insert value equal to the p+1st knot from the beginning or end into a periodic curve.\n";
                    return;
                }
            }
            /* // This way it would be nicer but it does not cover all the dubious cases.
            for( It it = begin; *it < blue1; ++it )
                *it += activeLength;
            for( It it = end-1; *it > blue2; --it )
                *it -= activeLength;*/

            std::sort( begin, end );

            // We copy some of the control points to pretend for a while that the basis is not periodic.
            //gsMatrix<T> trueCoefs = this->basis().perCoefs( this->coefs() );
            gsBoehmRefine( this->basis().knots(), this->coefs(), this->degree(), begin, end );
            //this->basis().enforceOuterKnotsPeriodic();
            //this->coefs() = trueCoefs; // Isn't this memory-leaking?
            // Periodic basis is restored using the so-called Student's Algorithm: forget as much as possible =).
            //this->coefs().conservativeResize( this->basis().size(), this->coefs().cols() );
        }
        else // non-periodic
            gsBoehmRefine( this->basis().knots(), this->coefs(), this->degree(), begin, end);
    }

    void elevate( int r )
        {
            // uses knot insertion
            //gsMatrix<T> * c = this->coefs();
            bspline::degreeElevateBSpline(this->basis(), this->coefs(), r);
        }    

    /// @brief Returns true iff the point p is contained (approximately) on
    /// the curve, with the given tolerance.
    // Under Construction..( TO DO )
    bool contains( gsMatrix<T> const & p, T const & tol = 1e-6 )
        {
            assert( p.cols()==1 );
            gsBSplineSolver<T> slv;
            std::vector<T> roots;
            int dim = this->geoDim();
            gsMatrix<T> ev, tmp(1,1);
            int i(1);

            slv.allRoots( *this, roots, 0, p(0,0), 1e-6, 100);
            for ( typename std::vector<T>::const_iterator it=roots.begin();
                  it != roots.end(); ++it)
            {
                tmp(0,0)= *it;
                this->eval_into(tmp,ev);
                for (i=1; i!=dim; ++i )
                    if ( math::abs( ev(i,0) - p(i,0) ) > tol ) 
                        break;
                if (i==dim) 
                    return true;
            }
            return false;
        }

    /// Sample \a npoints uniformly distributed (in parameter domain) points on the curve.
    typename gsMatrix<T>::uPtr sample(int npoints = 50) const
    {      
        typename gsMatrix<T>::uPtr images ( new gsMatrix<T>() );
        gsMatrix<T> interval = this->parameterRange();
        gsMatrix<T> pts = gsPointGrid( interval(0,0), interval(0,1), npoints );
        this->eval_into( pts, *images );
        return images;
    }

    /// Tries to convert the curve into periodic.
    void setPeriodic(bool flag = true)
    {
        this->basis().setPeriodic(flag);
        this->m_coefs = this->basis().perCoefs( this->m_coefs );
    }

    inline int numCoefs() const { return this->m_coefs.rows() - this->basis().numCrossingFunctions(); }

    inline bool isPeriodic() const { return this->basis().isPeriodic(); }


    /// Returns true if this curve is closed
    bool isClosed(T tol = 1e-10) const 
    { 
        return ( this->basis().isPeriodic() || 
                 (m_coefs.row(0) - m_coefs.row(m_coefs.rows()-1)).squaredNorm()<tol );  
    }

    /// \brief Return true if point \a u is on the curve with
    /// tolerance \a tol
    bool isOn(gsMatrix<T> const &u, T tol = 1e-3) const;

    /// \brief Return true if point \a u is an endpoint (corner) of
    /// the curve with tolerance \a tol
    bool isPatchCorner(gsMatrix<T> const &v, T tol = 1e-3) const;

    /// \brief Modifies the parameterization such that the point \a v
    /// is the starting point of the curve. Assumes that \a v is
    /// either the starting or the ending point of the curve
    void setOriginCorner(gsMatrix<T> const &v);

    /// \brief Modifies the parameterization such that the point \a v
    /// is the ending point of the curve. Assumes that \a v is
    /// either the starting or the ending point of the curve
    void setFurthestCorner(gsMatrix<T> const &v);
    
protected:
    
    using Base::m_coefs;
    
    // TODO Check function
    // check function: check the coefficient number, degree, knot vector ...

}; // class gsBSpline
    

} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBSpline.hpp)
#endif
