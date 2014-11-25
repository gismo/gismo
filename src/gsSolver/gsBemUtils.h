/** @file gsBemUtils.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Falini, A. Mantzaflaris
*/

#pragma once


#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsFunction.h>
//#include <gsCore/gsFunctionExpr.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsModeling/gsPlanarDomain.h>

namespace gismo
{


// Class for representing the solution of a Bem problem
template<class T>
class gsGreenFunction2d : public gsFunction<T>
{
public:

    gsGreenFunction2d()  { }
    
    ~gsGreenFunction2d() { }

public:
    // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:

  /// Clone function. Used to make a copy of the object
  gsGreenFunction2d * clone() const
    { return new gsGreenFunction2d; };


    /// Prints the object as a string.
  std::ostream &print(std::ostream &os) const
    {
      os << "gsGreenFunction is ....\n"; return os;

    }

    /// Evaluate the function at all columns of \a u into \a result
  void setSourcePoint(const gsMatrix<T>& u) const
    {
      source = u;
    }

    /// Evaluate the function at all columns of \a u into \a result
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( 1, u.cols() );
        
        for ( index_t i=0; i!= u.cols(); ++i )
        {
            // (y-w) / (2*pi*( (x-z)^2 + (y-w)^2 ) )
            result(0,i)= log( (u.col(i) - source).norm() )/ (2*3.141592653589793238462) ;
        }

    }

    /// Evaluate the function at all columns of \a u into \result
    void grad_into(const gsMatrix<T>& u, gsMatrix<T>& result) const

//            // (x-z) / (2*pi*( (x-z)^2 + (y-w)^2 ) )
//            //result(0,i)= ( u(0,i) - source(0,0) ) / tmp ;
//            // (y-w) / (2*pi*( (x-z)^2 + (y-w)^2 ) )
//            //result(1,i)= ( u(1,i) - source(1,0) )/ tmp;

    {
        result.resize( 2, u.cols() );
        
        for ( index_t i=0; i < u.cols(); ++i )
        {
            const T dst = ( u.col(i) - source ).norm();
            const T tmp = 2*3.141592653589793238462 * ( dst*dst );
            // (x-z) / (2*pi*( (x-z)^2 + (y-w)^2 ) )
            result.col(i) = (u.col(i) - source) / tmp ;
        }
    
        }


    /// Evaluate the green function wrt the collocation point
    void source_grad_into(const gsMatrix<T>& u, gsMatrix<T>& result) const

//        const T tmp = 2*3.141592653589793238462* (u-source).squaredNorm();
//        // - (x-z) / (2*pi*( (x-z)^2 + (y-w)^2 ) )
//        result(0,i)= ( source(0,0)-u(0,i) )/ tmp;
//        // - (y-w) / (2*pi*( (x-z)^2 + (y-w)^2 ) )
//        result(1,i)=  (source(1,0)- u(1,i)  )/ tmp;

    {
        result.resize( 2, u.cols() );

        for ( index_t i=0; i!= u.cols(); ++i )
        {
            const T tmp = 2*3.141592653589793238462* (u.col(i)-source).squaredNorm();
            // - (x-z) / (2*pi*( (x-z)^2 + (y-w)^2 ) )
            result.col(i) =  ( source - u.col(i) ) / tmp;
        }
    }


    /// Evaluate the function
    gsMatrix<T> * source_mderiv(const gsMatrix<T>& u) const

    {
        gsMatrix<T> * result = new gsMatrix<T>;
        this->source_mderiv_into(u, *result);
        return result;
    }
    
    /// Evaluate the green mixed derivatives wrt the collocation point
    /// result.row(0) is the mixed derivative of G respect y1x1
    /// result.row(1) is the mixed derivative of G respect y1x2
    /// result.row(2) is the mixed derivative of G respect y2x1
    /// result.row(3) is the mixed derivative of G respect y2x2
    void source_mderiv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        using std::pow;
        result.resize( 4, u.cols() );

        for ( index_t i=0; i!= u.cols(); ++i )
        {
            
            const T tmp= (u.col(i)-source).squaredNorm() ;
            
            result(0,i)= 1/3.141592653589793238462*(-1/(2*tmp) + ( (u(0,i)-source(0,0))*(u(0,i)-source(0,0)) ) /(tmp*tmp) );
            //            
            result(1,i)= ( (u(0,i)-source(0,0)) *( u(1,i)-source(1,0) ) )/(tmp*tmp*3.141592653589793238462) ;
            //
            result(2,i)= ( ( u(0,i)-source(0,0) )*( u(1,i)-source(1,0) ))/(tmp*tmp*3.141592653589793238462) ;
            //
            result(3,i)= 1/3.141592653589793238462*(-1/(2*tmp) +( (u(1,i)-source(1,0))*(u(1,i)-source(1,0)) )/(tmp*tmp) );
        }


        }


    const gsMatrix<T,2,1> & getSource() const {return source;}
    
protected:
    
    // Source point
    mutable gsMatrix<T,2,1> source;
};



// Class for representing the solution of a Bem problem
template<class T>
class gsBemSolution : public gsFunction<T>
{
public:

    gsBemSolution(gsCurveLoop<T> * loop, gsGeometry<T> * flux,
                  gsFunction<T> * boundary_fun,
                  std::vector<T> & breaks,
                  bool par_rhs= false)

    : m_pdomain(loop),  m_par_rhs(par_rhs)

    {
        m_flux.push_back(flux);

        m_boundary_fun.push_back(boundary_fun);

        m_green = new gsGreenFunction2d<T>;

        m_breaks.push_back(
                flux->basis().domain()->breaks()
                );
        
     //  m_breaks.swap(breaks);
  }

    gsBemSolution(gsPlanarDomain<T> * domain, std::vector<gsGeometry<T> *> const & flux,
                  std::vector<gsFunction<T> *> const & boundary_fun,
                  std::vector< std::vector<T> > & breaks,
                  bool par_rhs= false)
        : m_pdomain( domain->clone() ), m_par_rhs(par_rhs)
    {
        m_flux=flux;

        m_boundary_fun=boundary_fun;

        //m_breaks.swap(breaks);
        
        for ( std::size_t i= 0; i< flux.size(); ++i )
            m_breaks.push_back(
                flux[i]->basis().domain()->breaks()
                );

        m_green = new gsGreenFunction2d<T>;
    }
    

    ~gsBemSolution()
    {
        delete m_green;
        delete m_pdomain;
        freeAll(m_flux);
        freeAll(m_boundary_fun);
    }
    
    /// Clone function. Used to make a copy of the object
    gsBemSolution * clone() const
    {
        // TO DO
        return NULL;
    }
    
    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "gsBemSolution.\n";
        //os<< "Integration points: "<< m_ngrid->cols() <<"\n";
        //os<< "Integration points: "<< m_ngrid <<"\n";
        return os;
    }

public:


    std::vector<gsGeometry<T> *> flux() { return m_flux; }

    
    gsPlanarDomain<T> *getDomain ()  {return m_pdomain;}

//<<<<<<< .mine



//  }

//    gsBemSolution(gsPlanarDomain<T> * domain, std::vector<gsGeometry<T> *> const & flux,
//                  std::vector<gsFunction<T> *> const & boundary_fun,
//                  std::vector< std::vector<T> > & breaks,
//                  //std::vector<gsMatrix<T> *> const & ngrid, std::vector<gsVector<T> *> const & wgrid,
//                  bool par_rhs= false)
//        : m_pdomain( domain->clone() ), m_par_rhs(par_rhs)
//  {
//     m_flux=flux;
//     m_boundary_fun=boundary_fun;
//     //m_ngrid=ngrid;
//     //m_wgrid=wgrid;
//     m_breaks.swap(breaks);

//     m_green = new gsGreenFunction2d<T>;
//  }

//  ~gsBemSolution()
//  {
//      delete m_green;
//      delete m_pdomain;
//      freeAll(m_flux);
//      freeAll(m_boundary_fun);
//  }
//=======
    const gsGreenFunction2d<T> & getGreen () const { return *m_green; }

    const std::vector<gsFunction<T> *> & getBoundary_fun() const {return m_boundary_fun; }


//public:

//    std::vector<gsGeometry<T> *> flux() { return m_flux; }

//    gsPlanarDomain<T> *getDomain ()  {return m_pdomain;}

//    const gsGreenFunction2d<T> & getGreen () const { return *m_green; }

//    const std::vector<gsFunction<T> *>  getBoundary_fun() const {return m_boundary_fun; }

//    const std::vector<gsMatrix<T> *>  getNgrid() const {return m_ngrid;}

//    const std::vector<gsVector<T> *>  getWgrid() const {return m_wgrid;}

    bool  getPar_rhs() {return m_par_rhs;}

  int targetDim() const  { return 1;}

    /// Evaluate the function at all columns of \a u into \result
   virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {

        GISMO_ASSERT( u.rows()==2, "gsBemSolution accepts points with two coordinates." );
        //GISMO_ASSERT( m_par_rhs && ,  "" );

        std::vector<gsGeometry<T>*> all_loops;

        const index_t nloops = m_flux.size();
        result.setZero(1, u.cols());

        gsMatrix<T> quNodes  ;
        gsVector<T> quWeights;
        gsVector<index_t>  numNodes(1);
        gsGaussRule<T> QuRule;

        T res_iv;


        // variables for unit normal and physical point

        gsMatrix<T> fluxEv;

       // gsVector<T> unormal;

        for(index_t v=0; v != nloops; v++)
        {
            // Get a single curve representing every loop
            all_loops.push_back(  m_pdomain->loop(v).singleCurve() );
        }

        gsMatrix<T> green_val, green_grad, bd_val;
      //  gsInfo<<"\n what is the value of m_par_rhs ? "<<m_par_rhs<<"\n";


        for( index_t i=0; i< u.cols(); ++i )
        {
        
    // first time: computing integral over (0,1) and evaluate the solution u
            for(index_t v=0; v< nloops; v++)

            {
                T parValue;
                gsMatrix<T> MparValue(1,1);
                MparValue.setZero();
                //if ( m_pdomain->onBoundary( u.col(i) ) ) // (!) v
                if ( m_pdomain->loop(v).isOn( u.col(i),parValue ) ) // (!) v
                {
                  //m_boundary_fun[v]->eval_into( m_par_rhs ?
                  // *all_loops[v]->eval(u.col(i)) :  u.col(i), bd_val );
                  // result(0,i) = bd_val.value();
                    MparValue(0,0)=parValue;
                    m_boundary_fun[v]->eval_into( m_par_rhs ?
                    MparValue :  u.col(i), bd_val );
                    result(0,i) = bd_val.value();

                    //result(0,i) *= 2.0;
                    
                    break;
                }

                numNodes[0] = m_flux[v]->basis().minDegree()+1;
                QuRule.setNodes(numNodes);
                
                //taking breaks points :
               const std::vector<T> & breaks = m_breaks[v];

                // Evaluator on the current loop
                typename gsGeometry<T>::Evaluator loop_v_eval (
                    all_loops[v]->evaluator(NEED_VALUE   |
                                            NEED_JACOBIAN) );
                
                m_green->setSourcePoint( u.col(i) ) ;

                // Start iteration over sub-elements
                for (std::size_t e=1; e!= breaks.size(); ++e)
                {

                    // Map the Quadrature rule to the element
                    QuRule.mapTo(breaks[e-1], breaks[e], quNodes, quWeights );

                    res_iv = evalOnElement( quNodes, quWeights, v, loop_v_eval, 
                                            green_val, green_grad, bd_val, fluxEv);

// /*
                    // Adaptive step
                    std::stack<iInterval> iStack;
                    iStack.push(iInterval(breaks[e-1], breaks[e], res_iv));
                    T lVal, rVal;
                    res_iv = 0;
                    do

                    {

                        const iInterval curr = iStack.top();
                        iStack.pop();


                        const T mid = (curr.left +  curr.right ) /2.0;

                        // Map the Quadrature rule to the element
                        QuRule.mapTo(curr.left, mid, quNodes, quWeights );
                        
                        lVal = evalOnElement( quNodes, quWeights, v, loop_v_eval, 
                                              green_val, green_grad, bd_val, fluxEv);


                        // Map the Quadrature rule to the element
                        QuRule.mapTo(mid, curr.right, quNodes, quWeights );
                        
                        // Evaluate loop values and derivatives
                        loop_v_eval->evaluateAt(quNodes);


                        rVal = evalOnElement( quNodes, quWeights, v, loop_v_eval, 
                                              green_val, green_grad, bd_val, fluxEv);


                        if ( math::abs(lVal + rVal - curr.value) > 1e-3  )
                        {
                            iStack.push( iInterval(curr.left, mid , lVal) );
                            iStack.push( iInterval(mid, curr.right, rVal) );
                        }
                        else
                            res_iv += lVal + rVal;
                    }
                    while ( ! iStack.empty() ); // Adaptive quadrature loop
//*/
                result(0,i) += res_iv;


                } // Element loop
                
            }// Next point 
        }
        
        freeAll(all_loops);
    }
    

//          result(0,i)= result_new(0,i);
//*/

//  ???????????????????????????????????????????????????????????????
//          if ( m_pdomain->onBoundary( u.col(i) ) )
//              result(0,i) *= 2.0;

//          }//goes to the next point

//      }

//      freeAll(all_loops);
//  }


    /// Evaluate the derivative of the solution
    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.setZero( 1, 2*u.cols() );

        std::vector<gsGeometry<T>*> all_loops;
        gsMatrix<T> Gmderiv(4,1);
        gsMatrix<T> Gderiv(2,1);
        gsMatrix<T> quNodes, bd_val, fev  ;
        gsVector<T> quWeights, unormal;
        gsVector<index_t>  numNodes(1);
        gsGaussRule<T> QuRule;

        const index_t nloops = m_flux.size();

        gsVector<T,2> res_iv(2);

        for(index_t v=0; v != nloops; v++)

        {
            // Get a single curve representing every loop
            all_loops.push_back(  m_pdomain->loop(v).singleCurve() );
        }


        for(index_t v=0; v< nloops; v++)
        {

            numNodes[0] = m_flux[v]->basis().minDegree()+1; // num of nodes over which we want to compute the integral
            QuRule.setNodes(numNodes);

            //taking breaks points :
            const std::vector<T> & breaks = m_breaks[v];// taking breaks for every loop

            // Evaluator on the current loop
            typename gsGeometry<T>::Evaluator loop_v_eval (
                all_loops[v]->evaluator(NEED_VALUE   |
                                        NEED_GRAD_TRANSFORM) );


            for ( index_t i=0; i!= u.cols(); ++i )
            {
                // Set source point in the Green's function
                m_green->setSourcePoint( u.col(i) ) ;

                //loop over breaks for generating gauss point in every span
                for ( std::size_t e=1; e!= breaks.size(); ++e)
                {
                    // Map the Quadrature rule to the element
                    QuRule.mapTo(breaks[e-1], breaks[e], quNodes, quWeights );

                    res_iv = derivOnElement(quNodes, quWeights, v, loop_v_eval,
                                            Gmderiv, Gderiv, bd_val, fev);

// 	   else
// 	   {
// 	     result(0,2*k) +=  (*m_wgrid[v])(i) *
//	      ( -tmp1 *(unormal(0,0)*Gmderiv(0,0)+unormal(1,0)*Gmderiv(1,0) ) -
// 		flux(0,0)* Gderiv(0,0) ) * tmp2;
//
//	    result(0,2*k+1) += (*m_wgrid[v])(i) *
//	      ( -tmp1 *(unormal(0,0)*Gmderiv(2,0)+unormal(1,0)*Gmderiv(3,0) ) -
// 		flux(0,0) * Gderiv(1,0) ) * tmp2 ;
//
//	  }
//                  }


                    // Adaptive step
                    std::stack<iIntervalV> iStack;
                    iStack.push(iIntervalV(breaks[e-1], breaks[e], res_iv(0,0),0));
                    iStack.push(iIntervalV(breaks[e-1], breaks[e], res_iv(1,0),1));
                    gsMatrix<T> lVal, rVal;
                    res_iv.setZero();
                    do
                    {
                        const iIntervalV curr = iStack.top();
                        iStack.pop();

                        const T mid = (curr.left +  curr.right ) /2.0;

                        // Map the Quadrature rule to the element
                        QuRule.mapTo(curr.left, mid, quNodes, quWeights );

                        lVal = derivOnElement( quNodes, quWeights, v, loop_v_eval,
                                               Gmderiv, Gderiv, bd_val, fev);


                        // Map the Quadrature rule to the element
                        QuRule.mapTo(mid, curr.right, quNodes, quWeights );

                        // Evaluate loop values and derivatives
                        loop_v_eval->evaluateAt(quNodes);


                        rVal = derivOnElement( quNodes, quWeights, v, loop_v_eval,
                                               Gmderiv, Gderiv, bd_val, fev);

                        if ( math::abs(lVal.at(curr.index,0) + rVal.at(curr.index,0) - curr.value) > 1e-3  )
                        {
                            iStack.push( iIntervalV(curr.left, mid , lVal(curr.index,0), curr.index ));
                            iStack.push( iIntervalV(mid, curr.right, rVal(curr.index,0), curr.index) );
                        }
                        else
                            res_iv.at(curr.index) += lVal.at(curr.index,0) + rVal.at(curr.index,0);

                       }
                    while ( ! iStack.empty() ); // Adaptive quadrature loop
//*/
//                    result.col(i) += res_iv.col(0);

                    result.template block<1,2>(0,2*i) += res_iv.transpose();

                }

            }
        }

        freeAll(all_loops);

    }

    T distanceL2(const gsFunction<T> & f) const;

    T MRE(const gsFunction<T> & f) const;

    inline T evalOnElement( gsMatrix<T> quNodes,
                            gsVector<T> quWeights,
                            int v,
                           const typename gsGeometry<T>::Evaluator & loop_v_eval,
                           gsMatrix<T> & green_val, 
                           gsMatrix<T> & green_grad, 
                           gsMatrix<T>  & bd_val,
                           gsMatrix<T>  & fluxEv  ) const
    {
        T result(0);
        
        gsVector<T> unormal(2);

        // Evaluate loop values and derivatives
        loop_v_eval->evaluateAt(quNodes);
        
        // Evaluate the boundary function
        m_boundary_fun[v]->eval_into( m_par_rhs ? 
                                      quNodes : loop_v_eval->values(), bd_val );
        
        // Evaluate the Green function
        m_green->eval_into( loop_v_eval->values(), green_val );
        
        // Evaluate the Green normal derivative
        m_green->grad_into( loop_v_eval->values(), green_grad ); 
        
        // Evaluate the flux
        m_flux[v]->eval_into(quNodes, fluxEv);
        
        // For all quadrature points in this interval
        for (index_t k=0; k!= quNodes.cols(); ++k)
        {
            // Compute the outer normal
            loop_v_eval->normal(k, unormal);
            
            const T intElement = quWeights[k] * unormal.norm();
            unormal.normalize();
            
            // Evaluate the Green normal derivative
            const T greenNDeriv = ( green_grad.col(k).transpose() * unormal ).value() ;
            
            result += intElement * 
                ( bd_val(0,k) * greenNDeriv - fluxEv(0,k) * green_val(0,k) );
         
        }// Quadrature loop

        return result;
    }

    inline gsMatrix<T> derivOnElement(gsMatrix<T> quNodes,gsVector<T> quWeights,int v,
                            const typename gsGeometry<T>::Evaluator & loop_v_eval,
                             gsMatrix<T> & Gmderiv,gsMatrix<T> & Gderiv,
                             gsMatrix<T> & bd_val,gsMatrix<T> & fev) const
    {

        gsMatrix<T> result(2,1);
        result.setZero();

        gsVector<T> unormal(2);

        // Evaluate loop values and derivatives
        loop_v_eval->evaluateAt(quNodes);

        // Evaluate the boundary function
        m_boundary_fun[v]->eval_into( m_par_rhs ?
                                      quNodes : loop_v_eval->values(), bd_val );

        // Evaluate mixed derivatives
        m_green->source_mderiv_into( loop_v_eval->values(), Gmderiv );

        // Evaluate source derivatives
        m_green->source_grad_into( loop_v_eval->values(), Gderiv );

        // Evaluate the flux
        m_flux[v]->eval_into(quNodes, fev);

        for (int k=0; k!= quNodes.cols(); ++k) // For all quadrature points in this interval
            {
                       // Compute the outer normal
              loop_v_eval->normal(k, unormal);

               const T intElement = quWeights[k] *unormal.norm();
               unormal.normalize();

               result.at(0,0) += intElement *
                   ( bd_val(0,k) * ( Gmderiv.template block<2,1>(0,k).transpose() * unormal ).value() -
                                     fev(0,k)* Gderiv(0,k) ) ;

               result.at(1,0) += intElement *
                   ( bd_val(0,k) * ( Gmderiv.template block<2,1>(2,k).transpose() * unormal ).value() -
                                     fev(0,k) * Gderiv(1,k) ) ;
               }

return result ;
    }

protected:
    
    struct iInterval
    {
        iInterval(T l, T r, T v)
        : left(l), right(r), value(v)
        { }

        T left;
        T right;
        T value;
    };

    struct iIntervalV
    {
        iIntervalV(T l, T r, T v, int i)
            : left(l), right(r), value(v), index(i)
        { }

        T left;
        T right;
        T value;
        int index;
    };

protected:

    gsGreenFunction2d<T> * m_green;
    //gsCurveLoop<T>    * m_loop; // m_pdomain
    gsPlanarDomain<T> *m_pdomain;
    //gsGeometry<T>     * m_flux;  // std::vector -- same size as planar domain loops
    std::vector<gsGeometry<T>* > m_flux;
    std::vector<gsFunction<T>* > m_boundary_fun; // vector

    std::vector< std::vector<T> > m_breaks;
    
//    std::vector<gsMatrix<T>* >  m_ngrid; // vector
//    std::vector< gsVector<T>* > m_wgrid;
    bool              m_par_rhs;

private:
    // Temporary variables
    
};


/// Function @name harmonicConjugate
/// computes the period of the harmonic conjugate of the harmonic measure of a given curve.
/// Given a planar domain \param domain with inner components C_1, C_2, ... ,C_n (holes)
/// the harmonic measure of every component is defined as the integral over C_j
/// of the normal derivative of the green function \param green.
/// The period is stored in a matrix \param p whose entries are computed as integrals over every
/// C_k of the normal derivative of the harmonic measure.
template<class T>
void harmonicConjugate(gsPlanarDomain<T> const & domain, const gsGreenFunction2d<T> &green, gsMatrix<T> &p)

{

    /* first step compute normal vectors for inner components (holes)
     */
    gsBSpline<T> *out_param;
    std::vector<gsGeometry<T>*> all_loops;
    std::vector<gsBSplineBasis<T> > base;

    std::vector<gsMatrix<T> *> all_ngrid_j, all_ngrid_k, quad_j, quad_k;
    std::vector<gsVector<T>*> all_wgrid_j, all_wgrid_k;

    gsMatrix<T> gr_points, quad1,quad2,u,qu, xev, qxev,unormal_inter, unormal2_inter, gree_mix;
    gsVector<T> unormal(2), unormal2(2);
    T temp1, temp2, costante;
    costante = -1/(2*3.1415926535);

    for(index_t i= 1; i<=domain.numHoles(); i++) //loop over inner components
    {
        out_param=dynamic_cast<gsBSpline<T> * >(domain.loop(i).singleCurve());
        all_loops.push_back( out_param );
        base.push_back( out_param->basis() );
        // range.push_back( out_param.parameterRange() );
        base[i-1].anchors_into(gr_points) ;
        //  all_loops[i]->eval_into( gr_points, y_points ) ;

        std::vector<T> breaks;
        breaks = base[i-1].domain()->breaks();

        //computing quadrature points :

        for ( int v= 1; v<gr_points.size()-1; ++v )
            breaks.push_back( gr_points(0,v) );
        std::sort(breaks.begin(), breaks.end() ) ;
        typename std::vector<T>::iterator itr =
        unique( breaks.begin(), breaks.end(), math::template almostEqual<T> );
        breaks.resize( itr - breaks.begin() ) ;

        gsMatrix<T> * ngrid = new gsMatrix<T>;
        gsVector<T> * wgrid = new gsVector<T>;
        iteratedGaussRule(*ngrid, *wgrid,3*base[i-1].degree(), breaks ) ;
        all_ngrid_j.push_back(ngrid);
        all_wgrid_j.push_back(wgrid);

        base[i-1].eval_into(*ngrid,quad1);
        quad_j.push_back(&quad1);

        iteratedGaussRule(*ngrid, *wgrid,3*base[i-1].degree(), breaks ) ;
        all_ngrid_k.push_back(ngrid);
        all_wgrid_k.push_back(wgrid);

        base[i-1].eval_into(*ngrid,quad2);
        quad_k.push_back(&quad2); //stores quad points over every inner loops
    }
    //computing double integrals
    p.conservativeResize(domain.numHoles(),domain.numHoles());
    p.setZero();

    for(int j=1; j<=domain.numHoles();++j)
    {
        for(int k=1; k<=domain.numHoles();++k)
        {
            //for quad points in j curve
            for(int qj=0; qj!=quad_j[j-1]->cols(); ++qj)
            {   u=all_ngrid_j[j-1]->col(qj);
                all_loops[j-1]->eval_into(u,qu);
                //u=quad[j-1].col(qj);

                all_loops[j-1]->deriv_into(u,unormal_inter);
                temp1=unormal_inter.norm();

                std::swap(unormal_inter(0), unormal_inter(1));
                unormal_inter(0,0) *=-1.0; //inward normal
                unormal_inter.normalize();
                unormal[0]=unormal_inter.at(0,0);
                unormal[1]=unormal_inter.at(1,0);

                //for all quadrature points
                for(int qk=0; qk!=quad_k[k-1]->cols(); ++qk)
                {
                    xev= all_ngrid_k[k-1]->col(qk);
                    all_loops[k-1]->eval_into(xev,qxev);

                    //                             gsInfo<<"\n xev "<<xev(0,0)<<"\n u "<< u(0,0) <<"\n";
                    //                             gsInfo<<"\n xev-u "<< xev(0,0)-u(0,0)<<"\n";
                    if( std::abs(xev(0,0)-u(0,0)) < 0.000001)
                        continue;

                    green.setSourcePoint(qxev);
                    green.source_mderiv_into(qu,gree_mix);
                    //                         gsInfo<<"\n gree_mix : \n"<<gree_mix <<"  for "<< j<<"\n" ;
                    //                          gsInfo<<"\n green source point : \n"<<green.getSource();

                    all_loops[k-1]->deriv_into(xev,unormal2_inter);
                    temp2=unormal2_inter.norm();

                    std::swap( unormal2_inter(0), unormal2_inter(1) );
                    unormal2_inter(0,0) *=-1.0; //inward normal
                    unormal2_inter.normalize();
                    unormal2[0]=unormal2_inter(0,0);
                    unormal2[1]=unormal2_inter(1,0);

                    p(j-1 , k-1) +=all_wgrid_j[j-1]->at(qj)*all_wgrid_k[k-1]->at(qk)
                    * temp1 * temp2*costante *
                    ( unormal2_inter(0,0)*unormal.dot(gree_mix.col(0).segment(0,2))+
                    unormal2_inter(1,0) * unormal.dot(gree_mix.col(0).segment(2,2) ) );


                }
            }
        }
    }


};

template<class T>
T gsBemSolution<T>::distanceL2(const gsFunction<T> & f) const
{
    T result(0);

    const index_t nloops = m_flux.size();

    std::vector<gsGeometry<T>*> all_loops;

    gsMatrix<T> fderiv, fluxEval;

    gsMatrix<T> quNodes  ;
    gsVector<T> quWeights, unormal;
    gsVector<index_t>  numNodes(1);
    gsGaussRule<T> QuRule;
    
    for(index_t v=0; v != nloops; v++)
    {
        // Get a single curve representing every loop
        all_loops.push_back(  m_pdomain->loop(v).singleCurve() );
    }

    // first time: computing integral over (0,1) and evaluate the solution u
    for(index_t v=0; v< nloops; v++)
    {
        numNodes[0] = m_flux[v]->basis().minDegree()+1;
        QuRule.setNodes(numNodes);
        
        //taking breaks points :
        const std::vector<T> & breaks = m_breaks[v];
        
        // Evaluator on the current loop
        typename gsGeometry<T>::Evaluator loop_v_eval (
            all_loops[v]->evaluator(NEED_VALUE   |
                                    NEED_JACOBIAN) );
        
        // Start iteration over sub-elements
        for (std::size_t e = 1; e != breaks.size(); ++e)
        {
            // Map the Quadrature rule to the element
            QuRule.mapTo(breaks[e-1], breaks[e], quNodes, quWeights );
            
            // Evaluate loop values and derivatives
            loop_v_eval->evaluateAt(quNodes);

            f.deriv_into( loop_v_eval->values(), fderiv );
            
            m_flux[v]->eval_into(quNodes, fluxEval );
            
            // For all quadrature points in this interval
            for (int k=0; k!= numNodes[0]; ++k)
            {
                // Compute the outer normal
                loop_v_eval->normal(k, unormal);
                
                const T intElement = quWeights[k] * unormal.norm();
                unormal.normalize();
                
                // Take the difference between approximation and exact value
                const T diff = (fderiv.template block<1,2>(0,2*k) * unormal).value() - fluxEval(0,k);

                result += intElement * diff * diff ;
            }


        }

    }


    freeAll(all_loops);


    return math::sqrt(result);
}

template<class T>
T gsBemSolution<T>::MRE(const gsFunction<T> & f) const

{
    int nloops = m_pdomain->numLoops();
    T sum ;
    std::vector<gsGeometry<T>*>   all_loops ;
    gsMatrix<T> gr_points, col_points, fev, bemev;

   //computation of collocation points: greville abscissae
  for(int v=0; v != nloops; v++)
  {
      all_loops.push_back(  m_pdomain->loop(v).singleCurve() );

      all_loops[v]->basis().anchors_into(gr_points);
      all_loops[v]->eval_into( gr_points, col_points ) ;


      for(int i=0; i<col_points.cols();++i)
      {
          f.eval_into(col_points,fev);
          this->eval_into(col_points,bemev);
          sum +=( ( ( fev.col(i)-bemev.col(i) ).norm() )
                  /( fev.col( i ).norm() ) );

      }
      sum /= ( col_points.cols() );
  }
  freeAll(all_loops);
  gsInfo<<"\n number of collocation points "<< col_points.cols()<<"\n";
  return sum;

}

}; // end namespace gismo
