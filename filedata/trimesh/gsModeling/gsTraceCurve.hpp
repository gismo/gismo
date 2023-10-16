/** @file gsTraceCurve.hpp

    @brief Provides functions for finding the preimage of a curve via a map

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Falini, A. Mantzaflaris
*/

#pragma once

#include <iostream>
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFunction.h>
#include <gsCore/gsFunctionExpr.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsUtils/gsPointGrid.h>


namespace gismo
{

/// \param map is a planar map given by two component functions, that
/// maps onto a template
/// \param p is a 2x1 matrix storing a chosen point in our template
/// \param x is a 2x1 matrix storing a starting point in our computational domain
/// \param result is a 2xn matrix containing points on the pre-image
/// curve
///
/// \ingroup Modeling
template<class T>
void gsTraceLine( std::pair<gsFunction<T>*,gsFunction<T>*>  & map,
                       gsMatrix<T> const & x,gsMatrix<T> const & p, gsMatrix<T> & result)
{
    T h;
    int loops=10;
    int n_points = 10;
    T tolerance = 0.0001;

     gsMatrix<T> line(2,1), px, py;
     gsMatrix<T> next (2,1), delta;
     gsMatrix<T> tangent (2,1) ;
     gsMatrix<T> current (2,1);
     gsMatrix<T> Jmap (2,2);



     map.first->eval_into(x, px);
     map.second->eval_into(x, py);



  // matrix \line stores the evaluation of the \map at the current point
     line.row(0)= px;
     line.row(1)= py;
     next.row(0)= px;
     next.row(1)= py;

     tangent = p.col(0)-line.col(0);

      h = tangent.norm()/n_points ;

     tangent/=tangent.norm();


  // "current" stores the point we start with
     current = x.col(0);

   int count=0;

   for(int i=0; i!=n_points; i++)
   { next = next + h*tangent;
      loops=10;

  do
  {
      count++;
      //gsDebug<< "Point "<<count <<"= "<< current.transpose() <<"\n";

      Jmap.row(0) =  map.first->jacobian(current );
      Jmap.row(1) =  map.second->jacobian(current);

      if ( math::abs( Jmap.determinant() ) < 0.0001) //0.025 amoeba_hole,austria_hole
            {
               gsDebug<< "\n trace line: \n Jacobian vanished at : " << current.transpose() <<"\n"
                   "with Jacobian = "<< math::abs( Jmap.determinant() )<<"\n \n";

               loops =-1 ;
                break; // stop Newton iteration
            }

      delta = Jmap.inverse()*(next.col(0)-line.col(0));
      current += delta;

      map.first->eval_into(current,px);
      map.second->eval_into(current,py);

      line.row(0)=px;
      line.row(1)=py;

      // gsDebug<<"loops="<<loops<<"\n";

  }
  while((--loops > 0) && ( (line-next).norm() > tolerance  ) );


  result = current;

  // gsDebug<<"points traced : "<< result<< "\n";



    // example:
    // gsFunctionExpr<T> g1("x");
    // gsFunctionExpr<T> g2("y");
    // std::pair<gsFunction<T>*,gsFunction<T>*> map(g1,g2);
    // gsMatrix<T> res, p(2,2);
    // p<< 0,1,
    //     0,1 ;
    // gsTraceCurve(map, p, res);



}

}


template<class T>
void gsTraceCurve( std::pair<gsFunction<T>*,gsFunction<T>*>  & map,
                       gsMatrix<T> const & x, gsBSpline<T> * bs,gsMatrix<T> & result,const int n_points = 50,
                       const T tolerance = 0.0001)
{

    int loops=10;// usually are 100


    result.resize(2,n_points-1);

       gsMatrix<T> current (2,1);
       gsMatrix<T> next (2,1), delta;
       gsMatrix<T> line(2,1), px, py;
       gsMatrix<T> Jmap (2,2);

      // to do: change with bs->parameterRange()
      gsMatrix<T> pts =  gsPointGrid<T>(0,1,n_points);

      gsMatrix<T> ev =  bs->eval(pts);

      //next = ev.col(0); // ev.col(0)==0 point removed!
      next= ev.col(1);

      gsMatrix<T> rr; // rr stores the initial point P0 from which the tracing starts
      gsTraceLine( map, x, next, rr );



       current = rr.col(0);


       result.col(n_points-2)= current;


         next=ev.col(1);

     map.first->eval_into(current, px);
     map.second->eval_into(current, py);

  // matrix \line stores the evaluation of the \map at the current point
     line.row(0)= px;
     line.row(1)= py;

   int count=0;

   for(int i=2; i!=n_points; i++)
   {
      next = ev.col(i);
   loops=10;
  do
  {

  count++;


  // TO DO: transpose the result of  gsFunction::grad and add transpose()
  Jmap.row(0) =  map.first->jacobian(current);
  Jmap.row(1) =  map.second->jacobian(current);


   delta = Jmap.inverse()*(next.col(0)-line.col(0));

   current += delta;

   //   gsMatrix<T> * grad1 = map.first->deriv() ;
   //   Jmap.row(0) = *grad1;
   //   delete grad1;

   map.first->eval_into(current,px);
   map.second->eval_into(current,py);

   line.row(0)=px;
   line.row(1)=py;

   //gsDebug<<"\n loops = "<<loops<<"\n";


  }
  while((--loops > 0) && ( (line-next).norm() > tolerance  ) );


  result.col(n_points-i-1) = current;
// gsDebug<<"\n result.col("<<i<<") : \n" << result.col(i)<<"\n";




    // example:
    // gsFunctionExpr<T> g1("x");
    // gsFunctionExpr<T> g2("y");
    // std::pair<gsFunction<T>*,gsFunction<T>*> map(g1,g2);
    // gsMatrix<T> res, p(2,2);
    // p<< 0,1,
    //     0,1 ;
    // gsTraceCurve(map, p, res);



}


}

/// \param map is a planar map given by two component functions,
/// \param x: point in the planar domain, should correspond to pre image (via \em map) to the
/// middle point of \em bs.
/// \param bs is the curve you want to trace
/// \param t0 and...
/// \param t1 parametric values, they represent the interval you want to discretize
/// \param result stores resulting points
/// \param n_points is the number of points you want to trace per curve
/// \param tolerance
template<class T>
void gsTraceCurvePart(std::pair<gsFunction<T>*,gsFunction<T>*>  & map,
                      gsMatrix<T> const & x, gsBSpline<T> * bs, T t0, T t1,
                      gsMatrix<T> & result,const int n_points = 50,
                      const T tolerance = 0.0001)
{
    int loops=100;

    gsMatrix<T> current (2,1);
    current=x;
    gsMatrix<T> next (2,1), delta, rr;
    gsMatrix<T> line(2,1), px, py;
    gsMatrix<T> Jmap (2,2);
    gsMatrix<T> medium_point (1,1);
    medium_point<< 0.5 ;

    gsMatrix<T> discretI = gsPointGrid<T>(t0,t1, n_points );
    gsMatrix<T> ev = bs->eval(discretI);
    map.first->eval_into(current, px);
    map.second->eval_into(current, py);
    line.row(0)= px;
    line.row(1)= py;

    result.resize(2,ev.cols()+1);
    result.col(0)=x;
    bool cont = true;
    index_t i;
    for(i =1; i< ev.cols()&& cont; i++)
    {
        next= ev.col(i);
        loops=100;
        do
        {
            // TO DO: transpose the result of  gsFunction::grad and add transpose()
            Jmap.row(0) =  map.first->jacobian(current);
            Jmap.row(1) =  map.second->jacobian(current);

            //gsDebug<<"jacobian matrix is \n"<< Jmap<<"\n";
            if ( math::abs( Jmap.determinant() ) < 0.01) //0.025 amoeba_hole,austria_hole
            {
                gsDebug<< "Jacobian vanished at : " << current.transpose() <<"\n"
                "with Jacobian = "<< math::abs( Jmap.determinant() )<<"\n \n";

                cont = false;
                break; // stop Newton iteration
            }

            delta = Jmap.inverse()*(next.col(0)-line.col(0));

            current += delta;

            map.first->eval_into(current,px);
            map.second->eval_into(current,py);

            line.row(0)=px;
            line.row(1)=py;



        }
        while((--loops > 0) && ( (line-next).norm() > tolerance  ) );

        result.col(i) = current;
    }

  // gsDebug<< "Size expected:"<< result.cols() <<"\n";
    result.conservativeResize( Eigen::NoChange, i);
    //gsDebug<< "Size finally :"<< result.cols()  <<"\n";





}

} // namespace gismo
