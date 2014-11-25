/** @file gsPlanarDomain.hpp

    @brief Provides implementation  of the PlanarDomain class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Falini, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsConstantFunction.h>
#include <gsCore/gsCurve.h>
#include <gsSolver/gsBemUtils.h>

#include <gsNurbs/gsNurbsCreator.h>

#include <gsUtils/gsMesh/gsMesh.h>

#include <gsSolver/gsBemLaplace.h>
#include <gsSolver/gsBemUtils.h>

#include <gsUtils/gsQuadrature.h>
#include <gsModeling/gsTraceCurve.hpp>
#include <gsModeling/gsTemplate.h>

#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsKnotVector.h>

namespace gismo
{


template <class T>
gsPlanarDomain<T>::gsPlanarDomain( std::vector< gsCurveLoop<T> *> const & loops)
{
    if(!loops[0]->is_ccw())
        loops[0]->reverse();
    m_loops.push_back(loops[0]);
    for(std::size_t i=1; i<loops.size();i++ )
    {

        if( loops[i]->is_ccw() )
            loops[i]->reverse();
        m_loops.push_back(loops[i]);
    }
    updateBoundingBox();
}


template<class T>
void gsPlanarDomain<T>::segment( int n_points, T tolerance,bool circle, gsMesh<T> &segmentation)
 // gsMatrix<T> const & start, int n_points,  tolerance
{
    // ***************************   Consider a template
gsTemplate<T> * pmy_template;
    if(circle)
    pmy_template = new gsTemplate<T>(this->numHoles() ); // Convex template with no holes

    else
    {
       bool sq=true;
       pmy_template= new gsTemplate<T>(sq,1);
    }

    gsTemplate<T>& my_template (*pmy_template);

     // ***************************   first step: solve Laplace problems


    // call : std::pair< gsFunction<T> *, gsFunction<T> *, > S = this->mapto( template );
    // u.first, u.second are the components of the solution
    std::pair <gsFunction<T>* , gsFunction<T>*> u;

    u=this->mapto(my_template, 2); // 2 for amoeba_hole, 3 for austria_hole
                                   // 4 for Puzzle1 domain


    // construct indice points: points inside the planar domain at a suitable distance from
    // the outer boundary

    gsMatrix<T> tab(2,2);
    tab= *safe(this->boundingBox());
    gsMatrix<T> storage;
    storage= uniformPointGrid<T>( tab.col(0), tab.col(1),600);
    gsMatrix<T> indice (2,storage.cols()); // matrix indice stores points of the grid which are
                                           // inside the planar domain
    indice.setZero();
    int k=0;

    gsBSpline<T>  b =*safe(dynamic_cast<gsBSpline<T> * > ( this->outer().singleCurve() ));

    gsMatrix<T> mod = b.coefs();
    mod *=0.80;//0.85
    gsKnotVector<T> kv = b.basis().knots();

    gsBSpline<T> *interior = new gsBSpline<T>(kv, mod);


    gsCurveLoop<T> * inner = new gsCurveLoop<T>( interior );
    gsPlanarDomain<T> Idomain( inner);
    
    if(this->numHoles()!=0)
    {
       gsMatrix<T> CI (9,2);
       CI <<-0.7 ,0,
            -0.7 ,0.3,
            -1, 0.3,
            -1.3 ,0.3,
            -1.3 ,0 ,
            -1.3, -0.3,
            -1 ,-0.3 ,
            -0.7 ,-0.3,
            -0.7, 0;
       gsKnotVector<T> KV2 (0,1,3,3,2) ;
       gsBSpline<T> *holeI = new gsBSpline<T>(KV2,CI);
       holeI->reverse();
        //    gsBSpline<T> *holeI=gsNurbsCreator<>::BSplineFatCircle();
        //holeI->reverse();
    Idomain.insertHole(new gsCurveLoop<T>( holeI ));
    }


    index_t j;
    for(j=0;j<storage.cols();j++)
    {
        if(  Idomain.inDomain(storage.col(j),1) && Idomain.inDomain( storage.col(j),0 ) ) //point inside the domain
        {
            indice.col(k).noalias() = storage.col(j);
            // cout<<"indice (0,k) : " << indice(0,k) << "\n k is : "<< k<<"\n";
            k++;

        }
    }
    indice.conservativeResize(2,k);

//******************************* JAKA EXAMPLES
   /* new strategy for non convex sets */
 //  indice.conservativeResize(2,1);
  /* choose a point inside the planar domain : i.e. (2.5,2)for Puzzle1Scaled times 10, puzzle1 scaled times 5 is (1.5,1)
   and store it in indice.  (point p (0.25, 0.2 for Puzzle1)) // 0.8,0.5; //Curtain (1.2,0.0), (0.0, 0.2) Puzzle2*/
 //  indice<<0.4, 0.1;
//******************************************************************
    gsMatrix<T> Image, Image1, Image2;
    u.first->eval_into(indice, Image1);
    u.second->eval_into(indice, Image2);
//    gsInfo<<"\n x coordinates of points mapped : \n"<<Image1<<"\n";
//    gsInfo<<"\n y coordinates of points mapped : \n"<<Image2<<"\n";
    Image.resize(2, Image1.cols() );
    Image.row(0) = Image1;
    Image.row(1) = Image2;

    gsMatrix<T> medium_point (1,1), middle, x, leftPart, rightPart;
    medium_point<< 0.5 ;
    // ***************************     second step: trace boundary curves, starting from a point
    // ***************************     in the middle, then applying the algorithm
    // ***************************     in both directions

    gsMatrix<T> param;

    for(int s=0;s!=(my_template.skeletonSize());++s)
    {

        middle=my_template.skeleton(s)->eval(medium_point);
        param = my_template.skeleton(s)->parameterRange();
        j=0;
        T dist = (Image.col(j) - middle ).norm();
        for (index_t i=1; i!= Image.cols(); ++i)
        {
            if( ( Image.col(i) - middle ).norm() < dist )
                j = i;
        }
/************ For JAKA EXAMPLE
//        T endp1, endp2;
//        gsMatrix<T> divid (1,1);
//        divid<< 1.0/T(3*n_points);
//        endp1= param(0,0) + divid(0,0);
//        endp2 = param(1,0) - divid(0,0);
***********************************************/

        gsTraceLine<T>(u, indice.col(j), middle, x );
        gsTraceCurvePart<T>(u,x,my_template.skeleton(s), param(0,1)/2, param(0,0), leftPart, n_points, tolerance );
        leftPart = leftPart.rowwise().reverse().eval();

        gsTraceCurvePart<T>(u,x,my_template.skeleton(s), param(0,1)/2, param(0,1), rightPart, n_points, tolerance );
        leftPart.conservativeResize( 2, leftPart.cols()+ rightPart.cols() );
        leftPart.rightCols( rightPart.cols() ) = rightPart;

        segmentation.addLine(leftPart);
    }


    // third step: construct Coon's patches


    //Clean Memory
    delete u.first;
    delete u.second;
    delete pmy_template;
}


template<class T>
void gsPlanarDomain<T>::getLamdas(gsBemSolution<T> & f, gsVector<T> &lambdas)
{
    std::vector<gsCurve<T>*> all_loops ;
    std::vector<gsBSplineBasis<T> > base;
    std::vector<gsMatrix<T> *> all_ngrid, quad;
    std::vector<gsVector<T>*> all_wgrid;
    gsMatrix<T> unormal_inter, gr_points, quad1, u, xev,val, gamma,p;
    gsVector<T> aux,unormal;
    gsBSpline<T> *component = NULL;
    T temp1, aiuto1,aiuto2, help1, help2;   
    for(int i=1; i<=this->numHoles(); ++i)
    {
        component = dynamic_cast<gsBSpline<T> * >(this->loop(i).singleCurve() );
        all_loops.push_back(component);
             
         base.push_back( component->basis() );

         base[i-1].anchors_into(gr_points) ;
// gsInfo<<"\n in pdomain checking : "<<"\n";
// this->check();
        
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
        iteratedGaussRule(*ngrid, *wgrid,3* base[i-1].degree(), breaks ) ;
        all_ngrid.push_back(ngrid);
        all_wgrid.push_back(wgrid);
        
        base[i-1].eval_into(*ngrid,quad1);
        quad.push_back(&quad1);
        
    }
    aux.resize(2);
    aux.setZero();
    unormal.resize(2);
    unormal.setZero();
    lambdas.resize(this->numHoles());
    lambdas.setZero();
    gamma.resize(this->numHoles(),1);
    gamma.setZero();
    
    for(int k=1; k<=this->numHoles(); ++k)
    {
        for(int q=0; q!=quad[k-1]->cols();++q)
        {
           u= all_ngrid[k-1]->col(q); 
           all_loops[k-1]->eval_into(u,xev);
           all_loops[k-1]->deriv_into(u,unormal_inter);
           temp1= unormal_inter.norm();
           
           std::swap(unormal_inter(0), unormal_inter(1) );
           unormal_inter(1,0) *= -1.0;
           unormal_inter.normalize();
           help1 = unormal_inter.at(0,0);
           help2 = unormal_inter.at(1,0);
           unormal[0]= help1;
           unormal[1]= help2;
           
           f.deriv_into(xev, val);
        
           aiuto1 = val(0,0);
           aux[0] = aiuto1;
           aiuto2 = val(0,1);
           aux[1] = aiuto2;
         
           gamma(k-1, 0) += all_wgrid[k-1]->at(q)*temp1*(unormal.dot(aux)) ;
        }
        
    }
    gamma *= -1.0; 
    
    const gsGreenFunction2d<T> & Gr = f.getGreen();
    harmonicConjugate(*this,Gr,p);
    
    lambdas = p.transpose().colPivHouseholderQr().solve(gamma);
    gsInfo<<"\n matrix p : \n "<< p <<"\n gamma :\n "<<gamma<<"\n lambadas : "<<lambdas<<"\n";

 delete component;
}

template<class T>
std::pair<gsFunction<T>* , gsFunction<T>* >
gsPlanarDomain<T>::mapto(gsTemplate<T> & my_template, int numRefine)
{
  //  GISMO_ASSERT( my_template.domain().numLoops()== this->numLoops(), "Template has different topology.");

    //******* basis vector of all basis functions of all curves in pd
    std::vector < gsBasis<T> * > basis;

    std::vector<gsFunction<T>*> bc_0;
    std::vector<gsFunction<T>*> bc_1;
    gsMatrix<T> param(1,1);
    param<<1;
    gsBSpline<T> * tmpl_loop = dynamic_cast<gsBSpline<T> * >( my_template.loop(0).singleCurve());

    for(index_t v=0; v<=this->numHoles(); v++)
    {
        gsBSpline<T> * this_loop = dynamic_cast<gsBSpline<T> * >( this->loop(v).singleCurve() );

        basis.push_back( this_loop->basis().clone() );
        if(v==0)
        {
          /* for the square we need :
           *  gsKnotVector<T>sqKnots(0,1,3,2,1,1);
           *  in the place of
           *  tmpl_loop->basis().knots()*/
        bc_0.push_back( new gsBSpline<T>( tmpl_loop->basis().knots(),  tmpl_loop->coefs().col(0))  );
        bc_1.push_back( new gsBSpline<T>( tmpl_loop->basis().knots(),  tmpl_loop->coefs().col(1))  );
        }
        else
        {

            bc_0.push_back( new gsFunctionExpr<T> ("0.167231") );//-0.334477 CII
            bc_1.push_back( new gsFunctionExpr<T> ("0.0") );
            
/***************************************************************
 * please, do not delete the commented part             
//             if(this->numLoops()==1)
//             { */
//             std::vector<gsFunction<T>*> bc_map(2);
//             bc_map[0]=bc_0[0];
//             bc_map[1]=bc_1[0];
//            gsMatrix<T> av = averageValue(bc_map, tmpl_loop->basis().knots().breaks() );
// //             gsMatrix<T> av (2,1);
// //             av<< 0.5,
// //                  0;
//            bc_0.push_back( new gsConstantFunction<T> (av(0,0), 1 ) );
//             bc_1.push_back( new gsConstantFunction<T> (av(1,0), 1 ) ); 
//             }
//             else
//             {
//                 bc_0.push_back( new gsConstantFunction<T> (*my_template.skeleton(0)->eval(param), 1 ) );
//                 bc_1.push_back( new gsConstantFunction<T> (*my_template.skeleton(1)->eval(param), 1 ) ); 
//             }
//         }


//             std::vector<gsFunction<T>*> bc_map(2);
//             bc_map[0]=bc_0[0];
//             bc_map[1]=bc_1[0];
//             gsMatrix<T> av = averageValue(bc_map, tmpl_loop->basis().knots().breaks() );
//             bc_0.push_back( new gsConstantFunction<T> (av(0,0), 1 ) );
//             bc_1.push_back( new gsConstantFunction<T> (av(1,0), 1 )  );
        
//***************************************/
       }
        delete this_loop;
    }

    
    delete tmpl_loop;

    for(index_t i=0; i<this->numLoops() ; i++)
    {
        for (int k = 0; k < numRefine; ++k)
            basis[i]->uniformRefine();
    }

    int count = 0;
    for(index_t i=0; i< this->numLoops() ; i++)
        count += basis[i]->size();

    gsInfo<<"Number of DoFs: "<< count <<"\n";

    //  gsBemLaplace<T> Lsolver_1( m_loops[0], f_1 , basis );
    gsBemLaplace<T> Lsolver_1( this, bc_0, basis );

    // using f_2 we can get u_2
    // gsBemLaplace<T> Lsolver_2( m_loops[0], f_2 , basis );
    gsBemLaplace<T> Lsolver_2( this, bc_1, basis );

    //////////////////********************
    gsBemSolution<T> * sol1= Lsolver_1.solve(true);
    gsBemSolution<T> * sol2= Lsolver_2.solve(true);
    ////////////////////*******************

    ///********************
   

    // as lambdas get image of centers of holes
    // from mpdomain getting centers of holes : TO DO, so far : by hand 
    // using sol1 and sol2 for getting coordinates of point and initializationg lambdas vector
    // construct again a gsBemSolution_forHoles object and using eval_into function . 
    

/*
     gsVector<T> lamb, mis;
     if(this->numHoles()!=0)
     {
         std::vector<gsFunction<T> *> bound1;
         std::vector<gsFunction<T> *>  bound2;
         std::vector<gsMatrix<T> *>  grid1;
         std::vector<gsMatrix<T> *>  grid2;
         std::vector<gsVector<T> *>  weight1;
         std::vector<gsVector<T> *>  weight2;
         std::vector<gsGeometry<T> *> flusso1;
         std::vector<gsGeometry<T> *> flusso2;
         for(unsigned i=0;i!=sol1->getBoundary_fun().size();++i)
         {
            bound1.push_back(sol1->getBoundary_fun()[i]);
            bound2.push_back(sol2->getBoundary_fun()[i]);
         }
         for(unsigned i=0; i!= sol1->flux().size();++i)
         {     
            flusso1.push_back(sol1->flux()[i]);
            flusso2.push_back(sol2->flux()[i]);
         }
         for(unsigned i=0;i!=sol1->getNgrid().size();++i)
         {
             grid1.push_back(sol1->getNgrid()[i]);
             grid2.push_back(sol2->getNgrid()[i]);
         }
          for(unsigned i=0;i!=sol1->getWgrid().size();++i)
          {
              weight1.push_back(sol1->getWgrid()[i]);
              weight2.push_back(sol2->getWgrid()[i]);
          }

              
////          gsInfo<<"\n sol1->getBoundary_fun() :"<<sol1->getBoundary_fun()[0]<<"\n";
         
////       gsInfo<<"\n matrix p is : "<< p <<"\n";


//     // getLamdas(*sol1,lamb);
//     //getLambdas function replaced by new strategy of mapping centers by sol1 (sol2)
//     // and taking images
//     gsMatrix<T> cent( 2,this->numHoles() ), Imag_cent ;
     
//     for(int s=1;s<=this->numHoles();++s)
//     {
//         gsBSpline<T> * inner_component = dynamic_cast<gsBSpline<T> * >( this->loop(s).singleCurve() );
//        gsMatrix<T> C= inner_component->coefs();
//        cent.col(s-1)=( ( C.row(0)+C.row(4) )/T(2)).transpose();
     
//    }
    
//       sol1->eval_into(cent,Imag_cent);
//       for(int i=0; i< Imag_cent.cols();++i)
//           lamb[i]= Imag_cent(0,i);
       
//      gsBemSolution_forHoles<T> * first_component = new gsBemSolution_forHoles<T>(this,flusso1,
//                                                    bound1, lamb, grid1, weight1, true);
    
     
//   //   getLamdas(*sol2,mis);
//      sol2->eval_into(cent,Imag_cent);
//      for(int i=0;i<Imag_cent.cols();++i)
//          mis[i]= Imag_cent(0,i);
      
//      gsBemSolution_forHoles<T> * second_component = new gsBemSolution_forHoles<T>(this, flusso2,
//                                                     bound2, mis, grid2, weight2, true);

//      // Free the bases
//      freeAll(basis);
     
//      return std:: make_pair(first_component, second_component);
//    }
//*/
   // else{
        // Free the bases
        freeAll(basis);
        return std::make_pair(sol1, sol2);    
     //   }

}

template<class T>
gsMatrix<T> gsPlanarDomain<T>::averageValue( std::vector<gsFunction<T>*> const &f,
                                             std::vector<T> const & breaks)
{
    gsMatrix<T> ngrid;
    gsVector<T> wgrid;
    gsMatrix<T> xev,B (f.size(),1);
    B.setZero();
    iteratedGaussRule(ngrid, wgrid, 2, breaks);
    for(std::size_t i=0; i<f.size();i++)
    {
        for(int k=0; k<ngrid.cols();k++) //for all quadrature points
        {
            const gsMatrix<T> & uc = ngrid.col(k);

            f[i]->eval_into(uc,xev);
            B(i,0)  += wgrid[k]*xev(0,0);
        }
    }
    gsDebugVar( B.transpose() );
    return B;
}




template<class T>
bool gsPlanarDomain<T>::inDomain( gsMatrix<T> const & u, int direction)
{
    GISMO_ASSERT(u.cols() == 1, "Waiting for a single point");

    // compute intersections with line x=u(0,0);
    std::vector<T> tmp = m_loops[0]->lineIntersections(direction, u(direction,0) );

    //   gsDebug<< "intersections:\n";
    //   for( std::size_t i = 0; i!=tmp.size(); ++i)
    //     gsInfo<< " "<< tmp[i];
    //   gsInfo<<"\n";

    if ( tmp.empty() ) // point outside the outer loop
        return false;

    gsAsMatrix<T> xx(tmp);
    gsMatrix<T> e;
    m_loops[0]->curve(0).eval_into( xx, e );

    int count = 0; // count even means out of the domain

    //checking if abscissae of intersection are greater then point's abscissa
    for(index_t i=0; i!=e.cols(); i++)
    {
        //gsInfo<<" matrix e (1,i): "<< e(1,i) <<"\n";
        if( (e( !direction, i )> u(!direction,0))  )
            count++; //gsInfo<< "count "<< count <<"\n";
    }
    if ( (count % 2) == 0 ) //this means point is outside the outer loop
        return false;

    for(index_t v = 1; v<this->numLoops(); v++)
    {
        //gsInfo<<"m_loops [v] : "<< m_loops[v] <<"\n";
        //gsInfo<<"m_loops[v]->lineIntersections(0,u(0,0)) "<<m_loops[v]->lineIntersections(0,u(0,0)) <<"\n";
        count=0;
        tmp= m_loops[v]->lineIntersections( direction, u(direction,0) );

        //gsInfo<<"tmp size is : "<< tmp.size() <<"\n";
        if(tmp.size()!=0) // intersections detected with hole loop
        {
            gsAsMatrix<T> Ev(tmp);
            m_loops[v]->curve(0).eval_into( Ev, e );

            for(index_t i=0; i!=e.cols(); i++)
            {
                // gsInfo<<"e(1,i) "<< e(1,i)<<"\n";
                if( e(!direction,i) > u(!direction,0) )
                    count++;
                //gsInfo<<"count_holes : "<<count_holes<<"\n";
            }
        }
        if( (count % 2) != 0 ) //this means point is inside hole v
            return false;
    }

    return  true; // if we get here then the point is in the domain
}


template<class T>
bool gsPlanarDomain<T>::onBoundary(gsMatrix<T> const & u)
{
    for(index_t v=0; v< this->numLoops();v++)
    {
        // true iff point u on boundary
        T parValue;
        if( m_loops[v]->isOn(u, parValue,1e-5 ) )
        { //v=this->numLoops();
            return true;
            //gsInfo<<"\n in loop "<<v<< "\n";
        }
    }
    return false;
}


template <class T>
std::ostream& gsPlanarDomain<T>::print(std::ostream &os) const
{
    os << "Outer boundary: " << *m_loops[0];
    if ( m_loops.size()>1 )
    {
        os << "Holes: ";
        for ( typename std::vector< gsCurveLoop<T> *>::const_iterator it =
                m_loops.begin()+1;  it != m_loops.end(); ++it)
            os << **it;
    }
    os << "Bounding box: ["<< m_bbox.col(0).transpose()<< ", "
        << m_bbox.col(1).transpose() << "]";
    return os;
}


/// linearly discriti
template <class T>
void gsPlanarDomain<T>::sampleLoop_into( int loopID, int npoints, int numEndPoints, gsMatrix<T> & u )
{
    assert( (loopID>=0) && (loopID < numLoops()) );

    int np; // new number of points
    switch (numEndPoints)
    {
    case (0):
        np = npoints-2;
        break;
    case (1):
        np = npoints-1;
        break;
    case (2):
        np = npoints;
        break;
    default:
        np = 0;
        break;
    }

    u.resize(2, (m_loops[loopID]->curves()).size() * np);
    int i=0;
    typename std::vector< gsCurve<T> *>::const_iterator it;
    gsMatrix<T> pts;
    gsMatrix<T> uCols;

    int firstInd=0;
    int secondInd=np-1;
    if (numEndPoints==0) {firstInd=1;secondInd=npoints-2;}

    for ( it= (m_loops[loopID]->curves()).begin(); it!= (m_loops[loopID]->curves()).end(); ++it )
    {
        //gsMatrix<T> * interval = (*it)->parameterRange();
        //gsMatrix<T> *  pts = gsPointGrid( interval->at(0), interval->at(1), npoints );
        pts.resize(1,np);
        for (int ii=firstInd;ii<=secondInd;ii++) pts(0,ii-firstInd)= T(ii)/(npoints-1);
        uCols.resize(2,np);
        (*it)->eval_into( pts, uCols );
        u.middleCols( i * np,np ) = uCols;

        (*it)->eval_into( pts, uCols );
        u.middleCols( i * np,np ) = uCols;
        ++i;
    }
}
template <class T>
T getDistance(gsVertex<T>* v1,gsVertex<T>* v2)  // todo: move as member of gsVertex 
{
    T x1=(v1->coords.x());
    T x2=(v2->coords.x());
    T y1=(v1->coords.y());
    T y2=(v2->coords.y());
    T z1=(v1->coords.z());
    T z2=(v2->coords.z());
    T dist=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
    return  dist;
}

template <class T>
void gsPlanarDomain<T>::sampleCurve_into( int loopID, int curveID, int npoints, gsMatrix<T> & u )
{
    assert( (loopID>=0) && (loopID < numLoops()) );
    assert( (curveID>=0) && (curveID < loop(loopID).size() ) );
    u.resize(2,npoints);
    typename std::vector< gsCurve<T> *>::const_iterator it;
    for ( it= (m_loops[loopID]->curves()).begin()+curveID; it!= (m_loops[loopID]->curves()).begin()+curveID+1 ; ++it )
    {
        //gsMatrix<T> * interval = (*it)->parameterRange();
        //gsMatrix<T> *  pts = gsPointGrid( interval->at(0), interval->at(1), npoints );
        gsMatrix<T> pts(1,npoints);
        for (int ii=0;ii<=npoints-1;ii++)
        { pts(0,ii)= T(ii)/(npoints-1); }; // todo: use gsPointGrid
        (*it)->eval_into( pts, u );
    }
}

/// Return a triangulation of the planar domain
template <class T>
gsMesh<T> * gsPlanarDomain<T>::toMesh(int npoints) const     // FOR NOW ONLY ONE LOOP
{
    gsMesh<T> * m = new gsMesh<T>();
    // Typedefs
    typedef typename gsVertex<T>::gsVertexHandle VertexHandle;
    typedef std::deque<VertexHandle>            VertexList;
    typedef typename std::deque<VertexHandle>::iterator  VertexListIter;

    GISMO_ASSERT(npoints > 0, "Number of points must be positive.");

#if FALSE
    T h_x = (m_bbox(0,1) - m_bbox(0,0) ) / (npoints-1) ,
      h_y = (m_bbox(1,1) - m_bbox(1,0) ) / (npoints-1) ;

    //========== 1. Compute y-parallel segments
    VertexList x_seg_begin; // Starting points of y-parallel segments
    VertexList x_seg_end;   // Endpoints  of y-parallel segments

    T xi = m_bbox(0,1);
    std::vector<T> tmp = m_loops[0]->lineIntersections(0, xi );// segments
    std::cout<<"----- INTs:  "<<tmp.size()<<"( "<<xi<<") \n";
    gsAsMatrix<T> xxi(tmp);
    gsMatrix<T> xev;
    m_loops[0]->curve(0).eval_into(xxi, xev );// Push forward onto the curve: xev.row(1)==xi
    gsVector<T> x_seg = xev.row(1) ;// xev.row(0)==xi
    std::sort( x_seg.data(), x_seg.data()+x_seg.size() );

    for ( index_t j= 0; j!= x_seg.size(); j+=2 )
    {
        VertexHandle v = m->addVertex( xi, x_seg[j] );
        if ( x_seg[j] < x_seg[j+1] )
        {
            x_seg_begin.push_back( v );
            v = m->addVertex( xi, x_seg[j+1] );
            x_seg_end.push_back ( v );
        }
        else
        {
            x_seg_end.push_back( v );
            v = m->addVertex( xi, x_seg[j+1] );
            x_seg_begin.push_back ( v );
        }
    }

    xi = m_bbox(0,0);
    for ( int i = 0; i!= npoints-1; ++i )
    {
        tmp = m_loops[0]->lineIntersections(0, xi );
        std::cout<<"----- INTs:  "<<tmp.size()<<"( "<<xi<<") \n";
        //std::cout<<"----- INTs x:  "<<tmp.size()<<", i="<<i<<"\n";
        gsAsMatrix<T> xxi(tmp);
        m_loops[0]->curve(0).eval_into(xxi, xev );
        x_seg = xev.row(1) ;// xev.row(0)==xi
        std::sort( x_seg.data(), x_seg.data()+x_seg.size() );

        for ( index_t j= 0; j!= x_seg.size(); j+=2 )
        {
            VertexHandle v = m->addVertex( xi, x_seg[j] );
            if ( x_seg[j] < x_seg[j+1] )
            {
                x_seg_begin.push_back( v );
                v = m->addVertex( xi, x_seg[j+1] );
                x_seg_end.push_back ( v );
            }
            else
            {
                x_seg_end.push_back( v );
                v = m->addVertex( xi, x_seg[j+1] );
                x_seg_begin.push_back ( v );
            }
        }
        xi += h_x;
    }

    //========== 2. Sort the segments wrt both endpoints
    std::sort( x_seg_begin.begin(), x_seg_begin.end(), Yless<T> );
    std::sort( x_seg_end.begin()  , x_seg_end.end()  , Yless<T> );

    std::cout<<"x-segments: "<<x_seg_begin.size()<<", "<<x_seg_end.size()<<"\n";

    //========== 3. March on x-parallel segments
    bool SegStart(false), SegEnd(false);
    VertexList line0; // x-sorted

    VertexListIter ss = x_seg_begin.begin(),
                   es = x_seg_end.begin();
    T yi = m_bbox(1,0);
    VertexList Yprev;

    for ( int i = 0; i!= npoints; ++i )
    {
        tmp = m_loops[0]->lineIntersections(1, yi );
        //std::cout<<"----- INTs y: "<<tmp.size()<<"\n";

        if ( ! tmp.empty() ) // Intersection event
        {
            std::cout<<"---  Intersection event "<<i<<"( "<<yi<<" )\n";

            // Compute intersections with y=yi
            gsAsMatrix<T> yyi(tmp);
            gsMatrix<T> yev;
            m_loops[0]->curve(0).eval_into(yyi, yev );
            gsVector<T> y_seg = yev.row(0) ;// yev.row(1)==yi
            std::sort( y_seg.data(), y_seg.data()+y_seg.size() );
            VertexHandle v;
            VertexList line1;

            // Look for SegEnd events
            while( es != x_seg_end.end() &&  (*es)->y() < yi )
            {// Property: every line has unique points wrt x-coord
                std::cout<<"SegEnd: "<< **es ;
                SegEnd=true;

                for ( VertexListIter it= line0.begin(); it!= line0.end(); ++it )
                    if (  (*it)->x() == (*es)->x() )
                    {
                        std::cout<<"SegEnd: removed "<< **it ;
                        line0.erase( it );
                        break;
                    }
                es++;
            }

            // Mirror previous lines
            for ( VertexListIter it= line0.begin(); it!= line0.end(); ++it )
            {
                v = m->addVertex( (*it)->x(),  yi );
                line1.push_back(v);
            }

            // Look for SegStart events
            while( ss != x_seg_begin.end() &&  (*ss)->y() < yi )
            {
                SegStart=true;
                std::cout<<"SegStart: added "<< (*ss)->x() <<".\n";
                line0.push_back(*(ss) );
                v = m->addVertex( (*ss)->x(),  yi );
                line1.push_back(v);
                ss++;
            }

            std::sort( line0.begin(), line0.end(), Xless<T> );
            std::sort( line1.begin(), line1.end(), Xless<T> );

            std::cout<<"line0  has "<< line0.size()  <<" points\n";
            std::cout<<"line1  has "<< line1.size()  <<" points\n";

            // add faces
            if ( ! line1.empty() && line1.size() == line0.size() )
            {
                std::cout<<"Tiling row "<<i<<".\n";
                VertexListIter it0=line0.begin();
                for ( VertexListIter it1= line1.begin(); it1!= line1.end()-1; ++it1, ++it0 )
                    m->addFace( *it0, *it1,  *(it1+1), *(it0+1) );
            }
            else
            {
                std::cout<<"Trouble..\n";
                VertexListIter it0=line0.begin();
                for ( VertexListIter it1= line1.begin(); it1!= line1.end(); ++it1, ++it0 )
                    std::cout<< (*it0)->x() <<" - "<< (*it1)->x()  <<" \n";
                while ( it0 != line0.end() )
                    std::cout<< (*it0++)->x() <<" -     * " <<" \n";

            }

            std::cout<<"Connecting endpoints.\n";

            VertexList Yseg;
            // Make boundary vertices on y=yi
            Yseg.push_back( m->addVertex( y_seg[0]             , yi ) );
            Yseg.push_back( m->addVertex( y_seg[y_seg.size()-1], yi ) );

            if (SegStart )
            {
                // connect x-boundary point to level 0
                if ( Yprev.size() )
                {
                    m->addFace( Yprev[0],  line0.front(), line0[1] );
                    m->addFace( Yprev[1],  line0.back() , line0[line0.size()-2] );
                }
                // connect x-boundary point to level 1
                m->addFace( Yseg[0],  line1.front(), line0.front() );
                m->addFace( Yseg[1],  line1.back(), line0.back() );
            }
            else if ( Yprev.size() )
            {
                // Connect y-boundary points
                m->addFace( line1.front(),  line0.front(), Yprev[0], Yseg[0] );
                m->addFace( line1.back() ,  line0.back() , Yprev[1], Yseg[1] );
            }

            Yprev = Yseg;
            line0 = line1;
            SegStart=SegEnd=false;
        }
        else // No intersection event
        {
            Yprev.clear();
            line0.clear();
        }

        // next y-line
        yi += h_y;
    }

    // repeat for yi=m_bbox(1,1)

    return m;
#endif
    T bbox_length_y =m_bbox(1,1)-m_bbox(1,0);
    int yPoints = static_cast<int>( npoints*bbox_length_y) ;
    int lb_yPoints=25;
    if(yPoints<lb_yPoints)
    {
        if ( yPoints > 0 )
            npoints *= lb_yPoints / yPoints;
        else
            npoints *= lb_yPoints / 2;
        
        yPoints  = lb_yPoints;
    }

    // vector of y-coords of npoints lines parrallel to x-axis
    // IDEA: Also use these as x_sample_guides !!!
    gsMatrix<T> y_samples = gsPointGrid<T>( m_bbox(1,0), m_bbox(1,1), yPoints);

    // x-coords of sampling points on the line y=yi (if any)
    // int: yi index,
    //gsVector *: a pointer to a vector of npoints x-coords of sample points on the line y=yi
    //std::vector<index_t>: a vector of indices that indicate
    //the line segments of y=yi that intesect the domain

    std::vector<std::vector< VertexHandle > >samples;
    std::vector<std::vector< VertexHandle > >intersections;

//    std::map<int,
//            std::pair<
//            std::vector< VertexHandle >,
//            std::vector<T>
//            >
//            > x_samples;
    for ( int i = 0; i!= yPoints; ++i )
    {
        std::vector<T> x_all;
        //std::cout<<" --- before intersections  " << i <<", y="<< y_samples(0,i)  <<"\n";
        for (size_t j=0;j<m_loops.size();j++)
        {
            std::vector<T> x = m_loops[j]->lineIntersections(1, y_samples(0,i) );
            if ( ! x.empty() )
            {
                gsCurve<T> * curve = m_loops[j]->singleCurve() ; // TO BE REMOVED later

                if ( x.size() == 1 )
                {
                    x.push_back(x[0]);
                    gsWarn<<"lineIntersection yielded a tangent"<<'\n';
                }
                else if (x.size()%2==1)
                {
                    x.push_back(x[0]);
                    gsWarn<<"lineIntersection yielded a tangent and more intersections, try another number of points"<<'\n';
                }

                gsAsMatrix<T> xx(x);
                gsMatrix<T> e = curve->eval( xx );
                for (std::size_t k=0;k<x.size();k++)
                {
                    x_all.push_back(e(0,k));
                }

                delete curve;

            }
        }
        std::sort(x_all.begin(), x_all.end() );
        gsVector<int> numPoints( x_all.size() / 2 );
        int k(0);
        for ( typename std::vector<T>::const_iterator it = x_all.begin();
              it < x_all.end(); it += 2 )
            numPoints[k++] = static_cast<int>( (*(it+1) - *it)* npoints );

        std::vector<VertexHandle> x_line;
        std::vector<VertexHandle> intersection_vec;
        for(index_t j=0;j<numPoints.size();j++)
        {
            x_line.push_back(m->addVertex(x_all[j*2],y_samples(0,i)));
            intersection_vec.push_back(m->addVertex(x_all[j*2],y_samples(0,i)));

            for(int k=0;k<numPoints[j];k++)
            {
                x_line.push_back(m->addVertex(x_all[j*2]*(numPoints[j]-k)/(numPoints[j]+1)+x_all[j*2+1]*(k+1)/(numPoints[j]+1),y_samples(0,i)));
            }
            x_line.push_back(m->addVertex(x_all[j*2+1],y_samples(0,i)));
            intersection_vec.push_back(m->addVertex(x_all[j*2+1],y_samples(0,i)));
        }
        samples.push_back(x_line);
        intersections.push_back(intersection_vec);
    }
    T checkDist=(y_samples(0,1)-y_samples(0,0))*5;
    for(size_t i=0;i<intersections.size()-1;i++)
    {
        std::size_t currentTop=0;
        std::size_t currentBot=0;

        T topToBotDist=0;
        T botToTopDist=0;
        if(samples[i].size()>0&&samples[i+1].size()>0)
        {
            while (currentTop!=samples[i].size()-1||currentBot!=samples[i+1].size()-1)
            {
                if(currentTop==samples[i].size()-1)
                {
                    VertexHandle v1=samples[i][currentTop];
                    VertexHandle v2=samples[i+1][currentBot];
                    VertexHandle v3=samples[i+1][currentBot+1];
                    if(getDistance(v1,v2)<checkDist&&getDistance(v1,v3)<checkDist&&getDistance(v2,v3)<checkDist)
                        m->addFace(v1,v2,v3);
                    currentBot++;
                }
                else if(currentBot==samples[i+1].size()-1)
                {
                    VertexHandle v1=samples[i][currentTop];
                    VertexHandle v2=samples[i+1][currentBot];
                    VertexHandle v3=samples[i][currentTop+1];
                    if(getDistance(v1,v2)<checkDist&&getDistance(v1,v3)<checkDist&&getDistance(v2,v3)<checkDist)
                        m->addFace(v1,v2,v3);
                    currentTop++;
                }
                else
                {
                    topToBotDist=getDistance(samples[i][currentTop],samples[i+1][currentBot+1]);
                    botToTopDist=getDistance(samples[i][currentTop+1],samples[i+1][currentBot]);
                    if (topToBotDist<botToTopDist)
                    {
                        VertexHandle v1=samples[i][currentTop];
                        VertexHandle v2=samples[i+1][currentBot];
                        VertexHandle v3=samples[i+1][currentBot+1];
                        if(getDistance(v1,v2)<checkDist&&getDistance(v1,v3)<checkDist&&getDistance(v2,v3)<checkDist)
                            m->addFace(v1,v2,v3);
                        currentBot++;
                    }
                    else
                    {
                        VertexHandle v1=samples[i][currentTop];
                        VertexHandle v2=samples[i+1][currentBot];
                        VertexHandle v3=samples[i][currentTop+1];
                        if(getDistance(v1,v2)<checkDist&&getDistance(v1,v3)<checkDist&&getDistance(v2,v3)<checkDist)
                            m->addFace(v1,v2,v3);
                        currentTop++;
                    }
                }
            }
        }

    }


//            // generate the Mesh vertices.
//            k = npoints - x.size(); // points to distribute to the (interior of) the segments
//            std::vector< VertexHandle > x_line;
//            x_line.reserve(npoints);
//            for (size_t v=0; v< x.size(); v+=2 )
//            {
//                x_line.push_back( m->addVertex(x[v], y_samples(0,i)) );
//                // TO DO
//                // ( ll[v/2] / ltotal )*k samples here
//                for( int j = 1; j!= k+1; ++j)
//                    x_line.push_back( m->addVertex(x[v] + T(j) * ll[v/2] / T(k+1), y_samples(0,i)) );

//                x_line.push_back( m->addVertex(x[v+1], y_samples(0,i) ) );
//            }
//            //std::cout<<"i: "<< i<< ", sz="<< x_line.size() <<", ints="<< x.size() <<"\n";
//            //std::cout<<"x= "<< x[0] <<", "<< x[1] <<"(y="<< y_samples(0,i) <<"\n";

//            x_samples.insert( make_pair(i,make_pair( x_line, x ) ) ); // i-th sample line


//    }

    //std::cout<<" --- mesh DONE  "<< *m  <<"\n";
    // for ( int i = 0; i!= 25; ++i )
    // 	std::cout<<"--"<<* m->vertex[i] ;

    // Zig-zag connection
//    for ( int i = 0; i!= npoints-1; ++i )
//    {
//        typename std::map<int,std::pair<std::vector<VertexHandle>, std::vector<T> > >::const_iterator
//            it0 = x_samples.find(i);

//        if ( it0 != x_samples.end() ) // Check whether intersection with line 0 exists
//        {
//            typename std::map<int,std::pair<std::vector<VertexHandle>, std::vector<T> > >:: const_iterator
//                it1 = x_samples.find(i+1);
//            if ( it1 != x_samples.end() ) // Check whether intersection with line 1 exists
//            {
//                //typename std::vector<T>::const_iterator s0 = it0->second.second.begin();
//                //typename std::vector<T>::const_iterator s1 = it1->second.second.begin();
//                typename std::vector<VertexHandle>::const_iterator f1 = it1->second.first.begin();
//                // Run through the sample points at line 0
//                for ( typename std::vector<VertexHandle>::const_iterator f0 =
//                        it0->second.first.begin();
//                        f0 != it0->second.first.end()-1 ; ++f0 )
//                {
//                    // TO DO
//                    //if ( *f0 == *s0 )// start segment level 0
//                    // if ( *f1 == *s1 )// start segment level 1
//                    m->addFace( *f0    , *f1    , *(f1+1), *(f0+1) ) ;
//                    f1++;
//                }
//            }
//        }
//    }
    //gsInfo<<*m;
    return m;

}


}
