 
#pragma once

#include <iostream>
#include <vector>

#include <gsCore/gsGeometry.h>
#include <gsCore/gsField.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsFunction.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsUtils/gsPointGrid.h>

#include <cpplot/cpplot.hpp> // External file

//#include <boost/shared_ptr.hpp>
//#include <boost/enable_shared_from_this.hpp>
//#include <boost/thread.hpp>
//#include <boost/thread/mutex.hpp>


namespace gismo
{
  /** 
      Class for scientific plotting
  */
  
template<class T>
class gsPlot
{

public:
typedef std::vector<double> dvec;
typedef std::vector< std::vector<double> > dmat;
typedef std::vector< std::vector<float> > tcvec;
typedef std::vector< std::vector< std::vector<float> > > tcmat;

public:

  /// Default empty constructor
  gsPlot() { };

  ~gsPlot() { }; //destructor


public:

// MatPlot
//void DISPLAY();

/// Plots a function
void function( gsFunction<T>  const & f)
{
    if ( fork_me() ) return;

    using namespace cpplot;
//    int dim = 1;
    std::vector< dvec > data;   
    std::vector< dmat > Zdata;

    int n = 100;
    gsMatrix<double> pts = gsPointGrid(0.0, 1.0, n) ;
    gsMatrix<T> * ev  = f.eval( pts ) ;

    //std::cout<<"eval cols "<< ev->cols() <<"\n";

    // for ( int i = 0; i< dim; ++i )
    // {
    //     std::vector<double> x(n);
    //     data.push_back( x ) ;
    // }
    std::vector<double> x(n);
    data.push_back( x ) ;
    data.push_back( x ) ;
    
    for ( index_t j= 0; j < ev->cols() ; ++j )
    {
        data[0][j] = pts(0,j) ; // x-coordinate
        data[1][j] = ev->at(0,j);
    }
    
    delete ev;

    figure_start();

    for ( unsigned j= 1; j < data.size(); ++j )
        plot( data[0], data[j]);
    figure_end();

    // Exit child process
    exit(0);
};


/// Plots a basis, 1D, 2D.
void basis( gsBasis<T> const & basis)
{
    if ( ! fork_me() ) return; // for debugging the plotting process
    //if ( fork_me() ) return;

    using namespace cpplot;
    int dim = basis.dim() ;
    std::vector< dvec > data;   
    std::vector< dmat > Zdata;


    gsVector<double> a ; 
    a=  basis.parameterRange()->col(0) ;
    gsVector<double> b ;
    b= basis.parameterRange()->col(1) ;
    gsVector<unsigned> np(dim);
    for ( int i = 0; i< dim; ++i )
        //np(i) = 40;
        np(i) = 5;

    int n= np(0);
    gsMatrix<double> pts = gsPointGrid(a,b,np) ;
    gsMatrix<T> * ev  = basis.eval( pts ) ;
    gsMatrix<unsigned> * act = basis.active( pts ) ;
    
    std::cout<<"eval\n"<< *ev;
    std::cout<<"act \n"<< *act;


    for ( int i = 0; i< dim; ++i )
    {
        data.push_back( linspace(a(i),b(i),n) ) ;
    }
    
    for ( index_t i = 0; i< basis.size(); ++i )
    {
        std::vector<double> x(n);
        data.push_back( x ) ;
    }

    if ( dim== 1 )
    {
        for ( index_t j= 0; j < ev->cols() ; ++j )
        {
            for ( index_t k= 0; k < ev->rows() ; ++k )
                data[act->at(k,j)+dim][j] = ev->at(k,j);
        }
        delete ev;
        delete act;


        figure_start();

        for ( unsigned j= 1; j < data.size(); ++j )
            plot( data[0], data[j]);

        figure_end();
    }
    
    if ( dim== 2 )
    {
        std::cout<<"plot "<< basis.size() <<"\n";
        for ( index_t i = 0; i< basis.size(); ++i )
            Zdata.push_back( dmat(n,dvec(n)) ) ;
        
std::cout<<"plot\n";

        for ( int i= 0; i < n ; ++i )
            for ( int j= 0; j < n ; ++j )
                for ( index_t k= 0; k < ev->rows() ; ++k )
                        Zdata[act->at(k,i*n+j)][i][j] = ev->at(k,i*n+j);

std::cout<<"plot OK \n";

        delete ev;
        delete act;

        std::vector<std::string> cc;
        cc.push_back("r");
        cc.push_back("g");
        cc.push_back("b");
        cc.push_back("y");
        cc.push_back("c");
        cc.push_back("p");

std::cout<<"fig \n"<< *act;

        figure_start();

         for ( unsigned j= 0; j < Zdata.size(); ++j )
         {
             surface_t p = surf ( data[0], data[1], Zdata[j]);
             //p->set("EdgeColor",   "none"   );
             p->set("EdgeColor",   "w"   );
             //p->set("FaceColor", cc[j%6] ); //cc[ j%6 ]
         }
         figure_end();
    }   
std::cout<<"Exit \n"<< *act;
    // Exit child process
    exit(0);
};

/// Plots a basis, 1D, 2D.
void geometry( gsGeometry<T> const & geo)
{
    if ( fork_me() ) return; // debug..
    using namespace cpplot;
    int dim = geo.geoDim() ;
    int pdim= geo.parDim();

    gsVector<double> a ; 
    a=  geo.basis().parameterRange()->col(0) ;
    gsVector<double> b ;
    b= geo.basis().parameterRange()->col(1) ;
    gsVector<unsigned> np(pdim);
    for ( int i = 0; i< pdim; ++i )
        np(i) = 30;

    int n= np(0);
    //int n= np(0);
    gsMatrix<double> pts = gsPointGrid(a,b,np) ;
    gsMatrix<T> * ev  = geo.eval( pts ) ;  

    if ( dim== 3 )
    {
        dmat Xdata(n,dvec(n));
        dmat Ydata(n,dvec(n));
        dmat Zdata(n,dvec(n));
        for ( int i= 0; i < n ; ++i )
            for ( int j= 0; j < n ; ++j )
            {
                Xdata[i][j] = ev->at(0,i*n+j);
                Ydata[i][j] = ev->at(1,i*n+j);
                Zdata[i][j] = ev->at(2,i*n+j);
            }
        delete ev;
        

        figure_start();

        spring();

        surface_t p = surface ( Xdata, Ydata, Zdata);
        p->set("EdgeColor","none");
        //p->set("FaceColor","flat");

        figure_end();
    }

    if ( dim== 2 )
    {
        dmat Xdata(n,dvec(n));
        dmat Ydata(n,dvec(n));
        dmat Zdata(n,dvec(n));

        for ( int i= 0; i < n ; ++i )
            for ( int j= 0; j < n ; ++j )
            {
                Xdata[i][j] = ev->at(0,i*n+j);
                Ydata[i][j] = ev->at(1,i*n+j);
            }
        delete ev;
        
        figure_start();

        surface_t p = pcolor ( Xdata, Ydata,Zdata);
        p->set("FaceColor","b");
        p->set("EdgeColor","none");

        figure_end();
    }

    // Exit child process
    exit(0);
};



/// Plots a basis, 1D, 2D.
void field( gsField<T> const & fl)
{
    
    if ( fork_me() ) return;
    using namespace cpplot;
    int dim = fl.geometry()->geoDim() ;
    int pdim= fl.geometry()->parDim();

    gsVector<double> a ; 
    a=  fl.geometry()->basis().parameterRange()->col(0) ;
    gsVector<double> b ;
    b= fl.geometry()->basis().parameterRange()->col(1) ;
    gsVector<unsigned> np(pdim);
    for ( int i = 0; i< pdim; ++i )
        np(i) = 30;

    int n= np(0);
    //int n= np(0);
    gsMatrix<double> pts = gsPointGrid(a,b,np) ;
    gsMatrix<T> * ev  = fl.point( pts ) ;  
    gsMatrix<T> * sc  = fl.pvalue( *ev ) ;  

    if ( dim== 3 )
    {
        dmat Xdata(n,dvec(n));
        dmat Ydata(n,dvec(n));
        dmat Zdata(n,dvec(n));
        dmat Cdata(n,dvec(n));

        for ( int i= 0; i < n ; ++i )
            for ( int j= 0; j < n ; ++j )
            {
                Xdata[i][j] = ev->at(0,i*n+j);
                Ydata[i][j] = ev->at(1,i*n+j);
                Zdata[i][j] = ev->at(2,i*n+j);
                Cdata[i][j] = sc->at(0,i*n+j);
            }
        delete ev;
        
        figure_start();

        surface ( Xdata, Ydata, Zdata );
        //set("EdgeColor","k");
        //set("FaceColor","flat");

        figure_end();
    }

    // Exit child process
    exit(0);
};





private:

    bool fork_me()
        {
            pid_t pID = fork();
            if (pID < 0) 
            {
                // Throw exception
                std::cerr << "Failed to fork" << std::endl;
                exit(1);
            }
            else if (pID != 0)
            {
                // True returned to parent process
                return true; 
            }
            // False returned to child process
            return false; 
        }



}; // class gsPlot


//////////////////////////////////////////////////
//////////////////////////////////////////////////




}; // namespace gismo
