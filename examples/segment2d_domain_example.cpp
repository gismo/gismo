// example segmentation

#include <iostream>
#include <math.h>

#include <gismo.h>

#include <gsSolverDev/gsBemUtils.h>


template<class T>
void segment(gsPlanarDomain<T> & pdomain,  int n_points, T tolerance,bool circle, gsMesh<T> &segmentation)
 // gsMatrix<T> const & start, int n_points,  tolerance
{
    // ***************************   Consider a template
gsTemplate<T> * pmy_template;
    if(circle)
    pmy_template = new gsTemplate<T>(pdomain.numHoles() ); // Convex template with no holes

    else
    {
       bool sq=true;
       pmy_template= new gsTemplate<T>(sq,1);
    }

    gsTemplate<T>& my_template (*pmy_template);

     // ***************************   first step: solve Laplace problems


    // call : std::pair< gsFunction<T> *, gsFunction<T> *, > S = pdomain.mapto( template );
    // u.first, u.second are the components of the solution
    std::pair <gsFunction<T>* , gsFunction<T>*> u;

    pdomain.mapto(pdomain,my_template, 2); // 2 for amoeba_hole, 3 for austria_hole
                                   // 4 for Puzzle1 domain


    // construct indice points: points inside the planar domain at a suitable distance from
    // the outer boundary

    gsMatrix<T> tab(2,2);
    tab= *safe(pdomain.boundingBox());
    gsMatrix<T> storage;
    storage= uniformPointGrid<T>( tab.col(0), tab.col(1),600);
    gsMatrix<T> indice (2,storage.cols()); // matrix indice stores points of the grid which are
                                           // inside the planar domain
    indice.setZero();
    int k=0;

    gsBSpline<T>  b =*safe(dynamic_cast<gsBSpline<T> * > ( pdomain.outer().singleCurve() ));

    gsMatrix<T> mod = b.coefs();
    mod *=0.80;//0.85
    gsKnotVector<T> kv = b.basis().knots();

    gsBSpline<T> *interior = new gsBSpline<T>(kv, mod);


    gsCurveLoop<T> * inner = new gsCurveLoop<T>( interior );
    gsPlanarDomain<T> Idomain( inner);
    
    if(pdomain.numHoles()!=0)
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
std::pair<gsFunction<T>* , gsFunction<T>* >
mapto(gsPlanarDomain<T> & pdomain, gsTemplate<T> & my_template, int numRefine)
{
  //  GISMO_ASSERT( my_template.domain().numLoops()== pdomain.numLoops(), "Template has different topology.");

    //******* basis vector of all basis functions of all curves in pd
    std::vector < gsBasis<T> * > basis;

    std::vector<gsFunction<T>*> bc_0;
    std::vector<gsFunction<T>*> bc_1;
    gsMatrix<T> param(1,1);
    param<<1;
    gsBSpline<T> * tmpl_loop = dynamic_cast<gsBSpline<T> * >( my_template.loop(0).singleCurve());

    for(index_t v=0; v<=pdomain.numHoles(); v++)
    {
        gsBSpline<T> * this_loop = dynamic_cast<gsBSpline<T> * >( pdomain.loop(v).singleCurve() );

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
//             if(pdomain.numLoops()==1)
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

    for(index_t i=0; i<pdomain.numLoops() ; i++)
    {
        for (int k = 0; k < numRefine; ++k)
            basis[i]->uniformRefine();
    }

    int count = 0;
    for(index_t i=0; i< pdomain.numLoops() ; i++)
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
     if(pdomain.numHoles()!=0)
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
//     gsMatrix<T> cent( 2,pdomain.numHoles() ), Imag_cent ;
     
//     for(int s=1;s<=pdomain.numHoles();++s)
//     {
//         gsBSpline<T> * inner_component = dynamic_cast<gsBSpline<T> * >( pdomain.loop(s).singleCurve() );
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


using namespace std;
using namespace gismo;

int main(int argc, char *argv[])
{
    gsPlanarDomain<> * Pdomain = NULL;
    int n_points(10);
    double tolerance = 1e-4;
    unsigned smpl(30);
    bool plot = false;
    std::vector<real_t> start_v;
    std::string fn("");

    try
    {
        gsCmdLine cmd("Segmentation in quadrangular patches of a planar domain ");

        gsArgVal<std::string> a5("g","geometry","File containing Geometry (.axl, .txt)",
                                 false,"", "geometry file", cmd );
        gsArgSwitch a4("", "plot", "Plot result with ParaView ", cmd);
        gsArgMultiVal<real_t> a3("i","starting_point",
                                 "point inside the computational domain from which looking for the starting point of the tracing curves", 
                                 false,"fixed starting point", cmd );
        gsArgVal<int> a2("s","samples", 
                         "Number of samples", 
                         false,smpl, "samples", cmd );
        gsArgVal<int> a1("p","n_points", 
                         "Number of points traced per curve", 
                         false,n_points, "points traced", cmd );
        gsArgVal<> a0("t","tolerance", 
                            "Required accuracy", 
                            false,tolerance, "tolerance", cmd );
        
        cmd.parse(argc,argv);
        
        
        fn = a5.getValue();
        
        if (fn.empty() )
        {
            fn = GISMO_SOURCE_DIR;
            fn+="/filedata/planar/amoeba0_pdomain.xml";
        }
        
        // Pdomain =  gsReadFile<>( fn ) ;
        gsFileData<>  filedata(fn);
        
        if( filedata.has<gsPlanarDomain<> >() )
            Pdomain = filedata.getFirst< gsPlanarDomain<> >();
        
        else
        {
            
            if(filedata.has<gsCurve <> >() )
            {
                gsCurve<> * geo = filedata.getAnyFirst< gsCurve<> >();
                gsCurveLoop<> * cp = new gsCurveLoop<>( geo );
                Pdomain = new gsPlanarDomain<>( cp );
            }
            
            
            else
            {
                gsWarn<< "Did not find any planar domain or geometry in "<< fn<<", quitting.\n";
                return false;
            }
        }
        
        plot = a4.getValue();
        start_v = a3.getValue();
        if(start_v.empty())
        {
            start_v.push_back(1.5);
            start_v.push_back(2);
        }
        
        smpl = a2.getValue();
        n_points = a1.getValue();
        tolerance = a0.getValue();
    } catch ( gsArgException& e )
    { cout << "Error: " << e.error() << " " << e.argId() << endl; }
    
    
    std::cout << "------------------------------------------------------------"
                 "\nInput Arguments: \n\n"
                 "input geometry: " << fn <<
                 "\n how many samples:  "<<smpl<<
                 "\n number of points traced per curve: "<<n_points<<
                 "\n tolerance: "<<tolerance<<"\n"
                 "------------------------------------------------------------"
                 "\n\n";


    gsAsMatrix<> start (start_v,2,1);
    
    bool c= true;
    gsMesh<> res;
    
    segment(*Pdomain, n_points, tolerance, c, res );
    //  cout<<"\n res \n "<<res << "\n";

    int exitCommand = 0;
    if(plot)
    {
        
        gsWriteParaview<>( *Pdomain, "Pdomain", smpl) ;
        
        gsWriteParaview<>( res, "Segmentation");

        //run: paraview        
        exitCommand = system("paraview Segmentation.vtp &");        
    }

    delete Pdomain;
    
    return exitCommand;
}
