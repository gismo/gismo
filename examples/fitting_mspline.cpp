/// gsFitting_general.cpp
/// Author:Gabor Kiss
/// for testing the class gsFitting

#include <iostream>
#include <time.h>

#include <gismo.h>




using namespace gismo;


int main(int argc, char *argv[])
{
    bool plot = false;
    //measuring the computational time
    //int clo=clock();
    std::string filename = "face.xml";
    gsGeometry<>::uPtr hbs;

    int np = 500;
    int nd = 25;
    int  err_type = 1;
    int function = 1;
    int eps = 0.01;
    
    gsCmdLine cmd("Hi, give me a file (.xml) with some multipatch geometry and basis.");
    cmd.addPlainString("filename", "File containing mp geometry and basis (.xml).", filename);
    cmd.addInt("f", "function", "Number of the function", function);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
     
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // The file data
    gsFileData<> data( filename );
    gsFunctionExpr<> f;
    gsMultiPatch<> mp;
    gsMultiBasis<> mb;
    gsSparseMatrix<> cf;
    gsMappedBasis<2, real_t> mbasis;

    // Load data as multipatch structure

    data.getFirst(mp);
    data.getFirst(mb);
    data.getFirst(cf);

    mbasis.init(mb, cf);

    gsInfo << "Number of patches: " << mp.nPatches() << "\n";
    gsInfo << "Number of points per patch: " << nd << "\n";
    gsInfo << "Number of basis functions: " << mbasis.size() << "\n";
    //gsInfo << "Basis" << mb << "\n";


    switch (function) {
    case 1:
        f = gsFunctionExpr<>("(exp(52*sqrt((10*x-2)^2+(10*y-3)^2)))^(-1)", 3);//one peak
        gsInfo << "Source function: " << f << "\n";
        break;
    case 2:
        f = gsFunctionExpr<>("y/2*(cos(4*(x^2+y-1)))^4", 3);//waves
        gsInfo << "Source function: " << f << "\n";
        break;
    case 3:
        f = gsFunctionExpr<>("3/4*exp(-((9*x-2)^2 + (9*y-2)^2)/4)+3/4*exp(-((9*x+1)^2)/49 - (9*y+1)/10)+1/2*exp(-((9*x-7)^2 + (9*y-3)^2)/4)-1/5*exp(-(9*x-4)^2 - (9*y-7)^2)", 3);//hils
        gsInfo << "Source function: " << f << "\n";
        break;
    case 4:
        f = gsFunctionExpr<>("sin((x*x+y*y+2/(5*pi))^(-1))", 3);//crater
        gsInfo << "Source function: " << f << "\n";
        break;
    case 5:
        f = gsFunctionExpr<>("(exp(2*sqrt((x)^2+(y)^2)))^(-1)", 3);//cusp
        gsInfo << "Source function: " << f << "\n";
        break;
    case 6:
        f = gsFunctionExpr<>("(exp(2*sqrt((10*x+3)^2+(10*y-3)^2)))^(-1) + (exp(2*sqrt((10*x-3)^2+(10*y+3)^2)))^(-1)", 3);//2 peaks
        gsInfo << "Source function: " << f << "\n";
        break;
    case 7:
        f = gsFunctionExpr<>("(1.5*exp(sqrt((10*x-3)^2+(10*y-3)^2)))^(-1)+ (1.5*exp(sqrt((10*x+3)^2+(10*y+3)^2)))^(-1) + (1.5*exp(sqrt((10*x)^2+(10*y)^2)))^(-1)", 3);//3 peaks
        gsInfo << "Source function: " << f << "\n";
        break;
    case 8:
        f = gsFunctionExpr<>("(x+y)/2+sqrt(((x-y)/2)^2)", 3);//Rvachev functions (R-functions)
        gsInfo << "Source function: " << f << "\n";
        break;
    default:
        gsInfo << "Unknown function, please pick one of the functions 1 - 8" << "\n";
        return 0;
        break;
    }

    int nbp = mp.nPatches();

    //gsDebugVar(nbp);

    gsVector<index_t> offset(nbp+1);
    offset[0] = 0;
    for (size_t i = 1; i<nbp; ++i)
    {
        offset[i] = offset[i - 1] + nd;
    }
    offset[nbp] = offset[nbp - 1] + nd;

    //gsDebugVar(offset.transpose());

    gsMatrix<real_t> Mpar (2 , nbp * nd);
    gsMatrix<real_t> fval (f.targetDim(), nbp * nd);

    gsMultiPatch<> ptc;
    

    // loop on the nb of patches

    for (int i = 0; i < nbp ; i++) {

        gsDebugVar(i);
        gsMatrix<> para = mp.patch(i).parameterRange();
        //gsDebugVar(ptc);
        gsVector<> c0 = para.col(0);
        gsVector<> c1 = para.col(1);
        c0.array() += eps;
        c1.array() -= eps;
        gsDebugVar(c0);
        gsDebugVar(c1); 
        //the parameter values for the fitting
        gsMatrix<> pts = uniformPointGrid(c0, c1, nd);
       //  gsMatrix<> pts(2, nd);
        //gsInfo << "Parameter values used for fitting: " << "\n" << pts << "\n";
        Mpar.middleCols(offset[i], nd) = pts;
        
        //the evaluated values of the original surface
        gsMatrix<> mp_eval = mp.patch(i).eval(pts);
        //gsInfo<<"Original surface points: "<<"\n"<< mp_eval<<"\n";
        //gsDebugVar(mp_eval.rows());
        //gsDebugVar(mp_eval.cols());

        //gsDebugVar(f.eval(mp_eval));
        fval.middleCols(offset[i], nd) = f.eval(mp_eval);

    }

    //gsDebugVar(Mpar);
    //gsDebugVar(fval);

    gsFileData<> out;
    out << Mpar;
    out << fval;
    out << gsMatrix<index_t>(offset); //conversion
    out.save("fitting_test1");

    gsDebugVar(out.lastPath()); 

   // end loop
    
   
    
    //create the fitting object
    gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    gsInfo<<"Creating the multipatch fitting object"<<"\n";
    gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    gsFitting<> fitting(Mpar, fval, offset, mbasis);


    gsInfo<<"fit class created"<<"\n";
    //gsMatrix<> results(1,6);
    fitting.compute(0); //0



    //gsGeometry<> * test;
    //test = fitting.result();

    //gsDebugVar(test);

    //gsTHBSpline<2>  * hbs1 = static_cast< gsTHBSpline<2>  *> (test);
    //std::vector<real_t> errors;
    //fitting.get_Error(errors, err_type);
    
    real_t error;
    fitting.computeApproxError(error, 0);
    gsDebugVar(error);
    /*real_t min = 1000000;
    real_t max = -1000000;
    for(unsigned int j =0; j < errors.size();j++){
        if(errors[j]>max){
            max = errors[j];
        }
        if(errors[j]<min){
            min = errors[j];
        }
    } */
    /*results(0, 0) = hbs1->basis().maxLevel();
    results(0,1) = hbs1->basis().size();
    results(0,2) = min;
    results(0,3) = max;
    results(0,5) = error;*/

    /*

    real_t num = 0;
    for(unsigned int j = 0; j < errors.size(); j++){
        if(errors[j]< 0.00001){
            num++;
        }
    }
    results(0,4) = (num*100)/errors.size();
    gsFileData<> newdata;
    newdata << *test ;
    plot = true;
    if(plot){
        newdata.dump("gsThbs_loc_face_4");
        gsWriteParaview( *test , "gsThbs_loc_face_4", np);
    }
    gsMatrix<> hbs1_eval =hbs1->eval(pts);
    gsInfo<<"number of points"<<errors.size()<<"\n";
    plot_errors( hbs_eval, hbs1_eval,errors, "gsThbs_loc_face_error_4");
    gsInfo<<"results"<<results<<"\n"<<"\n";
    gsInfo<<results(0,0);
    gsInfo<<" & "<<results(0,1);
    gsInfo.setf(std::ios::scientific);
    for(int j = 2; j < results.cols();j++){
        gsInfo<<" & "<<results(0,j);
    }
    gsInfo<<"\n"<<"Finished"<<"\n";



    gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    gsInfo<<"Creating the Tensor fitting object"<<"\n";
    gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    gsMatrix<> results(1,6);
    gsFitting<> * t_fitting = new gsFitting<>(pts, hbs_eval, T_tbasis);

    t_fitting->compute(0.001);
    gsGeometry<> * t_test;
    t_test = t_fitting->result();
    gsInfo<<*t_test<<"\n";

        gsGeometry<> * test;
    test = t_fitting->result();
    std::vector<real_t> errors;
    t_fitting->get_Error(errors, err_type);
    real_t error;
    t_fitting->computeApproxError(error, 0);
    real_t min = 1000000;
    real_t max = -1000000;
    for(unsigned int j =0; j < errors.size();j++){
        if(errors[j]>max){
            max = errors[j];
        }
        if(errors[j]<min){
            min = errors[j];
        }
    }
    results(0,0) = 0;
    results(0,1) = test->basis().size();
    results(0,2) = min;
    results(0,3) = max;
    results(0,5) = error;

    real_t num = 0;
    for(unsigned int j = 0; j < errors.size(); j++){
        if(errors[j]< 0.00001){
            num++;
        }
    }
    results(0,4) = (num*100)/errors.size();
    gsInfo<<results(0,0);
    gsInfo<<" & "<<results(0,1);
    gsInfo.setf(std::ios::scientific);
    for(int j = 2; j < results.cols();j++){
        gsInfo<<" & "<<results(0,j);
    }

    gsFileData<> t_newdata;
    t_newdata << *t_test ;
    t_newdata.dump("t_fitting_1");
    gsWriteParaview( test , "t_fitting_1", np);
    gsInfo<<"\n"<<"Finished"<<"\n";
*/
    



  return 0;
}
