/** @file cahn-hilliard.cpp

    @brief Tutorial on how to use expression assembler to solve the Cahn-Hilliard equation

*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    bool plot_error = true;

    gsMultiPatch<> mp;
    mp.addPatch(*gsNurbsCreator<>:: BSplineSquare()); // create a multipatch with multibasis!???????
    
    gsFileData<> fd1, fd2;

    std::string fn1 = "finebasis.xml";
    std::string fn2 = "coarsebasis.xml";

    fd1.read(fn1);
    fd2.read(fn2);
 
    gsMultiBasis<> dbasis_fine, dbasis_coarse;

    gsMatrix<> Cfine, Ccoarse,Ccoarse_l2, Ccoarse_sch, Ccoarse_tay;
    

    fd1.getId(0,dbasis_fine);
    fd2.getId(0,dbasis_coarse);
    fd1.getId(1,Cfine);

    
    gsGeometry<>::uPtr Cfine_ = dbasis_fine.basis(0).makeGeometry(give(Cfine));

    gsQuasiInterpolate<real_t>::localIntpl(dbasis_coarse.basis(0),*Cfine_,Ccoarse);
    gsGeometry<>::uPtr Ccoarse_ = dbasis_coarse.basis(0).makeGeometry(give(Ccoarse));
    
    gsL2Projection<real_t>::projectFunction(dbasis_fine.basis(0), dbasis_coarse,*Cfine_,mp,Ccoarse_l2);
    gsGeometry<>::uPtr Ccoarse_l2_ = dbasis_coarse.basis(0).makeGeometry(give(Ccoarse_l2));

    gsQuasiInterpolate<real_t>::Schoenberg(dbasis_coarse.basis(0),*Cfine_,Ccoarse_sch);
    gsGeometry<>::uPtr Ccoarse_sch_ = dbasis_coarse.basis(0).makeGeometry(give(Ccoarse_sch));

    // // Taylor
    // gsQuasiInterpolate<real_t>::Taylor(dbasis_coarse.basis(0),*Cfine_,2,Ccoarse_tay);
    // // I need a gsBSplineBasis!!!!!!!
    // // gsBSplineBasis
    // gsGeometry<>::uPtr Ccoarse_tay_ = dbasis_coarse.basis(0).makeGeometry(give(Ccoarse_tay));

    // // Eval Based
    // gsQuasiInterpolate<real_t>::EvalBased(dbasis.basis(0),*Cnew_,false,CnewF);
    
    // gsWriteParaview(mp,*Ccoarse_,"Coarse"); // its a ptr
    // gsWriteParaview(mp,*Cfine_,"Fine"); // its a ptr

    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis_fine); // para hacer plot en coarse o en fine mesh!!!!! importantisimo
    gsExprEvaluator<> ev(A);

    auto G = ev.getMap(mp);

    auto cfine = ev.getVariable(*Cfine_,G);
    auto ccoarse = ev.getVariable(*Ccoarse_,G); //quasi interpolation
    auto ccoarse_l2 = ev.getVariable(*Ccoarse_l2_,G); // L2 interpolation
    auto ccoarse_sch = ev.getVariable(*Ccoarse_sch_,G); // quasi interpolation schoenberg
    //auto ccoarse_tay = ev.getVariable(*Ccoarse_tay_,G); // quasi interpolation schoenberg

    // ev.writeParaview(cfine,G,"Finesol");
    // ev.writeParaview(ccoarse,G,"Coarsesol");

    // ev.writeParaview((cfine-ccoarse).sqNorm(),G,"errordiff");


    gsDebugVar(ev.integralElWise(meas(G) * (ccoarse-cfine).sqNorm()));
    gsDebugVar(ev.integralElWise(meas(G) * (ccoarse_l2-cfine).sqNorm()));
    gsDebugVar(ev.integralElWise(meas(G) * (ccoarse_sch-cfine).sqNorm()));
    

    A.setIntegrationElements(dbasis_coarse); // to do plot in the coarse mesh!
    gsParaviewCollection error_collection("Error_Analysis", &ev);
    error_collection.options().setSwitch("plotElements", true);
    error_collection.options().setInt("plotElements.resolution", 4);
    error_collection.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 5000);



    if (plot_error)
    {
        error_collection.newTimeStep(&mp); // i need an updated mp!
        error_collection.addField((ccoarse-cfine).sqNorm(),"error QI");
        error_collection.addField((ccoarse_l2-cfine).sqNorm(),"error L2");
        error_collection.addField((ccoarse_sch-cfine).sqNorm(),"error schoenberg");
        //error_collection.addField((ccoarse_tay-cfine).sqNorm(),"error taylor");
        error_collection.saveTimeStep();
        error_collection.save();
    }

    //como se si es la malla fina o la gruesa???????????????
    //q comparo de todo
    //l2 proyeccion

}// end main
