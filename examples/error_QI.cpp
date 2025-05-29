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

    gsMatrix<> Cfine, Ccoarse,Ccoarse_l2, Ccoarse_sch, Ccoarse_tay, Ccoarse_l2_local, Ccoarse_taylor;
    

    fd1.getId(0,dbasis_fine);
    fd2.getId(0,dbasis_coarse);
    fd1.getId(1,Cfine);

    gsStopwatch clock;
    real_t taytime = 0;
    real_t l2time = 0;
    real_t localqi = 0;


    gsGeometry<>::uPtr Cfine_ = dbasis_fine.basis(0).makeGeometry(give(Cfine));
    gsInfo<<Cfine.size()<<"\n";
    gsInfo<<dbasis_fine.basis(0)<<"\n";


    clock.restart();
    gsQuasiInterpolate<real_t>::localIntpl(dbasis_coarse.basis(0),*Cfine_,Ccoarse);
    gsGeometry<>::uPtr Ccoarse_ = dbasis_coarse.basis(0).makeGeometry(give(Ccoarse));
    localqi += clock.stop();
    
    clock.restart();
    gsL2Projection<real_t>::projectFunction(dbasis_fine.basis(0), dbasis_coarse,*Cfine_,mp,Ccoarse_l2);
    gsGeometry<>::uPtr Ccoarse_l2_ = dbasis_coarse.basis(0).makeGeometry(give(Ccoarse_l2));
    l2time += clock.stop();

    gsQuasiInterpolate<real_t>::Schoenberg(dbasis_coarse.basis(0),*Cfine_,Ccoarse_sch);
    gsGeometry<>::uPtr Ccoarse_sch_ = dbasis_coarse.basis(0).makeGeometry(give(Ccoarse_sch));

    gsDebugVar(dbasis_coarse.basis(0).size());
    gsDebugVar(dbasis_fine.basis(0).size());

    gsQuasiInterpolate<real_t>::localL2(dbasis_coarse.basis(0),*Cfine_,Ccoarse_l2_local);
    gsGeometry<>::uPtr Ccoarse_local_L2 = dbasis_coarse.basis(0).makeGeometry(give(Ccoarse_l2_local));


    clock.restart();
    gsQuasiInterpolate<real_t>::localTaylor(dbasis_coarse.basis(0),*Cfine_,2,Ccoarse_taylor);
    gsGeometry<>::uPtr Ccoarse_taylor_geom = dbasis_coarse.basis(0).makeGeometry(give(Ccoarse_taylor));
    taytime += clock.stop();

    gsDebugVar(l2time);
    gsDebugVar(taytime);
    gsDebugVar(localqi);

    // GLOBAL INTERPOLATION
    gsMatrix<> anchors = dbasis_coarse.basis(0).anchors();
    gsSparseMatrix<> cmatrix = dbasis_coarse.basis(0).collocationMatrix(anchors);
    gsMatrix<> vals = Cfine_->eval(anchors);
    vals.transposeInPlace();
    gsSparseSolver<>::LU solver;
    solver.compute(cmatrix);
    gsMatrix<> coefs = solver.solve(vals);

    gsGeometry<>::uPtr Ccoarse_global_interpl = dbasis_coarse.basis(0).makeGeometry(give(coefs));


    // First basis for projection and then integration basis! (inverted)
    // static void projectFunctionLocal(      const gsBasis<T> &b,
    //                                             const gsMultiBasis<T>   &intbasis,
    //                                             const gsFunction<T>  &source,
    //                                             const gsMultiPatch<T>   &geometry,
    //                                             gsMatrix<T> & result);
                                                

    // // Taylor (only for 1D)
    // gsQuasiInterpolate<real_t>::Taylor(dbasis_coarse.basis(0),*Cfine_,2,Ccoarse_tay);
    // // I need a gsBSplineBasis!!!!!!!
    // // gsBSplineBasis
    // gsGeometry<>::uPtr Ccoarse_tay_ = dbasis_coarse.basis(0).makeGeometry(give(Ccoarse_tay));

    // // Eval Based
    // gsQuasiInterpolate<real_t>::EvalBased(dbasis.basis(0),*Cnew_,false,CnewF);
    
    // gsWriteParaview(mp,*Ccoarse_,"Coarse"); // its a ptr
    // gsWriteParaview(mp,*Cfine_,"Fine"); // its a ptr

    // gsExprAssembler<> A(1,1);
    // A.setIntegrationElements(dbasis_fine); // para hacer plot en coarse o en fine mesh!!!!! importantisimo
    // gsExprEvaluator<> ev(A);

    // auto G = ev.getMap(mp);

    // auto cfine = ev.getVariable(*Cfine_,G);
    // auto ccoarse = ev.getVariable(*Ccoarse_,G); //quasi interpolation
    // auto ccoarse_l2 = ev.getVariable(*Ccoarse_l2_,G); // L2 interpolation
    // auto ccoarse_sch = ev.getVariable(*Ccoarse_sch_,G); // quasi interpolation schoenberg
    // auto ccoarse_l2_local = ev.getVariable(*Ccoarse_local_L2,G); // quasi interpolation schoenberg
    // auto ccoarse_tay = ev.getVariable(*Ccoarse_taylor_geom,G); // quasi interpolation schoenberg
    // auto ccoarse_gl  = ev.getVariable(*Ccoarse_global_interpl,G); // global interpolation

    // // ev.writeParaview(cfine,G,"Finesol");
    // // ev.writeParaview(ccoarse,G,"Coarsesol");

    // // ev.writeParaview((cfine-ccoarse).sqNorm(),G,"errordiff");

    // gsMatrix<> vals_gl = Ccoarse_global_interpl->eval(anchors);
    // gsInfo<<"dvals = "<<(vals-vals_gl).norm()<<"\n";

    // gsInfo<<"Integral of the difference of the coefficients\n";
    // gsInfo<<"QI:  "<<(ev.integral(meas(G) * (ccoarse-cfine).sqNorm()))<<"\n";
    // gsInfo<<"L2:  "<<(ev.integral(meas(G) * (ccoarse_l2-cfine).sqNorm()))<<"\n";
    // gsInfo<<"Sch: "<<(ev.integral(meas(G) * (ccoarse_sch-cfine).sqNorm()))<<"\n";
    // gsInfo<<"L2L: "<<(ev.integral(meas(G) * (ccoarse_l2_local-cfine).sqNorm()))<<"\n";
    // gsInfo<<"T:   "<<(ev.integral(meas(G) * (ccoarse_tay-cfine).sqNorm()))<<"\n";
    // gsInfo<<"GL:   "<<(ev.integral(meas(G) * (ccoarse_gl-cfine).sqNorm()))<<"\n";
    
    // gsInfo<<"Difference of the integrals of the coefficients\n";
    // real_t cfine_int = ev.integral(cfine*meas(G));
    // gsInfo<<"QI:  "<<cfine_int-ev.integral(ccoarse*meas(G))<<"\n";
    // gsInfo<<"L2:  "<<cfine_int-ev.integral(ccoarse_l2*meas(G))<<"\n";
    // gsInfo<<"Sch: "<<cfine_int-ev.integral(ccoarse_sch*meas(G))<<"\n";
    // gsInfo<<"L2L: "<<cfine_int-ev.integral(ccoarse_l2_local*meas(G))<<"\n";
    // gsInfo<<"T:   "<<cfine_int-ev.integral(ccoarse_tay*meas(G))<<"\n";
    // gsInfo<<"GL:  "<<cfine_int-ev.integral(ccoarse_gl*meas(G))<<"\n";

    // A.setIntegrationElements(dbasis_coarse); // to do plot in the coarse mesh!
    // gsParaviewCollection error_collection("Error_Analysis", &ev);
    // error_collection.options().setSwitch("plotElements", true);
    // error_collection.options().setInt("plotElements.resolution", 4);
    // error_collection.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 10000);
    // error_collection.options().setInt("precision", 40); // 1e-18



    // if (plot_error)
    // {
    //     error_collection.newTimeStep(&mp); // i need an updated mp!
    //     error_collection.addField((ccoarse-cfine).sqNorm(),"error QI");
    //     error_collection.addField((ccoarse_l2-cfine).sqNorm(),"error L2");
    //     error_collection.addField((ccoarse_sch-cfine).sqNorm(),"error schoenberg");
    //     error_collection.addField((ccoarse_l2_local-cfine).sqNorm(),"error L2 local");
    //     error_collection.addField((ccoarse_l2_local-ccoarse).sqNorm(),"difference");
    //     // error_collection.addField((ccoarse_tay-cfine).sqNorm(),"error taylor");
    //     error_collection.addField((ccoarse_gl-cfine).sqNorm(),"error global interpl");
    //     error_collection.saveTimeStep();
    //     error_collection.save();
    // }

    //como se si es la malla fina o la gruesa???????????????
    //q comparo de todo
    //l2 proyeccion

}// end main
