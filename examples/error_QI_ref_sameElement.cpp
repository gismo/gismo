/** @file cahn-hilliard.cpp

    @brief Tutorial on how to use expression assembler to solve the Cahn-Hilliard equation

*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]


template <class T>
T computeNorm(const gsMultiBasis<T>  & integrationBasis,
                const gsFunctionSet<T> & currentBasis,
                const gsFunctionSet<T> & geometryMap,
                const gsFunctionSet<T> & sourceFunction1,
                const gsMatrix<T> & coefs,
                const gsOptionList     & options)
{
  
    // Create an assembler
    gsExprAssembler<T> A(1,1);
    A.options().update(options,gsOptionList::addIfUnknown); 

    // gsInfo << "Assembler options: " << A.options() << "\n";

    // Set the integration elements
    A.setIntegrationElements(integrationBasis);

    auto u = A.getSpace(currentBasis);
    u.setup();
    // gsMatrix<T> ccoarse = coefs; // coefs is const
    // auto sol= A.getSolution(u,ccoarse);

    // A.initSystem(); // it gives me seg fault?

    // Assign the geometry map
    typename gsExprAssembler<T>::geometryMap G = A.getMap(geometryMap);

    gsExprEvaluator<T> ev(A);
    // For the evaluator to use the options of the assembler (SameElement assumption!)
    ev.options().update(A.options(),gsOptionList::addIfUnknown);

    auto sol_before = ev.getVariable(sourceFunction1,G); // solution before projection      
    gsGeometry<>::uPtr sol_after_ptr  = currentBasis.basis(0).makeGeometry(coefs);
    auto sol_after = ev.getVariable(*sol_after_ptr,G);

    // gsInfo<<currentBasis.basis(0)<<"\n";
    // gsInfo << ev.options() << "\n";
    // gsInfo << A.options() << "\n";

    return ev.integral((sol_before-sol_after).sqNorm() * meas(G));

}

template <class T>
T basisUnion(const gsMultiBasis<T> & mb1, 
             const gsMultiBasis<T> & mb2,
             gsMatrix<T> & mbunion)
{
    return mbunion; 
}


int main(int argc, char *argv[])
{
    bool plot_error = true;

    gsMultiPatch<> mp;
    mp.addPatch(*gsNurbsCreator<>:: BSplineSquare()); // create a multipatch with multibasis!???????
 
    gsFileData<> fd1, fd2;

    std::string fn1 = "/home/lucasventa/gismo_random/basis_hugo_fine.xml"; // fine
    std::string fn2 = "/home/lucasventa/gismo_random/basis_hugo_crs.xml"; // coarse

    fd1.read(fn1);
    fd2.read(fn2);
 
    gsMultiBasis<> dbasis_fine, dbasis_coarse, dbasis_tensor;

    fd1.getId(0,dbasis_fine);
    fd2.getId(0,dbasis_coarse);
    // fd1.getId(1,Cfine);

    // Cast every basis of dbasis to a gsTHBSplineBasis
    for (size_t p=0; p!=dbasis_coarse.nBases(); p++)
    {
        GISMO_ASSERT((dynamic_cast<gsTHBSplineBasis<2,real_t>*>(&dbasis_coarse.basis(p))), "Basis is not a gsTHBSplineBasis");
        gsTHBSplineBasis<2,real_t> & b = static_cast<gsTHBSplineBasis<2,real_t> & >(dbasis_coarse.basis(p));
        dbasis_tensor.addBasis(b.tensorLevel(b.maxLevel()).clone());
    }

    gsMatrix<> Cfine, Ccoarse,Ccoarse_l2, Ccoarse_sch, Ccoarse_tay, Ccoarse_l2_local, Ccoarse_taylor, Cfine_l2;

    gsMatrix<> ccoarse = gsMatrix<>::Random(dbasis_coarse.size(),1);


    gsStopwatch clock;
    real_t taytime = 0;
    real_t l2time = 0;
    real_t localqi = 0;


    gsGeometry<>::uPtr Ccrs_ = dbasis_coarse.basis(0).makeGeometry(ccoarse);
    // gsInfo<<ccoarse.size()<<"\n";
    // gsInfo<<dbasis_coarse.basis(0)<<"\n";
    // gsInfo<<dbasis_fine.basis(0)<<"\n";


    clock.restart();
    gsQuasiInterpolate<real_t>::localIntpl(dbasis_fine.basis(0),*Ccrs_,Cfine);
    gsGeometry<>::uPtr Cfine_ = dbasis_fine.basis(0).makeGeometry(Cfine);
    localqi += clock.stop();
    
    // clock.restart();
    // real_t errorl2 = gsL2Projection<real_t>::projectFunction(dbasis_fine.basis(0), dbasis_fine,*Ccrs_,mp,Cfine_l2);
    // gsGeometry<>::uPtr Cfine_l2_ = dbasis_fine.basis(0).makeGeometry(give(Cfine_l2));
    // l2time += clock.stop();
    // gsInfo << "L2 error: " << errorl2 << "\n";


    // gsDebugVar(l2time);
    // gsDebugVar(taytime);
    // gsDebugVar(localqi);
    
    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis_fine); // para hacer plot en coarse o en fine mesh!!!!! importantisimo
    gsExprEvaluator<> ev(A);

    auto G = ev.getMap(mp);

    auto ccrs = ev.getVariable(*Ccrs_,G);
    auto cfine = ev.getVariable(*Cfine_,G); //quasi interpolation
    // auto cfine_l2 = ev.getVariable(*Cfine_l2_,G); // L2 interpolation
    // auto ccoarse_sch = ev.getVariable(*Ccoarse_sch_,G); // quasi interpolation schoenberg


    // gsInfo<<"Integral of the difference of the coefficients\n";
    // gsInfo<<"QI:  "<<(ev.integral(meas(G) * (ccrs-cfine).sqNorm()))<<"\n";

    // gsInfo<<"==================== ERRORS ====================\n";
    real_t error_1 = computeNorm(dbasis_tensor, dbasis_fine, mp, *Ccrs_, Cfine, A.options());
    real_t error_2 = computeNorm(dbasis_fine,   dbasis_fine, mp, *Ccrs_, Cfine, A.options());
    real_t error_3 = computeNorm(dbasis_coarse, dbasis_fine, mp, *Ccrs_, Cfine, A.options());

    // T computeNorm(const gsMultiBasis<T>  & integrationBasis,
    //             const gsFunctionSet<T> & currentBasis,
    //             const gsFunctionSet<T> & geometryMap,
    //             const gsFunctionSet<T> & sourceFunction1,
    //             const gsMatrix<T> & coefs,
    //             const gsOptionList     & options)

    gsInfo << "SameElement flag is: " 
       << (A.options().getSwitch("SameElement") ? "ON" : "OFF") << "\n";

    gsInfo<< "Error with integration on TPr basis (SE=ON): "<< error_1<<"\n";
    gsInfo<< "Error with integration on fine basis (SE=ON): "<< error_2<<"\n";
    gsInfo<< "Error with integration on crs basis (SE=ON): "<< error_3<<"\n";


    A.options().setSwitch("SameElement",false);

    gsInfo << "SameElement flag is: " 
       << (A.options().getSwitch("SameElement") ? "ON" : "OFF") << "\n";

    gsInfo<<"Fine basis:\n";
    gsInfo<<dbasis_fine.basis(0)<<"\n";
    gsInfo<<"Coarse basis:\n";
    gsInfo<<dbasis_coarse.basis(0)<<"\n";
    
    real_t error_4 = computeNorm(dbasis_tensor, dbasis_fine, mp, *Ccrs_, Cfine, A.options());
    real_t error_5 = computeNorm(dbasis_fine,   dbasis_fine, mp, *Ccrs_, Cfine, A.options());
    real_t error_6 = computeNorm(dbasis_coarse, dbasis_fine, mp, *Ccrs_, Cfine, A.options());

    gsInfo<< "Error with integration on TPr basis (SE=OFF): "<< error_4<<"\n";
    gsInfo<< "Error with integration on fine basis (SE=OFF): "<< error_5<<"\n";
    gsInfo<< "Error with integration on crs basis (SE=OFF): "<< error_6<<"\n";

    // gsWriteParaview(dbasis_fine.basis(0), "basis_hugo_fine");
    // gsWriteParaview(dbasis_coarse.basis(0), "basis_hugo_crs");

    gsMesh<> mesh_fine(dbasis_fine.basis(0));
    gsMesh<> mesh_crs(dbasis_coarse.basis(0));

    gsWriteParaview(mesh_fine, "mesh_hugo_fine");
    gsWriteParaview(mesh_crs, "mesh_hugo_crs");

}// end main
