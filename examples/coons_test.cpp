/** @file gsView.cpp

    @brief Produce Paraview file output from XML input, fo Visualizing  G+Smo objects

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>


using namespace gismo;

namespace expr
{

template<class T> class mixedDer_expr;
template<class T>
class mixedDer_expr : public gismo::expr::_expr<mixedDer_expr<T> >
{
    typename gismo::expr::gsFeVariable<T>::Nested_t u;

public:
    enum {ScalarValued = 0, ColBlocks = 0, Space = 0};


    typedef T Scalar;
    mixedDer_expr(const gismo::expr::gsFeVariable<T> & _u)
    : u(_u)
    { }

    gsMatrix<T> eval(const index_t k) const
    {
        //u.eval(k);        
        GISMO_ASSERT(0!=u.cSize(),"divide by 0");

        // numActive x parDim
        gsMatrix<T> a = u.data().values[2]
            .reshapeCol(k, u.data().values[2].rows()/u.cSize(), u.cSize() )
            .bottomRows(1).transpose(); //for 2D only
        //gsDebugVar(a);
        return a;
    }

    index_t rows() const { return u.rows();   }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        u.parse(evList);
        u.data().flags |= NEED_VALUE|NEED_ACTIVE|NEED_DERIV2;
    }

    const gismo::expr::gsFeSpace<Scalar> & rowVar() const { return u.rowVar(); }
    const gismo::expr::gsFeSpace<Scalar> & colVar() const
    {return gismo::expr::gsNullExpr<T>();}
};

}//namespace expr
template<class T>
::expr::mixedDer_expr<T> mixedDer(const gismo::expr::gsFeVariable<T> & u)
{ return ::expr::mixedDer_expr<T>(u); }

#define HierSpline gsTHBSpline
#define HierSplineBasis gsTHBSplineBasis

// input is the coons patch computed for many levels of global refinement
HierSpline<2> hierPatchFromLevels(std::vector<gsGeometry<>::uPtr> & geo)
{
    // make ad hoc hierarchical basis
    HierSplineBasis<2> thb(geo.front()->basis());
    std::vector<index_t> box(5);
    index_t l = 1;
    box[0]=l;
    box[1]=box[2]=1<<l;
    box[3]=box[4]=2<<l;
    thb.refineElements(box);
    box[0]= l= 2;
    box[1]=1<<l; box[2]=0;
    box[3]=2<<l; box[4]=1<<l;
    thb.refineElements(box);
    box[0]= l= 3;
    box[1]=0    ; box[2]=1<<l;
    box[3]=1<<l ; box[4]=2<<l;
    thb.refineElements(box);
    gsInfo << thb <<"\n";
    //gsWrite(thb, "result_thb");

    gsMatrix<> hcf(thb.size(),3);
    for(index_t i = 0; i < thb.size(); i++)
    {
        l = thb.levelOf(i);
        index_t fti = thb.flatTensorIndexOf(i);
        hcf.row(i) = geo[l]->coef(fti);
    }
            
    return HierSpline<2>(thb,hcf);
}


// A Hier. patch from levels with (ad hoc) full boundary refinement.
// input is the coons patch computed for many levels of global refinement
HierSpline<2> hierPatchFromLevels2(std::vector<gsGeometry<>::uPtr> & geo)
{
    HierSplineBasis<2> thb(geo.front()->basis());
    std::vector<index_t> box(20);
    index_t l = geo.size()-1;
    index_t u1 = thb.tensorLevel(l).component(0).size() - thb.degree(0);
    index_t v1 = thb.tensorLevel(l).component(1).size() - thb.degree(1);
    box[0]=box[5]=box[10]=box[15]=l;
    box[1]=box[2]=0;//
    box[3]=1; box[4]=v1;
    box[6]=box[7]=0;//
    box[8]=u1; box[9]=1;
    box[11]=u1-2; box[12]=0;//
    box[13]=u1-1; box[14]=v1-1;
    box[16]=0; box[17]=v1-2;//
    box[18]=u1-1; box[19]=v1-1;
    thb.refineElements(box);

    gsMatrix<> hcf(thb.size(),3);
    for(index_t i = 0; i < thb.size(); i++)
    {
        l = thb.levelOf(i);
        index_t fti = thb.flatTensorIndexOf(i);
        hcf.row(i) = geo[l]->coef(fti);
    }
            
    return HierSpline<2>(thb,hcf);
}


// Computes an H-Coons patch adaptively
// adaptivitity is driven from the exact TP solution
HierSpline<2> hierPatchAdaptive(std::vector<gsGeometry<>::uPtr> & geo,
                                HierSplineBasis<2> & thb)
{
    index_t refCriterion = PUCA, refExt = 1;
    real_t refParameter = 0.85;

    //First compute the h-coons on the input
    gsMatrix<> hcf(thb.size(),3);
    for(index_t i = 0; i < thb.size(); i++)
    {
        index_t l = thb.levelOf(i);
        index_t fti = thb.flatTensorIndexOf(i);
        hcf.row(i) = geo[l]->coef(fti);
    }
    HierSpline<2> hpatch(thb,hcf);
    gsMultiPatch<> mp(hpatch);
    gsMultiBasis<> mb;
    
    typedef gsExprEvaluator<real_t>::variable    variable;
    gsExprEvaluator<real_t> ev;
    //gsMultiBasis<> mb(geo.back()->basis());
    //ev.setIntegrationElements(mb);
    variable v = ev.getVariable(mp);
    variable w = ev.getVariable(*geo.back());

    index_t RefineLoopMax = geo.size();
    std::vector<bool> elMarked;
    // So, ready to start the adaptive refinement loop:
    gsInfo << "Basis: "<< thb <<"\n";

    mb = gsMultiBasis<>(geo.back()->basis());
    ev.setIntegrationElements(mb);
            
    for( int RefineLoop = 1; RefineLoop < RefineLoopMax ; RefineLoop++ )
    {
        gsInfo << "\n ====== Iteration " << RefineLoop << " of " << RefineLoopMax-1 << " ======" << "\n" << "\n";

        // mb = gsMultiBasis<>(mp);
        // ev.setIntegrationElements(mb);

        // The vector with element-wise local error estimates.
        ev.integralElWise( ( mixedDer(v) - mixedDer(w) ).sqNorm() );
        const std::vector<real_t> & elErrEst = ev.elementwise();

        // Mark elements for refinement, based on the computed local errors and
        // refCriterion and refParameter.
        //elMarked.resize( elErrEst.size() );
        gsMarkElementsForRef( elErrEst, refCriterion, refParameter, elMarked);
        gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true)<<"\n";
        gsInfo<<"J before: "<< ev.integral( mixedDer(v).sqNorm() ) <<"\n";
        gsRefineMarkedElements(mp, elMarked, refExt);

        // The integration mesh is changed
        // mb = gsMultiBasis<>(mp);
        // ev.setIntegrationElements(mb);

        gsInfo<<"J refined: "<< ev.integral( mixedDer(v).sqNorm() ) <<"\n";
        gsMatrix<> & cf = mp.patch(0).coefs();
        const HierSplineBasis<2> & bb =
            static_cast<const HierSplineBasis<2>&>(mp.patch(0).basis());
        gsInfo << "Basis: "<< bb <<"\n";
        for(index_t i = 0; i < bb.size(); i++)
        {
            index_t l = bb.levelOf(i);
            index_t fti = bb.flatTensorIndexOf(i);
            GISMO_ENSURE(l<(index_t)geo.size(),"level "<<l<<" too high");
            cf.row(i) = geo[l]->coef(fti);
        }
        
        gsInfo<<"J updated: "<< ev.integral( mixedDer(v).sqNorm() ) <<"\n";
        
        //if (dump)
        //{
        //    std::stringstream ss;
        //    ss << "adaptive_step_" << RefineLoop << ".xml";
        //    gsWrite(bases[0], ss.str());
        //// }
    }

    return static_cast<HierSpline<2>&>(mp.patch(0));
}


// test that the Coons patch has the permanence property
// by doing it on the coarse mesh and then refining to see that
// the result is actually the same
void testPermanence(std::vector<gsGeometry<>::uPtr> & geo)
{
    std::reverse(geo.begin(),geo.end());
    for (size_t i = 1; i<geo.size(); ++i)
    {
        for (size_t j = i; j<geo.size(); ++j)
            geo[j]->uniformRefine();
        
        gsInfo <<i<<". CF diff: "<< (geo.front()->coefs() - geo[i]->coefs()).norm() <<"\n";
    }
}

// returns the boundary curves of an example case
gsMultiPatch<> example01(index_t nref = 3)
{
    real_t tol = 1e-3;
    std::string fn = "curves3d/curve_boundary.xml";
    gsMultiPatch<> boundary;
    gsReadFile<>(fn, boundary);
    GISMO_ENSURE(!boundary.empty(), "The gsMultiPatch is empty - maybe file is missing or corrupt.");
    //gsInfo<<"Got "<< boundary <<"\n";
    boundary.computeTopology(tol);
    GISMO_ENSURE( boundary.isClosed(), "The boundary is not closed, adjust tolerance.");
    boundary.closeGaps(tol);
    
    for(index_t r = 0; r<nref; ++r)
        boundary.uniformRefine();
    return boundary;
}

// returns the boundary curves of an example case
gsMultiPatch<> example02(index_t nref = 3)
{
    gsMultiPatch<> res;
    
    gsBSplineBasis<> bsb(0,1,1,2);
    for(index_t r = 0; r<nref; ++r)
        bsb.uniformRefine();
    gsFunctionExpr<> patch("sin(2*pi*x)/3",1);
    gsMatrix<> tmp = patch.eval(bsb.anchors());
    auto scurve = bsb.interpolateAtAnchors(tmp);

    // refine but not fit
    //for(index_t r = 0; r<nref; ++r)
    //{
    //    bsb.uniformRefine();
    //    scurve->uniformRefine();
    // }
    
    gsTensorBSplineBasis<2> tbb (bsb.knots(), bsb.knots());
    gsMatrix<> cf = tbb.anchors().transpose();

    gsTensorBSpline<2> bs(tbb,cf);
    bs.embed3d();
    bs.coefs().col(2).head(bsb.size()) = scurve->coefs();
    bs.coefs().col(2).tail(bsb.size()) = - scurve->coefs();
    //gsWrite(bs, "result_bs");

    for (index_t i = 1; i!= 5; ++i)
        res.addPatch( bs.boundary(i) );

    //hcf.bottomRows(1).array() += 1;
    
    return res;
}

int main(int argc, char *argv[])
{
    std::string fn("coons.xml");
    index_t nRef(3);
    bool plot_mesh = false;
    bool plot_net = false;
    bool plot_boundary = false;
    bool get_basis = false;
    bool get_mesh = false;
    bool get_geo = false;

    //! [Parse Command line]
    gsCmdLine cmd("Hi, give me a file (eg: .xml) and I will try to draw it!");

    cmd.addSwitch("geometry", "Try to find and plot a geometry contained in the file", get_geo);
    cmd.addSwitch("mesh"    , "Try to find and plot a mesh contained in the file", get_mesh);
    cmd.addSwitch("basis"   , "Try to find and plot a basis contained in the file", get_basis);
    cmd.addInt   ("n", "nRef", "Number of samples to use for viewing", nRef);
    cmd.addSwitch("element"   , "Plot the element mesh (when applicable)", plot_mesh);
    cmd.addSwitch("controlNet", "Plot the control net (when applicable)", plot_net);
    cmd.addSwitch("boundary"  , "Plot the boundaries and interfaces of patches with colors", plot_boundary);
    cmd.addPlainString("filename", "File containing data to draw (.xml or third-party)", fn);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse Command line]

    if ( fn.empty() )
    {
        gsInfo<< cmd.getMessage();
        gsInfo<<"\nType "<< argv[0]<< " -h, to get the list of command line options.\n";
        return 0;
    }
    
    // /*
    // Generate examples:    
    std::vector<gsMultiPatch<> > bd;
    std::vector<gsGeometry<>::uPtr> coons, spring;
    bd.resize(nRef);
//    gsFileData<> fdc, fds, fdh;
    for (index_t i = 0; i<nRef; ++i)
    {
        bd[i] = example01(i);
        gsCoonsPatch<real_t> coonsgen(bd[i]);
        coonsgen.compute();
        gsWrite(coonsgen.compute(), "coons_t_" + util::to_string(i));
        coons.push_back( coonsgen.result().clone() );
        gsSpringPatch<real_t> springgen(bd[i]);
        springgen.compute();
        spring.push_back( springgen.result().clone() );
        //fdc.add( coons[i]);
    }

    gsInfo << "Coons patch: "<< *coons.back() <<"\n";
        
    gsWrite(bd.back(), "boudnary_02");
        
    HierSpline<2> thbc = hierPatchFromLevels2(coons);
    gsWrite(thbc, "coons_thb_02");

    
    HierSpline<2> thbs = hierPatchFromLevels2(spring);
    gsWrite(thbs, "spring_thb_02");
    gsWrite(*spring.back(), "spring_t_02");

    // Add more coons guys
    bd.resize(3*nRef);
    for (index_t i = nRef; i<2*nRef; ++i)
    {
        bd[i] = example01(i);
        gsCoonsPatch<real_t> coonsgen(bd[i]);
        coonsgen.compute();
        gsWrite(coonsgen.compute(), "coons_t_" + util::to_string(i));
        coons.push_back( coonsgen.result().clone() );
    }
    HierSpline<2> thba = hierPatchAdaptive(coons,thbs.basis());
    gsWrite(thba, "adaptive");
        
    //*/

    // Initiate the expression evaluator
    gsExprEvaluator<real_t> ev;
    gsMultiBasis<> mb(coons.back()->basis());
    ev.setIntegrationElements(mb);
    typedef gsExprEvaluator<real_t>::variable    variable;
    variable    u = ev.getVariable(*coons.back());
    variable    z = ev.getVariable(*spring.back());
    variable    v = ev.getVariable(thbc);
    variable    s = ev.getVariable(thbs);
    //ev.options().setInt("quB", 4); //increase quadrature accuracy
    //ev.options().setReal("quA", 2); //increase quadrature accuracy
    gsInfo<< "  Compare difference (Coons): "<< ev.max((u-v).sqNorm()) <<"\n";
    gsInfo<< "  Compare difference (Lapl): "<< ev.max((z-s).sqNorm()) <<"\n";

    //CLEANUP needed?
    
    for (index_t i = 0; i<nRef; ++i)
    {
        variable    w = ev.getVariable(*coons[i]);
        gsInfo<<"sz: "<<coons[i]->basis().size()<< "  Coons variational: "<< ev.integral( mixedDer(w).sqNorm() ) <<"\n";
    }
    gsInfo<< "H-Coons variational: "<< ev.integral( mixedDer(v).sqNorm() ) <<"\n";


    // /*
    gsInfo<<"Coons:\n";
    testPermanence(coons);
    gsInfo<<"Spring:\n";
    testPermanence(spring);
    //*/
    return EXIT_SUCCESS;
}

/*
  Permanence principle.


  Q: What is the Coons patch of a hierarchically refined patch?
  Q: What is the mesh required to obtain this patch?
  Q: Is H-Coons the same as th transfinite curve patch ??
  
  Assume a compatible boundary splines (same degree rep., knot interval)
  In the finest representation, there is a Coons patch
  Can we define a similar coons patch with THB ?

  For a fixed spline boundary: The Coons patch of a adaptively refined patch
  is the same as the one of a globally refined patch.

  The discrete laplacian does not satisfy the permanence principle,
  but it only APPROXIMATES it (ie. it is NOT a minimal surface). When
  we refine the approximation becomes better.
  
  Test 1: permanence principle (THB) very important: HB fails miserably

  Test 2: boundary approximation (THB)

  - Compare the H-Coons to to the finest grid Coons' patch. Check
    minimization! we would like them to be comparable!

    // Adaptive refinement:
    input: 4 boundary curves, with many knots, compatible or not compatible
    1. Start from a coarse mesh with these 4 curves (h-mesh with boundary refinement?, or not)
    2. Use a local-adaptive criterion for refinement (split element and compute functional in the smaller elements...)

    Properties:
    permanence
    bilinear precision
    twist-minimizing
    
 */
