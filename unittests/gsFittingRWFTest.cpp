/** TODO */

#include "gismo_unittest.h"
#include <gsModeling/gsFittingRWF.h>

SUITE(gsFittingRWFTest)
{
    TEST(stable_errguided_blade)
    {
        index_t iter = 9;
        std::string fn("/ya/ya135/ya13515/x/Testing/blade.xml");
        index_t numURef = 6;
        index_t numLRef = 6;
        real_t alpha = 1e-10;
        real_t toll = 5e-5;
        real_t tolerance = 5e-5;
        bool condcheck = false;
        index_t dreg = 2;

        // get data
        gsFileData<> fd_in(fn);
        gsMatrix<> uv, xyz;
        fd_in.getId<gsMatrix<> >(0, uv );
        fd_in.getId<gsMatrix<> >(1, xyz);

        // build fitting basis
        gsKnotVector<> u_knots (0, 1, 4, 4);
        gsKnotVector<> v_knots (0, 1, 0, 4);
        gsTensorBSplineBasis<2> T_tbasis( u_knots, v_knots );
        T_tbasis.uniformRefine( (1<<numURef)-1 );
        gsFittingRWFErrorGuided<2,real_t> ref( uv, xyz, T_tbasis);

        //build lambda
        gsFileData<> lbd_in("fitting/01lambda_3D-04.xml");
        gsTensorBSpline<2> lambda;
        lbd_in.getId<gsTensorBSpline<2> >(0, lambda);
        lambda.uniformRefine((1<<numLRef)-1);

        // compute alpha
        real_t h_rel  = lambda.basis().getMaxCellLength();
        alpha = std::pow(alpha, h_rel);

        for(index_t i=0; i < iter; i++)
            ref.nextIteration(alpha, toll, tolerance, lambda, false, true, condcheck, dreg, false);

        CHECK_CLOSE(ref.maxPointError(), 4.97e-05, 1e-7);
    }

    TEST(stable_suppguided_1D)
    {
        index_t iter = 9;
        std::string fn("/ya/ya135/ya13515/x/Testing/Example1D/BertsCurve/paperfinal/01Constraints/example_v_curve.xml");
        index_t numURef = 0;
        index_t numLRef = 6;
        real_t alpha = 1e-10;
        real_t toll = 1e-3;
        real_t tolerance = 1e-3;
        bool condcheck = false;
        index_t lknots = 4;

        // get data
        gsFileData<> fd_in(fn);
        gsMatrix<> uv, xyz;
        fd_in.getId<gsMatrix<> >(0, uv );
        fd_in.getId<gsMatrix<> >(1, xyz);

        // build fitting basis
        gsKnotVector<> u_knots (0, 1, 32, 4);
        gsBSplineBasis<>basis( u_knots );
        basis.uniformRefine( (1<<numURef)-1 );
        gsFittingRWF<1,real_t> ref( uv, xyz, basis);

        //build lambda
        gsBSpline<> lambda;
        gsFileData<> lbd_in("fitting/01lambda_1D-02.xml");
        lbd_in.getId<gsBSpline<> >(0, lambda );
        lambda.uniformRefine(lknots);
        lambda.uniformRefine((1<<numLRef)-1);

        ref.findLambda(lambda,0,1e-7);

        // compute alpha
        // const std::vector<real_t> & errors = ref.pointWiseErrors();
        alpha = math::pow(alpha,lambda.basis().getMaxCellLength()/(lambda.basis().support().coeff(1)-lambda.basis().support().coeff(0)));

        for(index_t i=0; i < iter; i++)
            ref.nextIteration(alpha, toll, tolerance, lambda, false, condcheck, 2, false);

        real_t maxerr = ref.maxPointError();

        CHECK_CLOSE(maxerr,0.001393879226060754,1e-10);
    }

    TEST(stable_errguided_1D)
    {
        index_t iter = 9;
        std::string fn("/ya/ya135/ya13515/x/Testing/Example1D/BertsCurve/paperfinal/01Constraints/example_v_curve.xml");
        index_t numURef = 0;
        index_t numLRef = 0;
        real_t alpha = 1e-10;
        real_t toll = 1e-3;
        real_t tolerance = 1e-3;
        bool condcheck = false;
        index_t lknots = 4;

        // get data
        gsFileData<> fd_in(fn);
        gsMatrix<> uv, xyz;
        fd_in.getId<gsMatrix<> >(0, uv );
        fd_in.getId<gsMatrix<> >(1, xyz);

        // build fitting basis
        gsKnotVector<> u_knots (0, 1, 32, 4);
        gsBSplineBasis<>basis( u_knots );
        basis.uniformRefine( (1<<numURef)-1 );
        gsFittingRWFErrorGuided<1,real_t> ref( uv, xyz, basis);

        //build lambda
        gsTensorBSpline<1,real_t> lambda;
        gsFileData<> lbd_in("fitting/01lambda_1D-02.xml");
        lbd_in.getId<gsTensorBSpline<1,real_t> >(0, lambda );
        lambda.uniformRefine(lknots);
        lambda.uniformRefine((1<<numLRef)-1);

        // compute alpha
        // const std::vector<real_t> & errors = ref.pointWiseErrors();
        alpha = math::pow(alpha,lambda.basis().getMaxCellLength()/(lambda.basis().support().coeff(1)-lambda.basis().support().coeff(0)));

        for(index_t i=0; i < iter; i++)
            ref.nextIteration(alpha, toll, tolerance, lambda, false, true, condcheck, 2, false);

        CHECK_CLOSE(ref.maxPointError(), 2.46e-04, 1e-6);
    }
}


