/** @file gsParametrization_test.cpp

    @brief Test parametrization functionality.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
**/

#include "gismo_unittest.h"

SUITE(gsParametrization_test)
{
    TEST(standardParametrization)
    {
	gsFileData<real_t> fd_mesh("stl/norm.stl");
	gsMesh<real_t>::uPtr mesh = fd_mesh.getFirst<gsMesh<real_t> >();

	real_t eps = 1e-5;

	gsParametrization<real_t> param(*mesh);
	
	param.compute();
	gsMatrix<real_t> uv  = param.createUVmatrix();
	gsMatrix<real_t> xyz = param.createXYZmatrix();

	CHECK_CLOSE(0.5, uv(0, 0), eps);
	CHECK_CLOSE(0, uv(0, 1), eps);
	CHECK_CLOSE(1, uv(0, 2), eps);
	CHECK_CLOSE(1, uv(0, 3), eps);
	CHECK_CLOSE(0, uv(0, 4), eps);

	CHECK_CLOSE(0.5, uv(1, 0), eps);
	CHECK_CLOSE(0, uv(1, 1), eps);
	CHECK_CLOSE(0, uv(1, 2), eps);
	CHECK_CLOSE(1, uv(1, 3), eps);
	CHECK_CLOSE(1, uv(1, 4), eps);

	CHECK_CLOSE(0, xyz(0, 0), eps);
	CHECK_CLOSE(-1, xyz(0, 1), eps);
	CHECK_CLOSE(0, xyz(0, 2), eps);
	CHECK_CLOSE(1, xyz(0, 3), eps);
	CHECK_CLOSE(0, xyz(0, 4), eps);

	CHECK_CLOSE(0, xyz(1, 0), eps);
	CHECK_CLOSE(0, xyz(1, 1), eps);
	CHECK_CLOSE(-1, xyz(1, 2), eps);
	CHECK_CLOSE(0, xyz(1, 3), eps);
	CHECK_CLOSE(1, xyz(1, 4), eps);

	CHECK_CLOSE(1, xyz(2, 0), eps);
	CHECK_CLOSE(0, xyz(2, 1), eps);
	CHECK_CLOSE(0, xyz(2, 2), eps);
	CHECK_CLOSE(0, xyz(2, 3), eps);
	CHECK_CLOSE(0, xyz(2, 4), eps);
    }

    struct inputs
    {
	gsMatrix<real_t> verticesV0, paramsV0, verticesV1, paramsV1;

	gsOptionList options;

	const real_t eps;

	//gsMesh<real_t>::uPtr mesh;
	
	inputs() : eps(1e-5)
	{
	    gsFileData<real_t> fd_v0("parametrization/powerplant-bottom.xml");
	    fd_v0.template getId<gsMatrix<real_t> >(0, paramsV0);
	    fd_v0.template getId<gsMatrix<real_t> >(1, verticesV0);

	    gsFileData<real_t> fd_v1("parametrization/powerplant-top.xml");
	    fd_v1.template getId<gsMatrix<real_t> >(0, paramsV1);
	    fd_v1.template getId<gsMatrix<real_t> >(1, verticesV1);

	    // TODO: Uncomment the following two lines and the last
	    // test starts failing, although it uses its own filedata
	    // and mesh. Why?

	    // gsFileData<real_t> fd_mesh("parametrization/powerplant-mesh.stl");
	    // mesh = fd_mesh.getFirst<gsMesh<real_t> >();

	    options.addInt("parametrizationMethod", "parametrizationMethod", 1);
	}
    };

    TEST_FIXTURE(inputs, gsPeriodicParametrizationOverlap)
    {
	// TODO: When I put mesh into inputs, strange things start happening.
	gsFileData<real_t> fd_mesh1("parametrization/powerplant-mesh.stl");
	gsMesh<real_t>::uPtr mesh1 = fd_mesh1.getFirst<gsMesh<real_t> >();

	gsFileData<real_t> fd_over("parametrization/powerplant-overlap.stl");
	gsMesh<real_t>::uPtr over = fd_over.getFirst<gsMesh<real_t> >();

	// Construct the parametrization.
	gsPeriodicParametrizationOverlap<real_t> param(*mesh1,
						       verticesV0, paramsV0,
						       verticesV1, paramsV1,
						       *over, options);

	param.compute();
	gsMatrix<real_t> uv  = param.createUVmatrix();
	gsMatrix<real_t> xyz = param.createXYZmatrix();
	param.restrictMatrices(uv, xyz);

	// Compare with pre-computed values.
	CHECK_CLOSE(0.420668, uv(0, 0), eps);
	CHECK_CLOSE(0.242996, uv(1, 0), eps);
	CHECK_CLOSE(0.25, xyz(0, 0), eps);
	CHECK_CLOSE(-0.203645, xyz(1, 0), eps);
	CHECK_CLOSE(0.199575, xyz(2, 0), eps);
	
	CHECK_CLOSE(0.300928, uv(0, 10), eps);
	CHECK_CLOSE(0.497183, uv(1, 10), eps);
	CHECK_CLOSE(0.0359812, xyz(0, 10), eps);
	CHECK_CLOSE(-0.255694, xyz(1, 10), eps);
	CHECK_CLOSE(0.389171, xyz(2, 10), eps);

	CHECK_CLOSE(0.988428, uv(0, 20), eps);
	CHECK_CLOSE(0.497183, uv(1, 20), eps);
	CHECK_CLOSE(-0.25, xyz(0, 20), eps);
	CHECK_CLOSE(0.0646077, xyz(1, 20), eps);
	CHECK_CLOSE(0.389171, xyz(2, 20), eps);

	CHECK_CLOSE(0.675928, uv(0, 30), eps);
	CHECK_CLOSE(0.497183, uv(1, 30), eps);
	CHECK_CLOSE(0.155361, xyz(0, 30), eps);
	CHECK_CLOSE(0.206246, xyz(1, 30), eps);
	CHECK_CLOSE(0.389171, xyz(2, 30), eps);

	CHECK_CLOSE(0.110855, uv(0, 40), eps);
	CHECK_CLOSE(0.7691, uv(1, 40), eps);
	CHECK_CLOSE(-0.202487, xyz(0, 40), eps);
	CHECK_CLOSE(-0.164435, xyz(1, 40), eps);
	CHECK_CLOSE(0.578768, xyz(2, 40), eps);

	CHECK_CLOSE(0.286153, uv(0, 50), eps);
	CHECK_CLOSE(0, uv(1, 50), eps);
	CHECK_CLOSE(0.0970307, xyz(0, 50), eps);
	CHECK_CLOSE(-0.419029, xyz(1, 50), eps);
	CHECK_CLOSE(0, xyz(2, 50), eps);

	CHECK_CLOSE(0.661153, uv(0, 60), eps);
	CHECK_CLOSE(0, uv(1, 60), eps);
	CHECK_CLOSE(0.227687, xyz(0, 60), eps);
	CHECK_CLOSE(0.364909, xyz(1, 60), eps);
	CHECK_CLOSE(0, xyz(2, 60), eps);

	CHECK_CLOSE(0.0447457, uv(0, 70), eps);
	CHECK_CLOSE(1, uv(1, 70), eps);
	CHECK_CLOSE(-0.307507, xyz(0, 70), eps);
	CHECK_CLOSE(-0.089105, xyz(1, 70), eps);
	CHECK_CLOSE(0.75, xyz(2, 70), eps);
    }

    TEST_FIXTURE(inputs, gsPeriodicParametrizationStitch)
    {
	// TODO: If any of the following two gsFileDatas gets removed,
	// ./bin/unittests -R gsParametrization_test
	// leads to a fail in this test.
	gsFileData<real_t> fd_v0b("parametrization/powerplant-bottom.xml");
	gsFileData<real_t> fd_v1b("parametrization/powerplant-top.xml");

	gsFileData<real_t> fd_mesh2("parametrization/powerplant-mesh.stl");
	gsMesh<real_t>::uPtr mesh2 = fd_mesh2.getFirst<gsMesh<real_t> >();

	gsMatrix<real_t> stitch;
	gsFileData<real_t> fd_stitch("parametrization/powerplant-stitch.xml");
	fd_stitch.getFirst<gsMatrix<real_t> >(stitch);

	// Construct the parametrization.
	gsPeriodicParametrizationStitch<real_t> param(*mesh2,
						      verticesV0, paramsV0,
						      verticesV1, paramsV1,
						      stitch, options);

	param.compute();
	gsMatrix<real_t> uv  = param.createUVmatrix();
	gsMatrix<real_t> xyz = param.createXYZmatrix();

	CHECK_CLOSE(0.420668, uv(0, 0), eps);
	CHECK_CLOSE(0.242996, uv(1, 0), eps);
	CHECK_CLOSE(0.25, xyz(0, 0), eps);
	CHECK_CLOSE(-0.203645, xyz(1, 0), eps);
	CHECK_CLOSE(0.199575, xyz(2, 0), eps);

	CHECK_CLOSE(0.300928, uv(0, 10), eps);
	CHECK_CLOSE(0.497183, uv(1, 10), eps);
	CHECK_CLOSE(0.0359812, xyz(0, 10), eps);
	CHECK_CLOSE(-0.255694, xyz(1, 10), eps);
	CHECK_CLOSE(0.389171, xyz(2, 10), eps);

	CHECK_CLOSE(0.988428, uv(0, 20), eps);
	CHECK_CLOSE(0.497183, uv(1, 20), eps);
	CHECK_CLOSE(-0.25, xyz(0, 20), eps);
	CHECK_CLOSE(0.0646077, xyz(1, 20), eps);
	CHECK_CLOSE(0.389171, xyz(2, 20), eps);

	CHECK_CLOSE(0.675928, uv(0, 30), eps);
	CHECK_CLOSE(0.497183, uv(1, 30), eps);
	CHECK_CLOSE(0.155361, xyz(0, 30), eps);
	CHECK_CLOSE(0.206246, xyz(1, 30), eps);
	CHECK_CLOSE(0.389171, xyz(2, 30), eps);

	CHECK_CLOSE(0.110855, uv(0, 40), eps);
	CHECK_CLOSE(0.7691, uv(1, 40), eps);
	CHECK_CLOSE(-0.202487, xyz(0, 40), eps);
	CHECK_CLOSE(-0.164435, xyz(1, 40), eps);
	CHECK_CLOSE(0.578768, xyz(2, 40), eps);

	CHECK_CLOSE(0.286153, uv(0, 50), eps);
	CHECK_CLOSE(0, uv(1, 50), eps);
	CHECK_CLOSE(0.0970307, xyz(0, 50), eps);
	CHECK_CLOSE(-0.419029, xyz(1, 50), eps);
	CHECK_CLOSE(0, xyz(2, 50), eps);

	CHECK_CLOSE(0.661153, uv(0, 60), eps);
	CHECK_CLOSE(0, uv(1, 60), eps);
	CHECK_CLOSE(0.227687, xyz(0, 60), eps);
	CHECK_CLOSE(0.364909, xyz(1, 60), eps);
	CHECK_CLOSE(0, xyz(2, 60), eps);

	CHECK_CLOSE(0.0447457, uv(0, 70), eps);
	CHECK_CLOSE(1, uv(1, 70), eps);
	CHECK_CLOSE(-0.307507, xyz(0, 70), eps);
	CHECK_CLOSE(-0.089105, xyz(1, 70), eps);
	CHECK_CLOSE(0.75, xyz(2, 70), eps);
    }
}
