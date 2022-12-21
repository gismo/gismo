#/** @file gsMappedSpline_example.cpp

    @brief Example using the gsMappedBasis and gsMappedSpline class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, P. Weinmueller, H. Verhelst
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]


template <class T>
class gsSingleBasis : public gismo::gsFunction<T>
{

protected:
    gsBasis<T> & _basis;
    mutable gsMapData<T> _tmp;
    index_t m_bfID;


public:
    /// Shared pointer for gsSingleBasis
    typedef memory::shared_ptr< gsSingleBasis > Ptr;

    /// Unique pointer for gsSingleBasis
    typedef memory::unique_ptr< gsSingleBasis > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsSingleBasis(gsBasis<T> & basis, index_t bfID) :
            _basis(basis), m_bfID(bfID), _basis_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN;
    }

    ~gsSingleBasis() { delete _basis_piece; }

GISMO_CLONE_FUNCTION(gsSingleBasis)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 1;}

    mutable gsSingleBasis<T> * _basis_piece; // why do we need this?

    const gsFunction<T> & piece(const index_t k) const
    {
        //delete _basis_piece;
        _basis_piece = new gsSingleBasis(_basis, m_bfID);
        return *_basis_piece;
    }

    // Input is parametric coordinates of 1-D \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( targetDim() , u.cols() );
        result = _basis.evalSingle(m_bfID, u);
    }
};


void createSplineBasisL2Projection(gsMultiBasis<> & mb_level1, gsMultiBasis<> & mb_level2, gsSparseMatrix<> & cf)
{
    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);

    A.setIntegrationElements(mb_level2);

    // Set the discretization space
    auto u = A.getSpace(mb_level2);

    gsMatrix<> mat_lvl1, mat_lvl2;
    mat_lvl1.setZero(mb_level1.basis(0).size(),mb_level2.basis(0).size());
    //mat_lvl2.setIdentity(mb_level2.basis(0).size(),mb_level2.basis(0).size());
    for (index_t bfID = 0; bfID < mb_level1.basis(0).size(); bfID ++) {
        gsSingleBasis<real_t> sb(mb_level1.basis(0), bfID);
        auto aa = A.getCoeff(sb);

        gsBoundaryConditions<> bc_empty;
        u.setup(bc_empty, dirichlet::homogeneous, 0);
        A.initSystem();

        A.assemble(u * u.tr(), u * aa);

        gsSparseSolver<>::CGDiagonal solver;
        solver.compute(A.matrix());
        gsMatrix<> solVector = solver.solve(A.rhs());

        auto u_sol = A.getSolution(u, solVector);
        gsMatrix<> sol;
        u_sol.extract(sol);
        mat_lvl1.row(bfID) = sol.transpose();
    }

    mat_lvl1 = mat_lvl1.transpose();
    cf = mat_lvl1.sparseView(1,1e-10);
}



int main(int argc, char *argv[])
{
    //! [Parse command line]
    index_t choice = 0;

    bool plot = false;
    bool xml = false;

    index_t discreteDegree = 3;
    index_t discreteRegularity = 2;
    index_t numRefine = 1;

    index_t discreteDegree2 = 3;
    index_t discreteRegularity2 = 2;
    index_t numRefine2 = 1;

    index_t geoDim = 2;
    index_t parDim = 2;

    //std::string fn("msplines/spline_test");
    std::string fn;

    gsCmdLine cmd("Example using mapped spline.");
    cmd.addInt("c", "choice", "Which example/case/choice do you want to run?", choice);
#/** EXAMPLE 0: Create a mapped spline from the geometry "-c 0" (default options)

The example run with the (given) geometry and create a multi-basis from the geometry.
Then the coeficient matrix is constructed with the identity matrix.

User options:
    -f  <string>    The filepath of the geometry.
                    If it is empty, the geometry will be a unit interval/square/cube. To choose which one, call "-D":
    -D  <int>       Define the dimension of the geometry domain: Interval(1), Square(2), Cube(3)

    -p  <int>       Discrete polynomial degree. The basis functions will set to degree "p"
    -r  <int>       Discrete regularity. The basis functions will set to regularity "r"
    -l  <int>       Discrete refinement. The basis functions will be uniformed refined "l"-times

    --xml           If you want to save the basis functions in a XML file.
    --plot          Create a ParaView visualization file with the mapped basis
**/

#/** EXAMPLE 1: Read a mapped spline from the xml file "-c 1"

The example read the basis function, coefficient sparse matrix and coordinates from the xml file. Use the flag "--plot"
to visualize the basis functions and the resulting geometry.

User options:
    -f  <string>    The filepath of the xml file.

    --plot          Create a ParaView visualization file with the mapped basis
**/

#/** EXAMPLE 2: Create the old basis functions as linear combination of the new basis "-c 2"

The example run with the (given) geometry and create a multi-basis from the geometry. We call the "old" space of that
basis functions as V(p,r,l). Then we set up the new space W(P,R,L) and do an L2-projection to obtain
the linear combination to represent the basis functions from V(p,r,l) with the basis of W(P,R,L), e.g.

v_j = \sum_i a_{i,j} w_i

which gives us the coeficient matrix.

Note: the transformation makes only sense iff V \subset W

User options:
    -f  <string>    The filepath of the geometry.
                    If it is empty, the geometry will be a unit interval/square/cube. To choose which one, call "-D":
    -D  <int>       Define the dimension of the geometry domain: Interval(1), Square(2), Cube(3)

    -p  <int>       Discrete polynomial degree. The basis functions will set to degree "p"
    -r  <int>       Discrete regularity. The basis functions will set to regularity "r"
    -l  <int>       Discrete refinement. The basis functions will be uniformed refined "l"-times

    -P  <int>       Discrete polynomial degree for space 2. The basis functions will set to degree "P"
    -R  <int>       Discrete regularity for space 2. The basis functions will set to regularity "R"
    -L  <int>       Discrete refinement for space 2. The basis functions will be uniformed refined "L"-times

    --xml           If you want to save the new basis functions in a XML file.
    --plot          Create a ParaView visualization file with the mapped basis
**/


    cmd.addString( "f", "file", "Input XML file prefix", fn );

    cmd.addInt("D", "geoDim", "Define the dimension of the geometry domain: Interval, Square, Cube.", geoDim);
    cmd.addInt("d", "parDim", "Define the dimension of the parameter domain.", parDim);

    cmd.addInt("p", "degree", "Set the degree for the basis", discreteDegree);
    cmd.addInt("r", "regularity", "Set the regularity for the basis", discreteRegularity);
    cmd.addInt("l", "loop", "Set the uniform refinement for the basis", numRefine);

    cmd.addInt("P", "degree2", "Set the degree for the second basis", discreteDegree2);
    cmd.addInt("R", "regularity2", "Set the regularity for the second basis", discreteRegularity2);
    cmd.addInt("L", "loop2", "Set the uniform refinement for the second basis", numRefine2);

    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("xml", "Save the XML file with the (special) basis functions", xml);
    cmd.getValues(argc,argv);
    //! [Parse command line]

    //! [Initialization]
    gsFileData<> fd;
    gsMultiPatch<real_t> mp;
    gsMultiBasis<real_t> mb;
    gsMultiBasis<real_t> mb_level2;
    gsSparseMatrix<real_t> cf;
    gsMatrix<real_t> coefs;
    //! [Initialization]

    //! [Set up the example]
    switch (choice) {
        // Example 0
        case 0:
            if (!fn.empty()) // Read geometry from XML file
            {
                fd.read(fn + ".xml");
                fd.getFirst(mp);
            }
            else
            {
                if (geoDim == 1)
                    mp.addPatch(gsNurbsCreator<>::BSplineUnitInterval(1));
                else if (geoDim == 2)
                    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1, 1, 1) );
                else if (geoDim == 3)
                    mp.addPatch( gsNurbsCreator<>::BSplineCube(1, 1, 1, 1) );
                else
                    gsInfo << "No geometry with dimension " << geoDim << "\n";
                mp.computeTopology();
            }

            mb = gsMultiBasis<>(mp);

            // Set degree and refinement
            mb.setDegree(discreteDegree);
            for (index_t i = 0; i < numRefine; i++)
                mb.uniformRefine(1, discreteDegree - discreteRegularity);

            // Set Trafo Matrix to identity for the mappedBasis
            cf.resize(mb.size(), mb.size());
            cf.setIdentity();

            break;
        // Example 1
        case 1:
            if (!fn.empty()) // Read data from XML file
            {
                fd.read(fn + ".xml");
                fd.getFirst(coefs); // Not really necessary
                fd.getFirst(mb);
                fd.getFirst(cf);
            }
            else
            {
                // replace it with default path/given example
                gsInfo <<"Please add the filepath of the xml file!\n";
            }
            break;

        case 2:
            if (!fn.empty()) // Read geometry from XML file
            {
                fd.read(fn + ".xml");
                fd.getFirst(mp);
            }
            else
            {
                if (geoDim == 1)
                    mp.addPatch(gsNurbsCreator<>::BSplineUnitInterval(1));
                else if (geoDim == 2)
                    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1, 1, 1) );
                else if (geoDim == 3)
                    mp.addPatch( gsNurbsCreator<>::BSplineCube(1, 1, 1, 1) );
                else
                    gsInfo << "No geometry with dimension " << geoDim << "\n";
                mp.computeTopology();
            }

            mb = gsMultiBasis<>(mp);
            // Set degree and refinement for the old basis functions
            mb.setDegree(discreteDegree);
            for (index_t i = 0; i < numRefine; i++)
                mb.uniformRefine(1, discreteDegree - discreteRegularity);

            mb_level2 = gsMultiBasis<>(mp);
            // Set degree and refinement for the new basis functions
            mb_level2.setDegree(discreteDegree2);
            for (index_t i = 0; i < numRefine2; i++)
                mb_level2.uniformRefine(1, discreteDegree2 - discreteRegularity2);

            createSplineBasisL2Projection(mb,mb_level2,cf);

            mb = mb_level2; // For the new cf
            break;

        default:
            gsInfo << "The choice doesn't exist. Please choose another example! \n";
            break;
    }
    //! [Set up the example]

    //! [Setup the Mapped Basis]
    gsMappedBasis<1,real_t> mbasis1;
    gsMappedBasis<2,real_t> mbasis2;
    gsMappedBasis<3,real_t> mbasis3;

    if (mb.dim() == 1) {
        mbasis1.init(mb,cf);
        gsInfo << "The MappedBasis has " << mbasis1.size() << " basis functions for all patches! \n";
    }
    else if (mb.dim() == 2) {
        mbasis2.init(mb,cf);
        gsInfo << "The MappedBasis has " << mbasis2.size() << " basis functions for all patches! \n";
    }
    else if (mb.dim() == 3) {
        mbasis3.init(mb,cf);
        gsInfo << "The MappedBasis has " << mbasis3.size() << " basis functions for all patches! \n";
    }
    //! [Setup the Mapped Basis]

    //! [Some computation on the basis]
/*
    gsMatrix<> points(2, 4), result_mspline;
    points << 0, 0, 0.2, 0, 0.5, 0.5, 1, 1;
    result_mspline = mbasis.basis(0).eval(points);
    result_mspline = mbasis.basis(0).deriv(points);
    result_mspline = mbasis.basis(0).deriv2(points);

    index_t bfID = 0;
    result_mspline = mbasis.basis(0).evalSingle(bfID, points);
 */
    //! [Some computation on the basis]


    //! [Setup the Mapped Spline]
    if (choice!=1) // If it is not example 2
    {
        gsMatrix<real_t> supp;
        if (mb.dim() == 1)
            supp = mbasis1.basis(0).support();
        if (mb.dim() == 2)
            supp = mbasis2.basis(0).support();
        if (mb.dim() == 3)
            supp = mbasis3.basis(0).support();

        gsVector<real_t> a = supp.col(0);
        gsVector<real_t> b = supp.col(1);
        unsigned npts = 0;
        if (mb.dim() == 1)
            npts = mbasis1.size();
        if (mb.dim() == 2)
            npts = mbasis2.size();
        if (mb.dim() == 3)
            npts = mbasis3.size();
        // Assume that the sqrt of mbasis1.size() is solvable
        gsVector<unsigned> np = uniformSampleCount(a,b, npts );
        gsMatrix<real_t> pts = gsPointGrid(a,b,np);

        coefs = pts.transpose();
    }

    gsMappedSpline<1,real_t> mspline1;
    gsMappedSpline<2,real_t> mspline2;
    gsMappedSpline<3,real_t> mspline3;
    if (mb.dim() == 1)
        mspline1.init(mbasis1,coefs);
    if (mb.dim() == 2)
        mspline2.init(mbasis2,coefs);
    if (mb.dim() == 3)
        mspline3.init(mbasis3,coefs);
    //! [Setup the Mapped Spline]

    //! [Export mapped spline to xml file]
    if (xml)
    {
        fd.clear();
        fd << coefs;
        fd << mb;
        fd << cf;
        fd.save( "MappedSpline.xml");
        gsInfo << "The mapped spline is export to MappedSpline.xml\n";
    }
    //! [Export mapped spline to xml file]

    //! [Export visualization in ParaView]
    if (plot && mb.dim() == 1)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( mbasis1.basis(0), "MappedBasis", 1000, false);

        std::string fileName;
        std::string basename = "MappedBasisSingle";

        gsParaviewCollection collection(basename);
        for (index_t i = 0; i < mbasis1.basis(0).size(); i++)
        {
            fileName = basename + "_0_" + util::to_string(i);
            gsWriteParaview_basisFnct(i, mbasis1.basis(0), fileName, 1000);
            collection.addPart(fileName + ".vts", -1, i);
        }
        collection.save();

        gsWriteParaview<>( mspline1, "MappedSpline", 1000);
    }
    else if (plot && mb.dim() == 2)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( mbasis2.basis(0), "MappedBasis", 1000, false);

        std::string fileName;
        std::string basename = "MappedBasisSingle";

        gsParaviewCollection collection(basename);
        for (index_t i = 0; i < mbasis2.basis(0).size(); i++) {
            fileName = basename + "_0_" + util::to_string(i);
            gsWriteParaview_basisFnct(i, mbasis2.basis(0), fileName, 1000);
            collection.addPart(fileName + ".vts", -1, i);
        }
        collection.save();

        gsWriteParaview<>( mspline2, "MappedSpline", 1000);
    }
    else if (plot && mb.dim() == 3)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( mbasis3.basis(0), "MappedBasis", 1000, false);

        std::string fileName;
        std::string basename = "MappedBasisSingle";

        gsParaviewCollection collection(basename);
        for (index_t i = 0; i < mbasis3.basis(0).size(); i++)
        {
            fileName = basename + "_0_" + util::to_string(i);
            gsWriteParaview_basisFnct(i, mbasis3.basis(0), fileName, 1000);
            collection.addPart(fileName + ".vts", -1, i);
        }
        collection.save();

        gsWriteParaview<>( mspline3, "MappedSpline", 1000);
    }
    if (plot) {
        gsInfo << "The mapped basis is plotted to MappedBasis.pvd\n";
        gsInfo << "The single mapped basis is plotted to MappedBasisSingle.pvd\n";
        gsInfo << "The mapped spline is plotted to MappedSpline.pvd\n";
    }

    //! [Export visualization in ParaView]
}
