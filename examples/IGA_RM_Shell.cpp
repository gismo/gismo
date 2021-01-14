// IGA_RM_Shell.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
/** @file IGA_RM_Shell.cpp

    @brief Provides main process for Reissner-Mindlin shells.

    Author(s): Y. Xia, HS. Wang
    
    Date:   2020-12-23
*/
#include <iostream>
#include <gismo.h>
#include "gsRMShellAssembler.hpp"
//#include "gsRMShellBase.h"
//#include "gsNURBSinfo.hpp"
//#include "gsRMShellBoundary.hpp"

using namespace std;
using namespace gismo;

int main(int argc, char* argv[])
{
	// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    // 1.初始化 ================================================================================
	// /////////////////////////////////////////////////////////////////////////////
	// Input options
	int		numElevate	    = 0;    //p细化，升阶
	int		numHrefine		= 0;	//h细化，插点
	int		plot		    = 1;    // 是否输出paraview图形
	int		testCase	    = 6;    // 算例选择
	bool	nonlinear	    = false;    // 非线性
	int		boolPatchTest   = false;    // 分片试验
	int		memIntegPtAdd   = 0;    // 添加的积分点的个数(默认gauss积分点数=基函数阶次)
    int     result          = 0;    
	int		dof_node		= 6;	// 节点自由度

	Material m_material;
	gsStopwatch clock;	// 计时器

    string input("../../TestCase/egg.xml");   // 模型数据文件
	string bc_file("../../TestCase/boundary.txt");      // 边界条件设置文件
    string output("../../TestResult/result.txt");    // 结果输出文件

	// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    // 2.数据读入 ================================================================================
	// /////////////////////////////////////////////////////////////////////////////

    // 文件读取方案A.命令行参数，在项目中右键属性》调试》命令行参数，可以修改，其主要是为cmd端操作提供方便
    // 关于该命令行的使用详见 https://gismo.github.io/commandLineArg_example.html
    gsCmdLine cmd("This program provides R-M shell numerical examples, give me your command!");
    cmd.addInt("t",  "testcase",     "Choose a test case. " ,               testCase);
    cmd.addInt("r",  "hRefine",      "Number of h-refinement steps",        numHrefine);
    cmd.addInt("e",  "degreeElevation", "Number of degree elevation steps", numElevate);
    cmd.addInt("i",  "integration", "Number of added integration point",     memIntegPtAdd);
    cmd.addInt("p",  "plot",         "Plot result in ParaView format",      plot);
    cmd.addSwitch("n", "nonlinear", "Nonlinear elasticity (otherwise linear)", nonlinear); //boolean
	cmd.addString("g", "geometry",  "InputFile containing Geometry (.xml, .axl, .txt)", input);
	cmd.addString("m", "name",      "Output file name",                     output);
	cmd.addString("c", "boundary",  "boundary condition",                   bc_file);
    // Read the arguments and update with the inputs, if given.
	try { cmd.getValues(argc, argv); }
	catch (int rv) { return rv; }
    // Print out the version information
    cmd.printVersion();
    gsInfo << "\nPrinting command line arguments:\n\n\n"
        << "testCase:\t\t"      << testCase     << "\n"
        << "numHref:\t\t"       << numHrefine   << "\n"
        << "numElevate:\t\t"    << numElevate   << "\n"
        << "memIntegPtAdd:\t\t" << memIntegPtAdd << "\n"
        << "plot:\t\t\t"          << plot         << "\n"
        << "nonlinear:\t\t"     << numElevate   << "\n"
        << "input:\t\t\t"         << input           << "\n"
        << "name_output:\t\t"   << output       << "\n"
        << "bc_file:\t\t"       << bc_file           << "\n";
	
	// 文件读取方案B.交互指令
	cout << "\n本 IGA RM Shell 程序可供测试的算例如下：\n";
	cout << "1.\t" << "矩形方板(简支&集中载荷)\n" ;
	cout << "2.\t" << "Pinched Cylinder\n";
	cout << "3.\t" << "Pinched Hemisphere\n";
	cout << "4.\t" << "Scordelis_Lo_roof\n";
	cout << "\n请输入要测试的算例(No.): ";
	cin >> testCase;
	cout << endl;
	cout << "请输入 H 细化次数(suggest:1-5): ";
	cin >> numHrefine;
	cout << endl;

	switch (testCase)
	{
	case 1:
	{
		// 固支/简支 集中/均布载荷 平板
		input ="../../TestCase/1-Rectangle_Plate.xml";
		//bc_file = "../../TestCase/1-Rectangle_Plate_Simple_pressure.txt"; //简支-均布载荷
		bc_file = "../../TestCase/1-Rectangle_Plate_Simple.txt"; //简支-集中载荷
		//bc_file = "../../TestCase/1-Rectangle_Plate_Fixed.txt"; //固支-集中载荷
		m_material.thickness = 1.0e-3; // F=2.5；原始的
		m_material.E_modulus = 2e11;
		m_material.poisson_ratio = 0.3;
		numElevate = 1; //输入基函数为1阶,升阶1次为2阶
		// 简支-集中载荷 位移参考解 w = 0.0063336
	}
		break;
	case 2:
	{   // 受压圆柱
		input ="../../TestCase/2-Pinched_cylinder.xml";
		bc_file = "../../TestCase/2-Pinched_cylinder.txt";
		m_material.thickness = 3.0E-2;
		m_material.E_modulus = 3e10;
		m_material.poisson_ratio = 0.3;
		// 位移参考解 w = -1.82E-07
	}
		break;
	case 3:
	{   // 受压半球
		input = "../../TestCase/3-Pinched_hemisphere.xml";//1/4模型
		bc_file = "../../TestCase/3-Pinched_hemisphere.txt"; // 1/4
		m_material.thickness = 0.04;
		m_material.E_modulus = 6.825e7;
		m_material.poisson_ratio = 0.3;
		//位移参考解 x = 9.40E-02
	}
		break;
	case 4:
	{	//Scordelis_Lo_roof 均布载荷
		input = "../../TestCase/4-Scordelis_Lo_roof.xml";	// 1/4
		bc_file = "../../TestCase/4-Scordelis_Lo_roof.txt"; //1/4
		m_material.thickness = 0.25;
		m_material.E_modulus = 4.32e8;
		m_material.poisson_ratio = 0.0;
		// 位移参考解 w = -3.02E-01
	}
		break;

	default:
		gsInfo << "Read data file ERROR!\n";
		break;
	}

    // 判断文件是否正确读入
    // 关于文件输入输出的使用详见 https://gismo.github.io/inputOutput_example.html
    if (!gsFileManager::fileExists(input))
    {
        gsWarn << "The file cannot be found!\n";
        return EXIT_FAILURE;
    }
    gsInfo << "\nRead file \"" << input << "\"\n";

	// 判断读入的文件中是否含有几何模型
    gsFileData<> fileData(input);
	gsGeometry<>::uPtr pGeom;
	if (fileData.has< gsGeometry<> >())
	{
		/* bool getFirst(Object & result)	const
        * Returns the first object of this type found in the XML data.
        * Doesn't look for nested objects. Writes it into the parameter.
        * https://gismo.github.io/classgismo_1_1gsFileData.html#a1d636af87cefa69eab119ac82b0d4632 */
		pGeom = fileData.getFirst< gsGeometry<> >();
	}
	else
	{
		gsWarn << "Input file doesn't have a geometry.\n";
		return EXIT_FAILURE;
	}
	gsInfo << "\nThe file contains: \n" << *pGeom << "\n";
	

	// 前处理--模型细化
	gsMultiPatch<> mul_patch;
	gsReadFile<> read_file(input, mul_patch);
	// p-refine
	if (numElevate != 0)
	{
		mul_patch.degreeElevate(numElevate);
	}
	// h-refine
	for (int r = 0; r < numHrefine; ++r)
	{
		mul_patch.uniformRefine();
	}
	// 输出细化后的图形
	gsWriteParaview(mul_patch, "../../TestResult/multi_patch", 1000, true);

	gsMultiBasis<> basis(mul_patch);
	gsInfo << "\nPatches: " << mul_patch.nPatches() << ", degree: " << basis.minCwiseDegree() << "\n";
	gsInfo << basis.basis(0) << "\n";

	// 获取节点向量信息
	gsNURBSinfo<real_t> nurbs_info(mul_patch, basis);

	// 读入边界条件txt文件
	if (!gsFileManager::fileExists(bc_file))
	{
		gsWarn << "The boundary file cannot be found!\n";
		return EXIT_FAILURE;
	}
	gsInfo << "\nRead file \"" << bc_file << "\"\n";
	ifstream bc_stream;
	bc_stream.open(bc_file);

	gsRMShellBoundary<real_t> BoundaryCondition(bc_stream, 
		basis, nurbs_info, dof_node);

	// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    // 3.数值计算 ================================================================================
	// /////////////////////////////////////////////////////////////////////////////
	
    /* 3.1 边界条件设置*/
    // 方案一 https://gismo.github.io/poisson_example.html
	// Define source function
	//gsFunctionExpr<> f("((pi*1)^2 + (pi*2)^2)*sin(pi*x*1)*sin(pi*y*2)",
		//"((pi*3)^2 + (pi*4)^2)*sin(pi*x*3)*sin(pi*y*4)", 2);
	// For homogeneous term, we can use this (last argument is the dimension of the domain)
	//gsConstantFunction<> f(0.0, 0.0, 2);
	
	// Print out source function and solution
	//gsInfo << "Source function " << f << "\n";

	gsBoundaryConditions<> bcInfo;
	// Every patch with a boundary need to be specified. In this
	// there are in total 8 sides (two for each patch)
	// Dirichlet Boundary conditions
	// First argument is the patch number
	//bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &g);

	// Neumann Boundary conditions
	//gsFunctionExpr<> hSouth("0.0", "0.0", 2);
	//bcInfo.addCondition(0, boundary::south, condition_type::neumann, &hSouth);
	gsFunctionExpr<> rhs;
	
	// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	// 进行中 =========================================================================================
	// /////////////////////////////////////////////////////////////////////////
	////////////// Setup solver and solve //////////////
	// Initialize Solver
	// Setup method for handling Dirichlet boundaries, options:
	//
	// * elimination: Eliminate the Dirichlet DoFs from the linear system.
	//
	// * nitsche: Keep the Dirichlet DoFs and enforce the boundary
	//
	// condition weakly by a penalty term.
	// Setup method for handling patch interfaces, options:
	//
	// * glue:Glue patches together by merging DoFs across an interface into one.
	//   This only works for conforming interfaces
	//
	// * dg: Use discontinuous Galerkin-like coupling between adjacent patches.
	//       (This option might not be available yet)
	// 参考 https://gismo.github.io/poisson_example.html
	//! [Assemble]
	// 3.1计算数据准备
	gsRMShellAssembler<real_t> assembler(mul_patch, basis, 
		m_material, BoundaryCondition, bc_stream,rhs, bcInfo,
		dirichlet::none, iFace::none, dof_node, memIntegPtAdd);

	// 3.2装配刚度阵 载荷阵
	// Generate system matrix and load vector
	gsInfo << "\nAssembling...";
	clock.restart();
	assembler.assemble();
	gsInfo << "done.\n";
	gsInfo << "!! Time to assemble : " << clock.stop() << "\n";
	gsInfo << "Have assembled a system (matrix and load vector) with "
		<< assembler.numDofs()*dof_node << " dofs.\n";
	
	// 3.3求解
	// Initialize the conjugate gradient solver
	
 	//gsInfo << "Solving...\n";
 	//gsSparseSolver<>::CGDiagonal solver(assembler.matrix());
 	//gsMatrix<> solVector = solver.solve(assembler.rhs());
 	//gsInfo << "Solved the system with CG solver.\n";
	
	// 另一个求解器 https://gismo.github.io/linearSolvers_example.html
 	
 	gsSparseSolver<>::LU solverLU;
 	gsInfo << "\nEigen's LU: Started solving... ";
 	clock.restart();
 	solverLU.compute(assembler.matrix());
 	gsMatrix<> solVector = solverLU.solve(assembler.rhs());
 	gsInfo << "done.\n";
 	gsInfo << "!! Time to solve : " << clock.stop() << "\n";

	
    // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    // 4.结果输出 ================================================================================
    // /////////////////////////////////////////////////////////////////////////////
	// 输出求解结果
	// 输出解 
	ofstream of_res(output.c_str());
	for (int i = 0; i < assembler.numDofs(); ++i)
	{
		of_res << setw(5) << i;
		for (int j = 0; j < 6; ++j)
		{
			of_res << setw(16) << solVector(i * 6 + j, 0);
		}
		of_res << endl;
	}
	of_res.close();

	// Construct the solution as a scalar field
	gsMultiPatch<> mpsol;
	assembler.constructSolution(solVector, mpsol);
	gsMultiPatch<> deformation = mpsol;
	gsMatrix<> d_before;
	gsMatrix<> d_after;
	for (index_t k = 0; k < mul_patch.nPatches(); ++k)
	{
		deformation.patch(k).coefs() -= mul_patch.patch(k).coefs();
	}
	gsField<> sol(mpsol, deformation);

	// 绘制位移结果 Paraview 图形
	if (plot)
	{
		// Write approximate and exact solution to paraview files
		gsInfo << "Plotting in Paraview...\n";
		gsWriteParaview<>(sol, "../../TestResult/RM_Shell", 1000);
		gsInfo << "Wrote paraview files: .. / .. / TestResult / RM_Shell.pvd\n";
		// Run paraview
		//gsFileManager::open("../../TestResult/RM_Shell.pvd");
	}
	else
	{
		gsInfo << "Done. No output created, re-run with --plot"
			"to get a ParaView file containing the solution.\n";
	}
	
    // 结果可视化模型
    // 关于文件输入输出的使用详见 https://gismo.github.io/inputOutput_example.html
    if (output.empty())
	{
		gsInfo << "Call program with option -o <basename> to write data to files\n";
		gsInfo << "<basename>Paraview.vtp, <basename>Paraview.pvd, <basename>.xml\n";
		return EXIT_SUCCESS;
	}

    // 绘制几何 paraview 图像
	// writing a paraview file
	const std::string out = output + "Paraview";
	gsWriteParaview(*pGeom, out);
	gsInfo << "Wrote paraview files: " << out << ".vtp\n";
	gsInfo << "Wrote paraview files: " << out << ".pvd\n";
	
    /* 生成几何模型的xml文件*/
    // writing a G+Smo .xml file
    gsFileData<> fd;
	fd << *pGeom;
	// output is a string. The extension .xml is added automatically
	fd.save(output);
	gsInfo << "Wrote G+Smo file:     " << output << ".xml \n";


	// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    // 5.完结撒花 ===============================================================================
	// /////////////////////////////////////////////////////////////////////////////
    
    gsInfo << "\nEnd of RM Shell program!\n";
	return EXIT_SUCCESS;
    return 0;
 }

