/** @file Reissner-Mindlin_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Yang Xia & Hugo Verhelst & HS. Wang
*/
#include <iostream>
#include <gismo.h>
#include "gsRMShell/gsRMShellAssembler.hpp"
#include "gsRMShell/gsPeridynamic.hpp"
#include "gsRMShell/gsHyperMeshOut.h"
#include "gsRMShell/gsUmfSolver.h"

using namespace std;
using namespace gismo;
using namespace DLUT::SAE::PERIDYNAMIC;

int LOAD_STEP = 1;
double RATIO_OF_HORIZON_MESHSIZE = 2.0;
double HORIZON = 3.0;
bool STRAINS = false; // 计算应力应变否？

int main(int argc, char* argv[])
{
	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	// 1.初始化 --------------------------------------------------------------------
	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	// Input options
	int		numElevate = 0;    //p 细化，升阶
	int		numHrefine = 0;	//h 细化，插点
	int		plot = 1;    // 是否输出paraview图形
	int		testCase = 6;    // 算例选择
	bool	nonlinear = false;    // 非线性
	int		boolPatchTest = false;    // 分片试验
	int		memIntegPtAdd = 0;    // 添加的积分点的个数(默认gauss积分点数=基函数阶次)
	int     result = 0;
	int		dof_node = 6;	// 节点自由度

	Material m_material;
	gsStopwatch clock;	// 计时器

	string input("../../TestCase/egg.xml");			// 模型数据文件
	string bc_file("../../TestCase/boundary.txt");  // 边界条件设置文件
	string output("../../TestResult/result.txt");   // 结果输出文件

	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	// 2.数据读入 -------------------------------------------------------------------
	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

	// 文件读取方案A.命令行参数，在项目中右键属性》调试》命令行参数，可以修改，其主要是为cmd端操作提供方便
	// 关于该命令行的使用详见 https://gismo.github.io/commandLineArg_example.html
	gsCmdLine cmd("This program provides R-M shell numerical examples, give me your command!");
	cmd.addInt("t", "testcase", "Choose a test case. ", testCase);
	cmd.addInt("r", "hRefine", "Number of h-refinement steps", numHrefine);
	cmd.addInt("e", "degreeElevation", "Number of degree elevation steps", numElevate);
	cmd.addInt("i", "integration", "Number of added integration point", memIntegPtAdd);
	cmd.addInt("p", "plot", "Plot result in ParaView format", plot);
	cmd.addSwitch("n", "nonlinear", "Nonlinear elasticity (otherwise linear)", nonlinear); //boolean
	cmd.addString("g", "geometry", "InputFile containing Geometry (.xml, .axl, .txt)", input);
	cmd.addString("m", "name", "Output file name", output);
	cmd.addString("c", "boundary", "boundary condition", bc_file);
	// Read the arguments and update with the inputs, if given.
	try { cmd.getValues(argc, argv); }
	catch (int rv) { return rv; }
	// Print out the version information
	cmd.printVersion();
	gsInfo << "\nPrinting command line arguments:\n\n\n"
		<< "testCase:\t\t" << testCase << "\n"
		<< "numHref:\t\t" << numHrefine << "\n"
		<< "numElevate:\t\t" << numElevate << "\n"
		<< "memIntegPtAdd:\t\t" << memIntegPtAdd << "\n"
		<< "plot:\t\t\t" << plot << "\n"
		<< "nonlinear:\t\t" << numElevate << "\n"
		<< "input:\t\t\t" << input << "\n"
		<< "name_output:\t\t" << output << "\n"
		<< "bc_file:\t\t" << bc_file << "\n";

	// 文件读取方案B.交互指令
	cout << "\n本 IGA RM Shell 程序可供测试的算例如下：\n";
	cout << "1.\t" << "矩形方板(简支&集中载荷)\n";
	cout << "2.\t" << "Scordelis_Lo_roof 均布载荷\n";
	cout << "3.\t" << "受压圆柱\n";
	cout << "4.\t" << "受压半球\n";
	cout << "5.\t" << "耦合用矩形板拉伸\n";
	cout << "6.\t" << "测试用任意算例（固支矩形板集中载荷）\n";
	cout << "7.\t" << "Patch Test\n";
	cout << "8.\t" << "平板拉伸\n";
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
		input = "../../TestCase/1-thin_square.xml";
		//bc_file = "../../TestCase/1-Rectangle_Plate_Simple_pressure.txt"; //简支-均布载荷
		bc_file = "../../TestCase/1-thin_square_SIDE_simple.txt"; //简支-集中载荷
		//bc_file = "../../TestCase/1-Rectangle_Plate_Fixed.txt"; //固支-集中载荷
		m_material.thickness = 1.0e-3; // F=2.5；原始的
		m_material.E_modulus = 2e11;
		m_material.poisson_ratio = 0.3;
		// 简支-集中载荷 位移参考解 w = 0.0063336
	}
	break;
	case 2:
	{
		//Scordelis_Lo_roof 均布载荷
		input = "../../TestCase/2-Scordelis_Lo_roof.xml";
		bc_file = "../../TestCase/2-Scordelis_Lo_roof.txt";
		m_material.thickness = 0.25;
		m_material.E_modulus = 4.32e8;
		m_material.poisson_ratio = 0.0;
		// 位移参考解 w = -3.02E-01
	}
	break;
	case 3:
	{
		// 受压圆柱
		input = "../../TestCase/3-pinched_cylinder.xml";//1/4模型
		bc_file = "../../TestCase/3-pinched_cylinder.txt"; // 1/4
		m_material.thickness = 3.0E-2;
		m_material.E_modulus = 3e10;
		m_material.poisson_ratio = 0.3;
		// 位移参考解 w = -1.82E-07
	}
	break;
	case 4:
	{
		// 受压半球
		input = "../../TestCase/4-Pinched_hemisphere.xml";	// 1/4
		bc_file = "../../TestCase/4-Pinched_hemisphere.txt"; //1/4
		m_material.thickness = 0.04;
		m_material.E_modulus = 6.825e7;
		m_material.poisson_ratio = 0.3;
		//位移参考解 x = 9.40E-02
	}
	break;
	case 5:
	{
		//耦合用矩形板拉伸
		input = "../../TestCase/5-Plate_square_L10.xml";	// 1/4
		bc_file = "../../TestCase/5-Plate_square_L10.txt"; //1/4
		m_material.thickness = 1.0;
		m_material.E_modulus = 1.0;
		m_material.poisson_ratio = 0.33;
		//位移参考解 v = -2.2152
	}
	break;
	case 6:
	{
		//测试用任意算例
		// 郑老师的6自由度PD壳论文中的平板算例
		input = "../../TestCase/6-Rectangle_Plate.xml";
		bc_file = "../../TestCase/6-Rectangle_Plate_Simple.txt";
		m_material.thickness = 0.5;
		m_material.E_modulus = 2e5;
		m_material.poisson_ratio = 0.3;
		//位移参考解0.5687(固支) 详见郑老师Micro beam bond的论文
	}
	break;
	case 7:
	{
		// patch test
		input = "../../TestCase/7-patch_test.xml";
		bc_file = "../../TestCase/7-patch_test.txt";
		m_material.thickness = 1.0;
		m_material.E_modulus = 1000;
		m_material.poisson_ratio = 0.3;
		boolPatchTest = 1;
		// numHrefine = 0;	// 不细化
	}
	break;
	case 8:
	{
		//孟师兄发表论文中的平板拉伸，带均布载荷
		input = "../../TestCase/8-square_stretch.xml";
		bc_file = "../../TestCase/8-square_stretch.txt";
		m_material.thickness = 1.0;
		m_material.E_modulus = 1.0;
		m_material.poisson_ratio = 0.33;
		//位移参考解 u = 10
	}
	break;
	case 9:
	{
		//检查刚度阵是否对称的
		input = "../../TestCase/square_six.xml";
		bc_file = "../../TestCase/square_six.txt";
		m_material.thickness = 1.0;
		m_material.E_modulus = 1.0;
		m_material.poisson_ratio = 0.33;
		//位移参考解 u = 10
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

	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	// 3.前处理 --------------------------------------------------------------------
	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	// 3.1 前处理--模型细化
	gsMultiPatch<> multipatch;
	gsReadFile<> read_file(input, multipatch);
	// p-refine
	if (numElevate != 0)
	{
		multipatch.degreeElevate(numElevate);
	}
	// h-refine
	for (int r = 0; r < numHrefine; ++r)
	{
		multipatch.uniformRefine();
	}
	// 输出细化后的图形
	//gsWriteParaview(multipatch, "../../TestResult/multi_patch", 1000, true);
	gsInfo << "Export paraview files: ../../TestResult/multi_patch.pvd\n";

	// 细化后的基函数信息
	gsMultiBasis<> multibasis(multipatch);
	gsInfo << "\nPatches: " << multipatch.nPatches() << "\n";
	gsInfo << "basis.basis(0): " << multibasis.basis(0) << "\n";
	gsInfo << "Patch 0, knot vector xi: \n" << multibasis[0].component(0).detail() << "\n";
	gsInfo << "Patch 0, knot vector eta:\n" << multibasis[0].component(1).detail() << "\n";

	// Map 信息
	gsDofMapper dofmapper(multibasis, multipatch.nPatches());
	multibasis.getMapper(true, dofmapper, true);

	// 3.2 获取节点向量,控制点坐标等信息
	gsNURBSinfo<real_t> nurbs_info(multipatch, multibasis, dofmapper);
	// gsNURBSinfo<real_t> nurbs_info(multipatch, multibasis); // 方案1 舍弃

	// 生成模型k文件
	// 方案1 可能不行
	/*index_t nofp = 0;
	gsGeometry<real_t>& Geo = multipatch[nofp];
	nurbs_info.OutputKfile(Geo, "../../TestResult/RMShell.k");*/
	// 方案2
	nurbs_info.ExportKfile("../../TestResult/RMShell.k");
	// 继续向k文件中添加材料参数等数据
	if (true)
	{
		// 显示指定app模式，防止已有数据被丢弃
		ofstream append("../../TestResult/RMShell.k", ofstream::app);
		// 输出材料参数
		append << "*MAT_ELASTIC\n";
		append << setw(10) << 1
			<< setw(10) << m_material.rho
			<< setw(10) << m_material.E_modulus
			<< setw(10) << m_material.poisson_ratio
			<< setw(10) << " " << setw(10) << " " << "\n";
		// 输出part信息
		append << "*PART\n";
		append << "auto1\n";
		append << setw(10) << 1
			<< setw(10) << 1
			<< setw(10) << 1 << "\n";
		// 输出壳的厚度
		append << "*SECTION_SHELL\n";
		append << setw(10) << 1
			<< setw(10) << 0 << "\n";
		append << setw(10) << m_material.thickness
			<< setw(10) << 1.0 << setw(10) << 1.0 << setw(10) << 1.0 << "\n";
		append.close();
	}

	// 3.3 读入边界条件txt文件
	if (!gsFileManager::fileExists(bc_file))
	{
		gsWarn << "The boundary file cannot be found!\n";
		return EXIT_FAILURE;
	}
	gsInfo << "\nRead boundary file \"" << bc_file << "\"\n";
	ifstream bc_stream;
	bc_stream.open(bc_file);

	// 3.4 将边界条件数据传入
	gsRMShellBoundary<real_t> BoundaryCondition(bc_stream,
		multibasis, nurbs_info, dof_node);

	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	// 4.数值计算 --------------------------------------------------------------------=
	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	// 参考 https://gismo.github.io/poisson_example.html

	// 4.1计算数据准备--等几何数据
	gsRMShellAssembler<real_t> assembler(multipatch, multibasis, dofmapper,
		nurbs_info, m_material, BoundaryCondition, bc_stream, dof_node, memIntegPtAdd);

	// 4.1计算数据准备--PD数据
	PDIGACouple<double> pdAss(&assembler, nurbs_info, m_material);

	// Generate system matrix and load vector
	// 4.2装配刚度阵 载荷阵
	// 4.2.1 装配 IGA 刚度阵 K
	gsInfo << "\nAssembling IGA K...";
	clock.restart();
	assembler.assemble();
	gsInfo << "done.\n";
	gsInfo << "=== >>> Time to assemble IGA global K : " << clock.stop() << "s.\n";
	gsInfo << "Have assembled a system (matrix and load vector) with "
		<< assembler.numDofs() * dof_node << " dofs.\n";

	// 4.2.2 装配 PD 刚度阵 K_PD
	//装配PD刚度阵
	gsInfo << "\nAssembling PD K...";
	clock.restart();
	pdAss.ImplicitAnalysis();
	gsInfo << "done.\n";
	gsInfo << "=== >>> Time to assemble PD global K : " << clock.stop() << "s.\n";

	// 4.2.3 合成耦合刚度阵
	gsInfo << "\nCoupling PD-IGA K... ";
	clock.restart();
	pdAss.combineStiffWithIGA();
	gsInfo << "done.\n";
	gsInfo << "=== >>> Time to Couple PD-IGA K : " << clock.stop() << "s.\n";

	// 4.2.4 添加边界条件
	gsInfo << "\nAdding boundary conditions... ";
	assembler.m_BoundaryCondition.setPressure(assembler.m_system);		// 均布载荷
	assembler.m_BoundaryCondition.setpLoad(assembler.m_system);			// 集中载荷
	assembler.m_BoundaryCondition.setDisp_constra(assembler.m_system);	// 位移约束
	gsInfo << "done.\n";

	// 4.2.5 分片试验
	if (boolPatchTest)
	{
		gsInfo << "\nPatch test ...";
		patch_test(assembler);
		gsInfo << "done.\n";
		return 0;
	}

	// 4.3求解
	/* Initialize the conjugate gradient solver */
	// 求解器1 CG 
	/*gsInfo << "Solving...\n";
	gsSparseSolver<>::CGDiagonal solver(assembler.matrix());
	gsMatrix<> solVector = solver.solve(assembler.rhs());
	gsInfo << "Solved the system with CG solver.\n";*/

	// 求解器2 LU
	/* 另一个求解器 https://gismo.github.io/linearSolvers_example.html
	gsSparseSolver<>::LU solverLU;
	gsInfo << "\nEigen's LU: Started solving... ";
	clock.restart();
	solverLU.compute(assembler.matrix());
	gsMatrix<> solVector = solverLU.solve(assembler.rhs());
	gsInfo << "done.\n";
	gsInfo << "=== >>> Time to solve : " << clock.stop() << "s.\n";*/

	// 求解器3 Umf
	/*Eigen::SparseMatrix<double> eiK;
	assembler.gsSparseToSparse(eiK);*/
	gsMatrix<> solVector;
	solVector.setZero(assembler.rhs().rows(), 1);
	gsInfo << "\nUmf Solver: Started solving... ";
	clock.restart();
	gsumf_solver(assembler.matrix(), solVector, assembler.rhs());
	gsInfo << "done.\n";
	gsInfo << "=== >>> Time to solve : " << clock.stop() << "s.\n";

	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	// 5.结果输出 --------------------------------------------------------------------=
	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	// 5.1 输出求解结果
	// 输出位移解 [u, v, w, thetax, thetay, thetaz]
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

	// 5.2 构造变形后的模型
	// Construct the solution as a scalar field
	gsMultiPatch<> mpsol; //gsMatrix<T>& coeffs = mpsol.patch(p).coefs();
	assembler.constructSolution(solVector, mpsol); // 提取出位移[u,v,w]，构造变形结果
	gsMultiPatch<> deformation = mpsol;
	for (index_t k = 0; k < multipatch.nPatches(); ++k)
	{
		deformation.patch(k).coefs() -= multipatch.patch(k).coefs();
	}
	gsField<> sol(mpsol, deformation);

	// 5.3 后处理 -- 结果可视化
	if (plot)
	{
		// 1. 绘制位移结果 Paraview 图形
		// Write approximate and exact solution to paraview files
		gsInfo << "\nPlotting in Paraview...\n";
		gsWriteParaview<>(sol, "../../TestResult/RM_Shell", 1000, true);
		gsInfo << "Export paraview files: \t.. / .. / TestResult / RM_Shell.pvd\n";
		// Run paraview
		//gsFileManager::open("../../TestResult/RM_Shell.pvd");

		// 2. 绘制位移应力应变结果 HyperView 图形
		gsInfo << "\nPlotting in HyperView...\n";
		gsMatrix<> solMat; // [n x 6]
		solMat.setZero(nurbs_info.nurbs_sumcp, dof_node);
		assembler.vector2Matrix(solVector, solMat);
		// 计算应力应变
		gsMatrix<> strainMat, stressMat;
		strainMat.setZero(nurbs_info.nurbs_sumcp, dof_node);
		stressMat.setZero(nurbs_info.nurbs_sumcp, dof_node);
		if (STRAINS)
		{
			assembler.StressStrain(solVector, strainMat, stressMat);
		}

		gsOutputHMRes("../../TestResult/RM_Shell.ascii",
			solMat, stressMat, strainMat);
		gsInfo << "Export HyperView files:\t../../TestResult/RM_Shell.ascii\n";
	}
	else
	{
		gsInfo << "Done. No output created, re-run with --plot"
			"to get a ParaView file containing the solution.\n";
	}

	// 5.3 结果可视化模型
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
	gsInfo << "Export paraview files: \t" << out << ".vtp\n";
	gsInfo << "Export paraview files: \t" << out << ".pvd\n";

	/* 生成几何模型的xml文件*/
	// writing a G+Smo .xml file
	gsFileData<> fd;
	fd << *pGeom;
	// output is a string. The extension .xml is added automatically
	fd.save(output);
	gsInfo << "Export G+Smo file:\t" << output << ".xml \n";


	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	// 6.完结撒花 --------------------------------------------------------------------
	// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	gsInfo << "\nEnd of RM Shell program!\n";
	//return EXIT_SUCCESS;
	return 0;
}



