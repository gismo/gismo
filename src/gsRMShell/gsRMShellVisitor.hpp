/** @file gsRMShellVisitor.h

    @brief RM shell element visitor.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Xia, HS. Wang

    Date:   2020-12-23
*/

#include "gsRMShellVisitor.h"
#include"gsMyBase/gsMyBase.h"

namespace gismo
{

    /** \brief Visitor for the RM shell equation.
     *
     * Assembles the bilinear terms
     */

    // 0.计算弹性矩阵 D
    template <class T, bool paramCoef>
    inline void gsRMShellVisitor<T, paramCoef>::computeDMatrix()
        {
            double E;
            double v;
            double k = 5.0 / 6.0;
            E = m_material.E_modulus;
            v = m_material.poisson_ratio;

            localD.setZero(dof, dof);
            localD(0, 0) = 1;   localD(0, 1) = v;
            localD(1, 0) = v;   localD(1, 1) = 1;
            localD(3, 3) = 0.5 - 0.5 * v;
            localD(4, 4) = k * (0.5 - 0.5 * v);
            localD(5, 5) = localD(4, 4);

            localD *= E / (1 - v * v);
        }

    // 1.设定积分规则
    template <class T, bool paramCoef>
    void gsRMShellVisitor<T, paramCoef>::initialize(const gsBasis<T>& basis,
        const index_t patchIndex,
        const gsOptionList& options,
        gsQuadRule<T>& rule)
        {
            m_patchcount = patchIndex;
            // Grab right-hand side for current patch
            //rhs_ptr = &pde_ptr->rhs()->piece(patchIndex);

            // Setup Quadrature
            //rule = gsQuadrature::get(basis, options); // harmless slicing occurs here

			gsVector<index_t> numNodes(basis.dim());
			unsigned digits = 0;
			for (index_t i = 0; i < basis.dim(); ++i)
			{// 积分点个数
				numNodes[i] = 2;
			}
			rule = gsGaussRule<T>(numNodes, digits);

            // Set Geometry evaluation flags
            //map_data.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
            map_data.flags = NEED_VALUE | NEED_JACOBIAN | NEED_GRAD_TRANSFORM;

        }

    // 2.1.compute greville abscissa points
    template <class T, bool paramCoef>
    inline void gsRMShellVisitor<T, paramCoef>::compute_greville()
        {
            //Returns the anchor points that represent the members of the basis in result. 
            //There is exactly one anchor point per basis function. 每个基函数只有一个锚点
            // for a B spline basis these are the Greville abscissae.
            gsMatrix<> greviAll;
            multi_patch.basis(m_patchcount).anchors_into(greviAll);
            //gsDebug << "m_mp : \n" << m_mp << endl;
            //gsDebug << "greviAll:\n" << greviAll << endl;

            // greviAll 中存放的是所有 Greville 横标点的值，当基函数是2阶时，
            // Greville 横标点与控制点是一一对应的，但是此处算出来的是[0,1],

            greville_pt.setZero(2, actives.rows());
            for (int i = 0; i < actives.rows(); ++i)
            {
                greville_pt.col(i) = greviAll.col(actives(i, 0));
            }

            //gsDebug << "greville used : my_greville\n" <<  my_greville << endl;
        }

    // 2.2.计算greville点的法向量
    template <class T, bool paramCoef>
    inline void gsRMShellVisitor<T, paramCoef>::comput_norm(gsGenericGeometryEvaluator< T, 2, 1 >& geoEval,
        const gsMatrix<T>& Nodes,
        gsMatrix<T>& norm)
        {
            //geo.outerNormal_into(Nodes, norm);

            geoEval.evaluateAt(Nodes);
            int count = Nodes.cols();
            norm.setZero(3, count);
            for (int i = 0; i < count; ++i)
            {
                geoEval.normal(i, m_normal);
                m_normal.normalize(); // 单位化
                norm(0, i) = m_normal(0);
                norm(1, i) = m_normal(1);
                norm(2, i) = m_normal(2);
            }
            //gsDebug << "greville norm:\n" << norm << endl;
        }

    // 2.3.计算坐标转换矩阵 T
    template <class T, bool paramCoef>
    inline void gsRMShellVisitor<T, paramCoef>::comput_gauss_tnb(gsMatrix<T> const& quNodes)
        {
            m_t1.setZero(3, quNodes.cols());
            m_t2.setZero(3, quNodes.cols());
            m_Gn.setZero(3, quNodes.cols());
            //gsDebug << "control points: \n" <<  control_point_coor << endl;
            //loop over every Gauss point i
            //gsDebug << " m_AllDersData[1] :\n " << m_AllDersData[1] << endl;

            for (int i = 0; i < quNodes.cols(); ++i)
            {
                for (int j = 0; j < (numActive); ++j)
                {
                    m_t1(0, i) += basisData[1](2 * j, i) * localCp(j, 0); //x
                    m_t1(1, i) += basisData[1](2 * j, i) * localCp(j, 1); //y
                    m_t1(2, i) += basisData[1](2 * j, i) * localCp(j, 2); //z

                    m_t2(0, i) += basisData[1](2 * j + 1, i) * localCp(j, 0);
                    m_t2(1, i) += basisData[1](2 * j + 1, i) * localCp(j, 1);
                    m_t2(2, i) += basisData[1](2 * j + 1, i) * localCp(j, 2);
                }

                // cross product to get n
                gsVector3d<real_t> temp_t1 = m_t1.col(i);
                gsVector3d<real_t> temp_t2 = m_t2.col(i);
                m_Gn.col(i) = temp_t1.cross(temp_t2); // 向量叉积 cross product
                //m_Gn.col(i) = m_t1.col(i).cross(m_t2.col(i)); // 这么写会报错
                m_Gn.col(i).normalize();
                m_t1.col(i).normalize();

                gsVector3d<real_t> temp_Gn = m_Gn.col(i);
                m_t2.col(i) = temp_Gn.cross(temp_t1);
                m_t2.col(i).normalize();

                //gsDebug << "t1 = \n" << m_t1 << endl;
                //gsDebug << "t2 = \n" << m_t2 << endl;
                //gsDebug << "Gn = \n" << m_Gn << endl;
            }
        }

    // 2.Evaluate on element.
    template <class T, bool paramCoef>
    inline void gsRMShellVisitor<T, paramCoef>::evaluate(const gsBasis<T>& basis,
        const gsGeometry<T>& geo,
        const gsMatrix<T>& quNodes)
        {
            // 高斯积分点
            map_data.points = quNodes;
            //gsDebug << "quNodes:\n" << map_data.points << "\n";

            // Compute the active basis functions
            // Assumes actives are the same for all quadrature points on the elements
            // 获取单元的控制点编号，存放在 actives
            basis.active_into(map_data.points.col(0), actives);
            numActive = actives.rows();
            //gsDebug << "No of cps:\n" << actives << "\n";

            // Initialize geometry evaluator -- TODO: Initialize inside visitor
            //typename gsGeometry<T>::Evaluator geoEval(multi_patch[m_patchcount].evaluator(map_data.flags)); // 引发异常
            //typename gsGeometry<T>::Evaluator geoEval(geo.evaluator(map_data.flags));
            globalCp.setZero(basis.size(), basis.dim() + 1);
            gsGenericGeometryEvaluator<T, 2, 1> geoEval(geo, map_data.flags);
            globalCp = geoEval.geometry().coefs();

            // 计算 Greville 横标点
            compute_greville();
            //gsDebug << "greville_pt:\n" << greville_pt << "\n";

            // Compute norm on Greville points 计算 Greville 点的法向量
            comput_norm(geoEval, greville_pt, greville_n);

            // 计算高斯点的非零基函数值和他们的1阶导数，并写入basisData
            // Evaluate basis functions on element
            basis.evalAllDers_into(map_data.points, 1, basisData);

            // Compute image of Gauss nodes under geometry mapping as well as Jacobians
            // geo.computeMap(map_data);

            // Evaluate right-hand side at the geometry points paramCoef
            // specifies whether the right hand side function should be
            // evaluated in parametric(true) or physical (false)
            //rhs_ptr->eval_into((paramCoef ? map_data.points : map_data.values[0]), rhsVals);

            // 计算单元的控制点坐标
            localCp.setZero(numActive, 3);
            //gsDebug << "global control point:s \n" << globalCp << "\n";
            for (int i = 0; i < numActive; ++i)
            {
                localCp.row(i) = globalCp.row(actives(i, 0));
            }
            //gsDebug << "local control points:\n " << localCp << "\n";

            // 计算坐标转换矩阵T
            comput_gauss_tnb(map_data.points);

            /*计算均布载荷的时候，计算表面积要用到雅可比矩阵，为了避免重复计算，
            在计算单元刚度阵的时候先把雅可比矩阵计算出来，放到rhs中，最后再统一乘载荷值*/
            if (rhs_cols == 0)
            { // 均布载荷可能没有，也可能不止一个
                rhs_cols = 1;
            }
            // Initialize local matrix/rhs
            localMat.setZero(dof * numActive, dof * numActive);
            localRhs.setZero(dof * numActive, rhs_cols);

            // 辅助矩阵
            B_aW.setZero(dof, 3 * numActive);
            B_tW.setZero(dof, 3 * numActive);
            globalB.setZero(dof, dof * numActive);
            Jacob.setZero(3, 3);
        }

    // 3.1.计算雅可比矩阵
    template <class T, bool paramCoef>
    inline void gsRMShellVisitor<T, paramCoef>::computeJacob(gsMatrix<T>& bGrads,
        gsMatrix<T>& bVals, 
        gsMatrix<T>& J,
        const real_t z,
        const index_t k)
        {

            J.setZero(3, 3);
            real_t thick = m_material.thickness;

            for (index_t j = 0; j < numActive; ++j)
            {
                // z is gauss points in [-1,1]
                J(0, 0) += bGrads(2 * j, k) * (localCp(j, 0) + (thick / 2.) * z * greville_n(0, j));
                J(0, 1) += bGrads(2 * j, k) * (localCp(j, 1) + (thick / 2.) * z * greville_n(1, j));
                J(0, 2) += bGrads(2 * j, k) * (localCp(j, 2) + (thick / 2.) * z * greville_n(2, j));
                /*gsDebug << J(0,0) << "  " << J(0,1) << "  " << J(0,2) << endl;*/
                J(1, 0) += bGrads(2 * j + 1, k) * (localCp(j, 0) + (thick / 2.) * z * greville_n(0, j));
                J(1, 1) += bGrads(2 * j + 1, k) * (localCp(j, 1) + (thick / 2.) * z * greville_n(1, j));
                J(1, 2) += bGrads(2 * j + 1, k) * (localCp(j, 2) + (thick / 2.) * z * greville_n(2, j));
                /*gsDebug << J(1,0) << "  " << J(1,1) << "  " << J(1,2) << endl;*/
                J(2, 0) += (thick / 2.) * bVals(j, k) * greville_n(0, j);
                J(2, 1) += (thick / 2.) * bVals(j, k) * greville_n(1, j);
                J(2, 2) += (thick / 2.) * bVals(j, k) * greville_n(2, j);
                /*gsDebug << J(2,0) << "  " << J(2,1) << "  " << J(2,2) << endl;*/
            }
        }

    /// 3.2.Computes the membrane and flexural strain first derivatives at greville point \em k
    template <class T, bool paramCoef>
    inline void gsRMShellVisitor<T, paramCoef>::computeStrainB(gsMatrix<T>& bVals,
        const gsMatrix<T>& bGrads,
        const real_t t,
        const index_t k)
        {

            gsVector<T, 3> m_v, n_der;
            B_aW.setZero(dof, 3 * numActive);
            B_tW.setZero(dof, 3 * numActive);
            globalB.setZero(dof, dof * numActive);

            gsMatrix<T> J_Inv = Jacob.inverse();

            for (index_t j = 0; j < numActive; ++j)
            {
                gsMatrix<T> tempB_aw(3, 1);
                tempB_aw(0, 0) = bGrads(j * 2, k);
                tempB_aw(1, 0) = bGrads(j * 2 + 1, k);
                tempB_aw(2, 0) = 0;

                gsMatrix<T> temp_B = J_Inv * tempB_aw;

                B_aW(0, j * 3) = temp_B(0, 0);
                B_aW(1, j * 3 + 1) = temp_B(1, 0);
                B_aW(2, j * 3 + 2) = temp_B(2, 0);
                B_aW(3, j * 3) = B_aW(1, j * 3 + 1);
                B_aW(3, j * 3 + 1) = B_aW(0, j * 3);
                B_aW(4, j * 3 + 1) = B_aW(2, j * 3 + 2);
                B_aW(4, j * 3 + 2) = B_aW(1, j * 3 + 1);
                B_aW(5, j * 3) = B_aW(2, j * 3 + 2);
                B_aW(5, j * 3 + 2) = B_aW(0, j * 3);
            }

            gsMatrix<T> pgBar;
            pgBar.setZero(3, numActive);
            real_t thick = m_material.thickness;

            for (index_t j = 0; j < numActive; ++j)
            {
                gsMatrix<T> tempB_tw(3, 1);
                tempB_tw(0, 0) = (thick / 2.0) * t * bGrads(j * 2, k);
                tempB_tw(1, 0) = (thick / 2.0) * t * bGrads(j * 2 + 1, k);
                tempB_tw(2, 0) = (thick / 2.0) * bVals(j, k);

                gsMatrix<T> temRax = J_Inv * tempB_tw;
                pgBar(0, j) = temRax(0, 0);
                pgBar(1, j) = temRax(1, 0);
                pgBar(2, j) = temRax(2, 0);
            }

            for (index_t j = 0; j < numActive; ++j)
            {
                B_tW(0, 3 * j + 1) = (1.0) * pgBar(0, j) * greville_n(2, j);
                B_tW(0, 3 * j + 2) = (-1.0) * pgBar(0, j) * greville_n(1, j);
                B_tW(1, 3 * j) = (-1.0) * pgBar(1, j) * greville_n(2, j);
                B_tW(1, 3 * j + 2) = (1.0) * pgBar(1, j) * greville_n(0, j);
                B_tW(2, 3 * j) = (1.0) * pgBar(2, j) * greville_n(1, j);
                B_tW(2, 3 * j + 1) = (-1.0) * pgBar(2, j) * greville_n(0, j);
                B_tW(3, 3 * j) = (-1.0) * pgBar(0, j) * greville_n(2, j);
                B_tW(3, 3 * j + 1) = (1.0) * pgBar(1, j) * greville_n(2, j);
                B_tW(3, 3 * j + 2) = (1.0) * pgBar(0, j) * greville_n(0, j)
                    - pgBar(1, j) * greville_n(1, j);
                B_tW(4, 3 * j) = (1.0) * pgBar(1, j) * greville_n(1, j)
                    - pgBar(2, j) * greville_n(2, j);
                B_tW(4, 3 * j + 1) = (-1.0) * pgBar(1, j) * greville_n(0, j);
                B_tW(4, 3 * j + 2) = (1.0) * pgBar(2, j) * greville_n(0, j);
                B_tW(5, 3 * j) = (1.0) * pgBar(0, j) * greville_n(1, j);
                B_tW(5, 3 * j + 1) = (1.0) * pgBar(2, j) * greville_n(2, j)
                    - pgBar(0, j) * greville_n(0, j);
                B_tW(5, 3 * j + 2) = (-1.0) * pgBar(2, j) * greville_n(1, j);
            }

            // 把两部分B矩阵（与位移相关的部分和与转动相关的部分）组合起来
            for (int j = 0; j < numActive; ++j)
            {
                globalB(0, dof * j + 0) = B_aW(0, 3 * j + 0);
                globalB(0, dof * j + 4) = B_tW(0, 3 * j + 1);
                globalB(0, dof * j + 5) = B_tW(0, 3 * j + 2);

                globalB(1, dof * j + 1) = B_aW(1, 3 * j + 1);
                globalB(1, dof * j + 3) = B_tW(1, 3 * j + 0);
                globalB(1, dof * j + 5) = B_tW(1, 3 * j + 2);

                globalB(2, dof * j + 2) = B_aW(2, 3 * j + 2);
                globalB(2, dof * j + 3) = B_tW(2, 3 * j + 0);
                globalB(2, dof * j + 4) = B_tW(2, 3 * j + 1);

                globalB(3, dof * j + 0) = B_aW(3, 3 * j + 0);
                globalB(3, dof * j + 1) = B_aW(3, 3 * j + 1);
                globalB(3, dof * j + 3) = B_tW(3, 3 * j + 0);
                globalB(3, dof * j + 4) = B_tW(3, 3 * j + 1);
                globalB(3, dof * j + 5) = B_tW(3, 3 * j + 2);

                globalB(4, dof * j + 1) = B_aW(4, 3 * j + 1);
                globalB(4, dof * j + 2) = B_aW(4, 3 * j + 2);
                globalB(4, dof * j + 3) = B_tW(4, 3 * j + 0);
                globalB(4, dof * j + 4) = B_tW(4, 3 * j + 1);
                globalB(4, dof * j + 5) = B_tW(4, 3 * j + 2);

                globalB(5, dof * j + 0) = B_aW(5, 3 * j + 0);
                globalB(5, dof * j + 2) = B_aW(5, 3 * j + 2);
                globalB(5, dof * j + 3) = B_tW(5, 3 * j + 0);
                globalB(5, dof * j + 4) = B_tW(5, 3 * j + 1);
                globalB(5, dof * j + 5) = B_tW(5, 3 * j + 2);
            }

        }

    // 3.3.全局弹性阵
    template <class T, bool paramCoef>
    inline void gsRMShellVisitor<T, paramCoef>::comput_D_global(int const k,
        gsMatrix<T>& D_global)
        {
            D_global.setZero(dof, dof);

            gsMatrix<> Tran;
            Tran.setZero(dof, dof);

            double a1 = m_t1(0, k); double b1 = m_t1(1, k); double c1 = m_t1(2, k);
            double a2 = m_t2(0, k); double b2 = m_t2(1, k); double c2 = m_t2(2, k);
            double a3 = m_Gn(0, k); double b3 = m_Gn(1, k); double c3 = m_Gn(2, k);

            //refer to : https://en.wikipedia.org/wiki/Orthotropic_material 
            // Cook FEM Chinese version P237
            Tran(0, 0) = a1 * a1;    Tran(0, 1) = b1 * b1;    Tran(0, 2) = c1 * c1;
            Tran(0, 3) = a1 * b1;    Tran(0, 4) = b1 * c1;    Tran(0, 5) = c1 * a1;

            Tran(1, 0) = a2 * a2;    Tran(1, 1) = b2 * b2;    Tran(1, 2) = c2 * c2;
            Tran(1, 3) = a2 * b2;    Tran(1, 4) = b2 * c2;    Tran(1, 5) = c2 * a2;

            Tran(2, 0) = a3 * a3;    Tran(2, 1) = b3 * b3;    Tran(2, 2) = c3 * c3;
            Tran(2, 3) = a3 * b3;    Tran(2, 4) = b3 * c3;    Tran(2, 5) = c3 * a3;

            Tran(3, 0) = 2. * a1 * a2;         Tran(3, 1) = 2. * b1 * b2;         Tran(3, 2) = 2. * c1 * c2;
            Tran(3, 3) = a1 * b2 + a2 * b1;    Tran(3, 4) = b1 * c2 + b2 * c1;    Tran(3, 5) = c1 * a2 + c2 * a1;

            Tran(4, 0) = 2. * a2 * a3;         Tran(4, 1) = 2. * b2 * b3;         Tran(4, 2) = 2. * c2 * c3;
            Tran(4, 3) = a2 * b3 + a3 * b2;    Tran(4, 4) = b2 * c3 + b3 * c2;    Tran(4, 5) = c2 * a3 + c3 * a2;

            Tran(5, 0) = 2. * a1 * a3;         Tran(5, 1) = 2. * b1 * b3;         Tran(5, 2) = 2. * c1 * c3;
            Tran(5, 3) = a1 * b3 + a3 * b1;    Tran(5, 4) = b1 * c3 + b3 * c1;    Tran(5, 5) = c1 * a3 + c3 * a1;

            D_global = Tran.transpose() * localD * Tran;
        }

    // 3.4.修正局部刚度矩阵的操作，使其非奇异
    template <class T, bool paramCoef>
    inline void gsRMShellVisitor<T, paramCoef>::fixSingularity()
        {
            gsMatrix<T, 3, 3> M;
            real_t s = 1.0E-5;

            std::vector<real_t> MAX;
            MAX.assign(numActive, 0.0);

            for (int i = 0; i < numActive; ++i)
            {
                if (fabs(localMat(3 + dof * i, 3 + dof * i)) > fabs(localMat(4 + dof * i, 4 + dof * i)))
                {
                    MAX[i] = fabs(localMat(3 + dof * i, 3 + dof * i));
                }
                else
                {
                    MAX[i] = fabs(localMat(4 + dof * i, 4 + dof * i));
                }
                if (MAX[i] < fabs(localMat(5 + dof * i, 5 + dof * i)))
                {
                    MAX[i] = fabs(localMat(5 + dof * i, 5 + dof * i));
                }
                //gsDebug << "MAX[i] = " << MAX[i] << endl;        
            }

            // 对于平板单元，对应转角自由度的刚度阵为0，这里使用乘小数的方法，将K_AA^WW进行修正
            for (int j = 0; j < numActive; ++j)
            {

                M(0, 0) = s * MAX[j] * greville_n(0, j) * greville_n(0, j);
                M(0, 1) = s * MAX[j] * greville_n(0, j) * greville_n(1, j);
                M(0, 2) = s * MAX[j] * greville_n(0, j) * greville_n(2, j);
                M(1, 0) = s * MAX[j] * greville_n(1, j) * greville_n(0, j);
                M(1, 1) = s * MAX[j] * greville_n(1, j) * greville_n(1, j);
                M(1, 2) = s * MAX[j] * greville_n(1, j) * greville_n(2, j);
                M(2, 0) = s * MAX[j] * greville_n(2, j) * greville_n(0, j);
                M(2, 1) = s * MAX[j] * greville_n(2, j) * greville_n(1, j);
                M(2, 2) = s * MAX[j] * greville_n(2, j) * greville_n(2, j);

                //gsDebug << "greville_n:\n" <<greville_n << endl;
                //gsDebug << "M :\n"<<M << endl;

                localMat(3 + dof * j, 3 + dof * j) += M(0, 0);
                localMat(3 + dof * j, 4 + dof * j) += M(0, 1);
                localMat(3 + dof * j, 5 + dof * j) += M(0, 2);
                localMat(4 + dof * j, 3 + dof * j) += M(1, 0);
                localMat(4 + dof * j, 4 + dof * j) += M(1, 1);
                localMat(4 + dof * j, 5 + dof * j) += M(1, 2);
                localMat(5 + dof * j, 3 + dof * j) += M(2, 0);
                localMat(5 + dof * j, 4 + dof * j) += M(2, 1);
                localMat(5 + dof * j, 5 + dof * j) += M(2, 2);
            }
        }

    // 3.5.均布载荷计算准备
	template <class T, bool paramCoef>
    inline void gsRMShellVisitor<T, paramCoef>::comput_rhs_press(gsDistriLoads<T>& m_pressure,
        gsMatrix<T>& Jacob,
        gsMatrix<T>& bVals,
        real_t z_thin,
        index_t rhs_cols,
        index_t dots,
        index_t k)
    {
        for (index_t rr = 0; rr < rhs_cols; ++rr)
        {// 可能不止一个均布载荷
            if (m_pressure[rr].section == "SUR")
            {
                if (dots == 0)
                {
                    // cross product 
                    gsVector3d<real_t> temp_a1 = Jacob.row(0).transpose();
                    gsVector3d<real_t> temp_a2 = Jacob.row(1).transpose();
                    gsVector3d<real_t> temp_area;
                    temp_area.setZero();
                    temp_area = temp_a1.cross(temp_a2);
                    // 表面积
                    real_t mapArea = 0.0;
                    mapArea = temp_area.norm();

                    /* 当前rhsVals中存储的是 bVals * mapArea的值，
                    * 但是两者的维度不一致，不能直接相加，需要进行处理
                    * bVals:[numActive x (quWeights.rows())]
                    * rhsVals:[(numActive*dof) x rhs_cols]*/
                    for (index_t ll = 0; ll < numActive; ++ll)
                    {
                        gsVector<real_t> temp_bvals;
                        temp_bvals.setZero(dof, 1);
                        temp_bvals.setConstant(bVals(ll, k) * mapArea);
                        rhsVals.col(rr).segment(ll * dof, dof) += temp_bvals;
                        /*vec1.segment(pos,n) : the n coeffs in the range [pos : pos + n - 1]
                        关于矩阵和向量的操作详见
                        http://eigen.tuxfamily.org/dox/group__QuickRefPage.html*/
                    }
                }
            }
            else if (m_pressure[rr].section == "UP")
            {
            }
            else if (m_pressure[rr].section == "DOWN")
            {
            }
            else if (m_pressure[rr].section == "LEFT")
            {
            }
            else if (m_pressure[rr].section == "RIGHT")
            {
            }
            else
            {
                gsInfo << "\nPressure section setting Error!\n";
            }
        }

    }
    

    // 3.算单刚
    template <class T, bool paramCoef>
    inline void gsRMShellVisitor<T, paramCoef>::assemble(gsDomainIterator<T>&,
        gsVector<T> const& quWeights)
        {
            gsMatrix<T>& bVals = basisData[0];
            gsMatrix<T>& bGrads = basisData[1];

            double detJ = 0.0;

            // z方向使用的2点高斯积分
            real_t GaussZ[2] = { -0.577350269189625764509148780502, 0.577350269189625764509148780502 };

            for (index_t i = 0; i < 2; ++i)
            { // z 方向循环

                for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
                {
                    //Computer Jacob matrix 计算Jacobi矩阵
                    computeJacob(bGrads, bVals, Jacob, GaussZ[i], k);   // 使用Greville点计算的
                    detJ = Jacob.determinant();

                    //Computer B matrix    
                    computeStrainB(bVals, bGrads, GaussZ[i], k);

                    // Compute elastic matrix on gauss point with coordinate transformation              
                    comput_D_global(k, globalD);

                    // 如果有均布载荷
                    rhsVals.setZero(dof * numActive, rhs_cols);
					if (m_pressure.numLoads() != 0)
					{
						// 为均布载荷计算做准备
						comput_rhs_press(m_pressure, Jacob, bVals, GaussZ[i], rhs_cols, i, k);
					}
                    
                    //z 方向积分点的权为 1.0;
                    gsMatrix<T> temMat = quWeights[k] *
                        (globalB.transpose() * globalD * globalB) * detJ * 1.0;

                    localMat += temMat;
                    localRhs += rhsVals * quWeights[k];

                }
            }
            // 为了避免刚度阵奇异而进行的修正
            //Add rotational stiffness                
            fixSingularity();
    }

    // 4.集总刚
    template <class T, bool paramCoef>
    inline void gsRMShellVisitor<T, paramCoef>::localToGlobal(const index_t patchIndex,
        const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseSystem<T>& system)
        {
            // Map patch-local DoFs to global DoFs
            //system.mapColIndices(actives, patchIndex, actives);

            // Add contributions to the system matrix and right-hand side
            //system.push(localMat, localRhs, actives, eliminatedDofs.front(), 0, 0);

            //Assert eliminatedDofs.rows() == rowMap.boundarySize()
            // from system.push (gsSparseSystem.h 972-)
            //gsDebug << "rows of m_matrix is: " << system.matrix().rows() << "\n";
            //gsDebug << "cols of m_matrix is: " << system.matrix().cols() << "\n";
            for (index_t i = 0; i != numActive; ++i)
            {

                const int ii = actives(i);
                for (index_t j = 0; j < numActive; ++j)
                {
                    const int jj = actives(j);

                    for (index_t m = 0; m < dof; ++m)
                    {
                        for (index_t n = 0; n < dof; ++n)
                        {
                            system.matrix().coeffRef(dof * ii + m, dof * jj + n) +=
                                localMat(dof * i + m, dof * j + n);
                        }
                    }

                }

				for (index_t l = 0; l < dof; ++l)
				{
					for (index_t k = 0; k < rhs_cols; ++k)
					{
						system.rhs().coeffRef(dof * ii + l, k) +=
							localRhs(dof * i + l, k);
					}
				}
            }
    }

} // namespace gismo


