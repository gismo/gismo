#pragma once
/** @file gsRMShellVisitor.h

    @brief RM shell element visitor.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Xia, HS. Wang

    Date:   2020-12-23
*/

#include <gsAssembler/gsQuadrature.h>
#include "gsGeometryEvaluator.hpp"
#include "gsRMShellBoundary.hpp"
#include "gsRMShellBase.h"

namespace gismo
{

    /** \brief Visitor for the RM shell equation.
     *
     * Assembles the bilinear terms
     */

    template <class T>
    class gsRMShellVisitor
    {
    public:

        /** \brief Constructor for gsRMShellVisitor.
         */
        gsRMShellVisitor(const gsPde<T>& pde)
        {
            pde_ptr = static_cast<const gsPoissonPde<T>*>(&pde);
        }
    public:
        // 构造函数
        gsRMShellVisitor(gsMultiPatch<T> const& patches,
            gsDofMapper         const& dofmapper,
            Material            material,
            gsDistriLoads<T>& pressure_,
            SSdata<T>& ssData,
            int                 dofs,
            int                 ptInt = 0) :
            multi_patch(patches),
            m_dofmap(dofmapper),
            m_material(material),
            m_pressure(pressure_),
            m_ssdata(ssData),
            dof(dofs),
            m_inte(ptInt)
        {
            rhs_cols = m_pressure.numLoads();
            computeDMatrix();
        }

        // 复制构造函数
        /*gsRMShellVisitor(gsRMShellVisitor<T>& C)
        {
            multi_patch = C.multi_patch;
            m_material = C.m_material;
            m_pressure = C.m_pressure;
            dof = C.dof;
            m_inte = C.m_inte;
            rhs_cols = C.rhs_cols;
        }*/

        ~gsRMShellVisitor()
        {

        }

        // 0.计算弹性矩阵 D
        inline void computeDMatrix();


        // 1.设定积分规则
        void initialize(const gsBasis<T>& basis,
            const index_t       patchIndex,
            const gsOptionList& options,
            gsQuadRule<T>& rule);


        // 2.1.compute greville abscissa points
        inline void compute_greville();


        // 2.2.计算greville点的法向量
        inline void comput_norm(gsGenericGeometryEvaluator< T, 2, 1 >& geoEval,
            const gsMatrix<T>& Nodes,
            gsMatrix<T>& norm);


        // 2.3.计算坐标转换矩阵 T
        inline void comput_gauss_tnb(gsMatrix<T> const& quNodes);

        // 2.Evaluate on element.
        inline void evaluate(const gsBasis<T>& basis,
            const gsGeometry<T>& geo,
            const gsMatrix<T>& quNodes);


        // 3.1.计算雅可比矩阵
        inline void computeJacob(gsMatrix<T>& bGrads,
            gsMatrix<T>& bVals,
            gsMatrix<T>& J,
            const real_t z,
            const index_t k);

        /// 3.2.Computes the membrane and flexural strain
        //  first derivatives at greville point \em k
        inline void computeStrainB(gsMatrix<T>& bVals,
            const gsMatrix<T>& bGrads,
            const real_t  t,
            const index_t k);

        // 3.3.全局弹性阵
        inline void comput_D_global(int const k,
            gsMatrix<T>& D_global);

        // 3.4.修正局部刚度矩阵的操作，使其非奇异
        inline void fixSingularity();

        // 3.5.均布载荷计算准备
        inline void comput_rhs_press(gsDistriLoads<T>& m_pressure,
            gsMatrix<T>& Jacob,
            gsMatrix<T>& bVals,
            index_t dots,
            index_t k,
            index_t rr);

        // 3.算单刚
        inline void assemble(gsDomainIterator<T>&,
            gsVector<T> const& quWeights);


        // 4.集总刚
        inline void localToGlobal(const index_t patchIndex,
            const std::vector<gsMatrix<T> >& eliminatedDofs,
            gsSparseSystem<T>& system);


    public:
        index_t          dof;           // 节点自由度
        index_t          m_patchcount;  // 片数
        index_t          m_inte;        // 增加积分点个数
        Material         m_material;    // 材料参数
        gsMatrix<real_t> globalCp;      // 全部控制点坐标
        gsMatrix<real_t> localCp;       // 单元对应的控制点坐标
        gsMatrix<real_t> greville_pt;   // 格雷维尔点
        gsMatrix<real_t> greville_n;    // 格雷维尔点处的法向量
        gsMatrix<real_t> m_t1;          // 坐标转换矩阵 t n b
        gsMatrix<real_t> m_t2;          // 坐标转换矩阵 t n b
        gsMatrix<real_t> m_Gn;          // 坐标转换矩阵 t n b
        gsVector<real_t> m_normal;      // 格雷维尔点处的法向量
        gsMultiPatch<T>  multi_patch;
        gsDofMapper      m_dofmap;
        gsDistriLoads<T> m_pressure;    // 均布载荷信息（实为局部变量）

    public:
        gsMatrix<real_t> localD;        // 局部坐标系下的D
        gsMatrix<real_t> globalD;       // 全局坐标系下的D
        gsMatrix<real_t> B_aW;          // 应变矩阵B
        gsMatrix<real_t> B_tW;          // 应变矩阵B
        gsMatrix<real_t> globalB;       // 单元的B阵
        gsMatrix<real_t> Jacob;         // 雅可比矩阵

        // 存储计算过程数据，以后计算应力应变，PS: SS 为 Stress Strain 的缩写
        SSdata<T>& m_ssdata;
        gsMatrix<real_t> SSDl;          // 用于计算应力应变的单元上的D
        gsMatrix<real_t> SSBl;          // 用于计算应力应变的单元上的B
        gsMatrix<real_t> SSNl;          // 用于计算应力应变的单元上的基函数N
        gsMatrix<index_t> SSPl;         // 用于计算应力应变的单元上的控制点编号P
    protected:
        // Pointer to the pde data
        const gsPoissonPde<T>* pde_ptr;

    protected:
        // Basis values
        std::vector<gsMatrix<T> >   basisData;  // 积分点的基函数值和基函数导数值
        std::vector < vector<gsMatrix<T> > >  basisData_p;// 积分点的基函数值和基函数导数值(均布载荷)
        gsMatrix<T>                 physGrad;   // 没有用到
        gsMatrix<index_t>           actives;    // 单元对应的控制点编号
        gsMatrix<index_t>           gloActives; // 单元对应的控制点编号的全局编号
        index_t                     numActive;  // 单元对应的控制点总数
        index_t                     rhs_cols;   // 均布载荷数
    protected:
        // Right hand side ptr for current patch
        const gsFunction<T>* rhs_ptr;    // 没有用到

        // Local values of the right hand side
        gsMatrix<T>                 rhsVals;    // 等于基函数值x表面积（均布载荷准备数据）

    protected:
        // Local matrices
        gsMatrix<T>                 localMat;   // 单元刚度阵
        gsMatrix<T>                 localRhs;   // 单元载荷阵

        gsMapData<T>                map_data;   // 积分规则信息
    };


} // namespace gismo


