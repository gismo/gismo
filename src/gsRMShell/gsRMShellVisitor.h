/** @file gsRMShellVisitor.h

    @brief RM shell element visitor.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Xia, HS. Wang

    Date:   2020-12-23
*/

#pragma once

#include"gsMyBase/gsMyBase.h"
#include <gsAssembler/gsQuadrature.h>
#include "gsCore/gsGeometryEvaluator.hpp"
#include "gsRMShellBoundary.hpp"
#include "gsRMShellBase.h"


namespace gismo
{

    /** \brief Visitor for the RM shell equation.
     *
     * Assembles the bilinear terms
     */

    template <class T, bool paramCoef = false>
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
            Material material,
            gsDistriLoads<T> & pressure_,
            int dofs,
            int ptInt = 0):
            multi_patch(patches),
            m_material(material),
            m_pressure(pressure_),
            dof(dofs),
            m_inte(ptInt)
		{
            rhs_cols = m_pressure.numLoads();
            computeDMatrix();
		}

        // 0.计算弹性矩阵 D
        inline void computeDMatrix();
        

        // 1.设定积分规则
        void initialize(const gsBasis<T>& basis,
            const index_t patchIndex,
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
            const index_t k );
        
        /// 3.2.Computes the membrane and flexural strain first derivatives at greville point \em k
        inline void computeStrainB(gsMatrix<T>& bVals,
            const gsMatrix<T>& bGrads,
            const real_t t,
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
            real_t  z_thin,
            index_t rhs_cols,
			index_t dots,
			index_t k);
		
        // 3.算单刚
        inline void assemble(gsDomainIterator<T>&,
            gsVector<T> const& quWeights);
       

        // 4.集总刚
        inline void localToGlobal(const index_t patchIndex,
            const std::vector<gsMatrix<T> >& eliminatedDofs,
            gsSparseSystem<T>& system);
       

    public:
        index_t dof;
        index_t m_patchcount;
        index_t m_inte;
        Material m_material;
        gsMatrix<real_t> globalCp;
        gsMatrix<real_t> localCp;
        gsMatrix<real_t> greville_pt;
        gsMatrix<real_t> greville_n;
        gsMatrix<real_t> m_t1;
        gsMatrix<real_t> m_t2;
        gsMatrix<real_t> m_Gn;
        gsVector<real_t> m_normal;
        gsMultiPatch<T> multi_patch;
        gsDistriLoads<T> m_pressure;
    
    public:
        gsMatrix<real_t> localD;
        gsMatrix<real_t> globalD;
        gsMatrix<real_t> B_aW;
        gsMatrix<real_t> B_tW;
        gsMatrix<real_t> globalB;
        gsMatrix<real_t> Jacob;

    protected:
        // Pointer to the pde data
        const gsPoissonPde<T>* pde_ptr;

    protected:
        // Basis values
        std::vector<gsMatrix<T> > basisData;
        gsMatrix<T>        physGrad;
        gsMatrix<index_t> actives;
        index_t numActive;
        index_t rhs_cols;
    protected:
        // Right hand side ptr for current patch
        const gsFunction<T>* rhs_ptr;

        // Local values of the right hand side
        gsMatrix<T> rhsVals;

    protected:
        // Local matrices
        gsMatrix<T> localMat;
        gsMatrix<T> localRhs;

        gsMapData<T> map_data;
    };


} // namespace gismo


