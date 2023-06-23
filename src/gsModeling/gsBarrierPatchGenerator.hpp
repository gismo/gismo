/** @file gsBarrierPatchGenerator.hpp

    @brief Provides patch construction from boundary data by using barrier method. It is a
    reference implementation of the following paper. If you make use of the code or the
    idea/algorithm in your work, please cite our paper:
	Ji, Y., Yu, Y. Y., Wang, M. Y., & Zhu, C. G. (2021).
	Constructing high-quality planar NURBS parameterization for
	isogeometric analysis by adjustment control points and weights.
	Journal of Computational and Applied Mathematics, 396, 113615.
	(https://www.sciencedirect.com/science/article/pii/S0377042721002375)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji, H.M. Verhelst
*/

#pragma once

#include <gsModeling/gsBarrierCore.h>

namespace gismo
{
/**
\brief Computes a patch parametrization given a set of
boundary geometries.  Parametrization is not guaranteed to be
non-singular. Works for surface, volumes.

\tparam d dimension of physical domain
\tparam T Coefficient type

\param bRep Input boundary representation
\param initialMethod Specify initialization method
\param filename Name of input data file

\ingroup Modeling
*/

    template<short_t d, typename T>
    gsBarrierPatchGenerator<d, T>::gsBarrierPatchGenerator(const gsMultiPatch<T> &bRep,
                                         const std::string &filename)
            :
            m_bRep(bRep)
    {
        // Input data Sanity Check
        this->defaultOptions();

        m_filename = filename;
        m_filename = m_filename.substr(
                m_filename.find_last_of("/\\") + 1); // file name without a path
        m_filename = m_filename.substr(0, m_filename.find_last_of(
                ".\\")); // file name without an extension
        m_filename = m_filename.substr(0, m_filename.find_last_of("_\\"));

        gsMatrix<T> boundingBox;
        m_bRep.boundingBox(boundingBox);

        m_boundingBoxLeftBottomCorner = boundingBox.col(0).transpose();
        gsMatrix<T> boundingBoxRightTopCorner = boundingBox.col(1).transpose();

//        real_t modelSize = (boundingBoxRightTopCorner -
//                m_boundingBoxLeftBottomCorner).maxCoeff();
        for (short_t idim = 0; idim != d; ++idim)
        {
            m_scalingVec(idim) = m_boxsize /
                                 (boundingBoxRightTopCorner(idim) -
                                  m_boundingBoxLeftBottomCorner(idim));
//            m_scalingVec(idim) = m_boxsize / modelSize;
            m_scalingFactor *= m_scalingVec(idim);
        }
    }

    template<short_t d, typename T>
    void gsBarrierPatchGenerator<d, T>::setMultiPatch(const gsMultiPatch<T> &multiPatch)
    {
        m_mp = multiPatch;
        m_mb = gsMultiBasis<T>(m_mp);
        gsBarrierCore<d, T>::scaling(m_mp, -m_boundingBoxLeftBottomCorner,
                                     m_scalingVec);

        gsDofMapper newMapper;
        newMapper.init(m_mb, m_mp.targetDim());
        for (size_t iptch = 0; iptch != m_mp.nPatches(); ++iptch)
        {
            gsMatrix<index_t> idx = m_mp.basis(iptch).allBoundary();
            for (index_t idim = 0; idim != m_mp.targetDim(); ++idim)
            {
                for (index_t idof = 0; idof != idx.size(); ++idof)
                {
                    newMapper.eliminateDof(idx(idof), iptch, idim);
                }
            }
        }
        newMapper.finalize();

        gsInfo << "#Numb of free  variables is " << newMapper.freeSize()
               << "\n";
        gsInfo << "#Numb of fixed variables is " << newMapper.boundarySize()
               << "\n";
        gsInfo << "#Numb of total variables is " << newMapper.size()
               << "\n\n\n";

        m_mapper = newMapper;
    }

    template<short_t d, typename T>
    void gsBarrierPatchGenerator<d, T>::defaultOptions()
    {
        m_options = gsBarrierCore<d, T>::defaultOptions();
        m_options.addInt("InitialMethod",
                         "Initialization Method: 0 Coons' patch (default), 1 Spring patch, 2: Cross-Ap. patch",
                         0);

        // TODO: (legacy) remove plotting from the hpp file
        // Set parameters and options for plot the results
        m_options.addInt("plot.npts", "Number of sampling points for plotting",
                         3000);
        m_options.addSwitch("plot.mesh",
                            "If true, plot the element mesh of the parameterization",
                            false);
        m_options.addSwitch("plot.net",
                            "If true, plot the control net of the parameterization",
                            false);
    }

    template<short_t d, typename T>
    gsMultiPatch<T> gsBarrierPatchGenerator<d, T>::getMultiPatch() { return m_mp; }

    template<short_t d, typename T>
    void gsBarrierPatchGenerator<d, T>::initialization()
    {
        // TODO: add some other initialization methods to test the robustness:
        // 1. discrete Coons 2. Smoothness energy 3. Spring model etc.
        switch (m_options.askInt("InitialMethod", 0))
        {
            case 1:
            {
                gsInfo << "Using Coons' patch construction.\n";
                gsCoonsPatch<T> coons(m_bRep);
                gsInfo << "Created a " << coons.compute() << "\n";
                m_mp.addPatch(coons.result());
                break;
            }
            case 2:
            {
//                // Cross Approximation patch method
//                // Question: only works for 3D surfaces? i.e., for mapping: x: R^2 --> R^3
//
//                gsInfo << "Using cross approximation construction.\n";
//                gsCrossApPatch<T> cross(m_bRep);
//                gsInfo << "Created a " << cross.compute() << "\n";
//                //if (save) gsWrite(spring.result(), "result_patch");
//                m_mp.addPatch(cross.result());
                break;
            }
            case 3:
            {
                // construct a parameterization with the inner control points all equal to (0, 0)
                gsInfo << "Set all the inner control points to a same point.\n";
                gsCoonsPatch<T> coons(m_bRep);
                coons.compute();
                m_mp.addPatch(coons.result());

                makeMapper();
                gsVector<T> initialZeros(m_mapper.freeSize());
                initialZeros.setConstant(m_boxsize / 2.0);
                convert_gsFreeVec_to_mp(initialZeros, m_mapper, m_mp);
                gsInfo << "Created a same point Patch." << "\n";
                break;
            }
            case 4:
            {
                // Smoothness energy method
                // However, the results seems similar with Spring model method?

                break;
            }
            case 0:
            default:
                // Spring model method
                gsInfo << "Using spring patch construction.\n";
                gsSpringPatch<T> spring(m_bRep);
                gsInfo << "Created a " << spring.compute() << "\n";
                m_mp.addPatch(spring.result());
                break;
        }
        m_mp.computeTopology();
    }

    template<short_t d, typename T>
    void gsBarrierPatchGenerator<d, T>::makeMapper()
    {
        ////! [Make mapper for the design DoFs]
        // Now, we set all the inner control points as optimization variables
        // It is also possible to set only a part of them as optimization variables
        m_mapper.init(m_mb, m_mp.targetDim());
        for (size_t iptch = 0; iptch != m_mp.nPatches(); ++iptch)
        {
            gsMatrix<index_t> idx = m_mp.basis(iptch).allBoundary();
            for (index_t idim = 0; idim != m_mp.targetDim(); ++idim)
            {
                m_mapper.markBoundary(iptch, idx, idim);
            }
        }
        m_mapper.finalize();

        gsInfo << "#Numb of free  variables is " << m_mapper.freeSize() << "\n";
        gsInfo << "#Numb of fixed variables is " << m_mapper.boundarySize()
               << "\n";
        gsInfo << "#Numb of total variables is " << m_mapper.size() << "\n\n\n";
    }

    template<short_t d, typename T>
    const gsGeometry<T> &gsBarrierPatchGenerator<d, T>::compute()
    {
        m_clock.restart();
        // STEP 0: scale the computational domain to [0,1]^d for better numerical stability
        gsBarrierCore<d, T>::scaling(m_bRep, -m_boundingBoxLeftBottomCorner,
                                     m_scalingVec);

        // STEP 1: Initial guess construction
        initialization();

        // write initialization to .xml file
//        gsBarrierCore<d,T>::scalingUndo(m_mp,-m_boundingBoxLeftBottomCorner,m_scalingVec);
//        gsWrite(m_mp, m_filename + "_init");
//        gsBarrierCore<d,T>::scaling(m_mp,-m_boundingBoxLeftBottomCorner,m_scalingVec);

        m_mb = gsMultiBasis<T>(m_mp);
        makeMapper();
        m_timer += m_clock.stop();

        m_mp = gsBarrierCore<d, T>::compute(m_mp, m_mapper, m_options);

        // STEP 4: undo scaling to the size of the input data
        gsBarrierCore<d, T>::scalingUndo(m_mp, -m_boundingBoxLeftBottomCorner,
                                         m_scalingVec);

        // STEP 5: output results and print parameterization quality information
        outputResult();

        return m_mp.patch(0);
    }

    template<short_t d, typename T>
    void gsBarrierPatchGenerator<d, T>::outputResult() const
    {

        m_evaluator.setIntegrationElements(m_mb);
        geometryMap G = m_evaluator.getMap(m_mp);
        GISMO_ENSURE(m_mp.nPatches() == 1, "Does not yet work for multi-patch, "
                                           "but multi-patch has "
                << m_mp.nPatches() << " patches");

        // resulting mesh
        m_evaluator.options() = m_options;
        if (d == 2) m_evaluator.options().setInt("plot.npts", 2000);
        else if (d == 3) m_evaluator.options().setInt("plot.npts", 10000);

        // Scaled Jacobian metric
        T minScaledJacobian, maxScaledJacobian, integralScaledJacobian;
        if (d == 2)
        {
            auto metric_ScaledJacobian =
                    jac(G).det() / (jac(G)[0].norm() * jac(G)[1].norm());
            m_evaluator.writeParaview(metric_ScaledJacobian, G,
                                      m_filename + "_metric_ScaledJacobian");
            minScaledJacobian = m_evaluator.template min(metric_ScaledJacobian);
            maxScaledJacobian = m_evaluator.template max(metric_ScaledJacobian);
            integralScaledJacobian = m_evaluator.template integral(
                    metric_ScaledJacobian);
        } else if (d == 3)
        {
            auto metric_ScaledJacobian = jac(G).det() /
                                         (jac(G)[0].norm() * jac(G)[1].norm() *
                                          jac(G)[2].norm());
            m_evaluator.writeParaview(metric_ScaledJacobian, G,
                                      m_filename + "_metric_ScaledJacobian");
            minScaledJacobian = m_evaluator.template min(metric_ScaledJacobian);
            maxScaledJacobian = m_evaluator.template max(metric_ScaledJacobian);
            integralScaledJacobian = m_evaluator.template integral(
                    metric_ScaledJacobian);
        }
        gsInfo << "Parameterization quality info (reference, not strictly):\n";
        gsInfo << " Scaled Jacobian:    min.: " << minScaledJacobian
               << "   max: " << maxScaledJacobian
               << "     integral: " << integralScaledJacobian << "\n";

        // Uniformity metric
        T minUniformity, maxUniformity, integralUniformity;
        T areaval = gsBarrierCore<d, T>::computeArea(m_mp);

        gsConstantFunction<T> areaConstFunc(areaval, d);
        auto area = m_evaluator.getVariable(areaConstFunc);
        auto metric_Uniformity = pow((jac(G).det() - area.val()) / area.val(),
                                     2);
        m_evaluator.writeParaview(metric_Uniformity, G,
                                  m_filename + "_metric_Uniformity");
        minUniformity = m_evaluator.template min(metric_Uniformity);
        maxUniformity = m_evaluator.template max(metric_Uniformity);
        integralUniformity = m_evaluator.template integral(metric_Uniformity);
        gsInfo << " Uniformity:    min.: " << minUniformity << "    max: "
               << maxUniformity << "    integral: "
               << integralUniformity << "\n";

//        // Frobenius condition number metric
//        T minFrobCondNum, maxFrobCondNum, integralFrobCondNum;
//        auto metric_FrobCondNum = jac(G).norm() * jac(G).inv().norm();
//        m_evaluator.writeParaview(metric_FrobCondNum, G, m_filename + "_metric_FrobCondNum");
//        minFrobCondNum = m_evaluator.template min(metric_FrobCondNum);
//        maxFrobCondNum = m_evaluator.template max(metric_FrobCondNum);
//        integralFrobCondNum = m_evaluator.template integral(metric_FrobCondNum);
//        gsInfo << " Frobenius condition number:    min.: " << minFrobCondNum << "   max: " << maxFrobCondNum
//               << "     integral: " << integralFrobCondNum << "\n";
    }
}// namespace gismo
