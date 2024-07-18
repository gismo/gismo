/** @file gsLofting.hpp

    @brief Implementation of gsLofting class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Imperatore
*/

#pragma once

# include <gsNurbs/gsKnotVector.h>
# include <gsNurbs/gsBSpline.h>
# include <gsNurbs/gsTensorBSpline.h>
# include <vector>


namespace gismo
{
    //destructor
    template<class T>
    gsLofting<T>::~gsLofting()
    {
        if ( m_result )
            delete m_result;
    }

    // constructors
    template<class T>
    gsLofting<T>::gsLofting(gsMultiPatch<T> container)
    {
        GISMO_ASSERT(container.nPatches() > 2, "Not enough lofting curves: #curves  = " << container.nPatches());
        m_container = container;
        m_result = nullptr;

    }

    template<class T>
    gsLofting<T>::gsLofting(gsMultiPatch<T> container, index_t deg_v)
    {
        GISMO_ASSERT(container.nPatches() > 2, "Not enough lofting curves: #curves  = " << container.nPatches());
        GISMO_ASSERT(container.nPatches() > deg_v, "Degree in v-direction must be lower, i.e., deg_v < " << container.nPatches());
        m_container = container;
        m_vdeg = deg_v;

        // make compatible curves
        m_result = nullptr;
    }

    template<class T>
    gsLofting<T>::gsLofting(gsMultiPatch<T> container, gsKnotVector<T> kv)
    {
        GISMO_ASSERT(container.nPatches() > 2, "Not enough lofting curves: #curves  = " << container.nPatches());
        GISMO_ASSERT(container.nPatches() == kv.size() - (kv.degree()+1), "Wrong dimension of univariate space v-direction: " << kv.size() - (kv.degree()+1) << "!=" << container.nPatches());
        m_container = container;
        m_N = container.nPatches(); // N number of section curves
        m_kv = kv;

        // make compatible curves
        const gsBSpline<T> & crv = static_cast<const gsBSpline<> &>(container.patch(0));
        m_ku = crv.basis().knots();

        gsTensorBSplineBasis<2, T> basis(m_ku, m_kv);
        m_basis = &basis;

        m_m = container.patch(0).coefs().rows(); //number of control points for each curve
        m_DIM = container.patch(0).coefs().cols(); //number of control points for each curve

        m_result = nullptr;
    }

    // functions

    // compute the iso-parametric values 
    // alpha : 1-D parameterization method = (1) chord-length; (0.5) centripetal; (0) uniform.
    template<class T>
    void gsLofting<T>::isovparameters(T alpha, gsMatrix<T> & v_parameters)
    {   
        
        //index_t N = m_container.nPatches(); // number of section curves
        //index_t m_m = m_container.patch(0).coefs().rows(); // number of control points in each section curve

        v_parameters.resize(1, m_N);
        v_parameters(0,0)   = 0;
        v_parameters(0,m_N-1) = 1;


        std::vector<T> distances(m_N-1);
        std::vector<T> L(m_m);
        
        for(index_t i = 0; i < m_m; i++) // number of control points in each section curve
        {
         for(index_t k = 1; k < m_N; k++) // number of section curves
         {
            distances.at(k-1) = math::pow((m_container.patch(k).coefs().row(i)-m_container.patch(k-1).coefs().row(i)).norm(),alpha);
         }
            L.at(i) = std::accumulate(distances.begin(),distances.end(), T(0));
        }

        
        for(index_t k = 1; k < m_N-1; k++)//k=2:N-1
        {
         T num = 0;

         for(index_t i = 0; i < m_m; i++) // number of control points in each section curve i=1:m
         {
            num += math::pow((m_container.patch(k).coefs().row(i)-m_container.patch(k-1).coefs().row(i)).norm(),alpha)/L[i];
         }  
            v_parameters(0,k) = v_parameters(0, k-1) + (1./m_m)*num;
        }
    }

    template<class T>
    void gsLofting<T>::compute()
    {
        // Wipe out previous result
        if ( m_result )
            delete m_result;


        gsBSplineBasis<T> bv(m_kv);
        gsMatrix<T> params;
        isovparameters(1., params);

        gsSparseMatrix<T> N_mat = bv.collocationMatrix(params);
        gsDebugVar(N_mat.cols());
        gsDebugVar(N_mat.rows());

        N_mat.makeCompressed();

        typename gsSparseSolver<T>::BiCGSTABILUT solver( N_mat );

        if ( solver.preconditioner().info() != gsEigen::Success )
        {
            gsWarn<<  "The preconditioner failed. Aborting.\n";
            
            return;
        }

        //gsMatrix<T> PX(m_m,m_N);, PY(m_m,m_N), PZ(m_m,m_N);

        gsMatrix<T> coefs(m_m * m_N, m_DIM);
        gsMatrix<T> c;
        coefs.setZero();

        std::vector<gsMatrix<T>> store_c;
        for(index_t dd = 0; dd < m_DIM; dd++)
        {
            gsMatrix<T> P(m_m,m_N), cc(m_N, m_m);
            for(index_t col = 0; col < m_N; col ++)
            {
                P.col(col) = m_container.patch(col).coefs().col(dd);
            }

            P.transposeInPlace();
            c = solver.solve(P); //toDense() 
            // store_c.push_back(c.transposeInPlace());
            c.transposeInPlace();
            store_c.push_back(c);
        }

        for (index_t dd = 0; dd < m_DIM; dd++)
        {
            for(index_t n = 0; n < m_N; n++)
            {
                for(index_t m = 0; m < m_m; m++)
                    coefs(m + n * m_m, dd) = store_c[dd](m,n);
            }
        }

        gsInfo << "Until here?\n";
        gsInfo << "Coefs:\n" << coefs << "\n";
        //const gsBasis<T> * bb = dynamic_cast<const gsBasis<T> *>(m_basis);
        // gsTensorBSplineBasis<2,T> * bb = dynamic_cast< gsTensorBSplineBasis<2,T> *>(m_basis);
        
        gsTensorBSplineBasis<2, T> bb(m_ku, m_kv);
        gsDebugVar(bb);
        m_result = bb.makeGeometry( give(coefs) ).release();

    }

    /// Make compatible curves to be lofted.
    /// compatible curves: same degree and knot vectors.
    template<class T>
    void gsLofting<T>::make_compatible()
    {
        index_t maxD = 0;
        // for(index_t dd = 0; dd < m_container.size(); dd++)
        // {
        //     if (maxD < m_container[dd].degree())
        //         maxD < m_container[dd].degree();
        // }

        // for(index_t dd = 0; dd < m_container.size(); dd++)
        // {
        //     index_t currD = m_container[dd].degree();
        //     m_container[dd].degreeElevate(maxD-currD);
        // }
    }




} // namespace gismo