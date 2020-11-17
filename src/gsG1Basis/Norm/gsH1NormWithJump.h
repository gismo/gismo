/** @file gsSeminormH1.h

    @brief Computes the H1 norm with the jump; for approx. g1 Basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/


#pragma once

namespace gismo
{

/** @brief The gsSeminormH1 class provides the functionality
 * to calculate the H1 - seminorm between a field and a function.
 *
 * \ingroup Assembler
*/
template <class T>
class gsH1NormWithJump
{

public:

    gsH1NormWithJump(const gsMultiPatch<> & multiPatch,
                     std::vector<gsMultiBasis<>> & multiBasis,
                     gsSparseMatrix<T> & _field1)
        : m_mp(multiPatch), m_mb(&multiBasis) , m_sparse(&_field1)
    {
        m_value.resize(m_mp.interfaces().size());
        m_value.setZero();
    }


public:

    void compute(gsG1System<T> g1System, bool isogeometric, std::string typeOfnorm)
    {

        for (size_t numInt = 0; numInt < m_mp.interfaces().size(); numInt++ )
        {
            T value = T(0.0);

            gsMatrix<T> quNodes_L; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights
            gsQuadRule<T> QuRule_L; // Reference Quadrature rule

            gsMatrix<T> quNodes_R; // Temp variable for mapped nodes
            gsQuadRule<T> QuRule_R; // Reference Quadrature rule

            // Evaluation flags for the Geometry map
            unsigned evFlags(0);

            boundaryInterface & iFace = m_mp.interfaces()[numInt];

            size_t L = iFace.second().patch;
            size_t R = iFace.first().patch;

            std::vector<size_t> numInt_vector;
            if (typeOfnorm == "edge")
                numInt_vector.push_back(numInt);
            else if (typeOfnorm == "vertex")
            {
                for(size_t numVer=0; numVer < m_mp.vertices().size(); numVer++)
                {
                    std::vector<patchCorner> allcornerLists = m_mp.vertices()[numVer];
                    std::vector<size_t> patchIndex;
                    std::vector<size_t> vertIndex;

                    bool patch_L = false, patch_R = false;
                    for (size_t j = 0; j < allcornerLists.size(); j++)
                    {
                        if (allcornerLists[j].patch == L)
                            patch_L = true;
                        if (allcornerLists[j].patch == R)
                            patch_R = true;
                    }
                    if (patch_L && patch_R)
                        numInt_vector.push_back(numVer);
                }
            }


            // Obtain an integration domain
            gsBasis<T> & dom_L = m_mb->at(isogeometric ? 0 : 1).basis(L);
            gsBasis<T> & dom_R = m_mb->at(isogeometric ? 0 : 1).basis(R);
/*
            if (dom_L.getMinCellLength() < dom_R.getMinCellLength()) // TODO MAYBE ONLY FOR TWO PATCH + FINER GRID LEFT
            {
                gsSparseMatrix<T> sol_sparse(m_sparse->rows(),dom_L.size()*2);
                for (index_t i = 0; i < m_sparse->rows(); i++)
                {
                    gsMatrix<> coefs_temp(numBasisFunctions[R+1] - numBasisFunctions[R],1);
                    for (index_t j = numBasisFunctions[R]; j < numBasisFunctions[R+1]; j++)
                        coefs_temp(j-numBasisFunctions[R],0) = m_sparse->at(i,j);

                    gsGeometry<>::uPtr geo_R = dom_R.makeGeometry(coefs_temp);
                    gsTensorBSpline<2, T> geo_R_sol = dynamic_cast<gsTensorBSpline<2, T> &> (*geo_R);
                    std::vector<int> mul={0, m_mp.patch(R).degree(1) - 1}; // TODO
                    geo_R_sol.uniformRefine(1,mul);

                    for (index_t j = numBasisFunctions[R]; j < numBasisFunctions[R+1]; j++)
                        sol_sparse.at(i,j) = geo_R_sol.coefs()(j - numBasisFunctions[R],0);

                    for (index_t j = numBasisFunctions[L]; j < numBasisFunctions[L+1]; j++)
                        sol_sparse.at(i,j) = m_sparse->at(i,j);
                }

                sol_sparse.prune(1e-14);
                sol_sparse.makeCompressed();

                dom_R = dom_L;
                //m_sparse->swap(sol_sparse);
            }
*/
            boxSide side_L = iFace.second();
            boxSide side_R = iFace.first();

            // Initialize visitor
            initializeb(dom_L, QuRule_L, evFlags, side_L);
            initializeb(dom_R, QuRule_R, evFlags, side_R);

            // Initialize geometry evaluator
            typename gsGeometryEvaluator<T>::uPtr geoEval_L(getEvaluator(evFlags, m_mp.patch(L)));
            // Initialize geometry evaluator
            typename gsGeometryEvaluator<T>::uPtr geoEval_R(getEvaluator(evFlags, m_mp.patch(R)));

            const gsAffineFunction<> ifaceMap(m_mp.getMapForInterface(iFace, dom_R.size() >= dom_L.size() ? 1 : -1));

            typename gsBasis<T>::domainIter domIt = dom_R.size() >= dom_L.size() ? dom_R.makeDomainIterator(side_R) : dom_L.makeDomainIterator(side_L);
            for (; domIt->good(); domIt->next())
            {
                // Map the Quadrature rule to the element
                gsMatrix<T> domItCorner(2,2);
                domItCorner.col(0) = domIt->lowerCorner();
                domItCorner.col(1) = domIt->upperCorner();

                dom_R.size() >= dom_L.size() ? QuRule_R.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes_R, quWeights) :
                    QuRule_L.mapTo(domItCorner.col(0), domItCorner.col(1), quNodes_L, quWeights);

                ifaceMap.eval_into(domItCorner,domItCorner);
                if (domItCorner(1-side_L.direction(),0) > domItCorner(1-side_L.direction(),1) && dom_R.size() >= dom_L.size()) // integral border switched
                {
                    gsMatrix<T> temp_domItCorner = domItCorner;
                    domItCorner.col(1) = temp_domItCorner.col(0);
                    domItCorner.col(0) = temp_domItCorner.col(1);
                }
                else if (domItCorner(1-side_R.direction(),0) > domItCorner(1-side_R.direction(),1) && dom_R.size() < dom_L.size()) // integral border switched
                {
                    gsMatrix<T> temp_domItCorner = domItCorner;
                    domItCorner.col(1) = temp_domItCorner.col(0);
                    domItCorner.col(0) = temp_domItCorner.col(1);
                }

                dom_R.size() >= dom_L.size() ? QuRule_L.mapTo(domItCorner.col(0), domItCorner.col(1), quNodes_L, quWeights):
                    QuRule_R.mapTo(domItCorner.col(0), domItCorner.col(1), quNodes_R, quWeights);

                // Evaluate on quadrature points
                evaluateb(*geoEval_L, quNodes_L, dom_L, *geoEval_R, quNodes_R, dom_R, g1System, numInt_vector, m_sparse, typeOfnorm);

                // Accumulate value from the current element (squared)
                computeb(*geoEval_L, side_L, *geoEval_R, quWeights, value);

            }
            m_value(numInt) = takeRoot(value);
        }

    }


    gsVector<T> value() const { return m_value; }

protected:

    void initializeb(const gsBasis<T> & basis,
                     gsQuadRule<T> & rule,
                     unsigned      & evFlags,
                     boxSide side) // replace with geoEval ?
    {
        // Setup Quadrature
        const unsigned d = basis.dim();
        const int dir = side.direction();
        gsVector<index_t> numQuadNodes( d );
        for (unsigned i = 0; i < d; ++i)
            numQuadNodes[i] = basis.degree(i) + 1;
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE|NEED_MEASURE|NEED_JACOBIAN|NEED_GRAD_TRANSFORM;
    }

    // Evaluate on element.
    void evaluateb(gsGeometryEvaluator<T> & geoEval,
                   gsMatrix<T>            & quNodes,
                   const gsBasis<T> & basis,
                   gsGeometryEvaluator<T> & geoEval_R,
                   gsMatrix<T>            & quNodes_R,
                   const gsBasis<T> & basis_R,
                   gsG1System<T> g1System,
                   std::vector<size_t> numInt,
                   const gsSparseMatrix<T> * sol_sparse,
                   std::string typeOfnorm)
    {

        // Evaluate first function
        gsMatrix<unsigned> actives;
        gsMatrix<T> bGrads;

        basis.active_into(quNodes.col(0), actives);

        // Evaluate basis functions on element
        basis.deriv_into(quNodes, bGrads);

        f1ders.setZero(2, bGrads.cols());
        if (typeOfnorm == "edge")
        {
            for (index_t i = g1System.get_numInterfaceFunctions()[numInt[0]]; i < g1System.get_numInterfaceFunctions()[numInt[0] + 1]; i++)
                for (index_t j = 0; j < actives.rows(); j++)
                    f1ders += sol_sparse->at(i, g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives.at(j))
                        * bGrads.block(2 * j, 0, 2, bGrads.cols());
        }
        else if (typeOfnorm == "vertex")
        {
            for (size_t num = 0; num < numInt.size(); num++)
                for (index_t i = g1System.get_numVertexFunctions()[numInt[num]]; i < g1System.get_numVertexFunctions()[numInt[num] + 1];
                     i++)
                    for (index_t j = 0; j < actives.rows(); j++)
                        f1ders += sol_sparse->at(i, g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives.at(j))
                            * bGrads.block(2 * j, 0, 2, bGrads.cols());
        }
        else if (typeOfnorm == "all")
        {
            for (index_t i = 0; i < sol_sparse->rows()-1; i++)
                for (index_t j = 0; j < actives.rows(); j++)
                    f1ders += sol_sparse->at(i, g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives.at(j))
                        * bGrads.block(2 * j, 0, 2, bGrads.cols());
        }

        geoEval.evaluateAt(quNodes);


        // Evaluate second function
        gsMatrix<unsigned> actives2;
        gsMatrix<T> bGrads2;

        basis_R.active_into(quNodes_R.col(0), actives2);

        // Evaluate basis functions on element
        basis_R.deriv_into(quNodes_R, bGrads2);

        f2ders.setZero(2, bGrads2.cols());
        if (typeOfnorm == "edge")
        {
            for (index_t i = g1System.get_numInterfaceFunctions()[numInt[0]]; i < g1System.get_numInterfaceFunctions()[numInt[0] + 1]; i++)
                for (index_t j = 0; j < actives2.rows(); j++)
                    f2ders += sol_sparse->at(i, g1System.get_numBasisFunctionsInterface()[geoEval_R.id()] + actives2.at(j))
                        * bGrads2.block(2 * j, 0, 2, bGrads2.cols());
        }
        else if (typeOfnorm == "vertex")
        {
            for (size_t num = 0; num < numInt.size(); num++)
                for (index_t i = g1System.get_numVertexFunctions()[numInt[num]]; i < g1System.get_numVertexFunctions()[numInt[num] + 1];
                     i++)
                    for (index_t j = 0; j < actives2.rows(); j++)
                        f2ders += sol_sparse->at(i, g1System.get_numBasisFunctionsInterface()[geoEval_R.id()] + actives2.at(j))
                            * bGrads2.block(2 * j, 0, 2, bGrads2.cols());
        }
        else if (typeOfnorm == "all")
        {
            for (index_t i = 0; i < sol_sparse->rows()-1; i++)
                for (index_t j = 0; j < actives2.rows(); j++)
                    f2ders += sol_sparse->at(i, g1System.get_numBasisFunctionsInterface()[geoEval_R.id()] + actives2.at(j))
                        * bGrads2.block(2 * j, 0, 2, bGrads2.cols());
        }

        geoEval_R.evaluateAt(quNodes_R);
/*
        // alpha
        index_t m_uv = 1;
        gsMatrix<> alpha;
        gsMatrix<> uv, ev;
        // alpha^S
        if (m_uv==1)
        {
            uv.setZero(2,quNodes.cols());
            uv.bottomRows(1) = quNodes.row(m_uv); // v
        }
        else if (m_uv==0)
        {
            uv.setZero(2,quNodes.cols());
            uv.topRows(1) = quNodes.row(m_uv); // u
        }

        // ======== Determine bar{alpha^(L)} == Patch 0 ========
        const gsGeometry<> & P0 = m_mp.patch(0); // iFace.second().patch = 0

        for (index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i), ev);
            uv(0, i) = 1 * ev.determinant();
            uv(1,i) = (ev.col(m_uv)).norm();
        }
        gsMatrix<> ones;
        ones.setOnes(1,quNodes.cols());
        alpha = ones.cwiseQuotient(uv.row(0).cwiseProduct(uv.row(1)));


        // beta
        // beta^S
        gsMatrix<> beta;
        if (m_uv==1)
        {
            uv.setZero(2,quNodes.cols());
            uv.bottomRows(1) = quNodes.row(m_uv); // v
        }
        else if (m_uv==0)
        {
            uv.setZero(2,quNodes.cols());
            uv.topRows(1) = quNodes.row(m_uv); // u
        }

        const index_t d = m_mp.parDim();
        gsVector<> D0(d);

        // ======== Determine bar{beta}^L ========
        for(index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i),ev);
            D0 = ev.col(m_uv);
            real_t D1 = 1/ D0.norm();
            uv(0,i) = - 1 * D1 * D1 * ev.col(1).transpose() * ev.col(0);

        }

        beta = uv.row(0);

        normal.setOnes(2, quNodes.cols());
        //normal.row(0) = alpha;
        normal.row(1) = - beta;//.cwiseProduct(alpha);
*/
    }

    // assemble on element
    inline T computeb(gsGeometryEvaluator<T> & geoEval,
                      boxSide side,
                      gsGeometryEvaluator<T> & geoEval_R,
                      gsVector<T> const      & quWeights,
                      T & accumulated)
    {

        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            //const T d = element.dim();
            // Transform the gradients
            geoEval.transformGradients(k, f1ders, f1pders);
            geoEval_R.transformGradients(k, f2ders, f2pders);

            // Compute the unit normal
            gsVector<T> unormal;

            geoEval.outerNormal(k, side, unormal);

            const T weight = quWeights[k] * unormal.norm() ;

            // f2ders : N X 1
            sum += weight * ( (f1pders - f2pders).transpose() * unormal ).squaredNorm() ;
            //sum += weight * ( (f1ders - f2ders).transpose() * normal ).squaredNorm() ;
        }
        accumulated += sum;

        return sum;
    }

    inline T takeRoot(const T v) { return math::sqrt(v);}



private:
    gsMatrix<T> f1ders, f2ders;
    gsMatrix<T> f1pders, f2pders; // f2pders only needed if f2param = true

    gsMatrix<T> normal;

    gsMultiPatch<T> m_mp;
    std::vector<gsMultiBasis<>> * m_mb;
    gsSparseMatrix<T> * m_sparse;

protected:
    gsVector<T> m_value;     // the total value of the norm


};


} // namespace gismo