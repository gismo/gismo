/** @file gsNormL2.h

    @brief Computes the L2 norm.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat
*/


#pragma once

#include <gsG1Basis/Norm/gsNorm.h>
#include <gsG1Basis/Norm/gsG1ASVisitorResidualNormL2.h>


namespace gismo
{


/** @brief The gsNormL2 class provides the functionality
 * to calculate the L2 - norm between a field and a function.
 *
 * \ingroup Assembler
*/
template <class T = real_t, class Visitor = gsG1ASVisitorResidualNormL2<T> >
class gsG1ASResidualNormL2
{

public:

    gsG1ASResidualNormL2(const gsMultiPatch<T> & multiPatch,
                         std::vector<gsSparseMatrix<>> & _field1,
                         std::vector<gsMultiBasis<>> & _func2,
                         bool _f2param = false)
                        : patchesPtr( &multiPatch ),
                          sparseMatrix(&_field1), basisVec(&_func2), f2param(_f2param)
    {

    }


public:

    void plotElWiseError(std::vector<real_t> &elWiseError)
    {

        std::string fn = "ElWiseError";
        index_t npts = elWiseError.size() / 2;
        gsParaviewCollection collection2(fn);
        std::string fileName2;

        for ( size_t pp =0; pp < patchesPtr->nPatches(); ++pp ) // Patches
        {
            fileName2 = fn + util::to_string(pp);
            //writeSinglePatchField( field, i, fileName, npts );

            const gsFunction<T> & geometry = patchesPtr->patch(pp);

            const int n = geometry.targetDim();
            const int d = geometry.domainDim();

            gsMatrix<T> ab = geometry.support();
            gsVector<T> a = ab.col(0);
            gsVector<T> b = ab.col(1);

            gsVector<unsigned> np = uniformSampleCount(a, b, npts);

            gsMatrix<T> pts = gsPointGrid(a, b, np);
            gsMatrix<T> eval_geo = geometry.eval(pts);//pts

            gsMatrix<T> eval_field(1,npts);
            for (index_t i = 0; i < npts; i++)
            {
                eval_field(0,i) = elWiseError[pp*npts + i];
            }

            // Here add g1 basis
            //eval_field.setZero();

            if ( 3 - d > 0 )
            {
                np.conservativeResize(3);
                np.bottomRows(3-d).setOnes();
            }
            else if (d > 3)
            {
                gsWarn<< "Cannot plot 4D data.\n";
                return;
            }

            if ( 3 - n > 0 )
            {
                eval_geo.conservativeResize(3,eval_geo.cols() );
                eval_geo.bottomRows(3-n).setZero();
            }
            else if (n > 3)
            {
                gsWarn<< "Data is more than 3 dimensions.\n";
            }

            if ( eval_field.rows() == 2)
            {
                eval_field.conservativeResize(3,eval_geo.cols() );
                eval_field.bottomRows(1).setZero(); // 3-field.dim()
            }
            if ( eval_field.rows() == 1)
            {
                eval_field.conservativeResize(3,eval_geo.cols() );
                eval_field.bottomRows(2).setZero(); // 3-field.dim()
            }

            gsWriteParaviewTPgrid(eval_geo, eval_field, np.template cast<index_t>(), fileName2);


            collection2.addPart(fileName2, ".vts");
        }
        collection2.save();
    }

    /// @brief Returns the computed norm value
    T value() const { return m_value; }

    std::vector<T> elWise_value() const { return m_elWise; }

    void compute(std::vector<gsG1System<real_t>> & g1SysVec, bool storeElWise = false)
    {
        boxSide side = boundary::none;

        if ( storeElWise )
            m_elWise.clear();

        m_value = T(0.0);

#pragma omp parallel
        {
#ifdef _OPENMP
            // Create thread-private visitor
        Visitor visitor;
        const int tid = omp_get_thread_num();
        const int nt  = omp_get_num_threads();
#else
            Visitor visitor;
#endif
            gsMatrix<T> quNodes; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights
            gsQuadRule<T> QuRule; // Reference Quadrature rule

            // Evaluation flags for the Geometry map
            unsigned evFlags(0);

            for (size_t pn = 0; pn < patchesPtr->nPatches(); ++pn)// for all patches
            {
                // Obtain an integration domain
                visitor.initialize(basisVec, QuRule, evFlags);

                // Initialize geometry evaluator
                typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, patchesPtr->patch(pn)));

                const gsBasis<T> & dom = basisVec->at(1).basis(pn);

                typename gsBasis<T>::domainIter domIt = dom.makeDomainIterator(side);

                // TODO: optimization of the assembling routine, it's too slow for now
                // Start iteration over elements
#ifdef _OPENMP
                for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
                for (; domIt->good(); domIt->next() )
#endif
                {

                    // Map the Quadrature rule to the element
                    QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

                    // Evaluate on quadrature points
                    visitor.evaluate(*geoEval, sparseMatrix, basisVec, g1SysVec, quNodes);


                    //visitor.compute(*domIt, *geoEval, quWeights, m_value);
                    T temp = 0.0;
                    const T result = visitor.compute(*domIt, *geoEval, quWeights, temp);
#pragma omp critical
                    {
                        m_value += result;
                        if (storeElWise)
                            m_elWise.push_back(takeRoot(result));
                    };
                }
            }

        }//omp parallel
        m_value = takeRoot(m_value);

    }

protected:

    inline T takeRoot(const T v)
    {

                return math::pow(v, static_cast<T>(1)/2 );
        
    }



protected:

    const gsMultiPatch<T> * patchesPtr;

    const std::vector<gsSparseMatrix<>> * sparseMatrix;

    const std::vector<gsMultiBasis<>> * basisVec;

private:

    bool f2param;

protected:
    std::vector< gsMultiPatch<>> m_G1Basis;


protected:

    std::vector<T> m_elWise;    // vector of the element-wise values of the norm
    T              m_value;     // the total value of the norm




};



};// namespace gismo





