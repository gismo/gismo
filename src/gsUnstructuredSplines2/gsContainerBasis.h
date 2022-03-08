/** @file gsContainerBasis.h

    @brief Provides declaration of Basis abstract interface. Similar to gsMultiBasis, but
    without topology.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/


#pragma once

#include<gsCore/gsBasis.h>
#include<gsSolver/gsBlockOp.h>


namespace gismo
{

    template<short_t d,class T>
    class gsContainerBasis  : public gsBasis<T>
    {

    public:
        /// Shared pointer for gsContainerBasis
        typedef memory::shared_ptr< gsContainerBasis > Ptr;

        /// Unique pointer for gsContainerBasis
        typedef memory::unique_ptr< gsContainerBasis > uPtr;

        gsContainerBasis(index_t numSubspaces)
        {
            basisContainer.clear();
            for (index_t i = 0; i != numSubspaces; i++)
                basisContainer.push_back(gsTensorBSplineBasis<d,T>());

        }

        ~gsContainerBasis() {};

        // boundary(boxSide const & side)

// implementations of gsBasis
    public:

        static uPtr make(   const gsContainerBasis& other)
        { return uPtr( new gsContainerBasis( other ) ); }


    public:

    GISMO_CLONE_FUNCTION(gsContainerBasis)

        short_t domainDim() const
        {
            return d;
        }

        void connectivity(const gsMatrix<T> & nodes,
                          gsMesh<T>   & mesh) const
        {
            GISMO_UNUSED(nodes); GISMO_UNUSED(mesh);
            GISMO_NO_IMPLEMENTATION;
        }

        memory::unique_ptr<gsGeometry<T> > makeGeometry( gsMatrix<T> coefs ) const
        {
            GISMO_UNUSED(coefs);
            GISMO_NO_IMPLEMENTATION;
        }

        std::ostream &print(std::ostream &os) const
        {
            GISMO_UNUSED(os);
            GISMO_NO_IMPLEMENTATION;
        }

        void uniformRefine()
        {
            for (size_t i=0; i< basisContainer.size(); ++i)
                basisContainer[i].uniformRefine();
        }

        // Returm max degree of all the spaces, otherwise i =
        short_t degree(short_t dir) const
        {
            short_t deg = 0;
            for (size_t i=0; i< basisContainer.size(); ++i)
                if (basisContainer[i].degree(dir) > deg)
                    deg = basisContainer[i].degree(dir);

            return deg;
        }

        index_t size() const {
            index_t sz = 0;
            for (size_t i=0; i< basisContainer.size(); ++i)
                sz += basisContainer[i].size();
            return sz;
        }

        void reverse()
        {
            for (size_t i=0; i< basisContainer.size(); ++i)
            {
                if (gsTensorBSplineBasis<d, T> * mb =
                        dynamic_cast<gsTensorBSplineBasis<d, T> *>(&basisContainer[i]) )
                {
                    gsTensorBSplineBasis<d, T> basis_temp = dynamic_cast<gsTensorBSplineBasis<d,T>&>(basisContainer[i]);
                    gsTensorBSplineBasis<d, T> newTensorBasis(basis_temp.knots(1),basis_temp.knots(0));
                    basisContainer[i] = newTensorBasis;
                }
                /*
                else if (gsTensorNurbsBasis<d, T> * mb =
                        dynamic_cast<gsTensorNurbsBasis<d, T> *>(&basisContainer[i]) )
                {
                    gsTensorNurbsBasis<d, T> basis_temp = dynamic_cast<gsTensorNurbsBasis<d,T>&>(basisContainer[i]);
                    basis_temp.swapDirections(1,0);
                    basisContainer[i] = basis_temp;
                }
                */
                else
                    gsInfo << "Works for now just with gsTensorBSplineBasis<d, T> \n";

            }
        }

        index_t nPieces() const {return basisContainer.size();}

        typename gsBasis<T>::domainIter makeDomainIterator(const boxSide & side) const
        {
            // Using the inner basis for iterating
            return basisContainer[0].makeDomainIterator(side);
        }

        typename gsBasis<T>::domainIter makeDomainIterator() const
        {
            // Using the inner basis for iterating
            return basisContainer[0].makeDomainIterator();
        }

        const gsBasis<T> & component(short_t i) const
        {
            return basisContainer[0].component(i);
        }

/*        void matchWith(const boundaryInterface & bi, const gsBasis<T> & other,
                       gsMatrix<index_t> & bndThis, gsMatrix<index_t> & bndOther) const;
*/


        gsMatrix<T> support() const
        {
            return basisContainer[0].support();
        }

        void uniformCoarsen_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots)
        {
//            gsMatrix<T> matrix(0,0);
//            gsDebugVar(basisContainer.size());
//            for (size_t i=0; i< basisContainer.size(); ++i)
//            {
//                gsDebugVar(i);
//                gsSparseMatrix<T,RowMajor> mat_temp;
//                basisContainer[i].uniformCoarsen_withTransfer(mat_temp, numKnots);
//                matrix.conservativeResize(matrix.rows()+mat_temp.rows(), matrix.cols() + mat_temp.cols());
//                matrix.bottomRightCorner(mat_temp.rows(), mat_temp.cols()) = mat_temp;
//                gsDebugVar(matrix.size());
//            }
//            transfer = matrix.sparseView();


//            index_t sz = basisContainer.size();
//            typename gsBlockOp<T>::Ptr block = gsBlockOp<T>::make(sz,sz);
//
//            for (size_t i=0; i< basisContainer.size(); ++i)
//            {
//                gsDebugVar(i);
//                gsSparseMatrix<T,RowMajor> mat_temp;
//                basisContainer[i].uniformCoarsen_withTransfer(mat_temp, numKnots);
//                block->addOperator(i,i,makeMatrixOp(mat_temp.moveToPtr()));
//            }
//            gsDebugVar(block->cols());
//            gsDebugVar(block->rows());
//            gsMatrix<> matrix;
//            block->toMatrix(matrix);
//            gsDebugVar(matrix.rows());
//            gsDebugVar(matrix.cols());
//            transfer = matrix.sparseView();

            basisContainer[0].uniformCoarsen_withTransfer(transfer, numKnots);
            gsDebugVar(transfer.rows());
            gsDebugVar(transfer.cols());
        }

        gsBasis<T>* boundaryBasis_impl(boxSide const & s) const
        {
            /*
            if (basisContainer.size() != 1)
            {
                typename gsBSplineTraits<d-1, T>::Basis* bBasis = new typename gsBSplineTraits<d-1, T>::Basis(*basisContainer[s.index()].boundaryBasis(s));
                return bBasis;
            }
            else {
                typename gsBSplineTraits<d-1, T>::Basis* bBasis = new typename gsBSplineTraits<d-1, T>::Basis(*basisContainer[0].boundaryBasis(s));
                return bBasis;
            }
             */

            gsDebugVar(s.index());
            // Maybe not working for approx C1 Basis functions
            typename gsBSplineTraits<d-1, T>::Basis* bBasis = new typename gsBSplineTraits<d-1, T>::Basis(*basisContainer[0].boundaryBasis(s));
            gsDebugVar(*bBasis);
            return bBasis;
        }

        void active_into(const gsMatrix<T> & u, gsMatrix<index_t> & result) const
        {
            GISMO_ASSERT(u.rows() == d, "Dimension of the points in active_into is wrong");
            //GISMO_ASSERT(u.cols() == 1, "Active_into is wrong");

            index_t nr = 0;
            std::vector<gsMatrix<index_t>> result_temp;
            result_temp.resize(basisContainer.size());
            for (size_t i=0; i< basisContainer.size(); ++i)
            {
                basisContainer[i].active_into(u, result_temp[i]);
                nr += result_temp[i].rows();
            }

            result.resize(nr,u.cols());

            index_t shift = 0;
            index_t shift_rows = 0;
            for (size_t i=0; i< basisContainer.size(); ++i)
            {
                result_temp[i].array() += shift;
                result.block(shift_rows, 0, result_temp[i].rows(), u.cols()) = result_temp[i];
                shift += basisContainer[i].size();
                shift_rows += result_temp[i].rows();

            }
        }

        void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
        {
            result.resize(0, u.cols());
            for (size_t i=0; i< basisContainer.size(); ++i)
            {
                gsMatrix<T> result_temp;
                basisContainer[i].eval_into(u, result_temp);
                result.conservativeResize(result.rows()+result_temp.rows(), result.cols());
                result.bottomRows(result_temp.rows()) = result_temp;
            }
        }

        void evalSingle_into(index_t k, const gsMatrix<T> & u, gsMatrix<T>& result) const
        {
            result.resize(0, u.cols());
            for (size_t i=0; i< basisContainer.size(); ++i)
            {
                gsMatrix<T> result_temp;
                basisContainer[i].evalSingle_into(k, u, result_temp);
                result.conservativeResize(result.rows()+result_temp.rows(), result.cols());
                result.bottomRows(result_temp.rows()) = result_temp;
            }
        }

        void deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
        {
            result.resize(0, u.cols());
            for (size_t i=0; i< basisContainer.size(); ++i)
            {
                gsMatrix<T> result_temp;
                basisContainer[i].deriv_into(u, result_temp);
                result.conservativeResize(result.rows()+result_temp.rows(), result.cols());
                result.bottomRows(result_temp.rows()) = result_temp;

            }
        }

        void deriv2_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
        {
            result.resize(0, u.cols());
            for (size_t i=0; i< basisContainer.size(); ++i)
            {
                gsMatrix<T> result_temp;
                basisContainer[i].deriv2_into(u, result_temp);
                result.conservativeResize(result.rows()+result_temp.rows(), result.cols());
                result.bottomRows(result_temp.rows()) = result_temp;

            }
        }

        void matchWith(const boundaryInterface &bi, const gsBasis<T> &other, gsMatrix<int> &bndThis,
                       gsMatrix<int> &bndOther) const {
            // First side
            // edge
            gsDebug << "matchWith() function is not implemented yet \n";

         /*
            gsMatrix<index_t> indizes1, indizes2;
            short_t side_id = bi.first().side().index();
            indizes1 = basisG1Container[side_id].boundaryOffset(boxSide(bi.first().side()), 0);
            indizes2 = basisG1Container[side_id].boundaryOffset(boxSide(bi.first().side()), 1);

            bndThis = indizes1;
            bndThis.resize(indizes1.rows() + indizes2.rows(), indizes1.cols());
            bndThis.block(0, 0, indizes1.rows(), indizes1.cols()) = indizes1;
            bndThis.block(indizes1.rows(), 0, indizes2.rows(), indizes1.cols()) = indizes2;

            index_t shift = 0;
            for (index_t i = 0; i < side_id; ++i)
                shift += basisG1Container[side_id].size();

            bndOther.array() += shift;

            bndOther = bndThis;
        */
        }


        gsMatrix<int> boundaryOffset(const boxSide & bside, int offset) const
        {
            gsMatrix<int> result;

            // TODO TEST WITH:
            /*
            index_t shift = 0;
            result.resize(0, 1);
            for (size_t i=0; i< basisContainer.size(); ++i)
            {
                gsMatrix<int> result_temp;
                result_temp = basisContainer[i].boundaryOffset(bside, offset);
                result_temp.array() += shift;
                result.conservativeResize(result.rows()+result_temp.rows(), 1);
                result.bottomRows(result_temp.rows()) = result_temp;
                shift += basisContainer[i].size();
            }
            */
            // Edges
            short_t side_id = bside.index();
            if (basisContainer.size() != 1)
            {
                index_t shift = 0;
                for (index_t i=0; i< side_id; ++i)
                    shift += basisContainer[i].size();
                result = basisContainer[side_id].boundaryOffset(bside, offset);
                result.array() += shift;
            }
            else if (basisContainer.size() == 1)
                result = basisContainer[0].boundaryOffset(bside, offset);


            // Vertices:
            if (basisContainer.size() != 1)
            {
                std::vector<boxCorner> containedCorners;
                bside.getContainedCorners(d, containedCorners);

                GISMO_ASSERT(containedCorners.size() != 0, "No contained corner");


                for (size_t nc = 0; nc < containedCorners.size(); nc++) {
                    index_t corner_id = containedCorners[nc].m_index + 4; // + 4 included bcs of 4 sides!

                    index_t shift = 0;
                    for (index_t i = 0; i < corner_id; ++i)
                        shift += basisContainer[i].size();

                    gsMatrix<int> result_temp;
                    result_temp = basisContainer[corner_id].boundaryOffset(bside, offset);
                    result_temp.array() += shift;

                    index_t r_rows = result.rows() + result_temp.rows();
                    result.conservativeResize(r_rows, 1);
                    result.bottomRows(result_temp.rows()) = result_temp;
                }
            }
            return result;
        }


// implementations of gsContainerBasis
    public:

        // basisContainer:
        void setBasis(index_t row, gsTensorBSplineBasis<d,T> basis) { basisContainer[row] = basis; }

        const gsBasis<T> & piece(const index_t k) const { return basisContainer[k]; }
        // basisContainer END

        //std::vector<gsTensorBSplineBasis<d,T>> & getBasisContainer() { return basisContainer; }

        // Data members
    protected:

        // Collection of the subspaces
        std::vector<gsTensorBSplineBasis<d, T>> basisContainer;
    };

}
