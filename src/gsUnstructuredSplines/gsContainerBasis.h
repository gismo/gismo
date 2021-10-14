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

        gsContainerBasis(index_t numberSubspaces)
        {
            basisContainer.clear();
            for (index_t i = 0; i != numberSubspaces; i++)
                basisContainer.addBasis(gsBasis<T>*() );
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

        void matchWith(const boundaryInterface & bi, const gsBasis<T> & other,
                       gsMatrix<index_t> & bndThis, gsMatrix<index_t> & bndOther) const;

        void active_into(const gsMatrix<T> & u, gsMatrix<index_t> & result) const
        {
            GISMO_ASSERT(u.rows() == d, "Dimension of the points in active_into is wrong");
            //GISMO_ASSERT(u.cols() == 1, "Active_into is wrong");

            index_t nr = 0;
            std::vector<gsMatrix<index_t>> result_temp;
            result_temp.resize(basisContainer.size());
            for (size_t i=0; i< basisContainer.size(); ++i)
                if (rowContainer[i] != 0)
                {
                    basisContainer[i].active_into(u, result_temp[i]);
                    nr += result_temp[i].rows();
                }

            result.resize(nr,u.cols());

            index_t shift = 0;
            index_t shift_rows = 0;
            for (size_t i=0; i< basisContainer.size(); ++i)
            {
                if (rowContainer[i] != 0)
                {
                    result_temp[i].array() += shift;
                    result.block(shift_rows, 0, result_temp[i].rows(), u.cols()) = result_temp[i];
                    shift += basisContainer[i].size();
                    shift_rows += result_temp[i].rows();
                }
            }
        }

        void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
        {
            result.resize(0, u.cols());
            for (size_t i=0; i< basisContainer.size(); ++i)
            {
                if (rowContainer[i] != 0)
                {
                    gsMatrix<T> result_temp;
                    basisContainer[i].eval_into(u, result_temp);
                    result.conservativeResize(result.rows()+result_temp.rows(), result.cols());
                    result.bottomRows(result_temp.rows()) = result_temp;
                }
            }
        }

        void deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
        {
            result.resize(0, u.cols());
            for (size_t i=0; i< basisContainer.size(); ++i)
            {
                if (rowContainer[i] != 0)
                {
                    gsMatrix<T> result_temp;
                    basisContainer[i].deriv_into(u, result_temp);
                    result.conservativeResize(result.rows()+result_temp.rows(), result.cols());
                    result.bottomRows(result_temp.rows()) = result_temp;
                }
            }
        }

        void deriv2_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
        {
            result.resize(0, u.cols());
            for (size_t i=0; i< basisContainer.size(); ++i)
            {
                if (rowContainer[i] != 0)
                {
                    gsMatrix<T> result_temp;
                    basisContainer[i].deriv2_into(u, result_temp);
                    result.conservativeResize(result.rows()+result_temp.rows(), result.cols());
                    result.bottomRows(result_temp.rows()) = result_temp;
                }
            }
        }


        gsMatrix<int> boundaryOffset(const boxSide &side, int offset) const
        {
            gsMatrix<int> result;

            // Edges
            short_t side_id = side.index();
            index_t shift = 0;
            for (index_t i=0; i< side_id; ++i)
                shift += basisContainer[i].size();

            result = basisContainer[side_id].boundaryOffset(boxSide(side_id), offset);
            result.array() += shift;

            // Vertices:
            std::vector<boxCorner> containedCorners;
            side.getContainedCorners(d, containedCorners);

            GISMO_ASSERT(containedCorners.size() != 0, "No contained corner");

            for (size_t nc = 0; nc < containedCorners.size(); nc++) {
                index_t corner_id = containedCorners[nc].m_index + 4; // + 4 included bcs of 4 sides!

                index_t shift = 0;
                for (index_t i=0; i< corner_id; ++i)
                    shift += basisContainer[i].size();

                gsMatrix<int> result_temp;
                result_temp = basisContainer[corner_id].boundaryOffset(boxSide(side_id), offset);
                result_temp.array() += shift;

                result.conservativeResize(result.rows()+result_temp.rows(), 1 );
                result.bottomRows(result_temp.rows()) = result_temp;

                if (offset == 1) // DIRTY AND QUICK
                {
                    result_temp = basisContainer[corner_id].boundaryOffset(boxSide(side_id), offset+1); // TODO
                    result_temp.array() += shift;

                    result.conservativeResize(result.rows()+result_temp.rows(), 1 );
                    result.bottomRows(result_temp.rows()) = result_temp;
                }

            }

            return result;
        }


// implementations of gsContainerBasis
    public:

        // basisContainer:
        void setBasis(index_t row, gsBasis<T> & innerBasis) { basisContainer[row] = innerBasis; }
        gsBasis<T> & getBasis(index_t row) { return basisContainer[row]; }
        // basisContainer END

        gsMultiBasis<T> & getBasisContainer() { return basisContainer; }

        // Data members
    protected:

        // Collection of the subspaces
        gsMultiBasis<T> basisContainer;
    };

}
