/** @file gsExprAssemblerSink_test.cpp

    @brief Checks that assembly into an external sink (assemble_into and
    friends) gives the same system as the internal fiber matrix, and that
    the patterns passed by computePattern*_into contain every assembled
    entry, including the couplings across interfaces.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "gismo_unittest.h"
#include <set>

namespace
{
using namespace gismo;

// Records everything a sink receives
struct RecordingSink
{
    std::vector<gsEigen::Triplet<real_t,index_t> > mat;
    std::vector<std::pair<index_t,real_t> > rhs;
    std::set<std::pair<index_t,index_t> > pattern;

    void addMatrix(const gsVector<index_t> & rows, const gsVector<index_t> & cols,
                   const gsMatrix<real_t> & block)
    {
        for (index_t j = 0; j != cols.size(); ++j)
            for (index_t i = 0; i != rows.size(); ++i)
                if (rows[i] >= 0 && cols[j] >= 0)
                    mat.push_back(gsEigen::Triplet<real_t,index_t>(rows[i], cols[j], block(i,j)));
    }

    void addRhs(const gsVector<index_t> & rows, const gsMatrix<real_t> & block)
    {
        for (index_t i = 0; i != rows.size(); ++i)
            if (rows[i] >= 0)
                rhs.push_back(std::make_pair(rows[i], block(i,0)));
    }

    void addPattern(const gsVector<index_t> & rows, const gsVector<index_t> & cols)
    {
        for (index_t j = 0; j != cols.size(); ++j)
            for (index_t i = 0; i != rows.size(); ++i)
                if (rows[i] >= 0 && cols[j] >= 0)
                    pattern.insert(std::make_pair(rows[i], cols[j]));
    }

    gsSparseMatrix<real_t> matrix(index_t n) const
    {
        gsSparseMatrix<real_t> K(n, n);
        K.setFromTriplets(mat.begin(), mat.end()); // sums duplicates
        return K;
    }

    gsMatrix<real_t> vector(index_t n) const
    {
        gsMatrix<real_t> b;
        b.setZero(n, 1);
        for (const auto & e : rhs) b(e.first, 0) += e.second;
        return b;
    }
};

}

SUITE(gsExprAssemblerSink_test)
{
    TEST(SinkMatchesFiberMatrix)
    {
        // two patches, side by side
        gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(2, 1, 1.0);
        gsMultiBasis<> mb(mp, true);
        mb.setDegree(2);
        mb.uniformRefine();
        mb.uniformRefine();

        gsFunctionExpr<> f("x*y+1", 2), g("sin(x)+y", 2);
        gsBoundaryConditions<> bc;
        bc.addCondition(0, boundary::west , condition_type::dirichlet, &g);
        bc.addCondition(0, boundary::south, condition_type::dirichlet, &g);
        bc.addCondition(1, boundary::east , condition_type::neumann  , &g);
        bc.addCondition(1, boundary::north, condition_type::neumann  , &g);
        bc.setGeoMap(mp);

        // glued (0) and discontinuous (-1) interfaces
        for (index_t icont = 0; icont >= -1; --icont)
        {
            gsExprAssembler<> A(1,1);
            A.setIntegrationElements(mb);
            gsExprAssembler<>::geometryMap G = A.getMap(mp);
            gsExprAssembler<>::space u = A.getSpace(mb);
            auto ff = A.getCoeff(f, G);
            auto gN = A.getBdrFunction(G);
            u.setup(bc, dirichlet::interpolation, icont);
            A.initSystem();
            const index_t n = A.numDofs();
            const real_t alpha = 10;

            // the same terms through the internal fiber matrix ...
            A.assemble(igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G));
            A.assembleBdr(bc.get("Neumann"), u * u.tr() * nv(G).norm(), u * gN.tr() * nv(G).norm());
            A.assembleIfc(mp.interfaces(),
                           alpha * u.left() * u.left().tr()  * nv(G).norm(),
                          -alpha * u.right()* u.left().tr()  * nv(G).norm(),
                          -alpha * u.left() * u.right().tr() * nv(G).norm(),
                           alpha * u.right()* u.right().tr() * nv(G).norm());
            const gsSparseMatrix<real_t> K = A.matrix();
            const gsMatrix<real_t> b = A.rhs();

            // ... and through a sink
            RecordingSink S;
            A.assemble_into(S, igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G));
            A.assembleBdr_into(S, bc.get("Neumann"), u * u.tr() * nv(G).norm(), u * gN.tr() * nv(G).norm());
            A.assembleIfc_into(S, mp.interfaces(),
                                alpha * u.left() * u.left().tr()  * nv(G).norm(),
                               -alpha * u.right()* u.left().tr()  * nv(G).norm(),
                               -alpha * u.left() * u.right().tr() * nv(G).norm(),
                                alpha * u.right()* u.right().tr() * nv(G).norm());
            const gsSparseMatrix<real_t> K2 = S.matrix(n);
            const gsMatrix<real_t> b2 = S.vector(n);

            CHECK( (K - K2).norm() <= 1e-12 * K.norm() );
            CHECK( (b - b2).norm() <= 1e-12 * b.norm() );

            // every assembled entry must be in the pattern
            A.computePattern_into(S, igrad(u, G) * igrad(u, G).tr());
            A.computePatternBdr_into(S, bc.get("Neumann"), u * u.tr());
            A.computePatternIfc_into(S, mp.interfaces(),
                                     u.left() * u.left().tr(),  u.right()* u.left().tr(),
                                     u.left() * u.right().tr(), u.right()* u.right().tr());
            index_t missing = 0;
            for (index_t k = 0; k != K2.outerSize(); ++k)
                for (gsSparseMatrix<real_t>::InnerIterator it(K2, k); it; ++it)
                    if (0 != it.value() && !S.pattern.count(std::make_pair(it.row(), it.col())))
                        ++missing;
            CHECK_EQUAL(0, missing);

            // the internal interface pattern covers the couplings across
            // the interface: the assembly does not add any entry to it
            gsExprAssembler<> B(1,1);
            B.setIntegrationElements(mb);
            gsExprAssembler<>::geometryMap G2 = B.getMap(mp);
            gsExprAssembler<>::space v = B.getSpace(mb);
            v.setup(bc, dirichlet::interpolation, icont);
            B.initSystem();
            B.computePatternIfc(mp.interfaces(),
                                v.left() * v.left().tr(),  v.right()* v.left().tr(),
                                v.left() * v.right().tr(), v.right()* v.right().tr());
            const index_t nnz = B.fiberMatrix().nonZeros();
            B.assembleIfc(mp.interfaces(),
                           alpha * v.left() * v.left().tr()  * nv(G2).norm(),
                          -alpha * v.right()* v.left().tr()  * nv(G2).norm(),
                          -alpha * v.left() * v.right().tr() * nv(G2).norm(),
                           alpha * v.right()* v.right().tr() * nv(G2).norm());
            CHECK_EQUAL(nnz, B.fiberMatrix().nonZeros());
        }
    }
}
