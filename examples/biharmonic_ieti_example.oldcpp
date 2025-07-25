/** @file biharmonic_ieti_example.cpp

    @brief Biharmonic example for an extremely simple multipatch domain

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <ctime>
#include <gismo.h>
#include <gsSolver/gsLowRankCorrectedOp.h>


using namespace gismo;

gsDofMapper setupTwoLayerDofMapper(const gsMultiPatch<>& mp, const gsMultiBasis<>& mb, unsigned bdyConds)
{
    gsVector<index_t> patchDofSizes(mp.nPatches());
    for (size_t k=0; k<mp.nPatches(); ++k)
        patchDofSizes[k] = mb[k].size();

    gsDofMapper dm(patchDofSizes);
    for (gsBoxTopology::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const index_t k1 = it->first().patch;
        const index_t k2 = it->second().patch;
        gsVector<index_t> s1 = mb.basis(k1).boundary(it->first().side()),
                          s2 = mb.basis(k2).boundary(it->second().side()),
                          s1o = mb.basis(k1).boundaryOffset(it->first().side(),1),
                          s2o = mb.basis(k2).boundaryOffset(it->second().side(),1);

        // TODO: We assume for now that the orientation matches!
        GISMO_ASSERT( s1.rows() == s2.rows() && s1.rows() == s1o.rows() && s2.rows() == s2o.rows(), "");
        for (index_t i=0;i<s1.rows();++i)
        {
            dm.matchDof(k1,s1[i],k2,s2[i]);
            dm.matchDof(k1,s1o[i],k2,s2o[i]);
        }
    }

    if (bdyConds)  // Eliminate boundary layer if bdyConds == 1 or == 2
    {
        for (gsBoxTopology::const_biterator it = mp.bBegin(); it != mp.bEnd(); ++it)
        {
            const index_t k = it->patchIndex();
            gsVector<index_t> s = mb.basis(k).boundary(it->side()),
                              so = mb.basis(k).boundaryOffset(it->side(),1);

            GISMO_ASSERT( s.rows() == so.rows(), "");
            for (index_t i=0;i<s.rows();++i)
            {
                dm.eliminateDof(s[i],k);
                if (bdyConds == 2)
                    dm.eliminateDof(so[i],k);
            }
        }
    }

    dm.finalize();
    return dm;
}

std::pair< gsMatrix<>, gsMatrix<> >
sampleBoundary(const gsGeometry<>& geo, boxSide s, bool inverted, index_t numberSamples)
{
    const short_t dim = 2;

    GISMO_ASSERT( s.index()>0 && s.index()<=4, "Invalid boxSide." );
    const short_t dir = (s.index()-1)/2;
    const short_t prm = (s.index()-1)%2;
    const gsMatrix<>& parameterRange = geo.parameterRange();

    gsMatrix<> sample(dim,numberSamples);
    for (index_t i=0; i<numberSamples; ++i)
    {
        index_t j = inverted ? (numberSamples-1-i) : i;
        sample(  dir, j) = parameterRange(  dir,0) + prm                     * (parameterRange(  dir,1)-parameterRange(  dir,0));
        sample(1-dir, j) = parameterRange(1-dir,0) + i/(numberSamples - 1.0) * (parameterRange(1-dir,1)-parameterRange(1-dir,0));
    }

    gsMatrix<> selector(2,4);
    selector.setZero();
    selector(0,dir)   = prm?1.0:-1.0;
    selector(1,dir+2) = prm?1.0:-1.0;

    return std::pair< gsMatrix<>, gsMatrix<> >( geo.eval(sample), selector*geo.deriv(sample) );
}


bool
checkC1(const gsMultiPatch<>& mp)
{
    bool isC1 = true;
    for (gsBoxTopology::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const index_t k1 = it->first().patch;
        const index_t k2 = it->second().patch;

        gsVector<index_t> cm;
        it->cornerMap(cm);
        GISMO_ASSERT (cm.size()==2 && cm[1] == 1-cm[0] && (cm[0] == 0 || cm[0] == 1), "Unexcpected corner map:\n"<<cm);
        const bool inverted = cm[0];

        std::pair< gsMatrix<>, gsMatrix<> > data1 = sampleBoundary(mp[k1], it->first(),  inverted, 11);
        std::pair< gsMatrix<>, gsMatrix<> > data2 = sampleBoundary(mp[k2], it->second(), false,    11);

        if ((data1.first-data2.first).cwiseAbs().maxCoeff()>1e-6)
        {
            gsDebug << "Interface between patches " << k1 << " and " << k2 << " is not C0.\n";
            isC1 = false;
        }
        if ((data1.second+data2.second).cwiseAbs().maxCoeff()>1e-6)
        {
            gsDebug << "Interface between patches " << k1 << " and " << k2 << " is not C1.\n";
            isC1 = false;
        }

    }
    return isC1;
}


std::vector<std::vector<std::vector<index_t>>>
setupTwoLayerSkeletonDofs(const gsMultiPatch<>& mp, const gsMultiBasis<>& mb, int bdyCond)
{
    GISMO_ENSURE ( bdyCond==0, "This function is not bounary-aware.");
    std::vector<std::vector<std::vector<index_t>>> result(mp.nPatches());
    for (size_t k=0; k<mp.nPatches(); ++k)
        result[k].resize(2);

    for (gsBoxTopology::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const index_t k1 = it->first().patch;
        const index_t k2 = it->second().patch;
        gsVector<index_t> s1 = mb.basis(k1).boundary(it->first().side()),
                          s2 = mb.basis(k2).boundary(it->second().side()),
                          s1o = mb.basis(k1).boundaryOffset(it->first().side(),1),
                          s2o = mb.basis(k2).boundaryOffset(it->second().side(),1);

        // We assume for now that the orientation matches!
        GISMO_ASSERT( s1.rows() == s2.rows() && s1.rows() == s1o.rows() && s2.rows() == s2o.rows(), "");
        for (index_t i=0;i<s1.rows();++i)
        {
            result[k1][0].push_back(s1 [i]);
            result[k1][1].push_back(s1o[i]);
            result[k2][0].push_back(s2 [i]);
            result[k2][1].push_back(s2o[i]);
        }
    }

    for (size_t k=0; k<mp.nPatches(); ++k)
    {
        for (size_t j=0; j<2; ++j)
        {
            std::sort(result[k][j].begin(),result[k][j].end());
            result[k][j].erase( std::unique(result[k][j].begin(),result[k][j].end()), result[k][j].end() );
        }
        /*
        std::vector<index_t> tmp;
        std::set_difference(result[k][1].begin(),result[k][1].end(),
                            result[k][0].begin(),result[k][0].end(),
                            std::inserter(tmp, tmp.begin()));
        tmp.swap(result[k][1]);
        //*/
    }

    return result;
}

struct corner_t {
    std::vector<std::pair<index_t,index_t>> data;
    void push_back(std::pair<index_t,index_t> p) { data.push_back(p); }
    bool operator==(const corner_t& other) const
    {
        if (data.size()!=other.data.size())
            return false;
        for (size_t i=0;i<data.size();++i)
        {
            int found = 0;
            for (size_t j=0;j<other.data.size();++j)
                if (data[i].first == other.data[j].first && data[i].second == other.data[j].second)
                    found += 1;
            if (found!=1)
                return false;
        }
        return true;
    }

};


std::vector<std::vector<std::pair<index_t,gsSparseVector<>>>>
cornersFromJumpMatrices( const std::vector<gsSparseMatrix<real_t,RowMajor>>& sms )
{
    // This function guesses the corners from the jump matrices; each dof belonging
    // to more than 2 patches is a corner.

    // lmultsof[k][i] ... which L-multipliers act on dof i of patch k?
    std::vector<std::vector<std::vector<index_t>>> lmultsof(sms.size());
    // dofsof[l][0 and 1] ... a pair (patch, index) of dofs connected by L-multiplier l
    std::vector<std::vector<std::pair<index_t,index_t>>> dofsof(sms[0].rows());
    for (size_t k=0; k<sms.size(); ++k)
    {
        lmultsof[k].resize(sms[k].cols());
        for (index_t i=0; i<sms[k].outerSize(); ++i)
            for (gsSparseMatrix<real_t,RowMajor>::InnerIterator it(sms[k],i); it; ++it)
            {
                lmultsof[k][it.col()].push_back(it.row());
                dofsof[it.row()].push_back(std::pair<index_t,index_t>(k, it.col()));
            }
    }
    /*gsInfo << "lmultsof: \n";
    for (size_t k=0;k<lmultsof.size(); ++k)
    {
        gsInfo << "[" << k << ":";
        for (size_t i=0;i<lmultsof[k].size(); ++i)
            if (lmultsof[k][i].size()>0)
            {
                gsInfo << "{" << i << ":";
                for (size_t j=0;j<lmultsof[k][i].size(); ++j)
                    gsInfo << " " << lmultsof[k][i][j];
                gsInfo << "}";
            }
        gsInfo << "]\n";
    }
    gsInfo << "dofsof: \n";
    for (size_t k=0;k<dofsof.size(); ++k)
    {
        gsInfo << "[" << k << ":";
        for (size_t i=0;i<dofsof[k].size(); ++i)
            gsInfo << " " << dofsof[k][i].first << "/" << dofsof[k][i].second;
        gsInfo << "]\n";
    }
    //*/


    std::vector<corner_t> result;
    for (size_t k=0; k<lmultsof.size(); ++k)
        for (size_t i=0; i<lmultsof[k].size(); ++i)
            if (lmultsof[k][i].size()>1)
            {
                corner_t corner;
                corner.push_back(std::pair<index_t,index_t>(k,i));
                for (size_t l=0; l<lmultsof[k][i].size(); ++l)
                {
                    index_t multiplier = lmultsof[k][i][l];
                    GISMO_ENSURE(dofsof[multiplier].size()==2,"");
                    if (dofsof[multiplier][0].first != (index_t)k)
                        corner.push_back(dofsof[multiplier][0]);
                    else if (dofsof[multiplier][1].first != (index_t)k)
                        corner.push_back(dofsof[multiplier][1]);
                    else
                    {
                        GISMO_ENSURE(0, "");
                    }
                }
                if (find(result.begin(), result.end(), corner) == result.end())
                    result.push_back(corner);
            }

    std::vector<std::vector<std::pair<index_t,gsSparseVector<>>>> finalResult;
    for (size_t i=0; i<result.size(); ++i)
    {
        std::vector<std::pair<index_t,gsSparseVector<>>> corner;
        for (size_t j=0; j<result[i].data.size(); ++j)
        {
            gsSparseVector<> sv(sms[result[i].data[j].first].cols());
            sv.setZero();
            //gsInfo << i << "%%" << j << "%%" <<
            //    result[i].data[j].first << "%%" << result[i].data[j].second << "%%"
            //    << sms[result[i].data[j].first].cols() << "%%" << result[i].data[j].second << std::endl;
            sv[result[i].data[j].second] = 1;
            corner.push_back(std::pair<index_t,gsSparseVector<>>(result[i].data[j].first,sv));
        }
        finalResult.push_back(corner);
    }
    return finalResult;

}


gsSparseMatrix<> makeTransformer(const gsBasis<>& basis)
{
    // This function creates a sparse matix that changes the sign of the
    // the n-1 st row and column (1D). This is tensorized. This function's
    // result is put into applyDofMapperTwoSided to respect the boundary
    // conditions.
    const index_t d = basis.dim();
    gsSparseMatrix<> result;
    for (index_t i=0; i<d; ++i)
    {
        const index_t ndofs1D = basis.component(d-1-i).size();
        GISMO_ASSERT( ndofs1D>3, "" );

        gsSparseMatrix<> transformer1D(ndofs1D,ndofs1D);
        transformer1D.setIdentity();

        transformer1D(0,0)=1;
        transformer1D(0,1)=1;

        transformer1D(1,1)=1;

        transformer1D(ndofs1D-2,ndofs1D-2)=-1;

        transformer1D(ndofs1D-1,ndofs1D-2)=1;
        transformer1D(ndofs1D-1,ndofs1D-1)=1;

        if (i==0)
            result = give(transformer1D);
        else
            result = result.kron(transformer1D);
    }
    return result;
}

gsSparseMatrix<>
applyDofMapperTwoSided(const gsSparseMatrix<>&sm, const gsDofMapper& dm)
{
    gsSparseEntries<> se;
    se.reserve(sm.nonZeros());
    for ( index_t i=0; i<sm.outerSize(); ++i )
        for ( gsSparseMatrix<>::InnerIterator it(sm,i); it; ++it )
            if (dm.is_free(it.row(), 0) && dm.is_free(it.col(), 0))
                se.add( dm.index(it.row(), 0), dm.index(it.col(), 0), it.value() );

    gsSparseMatrix<> result(dm.freeSize(), dm.freeSize());
    result.setFrom(se);
    return result;
}

class deluxePreconder {

    static size_t idx_of(std::vector<size_t>& vec, size_t what)
    {
        for (size_t i=0; i<vec.size(); ++i)
            if (vec[i]==what)
                return i;
        vec.push_back(what);
        return vec.size()-1;
    }

public:
    struct edge {
        size_t first;
        size_t second;
        std::vector<size_t> firstIndices;
        std::vector<size_t> secondIndices;
        std::vector<size_t> lagrangeIndices;
        edge(size_t _first, size_t _second) : first(_first), second(_second) {}
    };

    std::vector<edge> identifyEdges() const
    {
        // TODO: handle corner dofs!
        const size_t nPatches = m_localJumpMatrices.size();
        gsMatrix<index_t> dofslist;
        dofslist.setZero(m_localJumpMatrices[0].rows(), 2);
        index_t edgeCount = 0;

        std::vector<edge> edges;

        for (size_t k=0; k<nPatches; ++k)
        {
            std::vector<size_t> neighbors;

            for (index_t i=0; i<m_localJumpMatrices[k].outerSize(); ++i)
            {
                for (gsSparseMatrix<real_t,RowMajor>::InnerIterator it(m_localJumpMatrices[k],i); it; ++it)
                {
                    if (dofslist( it.row(), 0 )==0)
                    {
                        dofslist( it.row(), 0 ) = k+1;
                        dofslist( it.row(), 1 ) = it.col();
                        GISMO_ASSERT( (it.value()>.99 && it.value()<1.01) || (-it.value()>.99 && -it.value()<1.01), it.value());
                    }
                    else if (dofslist( it.row(), 0 )>0)
                    {
                        GISMO_ASSERT( (it.value()>.99 && it.value()<1.01) || (-it.value()>.99 && -it.value()<1.01), it.value());
                        size_t l = dofslist(it.row(),0)-1;
                        size_t e = idx_of(neighbors,l);
                        e+=edgeCount;
                        if (e==edges.size())
                            edges.emplace_back(l,k);
                        GISMO_ASSERT (e<edges.size(), "Internal error.");
                        edges[e].firstIndices.push_back(dofslist(it.row(),1));
                        edges[e].secondIndices.push_back(it.col());
                        edges[e].lagrangeIndices.push_back(it.row());
                        dofslist(it.row(),0) = -1; // mark as "done"
                    }
                    else
                        GISMO_ENSURE(0, "No L-multiplier should connect to 3 patches!");
                }
            }
            edgeCount += neighbors.size();
        }
        gsInfo << "[deluxe-pc: " << edgeCount << "edges]" << std::flush;
        /*
        for (size_t l=0; l<edges.size(); ++l)
        {
            gsInfo << "Edge " << l << ": ";
            gsInfo << "(" << edges[l].first << "," << edges[l].second << ") [ ";
            GISMO_ASSERT(edges[l].firstIndices.size()==edges[l].secondIndices.size(), "Internal error.");
            for (size_t i=0; i<edges[l].firstIndices.size(); ++i)
                gsInfo << "(" << edges[l].firstIndices[i] << "," << edges[l].secondIndices[i] << "," << edges[l].lagrangeIndices[i] << ") ";
            gsInfo << "]\n";
        }
        //*/


        return edges;

    }

    gsSparseMatrix<> setupLocalProblem(const edge& e) const
    {
        const size_t firstPatch = e.first;
        const size_t secondPatch = e.second;
        const index_t firstSize = m_localSystems[firstPatch].rows();
        const index_t secondSize = m_localSystems[secondPatch].rows();
        const index_t lagrangeSize = e.lagrangeIndices.size();
        const index_t localDirichletSize = m_localSkeletonDofs[e.first].size() + m_localSkeletonDofs[e.second].size() - 2*lagrangeSize;

        const index_t size = firstSize+secondSize+lagrangeSize+localDirichletSize;

        gsSparseEntries<> se;
        // first matrix
        for (index_t i=0; i<m_localSystems[firstPatch].outerSize(); ++i)
            for (gsSparseMatrix<real_t>::InnerIterator it(m_localSystems[firstPatch],i); it; ++it)
                se.add(it.row(), it.col(), it.value());
        // second matrix
        for (index_t i=0; i<m_localSystems[secondPatch].outerSize(); ++i)
            for (gsSparseMatrix<real_t>::InnerIterator it(m_localSystems[secondPatch],i); it; ++it)
                se.add(firstSize+it.row(), firstSize+it.col(), it.value());
        // jump blocks
        for (index_t i=0; i<lagrangeSize; ++i)
        {
            real_t value1 = m_localJumpMatrices[firstPatch ](e.lagrangeIndices[i], e.firstIndices [i]);
            real_t value2 = m_localJumpMatrices[secondPatch](e.lagrangeIndices[i], e.secondIndices[i]);
            GISMO_ASSERT( (value1>.99 && value1<1.01) || (-value1>.99 && -value1<1.01), "");
            GISMO_ASSERT( (value2>.99 && value2<1.01) || (-value2>.99 && -value2<1.01), "");
            se.add(firstSize+secondSize+i,             e.firstIndices[i] , value1 );
            se.add(firstSize+secondSize+i, firstSize + e.secondIndices[i], value2 );
            se.add(            e.firstIndices[i] , firstSize+secondSize+i, value1 );
            se.add(firstSize + e.secondIndices[i], firstSize+secondSize+i, value2 );
        }

        std::vector<index_t> outerSkeletonFirst;
        std::copy_if(m_localSkeletonDofs[e.first].begin(), m_localSkeletonDofs[e.first].end(), std::back_inserter(outerSkeletonFirst),
            [&](index_t arg)
            { return (std::find(e.firstIndices.begin(), e.firstIndices.end(), arg) == e.firstIndices.end());});

        std::vector<index_t> outerSkeletonSecond;
        std::copy_if(m_localSkeletonDofs[e.second].begin(), m_localSkeletonDofs[e.second].end(), std::back_inserter(outerSkeletonSecond),
            [&](index_t arg)
            { return (std::find(e.secondIndices.begin(), e.secondIndices.end(), arg) == e.secondIndices.end());});

        GISMO_ASSERT( outerSkeletonFirst.size()+outerSkeletonSecond.size() == localDirichletSize, "Internal error.");
        for (size_t i=0; i<outerSkeletonFirst.size(); ++i)
        {
            se.add( firstSize+secondSize+lagrangeSize+i, outerSkeletonFirst[i], 1 );
            se.add( outerSkeletonFirst[i], firstSize+secondSize+lagrangeSize+i, 1 );
        }
        for (size_t i=0; i<outerSkeletonSecond.size(); ++i)
        {
            se.add( firstSize+secondSize+lagrangeSize+outerSkeletonFirst.size()+i, firstSize+outerSkeletonSecond[i], 1 );
            se.add( firstSize+outerSkeletonSecond[i], firstSize+secondSize+lagrangeSize+outerSkeletonFirst.size()+i, 1 );
        }

        // setup matrix
        gsSparseMatrix<> result(size,size);
        result.setFrom(se);
        result *= -1; // Schur complement has form C - B^t A^{-1} B, so must be negative!

        //gsInfo << "***** local edge system *****\n\n" << result << "\n*****\n";

        return result;
    }

    gsSparseMatrix<> getLocalEmbedding(const edge& e)
    {
        const size_t firstPatch = e.first;
        const size_t secondPatch = e.second;
        const index_t firstSize = m_localSystems[firstPatch].rows();
        const index_t secondSize = m_localSystems[secondPatch].rows();
        const index_t lagrangeSize = e.lagrangeIndices.size();
        const index_t localDirichletSize = m_localSkeletonDofs[e.first].size() + m_localSkeletonDofs[e.second].size() - 2*lagrangeSize;
        const index_t localSize = firstSize+secondSize+lagrangeSize+localDirichletSize;
        const index_t globalLagrangeSize = m_localJumpMatrices[0].rows();

        gsSparseEntries<> se;
        se.reserve(lagrangeSize);

        for (index_t i=0; i<lagrangeSize; ++i)
            se.add( e.lagrangeIndices[i], firstSize+secondSize+i, 1 );

        // setup matrix
        gsSparseMatrix<> result(globalLagrangeSize,localSize);
        result.setFrom(se);
        return result;
    }


public:

    void addSubdomain(
            gsSparseMatrix<> localSystem,
            gsSparseMatrix<real_t,RowMajor> localJumpMatrix,
            std::vector<index_t> localSkeletonDofs
    )
    {
        m_localSystems.push_back(give(localSystem));
        m_localJumpMatrices.push_back(give(localJumpMatrix));
        m_localSkeletonDofs.push_back(give(localSkeletonDofs));
    }


    void setupPreconditioner()
    {
        m_preconditioner = gsSumOp<>::make();
        std::vector<edge> edges = identifyEdges();
        for (size_t e=0; e<edges.size(); ++e)
        {
            gsLinearOperator<>::Ptr solver = makeSparseLUSolver(setupLocalProblem(edges[e]));
            memory::shared_ptr<gsSparseMatrix<>> embedding = getLocalEmbedding(edges[e]).moveToPtr();
            gsLinearOperator<>::Ptr prolong = makeMatrixOp(embedding.get()->transpose());
            gsLinearOperator<>::Ptr restrict = makeMatrixOp(embedding);
            m_preconditioner->addOperator( gsProductOp<>::make(prolong, solver, restrict) );
        }
    }

    gsLinearOperator<>::Ptr preconditioner()
    {
        return m_preconditioner;
    }


private:
    std::vector<gsSparseMatrix<>> m_localSystems;
    std::vector<gsSparseMatrix<real_t,RowMajor>> m_localJumpMatrices;
    std::vector<std::vector<index_t>> m_localSkeletonDofs;
    gsSumOp<>::Ptr m_preconditioner;


};


int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    index_t rhsType = 2;
    real_t sigma = .5;
    index_t nPatchesX = 4;
    index_t nPatchesY = 4;
    index_t geoExample = 0;
    index_t degree = 2;
    index_t multiplicity = 1;
    index_t refinements = 2;
    real_t robin = 0;
    real_t alpha = 1;
    int bdyConds = 2;
    std::string primals("x");
    std::string dualPreconder("d");
    real_t extremelyDeluxeParameter = 1;
    std::string solverType("cg");
    real_t tolerance = 1.e-8;
    index_t maxIterations = 1000;
    std::string out;
    bool plot = false;

    gsCmdLine cmd("Biharmonic IETI example for an extremely simple multipatch domain.");
    cmd.addInt   ("t", "RhsType",               "Chosen right-hand side", rhsType);
    cmd.addReal  ("s", "Sigma",                 "Poisson ratio", sigma);
    cmd.addInt   ("x", "PatchesX",              "Number of patches (coordinate direction x)", nPatchesX);
    cmd.addInt   ("y", "PatchesY",              "Number of patches (coordinate direction y)", nPatchesY);
    cmd.addInt   ("",  "GeoExample",            "How to represent the first patch?", geoExample);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addInt   ("m", "Multiplicity",          "Multiplicity of knots for B-spline discretization space", multiplicity);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addReal  ("o", "Robin",                 "Penalty parameter for Robin boundary conditions", robin);
    cmd.addReal  ("a", "Alpha",                 "Scaling parameter for reaction term", alpha);
    cmd.addInt   ("b", "BdyConds",              "Bounday conditions: (0) \u0394u, \u2202n\u0394u; (1) u, \u0394u; (2) u, \u2202nu", bdyConds);
    cmd.addString("c", "Primals",               "Chosen primal dofs: (0) no, (c) classical corners, (x) eXtended cornerdofs", primals);
    cmd.addString("d", "DualPreconder",         "Use preconder: (s) stanard Dirichlet, (c) componentwise Dirichlet, (d) deluxe, (e) extremely deluxe, (b) bruteforce", dualPreconder);
    cmd.addReal  ("",  "ExtremelyDeluxeParameter", "", extremelyDeluxeParameter);
    cmd.addString("",  "Solver",                "Which solver to use: \"cg\" or \"gmres\".", solverType);
    cmd.addReal  ("",  "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Stopping criterion for linear solver", maxIterations);
    cmd.addString("",  "out",                   "Write solution and used options to file", out);
    cmd.addSwitch(     "plot",                  "Plot the result with Paraview", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Run biharmonic_ieti_example with options:\n" << cmd << "\n";

    /********************** Define rhs **********************/
    const index_t dim = 2;
    const char *rhsTypes[] =
        {
            "32*pi^4*sin(2*pi*x)*sin(2*pi*y)",
            "2*pi^4*sin(pi*x)*sin(pi*y)",
            "1/8*pi^4*sin(pi*x/2)*sin(pi*y/2)",
            "2/100*pi^4*sin(pi*x/10)*sin(pi*y/10)"
        };
    if (rhsType < 0 || (size_t)rhsType >= util::size(rhsTypes))
    {
        gsInfo << "Invalid choice for --RhsType (-t).\n";
        return -1;
    }
    if (bdyConds < 0 || bdyConds > 2)
    {
        gsInfo << "Invalid choice for --BdyCond (-b).\n";
        return -1;
    }

    gsInfo << "Rhs function is " << rhsTypes[rhsType] << "\n";
    gsFunctionExpr<> f(rhsTypes[rhsType],dim);

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;
    gsMultiPatch<> mp;
    /// cases
    if (geoExample==0)
        mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,1,1));
    else if (geoExample==1)
        mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,-1,1,0,90));
    else if (geoExample==2)
        mp.addPatch(gsNurbsCreator<>::BSplineRectangle(-1,-1,0,0,180));
    else if (geoExample==3)
        mp.addPatch(gsNurbsCreator<>::BSplineRectangle(-1,0,0,1,270));
    else if (geoExample==4)
        mp.addPatch(gsNurbsCreator<>::BSplineRectangle(1,0,0,1));
    else if (geoExample==5)
        mp.addPatch(gsNurbsCreator<>::BSplineRectangle(1,1,0,0));
    else
    {
        GISMO_ENSURE(false, "Invalid geoExample");
    }

    /// end
    for (index_t i=0; i<nPatchesX; ++i)
        for (index_t j=0; j<nPatchesY; ++j)
            if (i+j)
                mp.addPatch(gsNurbsCreator<>::BSplineRectangle(i,j,i+1,j+1));
    mp.computeTopology();

    const index_t nPatches = mp.nPatches();
    gsInfo << "done: " << nPatches << " patches, " << mp.interfaces().size() << " interfaces.\n";

    /************ Setup bases and adjust degree *************/

    gsMultiBasis<> mb(mp);
    gsInfo << "Setup bases and adjust degree... " << std::flush;
    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    for ( index_t i = 0; i < refinements; ++i )
        mb.uniformRefine();

    for ( size_t i = 0; i < mb.nBases(); ++ i )
    {
        if (multiplicity>1 && multiplicity<degree)
            mb[i].reduceContinuity(multiplicity-1);
        else if (multiplicity!=1)
        {
            gsInfo << "Multiplicity must be at least 1 and at most degree-1\n";
            return -1;
        }
    }

    gsInfo << "done.\n";

    /******************* Setup dofMapper ********************/
    if (!checkC1(mp))
    {
        gsInfo << "This is not a C1 geometry.\n";
        return -1;
    }

    gsInfo << "Setup dofMapper... " << std::flush;
    gsDofMapper dm = setupTwoLayerDofMapper(mp, mb, robin==0.?bdyConds:0);
    std::vector<std::vector<std::vector<index_t>>> skeletonDofs;
    if (dualPreconder=="c")
        skeletonDofs = setupTwoLayerSkeletonDofs(mp, mb, robin==0.?bdyConds:0);
    gsInfo << "done:\n" << dm << "\n";

    /****************** Setup ietimapper ********************/
    gsInfo << "Setup ietimapper... " << std::flush;
    gsMatrix<> fixedPart;
    fixedPart.setZero(dm.boundarySize(),1);
    gsIetiMapper<> ietiMapper(mb,dm,fixedPart);
    ietiMapper.computeJumpMatrices(/*fullyRedundant=*/false,/*excludeCorners=*/false,/*excludeDofsForSeveralPatches=*/true);
    gsIetiMapper<> ietiMapper2(mb,dm,fixedPart);
    ietiMapper2.computeJumpMatrices(/*fullyRedundant=*/true,/*excludeCorners=*/false,/*excludeDofsForSeveralPatches=*/false);

    if (primals == "c")
    {
        ietiMapper.cornersAsPrimals();
        gsInfo << "[" << ietiMapper.nPrimalDofs() << " classical cornerdofs added as primals]";
    }
    else if (primals == "x")
    {

        std::vector<gsSparseMatrix<real_t,RowMajor>> jm;
        jm.reserve(nPatches);
        for (index_t k=0; k<nPatches; ++k)
            jm.push_back(ietiMapper2.jumpMatrix(k));
        //gsInfo << "im2-jumps[0]:" << jm[0].rows() << "x" << jm[0].cols() << "\n";
        std::vector<std::vector<std::pair<index_t,gsSparseVector<>>>> data = cornersFromJumpMatrices(jm);

        if (! true )
        {
            gsInfo << "[";
            for (size_t i=0; i<data.size(); ++i)
            {
                gsInfo << "[";
                for (size_t j=0; j<data[i].size(); ++j)
                {
                    gsInfo << " {" << data[i][j].first << "," << gsSparseVector<>::InnerIterator(data[i][j].second,0).row() << "}";
                }
                gsInfo << " ]\n";
            }
            gsInfo << "]\n\n";
        }

        for (size_t i=0; i<data.size(); ++i)
            ietiMapper.customPrimalConstraints(data[i]);
        gsInfo << "[" << data.size() << " eXtended cornerdofs added as primals]";
    }
    else if (primals == "0")
        gsInfo << "[no primaldofs added]";
    else
    {
        gsInfo << "Invalid choice for --Primals.\n";
        return -1;
    }
    gsIetiSystem<> ieti;
    ieti.reserve(nPatches+1);

    gsScaledDirichletPrec<> prec;
    prec.reserve(nPatches);

    gsPrimalSystem<> primal(ietiMapper.nPrimalDofs());

    gsInfo << "done.\n";

    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrix... " << std::flush;

    std::vector<gsSparseMatrix<>> localBasisTransforms; localBasisTransforms.reserve(nPatches);
    std::vector<gsSparseMatrix<>> localStiffnessMatrices; localStiffnessMatrices.reserve(nPatches);
    std::vector<gsMatrix<>> localRhsVectors; localRhsVectors.reserve(nPatches);
    std::vector<gsLinearOperator<>::Ptr> localSchurs; localSchurs.reserve(nPatches);

    deluxePreconder deluxe;

    for (index_t k=0; k<nPatches; ++k)
    {
        gsInfo << "[" << k << "] " << std::flush;
        gsMultiPatch<> mp_local(mp[k]);
        gsMultiBasis<> mb_local(mb[k]);

        //! [Problem setup]
        gsExprAssembler<real_t> A(1,1);

        // Elements used for numerical integration
        A.setIntegrationElements(mb_local);
        gsExprEvaluator<real_t> ev(A);

        // Set the geometry map
        auto G = A.getMap(mp_local);
        auto ff = A.getCoeff(f, G);
        auto u = A.getSpace(mb_local);

        u.setup();
        ietiMapper.initFeSpace(u,k);

        A.initSystem();
        A.assemble(sigma *ilapl(u, G) * ilapl(u, G).tr() * meas(G), u * ff * meas(G));
        A.assemble((1-sigma) * ihess(u, G) * ihess(u, G).tr() * meas(G));
        A.assemble(alpha*u*u.tr()*meas(G));


        gsSparseMatrix<> transformer = applyDofMapperTwoSided(makeTransformer(mb[k]),ietiMapper.dofMapperLocal(k));

        // Fetch data
        gsSparseMatrix<real_t, RowMajor> jumpMatrix  = ietiMapper.jumpMatrix(k);
        gsSparseMatrix<>                 localMatrix = transformer*A.matrix()*transformer.transpose();
        gsMatrix<>                       localRhs    = transformer*A.rhs();

        //gsInfo << "\nlocalMatrix:"<<localMatrix.rows()<<"x"<<localMatrix.cols()<<";jumpMatrix:"<<jumpMatrix.rows()<<"x"<<jumpMatrix.cols();

        GISMO_ASSERT(jumpMatrix.cols() == localMatrix.rows(), "");

        // Penalize (Dirichlet) boundary
        if (robin)
        {
            gsInfo << "[robin]";
            for (gsBoxTopology::const_biterator it = mp.bBegin(); it != mp.bEnd(); ++it)
            {
                if (it->patchIndex()==k)
                {
                    // gsInfo << "Found biterator: " << *it << "\n";
                    gsVector<index_t> s1 = mb_local.basis(0).boundary(it->side());
                    for (index_t i=0; i<s1.rows(); ++i)
                        localMatrix(s1[i],s1[i]) += robin;
                }
            }
        }

        // Store
        localStiffnessMatrices.push_back(localMatrix);
        localRhsVectors.push_back(localRhs);
        localBasisTransforms.push_back(transformer);
        localSchurs.push_back(gsScaledDirichletPrec<>::schurComplement(localMatrix,ietiMapper.skeletonDofs(k)));

        if (dualPreconder=="s")
        {
            prec.addSubdomain(
                gsScaledDirichletPrec<>::restrictToSkeleton(
                    jumpMatrix,
                    localMatrix,
                    ietiMapper.skeletonDofs(k)
                )
            );
        }
        else if (dualPreconder=="c")
        {
            for (index_t j=0; j<2; ++j)
                prec.addSubdomain(
                    gsScaledDirichletPrec<>::restrictToSkeleton(
                        jumpMatrix,
                        localMatrix,
                        skeletonDofs[k][j]
                    )
                );
        }
        else if (dualPreconder=="d" || dualPreconder=="e" || dualPreconder=="b")
        {
#ifndef NDEBUG
            gsInfo << "Skeleton=[";
            for (size_t i=0; i<ietiMapper.skeletonDofs(k).size(); ++i)
                gsInfo  << ietiMapper.skeletonDofs(k)[i] << " ";
            gsInfo << "]";
#endif
            deluxe.addSubdomain(localMatrix, jumpMatrix, ietiMapper.skeletonDofs(k));
        }
        else
        {
            gsInfo << "\nUnknown dual preconder.\n";
            return -1;
        }

        // This function writes back to jumpMatrix, localMatrix, and localRhs,
        // so it must be called after prec.addSubdomain().
        //! [Patch to primals]
        primal.handleConstraints(
            ietiMapper.primalConstraints(k),
            ietiMapper.primalDofIndices(k),
            jumpMatrix,
            localMatrix,
            localRhs
        );
        //! [Patch to primals]

        // Add the patch to the Ieti system
        //! [Patch to system]
        ieti.addSubdomain(
            jumpMatrix.moveToPtr(),
            makeMatrixOp(localMatrix.moveToPtr()),
            give(localRhs)
        );
        //! [Patch to system]
    //! [End of assembling loop]
    } // end for
    //! [End of assembling loop]

    // Add the primal problem if there are primal constraints
    //! [Primal to system]
    if (ietiMapper.nPrimalDofs()>0)
    {
        gsInfo << "[P] " << std::flush;
        gsLinearOperator<>::Ptr localSolver = makeSparseLUSolver(primal.localMatrix());

        ieti.addSubdomain(
            gsSparseMatrix<real_t,RowMajor>(primal.jumpMatrix()).moveToPtr(),
            makeMatrixOp(gsSparseMatrix<>(primal.localMatrix()).moveToPtr()),
            give(primal.localRhs()),
            localSolver
        );
    }
    //! [Primal to system]

    gsInfo << "done. " << ietiMapper.nPrimalDofs() << " primal dofs.\n";


    /**************** Setup solver and solve ****************/

    gsInfo << "Setup solver and solve... \n";

    // Tell the preconditioner to set up the scaling
    //! [Setup scaling]
    gsLinearOperator<>::Ptr preconder;
    if (dualPreconder=="d")
    {
        gsInfo << "    Setup deluxe preconder... " << std::flush;
        deluxe.setupPreconditioner();
        preconder = deluxe.preconditioner();
    }
    else if (dualPreconder=="e")
    {
        gsInfo << "    Setup extremely deluxe preconder... " << std::flush;
        deluxe.setupPreconditioner();
        preconder = gsLowRankCorrectedOp<>::make( deluxe.preconditioner(), extremelyDeluxeParameter*primal.localMatrix(), primal.jumpMatrix(), primal.jumpMatrix() );
    }
    else if (dualPreconder=="b")
    {
        gsInfo << "    Setup bruteforce preconder... " << std::flush;
        gsSparseMatrix<> sm;
        {
            gsMatrix<> m;
            ieti.schurComplement()->toMatrix(m);
            sm = m.sparseView(m(0,0),1e-8);
        }
        gsSparseMatrix<> sm2(sm.rows(), sm.cols());
        std::vector<deluxePreconder::edge> edges = deluxe.identifyEdges();
        std::vector<bool> taken(sm.rows());
        for (size_t k=0; k<edges.size(); ++k)
        {
            gsSparseEntries<> se;
            se.reserve(edges[k].lagrangeIndices.size());
            for (size_t i=0; i<edges[k].lagrangeIndices.size(); ++i)
            {
                GISMO_ENSURE ( !taken[edges[k].lagrangeIndices[i]], "");
                taken[edges[k].lagrangeIndices[i]] = true;
                se.add(edges[k].lagrangeIndices[i],edges[k].lagrangeIndices[i],1.);
            }
            gsSparseMatrix<> trans(sm.rows(), sm.cols());
            trans.setFrom(se);
            sm2 += trans * sm * trans;
        }
        for (size_t i=0; i<taken.size(); ++i)
        {
            GISMO_ENSURE (taken[i],"i="<<i);
        }
        //preconder = makeSparseLUSolver(sm2);
        preconder = makeSparseCholeskySolver(sm2);
    }
    else if (dualPreconder=="c" || dualPreconder=="s")
    {
        gsInfo << "    Setup sclaled Dirichlet preconder... " << std::flush;
        prec.setupMultiplicityScaling();
        preconder = prec.preconditioner();
    }
    else
        GISMO_ENSURE (0, "Unknown dualPreconder.");
    //! [Setup scaling]

    gsInfo << "done.\n    Setup rhs... " << std::flush;
    // Compute the Schur-complement contribution for the right-hand-side
    //! [Setup rhs]
    gsMatrix<> rhsForSchur = ieti.rhsForSchurComplement();
    //! [Setup rhs]

    gsInfo << "done.\n    Setup cg solver for Lagrange multipliers and solve... " << std::flush;
    // Initial guess
    //! [Define initial guess]
    gsMatrix<> lambda;
    lambda.setRandom( ieti.nLagrangeMultipliers(), 1 );
    //! [Define initial guess]

    gsMatrix<> errorHistory;

    // This is the main cg iteration
    //! [Solve]

#ifndef NDEBUG
    {
        gsMatrix<> m;
        ieti.schurComplement()->toMatrix(m);
        m=m.unaryExpr([](double x){return (abs(x)<1e-4)?0.:x;});
        gsInfo << "\n\nProblem F:\n" << std::setprecision(3) <<  m << "\n\n";
    }
    {
        gsMatrix<> m;
        preconder->toMatrix(m);
        m=m.unaryExpr([](double x){return (abs(x)<1e-4)?0.:x;});
        gsInfo << "\n\nPreconder M:\n" << std::setprecision(3) <<  m << "\n\n";
    }
#endif

    real_t conditionNumber = -1;
    gsMatrix<real_t> eigenvalues;
    if (solverType == "cg")
    {
        gsConjugateGradient<> solver( ieti.schurComplement(), preconder );
        solver.setOptions( cmd.getGroup("Solver") );
        solver.setCalcEigenvalues(true);
        solver.solveDetailed( rhsForSchur, lambda, errorHistory );
        conditionNumber = solver.getConditionNumber();
        solver.getEigenvalues(eigenvalues);
    }
    else if (solverType == "gmres")
    {
        gsGMRes<> solver( ieti.schurComplement(), preconder );
        solver.setOptions( cmd.getGroup("Solver") );
        solver.solveDetailed( rhsForSchur, lambda, errorHistory );
    }
    else
    {
        gsInfo << "Unknown --solver.\n";
        return -1;
    }
    //! [Solve]

    gsInfo << "done.\n    Reconstruct solution from Lagrange multipliers... " << std::flush;
    // Now, we want to have the global solution for u
    //! [Recover]
    std::vector<gsMatrix<>> localSol = primal.distributePrimalSolution(
            ieti.constructSolutionFromLagrangeMultipliers(lambda)
        );
    gsMatrix<> globalSol = ietiMapper.constructGlobalSolutionFromLocalSolutions(localSol);
    //! [Recover]
    gsInfo << "done.\n\n";

    /******************** Print end Exit ********************/

    const index_t iter = errorHistory.rows()-1;
    const bool success = errorHistory(iter,0) < tolerance;
    if (success)
        gsInfo << "Reached desired tolerance after " << iter << " iterations:\n";
    else
        gsInfo << "Did not reach desired tolerance after " << iter << " iterations:\n";

    if (errorHistory.rows() < 20)
        gsInfo << errorHistory.transpose() << "\n\n";
    else
        gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

    if (solverType == "cg")
    {
        gsInfo << "Estimated condition number: " << conditionNumber << "\n";
        gsInfo << "Eigenvalues: " << eigenvalues.transpose() << "\n";
    }

    /********************** Output **************************/
    if (!out.empty())
    {
        /*gsFileData<> fd;
        std::time_t time = std::time(NULL);
        fd.add(cmd);
        for (index_t k=0; k<nPatches; ++k)
            fd.add(gsMatrix<>(localBasisTransforms[k].transpose() * localSol[k]));

        for (index_t k=0; k<nPatches; ++k)
            fd.add(gsSparseMatrix<>(ietiMapper.jumpMatrix(k)));

        for (index_t k=0; k<nPatches; ++k)
            fd.add(localStiffnessMatrices[k]);

        for (index_t k=0; k<nPatches; ++k)
        {
            gsMatrix<> tmp;
            localSchurs[k]->toMatrix(tmp);
            fd.add(tmp);
        }

        gsMatrix<> sc;
        ieti.schurComplement()->toMatrix(sc);
        fd.add(sc);

        gsMatrix<> pc;
        prec.preconditioner()->toMatrix(pc);
        fd.add(pc);

        fd.addComment(std::string("biharmonic_ieti_example   Timestamp:")+std::ctime(&time));
        fd.save(out);*/
        std::ofstream outfile (out, std::ios_base::app);
        outfile << "biharmonic_ieti_example\t"
                << degree << "\t"
                << refinements << "\t"
                << conditionNumber << "\t"
                << iter << "\t"
                << "***" << "\t"
                << nPatchesX << "\t"
                << nPatchesY << "\t"
                << rhsType << "\t"
                << robin << "\t"
                << alpha << "\t"
                << bdyConds << "\t"
                << primals << "\t"
                << dualPreconder << "\t"
                << extremelyDeluxeParameter << "\t"
                << solverType << "\n";

        gsInfo << "Write solution to file " << out << "\n";
    }

    if (plot)
    {
        gsInfo << "Write Paraview data to file ieti_result.pvd\n";
        gsMultiPatch<> mpsol;
        for (index_t k=0; k<nPatches; ++k)
            mpsol.addPatch( mb[k].makeGeometry( localBasisTransforms[k].transpose() * localSol[k] ) );
        gsWriteParaview<>( gsField<>( mp, mpsol ), "ieti_result", 1000);
    }
    if (!plot&&out.empty())
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution or --out to write solution to xml file.\n";
    }

    return 0;
}
