/** @file gsG1BiharmonicAssembler.h
    @brief Provides assembler for a homogenius Biharmonic equation.
    This file is part of the G+Smo library.
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    Author(s): P. Weinmueller
*/

#pragma once

#include <gsAssembler/gsAssembler.h>

#include <gsPde/gsBiharmonicPde.h>

#include <gsAssembler/gsVisitorBiharmonic.h>
#include <gsAssembler/gsVisitorNeumann.h>
#include <gsAssembler/gsVisitorNeumannBiharmonic.h>
//#include <gsAssembler/gsVisitorNitscheBiharmonic.h>

# include <gsG1Basis/gsG1Mapper_pascal.h>

# include <gsG1Basis/gsG1System.h>

namespace gismo
{

/** @brief
    Implementation of a homogeneous Biharmonic Assembler.
    It sets up an assembler and assembles the system patch wise and
    combines the patch-local stiffness matrices into a global system.
    Dirichlet boundary can only be enforced strongly (i.e Nitsche is
    not implemented).
*/
template <class T, class bhVisitor = gsVisitorBiharmonic<T> >
class gsG1BiharmonicAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
/** @brief
    Constructor of the assembler object.
    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] bconditions  is a gsBoundaryConditions object that holds boundary conditions on the form:
    \f[ \text{Dirichlet: } u = g \text{ on } \Gamma, \text{ and Neumann: } \nabla \Delta u \cdot \mathbf{n} = h \text{ on } \Gamma\f]
    \param[in] bconditions2 is a gsBoundaryConditions object that holds Neumann boundary conditions on the form:
    \f[\text{Neumann: } \nabla \Delta u \cdot \mathbf{n} = g\, \rightarrow \,(g,\nabla v \cdot \mathbf{n})_\Gamma, \f] where \f$ g \f$ is the Neumann data,
    \f$ v \f$ is the test function and \f$ \Gamma\f$ is the boundary side.
    \param[in] rhs is the right-hand side of the Biharmonic equation, \f$\mathbf{f}\f$.
    \param[in] dirStrategy option for the treatment of Dirichlet boundary in the \em bconditions object.
    \param[in] intStrategy option for the treatment of patch interfaces
*/
    gsG1BiharmonicAssembler(gsMultiPatch<T> const        & patches,
                            gsMultiBasis<T> const         & bases,
                            gsBoundaryConditions<T> const & bconditions,
                            gsBoundaryConditions<T> const & bconditions2,
                            const gsPiecewiseFunction<T>          & rhs)
        : m_ppde(patches,bconditions,bconditions2,rhs)
    {
        iFace::strategy intStrategy = iFace::none;
        dirichlet::strategy dirStrategy = dirichlet::none;

        m_options.setInt("DirichletStrategy", dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);

        Base::initialize(m_ppde, bases, m_options);
    }

    void refresh();

    void assemble();

    void constructSolution(const gsMatrix<T>& solVector,
                           gsMultiPatch<T>& result, int unk = 0);

    void computeDirichletDofsL2Proj(gsG1System<real_t> &  g1System);

    void plotParaview(gsField<> &solField_interior, std::vector<gsMultiPatch<T>> &result);



    gsMatrix<> get_bValue() { return m_g1_ddof; };

protected:

    // fixme: add constructor and remove this
    gsBiharmonicPde<T> m_ppde;

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

protected:

    gsMatrix<T> m_g1_ddof;

};

template<class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::constructSolution(const gsMatrix<T> &solVector, gsMultiPatch<T> &result, int unk)
{
    result.clear();
    const index_t dim = ( 0!=solVector.cols() ? solVector.cols() :  m_ddof[unk].cols() );

    gsMatrix<T> coeffs;
    index_t globalID = 0;
    for (size_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p)
    {
        // Reconstruct solution coefficients on patch p
        const int sz = m_bases[m_system.colBasis(unk)][p].size();

        coeffs.resize(sz, dim);
        coeffs = solVector.block(globalID,0,sz,1);

        globalID += sz;

        result.addPatch(m_bases[m_system.colBasis(unk)][p].makeGeometry(give(coeffs)));
    }
}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::plotParaview(gsField<> &solField_interior, std::vector<gsMultiPatch<T>> &g1Basis)
{

    std::string fn = "G1Biharmonic";
    index_t npts = 2000;
    gsParaviewCollection collection2(fn);
    std::string fileName2;

    for ( size_t pp =0; pp < m_pde_ptr->domain().nPatches(); ++pp ) // Patches
    {
        fileName2 = fn + util::to_string(pp);
        //writeSinglePatchField( field, i, fileName, npts );

        const gsFunction<T> & geometry = m_pde_ptr->domain().patch(pp);
        const gsFunction<T> & parField = solField_interior.function(pp);

        const int n = geometry.targetDim();
        const int d = geometry.domainDim();

        gsMatrix<T> ab = geometry.support();
        gsVector<T> a = ab.col(0);
        gsVector<T> b = ab.col(1);

        gsVector<unsigned> np = uniformSampleCount(a, b, npts);
        gsMatrix<T> pts = gsPointGrid(a, b, np);

        gsMatrix<T> eval_geo = geometry.eval(pts);//pts
        gsMatrix<T> eval_field = solField_interior.isParametric() ? parField.eval(pts) : parField.eval(eval_geo);

        // Here add g1 basis
        //eval_field.setZero();
        for (size_t i = 0; i < g1Basis[pp].nPatches(); i++)
        {
            gsField<> temp_field(m_pde_ptr->domain().patch(pp),g1Basis[pp].patch(i));
            eval_field += temp_field.value(pts);
        }


        gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
        //gsFunctionExpr<> solVal("y",2);
        gsField<> exact(m_pde_ptr->domain().patch(pp), solVal, false );
        //eval_field -= exact.value(pts);

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

        gsWriteParaviewTPgrid(eval_geo, eval_field, np.template cast<index_t>(), fileName2);


        collection2.addPart(fileName2, ".vts");
    }
    collection2.save();
}


template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::refresh()
{
    // We use predefined helper which initializes the system matrix
    // rows and columns using the same test and trial space
    Base::scalarProblemGalerkinRefresh();

}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::assemble()
{
    GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0], 2, 1, 0.333333);

    m_system.reserve(nz, this->pde().numRhs());

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    m_ddof.resize(m_system.numUnknowns());
    m_ddof[0].setZero(m_system.colMapper(0).boundarySize(), m_system.unkSize(0) * m_system.rhs().cols());


    // Assemble volume integrals
    Base::template push<bhVisitor>();

    // Neuman conditions of first kind
    //Base::template push<gsVisitorNeumann<T> >(
    //    m_ppde.bcFirstKind().neumannSides() );

    // Neuman conditions of second kind
    Base::template push<gsVisitorNeumannBiharmonic<T> >(m_ppde.bcSecondKind().neumannSides());


    if (m_options.getInt("InterfaceStrategy") == iFace::dg)
        gsWarn << "DG option ignored.\n";

    /*
    // If requested, force Dirichlet boundary conditions by Nitsche's method
    this->template push<gsVisitorNitscheBiharmonic<T> >(
    m_ppde.bcSecondKind().dirichletSides() );
    */

    // Assembly is done, compress the matrix
    Base::finalize();
}


template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::computeDirichletDofsL2Proj(gsG1System<real_t> &  g1System)
{
    gsVector<> numBoundaryVertexFunctions = g1System.get_numBoundaryVertexFunctions();
    gsVector<> numBoundaryEdgeFunctions = g1System.get_numBoundaryEdgeFunctions();

    size_t unk_ = 0;

    m_g1_ddof.resize( g1System.boundary_size(), m_system.unkSize(unk_)*m_system.rhs().cols());  //m_pde_ptr->numRhs() );
    m_g1_ddof.setZero();

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>        globProjRhs;
    globProjRhs.setZero( g1System.boundary_size(), m_system.unkSize(unk_)*m_system.rhs().cols() );

    // Temporaries
    gsVector<T> quWeights;

    gsMatrix<T> rhsVals;
    gsMatrix<unsigned> globIdxAct;
    gsMatrix<T> basisVals;

    gsMapData<T> md(NEED_MEASURE);

    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
              iter = m_pde_ptr->bc().dirichletBegin();
          iter != m_pde_ptr->bc().dirichletEnd(); ++iter )
    {
        if (iter->isHomogeneous() )
            continue;

        GISMO_ASSERT(iter->function()->targetDim() == m_system.unkSize(unk_)*m_system.rhs().cols(),
                     "Given Dirichlet boundary function does not match problem dimension."
                         <<iter->function()->targetDim()<<" != "<<m_system.unkSize(unk_)<<"\n");

        const size_t unk = iter->unknown();
        if(unk!=unk_)
            continue;
        const size_t patchIdx = iter->patch();
        const index_t sideIdx = iter->side();

        size_t row_Edge = 0;
        for (size_t numBdy = 0; numBdy < m_pde_ptr->domain().boundaries().size(); numBdy++ )
            if (m_pde_ptr->domain().boundaries()[numBdy].patch == patchIdx && m_pde_ptr->domain().boundaries()[numBdy].m_index == sideIdx)
                row_Edge = numBdy;

        gsMultiPatch<T> multiPatch_Edges;
        gsTensorBSplineBasis<2,real_t> temp_basis = dynamic_cast<gsTensorBSplineBasis<2,real_t>  &>(m_bases[unk_].basis(patchIdx));
        for (size_t i = 0; i < numBoundaryEdgeFunctions[row_Edge+1] - numBoundaryEdgeFunctions[row_Edge]; i++)
        {
            index_t ii = numBoundaryEdgeFunctions[row_Edge] + i;
            multiPatch_Edges.addPatch(temp_basis.makeGeometry(g1System.getSingleBasis(ii, patchIdx).transpose()));
        }


        // VERTEX
        std::pair<index_t,index_t> vertex_pair;
        switch (sideIdx)
        {
            case 1:
                vertex_pair = std::make_pair(1,3);
                break;
            case 2:
                vertex_pair = std::make_pair(2,4);
                break;
            case 3:
                vertex_pair = std::make_pair(1,2);
                break;
            case 4:
                vertex_pair = std::make_pair(3,4);
                break;
            default:
                break;
        }

        size_t row_Vertex_0 = 0, row_Vertex_1 = 0;
        for (size_t numVert = 0; numVert < m_pde_ptr->domain().vertices().size(); numVert++ )
            for(size_t j = 0; j < m_pde_ptr->domain().vertices()[numVert].size(); j++)
            {
                if (m_pde_ptr->domain().vertices()[numVert][j].patch == patchIdx && m_pde_ptr->domain().vertices()[numVert][j].m_index == vertex_pair.first)
                    row_Vertex_0 = numVert;
                if (m_pde_ptr->domain().vertices()[numVert][j].patch == patchIdx && m_pde_ptr->domain().vertices()[numVert][j].m_index == vertex_pair.second)
                    row_Vertex_1 = numVert;
            }



        gsMultiPatch<T> multiPatch_Vertex_0, multiPatch_Vertex_1;
        for (size_t i = 0; i < numBoundaryVertexFunctions[row_Vertex_0+1] - numBoundaryVertexFunctions[row_Vertex_0]; i++)
        {
            index_t ii =  numBoundaryVertexFunctions[row_Vertex_0] + i;
            multiPatch_Vertex_0.addPatch(temp_basis.makeGeometry(g1System.getSingleBasis(ii, patchIdx).transpose()));
        }


        for (size_t i = 0; i < numBoundaryVertexFunctions[row_Vertex_1+1] - numBoundaryVertexFunctions[row_Vertex_1]; i++)
        {
            index_t ii =  numBoundaryVertexFunctions[row_Vertex_1] + i;
            multiPatch_Vertex_1.addPatch(temp_basis.makeGeometry(g1System.getSingleBasis(ii, patchIdx).transpose()));
        }

//        if (patchIdx == 0 && sideIdx == 3)
//            gsWriteParaview(multiPatch_Vertex_0.patch(0),"test",5000);

        const gsBasis<T> & basis = m_bases[unk_].basis(patchIdx); // Assume that the basis is the same for all the basis functions

        const gsGeometry<T> & patch = m_pde_ptr->patches()[patchIdx];

        // Set up quadrature to degree+1 Gauss points per direction,
        // all lying on iter->side() except from the direction which
        // is NOT along the element

        gsGaussRule<T> bdQuRule(basis, 1.0, 1, iter->side().direction());

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(iter->side());

        for(; bdryIter->good(); bdryIter->next() )
        {
            bdQuRule.mapTo( bdryIter->lowerCorner(), bdryIter->upperCorner(),
                            md.points, quWeights);

            //geoEval->evaluateAt( md.points );
            patch.computeMap(md);

            // the values of the boundary condition are stored
            // to rhsVals. Here, "rhs" refers to the right-hand-side
            // of the L2-projection, not of the PDE.
            rhsVals = iter->function()->eval( m_pde_ptr->domain()[patchIdx].eval( md.points ) );

            basisVals.setZero(numBoundaryEdgeFunctions[row_Edge+1] - numBoundaryEdgeFunctions[row_Edge] +
                numBoundaryVertexFunctions[row_Vertex_0+1] - numBoundaryVertexFunctions[row_Vertex_0] +
                numBoundaryVertexFunctions[row_Vertex_1+1] - numBoundaryVertexFunctions[row_Vertex_1], md.points.dim().second);
            for (size_t i = 0; i < numBoundaryEdgeFunctions[row_Edge+1] - numBoundaryEdgeFunctions[row_Edge]; i++) // Edge
                basisVals.row(i) += multiPatch_Edges.patch(i).eval(md.points);

            for (size_t i = 0; i < numBoundaryVertexFunctions[row_Vertex_0+1] - numBoundaryVertexFunctions[row_Vertex_0]; i++) // Left vertex
                basisVals.row(numBoundaryEdgeFunctions[row_Edge+1] - numBoundaryEdgeFunctions[row_Edge] + i) += multiPatch_Vertex_0.patch(i).eval(md.points);

            for (size_t i = 0; i < numBoundaryVertexFunctions[row_Vertex_1+1] - numBoundaryVertexFunctions[row_Vertex_1]; i++) // Right vertex
                basisVals.row(numBoundaryEdgeFunctions[row_Edge+1] - numBoundaryEdgeFunctions[row_Edge] + numBoundaryVertexFunctions[row_Vertex_0+1] -
                      numBoundaryVertexFunctions[row_Vertex_0] + i) += multiPatch_Vertex_1.patch(i).eval(md.points);


            globIdxAct.setZero(multiPatch_Edges.nPatches() + multiPatch_Vertex_0.nPatches() + multiPatch_Vertex_1.nPatches(),1);
            gsVector<unsigned> vec;
            if (numBoundaryEdgeFunctions[row_Edge+1] - numBoundaryEdgeFunctions[row_Edge] == 2)
            {
                vec.resize(2);
                vec.at(0) = numBoundaryEdgeFunctions[row_Edge] - numBoundaryEdgeFunctions[0];
                vec.at(1) = numBoundaryEdgeFunctions[row_Edge] +1 - numBoundaryEdgeFunctions[0];
            }
            else
                vec.setLinSpaced(numBoundaryEdgeFunctions[row_Edge+1] - numBoundaryEdgeFunctions[row_Edge], numBoundaryEdgeFunctions[row_Edge] - numBoundaryEdgeFunctions[0],
                                 numBoundaryEdgeFunctions[row_Edge+1] - numBoundaryEdgeFunctions[0]);

            globIdxAct.block(0,0,multiPatch_Edges.nPatches(),1) = vec;

            if (numBoundaryVertexFunctions[row_Vertex_0+1] - numBoundaryVertexFunctions[row_Vertex_0] == 2)
            {
                vec.resize(2);
                vec.at(0) = numBoundaryVertexFunctions[row_Vertex_0] - numBoundaryEdgeFunctions[0];
                vec.at(1) = numBoundaryVertexFunctions[row_Vertex_0] - numBoundaryEdgeFunctions[0] + 1;
            }
            else
                vec.setLinSpaced(numBoundaryVertexFunctions[row_Vertex_0+1] - numBoundaryVertexFunctions[row_Vertex_0],
                                 numBoundaryVertexFunctions[row_Vertex_0] - numBoundaryEdgeFunctions[0],
                                 numBoundaryVertexFunctions[row_Vertex_0+1] - numBoundaryEdgeFunctions[0]);

            globIdxAct.block(multiPatch_Edges.nPatches(),0,multiPatch_Vertex_0.nPatches(),1) = vec;

            if (numBoundaryVertexFunctions[row_Vertex_1+1] - numBoundaryVertexFunctions[row_Vertex_1] == 2)
            {
                vec.resize(2);
                vec.at(0) = numBoundaryVertexFunctions[row_Vertex_1] - numBoundaryEdgeFunctions[0];
                vec.at(1) = numBoundaryVertexFunctions[row_Vertex_1] - numBoundaryEdgeFunctions[0] + 1;
            }
            else
                vec.setLinSpaced(numBoundaryVertexFunctions[row_Vertex_1+1] - numBoundaryVertexFunctions[row_Vertex_1],
                                 numBoundaryVertexFunctions[row_Vertex_1] - numBoundaryEdgeFunctions[0],
                                 numBoundaryVertexFunctions[row_Vertex_1+1] - numBoundaryEdgeFunctions[0]);

            globIdxAct.block(multiPatch_Edges.nPatches() + multiPatch_Vertex_0.nPatches(),0,multiPatch_Vertex_1.nPatches(),1) = vec;


            std::vector<index_t> eltBdryFcts;
            eltBdryFcts.reserve(g1System.boundary_size());
            for( size_t i=0; i < multiPatch_Edges.nPatches(); i++)
                eltBdryFcts.push_back( i );
            for( size_t i=0; i < multiPatch_Vertex_0.nPatches(); i++)
                eltBdryFcts.push_back( multiPatch_Edges.nPatches() + i );
            for( size_t i=0; i < multiPatch_Vertex_1.nPatches(); i++)
                eltBdryFcts.push_back( multiPatch_Edges.nPatches() + multiPatch_Vertex_0.nPatches() + i );


            // Do the actual assembly:
            for( index_t k=0; k < md.points.cols(); k++ )
            {
                const T weight_k = quWeights[k] * md.measure(k);

                // Only run through the active boundary functions on the element:
                for( size_t i0=0; i0 < eltBdryFcts.size(); i0++ )
                {
                    unsigned ii, jj;
                    // Each active boundary function/DOF in eltBdryFcts has...
                    // ...the above-mentioned "element-wise index"
                    const index_t i = eltBdryFcts[i0];
                    // ...the boundary index.
                    ii = globIdxAct.at(i);

                    for( size_t j0=0; j0 < eltBdryFcts.size(); j0++ )
                    {
                        const index_t j = eltBdryFcts[j0];
                        jj = globIdxAct.at(j);

                        // Use the "element-wise index" to get the needed
                        // function value.
                        // Use the boundary index to put the value in the proper
                        // place in the global projection matrix.
                        projMatEntries.add(ii, jj, weight_k * basisVals(i,k) * basisVals(j,k));
                    } // for j
                    globProjRhs.row(ii) += weight_k *  basisVals(i,k) * rhsVals.col(k).transpose();
                } // for i
            } // for k
        } // bdryIter

    } // boundaryConditions-Iterator
    gsSparseMatrix<T> globProjMat( g1System.boundary_size(),
                                   g1System.boundary_size());
    globProjMat.setFrom( projMatEntries );
    globProjMat.makeCompressed();

    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    typename gsSparseSolver<T>::CGDiagonal solver;
    m_g1_ddof = solver.compute( globProjMat ).solve ( globProjRhs );

}


} // namespace gismo
