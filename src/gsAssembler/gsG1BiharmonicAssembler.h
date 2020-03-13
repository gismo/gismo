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

    void constructG1Solution(const gsMatrix<T> &solVector, gsField<> &solField_interior, std::vector<gsMultiPatch<T>> &result, gsG1System<real_t> & g1System);



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
    for (index_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p)
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
void gsG1BiharmonicAssembler<T,bhVisitor>::constructG1Solution(const gsMatrix<T> &solVector, gsField<> &solField_interior, std::vector<gsMultiPatch<T>> &result, gsG1System<real_t> & g1System)
{

    gsMultiPatch<T> init_edge;
    std::vector<gsMultiPatch<T>> g1Basis(m_pde_ptr->domain().nPatches(), init_edge);
    for ( index_t rowEdge = 0; rowEdge < m_pde_ptr->domain().boundaries().size(); rowEdge++ )
    {
        index_t patchIdx = m_pde_ptr->domain().boundaries()[rowEdge].patch;
        gsTensorBSplineBasis<2,real_t> temp_basis = dynamic_cast<gsTensorBSplineBasis<2,real_t>  &>(m_bases[0].basis(patchIdx));
        for (index_t i = 0; i < g1System.numBoundaryEdgeFcts(rowEdge); i++) // each boundary edge
            g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(g1System.getSingleBasis(g1System.interfaceSize() + g1System.getAllBoundaryEdgeFunctions()[rowEdge] + i, patchIdx).transpose() *
                m_g1_ddof.at(g1System.getAllBoundaryEdgeFunctions()[rowEdge] + i)));
        for (index_t i = 0; i < g1System.numEdgeFcts(rowEdge); i++) // each edge
            g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(g1System.getSingleBasis(g1System.interfaceSize() + g1System.boundarySize_Edge() + g1System.getAllEdgeFunctions()[rowEdge] + i, patchIdx).transpose() *
                solVector.at(g1System.interfaceSize() + g1System.boundarySize_Edge() + g1System.getAllEdgeFunctions()[rowEdge] + i)));
    }

    for ( index_t rowInt = 0; rowInt < m_pde_ptr->domain().interfaces().size(); rowInt++ ) // each interface edge
    {
        index_t patchIdx = m_pde_ptr->domain().interfaces()[rowInt].first().patch;
        gsTensorBSplineBasis<2,real_t> temp_basis = dynamic_cast<gsTensorBSplineBasis<2,real_t>  &>(m_bases[0].basis(patchIdx));
        for (index_t i = 0; i < g1System.numInterfaceFcts(rowInt); i++)
            g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(g1System.getSingleBasis(g1System.getAllInterfaceFunctions()[rowInt] + i, patchIdx).transpose() *
                solVector.at(g1System.getAllInterfaceFunctions()[rowInt] + i)));

        patchIdx = m_pde_ptr->domain().interfaces()[rowInt].second().patch;
        temp_basis = dynamic_cast<gsTensorBSplineBasis<2,real_t>  &>(m_bases[0].basis(patchIdx));
        for (index_t i = 0; i < g1System.numInterfaceFcts(rowInt); i++)
            g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(g1System.getSingleBasis(g1System.getAllInterfaceFunctions()[rowInt] + i, patchIdx).transpose() *
                solVector.at(g1System.getAllInterfaceFunctions()[rowInt] + i)));
    }



    for ( index_t rowVertex = 0; rowVertex < m_pde_ptr->domain().vertices().size(); rowVertex++ )
    {
        std::vector<patchCorner> allcornerLists = m_pde_ptr->domain().vertices()[rowVertex];

        for(size_t j = 0; j < allcornerLists.size(); j++)
        {
            index_t patchIdx = m_pde_ptr->domain().vertices()[rowVertex][j].patch;
            gsTensorBSplineBasis<2, real_t>
                temp_basis = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(m_bases[0].basis(patchIdx));
            for (index_t i = 0; i < g1System.numBoundaryVertexFcts(rowVertex); i++) // each boundary vertex
                g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(g1System.getSingleBasis(
                    g1System.interfaceSize() + g1System.edge_size() + g1System.boundarySize_Edge() +
                        g1System.getAllBoundaryVertexFunctions()[rowVertex] + i, patchIdx).transpose() * m_g1_ddof
                    .at(g1System.boundarySize_Edge() + g1System.getAllBoundaryVertexFunctions()[rowVertex] + i)));
            for (index_t i = 0; i < g1System.numVertexFcts(rowVertex); i++) // each dofs vertex
                g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(g1System.getSingleBasis(
                    g1System.interfaceSize() + g1System.edge_size() + g1System.boundarySize_Edge() +
                        g1System.boundarySize_Vertex() + g1System.getAllVertexFunctions()[rowVertex] + i, patchIdx).transpose() * solVector
                    .at(g1System.interfaceSize() + g1System.edge_size() + g1System.boundarySize_Edge() + g1System.boundarySize_Vertex() + g1System.getAllVertexFunctions()[rowVertex] + i)));
        }
    }


    result = g1Basis;


    std::string fn = "G1Biharmonic";
    index_t npts = 5000;
    gsParaviewCollection collection2(fn);
    std::string fileName2;

    for ( unsigned pp =0; pp < m_pde_ptr->domain().nPatches(); ++pp ) // Patches
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



        //gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
        //gsField<> exact( field.patch(i), solVal, false );
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
    Base::template push<gsVisitorNeumannBiharmonic<T> >(
        m_ppde.bcSecondKind().neumannSides());

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
    size_t unk_ = 0;

    m_g1_ddof.resize( g1System.boundarySize_Edge() + g1System.boundarySize_Vertex(), m_system.unkSize(unk_)*m_system.rhs().cols());  //m_pde_ptr->numRhs() );

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>        globProjRhs;
    globProjRhs.setZero( g1System.boundarySize_Edge() + g1System.boundarySize_Vertex() , m_system.unkSize(unk_)*m_system.rhs().cols() );

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
        const index_t patchIdx = iter->patch();
        const index_t sideIdx = iter->side();

        size_t row_Edge = 0;
        for (size_t numBdy = 0; numBdy < m_pde_ptr->domain().boundaries().size(); numBdy++ )
            if (m_pde_ptr->domain().boundaries()[numBdy].patch == patchIdx && m_pde_ptr->domain().boundaries()[numBdy].m_index == sideIdx)
                row_Edge = numBdy;

        gsMultiPatch<T> multiPatch_Edges;
        gsTensorBSplineBasis<2,real_t> temp_basis = dynamic_cast<gsTensorBSplineBasis<2,real_t>  &>(m_bases[unk_].basis(patchIdx));
        for (index_t i = 0; i < g1System.numBoundaryEdgeFcts(row_Edge); i++)
            multiPatch_Edges.addPatch(temp_basis.makeGeometry(g1System.getSingleBasis(g1System.interfaceSize() +
                       g1System.getAllBoundaryEdgeFunctions()[row_Edge] + i, patchIdx).transpose()));



        if (patchIdx == 0 && sideIdx == 2)
        {
            gsField<> fiel(m_pde_ptr->domain().patch(0),multiPatch_Edges.patch(0));
            gsWriteParaview(fiel,"test",5000);
            gsField<> field(m_pde_ptr->domain().patch(0),multiPatch_Edges.patch(1));
            gsWriteParaview(field,"test1",5000);
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
        for (index_t i = 0; i < g1System.numBoundaryVertexFcts(row_Vertex_0); i++)
            multiPatch_Vertex_0.addPatch(temp_basis.makeGeometry(g1System.getSingleBasis(g1System.interfaceSize() +
                          g1System.boundarySize_Edge() + g1System.edge_size() + g1System.getAllBoundaryVertexFunctions()[row_Vertex_0] + i, patchIdx).transpose()));

        for (index_t i = 0; i < g1System.numBoundaryVertexFcts(row_Vertex_1); i++)
            multiPatch_Vertex_1.addPatch(temp_basis.makeGeometry(g1System.getSingleBasis(g1System.interfaceSize() +
                g1System.boundarySize_Edge() + g1System.edge_size() + g1System.getAllBoundaryVertexFunctions()[row_Vertex_1] + i, patchIdx).transpose()));


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

            basisVals.setZero(g1System.numBoundaryEdgeFcts(row_Edge) + g1System.numBoundaryVertexFcts(row_Vertex_0) +
                g1System.numBoundaryVertexFcts(row_Vertex_1),md.points.dim().second);
            for (size_t i = 0; i < g1System.numBoundaryEdgeFcts(row_Edge); i++) // Edge
                basisVals.row(i) += multiPatch_Edges.patch(i).eval(md.points);

            for (size_t i = 0; i < g1System.numBoundaryVertexFcts(row_Vertex_0); i++) // Left vertex
                basisVals.row(g1System.numBoundaryEdgeFcts(row_Edge) + i) += multiPatch_Vertex_0.patch(i).eval(md.points);

            for (size_t i = 0; i < g1System.numBoundaryVertexFcts(row_Vertex_1); i++) // Right vertex
                basisVals.row(g1System.numBoundaryEdgeFcts(row_Edge) + g1System.numBoundaryVertexFcts(row_Vertex_0) + i) += multiPatch_Vertex_1.patch(i).eval(md.points);

            // Indices involved here:
            // --- Local index:
            // Index of the basis function/DOF on the patch.
            // Does not take into account any boundary or interface conditions.
            // --- Global Index:
            // Each DOF has a unique global index that runs over all patches.
            // This global index includes a re-ordering such that all eliminated
            // DOFs come at the end.
            // The global index also takes care of glued interface, i.e., corresponding
            // DOFs on different patches will have the same global index, if they are
            // glued together.
            // --- Boundary Index (actually, it's a "Dirichlet Boundary Index"):
            // The eliminated DOFs, which come last in the global indexing,
            // have their own numbering starting from zero.

            // Get the global indices (second line) of the local
            // active basis (first line) functions/DOFs:
            globIdxAct.setZero(multiPatch_Edges.nPatches() + multiPatch_Vertex_0.nPatches() + multiPatch_Vertex_1.nPatches(),1);
            gsVector<unsigned> vec;
            if (g1System.numBoundaryEdgeFcts(row_Edge) == 2)
            {
                vec.resize(2);
                vec.at(0) = g1System.getAllBoundaryEdgeFunctions()[row_Edge];
                vec.at(1) = g1System.getAllBoundaryEdgeFunctions()[row_Edge] +1;
            }
            else
                vec.setLinSpaced(g1System.numBoundaryEdgeFcts(row_Edge),g1System.getAllBoundaryEdgeFunctions()[row_Edge],g1System.getAllBoundaryEdgeFunctions()[row_Edge+1]);

            globIdxAct.block(0,0,multiPatch_Edges.nPatches(),1) = vec;

            if (g1System.numBoundaryVertexFcts(row_Vertex_0) == 2)
            {
                vec.resize(2);
                vec.at(0) = g1System.getAllBoundaryVertexFunctions()[row_Vertex_0] + g1System.boundarySize_Edge();
                vec.at(1) = g1System.getAllBoundaryVertexFunctions()[row_Vertex_0] +1 + g1System.boundarySize_Edge();
            }
            else
                vec.setLinSpaced(g1System.numBoundaryVertexFcts(row_Vertex_0),g1System.getAllBoundaryVertexFunctions()[row_Vertex_0] + g1System.boundarySize_Edge(),g1System.getAllBoundaryVertexFunctions()[row_Vertex_0+1] + g1System.boundarySize_Edge());
            globIdxAct.block(multiPatch_Edges.nPatches(),0,multiPatch_Vertex_0.nPatches(),1) = vec;

            if (g1System.numBoundaryVertexFcts(row_Vertex_1) == 2)
            {
                vec.resize(2);
                vec.at(0) = g1System.getAllBoundaryVertexFunctions()[row_Vertex_1] + g1System.boundarySize_Edge();
                vec.at(1) = g1System.getAllBoundaryVertexFunctions()[row_Vertex_1] +1 + g1System.boundarySize_Edge();
            }
            else
                vec.setLinSpaced(g1System.numBoundaryVertexFcts(row_Vertex_1),g1System.getAllBoundaryVertexFunctions()[row_Vertex_1] + g1System.boundarySize_Edge(),g1System.getAllBoundaryVertexFunctions()[row_Vertex_1+1] + g1System.boundarySize_Edge());
            globIdxAct.block(multiPatch_Edges.nPatches() + multiPatch_Vertex_0.nPatches(),0,multiPatch_Vertex_1.nPatches(),1) = vec;

            // Out of the active functions/DOFs on this element, collect all those
            // which correspond to a boundary DOF.
            // This is checked by calling mapper.is_boundary_index( global Index )

            // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
            // something like a "element-wise index"
            std::vector<index_t> eltBdryFcts;
            eltBdryFcts.reserve(g1System.boundarySize_Edge() + g1System.boundarySize_Vertex());
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
    gsSparseMatrix<T> globProjMat( g1System.boundarySize_Edge() + g1System.boundarySize_Vertex() ,
                                   g1System.boundarySize_Edge() + g1System.boundarySize_Vertex());
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
