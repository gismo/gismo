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
//#include <gsAssembler/gsVisitorBiharmonic.h>
#include <gsAssembler/gsG1ASVisitorBiharmonic.h>
#include <gsAssembler/gsVisitorNeumann.h>
#include <gsAssembler/gsVisitorLaplaceBoundaryBiharmonic.h>
//#include <gsAssembler/gsVisitorNitscheBiharmonic.h>
#include <gsAssembler/gsVisitorMixed.h>

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
template <class T, class bhVisitor = gsG1ASVisitorBiharmonic<T> >
//template <class T, class bhVisitor = gsVisitorBiharmonic<T> >
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

    void assemble(bool isogeometric = true, index_t patchIdx = 0);

    void applyMixed(gsVisitorMixed<T> & visitor, index_t patchIdx = 0);

    void constructSolution(const gsMatrix<T>& solVector,
                           gsMultiPatch<T>& result, int unk = 0);

    void computeDirichletDofsL2Proj(std::vector<gsMultiBasis<>> & mb, gsG1System<real_t> &  g1System, bool isogeometric);
    void computeDirichletAndNeumannDofsL2Proj(gsG1System<real_t> &  g1System);


    void plotParaview(gsField<> &solField_interior, std::vector<gsMultiPatch<T>> &result);

    gsMatrix<> get_bValue() { return m_g1_ddof; };

    void constructSystem(const gsSparseMatrix<> & mat22, const gsMatrix<> &f2, const gsSparseMatrix<> & mat12, const gsSparseMatrix<> & mat21, const gsMultiBasis<> & multiBasis);

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
void gsG1BiharmonicAssembler<T,bhVisitor>::constructSystem(const gsSparseMatrix<> &mat22, const gsMatrix<> &f2, const gsSparseMatrix<> &mat12, const gsSparseMatrix<> &mat21, const gsMultiBasis<> & multiBasis)
{


    //gsInfo << "DIMENSION: " << m_system.matrix().dim() << "\n" << mat22.dim() << "\n" << mat12.dim() << "\n" << mat21.dim() << "\n";


    gsStopwatch clock;

    index_t rows_K = m_system.matrix().rows() + mat22.rows();
    index_t cols_K = m_system.matrix().cols() + mat22.cols();

    index_t rows_f = m_system.rhs().rows() + f2.rows();
    index_t cols_f = m_system.rhs().cols();

    int nnz = (int) 1.2 * (m_system.matrix().nonZeros() + mat22.nonZeros() + mat12.nonZeros() + mat21.nonZeros());

    //K.conservativeResize(rows_K*2,cols_K*2);
    gsSparseMatrix<real_t> K_temp(rows_K,cols_K);

    clock.restart();
    index_t rows_f_old = m_system.rhs().rows();
    m_system.rhs().conservativeResize(rows_f,cols_f);
    m_system.rhs().block(rows_f_old,0,f2.rows(),cols_f) = f2;
    //gsInfo << "Matrix manipulation: " << clock.stop() << "\n";

    clock.restart();
    K_temp.reserve(nnz);

    typedef Eigen::Triplet<real_t> TT;
    std::vector<TT> tripletList;
    tripletList.reserve(nnz);
    for (int k = 0; k < m_system.matrix().outerSize(); ++k)
        for (gsSparseMatrix<real_t>::InnerIterator it(m_system.matrix(), k); it; ++it)
        {
            tripletList.push_back(TT(it.row(), it.col(), it.value()));
            //gsInfo << it.value() << " : " << mat12.at(it.row(),it.col()) << "\n";
        }


    for (int k = 0; k < mat22.outerSize(); ++k)
        for (gsSparseMatrix<real_t>::InnerIterator it(mat22, k); it; ++it)
            tripletList.push_back(TT(it.row()+  m_system.matrix().rows(), it.col()+  m_system.matrix().cols(), it.value()));
/*
    for (int k = 0; k < m_system.matrix().outerSize(); ++k)
        for (gsSparseMatrix<real_t>::InnerIterator it(m_system.matrix(), k); it; ++it)
            tripletList.push_back(TT(it.row()+  m_system.matrix().cols(), it.col() , it.value()));

    for (int k = 0; k < m_system.matrix().outerSize(); ++k)
        for (gsSparseMatrix<real_t>::InnerIterator it(m_system.matrix(), k); it; ++it)
            tripletList.push_back(TT(it.row()+ m_system.matrix().rows(), it.col(), it.value()));
*/

    for (int k = 0; k < mat12.outerSize(); ++k)
        for (gsSparseMatrix<real_t>::InnerIterator it(mat12, k); it; ++it)
        {
            if (it.row() < multiBasis.basis(0).size())
            {
                //gsInfo << it.row() << " : " << it.col() << " : " << it.value() << " : " << m_system.matrix().at(it.row(), it.col()-m_system.matrix().cols()/2) << "\n";
                tripletList.push_back(TT(it.row(), it.col()+multiBasis.basis(0).size(), it.value()));
            }

            if (it.row() >= multiBasis.basis(0).size())
            {
                tripletList.push_back(TT(it.row()+multiBasis.basis(0).size(), it.col(), it.value()));
            }
        }


    for (int k = 0; k < mat21.outerSize(); ++k)
        for (gsSparseMatrix<real_t>::InnerIterator it(mat21, k); it; ++it)
        {
            if (it.row() < multiBasis.basis(0).size())
            {
                //gsInfo << it.row() << " : " << it.col() << " : " << it.value() << " : " << m_system.matrix().at(it.row()+m_system.matrix().cols()/2, it.col()) << "\n";
                tripletList.push_back(TT(it.row()+multiBasis.basis(0).size(), it.col()+multiBasis.basis(1).size()+multiBasis.basis(0).size(), it.value()));
            }

            if (it.row() >= multiBasis.basis(0).size())
            {
                tripletList.push_back(TT(it.row()+multiBasis.basis(1).size()+multiBasis.basis(0).size(), it.col()+multiBasis.basis(0).size(), it.value()));
            }

        }



    K_temp.setFromTriplets(tripletList.begin(), tripletList.end());
    m_system.matrix().swap(K_temp);
    m_system.matrix().makeCompressed();
    //gsInfo << "Sparse manipulation: " << clock.stop() << " with " << K_temp.nonZeros() << "\n";
}


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
    index_t npts = 10000;
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
    //Base::scalarProblemGalerkinRefresh();

    gsDofMapper map(m_bases[0]);
    map.finalize();

    // 2. Create the sparse system
    m_system = gsSparseSystem<T>(map);//1,1

}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::assemble(bool isogeometric, index_t patchIdx)
{
    GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][1], 2, 1, 0.333333);

    m_system.reserve(nz, this->pde().numRhs());

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    m_ddof.resize(m_system.numUnknowns());
    m_ddof[0].setZero(m_system.colMapper(0).boundarySize(), m_system.unkSize(0) * m_system.rhs().cols());

    if (isogeometric)
    {
        // Assemble volume integrals
        Base::template push<bhVisitor>();
    }

    // Neuman conditions of first kind
    Base::template push<gsVisitorNeumann<T> >(
        m_ppde.bcFirstKind().neumannSides() );

    // Laplace conditions of second kind
    Base::template push<gsVisitorLaplaceBoundaryBiharmonic<T> >(
        m_ppde.bcSecondKind().laplaceSides() );


    if (m_options.getInt("InterfaceStrategy") == iFace::dg)
        gsWarn << "DG option ignored.\n";

    if (!isogeometric)
    {
        gsVisitorMixed<T> visitorMixed;
        applyMixed(visitorMixed, patchIdx);

    }

    /*
    // If requested, force Dirichlet boundary conditions by Nitsche's method
    this->template push<gsVisitorNitscheBiharmonic<T> >(
    m_ppde.bcSecondKind().dirichletSides() );
    */

    // Assembly is done, compress the matrix
    Base::finalize();
}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::applyMixed(gsVisitorMixed<T> & visitor, index_t patchIdx)
{
    const int patchIndex1      = patchIdx;
    const int patchIndex2      = patchIdx == 0 ? 1 : 0;

    const gsBasis<T> & B1 = m_bases[0][patchIndex1];// (!) unknown 0
    const gsBasis<T> & B2 = m_bases[0][patchIndex2];

#pragma omp parallel
{
    gsQuadRule<T> quRule ; // Quadrature rule
    gsMatrix<T> quNodes1;// Mapped nodes
    gsVector<T> quWeights;         // Mapped weights

    gsVisitorMixed<T>
#ifdef _OPENMP
        // Create thread-private visitor
    visitor_(visitor);
    const int tid = omp_get_thread_num();
    const int nt  = omp_get_num_threads();
#else
        &visitor_ = visitor;
#endif
    // Initialize
    visitor_.initialize(B1, quRule);

    const gsGeometry<T> & patch1 = m_pde_ptr->patches()[patchIndex1];

    // Initialize domain element iterators
    typename gsBasis<T>::domainIter domIt = B1.makeDomainIterator();

    // Start iteration over elements
#ifdef _OPENMP
        for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
        for (; domIt->good(); domIt->next() )
#endif
        {
            // Compute the quadrature rule on both sides
            quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes1, quWeights);


            // Perform required evaluations on the quadrature nodes
            visitor_.evaluate(B1, patch1, B2, quNodes1);

            // Assemble on element
            visitor_.assemble(*domIt,*domIt, quWeights);

            // Push to global patch matrix (m_rhs is filled in place)
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(patchIndex1, patchIndex2, m_ddof, m_system);
         }
}//omp parallel
}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::computeDirichletDofsL2Proj(std::vector<gsMultiBasis<>> & mb, gsG1System<real_t> &  g1System, bool isogeometric)
{
    gsVector<> numBoundaryVertexFunctions = g1System.get_numBoundaryVertexFunctions();
    gsVector<> numBoundaryEdgeFunctions = g1System.get_numBoundaryEdgeFunctions();

    size_t unk_ = 0;

    m_g1_ddof.resize( g1System.boundary_size(), m_system.unkSize(unk_)*m_system.rhs().cols());  //m_pde_ptr->numRhs() );
    m_g1_ddof.setZero();

/*
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
        gsTensorBSplineBasis<2,real_t> temp_basis = dynamic_cast<gsTensorBSplineBasis<2,real_t>  &>(mb[0].basis(patchIdx)); // TODO general
        gsTensorBSplineBasis<2,real_t> temp_basis2 = dynamic_cast<gsTensorBSplineBasis<2,real_t>  &>(mb[isogeometric ? 0 : 1].basis(patchIdx)); // TODO general

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
            if (row_Vertex_0 == 1 || row_Vertex_0 == 3)
                multiPatch_Vertex_0.addPatch(temp_basis2.makeGeometry(g1System.getSingleInterfaceBasis(ii, patchIdx).transpose()));
            else
                multiPatch_Vertex_0.addPatch(temp_basis.makeGeometry(g1System.getSingleBasis(ii, patchIdx).transpose()));
        }


        for (size_t i = 0; i < numBoundaryVertexFunctions[row_Vertex_1+1] - numBoundaryVertexFunctions[row_Vertex_1]; i++)
        {
            index_t ii =  numBoundaryVertexFunctions[row_Vertex_1] + i;
            if (row_Vertex_1 == 1 || row_Vertex_1 == 3)
                multiPatch_Vertex_1.addPatch(temp_basis2.makeGeometry(g1System.getSingleInterfaceBasis(ii, patchIdx).transpose()));
            else
                multiPatch_Vertex_1.addPatch(temp_basis.makeGeometry(g1System.getSingleBasis(ii, patchIdx).transpose()));
        }

        if (patchIdx == 0 && row_Vertex_1 == 3)
            gsWriteParaview(multiPatch_Vertex_1,"test",5000);


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

                        //gsInfo << "HIER: " << jj << " : " << ii << " : " << i <<" : " << j << " : " << k << "\n";
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

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::computeDirichletAndNeumannDofsL2Proj(gsG1System<real_t> &  g1System)
{
    gsDofMapper mapper(m_bases[0]);

    gsMatrix<unsigned> act;
    for (typename gsBoundaryConditions<T>::bcContainer::const_iterator it
        = m_ppde.bcFirstKind().dirichletSides().begin(); it!= m_ppde.bcFirstKind().dirichletSides().end(); ++it)
    {
        act = m_bases[0].basis(it->patch()).boundaryOffset(it->side().index(), 0); // First
        mapper.markBoundary(it->patch(), act);
    }

    for (typename gsBoundaryConditions<T>::bcContainer::const_iterator it
        = m_ppde.bcSecondKind().neumannSides().begin(); it!= m_ppde.bcSecondKind().neumannSides().end(); ++it)
    {
        act = m_bases[0].basis(it->patch()).boundaryOffset(it->side().index(), 1); // Second
        mapper.markBoundary(it->patch(), act);

    }

    // Corner boundary dofs
    for(size_t numVer=0; numVer < m_ppde.domain().vertices().size(); numVer++)
    {
        std::vector<patchCorner> allcornerLists = m_ppde.domain().vertices()[numVer];
        std::vector<size_t> patchIndex;
        std::vector<size_t> vertIndex;
        for(size_t j = 0; j < allcornerLists.size(); j++)
        {
            patchIndex.push_back(allcornerLists[j].patch);
            vertIndex.push_back(allcornerLists[j].m_index);
        }

        if(g1System.get_kindOfVertex()[numVer] == 0) // Internal vertex
            continue;

        for(size_t numPat=0; numPat < patchIndex.size(); numPat++)
        {
            size_t pID = patchIndex[numPat];
            size_t vID = vertIndex[numPat];
            index_t dim_u = m_bases[0].basis(pID).component(0).size();
            index_t dim_v = m_bases[0].basis(pID).component(1).size();

            index_t supp_size = 3; // TODO FIX SUPPORT SIZE

            index_t start_u = 0, start_v = 0, end_u = 0, end_v = 0;
            switch (vID)
            {
                case 1:
                    start_u = 0;
                    start_v = 0;
                    end_u = supp_size;
                    end_v = supp_size;
                    break;
                case 2:
                    start_u = dim_u - supp_size;
                    start_v = 0;
                    end_u = dim_u;
                    end_v = supp_size;
                    break;
                case 3:
                    start_u = 0;
                    start_v = dim_v - supp_size;
                    end_u = supp_size;
                    end_v = dim_v;
                    break;
                case 4:
                    start_u = dim_u - supp_size;
                    start_v = dim_v - supp_size;
                    end_u = dim_u;
                    end_v = dim_v;
                    break;
                default:
                    break;
            }

            gsMatrix<unsigned> bfID((end_u - start_u) * (end_v - start_v),1);
            bfID.setZero();

            index_t ii = 0;
            for(index_t j = start_v; j < end_v; j++)
                for(index_t i = start_u; i < end_u; i++)
                {
                    bfID(ii,0) = j*dim_u + i;
                    ii += 1;
                }

            mapper.markBoundary(pID,bfID);
        }


    }

    mapper.finalize();

    m_g1_ddof.resize( g1System.boundary_size(), m_system.unkSize(0)*m_system.rhs().cols());  //m_pde_ptr->numRhs() );

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>        globProjRhs;
    globProjRhs.setZero( mapper.boundarySize(), m_system.unkSize(0)*m_system.rhs().cols() );

    // Temporaries
    gsVector<T> quWeights;

    gsMatrix<T> rhsVals, rhsVals2;
    gsMatrix<unsigned> globIdxAct;
    gsMatrix<T> basisVals, basisGrads, physBasisGrad;

    gsVector<T> unormal;

    real_t lambda = 0.01;

    gsMapData<T> md(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM);


    typename gsBoundaryConditions<T>::const_iterator
        iter_dir = m_pde_ptr->bc().dirichletBegin();

    for ( typename gsBoundaryConditions<T>::const_iterator
              iter = m_ppde.bcSecondKind().neumannBegin();
          iter != m_ppde.bcSecondKind().neumannEnd(); ++iter )
    {

        GISMO_ASSERT(iter->function()->targetDim() == m_system.unkSize(0)*m_system.rhs().cols(),
                     "Given Dirichlet boundary function does not match problem dimension."
                         <<iter->function()->targetDim()<<" != "<<m_system.unkSize(0)<<"\n");

        const int unk = iter->unknown();

        const int patchIdx   = iter->patch();
        const gsBasis<T> & basis = (m_bases[unk])[patchIdx];

        GISMO_ASSERT(iter_dir->patch() == patchIdx && iter_dir->side().index() == iter->side().index(),
                     "Given Dirichlet boundary edge does not match the neumann edge."
                         <<iter_dir->patch()<<" != "<<patchIdx<<" and "
                         <<iter_dir->side().index()<<" != "<<iter->side().index()<<"\n");

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
            rhsVals = iter_dir->function()->eval( m_pde_ptr->domain()[patchIdx].eval( md.points ) );
            rhsVals2 = iter->function()->eval( m_pde_ptr->domain()[patchIdx].eval( md.points ) );

            basis.eval_into( md.points, basisVals);
            basis.deriv_into( md.points, basisGrads);

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
            basis.active_into(md.points.col(0), globIdxAct );
            mapper.localToGlobal( globIdxAct, patchIdx, globIdxAct);

            // Out of the active functions/DOFs on this element, collect all those
            // which correspond to a boundary DOF.
            // This is checked by calling mapper.is_boundary_index( global Index )

            // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
            // something like a "element-wise index"
            std::vector<index_t> eltBdryFcts;
            eltBdryFcts.reserve(mapper.boundarySize());
            for( index_t i=0; i < globIdxAct.rows(); i++)
                if( mapper.is_boundary_index(globIdxAct(i,0)) )
                    eltBdryFcts.push_back( i );

            // Do the actual assembly:
            for( index_t k=0; k < md.points.cols(); k++ )
            {
                // Compute the outer normal vector on the side
                outerNormal(md, k, iter->side(), unormal);

                // Multiply quadrature weight by the measure of normal
                //const T weight_k = quWeights[k] * md.measure(k);
                const T weight_k = quWeights[k] * unormal.norm();

                unormal.normalize();

                transformGradients(md, k, basisGrads, physBasisGrad);

                // Only run through the active boundary functions on the element:
                for( size_t i0=0; i0 < eltBdryFcts.size(); i0++ )
                {
                    // Each active boundary function/DOF in eltBdryFcts has...
                    // ...the above-mentioned "element-wise index"
                    const unsigned i = eltBdryFcts[i0];
                    // ...the boundary index.
                    const unsigned ii = mapper.global_to_bindex( globIdxAct( i ));

                    for( size_t j0=0; j0 < eltBdryFcts.size(); j0++ )
                    {
                        const unsigned j = eltBdryFcts[j0];
                        const unsigned jj = mapper.global_to_bindex( globIdxAct( j ));

                        // Use the "element-wise index" to get the needed
                        // function value.
                        // Use the boundary index to put the value in the proper
                        // place in the global projection matrix.
                        projMatEntries.add(ii, jj, weight_k * (basisVals(i,k) * basisVals(j,k) + lambda *
                            ((physBasisGrad.col(i).transpose() * unormal)(0,0) * (physBasisGrad.col(j).transpose() * unormal )(0,0))));
                    } // for j

                    globProjRhs.row(ii) += weight_k * ( basisVals(i,k) * rhsVals.col(k).transpose() + lambda *
                        ( physBasisGrad.col(i).transpose() * unormal ) * (rhsVals2.col(k).transpose() * unormal));

                } // for i
            } // for k
        } // bdryIter
        iter_dir++;
    } // boundaryConditions-Iterator

    gsSparseMatrix<T> globProjMat( mapper.boundarySize(), mapper.boundarySize() );
    globProjMat.setFrom( projMatEntries );
    globProjMat.makeCompressed();

    gsSparseMatrix<T> B_0_sparse;
    B_0_sparse.resize(g1System.boundary_size(), mapper.boundarySize());
    B_0_sparse.setZero();

    for(size_t i = 0; i < g1System.boundary_size(); i++)
    {
        for (size_t patchIdx = 0; patchIdx < m_ppde.domain().nPatches(); patchIdx++)
            for(index_t j = 0; j < (m_bases[0])[patchIdx].size(); j++)
            {
               if (mapper.is_boundary(j,patchIdx))
               {
                   index_t jj = mapper.bindex(j,patchIdx);
                   if ((g1System.getSingleBoundaryBasis(i,patchIdx))(0,j) * (g1System.getSingleBoundaryBasis(i,patchIdx))(0,j) > 10e-25)
                       B_0_sparse.insert(i,jj) = (g1System.getSingleBoundaryBasis(i,patchIdx))(0,j);
               }
            }
    }
    B_0_sparse.makeCompressed();

    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    typename gsSparseSolver<T>::CGDiagonal solver;
    m_g1_ddof = solver.compute( B_0_sparse * globProjMat * B_0_sparse.transpose() ).solve ( B_0_sparse * globProjRhs );
    */
}

} // namespace gismo
