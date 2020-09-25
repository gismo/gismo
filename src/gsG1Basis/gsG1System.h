/** @file gsG1System.h

    @brief Create a G1-System for a Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once


namespace gismo
{
template<class T>
class gsG1System
{
public:
/*
  The mapper/structure for the G1System is as follow:

  - std::vector<gsVector<>> numBasisFunctions:

    0) Interface basis functions
    1) Edge basis functions
    2) Vertex basis functions

    3) Boundary edge basis functions
    4) Boundary vertex basis functions

    5) Interior basis functions

  Each entry in numBasisFunctions has a vector which has the size of interfaces/edges/vertices/patches + 1,
  e.g. for interface basis functions:
  Assume we have k interfaces, where the i-th interface has n_i basis functions. So the vector is

  numBasisFunctions[0] =
   [0,
   n_1,
   n_1 + n_2,
   ...
   n_1 + n_2 + ... + n_k]

   That means that we have in total numBasisFunctions[0][k] interface basis functions and the basis functions are
   stored with a global number. Furthermore, we obtain the number of basis functions n_i in a single interface i with
   n_i = numBasisFunctions[0][i+1] - numBasisFunctions[0][i]

   !!! The number is always given in the global sense except the interior basis functions !!! E.g.
   The counting for the edge basis functions continue from the last number in the interface basis functions:
   numBasisFunctions[1] =
   [n_1 + n_2 + ... + n_k,
   ...]

   It follows that the last entry in the vector is the same as the first in the next vector, e.g.
   numBasisFunctions[0].last() == numBasisFunctions[1][0]

*/
    gsG1System(gsMultiPatch<> & mp,
               gsMultiBasis<> & mb,
               bool neumannBdy = false,
               bool twoPatch = false)
               : m_twoPatch(twoPatch), m_neumannBdy(neumannBdy)
    {
        numBasisFunctions.resize(6);

        if (m_twoPatch)
            initialize_twoPatch(mp,mb);
        else
            initialize(mp, mb);
    }

    void initialize_twoPatch(gsMultiPatch<> & mp, gsMultiBasis<> mb);
    void initialize(gsMultiPatch<> & mp, gsMultiBasis<> mb);
    void finalize(gsMultiPatch<> & mp, gsMultiBasis<> & mb, gsMatrix<> g1);
    gsMatrix<> solve(gsSparseMatrix<real_t> K, gsMatrix<> f);

    void insertInterfaceEdge(gsMultiPatch<> & mp, boundaryInterface item, index_t iID ,index_t bfID);
    void insertBoundaryEdge(gsMultiPatch<> & mp, patchSide item, index_t bID ,size_t bfID);
    void insertVertex(gsMultiPatch<> & mp, std::vector<size_t> patchIndex, index_t vID, index_t nDofs, index_t bfID);

    void constructG1Solution(const gsMatrix<T> &solVector, std::vector<gsMultiPatch<>> & result, const gsMultiPatch<> & geo);
    void constructSparseG1Solution(const gsMatrix<T> &solVector, gsSparseMatrix<T> & result);
    gsVector<> get_numBasisFunctions() { return numBasisFunctions[5]; }
    gsVector<> get_numInterfaceFunctions() { return numBasisFunctions[0]; }
    gsVector<> get_numBoundaryEdgeFunctions() {return numBasisFunctions[3]; };
    gsVector<> get_numBoundaryVertexFunctions() {return numBasisFunctions[4]; };
    gsVector<> get_numVertexFunctions() {return numBasisFunctions[2]; };

    gsVector<> get_kindOfVertex() {return kindOfVertex; };



    size_t boundary_size() { return numBasisFunctions[4].last() - numBasisFunctions[3][0]; }

    size_t sizePlusInterface(index_t i) { return  sizePlusInt[i]; };
    size_t sizePlusBoundary(index_t i) { return  sizePlusBdy[i]; };


    gsMatrix<> getSingleBasis(index_t global_row, index_t patchIdx) {
        return D_sparse.block(global_row, numBasisFunctions[5][patchIdx], 1, numBasisFunctions[5][patchIdx+1] - numBasisFunctions[5][patchIdx]);
    };

    gsMatrix<> getSingleBoundaryBasis(index_t boundary_row, index_t patchIdx) {
        return D_sparse.block(dim_G1_Dofs + boundary_row, numBasisFunctions[5][patchIdx], 1, numBasisFunctions[5][patchIdx+1] - numBasisFunctions[5][patchIdx]);
    };

protected:
    bool m_twoPatch, m_neumannBdy;

    size_t dim_K, dim_G1_Dofs, dim_G1_Bdy;

    std::vector<gsVector<>> numBasisFunctions;

    gsVector<> kindOfVertex;
    gsVector<size_t> sizePlusInt, sizePlusBdy;

    gsSparseMatrix<T> D_sparse, D_0_sparse, D_boundary_sparse;
    gsMatrix<> m_g1;

}; // class gsG1System

template<class T>
void gsG1System<T>::initialize_twoPatch(gsMultiPatch<> & mp, gsMultiBasis<> mb)
{
    // Number of the patches
    index_t numPatches = mp.nPatches();

    // Dimension of Basis plus at the edge
    sizePlusInt.setZero(mp.interfaces().size());
    sizePlusBdy.setZero(mp.boundaries().size());

    // Kind of vertex
    // -1 Boundary vertex
    // 0 Internal vertex
    // 1 Interface boundary vertey
    kindOfVertex.setZero(mp.vertices().size());

    // Get the dimension of the basis functions for each patch
    numBasisFunctions[5].setZero(numPatches+1);

    numBasisFunctions[0].setZero(mp.interfaces().size()+1);

    numBasisFunctions[1].setZero(mp.boundaries().size()+1);
    numBasisFunctions[3].setZero(mp.boundaries().size()+1);

    numBasisFunctions[2].setZero(mp.vertices().size()+1);
    numBasisFunctions[4].setZero(mp.vertices().size()+1);

    for (size_t i = 0; i < mb.nBases(); i++ )
        numBasisFunctions[5][i+1] = numBasisFunctions[5][i] + mb.basis(i).size();

    for (size_t i = 0; i < mp.interfaces().size(); i++)
    {
        // Get the dimension for the spaces at the edges
        index_t dir = mp.interfaces()[i].first().m_index < 3 ? 1 : 0;
        gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(mb.basis(mp.interfaces()[i].first().patch).component(dir)); // If the interface matches!!!
        index_t m_p = basis_edge.maxDegree();
        index_t m_r = 1; // Here fixed to 1 TODO MORE GENERAL
        index_t m_n = basis_edge.numElements();

        index_t numIntBdy;
        if (m_neumannBdy)
            numIntBdy = 8;
        else
            numIntBdy = 4;

        numBasisFunctions[0][i+1] = numBasisFunctions[0][i] + 2 * (m_p - m_r - 1) * (m_n - 1) + 2 * m_p + 1 - numIntBdy; // 1+ and 1- times 2
        sizePlusInt[i] = (m_p - m_r - 1) * (m_n - 1) + m_p + 1;
    }
    for (size_t i = 0; i < mp.boundaries().size(); i++)
    {
        // Get the dimension for the spaces at the edges
        index_t dir = mp.boundaries()[i].m_index < 3 ? 1 : 0;
        gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(mb.basis(mp.boundaries()[i].patch).component(dir)); // 0 -> u, 1 -> v

        if (m_neumannBdy)
        {
            numBasisFunctions[3][i+1] = numBasisFunctions[3][i] + 2*basis_edge.size() - 8; // boundary edge functions
            numBasisFunctions[1][i+1] = numBasisFunctions[1][i]; // edge functions
        }
        else
        {
            numBasisFunctions[3][i+1] = numBasisFunctions[3][i] + basis_edge.size() - 4; // boundary edge functions
            numBasisFunctions[1][i+1] = numBasisFunctions[1][i] + basis_edge.size() - 4; // edge functions
        }

        sizePlusBdy[i] = basis_edge.size() - 4;
    }

    for (size_t i = 0; i < mp.vertices().size(); i++)
    {
        if (mp.vertices()[i].size() == 1)
        {
            // |  o  o  o
            // |  x  x  o
            // |  x  x  o
            // |__________
            kindOfVertex[i] = -1; // Boundary vertex
            if (m_neumannBdy)
            {
                numBasisFunctions[2][i+1] = numBasisFunctions[2][i] + 0; // vertex functions // TODO
                numBasisFunctions[4][i+1] = numBasisFunctions[4][i] + 4; // TODO
            }
            else
            {
                numBasisFunctions[2][i+1] = numBasisFunctions[2][i] + 1; // vertex functions // TODO
                numBasisFunctions[4][i+1] = numBasisFunctions[4][i] + 3; // TODO
            }
        }
        else
        {
            gsMultiPatch<> temp_mp;
            for (size_t j = 0; j < mp.vertices()[i].size(); j++)
                temp_mp.addPatch(mp.patch(mp.vertices()[i][j].patch));

            temp_mp.computeTopology();
            if (mp.vertices()[i].size() == temp_mp.interfaces().size())
            {
                kindOfVertex[i] = 0; // Internal vertex
                numBasisFunctions[2][i+1] = numBasisFunctions[2][i] + 0; // vertex functions
                numBasisFunctions[4][i+1] = numBasisFunctions[4][i];
            }

            else
            {
                kindOfVertex[i] = 1; // Interface-Boundary vertex

                if (m_neumannBdy)
                {
                    numBasisFunctions[2][i+1] = numBasisFunctions[2][i] + 0; // TODO
                    numBasisFunctions[4][i+1] = numBasisFunctions[4][i] + 4; // TODO
                }
                else
                {
                    numBasisFunctions[2][i+1] = numBasisFunctions[2][i] + 0; // TODO
                    numBasisFunctions[4][i+1] = numBasisFunctions[4][i] + 2; // TODO
                }
            }
        }
    }
    numBasisFunctions[1] = numBasisFunctions[1].array() + numBasisFunctions[0].last();
    numBasisFunctions[2] = numBasisFunctions[2].array() + numBasisFunctions[1].last();
    numBasisFunctions[3] = numBasisFunctions[3].array() + numBasisFunctions[2].last();
    numBasisFunctions[4] = numBasisFunctions[4].array() + numBasisFunctions[3].last();

//    gsInfo << "Num Basis Functions " << numBasisFunctions[5] << "\n";
//    gsInfo << "Num Interface Functions " << numBasisFunctions[0] << "\n";
//    gsInfo << "Num Edges Functions " << numBasisFunctions[1] << "\n";
//    gsInfo << "Num Boundary Edges Functions " << numBasisFunctions[3] << "\n";
//    gsInfo << "Num Vertex Functions " << numBasisFunctions[2] << "\n";
//    gsInfo << "Num Boundary Vertex Functions " << numBasisFunctions[4] << "\n";
//    gsInfo << "Kind of Vertex Functions " << kindOfVertex << "\n";
//    gsInfo << "Size of plus space Bdy  " << sizePlusBdy << "\n";
//    gsInfo << "Size of plus space Int  " << sizePlusInt << "\n";

    // Setting the final matrix
    dim_K = numBasisFunctions[5].last(); // interior basis dimension
    dim_G1_Dofs = numBasisFunctions[2].last() ; // edges basis dimension
    dim_G1_Bdy = numBasisFunctions[4].last() - numBasisFunctions[3][0]; // vertex basis dimension

    // Full matrix
    D_sparse.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K, dim_K);
    D_sparse.reserve(3*dim_K);
    D_sparse.setZero();

    // Without boundary
    D_0_sparse.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K, dim_K);
    D_0_sparse.reserve(3*dim_K);
    D_0_sparse.setZero();

    // Only boundary
    D_boundary_sparse.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K , dim_K);
    D_boundary_sparse.reserve(3*dim_K);
    D_boundary_sparse.setZero();

    // Boundary values
    m_g1.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K, 1);
    m_g1.setZero();


}

template<class T>
void gsG1System<T>::initialize(gsMultiPatch<> & mp, gsMultiBasis<> mb)
{
    // Number of the patches
    index_t numPatches = mp.nPatches();

    // Dimension of Basis plus at the edge
    sizePlusInt.setZero(mp.interfaces().size());
    sizePlusBdy.setZero(mp.boundaries().size());

    // Kind of vertex
    // -1 Boundary vertex
    // 0 Internal vertex
    // 1 Interface boundary vertey
    kindOfVertex.setZero(mp.vertices().size());

    // Get the dimension of the basis functions for each patch
    numBasisFunctions[5].setZero(numPatches+1);

    numBasisFunctions[0].setZero(mp.interfaces().size()+1);

    numBasisFunctions[1].setZero(mp.boundaries().size()+1);
    numBasisFunctions[3].setZero(mp.boundaries().size()+1);

    numBasisFunctions[2].setZero(mp.vertices().size()+1);
    numBasisFunctions[4].setZero(mp.vertices().size()+1);

    for (size_t i = 0; i < mb.nBases(); i++ )
        numBasisFunctions[5][i+1] = numBasisFunctions[5][i] + mb.basis(i).size();

    for (size_t i = 0; i < mp.interfaces().size(); i++)
    {
        // Get the dimension for the spaces at the edges
        index_t dir = mp.interfaces()[i].first().m_index < 3 ? 1 : 0;
        gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(mb.basis(mp.interfaces()[i].first().patch).component(dir)); // If the interface matches!!!
        index_t m_p = basis_edge.maxDegree();
        index_t m_r = 1; // Here fixed to 1 TODO MORE GENERAL
        index_t m_n = basis_edge.numElements();

        numBasisFunctions[0][i+1] = numBasisFunctions[0][i] + 2 * (m_p - m_r - 1) * (m_n - 1) + 2 * m_p - 9;
        sizePlusInt[i] = (m_p - m_r - 1) * (m_n - 1) + m_p + 1;
    }
    for (size_t i = 0; i < mp.boundaries().size(); i++)
    {
        // Get the dimension for the spaces at the edges
        index_t dir = mp.boundaries()[i].m_index < 3 ? 1 : 0;
        gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(mb.basis(mp.boundaries()[i].patch).component(dir)); // 0 -> u, 1 -> v
        index_t m_p = basis_edge.maxDegree();
        index_t m_r = 1; // Here fixed to 1 TODO MORE GENERAL
        index_t m_n = basis_edge.numElements();

        if (m_neumannBdy)
        {
            numBasisFunctions[3][i+1] = numBasisFunctions[3][i] + 2*(m_p - m_r - 1) * (m_n - 1) + 2*m_p + 1 - 10; // Boundary edge
            numBasisFunctions[1][i+1] = numBasisFunctions[1][i]; // Edge
        }
        else
        {
            numBasisFunctions[3][i+1] = numBasisFunctions[3][i] + (m_p - m_r - 1) * (m_n - 1) + m_p + 1 - 6; // Boundary edge
            numBasisFunctions[1][i+1] = numBasisFunctions[1][i] + (m_p - m_r - 1) * (m_n - 1) + m_p - 4; // Edge
        }

        sizePlusBdy[i] = (m_p - m_r - 1) * (m_n - 1) + m_p + 1;
    }

    for (size_t i = 0; i < mp.vertices().size(); i++)
    {
        if (mp.vertices()[i].size() == 1)
        {
            kindOfVertex[i] = -1; // Boundary vertex
            numBasisFunctions[2][i+1] = numBasisFunctions[2][i] + 1; // TODO
            numBasisFunctions[4][i+1] = numBasisFunctions[4][i] + 6; // TODO
        }
        else
        {
            gsMultiPatch<> temp_mp;
            for (size_t j = 0; j < mp.vertices()[i].size(); j++)
                temp_mp.addPatch(mp.patch(mp.vertices()[i][j].patch));

            temp_mp.computeTopology();
            if (mp.vertices()[i].size() == temp_mp.interfaces().size())
            {
                kindOfVertex[i] = 0; // Internal vertex
                numBasisFunctions[2][i+1] = numBasisFunctions[2][i] + 6;
                numBasisFunctions[4][i+1] = numBasisFunctions[4][i];
            }

            else
            {
                kindOfVertex[i] = 1; // Interface-Boundary vertex
                numBasisFunctions[2][i+1] = numBasisFunctions[2][i] + 3; // TODO
                numBasisFunctions[4][i+1] = numBasisFunctions[4][i] + 6; // TODO
            }
        }
    }
    numBasisFunctions[1] = numBasisFunctions[1].array() + numBasisFunctions[0].last();
    numBasisFunctions[2] = numBasisFunctions[2].array() + numBasisFunctions[1].last();
    numBasisFunctions[3] = numBasisFunctions[3].array() + numBasisFunctions[2].last();
    numBasisFunctions[4] = numBasisFunctions[4].array() + numBasisFunctions[3].last();


//    gsInfo << "Num Basis Functions " << numBasisFunctions[5] << "\n";
//    gsInfo << "Num Interface Functions " << numBasisFunctions[0] << "\n";
//    gsInfo << "Num Edges Functions " << numBasisFunctions[1] << "\n";
//    gsInfo << "Num Boundary Edges Functions " << numBasisFunctions[3] << "\n";
//    gsInfo << "Num Vertex Functions " << numBasisFunctions[2] << "\n";
//    gsInfo << "Num Boundary Vertex Functions " << numBasisFunctions[4] << "\n";
//    gsInfo << "Kind of Vertex Functions " << kindOfVertex << "\n";
//    gsInfo << "Size of plus space Bdy  " << sizePlusBdy << "\n";
//    gsInfo << "Size of plus space Int  " << sizePlusInt << "\n";


    // Setting the final matrix
    dim_K = numBasisFunctions[5].last(); // interior basis dimension
    dim_G1_Dofs = numBasisFunctions[2].last() ; // edges basis dimension
    dim_G1_Bdy = numBasisFunctions[4].last() - numBasisFunctions[3][0]; // vertex basis dimension

    // Full matrix
    D_sparse.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K, dim_K);
    D_sparse.reserve(3*dim_K);
    D_sparse.setZero();

    // Without boundary
    D_0_sparse.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K, dim_K);
    D_0_sparse.reserve(3*dim_K);
    D_0_sparse.setZero();

    // Only boundary
    D_boundary_sparse.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K , dim_K);
    D_boundary_sparse.reserve(3*dim_K);
    D_boundary_sparse.setZero();

    // Boundary values
    m_g1.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K, 1);
    m_g1.setZero();


}

template<class T>
void gsG1System<T>::constructSparseG1Solution(const gsMatrix<T> & solVector,
                                              gsSparseMatrix<T> & result)
{
    result.clear();
    result = D_sparse.block(0,0,dim_G1_Dofs + dim_G1_Bdy + 1, dim_K); // + 1 for the interior solution

    // G1 solution
    for (size_t i = 0; i < dim_G1_Dofs; i++)
        result.row(i) *= solVector.at(i);
    for (size_t i = dim_G1_Dofs; i < dim_G1_Dofs + dim_G1_Bdy; i++)
        result.row(i) *= m_g1.at(i);
    for (size_t i = 0; i < dim_K; i++)
        result.insert(dim_G1_Dofs + dim_G1_Bdy,i) = solVector.at(dim_G1_Dofs + dim_G1_Bdy + i); // Interior solution


    result.makeCompressed();
}

template<class T>
void gsG1System<T>::constructG1Solution(const gsMatrix<T> & solVector, std::vector<gsMultiPatch<>> & result, const gsMultiPatch<> & geo)
{
    gsMultiPatch<T> init_mp;
    std::vector<gsMultiPatch<T>> g1Basis(geo.nPatches(), init_mp);

    for ( size_t rowInt = 0; rowInt < geo.interfaces().size(); rowInt++ ) // each interface edge
    {
        index_t patchIdx = geo.interfaces()[rowInt].first().patch;
        gsTensorBSplineBasis<2,real_t> temp_basis = dynamic_cast<gsTensorBSplineBasis<2,real_t>  &>(geo.basis(patchIdx));
        for (size_t i = 0; i < numBasisFunctions[0][rowInt+1] - numBasisFunctions[0][rowInt]; i++)
        {
            index_t ii = numBasisFunctions[0][rowInt] + i;
            g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(D_sparse.block(ii,numBasisFunctions[5][patchIdx],1,temp_basis.size()).transpose() *
                solVector.at(ii)));
        }

        patchIdx = geo.interfaces()[rowInt].second().patch;
        temp_basis = dynamic_cast<gsTensorBSplineBasis<2,real_t>  &>(geo.basis(patchIdx));
        for (size_t i = 0; i < numBasisFunctions[0][rowInt+1] - numBasisFunctions[0][rowInt]; i++)
        {
            index_t ii = numBasisFunctions[0][rowInt] + i;
            g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(D_sparse.block(ii,numBasisFunctions[5][patchIdx],1,temp_basis.size()).transpose() *
                solVector.at(ii)));
        }
    }


    for ( size_t rowEdge = 0; rowEdge < geo.boundaries().size(); rowEdge++ )
    {
        index_t patchIdx = geo.boundaries()[rowEdge].patch;
        gsTensorBSplineBasis<2,real_t> temp_basis = dynamic_cast<gsTensorBSplineBasis<2,real_t>  &>(geo.basis(patchIdx));
        for (size_t i = 0; i < numBasisFunctions[3][rowEdge+1] - numBasisFunctions[3][rowEdge]; i++) // each boundary edge
        {
            index_t ii = numBasisFunctions[3][rowEdge] + i;
            g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(D_sparse.block(ii,numBasisFunctions[5][patchIdx],1,temp_basis.size()).transpose() *
                m_g1.at(ii)));
        }
        for (size_t i = 0; i < numBasisFunctions[1][rowEdge+1] - numBasisFunctions[1][rowEdge]; i++) // each edge
        {
            index_t ii = numBasisFunctions[1][rowEdge] + i;
            g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(D_sparse.block(ii,numBasisFunctions[5][patchIdx],1,temp_basis.size()).transpose() *
                solVector.at(ii)));
        }
    }


    for ( size_t rowVertex = 0; rowVertex < geo.vertices().size(); rowVertex++ )
    {
        std::vector<patchCorner> allcornerLists = geo.vertices()[rowVertex];

        for(size_t j = 0; j < allcornerLists.size(); j++)
        {
            index_t patchIdx = geo.vertices()[rowVertex][j].patch;
            gsTensorBSplineBasis<2, real_t>
                temp_basis = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(geo.basis(patchIdx));
            for (size_t i = 0; i < numBasisFunctions[4][rowVertex+1] - numBasisFunctions[4][rowVertex]; i++) // each boundary vertex
            {
                index_t ii = numBasisFunctions[4][rowVertex] + i;
                g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(D_sparse.block(ii,numBasisFunctions[5][patchIdx],1,temp_basis.size()).transpose() *
                    m_g1.at(ii)));
            }
            for (size_t i = 0; i < numBasisFunctions[2][rowVertex+1] - numBasisFunctions[2][rowVertex]; i++) // each dofs vertex
            {
                index_t ii = numBasisFunctions[2][rowVertex] + i;
                g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(D_sparse.block(ii,numBasisFunctions[5][patchIdx],1,temp_basis.size()).transpose() *
                    solVector.at(ii)));
            }
        }
    }

    result = g1Basis;
}


template<class T>
void gsG1System<T>::insertInterfaceEdge(gsMultiPatch<> & mp, boundaryInterface item, index_t iID ,index_t bfID)
{
    // Insert all coefficients of the g1 Basis at the interface
    for (size_t np = 0; np < 2; ++np) // two interface patches
        for (index_t j = 0; j < mp.patch(np).coefs().size(); j++) // all the coefs
            if (m_twoPatch && !m_neumannBdy)
            {
                index_t plusInt = sizePlusInt[0];
                if(bfID == 0 || bfID == plusInt-1 || bfID == plusInt || bfID == 2*plusInt - 2) // first and last of +/-
                {
                    if (bfID == 0 || bfID == plusInt)
                        if (mp.patch(np).coefs().at(j) * mp.patch(np).coefs().at(j)  > 10e-25)
                        {
                            index_t jj, ii;
                            ii = numBasisFunctions[4][1] + (bfID == 0 ? 0 : 1); // Boundary
                            jj = numBasisFunctions[5][np == 0 ? item.first().patch : item.second().patch] + j;
                            D_sparse.insert(ii,jj) = mp.patch(np).coefs().at(j);
                        }
                    if (bfID == plusInt-1 || bfID == 2*plusInt - 2)
                        if (mp.patch(np).coefs().at(j) * mp.patch(np).coefs().at(j)  > 10e-25)
                        {
                            index_t jj, ii;
                            ii = numBasisFunctions[4][3] + (bfID == plusInt-1 ? 0 : 1); // Boundary
                            jj = numBasisFunctions[5][np == 0 ? item.first().patch : item.second().patch] + j;
                            D_sparse.insert(ii,jj) = mp.patch(np).coefs().at(j);
                        }
                }
                else
                {
                    index_t plusInt = sizePlusInt[0];
                    if (mp.patch(np).coefs().at(j) * mp.patch(np).coefs().at(j)  > 10e-25)
                    {
                        index_t jj, ii, bfID_shift;
                        if (bfID < plusInt - 1)
                            bfID_shift = 1;
                        else
                            bfID_shift = 3;

                        ii = numBasisFunctions[0][iID] + bfID - bfID_shift;
                        jj = numBasisFunctions[5][np == 0 ? item.first().patch : item.second().patch] + j;
                        D_sparse.insert(ii,jj) = mp.patch(np).coefs().at(j);
                    }
                }
            } // two patch
            else
            {
                if (mp.patch(np).coefs().at(j) * mp.patch(np).coefs().at(j)  > 10e-25)
                {
                    index_t jj, ii;
                    ii = numBasisFunctions[0][iID] + bfID;
                    jj = numBasisFunctions[5][np == 0 ? item.first().patch : item.second().patch] + j;
                    D_sparse.insert(ii,jj) = mp.patch(np).coefs().at(j);
                }
            }
}

template<class T>
void gsG1System<T>::insertBoundaryEdge(gsMultiPatch<> & mp, patchSide item, index_t bID, size_t bfID)
{
    // Insert all coefficients of the g1 Basis at the interface
    for (index_t j = 0; j < mp.patch(0).coefs().size(); j++) // all the coefs
        if (mp.patch(0).coefs().at(j) * mp.patch(0).coefs().at(j)  > 10e-25)
        {

            index_t jj, ii;
            if (m_twoPatch && !m_neumannBdy)
            {
                if (bfID < sizePlusBdy[bID])
                    ii = numBasisFunctions[3][bID] + bfID;

                else
                    ii = numBasisFunctions[1][bID] + bfID - sizePlusBdy[bID];
            }
            else if (!m_neumannBdy)
            {
                if (bfID < sizePlusBdy[bID] - 6)
                    ii = numBasisFunctions[3][bID] + bfID;

                else
                    ii = numBasisFunctions[1][bID] + bfID - sizePlusBdy[bID] + 6;
            }
            else if (m_neumannBdy)
            {
                ii = numBasisFunctions[3][bID] + bfID; // all edges belongs to the boundary
            }


            jj = numBasisFunctions[5][item.patch] + j;
            D_sparse.insert(ii,jj) = mp.patch(0).coefs().at(j);
        }
}

template<class T>
void gsG1System<T>::insertVertex(gsMultiPatch<> & mp, std::vector<size_t> patchIndex, index_t vID, index_t nDofs, index_t bfID)
{
    // Insert all coefficients of the g1 Basis at the interface
    for (size_t np = 0; np < mp.nPatches(); ++np) // for each patch which has the vertex
        for (index_t j = 0; j < mp.patch(np).coefs().size(); j++) // all the coefs
            if (mp.patch(np).coefs().at(j) * mp.patch(np).coefs().at(j)  > 10e-25)
            {
                index_t jj, ii = -1;
                if (kindOfVertex[vID] == 0) // interior vertex
                    ii = numBasisFunctions[2][vID] + bfID; // all six belongs to Dofs
                else if (kindOfVertex[vID] == -1) // boundary
                {
                    if (bfID < nDofs)
                        ii = numBasisFunctions[2][vID] + bfID; // Dofs
                    else
                        ii = numBasisFunctions[4][vID] + bfID - nDofs; // Boundary
                }
                else if (kindOfVertex[vID] == 1) // interface boundary
                {
                    if (bfID < nDofs)
                        ii = numBasisFunctions[2][vID] + bfID; // Dofs
                    else
                        ii = numBasisFunctions[4][vID] + bfID - nDofs; // Boundary
                }
                jj = numBasisFunctions[5][patchIndex[np]] + j;
                D_sparse.insert(ii,jj) = mp.patch(np).coefs().at(j);
            }
}

template<class T>
void gsG1System<T>::finalize(gsMultiPatch<> & mp, gsMultiBasis<> & mb, gsMatrix<> g1)
{

    gsSparseMatrix<> B_0_sparse, B_boundary_sparse, temp;
    B_0_sparse.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K, dim_G1_Dofs + dim_G1_Bdy + dim_K);
    B_0_sparse.reserve(dim_G1_Dofs + dim_K);
    B_0_sparse.setZero();

    B_boundary_sparse.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K, dim_G1_Dofs + dim_G1_Bdy + dim_K);
    B_boundary_sparse.reserve(dim_G1_Bdy);
    B_boundary_sparse.setZero();

    // Add for the G1 Dofs to 1
    //for(size_t i = 0; i < dim_G1_Dofs; i++) PASCAL
    for(size_t i = 0; i < dim_G1_Dofs+dim_G1_Bdy; i++)
        B_0_sparse.insert(i,i) = 1;

    // Add for the Boundary the edges to 1
    for(size_t i = 0; i < dim_G1_Bdy; i++)
    {
        index_t ii = dim_G1_Dofs + i;
        //B_boundary_sparse.insert(ii,ii) = 1; PASCAL
    }

    // Add the identity to the end of D
    // Construct the internal matrix
    for(size_t np = 0; np < mb.nBases(); np++)
    {
        index_t dim_u = mb.basis(np).component(0).size();
        index_t dim_v = mb.basis(np).component(1).size();
        for(index_t j = 2; j < dim_v - 2; j++)
            for(index_t i = 2; i < dim_u - 2; i++)
            {
                index_t ii = dim_G1_Dofs + dim_G1_Bdy + numBasisFunctions[5][np] + j*dim_u + i;
                index_t jj = numBasisFunctions[5][np] + j*dim_u + i;
                B_0_sparse.insert(ii,ii) = 1;
                D_sparse.insert(ii,jj) = 1;
            }
    }
    D_sparse.makeCompressed();

    // Getting D_0
    D_0_sparse = B_0_sparse * D_sparse;
    D_0_sparse.makeCompressed();

    // Compute D_Boundary
    D_boundary_sparse = B_boundary_sparse * D_sparse;
    D_boundary_sparse.makeCompressed();

    // Set up the boundary vector
    m_g1.block(dim_G1_Dofs,0,dim_G1_Bdy,1) = g1;

}

template<class T>
gsMatrix<> gsG1System<T>::solve(gsSparseMatrix<real_t> K, gsMatrix<> f)
{

    gsSparseMatrix<real_t> A = D_0_sparse * K * D_0_sparse.transpose();

    gsVector<real_t> F = D_0_sparse * f;  //- D_0_sparse * K * D_boundary_sparse.transpose() * m_g1; PASCAL
//    gsInfo << "F: " << F << "\n";
//    gsInfo << "m_g1: " << m_g1 << "\n";

    gsSparseSolver<real_t>::CGDiagonal solver;
//    gsSparseSolver<real_t>::BiCGSTABILUT solver;

    solver.compute(A);
    gsMatrix<> solVector = solver.solve(F);


//    gsInfo << "Sol: " << solVector << "\n";
    return solVector;
}

} // namespace
