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

    gsG1System(gsMultiBasis<> & interior,
               std::vector<gsG1AuxiliaryPatch> & edges,
               std::vector<gsG1AuxiliaryPatch> & vertices,
               gsG1Mapper_pascal<T> & g1Mapper):
        m_edges(edges), m_vertices(vertices), m_g1Mapper(g1Mapper)
    {
        m_map_edge = m_g1Mapper.getMapper_Edges();
        m_map_vertex = m_g1Mapper.getMapper_Edges();

        dim_K = interior.size();
        dim_I = m_map_edge.freeSize() + m_map_edge.boundarySize();
        dim_D = m_map_edge.freeSize() + m_map_edge.boundarySize() + m_map_vertex.freeSize() + m_map_vertex.boundarySize();
    }

    void setGlobalMapper(gsDofMapper global) { map_global = global; };

    void initialize();
    void assemble();

    void plotParaview(gsMultiPatch<T> & mp,
                      std::vector<gsG1AuxiliaryPatch> & interface,
                      std::vector<gsG1AuxiliaryPatch> & boundaries,
                      std::vector<gsG1AuxiliaryPatch> & vertices,
                      std::string basename);

    void plotParaview(gsMultiPatch<T> & mp,
                      std::vector<gsG1AuxiliaryPatch> & edges,
                      std::vector<gsG1AuxiliaryPatch> & vertices,
                      std::string basename);

protected:
    std::vector<gsG1AuxiliaryPatch> & m_edges, m_vertices;

    gsG1Mapper_pascal<T> m_g1Mapper;
    gsDofMapper m_map_edge;
    gsDofMapper m_map_vertex;

    index_t dim_K, dim_I, dim_D;

    gsDofMapper map_global;

    gsSparseMatrix<T> D_sparse, D_0_sparse, D_boundary_sparse;

    gsMatrix<T> m_dofs, m_g;

}; // class gsG1System

template<class T>
void gsG1System<T>::initialize()
{
    // Full matrix
    D_sparse.resize(dim_D, dim_K);
    D_sparse.reserve(3*dim_K);
    D_sparse.setZero();

    // Without boundary
    D_0_sparse.resize(dim_D, dim_K);
    D_0_sparse.reserve(3*dim_K);
    D_0_sparse.setZero();

    // Only boundary
    D_boundary_sparse.resize(dim_D , dim_K);
    D_boundary_sparse.reserve(3*dim_K);
    D_boundary_sparse.setZero();

    m_g.resize(dim_D,1);
    m_g.setZero();
}

template<class T>
void gsG1System<T>::assemble()
{


    // First insert all coefficients of the g1 Basis
    for (size_t i = 0; i < m_edges.size(); ++i) // which edge
    {
        for (size_t j = 0; j < m_edges.at(i).getG1Basis().nPatches(); j++) // which basis functions
        {
            gsMatrix<unsigned > localDof(1,1), globalDof(1,1);
            localDof << j;
            m_map_edge.localToGlobal(localDof,i,globalDof);

            for (size_t k = 0; k < m_edges.at(i).getG1Basis().basis(j).size(); k++)
            {
                if (m_edges.at(i).getG1Basis().patch(j).coefs().at(k) * m_edges.at(i).getG1Basis().patch(j).coefs().at(k) > 10e-25)
                {
                    gsMatrix<unsigned > localDof_BF(1,1), globalDof_BF;
                    localDof_BF << k;
                    map_global.localToGlobal(localDof_BF,m_edges.at(i).getGlobalPatchIndex(),globalDof_BF);
                    D_sparse.insert(globalDof.at(0),globalDof_BF.at(0)) = m_edges.at(i).getG1Basis().patch(j).coefs().at(k);
                }

            }

        }
    } // edges

    for (size_t i = 0; i < m_vertices.size(); ++i) // which vertex
    {
        for (size_t j = 0; j < m_vertices.at(i).getG1Basis().nPatches(); j++) // which basis functions
        {
            gsMatrix<unsigned > localDof(1,1), globalDof(1,1);
            localDof << j;
            m_map_vertex.localToGlobal(localDof,i,globalDof);
            gsInfo << globalDof << " : ";
            for (size_t k = 0; k < m_vertices.at(i).getG1Basis().basis(j).size(); k++)
            {
                if (m_vertices.at(i).getG1Basis().patch(j).coefs().at(k) * m_vertices.at(i).getG1Basis().patch(j).coefs().at(k) > 10e-30)
                {
                    gsMatrix<unsigned > localDof_BF(1,1), globalDof_BF;
                    localDof_BF << k;
                    map_global.localToGlobal(localDof_BF,m_vertices.at(i).getGlobalPatchIndex(),globalDof_BF);
                    //gsInfo << globalDof_BF << " \n";
                    D_sparse.insert(dim_I + globalDof.at(0),globalDof_BF.at(0)) = m_vertices.at(i).getG1Basis().patch(j).coefs().at(k);
                }

            }
            gsInfo << "\n";

        }
    } // vertex



}

template<class T>
void gsG1System<T>::plotParaview(gsMultiPatch<T> & mp,
                                 std::vector<gsG1AuxiliaryPatch> & interface,
                                 std::vector<gsG1AuxiliaryPatch> & boundaries,
                                 std::vector<gsG1AuxiliaryPatch> & vertices,
                                 std::string basename)
{
    gsParaviewCollection collection(basename);
    std::string fileName;
    index_t iter = 0;
    for (std::vector<gsG1AuxiliaryPatch>::iterator auxGeo = interface.begin(); auxGeo != interface.end(); ++auxGeo)
    {
        for (size_t i = 2; i < auxGeo->getG1Basis().nPatches()-2; i++)
        {
            fileName = basename + "_" + util::to_string(iter);
            gsField<> temp_field(mp.patch(auxGeo->getGlobalPatchIndex()),auxGeo->getG1Basis().patch(i));
            gsWriteParaview(temp_field,fileName,5000);
            collection.addTimestep(fileName,iter,"0.vts");
            iter ++;
        }
    }
    for (std::vector<gsG1AuxiliaryPatch>::iterator auxGeo = boundaries.begin(); auxGeo != boundaries.end(); ++auxGeo)
    {
        for (size_t i = 0; i < auxGeo->getG1Basis().nPatches(); i++)
        {
            fileName = basename + "_" + util::to_string(iter);
            gsField<> temp_field(mp.patch(auxGeo->getGlobalPatchIndex()),auxGeo->getG1Basis().patch(i));
            gsWriteParaview(temp_field,fileName,5000);
            collection.addTimestep(fileName,iter,"0.vts");
            iter ++;
        }
    }
    for (std::vector<gsG1AuxiliaryPatch>::iterator auxGeo = vertices.begin(); auxGeo != vertices.end(); ++auxGeo)
    {
        for (size_t i = 0; i < auxGeo->getG1Basis().nPatches(); i++)
        {
            fileName = basename + "_" + util::to_string(iter);
            gsField<> temp_field(mp.patch(auxGeo->getGlobalPatchIndex()),auxGeo->getG1Basis().patch(i));
            gsWriteParaview(temp_field,fileName,5000);
            collection.addTimestep(fileName,iter,"0.vts");
            iter ++;
        }
    }

    collection.save();
}

template<class T>
void gsG1System<T>::plotParaview(gsMultiPatch<T> & mp,
                                 std::vector<gsG1AuxiliaryPatch> & edges,
                                 std::vector<gsG1AuxiliaryPatch> & vertices,
                                 std::string basename)
{
    gsParaviewCollection collection(basename);
    std::string fileName;
    index_t iter = 0;
    /*
    for (std::vector<gsG1AuxiliaryPatch>::iterator auxGeo = edges.begin(); auxGeo != edges.end(); ++auxGeo)
    {
        for (size_t i = 0; i < auxGeo->getG1Basis().nPatches(); i++)
        {
            fileName = basename + "_" + util::to_string(iter);
            gsField<> temp_field(mp.patch(auxGeo->getGlobalPatchIndex()),auxGeo->getG1Basis().patch(i));
            gsWriteParaview(temp_field,fileName,5000);
            collection.addTimestep(fileName,iter,"0.vts");
            iter ++;
        }
    }
     */
    for (std::vector<gsG1AuxiliaryPatch>::iterator auxGeo = vertices.begin(); auxGeo != vertices.end(); ++auxGeo)
    {
        for (size_t i = 0; i < auxGeo->getG1Basis().nPatches(); i++)
        {
            fileName = basename + "_" + util::to_string(iter);
            gsField<> temp_field(mp.patch(auxGeo->getGlobalPatchIndex()),auxGeo->getG1Basis().patch(i));
            gsWriteParaview(temp_field,fileName,5000);
            collection.addTimestep(fileName,iter,"0.vts");
            iter ++;
        }
    }

    collection.save();
}

} // namespace
