/** @file gsG1Mapper.h

    @brief Provides the mapper for the G1 Basis

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
class gsG1Mapper_pascal
{
public:
    gsG1Mapper_pascal(const gsMultiPatch<T> & mp)
    {
        m_4np = 4 * mp.nPatches();
        m_global_edge.setZero(m_4np);
        m_global_vertex.setZero(m_4np);
        m_isBoundary_edge.setZero(m_4np);
        m_isBoundary_vertex.setZero(m_4np);
        m_isCoupled_vertex.setZero(m_4np);
    }

    gsDofMapper getMapper_Vertex() {return map_vertex; }
    gsDofMapper getMapper_Edges() {return map_edges; }

    void markGlobalIndex_Edge(index_t row,
                         index_t local_side,
                         index_t patchIndex,
                         bool isBoundary)
    {
        m_global_edge.at(row) = patchIndex*4 + local_side;
        m_isBoundary_edge.at(row) = isBoundary;
    }

    index_t localToGlobal_Edge(index_t local_side,
                          index_t patchIndex)
    {
        for (index_t i = 0; i < m_4np; i++)
            if (m_global_edge.at(i) == patchIndex*4 + local_side)
                return i;

        gsInfo << "ERROR IN G1Mapper\n";
        return -1;
    }

    bool isBoundary_Edge(index_t local_side,
                         index_t patchIndex)
    {
        for (index_t i = 0; i < m_4np; i++)
            if (m_global_edge.at(i) == patchIndex*4 + local_side)
                return m_isBoundary_edge.at(i);

        gsInfo << "ERROR IN G1Mapper\n";
        return -1;
    }

    std::pair<index_t,index_t> findVertexOfEdge(index_t local_side)
    {
        std::pair<index_t,index_t> vertex_pair;
        switch (local_side)
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
        return vertex_pair;
    }

    std::pair<index_t,index_t> findEdgesOfVertex(index_t local_vertex)
    {
        std::pair<index_t,index_t> vertex_pair;
        switch (local_vertex)
        {
            case 1:
                vertex_pair = std::make_pair(3,1);
                break;
            case 2:
                vertex_pair = std::make_pair(3,2);
                break;
            case 3:
                vertex_pair = std::make_pair(4,1);
                break;
            case 0:
                vertex_pair = std::make_pair(4,2);
                break;
            default:
                break;
        }
        return vertex_pair;
    }

    void markGlobalIndex_Vertex(index_t row,
                              index_t local_side,
                              index_t patchIndex,
                              index_t isBoundary)
    {
        m_global_vertex.at(row) = patchIndex*4 + local_side;
        m_isBoundary_vertex.at(row) = isBoundary;
    }

    void coupleVertex(std::vector<size_t> patchIndex,
                      std::vector<size_t> vertIndex)
    {
        index_t val = localToGlobal_Vertex(vertIndex.at(0),patchIndex.at(0));
        for (index_t i = 0; i < patchIndex.size(); i++)
            m_isCoupled_vertex.at(localToGlobal_Vertex(vertIndex.at(i),patchIndex.at(i))) = val;
    }

    std::vector<index_t> findAllCoupledVertex(index_t val)
    {
        std::vector<index_t> coupled;
        for (index_t i = 0; i < m_isCoupled_vertex.size(); i++)
        {
            if (m_isCoupled_vertex.at(i) == val)
                coupled.push_back(i);
        }
        return coupled;
    }

    index_t localToGlobal_Vertex(index_t local_side,
                               index_t patchIndex)
    {
        for (index_t i = 0; i < m_4np; i++)
            if (m_global_vertex.at(i) == patchIndex*4 + local_side)
                return i;

        gsInfo << "ERROR IN G1Mapper\n";
        return -1;
    }

    index_t whereGlobalIndex_Vertex(index_t global)
    {
        for (index_t i = 0; i < m_4np; i++)
            if (m_global_vertex.at(i) == global)
                return i;

        gsInfo << "ERROR IN G1Mapper\n";
        return -1;
    }

    index_t whichBoundary_Vertex(index_t local_side,
                         index_t patchIndex)
    {
        for (index_t i = 0; i < m_4np; i++)
            if (m_global_vertex.at(i) == patchIndex*4 + local_side)
                return m_isBoundary_vertex.at(i);

        gsInfo << "ERROR IN G1Mapper\n";
        return -1;
    }

    void markBoundary_Edge(std::vector<gsG1AuxiliaryPatch> & g1_edges, gsMultiPatch<> & mp)
    {
        gsVector<index_t> sz;
        sz.resize(m_4np);
        for (index_t i = 0; i < m_4np; i++)
            sz.at(i) = g1_edges.at(i).getG1Basis().nPatches();

        gsDofMapper map(sz);

        for (index_t i = 0; i < m_4np; i++)
            if (m_isBoundary_edge.at(i))
            {
                gsVector<unsigned> vec;
                if (g1_edges.at(i).get_n_plus() == 8)
                    vec.setLinSpaced(g1_edges.at(i).get_n_plus()-6,3,g1_edges.at(i).get_n_plus()-4);
                else
                    vec.setLinSpaced(g1_edges.at(i).get_n_plus()-6,3,g1_edges.at(i).get_n_plus()-3);
                gsInfo << vec << " Plus " << g1_edges.at(i).get_n_plus() << " Minus " << g1_edges.at(i).get_n_minus() << "\n";
                gsMatrix<unsigned> bdr(g1_edges.at(i).get_n_plus()-6,1);
                bdr.col(0) = vec;
                map.markBoundary(i,bdr);
            }

        for (index_t i = 0; i < mp.interfaces().size(); i++)
        {
            gsInfo << " Interface " << mp.interfaces()[i].first().patch << " : " <<
               mp.interfaces()[i].second().patch << "\n";
            index_t global_patch_0 = mp.interfaces()[i].first().patch;
            index_t local_side_0 = mp.interfaces()[i].first().side();
            index_t global_patch_1 = mp.interfaces()[i].second().patch;
            index_t local_side_1 = mp.interfaces()[i].second().side();

            gsVector<unsigned> vec;
            vec.setLinSpaced(g1_edges.at(localToGlobal_Edge(local_side_0,global_patch_0)).getG1Basis().nPatches(),0,g1_edges.at(localToGlobal_Edge(local_side_0,global_patch_0)).getG1Basis().nPatches());
            gsMatrix<unsigned> bdr(g1_edges.at(localToGlobal_Edge(local_side_0,global_patch_0)).getG1Basis().nPatches(),1);
            bdr.col(0) = vec;
            map.matchDofs(localToGlobal_Edge(local_side_1,global_patch_1), bdr,localToGlobal_Edge(local_side_0,global_patch_0), bdr);
        }

        map.finalize();
        map_edges = map;
        gsInfo << "Mapper Edges \n";
        map.print();
    }

    void markBoundary_Vertex(gsMultiPatch<> & mp, std::vector<gsG1AuxiliaryPatch> & g1_vertex)
    {
        gsInfo << " PRINT : " << m_global_vertex << "\n";

        gsVector<index_t> sz;
        sz.resize(m_4np);
        for (index_t i = 0; i < m_4np; i++)
            sz.at(i) = g1_vertex.at(i).getG1Basis().nPatches();

        gsDofMapper map(sz);

        for (index_t i = 0; i < m_4np; i++)
        {
            if (m_isBoundary_vertex.at(i) == -1) // Boundary vertex
            {
                gsVector<unsigned> vec;
                vec.setLinSpaced(5,0,5);
                gsMatrix<unsigned> bdr(5,1);
                //bdr.col(0) = vec;
                bdr << 0,1,2,3,5;
                map.markBoundary(i,bdr);
            }
            else if (m_isBoundary_vertex.at(i) == 0) // Internal vertex
            {

            }
            else if (m_isBoundary_vertex.at(i) == 1) // Interface-Boundary vertex
            {
                std::vector<bool> tmp;
                gsMatrix<unsigned> bdr(3,1);
                index_t global_patch = (m_global_vertex.at(i)-1) /4;
                index_t local_vertIdx = m_global_vertex.at(i)%4;
                std::pair<index_t,index_t> edges = findEdgesOfVertex(local_vertIdx);

                gsInfo << "Patch " << global_patch << " : " << local_vertIdx << " : " << m_isBoundary_edge.at(localToGlobal_Edge(edges.second,global_patch)) << "\n";
                if (m_isBoundary_edge.at(localToGlobal_Edge(edges.second,global_patch))) // West or East
                    bdr << 0,2,5;
                else
                    bdr << 0,1,3;

                map.markBoundary(i,bdr);
            }
        }

        gsInfo << m_global_vertex << "\n";
        for (index_t val = 0; val < m_isCoupled_vertex.maxCoeff(); val ++)
        {
            std::vector<index_t> coupled = findAllCoupledVertex(val);
            if (coupled.size() > 1)
            {
                gsInfo << whereGlobalIndex_Vertex(m_global_vertex.at(coupled.at(0))) << "\n";
                gsVector<unsigned> vec;
                vec.setLinSpaced(6,0,6);
                gsMatrix<unsigned> bdr(6,1);
                bdr.col(0) = vec;

                index_t initial_patchID = whereGlobalIndex_Vertex(m_global_vertex.at(coupled.at(0)));
                for (index_t i = 1; i < coupled.size(); i++)
                {
                    gsInfo << whereGlobalIndex_Vertex(m_global_vertex.at(coupled.at(i))) << "\n";
                    map.matchDofs(initial_patchID, bdr, whereGlobalIndex_Vertex(m_global_vertex.at(coupled.at(i))), bdr);
                }
            }
        }

        map.finalize();
        map_vertex = map;
        map.print();
    }

protected:
    index_t m_4np;
    gsVector<index_t> m_global_edge, m_global_vertex;
    gsVector<index_t> m_isBoundary_edge, m_isBoundary_vertex;
    gsVector<index_t> m_isCoupled_vertex;

    gsDofMapper map_edges, map_vertex;

}; // Class gsG1Mapper
} // namespace gismo
