/** @file gsSpectra.h

    @brief Header file for using Spectra extension

    https://spectralib.org/doc

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsConfig.h>

#include <precice/SolverInterface.hpp>

namespace gismo {

template<class T>
class gsPreCICE //gsPreciceAdapter
{
public:

    /// Check these functions
    /// write_block_scalar_data
    /// set_mesh_vertex
    /// set_mesh_vertices



    gsPreCICE()
    :
    m_interface("participantName", "configurationFileName", 0, 1)
    {
        // #ifdef GISMO_WITH_MPI --> rank and size
        // m_interface = precice::SolverInterface( "participantName", "configurationFileName", 0, 1 ); // last arguments rank, size
    }

    gsPreCICE(std::string participantName, std::string configurationFileName)
    :
    m_interface(participantName, configurationFileName, 0, 1)
    {
        // #ifdef GISMO_WITH_MPI --> rank and size
        // m_interface = precice::SolverInterface( "participantName", "configurationFileName", 0, 1 ); // last arguments rank, size
    }

    bool isCouplingOngoing() const { return m_interface.isCouplingOngoing(); }
    bool isReadDataAvailable() const { return m_interface.isReadDataAvailable(); }
    bool isWriteDataRequired(real_t computedTimestepLength) const { return m_interface.isWriteDataRequired(computedTimestepLength); }
    bool isTimeWindowComplete() const { return m_interface.isTimeWindowComplete(); }
    bool isActionRequired(const std::string &action) const { return m_interface.isActionRequired(action); }
    void markActionFulfilled(const std::string &action) { m_interface.markActionFulfilled(action); }

    // TODO: These functions are precice constants and are preferably called outside of the class
    const std::string actionWriteInitialData() { return precice::constants::actionWriteInitialData(); }
    const std::string actionWriteIterationCheckpoint() { return precice::constants::actionWriteIterationCheckpoint(); }
    const std::string actionReadIterationCheckpoint() { return precice::constants::actionReadIterationCheckpoint(); }

    T initialize()
    {
        // gsMatrix<T> pointsTranspose = points;
        // pointsTranspose.blockTransposeInPlace(1);
        // // iterate over mesh to obtain point from basis
        // // collect everything to precice format
        // m_meshID = m_interface.getMeshID(m_meshName);

        // m_sizes = points.cols();
        // m_positions = pointsTranspose.data();

        // m_interface.setMeshVertices(m_meshID,m_sizes,m_positions,m_vertexIDs);

        m_precicedt = m_interface.initialize();
        //read is already available after initialize


        // T * values;
        // if (points.rows()==1)
        //     m_interface.readBlockScalarData(m_meshID,m_sizes,m_vertexIDs,values)

        return m_precicedt;
    }

    void addMesh(const std::string & meshName, const gsMatrix<T> & points)
    {
        m_meshNames.push_back(meshName);
        m_meshIDs.push_back(m_interface.getMeshID(m_meshNames.back()));

        gsMatrix<T> pointsTranspose = points;
        pointsTranspose.blockTransposeInPlace(1);
        m_sizes.push_back(points.cols());
        m_positions.push_back(pointsTranspose.data());

        index_t * vertexIDs;
        m_interface.setMeshVertices(m_meshIDs.back(),m_sizes.back(),m_positions.back(),vertexIDs);
        m_vertexIDs.push_back(vertexIDs);
    }

    real_t advance(T dt)
    {
        /// convert values to paramettric space??

        // m_interface.writeBlockVectorData(m_meshID,m_vertexIDs,values)
        // can also do the reading and the writing
        m_precicedt = m_interface.advance(dt);
        return m_precicedt;

        // values = m_interface.readBlockScalarData(....)
        // make a gismo function object and return



        //  // optional: checkpoints could be in here, but for very complicated cases it would be better to put the checkpoints outside

    }

    void finalize()
    {
        m_interface.finalize();
    }

    void readBlockScalarData(const index_t & dataID, const index_t & size, const index_t * IDs, gsMatrix<T> & values) const
    {
        T * values_ptr = new T[size];
        m_interface.readBlockScalarData(dataID,size,IDs,values_ptr);

        values = Eigen::Map<typename gsMatrix<T>::Base>(values_ptr,1,size);

        // CAST VALUES TO gsMATRIX
    }

    void readBlockScalarData(const index_t & meshID, const index_t & dataID, const gsMatrix<T> & coords, gsMatrix<T> & values) const
    {
        index_t * IDs = new index_t[coords.cols()];
        this->getMeshVertexIDsFromPositions(meshID,coords,IDs);

        this->readBlockScalarData(dataID,coords.cols(),IDs,values);
    }

    // void readBlockVectorData(const index_t & dataID, const index_t & size, const index_t * IDs, gsMatrix<T> & values) const
    // {
    //     T * values_ptr = new T[size];
    //     m_interface.readBlockVectorData(dataID,size,IDs,values_ptr);
    //     index_t rows =
    //     values = Eigen::Map<typename gsMatrix<T>::Base>(values_ptr,1,size);

    //     // CAST VALUES TO gsMATRIX
    // }

    // void readBlockVectorData(const index_t & meshID, const index_t & dataID, const gsMatrix<T> & coords, gsMatrix<T> & values) const
    // {
    //     index_t * IDs = new index_t[coords.cols()];
    //     this->getMeshVertexIDsFromPositions(meshID,coords,IDs);

    //     this->readBlockVectorData(dataID,coords.cols(),IDs,values);
    // }

    void writeBlockScalarData(const index_t & dataID, const index_t & size, const index_t * IDs, const gsMatrix<T> & values)
    {
        m_interface.writeBlockScalarData(dataID,size,IDs,values.data());
    }

    void writeBlockScalarData(const index_t & meshID, const index_t & dataID, const gsMatrix<T> & coords, const gsMatrix<T> & values)
    {
        index_t * IDs = new index_t[coords.cols()];
        this->getMeshVertexIDsFromPositions(meshID,coords,IDs);

        this->writeBlockScalarData(dataID,coords.cols(),IDs,values);
    }

    // void writeBlockVectorData(const index_t & dataID, const index_t & size, const index_t * IDs, const gsMatrix<T> & values) const
    // {
    //     m_interface.writeBlockVectorData(dataID,size,IDs,values.data());
    // }

    // void writeBlockVectorData(const index_t & meshID, const index_t & dataID, const gsMatrix<T> & coords, const gsMatrix<T> & values) const
    // {
    //     index_t * IDs = new index_t[coords.cols()];
    //     this->getMeshVertexIDsFromPositions(meshID,coords,IDs);

    //     this->writeBlockVectorData(dataID,coords.cols(),IDs,values);
    // }

    void getMeshVertexIDsFromPositions(const index_t & meshID, const gsMatrix<T> & coords, gsMatrix<index_t> & IDs) const
    {
        index_t * ID_ptr;
        this->getMeshVertexIDsFromPositions(meshID,coords,ID_ptr);
        for (index_t k=0; k!=coords.cols(); k++)
            gsDebugVar(ID_ptr[k]);

        // IDs = Eigen::Map<typename gsMatrix<index_t>::Base >{ID_ptr,1,static_cast<EIGEN_DEFAULT_DENSE_INDEX_TYPE>(coords.cols())};
        IDs = Eigen::Map<typename gsMatrix<index_t>::Base>(ID_ptr,1,coords.cols());
    }

    void getMeshVertexIDsFromPositions(const index_t & meshID, const gsMatrix<T> & coords, index_t * ID_ptr) const
    {
        gsMatrix<T> coordsTranspose = coords;
        coordsTranspose.blockTransposeInPlace(1);
        T * positions = coordsTranspose.data();
        m_interface.getMeshVertexIDsFromPositions(meshID,coords.cols(),positions,ID_ptr);
    }

    index_t getMeshID(const std::string &dataName) const
    {
        return m_interface.getMeshID(dataName);
    }

    index_t getDataID(const std::string &dataName, int meshID) const
    {
        return m_interface.getDataID(dataName,meshID);
    }

private:
    std::vector<std::string> m_meshNames;
    std::vector<index_t> m_meshIDs;
    std::vector<index_t> m_sizes;
    std::vector<T *> m_positions;
    std::vector<index_t *> m_vertexIDs;
    std::string m_meshName;
    index_t m_meshID;
    index_t m_meshSize;
    // T* m_positions;
    // index_t* m_meshIDs;

    precice::SolverInterface m_interface;
    T m_precicedt;
};

} //namespace gismo
