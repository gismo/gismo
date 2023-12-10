/** @file gsPreCICE.h

    @brief Header file for using gsPreCICE extension

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-...)
*/

#pragma once

#include <gsCore/gsConfig.h>
#include <precice/SolverInterface.hpp>

namespace gismo {

template<class T>
class gsPreCICE
{
public:

    /**
     * @brief      Default constructor
     */
    gsPreCICE()
    :
    m_interface("participantName", "configurationFileName", 0, 1)
    {
        // #ifdef GISMO_WITH_MPI --> rank and size
        // m_interface = precice::SolverInterface( "participantName", "configurationFileName", 0, 1 ); // last arguments rank, size
    }

    /**
     * @brief      Constructor
     *
     * @param[in]  participantName        The participant name
     * @param[in]  configurationFileName  The configuration file name
     *
     * todo: add GISMO_WITH_MPI
     */
    gsPreCICE(std::string participantName, std::string configurationFileName)
    :
    m_interface(participantName, configurationFileName, 0, 1) // last arguments rank, size
    {
    }

    /// See precice::SolverInterface::isCouplingOngoing
    bool isCouplingOngoing() const { return m_interface.isCouplingOngoing(); }
    /// See precice::SolverInterface::isReadDataAvailable
    bool isReadDataAvailable() const { return m_interface.isReadDataAvailable(); }
    /// See precice::SolverInterface::isWriteDataRequired
    bool isWriteDataRequired(real_t computedTimestepLength) const { return m_interface.isWriteDataRequired(computedTimestepLength); }
    /// See precice::SolverInterface::isTimeWindowComplete
    bool isTimeWindowComplete() const { return m_interface.isTimeWindowComplete(); }
    /// See precice::SolverInterface::isActionRequired
    bool isActionRequired(const std::string &action) const { return m_interface.isActionRequired(action); }
    /// See precice::SolverInterface::markActionFulfilled
    void markActionFulfilled(const std::string &action) { m_interface.markActionFulfilled(action); }

    // TODO: These functions are precice constants and are preferably called outside of the class
    /// See precice::SolverInterface::actionWriteInitialData
    const std::string actionWriteInitialData() { return precice::constants::actionWriteInitialData(); }
    /// See precice::SolverInterface::actionWriteIterationCheckpoint
    const std::string actionWriteIterationCheckpoint() { return precice::constants::actionWriteIterationCheckpoint(); }
    /// See precice::SolverInterface::actionReadIterationCheckpoint
    const std::string actionReadIterationCheckpoint() { return precice::constants::actionReadIterationCheckpoint(); }


    /**
     * @brief      Initializes the precice::SolverInterface
     *
     * @note        This function can be expanded with more initialization actions (i.e. an initial write/read)
     *
     * @return     the precice time-step
     */
    T initialize()
    {
        m_precicedt = m_interface.initialize();
        return m_precicedt;
    }

    void initialize_data()
    {
        m_interface.initializeData();
    }



    /**
     * @brief      Adds a mesh to the precice interface
     *
     * @param[in]  meshName  The mesh name
     * @param[in]  points    The points stored in columns
     */
    void addMesh(const std::string & meshName, const gsMatrix<T> & points)
    {
        m_meshNames.push_back(meshName);
        m_meshIDs.push_back(m_interface.getMeshID(m_meshNames.back()));

        // Makes a hard-copy? Can we prevent it?
        gsMatrix<T> pointsTranspose = points;
        pointsTranspose.blockTransposeInPlace(1);
        m_sizes.push_back(points.cols());
        m_positions.push_back(pointsTranspose);

        gsVector<index_t> vertexIDs(m_sizes.back());
        m_interface.setMeshVertices(m_meshIDs.back(),m_sizes.back(),m_positions.back().data(),vertexIDs.data());
        m_vertexIDs.push_back(vertexIDs);
    }

    /**
     * @brief      Advances the precice::SolverInterface
     *
     * @param[in]  dt    time step
     *
     * @note        This function can be expanded with a read/write action inside.
     *              Checkpoints can also be put inside, but it might not be useful in complicated cases
     *
     * @return     the precice timestep
     */
    T advance(T dt)
    {
        m_precicedt = m_interface.advance(dt);
        return m_precicedt;
    }

    /// See precice::SolverInterface::finalize for details
    void finalize()
    {
        m_interface.finalize();
    }

    /**
     * @brief      Reads a block of scalar data.
     *
     * @param[in]  dataID  The data ID
     * @param[in]  size    The number of points
     * @param[in]  IDs     The IDs of the points
     * @param      values  The values (column-wise)
     */
    void readBlockScalarData(const index_t & dataID, const index_t & size, const gsVector<index_t> & IDs, gsMatrix<T> & values) const
    {
        values.resize(1,size);
        m_interface.readBlockScalarData(dataID,size,IDs.data(),values.data());
    }

    /**
     * @brief      Reads a block of scalar data.
     *
     * @param[in]  meshID  The mesh ID
     * @param[in]  dataID  The data ID
     * @param[in]  coords  The coordinates of the points (column-wise)
     * @param      values  The values (column-wise)
     */
    void readBlockScalarData(const index_t & meshID, const index_t & dataID, const gsMatrix<T> & coords, gsMatrix<T> & values) const
    {
        gsVector<index_t> IDs;
        this->getMeshVertexIDsFromPositions(meshID,coords,IDs);
        this->readBlockScalarData(dataID,coords.cols(),IDs,values);
    }

    void readBlockVectorData(const index_t & dataID, const index_t & size, const gsVector<index_t> & IDs, gsMatrix<T> & values) const
    {
        int d = m_interface.getDimensions();
        values.resize(d,size);
        m_interface.readBlockVectorData(dataID,size,IDs.data(),values.data());
    }

    gsMatrix<T>  readBlockVectorData(const index_t & meshID, const index_t & dataID, const gsMatrix<T> & coords) const
    {
        gsMatrix<T> result;
        this->readBlockVectorData(meshID,dataID,coords,result);
        return result;
    }

    void readBlockVectorData(const index_t & meshID, const index_t & dataID, const gsMatrix<T> & coords, gsMatrix<T> & values) const
    {
        gsVector<index_t> IDs;
        this->getMeshVertexIDsFromPositions(meshID,coords,IDs);
        this->readBlockVectorData(dataID,coords.cols(),IDs,values);
    }

    /**
     * @brief      Writes a block of scalar data.
     *
     * @param[in]  dataID  The data ID
     * @param[in]  size    The number of points
     * @param[in]  IDs     The IDs of the points
     * @param[in]  values  The values (column-wise)
     */
    void writeBlockScalarData(const index_t & dataID, const index_t & size, const gsVector<index_t> & IDs, const gsMatrix<T> & values)
    {
        m_interface.writeBlockScalarData(dataID,size,IDs.data(),values.data());
    }

    /**
     * @brief      Writes a block of scalar data.
     *
     * @param[in]  meshID  The mesh ID
     * @param[in]  dataID  The data ID
     * @param[in]  coords  The coordinates of the points (column-wise)
     * @param      values  The values (column-wise)
     */
    void writeBlockScalarData(const index_t & meshID, const index_t & dataID, const gsMatrix<T> & coords, const gsMatrix<T> & values)
    {
        gsVector<index_t> IDs;
        this->getMeshVertexIDsFromPositions(meshID,coords,IDs);
        this->writeBlockScalarData(dataID,coords.cols(),IDs,values);
    }

    void writeBlockVectorData(const index_t & dataID, const index_t & size, const gsVector<index_t> & IDs, const gsMatrix<T> & values)
    {
        m_interface.writeBlockVectorData(dataID,size,IDs.data(),values.data());
    }

    void writeBlockVectorData(const index_t & meshID, const index_t & dataID, const gsMatrix<T> & coords, const gsMatrix<T> & values)
    {
        gsVector<index_t> IDs;
        this->getMeshVertexIDsFromPositions(meshID,coords,IDs);
        this->writeBlockVectorData(dataID,coords.cols(),IDs,values);
    }

    /**
     * @brief      Gets the mesh vertex IDs from positions of the points.
     *
     * @param[in]  meshID  The mesh ID
     * @param[in]  coords  The coordinates of the points
     * @param      IDs     The IDs of the points
     */
    void getMeshVertexIDsFromPositions(const index_t & meshID, const gsMatrix<T> & coords, gsVector<index_t> & IDs) const
    {
        IDs.resize(coords.cols());
        gsMatrix<T> coordsTranspose = coords;
        coordsTranspose.blockTransposeInPlace(1);
        m_interface.getMeshVertexIDsFromPositions(meshID,coords.cols(),coordsTranspose.data(),IDs.data());
    }

    /**
     * @brief      Returns the ID of the mesh
     *
     * @param[in]  dataName  The name of the data
     *
     * @return     The mesh ID.
     */
    index_t getMeshID(const std::string &dataName) const
    {
        return m_interface.getMeshID(dataName);
    }

    /**
     * @brief      Returns the ID of the data
     *
     * @param[in]  dataName  The name of the data
     * @param[in]  meshID    The mesh ID
     *
     * @return     The data ID.
     */
    index_t getDataID(const std::string &dataName, int meshID) const
    {
        return m_interface.getDataID(dataName,meshID);
    }

    /// See \a gsFunction
    virtual std::ostream &print(std::ostream &os) const
    {
        os << precice::getVersionInformation()<<"\n\n";
        os << "Interface has the following meshes:\n";
        for (size_t k=0; k!=m_meshNames.size(); k++)
            gsInfo<<m_meshNames[k]<<"\t (ID="<<m_meshIDs[k]<<")\n";
        return os;
    }

private:
    /// Stores all mesh names (might be useful later)
    std::vector<std::string> m_meshNames;
    /// Stores all mesh IDs (might be useful later)
    std::vector<index_t> m_meshIDs;
    /// Stores all mesh sizes (might be useful later)
    std::vector<index_t> m_sizes;
    /// Stores all mesh positions (might be useful later)
    std::vector<gsMatrix<T>> m_positions;
    /// Stores all mesh vertex IDs (might be useful later)
    std::vector<gsVector<index_t>> m_vertexIDs;
    /// Stores the precice::SolverInterface (see the precice Doxygen for details about this class)
    precice::SolverInterface m_interface;
    /// Stores the precice timestep
    T m_precicedt;
};

} //namespace gismo
