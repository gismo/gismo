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
#include <precice/precice.hpp>

namespace gismo {

template<class T>
class gsPreCICE
{

struct mapCompare
{
    bool operator()(const gsVector<T> &a, const gsVector<T>& b) const
    {
        return std::lexicographical_compare(a.begin(),a.end(),b.begin(),b.end());
    }
};

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

    // preCICE v2.5.0
    // /// See precice::SolverInterface::isReadDataAvailable
    // bool isReadDataAvailable() const { return m_interface.isReadDataAvailable(); }

    // /// See precice::SolverInterface::requiresInitialData
    // bool requiresInitialData() const { return m_interface.requiresInitialData(); }
    // / See precice::SolverInterface::isTimeWindowComplete
    bool isTimeWindowComplete() const { return m_interface.isTimeWindowComplete(); }

    // preCICE v3.0.0
    /// See precice::SolverInterface::requiresInitialData
    bool requiresInitialData() { return m_interface.requiresInitialData(); }
    /// See precice::SolverInterface::requiresReadingCheckpoint
    bool requiresReadingCheckpoint() { return m_interface.requiresReadingCheckpoint(); }
    /// See precice::SolverInterface::requiresWritingCheckpoint
    bool requiresWritingCheckpoint() { return m_interface.requiresWritingCheckpoint(); }

    // preCICE v2.5.0
    // /// See precice::SolverInterface::isActionRequired
    // bool isActionRequired(const std::string &action) const { return m_interface.isActionRequired(action); }
    // /// See precice::SolverInterface::markActionFulfilled
    // void markActionFulfilled(const std::string &action) { m_interface.markActionFulfilled(action); }

    // TODO: These functions are precice constants and are preferably called outside of the class
    /// See precice::SolverInterface::actionWriteInitialData
    // const std::string actionWriteInitialData() { return precice::constants::actionWriteInitialData(); }
    // /// See precice::SolverInterface::actionWriteIterationCheckpoint
    // const std::string actionWriteIterationCheckpoint() { return precice::constants::actionWriteIterationCheckpoint(); }
    // /// See precice::SolverInterface::actionReadIterationCheckpoint
    // const std::string actionReadIterationCheckpoint() { return precice::constants::actionReadIterationCheckpoint(); }


    /**
     * @brief      Initializes the precice::SolverInterface
     *
     * @note        This function can be expanded with more initialization actions (i.e. an initial write/read)
     *
     * @return     the precice time-step
     */
    T initialize()
    {
        m_interface.initialize();
        m_precicedt = m_interface.getMaxTimeStepSize();
        return m_precicedt;
    }

/*
    precice v2.5.0
    void initialize_data()
    {
        m_interface.initializeData();
    }
*/


    /**
     * @brief      Adds a mesh to the precice interface
     *
     * @param[in]  meshName  The mesh name
     * @param[in]  points    The points stored in columns
     */
    void addMesh(const std::string & meshName, const gsMatrix<T> & points, gsVector<index_t> & vertexIDs)
    {
        index_t index = m_meshNames.size();
        m_meshNames[meshName] = index;

        vertexIDs.resize(points.cols());
        m_interface.setMeshVertices(meshName,points,vertexIDs);

        // Create a look-up table from points to IDs
        std::map<gsVector<T>,index_t,mapCompare> map;
        for (index_t k=0; k!=points.cols(); k++)
            map[points.col(k)] = vertexIDs.at(k);

        m_maps.push_back(map);
    }

    void addMesh(const std::string & meshName, const gsMatrix<T> & points)
    {
        gsVector<index_t> vertexIDs;
        this->addMesh(meshName,points,vertexIDs);
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
        m_interface.advance(dt);
        m_precicedt = m_interface.getMaxTimeStepSize();
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
     * @param[in]  meshName The mesh name
     * @param[in]  dataID   The data name
     * @param[in]  coords   The coordinates of the points (column-wise)
     * @param      values   The values (column-wise)
     */
    void readData(const std::string & meshName, const std::string & dataName, const gsMatrix<T> & coords, gsMatrix<T> & values) const
    {
        this->readData(meshName,dataName,coords,m_interface.getMaxTimeStepSize(),values);
    }

    void readData(const std::string & meshName, const std::string & dataName, const gsMatrix<T> & coords, const T & timestep, gsMatrix<T> & values) const
    {
        gsVector<index_t> IDs;
        this->getMeshVertexIDsFromPositions(meshName,coords,IDs);
        this->readData(meshName,dataName,IDs,timestep,values);
    }

    void readData(const std::string & meshName, const std::string & dataName, const gsVector<index_t> & IDs, gsMatrix<T> & values) const
    {
        this->readData(meshName,dataName,IDs,m_interface.getMaxTimeStepSize(),values);
    }

    void readData(const std::string & meshName, const std::string & dataName, const gsVector<index_t> & IDs, const T & timestep, gsMatrix<T> & values) const
    {
        int d = m_interface.getDataDimensions(meshName,dataName);
        values.resize(1,d*IDs.size());
        m_interface.readData(meshName,dataName,IDs,timestep,values);
        values.resize(d,IDs.size());
    }

    gsMatrix<T>  readData(const std::string & meshName, const std::string & dataName, const gsMatrix<T> & coords) const
    {
        gsMatrix<T> result;
        this->readData(meshName,dataName,coords,result);
        return result;
    }

    /**
     * @brief      Writes a block of scalar data.
     *
     * @param[in]  meshID  The mesh ID
     * @param[in]  dataID  The data ID
     * @param[in]  coords  The coordinates of the points (column-wise)
     * @param      values  The values (column-wise)
     */
    void writeData(const std::string & meshName, const std::string & dataName, const gsMatrix<T> & coords, const gsMatrix<T> & values)
    {
        gsVector<index_t> IDs;
        this->getMeshVertexIDsFromPositions(meshName,coords,IDs);
        m_interface.writeData(meshName,dataName,IDs,values);
    }

    /**
     * @brief      Writes a block of scalar data.
     *
     * @param[in]  meshID  The mesh ID
     * @param[in]  dataID  The data ID
     * @param[in]  IDs     The IDs of the points to write on
     * @param      values  The values (column-wise)
     */
    void writeData(const std::string & meshName, const std::string & dataName, const gsVector<index_t> & IDs, const gsMatrix<T> & values)
    {
        m_interface.writeData(meshName,dataName,IDs,values);
    }

    /**
     * @brief      Gets the mesh vertex IDs from positions of the points. THIS IS SLOW!!
     *
     * @param[in]  meshName The mesh name
     * @param[in]  coords   The coordinates of the points
     * @param      IDs      The IDs of the points
     */
    void getMeshVertexIDsFromPositions(const std::string & meshName, const gsMatrix<T> & coords, gsVector<index_t> & IDs) const
    {
        IDs.resize(coords.cols());
        for (index_t k=0; k!=coords.cols(); k++)
            IDs.at(k) = m_maps.at(this->getMeshID(meshName)).at(coords.col(k));
    }

    /**
     * @brief      Sets the mesh access region.
     *
     * @param[in]  meshName  The mesh name to access
     * @param[in]  bbox      The bounding box of the region
     */
    void setMeshAccessRegion(const std::string & meshName, const gsMatrix<T> & bbox)
    {
        // std::vector<T> region(bbox.data(), bbox.data() + bbox.rows() * bbox.cols());
        m_interface.setMeshAccessRegion(meshName,bbox);
    }

    /**
     * [FROM PRECICE MANUAL]
     * @brief getMeshVertexIDsAndCoordinates Iterates over the region of
     *        interest defined by bounding boxes and reads the corresponding
     *        coordinates omitting the mapping.
     *
     * @param[in]  meshName corresponding mesh name
     * @param[out] ids ids corresponding to the coordinates
     * @param[out] coordinates the coordinates associated to the \p ids and
     *             corresponding data values
     *
     * @pre ids.size() == getMeshVertexSize(meshName)
     * @pre coordinates.size() == getMeshVertexSize(meshName) * getMeshDimensions(meshName)
     *
     * @pre This function can be called on received meshes as well as provided
     * meshes. However, you need to call this function after @p initialize(),
     * if the \p meshName corresponds to a received mesh, since the relevant mesh data
     * is exchanged during the @p initialize() call.
     *
     * @see getMeshVertexSize() to get the amount of vertices in the mesh
     * @see getMeshDimensions() to get the spacial dimensionality of the mesh
     */
    void getMeshVertexIDsAndCoordinates(const std::string & meshName, gsVector<index_t> & IDs, gsMatrix<T> & coords) const
    {
        const index_t nPoints = m_interface.getMeshVertexSize(meshName);
        index_t d = m_interface.getMeshDimensions(meshName);

        IDs.resize(nPoints);
        coords.resize(d,nPoints);
        m_interface.getMeshVertexIDsAndCoordinates(meshName,IDs,coords);
    }

    /**
     * @brief      Returns the ID of the mesh
     *
     * @param[in]  dataName  The name of the data
     *
     * @return     The mesh ID.
     */
    index_t getMeshID(const std::string &meshName) const
    {
        return m_meshNames.at(meshName);
    }

    /// See \a gsFunction
    virtual std::ostream &print(std::ostream &os) const
    {
        // os << precice::getVersionInformation()<<"\n\n";
        os << "Interface has the following meshes:\n";
        for (std::map<std::string,index_t>::const_iterator name=m_meshNames.cbegin(); name!=m_meshNames.cend(); name++)
            gsInfo<<name->first<<"\t (ID="<<name->second<<")\n";
        return os;
    }

private:
    /// Stores all mesh names (might be useful later)
    std::map<std::string,index_t> m_meshNames;
    /// Stores all mesh positions (might be useful later)
    std::vector<std::vector<T>> m_positions;
    /// Stores all mesh vertex IDs (might be useful later)
    std::vector<std::vector<index_t>> m_vertexIDs;
    /// Stores a map between mesh positions and vertex IDs
    std::vector<std::map<gsVector<T>,index_t,mapCompare>> m_maps;
    /// Stores the precice::SolverInterface (see the precice Doxygen for details about this class)
    precice::Participant m_interface;
    /// Stores the precice timestep
    T m_precicedt;
};

} //namespace gismo
