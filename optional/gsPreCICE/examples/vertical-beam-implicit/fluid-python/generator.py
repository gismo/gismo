import numpy as np
import time
import precice


n = 20
dn = 1 / n
timestepsize = 0.1

# generate mesh
y = np.linspace(0, 1, n + 1)

# preCICE setup
configFileName = "../precice-config.xml"
participantName = "Fluid"
solverProcessIndex = 0
solverProcessSize = 1
participant = precice.Participant(participantName, configFileName, solverProcessIndex, solverProcessSize)

# define coupling meshes
meshName = "Fluid-Mesh"
# meshID = participant.get_mesh_id(meshName)

positions = [[-0.5, y0] for y0 in y[:-1]]
vertex_ids = participant.set_mesh_vertices(meshName, positions)

print(positions)

# coupling data
writeData = "Stress"
readData = "Displacement"
# writedataID = participant.get_data_id(writeData, meshID)
# readdataID = participant.get_data_id(readData, meshID)

# initialize preCICE
participant.initialize()
precice_dt = participant.get_max_time_step_size()
dt = min(precice_dt, timestepsize)

timestep = 0
t = 0

while participant.is_coupling_ongoing():
    # read displacements from participant
    v = participant.read_data(meshName, readData, vertex_ids, dt)

    # save checkpoint
    if participant.requires_writing_checkpoint():
        print("DUMMY: Writing iteration checkpoint")

    # solve fluid equations
    print('t = ',t)
    # u = [[-1e4 * np.sin(3.1415923565*t), 0] for y0 in y[:-1]]
    u = [[-1e4, 0] for y0 in y[:-1]]

    # if (t==0):
    #     u = [[-1e4, 0] for y0 in y[:-1]]
    # else:
    #     u = [[0, 0] for y0 in y[:-1]]

    # write forces to participant
    participant.write_data(meshName, writeData, vertex_ids, u)

    # do the coupling
    participant.advance(dt)
    precice_dt = participant.get_max_time_step_size()

    dt = min(precice_dt, timestepsize)

    # read checkpoint if required
    if participant.requires_reading_checkpoint():
        print("DUMMY: Reading iteration checkpoint")

    if participant.is_time_window_complete():
        # advance variables
        timestep += 1
        t += dt

participant.finalize()



# # preCICE setup
# participant_name = "Fluid"
# config_file_name = "../precice-config.xml"
# solver_process_index = 0
# solver_process_size = 1
# interface = precice.Interface(participant_name, config_file_name, solver_process_index, solver_process_size)

# mesh_name = "Fluid-Mesh"
# mesh_id = interface.get_mesh_id(mesh_name)

# data_name = "Data"
# data_id = interface.get_data_id(data_name, mesh_id)

# positions = [[1, y0] for y0 in y[:-1]]

# vertex_ids = interface.set_mesh_vertices(mesh_id, positions)

# precice_dt = interface.initialize()

# t = 0

# while interface.is_coupling_ongoing():

#     dt = 0.01

#     print("Generating data")
#     dt = np.minimum(dt, precice_dt)
#     time.sleep(0.2)
#     u = 1 - 2 * np.random.rand(n)

#     interface.write_block_scalar_data(data_id, vertex_ids, u)

#     precice_dt = interface.advance(dt)

#     t = t + dt

# interface.finalize()