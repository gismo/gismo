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
interface = precice.Interface(participantName, configFileName, solverProcessIndex, solverProcessSize)

# define coupling meshes
meshName = "Fluid-Mesh"
meshID = interface.get_mesh_id(meshName)

positions = [[-0.5, y0] for y0 in y[:-1]]
vertex_ids = interface.set_mesh_vertices(meshID, positions)

print(positions)

# coupling data
writeData = "Force"
readData = "Displacement"
writedataID = interface.get_data_id(writeData, meshID)
readdataID = interface.get_data_id(readData, meshID)

# initialize preCICE
precice_dt = interface.initialize()
dt = min(precice_dt, timestepsize)

timestep = 0
t = 0

while interface.is_coupling_ongoing():
    # read displacements from interface
    if interface.is_read_data_available():
        v = interface.read_block_vector_data(readdataID, vertex_ids)
        print(v)

    # save checkpoint
    if interface.is_action_required(precice.action_write_iteration_checkpoint()):
        interface.mark_action_fulfilled(precice.action_write_iteration_checkpoint())

    # solve fluid equations
    print(t)
    u = [[-1e2 * np.sin(3.1415923565*t), 0] for y0 in y[:-1]]

    # if (t==0):
    #     u = [[-1e4, 0] for y0 in y[:-1]]
    # else:
    #     u = [[0, 0] for y0 in y[:-1]]

    # write forces to interface
    if interface.is_write_data_required(dt):
        interface.write_block_vector_data(writedataID, vertex_ids, u)

    # do the coupling
    precice_dt = interface.advance(dt)
    dt = min(precice_dt, timestepsize)

    # read checkpoint if required
    if interface.is_action_required(precice.action_read_iteration_checkpoint()):
        interface.mark_action_fulfilled(precice.action_read_iteration_checkpoint())

    if interface.is_time_window_complete():
        # advance variables
        timestep += 1
        t += dt

interface.finalize()



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