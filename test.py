import numpy as np
import smd_3d as smd

outer_radius = 1
ejection_length = 0.1
desorption_time = 10
boundary_coefficient = 1
bulk_coefficient = 0.01
theta = 2*np.pi/3
traps = np.array([
                     [0.,0.,0.,0.1]
                 ])
r = 0.7
initial_position = [r,0,0]
nb_simulations = 1000
tolerance = 0.001

# split = smd.splitting_boundary_targets_python(
#     outer_radius,
#     ejection_length,
#     desorption_time,
#     boundary_coefficient,
#     traps,
#     initial_position,
#     nb_simulations,
#     tolerance)

mfpt = smd.mfpt_python(
    outer_radius,
    bulk_coefficient,
    boundary_coefficient,
    ejection_length,
    desorption_time,
    traps, 
    initial_position, 
    nb_simulations,
    tolerance)

#print(f"split : {split}")
print(f"mfpt : {mfpt}")
