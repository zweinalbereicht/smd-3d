import numpy as np
import smd_3d as smd

outer_radius = 1
ejection_length = 0.6
desorption_time = 0.1
boundary_coefficient = 1
bulk_coefficient = 1
theta = 2*np.pi/3
traps = np.array([
                     [0.,0.,0.,0.1]
                 ])
r = 0.7

initial_position = [r,0,0]
nb_simulations = 500
tolerance = 0.001
evolve_time=1000

print('starting')
a=smd.stationnary_state_python(
    outer_radius,
    bulk_coefficient,
    boundary_coefficient,
    ejection_length,
    desorption_time,
    initial_position,
    evolve_time,
    nb_simulations,
    tolerance,)
print('done')

a=np.array(a)
print(len(a[a=='boundary'])/len(a))

# split = smd.splitting_boundary_targets_python(
#     outer_radius,
#     ejection_length,
#     desorption_time,
#     boundary_coefficient,
#     traps,
#     initial_position,
#     nb_simulations,
#     tolerance)

# mfpt = smd.mfpt_python(
#     outer_radius,
#     bulk_coefficient,
#     boundary_coefficient,
#     ejection_length,
#     desorption_time,
#     traps, 
#     initial_position, 
#     nb_simulations,
#     tolerance)

