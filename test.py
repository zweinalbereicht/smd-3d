import numpy as np
import smd_3d as smd

outer_radius = 3
ejection_length = 0.05
desorption_time =100
boundary_coefficient = 1
bulk_coefficient = 1
theta = 2*np.pi/3
traps = np.array([[-2,0.,0.,0.1]])
stickytraps = np.array([[0.1,0.,0.,2,desorption_time, boundary_coefficient]])
initial_position = [-2,1,0] #dedans, on est bon!
nb_simulations = 500
tolerance = 0.01
dt=0.00001
evolve_time=30

# a=smd.fpt_distribution_python(
#     outer_radius,
#     bulk_coefficient,
#     boundary_coefficient,
#     ejection_length,
#     desorption_time,
#     traps,
#     initial_position,
#     nb_simulations,
#     dt)
#
# print(np.mean(a))

#
# print('starting')

# a=smd.stationnary_state_python(
#     outer_radius,
#     bulk_coefficient,
#     boundary_coefficient,
#     ejection_length,
#     desorption_time,
#     stickytraps,
#     initial_position,
#     evolve_time,
#     nb_simulations,
#     tolerance,)
#
# a=np.array(a)
# print(len(a[a=='boundary'])/len(a))

# print('done')
#

# split = smd.splitting_boundary_targets_python(
#     outer_radius,
#     ejection_length,
#     desorption_time,
#     boundary_coefficient,
#     traps,
#     initial_position,
#     nb_simulations,
#     tolerance)

mfpt = smd.mfpt_stickytraps_python(
    outer_radius,
    bulk_coefficient,
    boundary_coefficient,
    ejection_length,
    desorption_time,
    stickytraps,
    traps,
    initial_position,
    nb_simulations,
    tolerance)

print(f"mfpt : {mfpt}")

