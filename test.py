import numpy as np
import smd_3d as smd

outer_radius = 1
ejection_length = 0.1
desorption_time = 1
boundary_coefficient = 1
theta = 2*np.pi/3
traps = np.array([
                     [np.cos(theta),np.sin(theta),-0.2,0.1]
                     ,[np.cos(2*theta),np.sin(2*theta),-0.2,0.1]
                     ,[np.cos(3*theta),np.sin(3*theta),-0.2,0.1]
                 ])

initial_position = [0,0,0]
target = 1
nb_simulations = 1000
tolerance = 0.01

split = smd.splitting_targets_python(
        outer_radius,
        ejection_length,
        desorption_time,
        boundary_coefficient,
        traps, 
        initial_position, 
        target,
        nb_simulations,
        tolerance)
print(f"should be close to 1/3 : {split}")
