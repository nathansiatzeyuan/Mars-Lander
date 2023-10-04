# uncomment the next line if running in a notebook
# %matplotlib inline
import numpy as np
import matplotlib.pyplot as plt

# mass, spring constant, initial position and velocity (Euler method)
m = 1
k = 1
x = 0
v = 1

# mass, spring constant, initial position and velocity (Verlet integrator)
m1 = 1
k1 = 1
x1 = 0
v1 = 1

# simulation time, timestep and time
t_max = 100
dt = 0.1
t_array = np.arange(0, t_max, dt)

# initialise empty lists to record trajectories
x_list = []
v_list = []
x1_list = []
v1_list = []

i=0
# Euler integration
for t in t_array:

    # append current state to trajectories （Euler method)
    x_list.append(x)
    v_list.append(v)

    # calculate new position and velocity （Euler method)
    a = -k * x / m
    x = x + dt * v
    v = v + dt * a
 
    # append current state to trajectories (Verlet integrator)
    if i < 2:
        x1_list.append(x)
        v1_list.append(v) 
    else:
        # calculate new position and velocity (Verlet integrator)
        x1_position = len(x1_list)-1
        x1 = x1_list[x1_position]
        a1 = -k1 * x1 / m1
        x1 = (2*x1) - (x1_list[x1_position-1]) + (dt **2)*a1
        v1 = (x1 - (x1_list[x1_position]))/dt
        x1_list.append(x1)
        v1_list.append(v1) 
    
    i+=1

# convert trajectory lists into arrays, so they can be sliced (useful for Assignment 2) (Euler Method)
x_array = np.array(x_list)
v_array = np.array(v_list)

# convert trajectory lists into arrays, so they can be sliced (useful for Assignment 2) (Verlet integrator)
x1_array = np.array(x1_list)
v1_array = np.array(v1_list)
print(x1_list)
print(v1_list)
# plot the position-time graph
plt.figure(1)
plt.clf()
plt.xlabel('time (s)')
plt.grid()
# plt.plot(t_array, x_array, label='x (m)')
# plt.plot(t_array, v_array, label='v (m/s)')
plt.plot(t_array, x1_array, label='x1 (m)')
plt.plot(t_array, v1_array, label='v1 (m/s)')
plt.legend()
plt.show()

