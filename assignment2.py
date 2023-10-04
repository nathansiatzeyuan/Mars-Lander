import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11  # Universal Gravitational Constant (m^3/kg/s^2)
M = 6.42e23      # Mass of Mars (kg)
m = 1000         # Mass of the moving body (kg)

#Intial Position
altitude_circular_orbit = 1000000
mars_radius = 3397000
radius_circular_orbit = altitude_circular_orbit + mars_radius  # Radius of Mars + Altitude

#Initial velocity for scenario 2,3,4
# Scenario 2: Circular Orbit (Orbital Trajectory)
# Calculate the required circular orbital speed
v_circular_orbit = np.sqrt(G * M / radius_circular_orbit)
print(v_circular_orbit)

# Scenario 3: Elliptical Orbit (Orbital Trajectory)
# Calculate the initial velocity (v_0) using the vis-viva equation
v_elliptical_orbit = 4000

# Scenario 4: Hyperbolic Escape (Orbital Trajectory)
v_escape = np.sqrt(2 * G * M / radius_circular_orbit)
print(v_escape)

# Common parameters for all scenarios
t_max = 10000
dt = 0.1  # Time step (seconds)
t_array = np.arange(0, t_max, dt)

# Lists to store position and velocity at each time step for each scenario
x_lists = [[] for _ in range(4)]
v_lists = [[] for _ in range(4)]

# Scenarios 2, 3, 4: Initial positions and velocities
x_values = [
    np.array([0, radius_circular_orbit, 0]),  # Descend
    np.array([0, radius_circular_orbit, 0]),  # Circular Orbit
    np.array([0, radius_circular_orbit, 0]),  # Elliptical Orbit
    np.array([0, radius_circular_orbit, 0])   # Hyperbolic Escape
]

v_values = [
    np.array([0, 0, 0]),                      # Descend
    np.array([v_circular_orbit, 0, 0]),       # Circular Orbit
    np.array([v_elliptical_orbit, 0, 0]),    # Elliptical Orbit
    np.array([v_escape, 0, 0])               # Hyperbolic Escape
]

# Numerical integration (Euler's method) for each scenario
for scenario in range(4):
    x = x_values[scenario]
    v = v_values[scenario]

    x_list = []
    v_list = []

    for t in t_array:
        # Append current state to trajectories (Euler method)
        x_list.append(x.copy())
        v_list.append(v.copy())

        # Calculate the gravitational force
        r = np.array(x_list[-1])  # Convert to NumPy array
        r_magnitude = np.linalg.norm(r)
        F_magnitude = (-G * M * m) / (r_magnitude ** 3)
        F = F_magnitude * r
        
        # Calculate new position and velocity (Euler method)
        a = F / m
        v = v + dt * a
        x = x + dt * v

        if scenario == 0 and x[1] < mars_radius:
            x[1] = mars_radius
            v[1] = 0  # Velocity becomes zero when reaching the surface

    # Store trajectory lists for each scenario
    x_lists[scenario] = np.array(x_list)
    v_lists[scenario] = np.array(v_list)

# Plot all trajectories on the same graph
plt.figure(1)
plt.clf()
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.grid()

# Plot Mars as a circle (orbital plane)
mars_circle = plt.Circle((0, 0), mars_radius, color='red', alpha=0.7)
plt.gca().add_patch(mars_circle)

# Scenarios 2, 3, 4: Plot Trajectories in the x-y plane
label = ['Descent','Circular Orbit', 'Elliptical Orbit', 'Escape Velocity']
for scenario in range(4):
    plt.plot(x_lists[scenario][:, 0], x_lists[scenario][:, 1], label=label[scenario])

plt.legend()
plt.show()


## Verlet integrator

# Lists to store position and velocity at each time step for each scenario
x_lists_verlet = [[] for _ in range(4)]
v_lists_verlet = [[] for _ in range(4)]

# Verlet Integrator for each scenario
for scenario in range(4):
    x = x_values[scenario].astype(float)
    v = v_values[scenario].astype(float)

    x_list = []
    v_list = []

    for t in t_array:
        # Append current state to trajectories (Verlet method)
        x_list.append(x.copy())
        v_list.append(v.copy())

        # Calculate the gravitational force
        r = np.array(x_list[-1])  # Convert to NumPy array
        r_magnitude = np.linalg.norm(r)
        F_magnitude = (-G * M * m) / (r_magnitude ** 3)
        F = F_magnitude * r

        # Calculate new position (Verlet method)
        a = F / m
        x += v * dt + 0.5 * a * dt ** 2

        # Calculate new velocity (Verlet method)
        r_new = np.array(x)
        r_magnitude_new = np.linalg.norm(r_new)
        F_magnitude_new = (-G * M * m) / (r_magnitude_new ** 3)
        F_new = F_magnitude_new * r_new
        a_new = F_new / m
        v += 0.5 * (a + a_new) * dt

        if scenario == 0 and x[1] < mars_radius:
            x[1] = mars_radius
            v[1] = 0  # Velocity becomes zero when reaching the surface

    # Store trajectory lists for each scenario
    x_lists_verlet[scenario] = np.array(x_list)
    v_lists_verlet[scenario] = np.array(v_list)

plt.figure(2)
plt.clf()
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.grid()

# Plot Mars as a circle (orbital plane)
mars_circle = plt.Circle((0, 0), mars_radius, color='red', alpha=0.7)
plt.gca().add_patch(mars_circle)
    
label = ['Descent', 'Circular Orbit', 'Elliptical Orbit', 'Escape Velocity']
# Plot each scenario separately using Verlet integrator
for scenario in range(4):
    plt.plot(x_lists_verlet[scenario][:, 0], x_lists_verlet[scenario][:, 1], label=label[scenario])

plt.legend()
plt.show()

