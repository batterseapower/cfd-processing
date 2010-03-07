#!/usr/bin/env python

from glob import glob
import sys
import os.path
import re
from math import sqrt, log, cos, cosh, tanh, sinh, pi

from common import *


directory = sys.argv[1]


# Parse out width
if "square" in directory:
    width = 0.05
else:
    m = re.search(r"width([\d\.]+)", directory)
    if m is None:
        print "Count not determine duct width: add widthN.M to the file name"
        sys.exit(1)
    else:
        width = float(m.group(1))

wall_to_wall_distance = 0.05
height = 0.05
length = 0.5
d_h = (width * height) / (2 * (width + height)) # Hydralic diameter

density = 1.18415
mu = 0.00001855
pressure = 0.05

c = width / 2
b = height / 2
pressure_gradient = -pressure / length


q_infinite_sum = infinite_sum(lambda n: tanh(0.5 * n * pi) / pow(n, 5))
q = pow(b, 4) / (6 * mu) * pressure_gradient * -(1 - 192 / pow(pi, 5) * q_infinite_sum)
vm = q / pow(b, 2)
z = 0

def velocity_theory(position):
    pres_infinite_sum = infinite_sum(lambda n: pow(-1, (n - 1) / 2.0) * (1 - cosh(n * pi * z / (2 * b)) / cosh(n * pi * c / (2 * b))) * cos(n * pi * (position - (wall_to_wall_distance / 2)) / (2 * b)) / pow(float(n), 3))
    return -16 * pow(b, 2) / (mu * pow(pi, 3)) * pressure_gradient * pres_infinite_sum


# Print some information about the model
printconstants("Model constants:", [("width", width), ("height", height), ("length", length), ("d_h", d_h), ("density", density), ("mu", mu), ("c", c), ("b", b), ("pressure gradient", pressure_gradient), ("q", q), ("vm", vm)])


# Pull in data
sim_velocity_profiles = {}
sim_scaled_velocity_profiles = {}
sim_scaled_nonhydralic_velocity_profiles = {}
for file in glob(os.path.join(directory, "*.csv")):
    print ""
    print "#", file
    
    # pl29sf1.13base0.1distance.csv
    # r"pl(\d+)sf([\d\.]+)base([\d\.]+)distance"
    file_base = os.path.basename(os.path.splitext(file)[0])
    
    # Extract data from CSV
    velocity_average, velocity_profile = extractdata(file)
    velocity_max = max([velocity for position, velocity in velocity_profile])
    
    # Trim data down to size (especially important if the velocity profile was built from position-indexed information)
    wall_to_center_distance = wall_to_wall_distance / 2.0
    half_velocity_profile = [(position, velocity) for position, velocity in velocity_profile if 0.0 < position <= wall_to_center_distance]
    
    # Build data for graphing
    sim_velocity_profiles[file_base] = velocity_profile
    sim_scaled_velocity_profiles[file_base] = [(velocity / velocity_average, 2 * position / d_h) for position, velocity in half_velocity_profile]
    sim_scaled_nonhydralic_velocity_profiles[file_base] = [(1 - (position / wall_to_center_distance), velocity / velocity_max) for position, velocity in half_velocity_profile]


# Find the theoretical velocity profile
all_positions = unions([set([position for position, _ in velocity_profile]) for velocity_profile in sim_velocity_profiles.values()])
theory_velocity_profile = sorted([(position, velocity_theory(position)) for position in all_positions])


# Guess at best model
print ""
print "# Simulation average squared error"
for sim, velocity_profile in sim_velocity_profiles.items():
    print "  ", sim, "=", average_error(velocity_profile, velocity_theory)


# OK, lets go:
# http://matplotlib.sourceforge.net/users/pyplot_tutorial.html
import matplotlib.pyplot as plt

# profile:
plt.clf()
plt.figure(figsize=(10, 8), dpi=80)

plt.xlabel("x")
plt.ylabel("u")

for sim, velocity_profile in [("theory", theory_velocity_profile)] + sim_velocity_profiles.items():
    positions, velocities = zip(*velocity_profile)
    plt.plot(positions, velocities, label=sim)

plt.legend(loc="lower center")
plt.savefig(os.path.join(directory, "velocity-profile"))

# scaled velocity:
plt.clf()
plt.figure(figsize=(10, 8), dpi=80)

plt.xlabel("u/u_avg")
plt.ylabel("2y/D_H")

for sim, scaled_velocity_profile in sim_scaled_velocity_profiles.items():
    scaled_velocity, scaled_y = zip(*scaled_velocity_profile)
    plt.plot(scaled_velocity, scaled_y, label=sim)

plt.legend(loc="upper left")
plt.savefig(os.path.join(directory, "scaled-velocity-profile"))

# scaled non-hydralic velocity:
plt.clf()
plt.figure(figsize=(10, 8), dpi=80)

plt.xlabel("2y/h")
plt.ylabel("v/v_max")

for sim, scaled_nonhydralic_velocity_profiles in sim_scaled_nonhydralic_velocity_profiles.items():
    scaled_position, scaled_velocity = zip(*scaled_nonhydralic_velocity_profiles)
    plt.plot(scaled_position, scaled_velocity, label=sim)

plt.legend(loc="lower left")
plt.savefig(os.path.join(directory, "scaled-nonhydralic-velocity-profile"))
