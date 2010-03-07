#!/usr/bin/env python

from glob import glob
import sys
import os.path
import re
from math import sqrt, log, cos, cosh, tanh, sinh, pi

from common import *


# NB: fixed for the model
wall_to_wall_distance = 0.05
width = 0.05
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
for file in concat(map(glob, sys.argv[1:])):
    print ""
    print "#", file
    
    # pl29sf1.13base0.1distance.csv
    # r"pl(\d+)sf([\d\.]+)base([\d\.]+)distance"
    file_base = os.path.splitext(file)[0]
    
    # Extract data from CSV
    velocity_average, velocity_profile = extractdata(file)
    
    # Build data for graphing
    sim_velocity_profiles[file_base] = velocity_profile


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
plt.savefig("velocity-profile")
