#!/usr/bin/env python

import csv
import sys
import re
from math import sqrt, log


def drop(n, xs):
    for x in xs:
        if n > 0:
            n = n - 1
        else:
            yield x

def unions(sets):
    allset = set()
    for aset in sets:
        allset.update(aset)
    return allset


# NB: fixed for the model
wall_to_wall_distance = 0.05
width = 0.05
height = 0.05


# NB: these are for base = 0.01
wall_shear_stress = 111.9
density = 1.18415
mu = 0.00001855

# These are base-independent:
ustar = sqrt(wall_shear_stress / density)
d_v = mu / ustar
d_h = (width * height) / (2 * (width + height)) # Hydralic diameter


def yplus(wall_distance):
    return wall_distance / d_v

def uplus(velocity):
    return velocity / ustar

def uplus_theory(yplus):
    # Computation Methods for Fluid Dynamics, p298
    if yplus <= 5:
        return yplus
    else:
        k = 0.14
        b = 5
        return (1 / k) * log(yplus) + b


# Show computed values
print "Computed constants:"
for name, value in { "ustar" : ustar, "d_v" : d_v, "d_h" : d_h }.items():
    print name, "=", value


# Pull in data
pl_wall_functions = {}
pl_scaled_velocity_profiles = {}
for file in sys.argv[1:]:
    # pl29sf1.13base0.1distance.csv
    # r"pl(\d+)sf([\d\.]+)base([\d\.]+)distance"
    m = re.search(r"pl(\d+)", file)
    pl = int(m.group(1))
    print "Prism layer count", pl, "for", file
    
    # "Position [0.0, 1.0, 0.0] (m)-point 1 (m)","Velocity: Magnitude-point 1 (m/s)"
    # 0.0,35.15489668712291
    # 2.499999993688107E-7,35.15489668712291
    reject_count = 0
    raw_velocity_positions = []
    for row in drop(1, csv.reader(open(file))):
        # Verify record is correct length
        if len(row) != 2:
            reject_count = reject_count + 1
            continue
        
        # Verify we can parse as numbers
        try:
            raw_velocity_positions.append((float(row[0]), float(row[1])))
        except ValueError:
            reject_count = reject_count + 1
    
    # Tell user about problems with the input data
    if reject_count != 0:
        print "Rejected", reject_count, "malformed row(s)"

    # Build the velocity profile
    velocity_positions = {}
    for position, velocity in raw_velocity_positions:
        # Record the minimum position that experienced a given velocity
        velocity_positions[velocity] = min(velocity_positions.get(velocity, sys.maxint), position)
    
    velocity_profile = sorted([(position, velocity) for velocity, position in list(velocity_positions.items())])
    
    # Determine average velocity
    distance_accumulator = 0.0
    velocity_average_accumulator = 0.0
    last_position = 0.0
    for position, velocity in sorted(raw_velocity_positions):
        distance = abs(position - last_position)
        distance_accumulator = distance_accumulator + distance
        velocity_average_accumulator = velocity_average_accumulator + velocity * distance
        
        last_position = position
    
    velocity_average = velocity_average_accumulator / distance_accumulator

    # Intermediate results
    print "Average velocity", velocity_average

    # Build data for graphing
    half_velocity_profile = [(position, velocity) for position, velocity in velocity_profile if position <= (wall_to_wall_distance / 2)]
    pl_wall_functions[pl] = [(yplus(position), uplus(velocity)) for position, velocity in half_velocity_profile]
    pl_scaled_velocity_profiles[pl] = [(velocity / velocity_average, 2 * position / d_h) for position, velocity in half_velocity_profile]


# Find the theoretical wall function
all_yplus = unions([set([yplus for yplus, _ in wall_function]) for wall_function in pl_wall_functions.values()])
theory_wall_function = [(yplus, uplus_theory(yplus)) for yplus in all_yplus]


# OK, lets go:
# http://matplotlib.sourceforge.net/users/pyplot_tutorial.html
import matplotlib.pyplot as plt

# y+ vs u+:
plt.clf()

plt.xlabel("y+")
plt.ylabel("u+")

for label, wall_function in [("theory", theory_wall_function)] + [("pl=" + str(pl), wall_function) for pl, wall_function in pl_wall_functions.items()]:
    yplus, uplus = zip(*wall_function)
    plt.semilogx(yplus, uplus, label=label)

plt.legend()
plt.savefig("yplus-vs-uplus")

# scaled velocity:
plt.clf()

plt.xlabel("u/u_avg")
plt.ylabel("2y/D_H")

for pl, scaled_velocity_profile in pl_scaled_velocity_profiles.items():
    scaled_velocity, scaled_y = zip(*scaled_velocity_profile)
    plt.plot(scaled_velocity, scaled_y, label="pl=" + str(pl))

plt.legend()
plt.savefig("scaled-velocity-profile")
