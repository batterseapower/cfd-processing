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

def printconstants(title, ks):
    print title
    for name, value in ks.items():
        print "  ", name, "=", value


# NB: fixed for the model
wall_to_wall_distance = 0.05
width = 0.05
height = 0.05
d_h = (width * height) / (2 * (width + height)) # Hydralic diameter

density = 1.18415
mu = 0.00001855


def uplus_theory(yplus):
    # Computation Methods for Fluid Dynamics, p298
    if yplus <= 5:
        return yplus
    else:
        k = 0.41
        b = 5
        return (1 / k) * log(yplus) + b


# Print some information about the model
printconstants("Model constants:", { "width" : width, "height" : height, "d_h" : d_h, "density" : density, "mu" : mu })


# Pull in data
pl_wall_functions = {}
pl_scaled_velocity_profiles = {}
for file in sys.argv[1:]:
    print ""
    print "#", file
    
    # pl29sf1.13base0.1distance.csv
    # r"pl(\d+)sf([\d\.]+)base([\d\.]+)distance"
    
    # Parse our prism layer count
    m = re.search(r"pl(\d+)", file)
    if m is None:
        print "Could not determine prism layer count: add plN to the file name"
        sys.exit(1)
    else:
        pl = int(m.group(1))
    
    # Parse out wall shear stress
    m = re.search(r"tw([\d\.]+)", file)
    if m is None:
        print "Count not determine wall shear stress: add twN.M to the file name"
        sys.exit(1)
    else:
        #wall_shear_stress = 111.9
        wall_shear_stress = float(m.group(1))
    
    # Setup physical constants for this model
    ustar = sqrt(wall_shear_stress / density)
    d_v = mu / ustar
    
    def yplus(wall_distance):
        return wall_distance / d_v

    def uplus(velocity):
        return velocity / ustar

    # Print some information about the simulation
    printconstants("Simulation constants:", { "pl" : pl, "wall shear stress" : wall_shear_stress, "ustar" : ustar, "d_v" : d_v })
    
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
    printconstants("Simulation summary:", { "average velocity" : velocity_average, "maximum velocity" : max([velocity for position, velocity in raw_velocity_positions]) })

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

plt.legend(loc="upper left")
plt.savefig("yplus-vs-uplus")

# scaled velocity:
plt.clf()

plt.xlabel("u/u_avg")
plt.ylabel("2y/D_H")

for pl, scaled_velocity_profile in pl_scaled_velocity_profiles.items():
    scaled_velocity, scaled_y = zip(*scaled_velocity_profile)
    plt.plot(scaled_velocity, scaled_y, label="pl=" + str(pl))

plt.legend(loc="upper left")
plt.savefig("scaled-velocity-profile")
