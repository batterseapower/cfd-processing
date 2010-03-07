#!/usr/bin/env python

from glob import glob
import sys
import os.path
import re
from math import sqrt, log

from common import *


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
printconstants("Model constants:", [("width", width), ("height", height), ("d_h", d_h), ("density", density), ("mu", mu)])


# Pull in data
sim_wall_functions = {}
sim_scaled_velocity_profiles = {}
for file in concat(map(glob, sys.argv[1:])):
    print ""
    print "#", file
    
    # pl29sf1.13base0.1distance.csv
    # r"pl(\d+)sf([\d\.]+)base([\d\.]+)distance"
    file_base = os.path.splitext(file)[0]
    
    # Parse our prism layer count
    m = re.search(r"pl(\d+)", file_base)
    if m is None:
        print "Could not determine prism layer count: add plN to the file name"
        sys.exit(1)
    else:
        pl = int(m.group(1))
    
    # Parse out stretch factor
    m = re.search(r"sf([\d\.]+)", file_base)
    sf = m is None and "unknown" or m.group(1)
    
    # Parse out wall shear stress
    m = re.search(r"tw([\d\.]+)", file_base)
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
    printconstants("Simulation constants:", [("pl", pl), ("wall shear stress", wall_shear_stress), ("ustar", ustar), ("d_v", d_v)])
    
    # Extract data from CSV
    velocity_average, velocity_profile = extractdata(file, extra_stats = lambda positions: [("minimum non-zero y+", yplus(min([position for position in positions if position > 0.0])))])
    
    # Trim data down to size (especially important if the velocity profile was built from position-indexed information)
    half_velocity_profile = [(position, velocity) for position, velocity in velocity_profile if 0.0 < position <= (wall_to_wall_distance / 2)]
    
    # Build data for graphing
    sim = "pl=" + str(pl) + ",sf=" + sf
    sim_wall_functions[sim] = [(yplus(position), uplus(velocity)) for position, velocity in half_velocity_profile]
    sim_scaled_velocity_profiles[sim] = [(velocity / velocity_average, 2 * position / d_h) for position, velocity in half_velocity_profile]


# Find the theoretical wall function
all_yplus = unions([set([yplus for yplus, _ in wall_function]) for wall_function in sim_wall_functions.values()])
theory_wall_function = sorted([(yplus, uplus_theory(yplus)) for yplus in all_yplus])


# Compute best fit to the theory
print ""
print "# Simulation average squared error"
for sim, wall_function in sim_wall_functions.items():
    print "  ", sim, "=", average_error(wall_function, uplus_theory)


# OK, lets go:
# http://matplotlib.sourceforge.net/users/pyplot_tutorial.html
import matplotlib.pyplot as plt

# y+ vs u+:
plt.clf()
plt.figure(figsize=(10, 8), dpi=80)

plt.xlabel("y+")
plt.ylabel("u+")

for sim, wall_function in [("theory", theory_wall_function)] + sim_wall_functions.items():
    yplus, uplus = zip(*wall_function)
    plt.semilogx(yplus, uplus, label=sim)

plt.legend(loc="upper left")
plt.savefig("yplus-vs-uplus")

plt.axis([0.1, 110, 0, 25])
plt.savefig("yplus-vs-uplus-experimental")

# scaled velocity:
plt.clf()
plt.figure(figsize=(10, 8), dpi=80)

plt.xlabel("u/u_avg")
plt.ylabel("2y/D_H")

for sim, scaled_velocity_profile in sim_scaled_velocity_profiles.items():
    scaled_velocity, scaled_y = zip(*scaled_velocity_profile)
    plt.plot(scaled_velocity, scaled_y, label=sim)

plt.legend(loc="upper left")
plt.savefig("scaled-velocity-profile")


# Build the superimposed data image
from PIL import *
from PIL import Image

lift_point = lambda f: lambda p1, p2: (f(p1[0], p2[0]), f(p1[1], p2[1]))
minus_point = lift_point(lambda x1, x2: x1 - x2)
divide_point = lift_point(lambda x1, x2: x1 / float(x2))
multiply_point = lift_point(lambda x1, x2: x1 * x2)

# Discover the boxes in both images that we want to superimpose
experimental = Image.open("experimental.png")
(experimental_width, experimental_height) = experimental.size
experimental_bottom_left = (88, 108) # x = 0.1, y = 0
experimental_top_right = (experimental_width - 82, experimental_height - 16)  # x = 100, y = 25

mine = Image.open("yplus-vs-uplus-experimental.png")
(mine_width, mine_height) = mine.size
mine_bottom_left = (127, 81)
mine_top_right = (mine_width - 111, mine_height - 81)

# Begin by scaling the experimental data correctly so the axis-enclosed boxes look the same
scale = divide_point(minus_point(mine_top_right, mine_bottom_left), minus_point(experimental_top_right, experimental_bottom_left))
experimental = experimental.resize(multiply_point(scale, experimental.size))
#experimental.save("yplus-vs-uplus-experimental-SCALED.png", "PNG")

# Now work out the transform to shift the experimental data so the origins match
#print mine_bottom_left, "-", multiply_point(scale, experimental_bottom_left)
transform = minus_point(mine_bottom_left, multiply_point(scale, experimental_bottom_left))
#print transform
experimental = experimental.offset(int(round(transform[0])), -int(round(transform[1])) + 8) # FIXME: fudge factor!

blended = Image.blend(mine, experimental.crop((0, 0, mine.size[0], mine.size[1])), 0.3)
blended.save("yplus-vs-uplus-experimental-combined.png", "PNG")