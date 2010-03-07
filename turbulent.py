#!/usr/bin/env python

import csv
import sys
import os.path
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
    for name, value in ks:
        if type(value) == type(0.0):
            value = "%.6g" % value
        
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
printconstants("Model constants:", [("width", width), ("height", height), ("d_h", d_h), ("density", density), ("mu", mu)])


# Pull in data
sim_wall_functions = {}
sim_scaled_velocity_profiles = {}
for file in sys.argv[1:]:
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
    
    # "Position [0.0, 1.0, 0.0] (m)-point 1 (m)","Velocity: Magnitude-point 1 (m/s)"
    # 0.0,35.15489668712291
    # 2.499999993688107E-7,35.15489668712291
    reject_count = 0
    raw_velocity_positions = []
    
    # Skip header row
    reader = csv.reader(open(file))
    header = reader.next()
    is_wall_distance = "Wall Distance" in header[0]

    if is_wall_distance:
        print "WARNING: using wall-distance position information, results may be inaccurate"
    
    # Input body of data
    all_rows = list(reader)
    last_position = None
    for row in all_rows:
        # Verify record is correct length
        if len(row) != 2:
            reject_count = reject_count + 1
            continue
        
        # Verify we can parse as numbers
        try:
            position, velocity = float(row[0]), float(row[1])
            
            # If reading a file containing wall distance, we only want half of the profile
            # or we get weirdy discontinituites in the output
            if is_wall_distance and position < last_position:
                # Diagnostics for the user to help them figure out disontinutities
                broke_at = len(raw_velocity_positions)
                expected_break_at = len(all_rows) / 2
                row_different = expected_break_at - broke_at
                
                print "Breaking wall-distance series at", position, "<", last_position
                print "Broke after", broke_at, "rows, expecting", expected_break_at, "rows"
                if row_different > 1:
                    print "CRITICAL WARNING: the number of rows got is less than our expectation by", row_different, "so the series is probably truncated"
                elif row_different < -1:
                    print "CRITICAL WARNING: the number of rows got exceeds our expectation by", abs(row_different), "so the series may contain false readings"
                
                break
            else:
                last_position = position
            
            raw_velocity_positions.append((position, velocity))
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
    
    velocity_profile = sorted([(position, velocity) for velocity, position in velocity_positions.items()])
    
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
    positions, velocities = zip(*raw_velocity_positions)
    printconstants("Simulation summary:", [("average velocity", velocity_average),
                                           ("minimum velocity", min(velocities)),
                                           ("maximum velocity", max(velocities)),
                                           ("minimum position", min(positions)),
                                           ("maximum position", max(positions)),
                                           ("minimum non-zero y+", yplus(min([position for position in positions if position > 0.0])))])
    
    # Trim data down to size (especially important if the velocity profile was built from position-indexed information)
    half_velocity_profile = [(position, velocity) for position, velocity in velocity_profile if position <= (wall_to_wall_distance / 2) and position > 0.0]
    
    # Build data for graphing
    sim = "pl=" + str(pl) + ",sf=" + sf
    sim_wall_functions[sim] = [(yplus(position), uplus(velocity)) for position, velocity in half_velocity_profile]
    sim_scaled_velocity_profiles[sim] = [(velocity / velocity_average, 2 * position / d_h) for position, velocity in half_velocity_profile]


# Find the theoretical wall function
all_yplus = unions([set([yplus for yplus, _ in wall_function]) for wall_function in sim_wall_functions.values()])
theory_wall_function = sorted([(yplus, uplus_theory(yplus)) for yplus in all_yplus])


# Compute best fit to the theory
# import bisect
# 
# print ""
# print "# Distance from theory:"
# best_match, best_match_distance = None, sys.maxint
# for sim, wall_function in sim_wall_functions.items():
#     wall_function_ypluses = [yplus for yplus, _ in wall_function]
#     
#     def wall_function_at(yplus):
#         # Linearly interpolate to approximate the wall function at this point
#         i = bisect.bisect(wall_function_ypluses, yplus)
#         yplus1, uplus1 = wall_function[i]
#         yplus2, uplus2 = wall_function[i - 1]
#         return uplus1 + (yplus - yplus1) * (uplus2 - uplus1) / (yplus2 - yplus1)
#     
#     sum_sqs_distance = sum([pow(wall_function_at(yplus) - uplus_theory, 2) for yplus, uplus_theory in theory_wall_function])
#     print "  ", sim + ":", sum_sqs_distance
#     if sum_sqs_distance < best_match_distance:
#         best_match, best_match_distance = sim, sum_sqs_distance
# 
# print "# Best match:", best_match


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

#plt.axis([0.1, 110, 0, 25])
plt.legend(loc="upper left")
plt.savefig("yplus-vs-uplus")

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
