import csv
import sys

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

def infinite_sum(term):
    n = 1
    result = 0.0
    contribution = sys.maxint
    while abs(contribution) > 0.000001 and n < 40:
        contribution = term(n)
        #print n, "=", "%.6g" % contribution
        result += contribution
        n += 2 # Only want odd terms

    return result

def average_error(observed, model):
    return sum([pow(model(x) - y, 2) for x, y in observed]) / len(observed)

def extractdata(file, extra_stats=lambda _: []):
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

    if is_wall_distance:
        # Build the velocity profile
        min_velocity_positions = {}
        max_velocity_positions = {}
        for position, velocity in raw_velocity_positions:
            # Record the positions that experienced a given velocity
            min_velocity_positions[velocity] = min(min_velocity_positions.get(velocity, sys.maxint), position)
            max_velocity_positions[velocity] = max(min_velocity_positions.get(velocity, -sys.maxint), position)
    
        # Use the average position
        velocity_profile = sorted([((min_velocity_positions[velocity] + max_velocity_positions[velocity]) / 2, velocity) for velocity in min_velocity_positions.keys()])
    else:
        velocity_profile = []
        
        # Assume the data for adjacent positions are adjacent in a position-indexed file
        velocity_in_run = None
        first_position_in_run = None
        for position, velocity in raw_velocity_positions:
            if velocity != velocity_in_run:
                if velocity_in_run is not None:
                    velocity_profile.append(((position + first_position_in_run) / 2 , velocity_in_run))
                
                first_position_in_run = position
                velocity_in_run = velocity
        
        if velocity_in_run is not None:
            velocity_profile.append(((raw_velocity_positions[-1][0] + first_position_in_run) / 2 , velocity_in_run))
        
        velocity_profile.sort()
    
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
                                           ("maximum position", max(positions))] + extra_stats(positions))
    
    return velocity_average, velocity_profile
