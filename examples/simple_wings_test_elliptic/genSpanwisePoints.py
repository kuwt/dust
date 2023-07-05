import numpy as np

####### input section ##########
rawinput = """
-0.383525000000000E+001 -0.350175000000000E+001 -0.316825000000000E+001 -0.283475000000000E+001 -0.250125000000000E+001 -0.216775000000000E+001 -0.183425000000000E+001 -0.150075000000000E+001 -0.116725000000000E+001 -0.833750000000000E+000 -0.500250000000000E+000 -0.166750000000000E+000  0.166750000000000E+000  0.500250000000000E+000  0.833750000000000E+000  0.116725000000000E+001  0.150075000000000E+001  0.183425000000000E+001  0.216775000000000E+001  0.250125000000000E+001  0.283475000000000E+001  0.316825000000000E+001  0.350175000000000E+001  0.383525000000000E+001
"""

root_chord_length = 1
##############################

# get x positions (1/4 chord length)
x = root_chord_length * 1.0/4.0

# get y positions
span_y_position = rawinput.split(" ")
span_y_position_temp = []
for elem in span_y_position:
    if elem != '':
        span_y_position_temp.append(float(elem))
        
#print(span_y_position_temp)
span_y_position = span_y_position_temp
print("There are {} elements".format(len(span_y_position)))

# get z positions 
z = 0

for elem in span_y_position:
    print("point = (/ {}, {}, {} /)".format(x,elem,z))

