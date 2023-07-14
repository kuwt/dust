import numpy as np

####### input section ##########

center_y = """
 -0.383525000000000E+001 -0.350175000000000E+001 -0.316825000000000E+001 -0.283475000000000E+001 -0.250125000000000E+001 -0.216775000000000E+001 -0.183425000000000E+001 -0.150075000000000E+001 -0.116725000000000E+001 -0.833750000000000E+000 -0.500250000000000E+000 -0.166750000000000E+000  0.166750000000000E+000  0.500250000000000E+000  0.833750000000000E+000  0.116725000000000E+001  0.150075000000000E+001  0.183425000000000E+001  0.216775000000000E+001  0.250125000000000E+001  0.283475000000000E+001  0.316825000000000E+001  0.350175000000000E+001  0.383525000000000E+001
"""

chord_raw_input = """
  0.775000000000000E+000  0.110600000000000E+001  0.129850000000000E+001  0.149100000000000E+001  0.161150000000000E+001  0.173200000000000E+001  0.180900000000000E+001  0.188600000000000E+001  0.192900000000000E+001  0.197200000000000E+001  0.198600000000000E+001  0.200000000000000E+001  0.200000000000000E+001  0.198600000000000E+001  0.197200000000000E+001  0.192900000000000E+001  0.188600000000000E+001  0.180900000000000E+001  0.173200000000000E+001  0.161150000000000E+001  0.149100000000000E+001  0.129850000000000E+001  0.110600000000000E+001  0.775000000000000E+000
"""

angle = 10.0 * np.pi/180.0
ratioOfChord = 0.5 # which position along the chord: e.g. quarter chord, 3 quarter chord 
##############################

# get y positions
y_position = center_y.split()

y_position_temp = []
for elem in y_position:
    y_position_temp.append(float(elem))
        
y_position_temp = y_position_temp
print("There are {} elements".format(len(y_position)))

# get x and z positions
chords = chord_raw_input.split()
chords_temp = []
for elem in chords:
    chords_temp.append(float(elem))
chords = chords_temp
print("There are {} elements".format(len(chords)))

x_position = []
z_position = []
for elem in chords:
    x = elem * np.cos(-angle) * ratioOfChord
    x_position.append(x)
    z = elem * np.sin(-angle) * ratioOfChord
    z_position.append(z)


for i, x in enumerate(x_position):
    print("point = (/ {}, {}, {} /)".format(x,y_position[i],z_position[i]))

