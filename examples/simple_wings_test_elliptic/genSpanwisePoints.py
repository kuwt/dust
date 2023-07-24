import numpy as np

####### input section ##########

center_y = """
 -0.764750000000000E+000 -0.698250000000000E+000 -0.631750000000000E+000 -0.565250000000000E+000 -0.498750000000000E+000 -0.432250000000000E+000 -0.365750000000000E+000 -0.299250000000000E+000 -0.232750000000000E+000 -0.166250000000000E+000 -0.997500000000000E-001 -0.332500000000000E-001  0.332500000000000E-001  0.997500000000000E-001  0.166250000000000E+000  0.232750000000000E+000  0.299250000000000E+000  0.365750000000000E+000  0.432250000000000E+000  0.498750000000000E+000  0.565250000000000E+000  0.631750000000000E+000  0.698250000000000E+000  0.764750000000000E+000
"""

chord_raw_input = """
  0.153000000000000E+000  0.188000000000000E+000  0.220500000000000E+000  0.253000000000000E+000  0.273500000000000E+000  0.294000000000000E+000  0.307500000000000E+000  0.321000000000000E+000  0.328000000000000E+000  0.335000000000000E+000  0.337500000000000E+000  0.340000000000000E+000  0.340000000000000E+000  0.337500000000000E+000  0.335000000000000E+000  0.328000000000000E+000  0.321000000000000E+000  0.307500000000000E+000  0.294000000000000E+000  0.273500000000000E+000  0.253000000000000E+000  0.220500000000000E+000  0.188000000000000E+000  0.153000000000000E+000
"""

angle = 5.0 * np.pi/180.0
ratioOfChord = 0 # which position along the chord: e.g. quarter chord, 3 quarter chord 
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

