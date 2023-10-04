import os


input_filename = "coupling_nodes.in"
output_filename = "coupling_nodes3.in"
mod_x = 0.003
mod_y = 0
mod_z = 0

################################
f = open(input_filename)
readcontent = f.read().split()
print(readcontent)
f.close()

new_x = []
new_y = []
new_z = []
for i,elem in enumerate(readcontent):
    if (i % 3 == 0):
        new_x.append( float(elem) + mod_x)
    if (i % 3 == 1):
        new_y.append( float(elem) + mod_y)
    if (i % 3 == 2):
        new_z.append( float(elem) + mod_z)


w = open(output_filename,"w")
for i, elem in enumerate(new_x):
    w.write("{} {} {}\n".format(new_x[i],new_y[i],new_z[i]))
w.close()
