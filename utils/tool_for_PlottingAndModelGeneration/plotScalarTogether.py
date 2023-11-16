import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import argparse
#import pdb; pdb.set_trace()
matplotlib.rcParams.update({'font.size': 32})
############ input ##################
timestep = 30
legend = "inflow"
files = [\
"/home/k/Share/code/dust/examples/simple_wings_test_elliptic/inflow.dat",\
#"/home/k/Share/code/dust/examples/simple_wings_test_elliptic/inflow_mid.dat",\
#"/home/k/Share/code/dust/examples/simple_wings_test_elliptic/inflow_cor.dat",\
]
#legends = ["all wake", "end vortex + particles","particles"]
legends = ["inflow"]
###################################

############# utility ###############
def readlineIntoFloats(line):
    content = line.split()
    content_temp = []
    for elem in content:
        content_temp.append(float(elem))
    content = content_temp
    return content
##############################
#print(os.listdir())
with open(files[0],'r') as file:
    listOfLines = file.readlines()
    
    ########## read the position information  #########
    line = listOfLines[0]
    content = readlineIntoFloats(line)
    x_table = []
    y_table = []
    z_table = []

    while len(content) != 0:
        x = content.pop(0)
        y = content.pop(0)
        z = content.pop(0)
        x_table.append(x)
        y_table.append(y)
        z_table.append(z)

    xdiff = np.diff(x_table)
    ydiff = np.diff(y_table)
    zdiff = np.diff(z_table)

    posdiff = np.sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff)
    pos = np.add.accumulate(posdiff)
    pos = np.insert(pos, 0, 0)

u_tables = []
for file in files:
    with open(file,'r') as file:
        listOfLines = file.readlines()
        ########## read the information at the selected time step #########
        line = listOfLines[timestep]
        content = readlineIntoFloats(line)
        
        # obtain the time 
        time = content.pop(0)
        # obtain the inflow 
        u_table = []
        while len(content) != 0:
            u = content.pop(0)
            u_table.append(u)
        u_tables.append(u_table)

   
####### ploting #####################
#plt.figure(0)
#plt.plot(pos,u_table)
#plt.legend([legend])
#plt.show()

########### plot #####################
plt.figure(1)
for i,file in enumerate(files):
    plt.plot(pos,u_tables[i])
plt.legend(legends)
plt.xlabel("span")
plt.show()
