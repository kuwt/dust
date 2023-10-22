import numpy as np
import os
import matplotlib.pyplot as plt
############ input ##################
timestep = 10

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
with open('inflow.dat','r') as file:
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
    
    ########## read the information at the selected time step #########
    line = listOfLines[timestep]
    content = readlineIntoFloats(line)
    
    # obtain the time 
    time = content.pop(0)
    # obtain the inflow 
    ux_table = []
    uy_table = []
    uz_table = []
    while len(content) != 0:
        ux = content.pop(0)
        uy = content.pop(0)
        uz = content.pop(0)
        ux_table.append(ux)
        uy_table.append(uy)
        uz_table.append(uz)
    
    #########processing ##########################
    xdiff = np.diff(x_table)
    ydiff = np.diff(y_table)
    zdiff = np.diff(z_table)

    posdiff = np.sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff)
    pos = np.add.accumulate(posdiff)
    pos = np.insert(pos, 0, 0)
   
    u_mag = np.sqrt(np.array(ux_table)* np.array(ux_table)+ np.array(uy_table)*np.array(uy_table) + np.array(uz_table)*np.array(uz_table))
    ####### ploting #####################
    plt.figure(0)
    plt.plot(pos,ux_table)
    plt.plot(pos,uy_table)
    plt.plot(pos,uz_table)
    plt.legend(["inflow_x","inflow_y", "inflow_z"])
    plt.savefig("all.png")
    plt.show()

    plt.figure(1)
    plt.plot(pos,u_mag)
    plt.legend(["inflow_mag"])
    plt.savefig("inflow_mag.png")
    plt.show()