import numpy as np
import os
import matplotlib.pyplot as plt
############ input ##################
pointprobes = 24
timestep = 0


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
with open('Postprocessing/post_prb01.dat','r') as file:
    listOfLines = file.readlines()
    
    ########## read the position information  #########
    line = listOfLines[2]
    content = readlineIntoFloats(line)
    y_position = content
    
    ########## read the information at the selected time step #########
    timestep = timestep + 7 # offset because of the organization of the file post_prbxx.dat
    line = listOfLines[timestep]
    content = readlineIntoFloats(line)
    
    # obtain the time 
    time = content.pop(0)
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
    
    plt.plot(y_position,ux_table)
    plt.plot(y_position,uy_table)
    plt.plot(y_position,uz_table)
    plt.legend(["ux","uy", "uz"])
    plt.show()
    
    plt.plot(y_position,uz_table)
    plt.legend(["uz"])
    plt.show()
