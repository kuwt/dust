import numpy as np
import os
import matplotlib.pyplot as plt
import argparse
############ input ##################
linerangel = 1
linerangeu = 1500
rmin = 0
rmax = 39

###################################
parser = argparse.ArgumentParser()
parser.add_argument("-i",dest='filename',required=True)
args = parser.parse_args()
print(args.filename)
filename = args.filename


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
with open(filename,'r') as file:
    listOfLines = file.readlines()
    
    timelist = np.zeros((linerangeu-linerangel)+1)

    ########## read the position information  #########
    line = listOfLines[0]
    content = readlineIntoFloats(line)
    x_table = []
    y_table = []
    z_table = []
    cur_r = 0
    while len(content) != 0:
        x = content.pop(0)
        y = content.pop(0)
        z = content.pop(0)
        if cur_r >= rmin and cur_r <= rmax:
            x_table.append(x)
            y_table.append(y)
            z_table.append(z)
        cur_r = cur_r + 1
    
    xdiff = np.diff(x_table)
    ydiff = np.diff(y_table)
    zdiff = np.diff(z_table)
    posdiff = np.sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff)
    pos = np.add.accumulate(posdiff)
    pos = np.insert(pos, 0, 0)
    rlist = np.array(pos)
    rlist = rlist + 1.08 # rotor starts from 1.08
    rlist = rlist/max(rlist)
    print("rlist = {}".format(rlist))
    print("rlist len = {}".format(len(rlist)))

    ########## read the information at the selected time step #########
    Zlist  = np.zeros((len(timelist),len(rlist)))

    cur_time_step = 0
    for linenum,line in enumerate(listOfLines):
        if linenum >= linerangel and linenum <= linerangeu:
            content = readlineIntoFloats(line)
            # obtain the time 
            time = content.pop(0)
            timelist[cur_time_step] = cur_time_step
            
            cur_r_step = 0
            # obtain the inflow 
            while len(content) != 0:
                value = content.pop(0)
                Zlist[cur_time_step][cur_r_step] = value
                cur_r_step = cur_r_step + 1
                if cur_r_step >rmax:
                    break
            cur_time_step = cur_time_step + 1

    print("Zlist = {}".format(Zlist))
    print("Zlist shape= {}".format(Zlist.shape))
    ########### plot #####################
    rmesh, timemesh = np.meshgrid(rlist,timelist)
    fig, ax = plt.subplots(dpi=120)
    cs = ax.contourf(timemesh, rmesh, Zlist, 300)
    cbar = fig.colorbar(cs)
    plt.show()
