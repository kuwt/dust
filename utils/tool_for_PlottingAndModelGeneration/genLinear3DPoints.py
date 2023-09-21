import numpy as np

####### input section ##########
startpt = np.array([0,4,0]) 
endpt = np.array([1.09,4,-0.1])
ptNum = 24
##############################
for i in range(ptNum):
    pt = startpt + (endpt - startpt) * i/ptNum
    print("point = (/ {}, {}, {} /)".format(pt[0],pt[1],pt[2]))

    

