import numpy
import sys
import vtk
import glob
import re
import math
from vtk.util import numpy_support as VN

input_filename = './Postpro/post_ff01_0005.vtr'
print ('Loading', input_filename)

######## read file ############
#reader = vtk.vtkXMLUnstructuredGridReader()
reader = vtk.vtkXMLRectilinearGridReader()
reader.SetFileName(input_filename)
reader.Update()
data = reader.GetOutput()

######## process  ############
gradientFilter = vtk.vtkGradientFilter()
gradientFilter.SetInputData(data)
gradientFilter.SetInputArrayToProcess(0,0,0,0,'velocity')
gradientFilter.SetResultArrayName('gradu')
gradientFilter.ComputeVorticityOn()  #Gets us the curl
gradientFilter.Update()
data_grad = gradientFilter.GetOutput()


#If you want to process the gradient tensor:
#wss_grad_vector = VN.vtk_to_numpy(data_grad .GetPointData().GetArray('gradu'))
#wss_grad_vector[i,0] #This is \partial u / \partial x
#wss_grad_vector[i,1] #This is \partial u / \partial y
#wss_grad_vector[i,2] #This is \partial u / \partial z
#wss_grad_vector[i,3] #This is \partial v / \partial x
#wss_grad_vector[i,4] #This is \partial v / \partial y
#wss_grad_vector[i,5] #This is \partial v / \partial z
#wss_grad_vector[i,6] #This is \partial w / \partial x
#wss_grad_vector[i,7] #This is \partial w / \partial y
#wss_grad_vector[i,8] #This is \partial w / \partial z

######## write results  ############
myoutput = vtk.vtkXMLDataSetWriter() 
myoutput.SetInputData(data_grad)
output_filename ='a.vtr'
myoutput.SetFileName(output_filename)
myoutput.Write()

print ('Done!')
