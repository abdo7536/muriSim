################# This python script contains a list of custom scripted python functions 
### that simplify user interactions with paraview.simple module to help create 
### effective pvbatch scripts for automated visualizations using paraview

################# All the dependencies to run inbuilt or custom functions that are called within the listed functions below need to be imported beforehand
### This section imports all the required function

try:
    import numpy as np
    print("numpy module is loaded...")
except ImportError:
    print("numpy module not installed...install numpy and try again")
try:
    import math
    print("math module is loaded...")
except ImportError:
    print("math module not installed...install math and try again")
try:
    import struct as st
    print("struct is installed...")
except ImportError:
    print("struct is not installed")  
try:
    import time                             
    print("time module is loaded...")
except ImportError:
    print("time module not installed...install numpy and try again")
try:
    import glob
    print("glob module is loaded...")
except ImportError:
    print("glob module not installed...install glob and try again")
try:
    import random
except ImportError:
    print("import module not installed...install random and try again")

##########################################################################################################################################################################
################# function to generate a list of file name groups
### Input(s): All inputs are mandatory
### str = The string name whose file name group needs to be created (including the path to the directory) [type - string]
### form = the file format (input exactly as it appears in the filenames) [type - string]

### Output(s):
### onlyfiles = A sorted list of file names of the user defined file extension [type - list of string]

################# Function Definition
def plname(str,form):
################# Returns list of path names that match the string provided
    onlyfiles = sorted(glob.glob( str + '*.' + form))
################# Return the function's output variables
    return (onlyfiles)

##########################################################################################################################################################################
################# function to generate X, Y and Z index points associated to sampling strategy chosen by the user
### Input(s): All inputs are mandatory
### trajTyp = Synthetic Observation trajectory type; 1 = balloon-like vertical trajectory (descending); 2 = Helical trajectory (descending); 3 = Other trajectory/sampling strategy (descending)
 
### Output(s):
### x_sampl = X indices for all points on the synthetic trajectory
### y_sampl = Y indices for all points on the synthetic trajectory
### z_sampl = Z indices for all points on the synthetic trajectory

################# Function Definition
def trajInd(trajTyp,numTraj,Xscal,Yscal,Zscal):
#################
    if trajTyp[0] == 1:         # Balloon-like descent
        ### Calculate the number of points to sample on each trajectory
        numSmpls = Zscal/trajTyp[1]
        ptsTraj = math.floor(Nz/numSmpls)
        
    elif trajTyp[0] == 2:       # DH helical descent

    elif trajTyp[0] == 3:       # Other trajectories

#################    

################# Return the function's output variables
    return(x_sampl, y_sampl, z_sampl)