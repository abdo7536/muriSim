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
### nm = The string name whose file name group needs to be created (including the path to the directory) [type - string]
### form = the file format (input exactly as it appears in the filenames) [type - string]

### Output(s):
### onlyfiles = A sorted list of file names of the user defined file extension [type - list of string]

################# Function Definition
def plname(nm,form):
################# Returns list of path names that match the string provided
    onlyfiles = sorted(glob.glob( nm + '*.' + form))
################# Return the function's output variables
    return (onlyfiles)

##########################################################################################################################################################################
################# function to generate X, Y and Z index points associated to sampling strategy chosen by the user
### Input(s): All inputs are mandatory
### trajTyp = Synthetic Observation trajectory type; 1 = balloon-like vertical trajectory (descending); 2 = Helical trajectory (descending); 3 = Other trajectory/sampling strategy (descending)
### Xref = User defined reference points in X
### Yref = User defined reference points in Y
### Nx = Number of wavenumbers/grid points in X
### Ny = Number of wavenumbers/grid points in Y
### Nz = Number of wavenumbers/grid points in Z

### Output(s):
### x_sampl = X indices for all points on the synthetic trajectory
### y_sampl = Y indices for all points on the synthetic trajectory
### z_sampl = Z indices for all points on the synthetic trajectory

################# Function Definition
def trajInd(trajTyp,Xref,Yref,Nx,Ny,Nz):
################# Check what trajectory is defined 
    if trajTyp[0] == 1:         # Balloon-like descent
        x_sampl = np.ones(Nz)*Xref
        y_sampl = np.ones(Nz)*Yref
        z_sampl = np.linspace(Nz,1,Nz)
    elif trajTyp[0] == 2:       # DH helical descent (NOT WORKING CURRENTLY)
        x_sampl = np.ones(Nz)
        y_sampl = np.ones(Nz)
        z_sampl = np.ones(Nz)
    elif trajTyp[0] == 3:       # Other trajectories (NOT WORKING CURRENTLY)
        x_sampl = np.ones(Nz)
        y_sampl = np.ones(Nz)
        z_sampl = np.ones(Nz)
################# Return the function's output variables
    return(x_sampl, y_sampl, z_sampl)

##########################################################################################################################################################################
################# function to get bytes of data corresponding to user specified data point in a 3D DNS data field
### Input(s): All inputs are mandatory
### TrX = The X index points on the trajectory
### TrY = The Y index points on the trajectory
### TrZ = The Z index points on the trajectory
### Nx = Number of wavenumbers/grid points in X
### Ny = Number of wavenumbers/grid points in Y
### Nz = Number of wavenumbers/grid points in Z
### fl = The structrue which holds all the data in a file

### Output(s):
### dat = The output data point (signed float)

def datGet(TrX,TrY,TrZ,Nx,Ny,Nz,fl):
    fldIdent = 1
    pts = Nx*Ny*Nz                  # The number of data points in each 3D dataset
    pt_jmp = ((TrZ-1)*Nx*Ny)+((TrY-1)*Nx)+(TrX-1)+((fldIdent-1)*pts)
    byt_jmp = int(pt_jmp*4)
    fl.seek(byt_jmp)
    tl_u = fl.tell()
    b=fl.read(4)
    dat = np.array(list(st.unpack('f', b)))
################# return the output data point
    return(dat)
