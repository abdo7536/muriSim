################# Script to sample data on a user defined trajectory from a subvol***.dat file
################# Author: Dr. Abhiram Doddi
################# Date Created: 1st April 2022
################# Last updated: 1st April 2022

##########################################################################################################################################################################
################# General Inputs: All inputs are mandatory
### Xl = Domain X dimension
Xl = 3
### Yl = Domain Y dimension
Yl = 2
### Zl = Domain Z dimension
Zl = 1
### Nx = Number of points in Domain X dimension
Nx = 2880
### Ny = Number of points in Domain Y dimension
Ny = 1920
### Nz = Number of points in Domain Z dimension
Nz = 960
### numTraj = Number of trajectories to sample
numTraj = 10
### locFlg = Flag to control the (X,Y) location for sampling; 1 = choose random (x,y) location for each trajectory; 2 = Specific (x,y) location (DEFAULT set to X/2, Y/2 location)
locFlg = 1
### datVar = The data variable to be stored; NOTE: This will only be used to name the output filename
datVar = 'U'

################# Data Inputs: All inputs are mandatory
### NOTE: These inputs are required to scale the scale-normalized DNS datasets  
### epsMeas = TKE dissipation rate [m^3s^-2]
epsMeas = 4.0e-3
### nu = Kinematic viscosity [m^2s^-1]
nu = 4.0e-4
### resMet = DNS resolution metric ( = grid spacing/ kolmogorov length scale)
resMet = 1.4413
### lxGW = GW X component wavelength (unitless/normalized) 
lxGw = 3
### lzGW = GW Z component wavelength (unitless/normalized) 
lzGw = 1
### aGW = GW amplitude relative to overturning wave amplitude (unitless)
a = 0.9

################# Function Specific Inputs: All inputs are mandatory
########## function name: 'lodatsin'
### flNm = The string name whose file name group needs to be created (including the path to the directory) [type - string]
flNm = '/Users/script_away/Projects/Documents/MURI_modeling/GWBData/subvol' + datVar + '*'
### form = the file format (input exactly as it appears in the filenames) [type - string]
form = 'dat'
### flNm = Input file name identifying character string (including the directory location)

################# Function Specific Inputs: All inputs are mandatory
########## function name: 'trajInd'
### trajTyp = Synthetic Observation trajectory type; 1 = balloon-like vertical trajectory (descending); 2 = Helical trajectory (descending); 3 = Other trajectory/sampling strategy (descending)
trajTyp = [1]                   # Example input for Balloon-like Trajectory

##########################################################################################################################################################################
################# Call custom function module developed to execute this program
### NOTE: This module inturn loads the python inbuilt modules as needed
try:
    from userFns import *
    print("Program required custom and python module is loaded...")
except ImportError:
    print("Custom library not installed.... Install and try again")

##########################################################################################################################################################################
################# Calculation Block: Calculate the synthetic observation trajectory based on user inputs
### Create required arrays (empty) to store synthetic trajectory data points
TrX = []
TrY = []
TrZ = []
### Calculate the scale parameters using measurements
eta = (nu**3/epsMeas)**(1/4)        # Kolmogorov length scale [m]
Zscal = Nz*eta*resMet               # Scaled Z DNS domain dimension [m]
Xscal = Zscal*(Xl/Zl)               # Scaled X DNS domain dimension [m]
Yscal = Zscal*(Yl/Zl)               # Scaled Y DNS domain dimension [m]

### DNS resolution calculations
dx = Xl/Nx             # Grid Resolution in X (normalized)
dy = Yl/Ny             # Grid Resolution in Y (normalized)
dz = Zl/Nz             # Grid Resolution in Z (normalized)
dxscal = Xscal/Nx      # Grid Resolution in X (scaled) [m]
dyscal = Yscal/Nx      # Grid Resolution in Y (scaled) [m]
dzscal = Zscal/Nx      # Grid Resolution in Z (scaled) [m]

### Grid Co-ordinate calculations
# NOTE: The data is written at the grid edges i.e. Number of data points = Nx*Ny*Nz specified in sam.inp file
grid_x = np.zeros(Nx)
grid_y = np.zeros(Ny)
grid_z = np.zeros(Nz)
grid_xScal = np.zeros(Nx)
grid_yScal = np.zeros(Ny)
grid_zScal = np.zeros(Nz)
grid_x[0:Nx] = (-Xl/2)+(np.linspace(1,Nx,Nx)-1)*dx         # Grid point locations in X direction (normalized)
grid_y[0:Ny] = (-Yl/2)+(np.linspace(1,Ny,Ny)-1)*dy         # Grid point locations in Y direction (normalized)
grid_z[0:Nz] = (np.linspace(1,Nz,Nz)-1)*dz                 # Grid point locations in Z direction (normalized)
grid_xScal[0:Nx] = (-Xscal/2)+(np.linspace(1,Nx,Nx)-1)*dxscal         # Grid point locations in X direction (scaled) [m]
grid_yScal[0:Ny] = (-Yscal/2)+(np.linspace(1,Ny,Ny)-1)*dyscal         # Grid point locations in Y direction (scaled) [m]
grid_zScal[0:Nz] = (np.linspace(1,Nz,Nz)-1)*dzscal                    # Grid point locations in Z direction (scaled) [m]

### Calculate the number datapoints to sample - based on User inputs for the defined observation strategy 
### Get X,Y reference points
if locFlg == 1:         # Pick random, non-repeating locations
    Xref = random.sample(range(Nx),numTraj)
    Yref = random.sample(range(Ny),numTraj)
elif locFlg == 2:         # Manually chosen (X,Y) reference co-ordinates at Domain X,Y centre
    Xref = np.ones(numTraj)*Nx/2
    Yref = np.ones(numTraj)*Ny/2
### Loop over the number of trajectories
for i in range(0,numTraj):
    ### Call the function to generate each trajectory
    [Xtr,Ytr,Ztr] = trajInd(trajTyp,Xref[i],Yref[i],Nx,Ny,Nz)
    TrX.append(Xtr)
    TrY.append(Ytr)
    TrZ.append(Ztr)

##########################################################################################################################################################################
################# Data Extraction block
### List the datafiles in the operation directory
filNms = plname(flNm,form)
pts = Nx*Ny*Nz                  # The number of data points in each 3D dataset
### Read required points from file: Loop over the number of files
for i in range(0,len(filNms)):
    savenm = filNms[i][:-18] + datVar + filNms[i][-11:-4] + '.txt'
    tmpVar = np.zeros(numTraj,np.shape(TrX)[1])
    for j in range(0,len(filNms)):
        for k in range(0,len(filNms)):

    with open(name_list[n],'br') as file:
        for i in range(0,len(x_sampl)):
            # read u
            fld_u = 1
            pt_jmp = ((z_sampl[i]-1)*Nx*Ny)+((y_sampl[i]-1)*Nx)+(x_sampl[i]-1)+((fld_u-1)*pts)
            byt_jmp = int(pt_jmp*4)
            file.seek(byt_jmp)
            tl_u = file.tell()
            b=file.read(4)  
            u[i] = np.array(list(st.unpack('f', b)))

    data = np.column_stack([u,v,w,T,p,t_diss,eps])

    np.savetxt(save_fil,data,delimiter=',')

        '''