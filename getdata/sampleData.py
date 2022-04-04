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
numTraj = 1
### locFlg = Flag to control the (X,Y) location for sampling; 1 = choose random (x,y) location for each trajectory; 2 = Uniformly spaced (x,y) location; 3 = Specific (x,y) location
locFlg = 1
# NOTE: If locFlg = 3 then manually set the (x,y) location inside the code!!

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
flNm = '/Users/script_away/Projects/Documents/MURI_modeling/GWBData/subvol*'
### form = the file format (input exactly as it appears in the filenames) [type - string]
form = 'dat'
### flNm = Input file name identifying character string (including the directory location)

################# Function Specific Inputs: All inputs are mandatory
########## function name: 'trajInd'
### trajTyp = Synthetic Observation trajectory type; 1 = balloon-like vertical trajectory (descending); 2 = Helical trajectory (descending); 3 = Other trajectory/sampling strategy (descending)
trajTyp = 1

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



### Calculate vertical sampling
x_sampl = np.ones(Nz)*(Nx/2)+1
y_sampl = np.ones(Nz)*(Ny/2)+1
z_sampl = np.linspace(1,Nz,Nz)

##########################################################################################################################################################################
################# Data Extraction block
## Read required points from file
pts = Nx*Ny*Nz
name_list = plname(flNm,form)
for n in range(0,9):
    tp = name_list[n]
    save_fil = tp[0:32] + 'txt_files_tru_vert/' + tp[32:-4] + '.txt'
    u = np.zeros(len(x_sampl))
    v = np.zeros(len(x_sampl))
    w = np.zeros(len(x_sampl))
    T = np.zeros(len(x_sampl))
    p = np.zeros(len(x_sampl))
    t_diss = np.zeros(len(x_sampl))
    eps = np.zeros(len(x_sampl))

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

    with open(name_list[n],'br') as file:
        for i in range(0,len(x_sampl)):
            # read v
            fld_v = 2
            pt_jmp = ((z_sampl[i]-1)*Nx*Ny)+((y_sampl[i]-1)*Nx)+(x_sampl[i]-1)+((fld_v-1)*pts)
            byt_jmp = int(pt_jmp*4)
            file.seek(byt_jmp)
            tl_v = file.tell()
            b=file.read(4)  
            v[i] = np.array(list(st.unpack('f', b)))

    with open(name_list[n],'br') as file:
        for i in range(0,len(x_sampl)):
            # read w
            fld_w = 3
            pt_jmp = ((z_sampl[i]-1)*Nx*Ny)+((y_sampl[i]-1)*Nx)+(x_sampl[i]-1)+((fld_w-1)*pts)
            byt_jmp = int(pt_jmp*4)
            file.seek(byt_jmp)
            tl_w = file.tell()
            b=file.read(4)  
            w[i] = np.array(list(st.unpack('f', b)))

    with open(name_list[n],'br') as file:
        for i in range(0,len(x_sampl)):
            # read T
            fld_T = 4
            pt_jmp = ((z_sampl[i]-1)*Nx*Ny)+((y_sampl[i]-1)*Nx)+(x_sampl[i]-1)+((fld_T-1)*pts)
            byt_jmp = int(pt_jmp*4)
            file.seek(byt_jmp)
            tl_T = file.tell()
            b=file.read(4)  
            T[i] = np.array(list(st.unpack('f', b)))

    with open(name_list[n],'br') as file:
        for i in range(0,len(x_sampl)):
            # read p
            fld_p = 5
            pt_jmp = ((z_sampl[i]-1)*Nx*Ny)+((y_sampl[i]-1)*Nx)+(x_sampl[i]-1)+((fld_p-1)*pts)
            byt_jmp = int(pt_jmp*4)
            file.seek(byt_jmp)
            tl_p = file.tell()
            b=file.read(4)  
            p[i] = np.array(list(st.unpack('f', b)))

    with open(name_list[n],'br') as file:
        for i in range(0,len(x_sampl)):
            # read t_diss
            fld_tdiss = 6
            pt_jmp = ((z_sampl[i]-1)*Nx*Ny)+((y_sampl[i]-1)*Nx)+(x_sampl[i]-1)+((fld_tdiss-1)*pts)
            byt_jmp = int(pt_jmp*4)
            file.seek(byt_jmp)
            tl_tdiss = file.tell()
            b=file.read(4)  
            t_diss[i] = np.array(list(st.unpack('f', b)))

    with open(name_list[n],'br') as file:
        for i in range(0,len(x_sampl)):
            # read eps
            fld_eps = 7
            pt_jmp = ((z_sampl[i]-1)*Nx*Ny)+((y_sampl[i]-1)*Nx)+(x_sampl[i]-1)+((fld_eps-1)*pts)
            byt_jmp = int(pt_jmp*4)
            file.seek(byt_jmp)
            tl_wps = file.tell()
            b=file.read(4)  
            eps[i] = np.array(list(st.unpack('f', b)))

    data = np.column_stack([u,v,w,T,p,t_diss,eps])

    np.savetxt(save_fil,data,delimiter=',')

