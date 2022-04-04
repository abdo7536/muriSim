################# Script to sample data on a user defined trajectory from a subvol***.dat file
################# Author: Dr. Abhiram Doddi
################# Date Created: 1st April 2022
################# Last updated: 1st April 2022

##########################################################################################################################################################################
################# Inputs: All inputs are mandatory
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

##########################################################################################################################################################################
################# All the dependencies to run inbuilt or custom functions that are called within the listed functions below need to be imported beforehand
### This section imports all the required function
try:
    import numpy as np
    print("numpy module is loaded...")
except ImportError:
    print("numpy module not installed...install numpy and try again")
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

##########################################################################################################################################################################
################# Calculation Block: Calculate the synthetic observation trajectory based on user inputs
### Scale parameters
h = 1                  # Shear depth
u = 1                  # Shear velocity
t = h/u                # Shear time scale

## UAS parameters
uas_srate = 800            # sampling rate of the UAS [Hz]
uas_speed = 14.5           # How fast the UAS travels horizontally [m/s]
uas_climb_rate = 2.5       # The rate at which the UAS nominally ascends or descends [m/s]
dhelix = 100               # choose between 75, 100, 150 [m]

## Guidance from Measurement
ldepth = 500                # choice between 400-600 [m]
shear_vel = 8               # no more than 10 [m/s]
t_meas = ldepth/shear_vel   # calculated based on nominally observed structure characteristics [t]

## Convert UAS parameters in scales of Simulation 
l_star = (6*h)/ldepth
vel_star = u/shear_vel

uasinsim_speed = vel_star * uas_speed
uasinsim_climb_rate = vel_star * uas_climb_rate
uasinsim_dhelix = l_star * dhelix
uas_dt = 1/uas_srate

## Input variables
Xl = 13*h                  # Box Size in X
Yl = 13*h                  # Box Size in Y
Zl = 39*h                  # Box Size in Z
Nx = 576                   # Grid Points in X
Ny = 576                   # Grid Points in Y
Nz = 1728                  # Grid Points in Z

## Fundamental Calculations
dx = Xl/Nx             # Grid Resolution in X
dy = Yl/Ny             # Grid Resolution in Y
dz = Zl/Nz             # Grid Resolution in Z

## Calculate the Grid Co-ordinates
grid_x = np.zeros(Nx)
grid_y = np.zeros(Ny)
grid_z = np.zeros(Nz)
grid_x[0:Nx] = (-Xl/2)+(np.linspace(1,Nx,Nx)-1)*dx         # Grid point locations in X direction
grid_y[0:Ny] = (-Yl/2)+(np.linspace(1,Ny,Ny)-1)*dy         # Grid point locations in Y direction
grid_z[0:Nz] = (np.linspace(1,Nz,Nz)-1)*dz                 # Grid point locations in Z direction

grid_cent_x = grid_x[0:-1] + dx/2
grid_cent_y = grid_y[0:-1] + dy
grid_cent_z = grid_z[0:-1] + dz

## Calculate tru vertical sampling
x_sampl = np.ones(Nz)*(Nx/2)+1
y_sampl = np.ones(Nz)*(Ny/2)+1
z_sampl = np.linspace(1,Nz,Nz)

##########################################################################################################################################################################
################# Extraction block
## Read required points from file
pts = Nx*Ny*Nz
name_list = plname('/p/work1/abdo7536/thesis2b/dir1/subvol1_0','dat')
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

