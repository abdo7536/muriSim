################# Script to sample data on a user defined trajectory from a subvol***.dat file
################# Author: Dr. Abhiram Doddi
################# Date Created: 1st April 2022
################# Last updated: 6th April 2022

##########################################################################################################################################################################
################# General Inputs: All inputs are mandatory
### Xl = Domain X dimension
Xl = 16
### Yl = Domain Y dimension
Yl = 16 
### Zl = Domain Z dimension
Zl = 16 
### Nx = Number of points in Domain X dimension
Nx = 1536 
### Ny = Number of points in Domain Y dimension
Ny = 1536
### Nz = Number of points in Domain Z dimension
Nz = 1536 
### numTraj = Number of trajectories to sample
numTraj = 1
### locFlg = Flag to control the (X,Y) location for sampling; 1 = choose random (x,y) location for each trajectory; 21 = randomly chosen Y reference co-ordinates at Domain X centre ; 22 = randomly chosen X reference co-ordinates at Domain Y centre; 3 = Read X,Y reference co-ordinates for each trajectory from a .txt file
locFlg = 1 
### allZ = Flag to sample every point in Z of the domain as 1 sample
allZ = 1
### datVar = The data variable to be stored; NOTE: This is the variable extracted from the DNS 3D datafield and stored at the user defined destination/directory
datVar = 'E'
### trajDir = The direction in which the synthetic observer traverses; NOTE: choice of X, Y and Z directions of synthetic observation traverse; NOTE: ONLY USED FOR NAMING THE FILES
trajDir = 'y'

################# Data Inputs: All inputs are mandatory
### NOTE: These inputs are required to scale the scale-normalized DNS datasets  
### epsMeas = TKE dissipation rate [m^3s^-2]
epsMeas = 3.738
### nu = Kinematic viscosity [m^2s^-1]
nu = 4.0e-4
### resMet = DNS resolution metric ( = grid spacing/ kolmogorov length scale)
resMet = 0.64016 
### balRt = HYFLITS balloon descent rate [m/s]
balRt = 2

################# Function Specific Inputs: All inputs are mandatory
########## function name: 'lodatsin'
### flNm = The string name whose file name group needs to be created (including the path to the directory) [type - string]
flNm = '/p/work1/abdo7536/SHIT04/analysis/getData/subvol' + datVar + '*'
### form = the file format (input exactly as it appears in the filenames) [type - string]
form = 'dat'
### flNm = Input file name identifying character string (including the directory location)

################# Function Specific Inputs: All inputs are mandatory
########## function name: 'trajInd'
### trajTyp = Synthetic Observation trajectory type; 1 = balloon-like vertical trajectory (descending); 2 = Horizontal Sampling in X; 3 = Horizontal Sampling in Y
trajTyp = [3]

##########################################################################################################################################################################
################# Call custom function module developed to execute this program
### NOTE: This module inturn loads the python inbuilt modules as needed
try:
    from userFns import *
    print("Program required custom and python module is loaded...")
except ImportError:
    print("Custom library not installed.... Install and try again")

##########################################################################################################################################################################
################# Calculation Block: Calculate the DNS datafield scaling and extract the 3D meshgrid
### Create required arrays (empty) to store synthetic trajectory data points and the meshgrid data points
TrX = []
TrY = []
TrZ = []
### Calculate the scale parameters using measurements
eta = (nu**3/epsMeas)**(1/4)        # Kolmogorov length scale [m]
Zscal = Nz*eta*resMet               # Scaled Z DNS domain dimension [m]
Xscal = Zscal*(Xl/Zl)               # Scaled X DNS domain dimension [m]
Yscal = Zscal*(Yl/Zl)               # Scaled Y DNS domain dimension [m]
print('Xl = ',Xscal,' [m]\n','Yl = ',Yscal,' [m]\n','Zl = ',Zscal,' [m]')

### DNS resolution calculations
dx = Xl/Nx             # Grid Resolution in X (normalized)
dy = Yl/Ny             # Grid Resolution in Y (normalized)
dz = Zl/Nz             # Grid Resolution in Z (normalized)
dxscal = Xscal/Nx      # Grid Resolution in X (scaled) [m]
dyscal = Yscal/Ny      # Grid Resolution in Y (scaled) [m]
dzscal = Zscal/Nz      # Grid Resolution in Z (scaled) [m]

### Grid Co-ordinate calculations
# NOTE: The data is written at the grid edges i.e. Number of data points = Nx*Ny*Nz specified in sam.inp file
grid_x = np.zeros((1,Nx))
grid_y = np.zeros((1,Ny))
grid_z = np.zeros((1,Nz))
grid_xScal = np.zeros((1,Nx))
grid_yScal = np.zeros((1,Ny))
grid_zScal = np.zeros((1,Nz))
grid_x[0][0:Nx] = (-Xl/2)+(np.linspace(1,Nx,Nx)-1)*dx         # Grid point locations in X direction (normalized)
grid_y[0][0:Ny] = (-Yl/2)+(np.linspace(1,Ny,Ny)-1)*dy         # Grid point locations in Y direction (normalized)
grid_z[0][0:Nz] = (np.linspace(1,Nz,Nz)-1)*dz                 # Grid point locations in Z direction (normalized)
grid_xScal[0][0:Nx] = (-Xscal/2)+(np.linspace(1,Nx,Nx)-1)*dxscal         # Grid point locations in X direction (scaled) [m]
grid_yScal[0][0:Ny] = (-Yscal/2)+(np.linspace(1,Ny,Ny)-1)*dyscal         # Grid point locations in Y direction (scaled) [m]
grid_zScal[0][0:Nz] = (np.linspace(1,Nz,Nz)-1)*dzscal                    # Grid point locations in Z direction (scaled) [m]
### Create a meshgrid (The meshgrid - [X,Y,Z] grid points - should be exported in the File)
[Ygrid,Xgrid,Zgrid] = np.meshgrid(grid_yScal,grid_xScal,grid_zScal)

##########################################################################################################################################################################
################# Trajectory points calculation block
### List the datafiles in the operation directory
filNms = plname(flNm,form)
dirTry = filNms[0][:-18]

### Calculate the number datapoints to sample - based on User inputs for the defined observation strategy 
### Get X,Y reference points
if locFlg == 1:         # Pick random, non-repeating locations
    tmpXref = random.sample(range(1,Nx),numTraj)
    tmpYref = random.sample(range(1,Ny),numTraj)
    Xref = np.transpose(np.reshape(tmpXref,(1,len(tmpXref))))
    Yref = np.transpose(np.reshape(tmpYref,(1,len(tmpYref))))
    strctRef = np.column_stack((Xref,Yref))
    refFlNm = dirTry + 'refPts.txt'
    np.savetxt(refFlNm,strctRef,delimiter=',')
elif locFlg == 21:        # randomly chosen Y reference co-ordinates at Domain X centre
    tmpXref = np.ones(numTraj)*Nx/2
    tmpYref = random.sample(range(1,Ny),numTraj)
    Xref = np.transpose(np.reshape(tmpXref,(1,len(tmpXref))))
    Yref = np.transpose(np.reshape(tmpYref,(1,len(tmpYref))))
    strctRef = np.column_stack((Xref,Yref))
    refFlNm = dirTry + 'refPts.txt'
    np.savetxt(refFlNm,strctRef,delimiter=',')
elif locFlg == 22:        # randomly chosen X reference co-ordinates at Domain Y centre
    tmpXref = random.sample(range(1,Nx),numTraj)
    tmpYref = np.ones(numTraj)*Ny/2
    Xref = np.transpose(np.reshape(tmpXref,(1,len(tmpXref))))
    Yref = np.transpose(np.reshape(tmpYref,(1,len(tmpYref))))
    strctRef = np.column_stack((Xref,Yref))
    refFlNm = dirTry + 'refPts.txt'
    np.savetxt(refFlNm,strctRef,delimiter=',')
elif locFlg == 3:         # Read (X,Y) reference co-ordinates from a text file
    refFlNm = dirTry + 'refPts.txt'
    with open(refFlNm) as da:
        refPts = da.readlines()
        Xref = np.zeros((len(refPts),1))
       	Yref = np.zeros((len(refPts),1))
        for i in range(0,len(refPts)):
            lnRef = refPts[i]
            Xref[i] = int(float(lnRef[0:24]))
            Yref[i] = int(float(lnRef[26:]))

### Compute interval sets in case of horizontal sampling trajectories
if allZ == 1:
    trajTyp.append(int(1))
    trajTyp.append(int(Nz))
elif allZ == 0:
    trajTyp.append(math.floor(Zscal/balRt))         # calculate the maximum number of (full) intervals in one single descent through the DNS datafield [integer]
    trajTyp.append(math.floor(balRt/dzscal))        # calculate the maximum number of grid points per each interval

smplPts = trajTyp[1]*trajTyp[2]                 # Compute the number of points to sample from the (scaled) DNS dataset

print('The number of sample intervals and sample points per interval = ', trajTyp)
print('The number of samples gathered per trajectory = ', smplPts)

### Do not proceed if the total number of points per interval or total points sampled is odd
if (trajTyp[2] % 2) != 0 or (smplPts % 2) != 0:
    print('ERROR:The number of points per interval is ODD or the total number of points sampled is ODD')
    print('WARNING: These need to be EVEN')
    sys.exit()

### Loop over the number of trajectories
for i in range(0,numTraj):
    ### Call the function to generate each trajectory
    [Xtr,Ytr,Ztr] = trajInd(trajTyp,Xref[i],Yref[i],Nx,Ny,Nz,smplPts)
    TrX.append(Xtr[0])
    TrY.append(Ytr[0])
    TrZ.append(Ztr[0])
    
##########################################################################################################################################################################
################# Data Extraction block
### Read required points from file: Loop over the number of files
for i in range(0,len(filNms)):
    ### Initialize arrays to store the trajectory data and grid data
    tmpVar = np.zeros((numTraj,np.shape(TrX)[1]))
    GrdX = np.zeros((numTraj,np.shape(TrX)[1]))
    GrdY = np.zeros((numTraj,np.shape(TrY)[1]))
    GrdZ = np.zeros((numTraj,np.shape(TrZ)[1]))
    ### Create a temporary variable to store all extracted data points on the trajectory
    ### Open the .dat data file
    with open(filNms[i],'br') as fl:
        for j in range(0,numTraj):
            for k in range(0,np.shape(TrX)[1]):    
                tmpVar[j][k] = datGet(TrX[j][k],TrY[j][k],TrZ[j][k],Nx,Ny,Nz,fl)
                GrdX[j][k] = Xgrid[int(TrX[j][k])-1][int(TrY[j][k])-1][int(TrZ[j][k])-1]
                GrdY[j][k] = Ygrid[int(TrX[j][k])-1][int(TrY[j][k])-1][int(TrZ[j][k])-1]
                GrdZ[j][k] = Zgrid[int(TrX[j][k])-1][int(TrY[j][k])-1][int(TrZ[j][k])-1]
            print('Gathered data for Trajectory = ',j)
        ### Create a table containing the extracted data and the grid centre co-ordinates in Z
        strct1 = np.transpose(tmpVar)
        strct2 = np.transpose(GrdX)
        strct3 = np.transpose(GrdY)
        strct4 = np.transpose(GrdZ)
        ### create a savename
        saveDat = filNms[i][:-18] + datVar + trajDir + "{:06d}".format(numTraj) + filNms[i][-11:-4] + '.txt'
        saveGX = filNms[i][:-18] + 'GridX' + datVar + trajDir + "{:06d}".format(numTraj) + filNms[i][-11:-4] + '.txt'
        saveGY = filNms[i][:-18] + 'GridY' + datVar + trajDir + "{:06d}".format(numTraj) + filNms[i][-11:-4] + '.txt'
        saveGZ = filNms[i][:-18] + 'GridZ' + datVar + trajDir + "{:06d}".format(numTraj) + filNms[i][-11:-4] + '.txt'
        ### save as textfile
        np.savetxt(saveDat,strct1,delimiter=',')
        np.savetxt(saveGX,strct2,delimiter=',')
        np.savetxt(saveGY,strct3,delimiter=',')
        np.savetxt(saveGZ,strct4,delimiter=',')
