%% Script to plot plane data from DNS (of epsilon, u, v, and w field) and show the vertical sampled trajectories on it
%% Author: Abhiram Doddi
%% Date Composed: 18th July 2022

%% House-cleaning
clear
close all
clc

%% Inputs
% General Inputs: All inputs are mandatory
Xl = 16;         % Xl = Domain X dimension
Yl = 16;         % Yl = Domain Y dimension
Zl = 16;         % Zl = Domain Z dimension
Nx = 1536;       % Nx = Number of points in Domain X dimension
Ny = 1536;       % Ny = Number of points in Domain Y dimension
Nz = 1536;       % Nz = Number of points in Domain Z dimension
% Data Inputs: All inputs are mandatory
% NOTE: These inputs are required to scale the scale-normalized DNS datasets  
epsMeas = 2.84e-3;        % epsMeas = TKE dissipation rate [m^3s^-2]
nu = 4.0e-4;              % nu = Kinematic viscosity [m^2s^-1]
nuDNS = 2.0e-3;           % The kunematic viscosity set in the DNS (unscaled)
resMet = 0.64016;         % resMet = DNS resolution metric ( = grid spacing/ kolmogorov length scale)
balRt = 2.0;              % balRt = HYFLITS balloon descent rate [m/s]
dirtry = '/Users/script_away/Projects/Documents/MURI_modeling/SHIT/run04/analysis2m_smPln/';
plnFlg = [0 0 1 1];                 % binary switches to turn on plane viz for u, v, w, and epsilon
exectry = pwd;
saSp = 1;
% Plot Inputs: All inputs are mandatory
ftsz = 22;

%% Load Data
if plnFlg(1) == 1               % for u
    flPlnU = strcat(dirtry,'yz1_0008_U.txt');
    temp = table2array(readtable(flPlnU));
    PlnU = reshape(temp,[Nx+1,Ny+1]);
    clear temp
end

if plnFlg(2) == 1               % for v
    flPlnV = strcat(dirtry,'yz1_0008_V.txt');
    temp = table2array(readtable(flPlnV));
    PlnV = reshape(temp,[Nx+1,Ny+1]);
    clear temp
end

if plnFlg(3) == 1               % for w
    flPlnW = strcat(dirtry,'yz1_0008_W.txt');
    temp = table2array(readtable(flPlnW));
    PlnW = reshape(temp,[Nx+1,Ny+1]);
    clear temp
end

if plnFlg(4) == 1               % for epsilon
    flPlnE = strcat(dirtry,'yz1_0008_E.txt');
    temp = table2array(readtable(flPlnE));
    PlnE = reshape(temp,[Nx+1,Ny+1]);
    clear temp
end
% also load the trajectory data
flrfPt = strcat(dirtry,'refPts.txt');
refXY = table2array(readtable(flrfPt));         % array containing reference grid indices (not co-ordinates) for each trajectory

%% calculate the grid
% calculate DNS resolution parameters
dx = Xl/Nx;             % Grid Resolution in X (normalized)
dy = Yl/Ny;             % Grid Resolution in Y (normalized)
dz = Zl/Nz;             % Grid Resolution in Z (normalized)
[GridXYx,GridXYy] = meshgrid((-Xl/2):dx:(Xl/2),(-Yl/2):dy:(Yl/2));
[GridXZx,GridXZz] = meshgrid((-Xl/2):dx:(Xl/2),0:dz:Zl);
[GridYZy,GridYZz] = meshgrid((-Yl/2):dy:(Yl/2),0:dz:Zl);

%% create trajectories
trajX = refXY(:,1)*dx + (-Xl/2);
trajY = refXY(:,2)*dy + (-Yl/2);
% create Z coordinates for each trajectory
trajZ = 0:dz:Zl;

%% Plot the planes and trajectories
figure(1)
clf
s = surface(GridYZy,GridYZz,PlnE);
hold on
xlabel('Y')
ylabel('Z')
s.EdgeColor = 'none';
colorbar
caxis([0 1])
for i = 1:1:length(trajY)
    Ycrd = ones([length(trajZ),1]).*trajY(i);
    Zcrd = ones([length(trajZ),1]).*max(max(PlnE));
    figure(1)
    plot3(Ycrd,trajZ,Zcrd,'r','LineWidth',0.5)
    clear Ycrd Zcrd
end