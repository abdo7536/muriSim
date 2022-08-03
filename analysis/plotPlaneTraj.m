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
saSp = 0;
% Plot Inputs: All inputs are mandatory
ftsz = 18;

%% Setup calculations
% Calculate scale parameters for DNS data
eta = (nu^3/epsMeas)^(1/4);         % Kolmogorov length scale [m]
tau = (nu/epsMeas)^(1/2);           % Kolmogorov time scale [s]
mu = eta/tau;                       % Kolmogorov velocity scale [m/s]
epsStar = nu*(mu/eta)^2;            % Scaling to be applied for TKE dissipation rate [m^3/s^2]

%% Load Data
if plnFlg(1) == 1               % for u
    flPlnU = strcat(dirtry,'yz1_0008_U.txt');
    temp = table2array(readtable(flPlnU));
    PlnU = transpose(reshape(temp,[Nx+1,Ny+1])).*mu;
    clear temp
end

if plnFlg(2) == 1               % for v
    flPlnV = strcat(dirtry,'yz1_0008_V.txt');
    temp = table2array(readtable(flPlnV));
    PlnV = transpose(reshape(temp,[Nx+1,Ny+1])).*mu;
    clear temp
end

if plnFlg(3) == 1               % for w
    flPlnW = strcat(dirtry,'yz1_0008_W.txt');
    temp = table2array(readtable(flPlnW));
    PlnW = transpose(reshape(temp,[Nx+1,Ny+1])).*mu;
    clear temp
end

if plnFlg(4) == 1               % for epsilon
    flPlnE = strcat(dirtry,'yz1_0008_E.txt');
    temp = table2array(readtable(flPlnE));
    PlnE = transpose(reshape(temp,[Nx+1,Ny+1])).*epsStar;
    clear temp
end
% also load the trajectory data
flrfPt = strcat(dirtry,'refPts.txt');
refXY = table2array(readtable(flrfPt));         % array containing reference grid indices (not co-ordinates) for each trajectory
% load 
flPDF = strcat(dirtry,'subvol1_004000_KE-diss-rate_pdf.txt');
dataPDF = table2array(readtable(flPDF));
bin = dataPDF(6:end,3);
PDF = dataPDF(6:end,4);
samples = dataPDF(6:end,5);

%% calculate the grid
% Calculate the DNS domain extents using Kolmogorov scaling
Zscal = Nz*eta*resMet;              % Scaled Z DNS domain dimension [m]
Xscal = Zscal*(Xl/Zl);              % Scaled X DNS domain dimension [m]
Yscal = Zscal*(Yl/Zl);              % Scaled Y DNS domain dimension [m]
% calculate DNS resolution parameters
dx = Xl/Nx;             % Grid Resolution in X (normalized)
dy = Yl/Ny;             % Grid Resolution in Y (normalized)
dz = Zl/Nz;             % Grid Resolution in Z (normalized)
dxscal = Xscal/Nx;      % Grid Resolution in X (scaled) [m]
dyscal = Yscal/Ny;      % Grid Resolution in Y (scaled) [m]
dzscal = Zscal/Nz;      % Grid Resolution in Z (scaled) [m]
[GridXYx,GridXYy] = meshgrid((-Xscal/2):dxscal:(Xscal/2),(-Yscal/2):dyscal:(Yscal/2));
[GridXZx,GridXZz] = meshgrid((-Xscal/2):dxscal:(Xscal/2),0:dzscal:Zscal);
[GridYZy,GridYZz] = meshgrid((-Yscal/2):dyscal:(Yscal/2),0:dzscal:Zscal);

%% create trajectories
trajX = refXY(:,1)*dxscal + (-Xscal/2);
trajY = refXY(:,2)*dyscal + (-Yscal/2);
% create Z coordinates for each trajectory
trajZ = 0:dzscal:Zscal;

%% Plot the planes and trajectories
if plnFlg(1) == 1               % for u
    figure(1)
    clf
    s = surface(GridYZy,GridYZz,PlnU);
    hold on
    xlabel('Y','FontSize', ftsz)
    ylabel('Z','FontSize', ftsz)
    s.EdgeColor = 'none';
    colorbar
    %caxis([0 1])
    title('U on YZ plane at X=0','FontSize', ftsz)
    for i = 1:1:length(trajY)
        Ycrd = ones([length(trajZ),1]).*trajY(i);
        Zcrd = ones([length(trajZ),1]).*max(max(PlnU));
        figure(1)
        plot3(Ycrd,trajZ,Zcrd,'r','LineWidth',0.25)
        text(trajY(i),trajZ(end-(15*i)),num2str(i),'FontSize', ftsz-4)
        clear Ycrd Zcrd
    end
    xlim([min(min(GridYZy)) max(max(GridYZy))])
    ylim([min(min(GridYZz)) max(max(GridYZz))])
    if saSp == 1
        cd(dirtry)
        set(gcf, 'Position',[100, 100, 1300, 1000])
        set(gcf, 'Position',[100, 100, 600, 450])
        savefig('U_YZpln.fig')
        cd(exectry)
    end
end

if plnFlg(2) == 1               % for v
    figure(2)
    clf
    s = surface(GridYZy,GridYZz,PlnV);
    hold on
    xlabel('Y','FontSize', ftsz)
    ylabel('Z','FontSize', ftsz)
    s.EdgeColor = 'none';
    colorbar
    %caxis([0 1])
    title('V on YZ plane at X=0','FontSize', ftsz)
    for i = 1:1:length(trajY)
        Ycrd = ones([length(trajZ),1]).*trajY(i);
        Zcrd = ones([length(trajZ),1]).*max(max(PlnV));
        figure(2)
        plot3(Ycrd,trajZ,Zcrd,'r','LineWidth',0.25)
        text(trajY(i),trajZ(end-(15*i)),num2str(i),'FontSize', ftsz-4)
        clear Ycrd Zcrd
    end
    xlim([min(min(GridYZy)) max(max(GridYZy))])
    ylim([min(min(GridYZz)) max(max(GridYZz))])
    if saSp == 1
        cd(dirtry) 
        set(gcf, 'Position',[100, 100, 1300, 1000])
        savefig('V_YZpln.fig')
        cd(exectry)
    end
end

if plnFlg(3) == 1               % for w
    figure(3)
    clf
    s = surface(GridYZy,GridYZz,PlnW);
    hold on
    xlabel('Y','FontSize', ftsz)
    ylabel('Z','FontSize', ftsz)
    s.EdgeColor = 'none';
    colorbar
    %caxis([0 1])
    title('W on YZ plane at X=0','FontSize', ftsz)
    for i = 1:1:length(trajY)
        Ycrd = ones([length(trajZ),1]).*trajY(i);
        Zcrd = ones([length(trajZ),1]).*max(max(PlnW));
        figure(3)
        plot3(Ycrd,trajZ,Zcrd,'r','LineWidth',0.25)
        text(trajY(i),trajZ(end-(15*i)),num2str(i),'FontSize', ftsz-4)
        clear Ycrd Zcrd
    end
    xlim([min(min(GridYZy)) max(max(GridYZy))])
    ylim([min(min(GridYZz)) max(max(GridYZz))])
    if saSp == 1
        cd(dirtry) 
        set(gcf, 'Position',[100, 100, 1300, 1000])
        savefig('W_YZpln.fig')
        cd(exectry)
    end
end

if plnFlg(4) == 1               % for epsilon
    figure(4)
    clf
    s = surface(GridYZy,GridYZz,log10(PlnE));
    hold on
    xlabel('Y','FontSize', ftsz)
    ylabel('Z','FontSize', ftsz)
    s.EdgeColor = 'none';
    colorbar
    caxis([-6 -1])
    title('\epsilon_{DNS} on YZ plane at X=0','FontSize', ftsz)
    for i = 1:1:length(trajY)
        Ycrd = ones([length(trajZ),1]).*trajY(i);
        Zcrd = ones([length(trajZ),1]).*max(max(PlnE));
        figure(4)
        plot3(Ycrd,trajZ,Zcrd,'r','LineWidth',0.25)
        text(trajY(i),trajZ(end-(15*i)),num2str(i),'FontSize', ftsz-4)
        clear Ycrd Zcrd
    end
    xlim([min(min(GridYZy)) max(max(GridYZy))])
    ylim([min(min(GridYZz)) max(max(GridYZz))])
    if saSp == 1
        cd(dirtry)
        set(gcf, 'Position',[100, 100, 1300, 1000])
        savefig('E_YZpln.fig')
        cd(exectry)
    end
end