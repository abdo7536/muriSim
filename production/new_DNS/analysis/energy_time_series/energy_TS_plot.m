%% Script to plot timeseries of energy and dissipation rate for SHITDNS run segments

%% House cleaning
clear
close all
clc

%% load data
load energy.mat
load resolution.mat

%% Plot data
figure(1)
clf
plot(table2array(energy01(:,2)),table2array(energy01(:,4)),'k','LineWidth',1)
hold on
plot(table2array(energy02(:,2)),table2array(energy02(:,4)),'k','LineWidth',1)
plot(table2array(energy03(:,2)),table2array(energy03(:,4)),'k','LineWidth',1)
plot(table2array(energy04(:,2)),table2array(energy04(:,4)),'k','LineWidth',1)
plot(table2array(energy05(:,2)),table2array(energy05(:,4)),'k','LineWidth',1)
plot(table2array(energy06(:,2)),table2array(energy06(:,4)),'k','LineWidth',1)
plot(table2array(energy07(:,2)),table2array(energy07(:,4)),'k','LineWidth',1)
plot(table2array(energy08(:,2)),table2array(energy08(:,4)),'k','LineWidth',1)
plot(table2array(energy09(:,2)),table2array(energy09(:,4)),'k','LineWidth',1)
plot(table2array(energy10(:,2)),table2array(energy10(:,4)),'k','LineWidth',1)
plot(table2array(energy11(:,2)),table2array(energy11(:,4)),'k','LineWidth',1)
plot(table2array(energy12(:,2)),table2array(energy12(:,4)),'k','LineWidth',1)
plot(table2array(energy13(:,2)),table2array(energy13(:,4)),'k','LineWidth',1)
plot(table2array(energy14(:,2)),table2array(energy14(:,4)),'k','LineWidth',1)
grid on
grid Minor
xlabel('time, [unitless]')
ylabel('energy dissipation rate, [unitless]')

figure(2)
clf
plot(table2array(energy01(:,2)),table2array(energy01(:,3)),'k','LineWidth',1)
hold on
plot(table2array(energy02(:,2)),table2array(energy02(:,3)),'k','LineWidth',1)
plot(table2array(energy03(:,2)),table2array(energy03(:,3)),'k','LineWidth',1)
plot(table2array(energy04(:,2)),table2array(energy04(:,3)),'k','LineWidth',1)
plot(table2array(energy05(:,2)),table2array(energy05(:,3)),'k','LineWidth',1)
plot(table2array(energy06(:,2)),table2array(energy06(:,3)),'k','LineWidth',1)
plot(table2array(energy07(:,2)),table2array(energy07(:,3)),'k','LineWidth',1)
plot(table2array(energy08(:,2)),table2array(energy08(:,3)),'k','LineWidth',1)
plot(table2array(energy09(:,2)),table2array(energy09(:,3)),'k','LineWidth',1)
plot(table2array(energy10(:,2)),table2array(energy10(:,3)),'k','LineWidth',1)
plot(table2array(energy11(:,2)),table2array(energy11(:,3)),'k','LineWidth',1)
plot(table2array(energy12(:,2)),table2array(energy12(:,3)),'k','LineWidth',1)
plot(table2array(energy13(:,2)),table2array(energy13(:,3)),'k','LineWidth',1)
plot(table2array(energy14(:,2)),table2array(energy14(:,3)),'k','LineWidth',1)
grid on
grid Minor
xlabel('time, [unitless]')
ylabel('energy, [unitless]')

figure(3)
clf
plot(table2array(resolution01(:,2)),table2array(resolution01(:,3)),'*k','LineWidth',1)
hold on
plot(table2array(resolution01(:,2)),table2array(resolution01(:,4)),'*b','LineWidth',1)
plot(table2array(resolution01(:,2)),table2array(resolution01(:,5)),'*r','LineWidth',1)
plot(table2array(resolution02(:,2)),table2array(resolution02(:,3)),'*k','LineWidth',1)
plot(table2array(resolution03(:,2)),table2array(resolution03(:,3)),'*k','LineWidth',1)
plot(table2array(resolution04(:,2)),table2array(resolution04(:,3)),'*k','LineWidth',1)
plot(table2array(resolution05(:,2)),table2array(resolution05(:,3)),'*k','LineWidth',1)
plot(table2array(resolution06(:,2)),table2array(resolution06(:,3)),'*k','LineWidth',1)
plot(table2array(resolution07(:,2)),table2array(resolution07(:,3)),'*k','LineWidth',1)
plot(table2array(resolution08(:,2)),table2array(resolution08(:,3)),'*k','LineWidth',1)
plot(table2array(resolution09(:,2)),table2array(resolution09(:,3)),'*k','LineWidth',1)
plot(table2array(resolution10(:,2)),table2array(resolution10(:,3)),'*k','LineWidth',1)
plot(table2array(resolution11(:,2)),table2array(resolution11(:,3)),'*k','LineWidth',1)
plot(table2array(resolution12(:,2)),table2array(resolution12(:,3)),'*k','LineWidth',1)
plot(table2array(resolution13(:,2)),table2array(resolution13(:,3)),'*k','LineWidth',1)
plot(table2array(resolution14(:,2)),table2array(resolution14(:,3)),'*k','LineWidth',1)
plot(table2array(resolution02(:,2)),table2array(resolution02(:,4)),'*b','LineWidth',1)
plot(table2array(resolution03(:,2)),table2array(resolution03(:,4)),'*b','LineWidth',1)
plot(table2array(resolution04(:,2)),table2array(resolution04(:,4)),'*b','LineWidth',1)
plot(table2array(resolution05(:,2)),table2array(resolution05(:,4)),'*b','LineWidth',1)
plot(table2array(resolution06(:,2)),table2array(resolution06(:,4)),'*b','LineWidth',1)
plot(table2array(resolution07(:,2)),table2array(resolution07(:,4)),'*b','LineWidth',1)
plot(table2array(resolution08(:,2)),table2array(resolution08(:,4)),'*b','LineWidth',1)
plot(table2array(resolution09(:,2)),table2array(resolution09(:,4)),'*b','LineWidth',1)
plot(table2array(resolution10(:,2)),table2array(resolution10(:,4)),'*b','LineWidth',1)
plot(table2array(resolution11(:,2)),table2array(resolution11(:,4)),'*b','LineWidth',1)
plot(table2array(resolution12(:,2)),table2array(resolution12(:,4)),'*b','LineWidth',1)
plot(table2array(resolution13(:,2)),table2array(resolution13(:,4)),'*b','LineWidth',1)
plot(table2array(resolution14(:,2)),table2array(resolution14(:,4)),'*b','LineWidth',1)
plot(table2array(resolution02(:,2)),table2array(resolution02(:,5)),'*r','LineWidth',1)
plot(table2array(resolution03(:,2)),table2array(resolution03(:,5)),'*r','LineWidth',1)
plot(table2array(resolution04(:,2)),table2array(resolution04(:,5)),'*r','LineWidth',1)
plot(table2array(resolution05(:,2)),table2array(resolution05(:,5)),'*r','LineWidth',1)
plot(table2array(resolution06(:,2)),table2array(resolution06(:,5)),'*r','LineWidth',1)
plot(table2array(resolution07(:,2)),table2array(resolution07(:,5)),'*r','LineWidth',1)
plot(table2array(resolution08(:,2)),table2array(resolution08(:,5)),'*r','LineWidth',1)
plot(table2array(resolution09(:,2)),table2array(resolution09(:,5)),'*r','LineWidth',1)
plot(table2array(resolution10(:,2)),table2array(resolution10(:,5)),'*r','LineWidth',1)
plot(table2array(resolution11(:,2)),table2array(resolution11(:,5)),'*r','LineWidth',1)
plot(table2array(resolution12(:,2)),table2array(resolution12(:,5)),'*r','LineWidth',1)
plot(table2array(resolution13(:,2)),table2array(resolution13(:,5)),'*r','LineWidth',1)
plot(table2array(resolution14(:,2)),table2array(resolution14(:,5)),'*r','LineWidth',1)
grid on
grid Minor
xlabel('time, [unitless]')
ylabel('[\delta_{x}/\eta_{x}]; [\delta_{y}/\eta_{y}]; [\delta_{z}/\eta_{z}], [unitless]')
legend('\delta_{x}/\eta_{x}', '\delta_{y}/\eta_{y}', '\delta_{z}/\eta_{z}','Location','SouthEast')