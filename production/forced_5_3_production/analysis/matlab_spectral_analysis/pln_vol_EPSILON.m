%% script to compute spectra, subsequently epsilon using the plane data
% sampled from the DNS
%% NOTE: This code applies scaling to all DNS data using HYFLITS measurements of TKE Dissipation rate and viscosity for 30 km AGL
%% NOTE: Wavenumbers extracted from SAM DNS dataset(s) are defined in angular units [multiplied by 2*PI] -- Accounting for this is necessary
%% NOTE: The synthetic observation/sampling of u', v', and w' fields are scaled and are in physical units referenced to HYFLITS data used for scaling
%% NOTE: The spectral analysis following the scaled data is conducted assuming linear frequency and NOT angular frequency
%% NOTE: However, the DNS spectra for E(k), E11(K3), E22(K3), and E33(K3) are all comuted assuming angular frequency
%% NOTE: Finally, for fair comparison betweeen DNS spectra and synthetic observations, the DNS spectrum derived TKEDR is scaled to HYFLITS measured TKEDR

%% House cleaning
clear
close all
clc

%% General Inputs
diag = 1;                               % switch to turn on program diagnostics
lean = 0;                               % unwanted figures are not plotted (saves time) 
Trec = 6;                               % Time record interval [s]
alp = 1.50;                             % Kolmogorov constant: 1.5 for 3D spectrum; 0.49-0.55 for longitudinal spectrum; 0.65 for transverse spectrum
alpL = 0.55;                            % Kolmogorov constant: 1.5 for 3D spectrum; 0.49-0.55 for longitudinal spectrum; 0.65 for transverse spectrum
cf = alpL*(2*pi)^(-2/3);                % Constant used by Frehlich et al 2003 [if derived in cycles/m]
cf_rad = cf*(2*pi)^(2/3);               % Constant used by Frehlich et al 2003 [if derived in rad/m]
alpT = 0.65;                            % Kolmogorov constant: 1.5 for 3D spectrum; 0.49-0.55 for longitudinal spectrum; 0.65 for transverse spectrum
winOn = 0;                              % switch to turn on or off the windowing
saSp = 1;                               % scitch to turn on or off saving figures and files
edges = -6:0.05:0;                      % histogram edges and bin widths
edgesSlp = -6:0.05:0;                   % histogram edges and bin widths for slopes
ftSz = 16;                              % figure text font size

% Simulation parameters set in sam.inp [ALL THESE QUANTITIES ARE SCALE NORMALIZED]
Xl = 01;                                % Xl = Domain X dimension
Yl = 01;                                % Yl = Domain Y dimension
Zl = 01;                                % Zl = Domain Z dimension
Nx = 2816;                              % Nx = Number of points in Domain X dimension
Ny = 2816;                              % Ny = Number of points in Domain Y dimension
Nz = 2816;                              % Nz = Number of points in Domain Z dimension
nuDNS = 1.5e-5;                         % The kinematic viscosity set in the DNS

% calculate DNS resolution parameters [ALL THESE QUANTITIES ARE SCALE NORMALIZED]
dx = Xl/Nx;                             % Grid Resolution in X
dy = Yl/Ny;                             % Grid Resolution in Y
dz = Zl/Nz;                             % Grid Resolution in Z

% Simulation outputs gathered from .out files [ALL THESE QUANTITIES ARE SCALE NORMALIZED]
resMet = mean([1.4437 1.4426 1.4418]);	% resMet = DNS resolution metric ( = grid spacing/ kolmogorov length scale)
DNSTKE = 7.8742E-01;                    % DNS mean Energy at the sampled time step
epsDNS = 8.2154E-01;                    % DNS mean TKEDR at the sampled time step
etaDNS = (nuDNS^3/epsDNS)^(1/4);        
tauDNS = (nuDNS/epsDNS)^(1/2);
muDNS = etaDNS/tauDNS;

% Scaling parameters applied from measurement [HYFLITS measurement at ~30 km AGL]
epsMeas = 0.98142e-3;                   % epsMeas = TKE dissipation rate [m^3s^-2]
nuMeas = 4.0e-4;                        % nu = Kinematic viscosity [m^2s^-1]
balRt = 002.5;                          % balRt = HYFLITS balloon descent rate [m/s]
etaMEas = (nuMeas^3/epsMeas)^(1/4);     % Kolmogorov length scale [m]
tauMeas = (nuMeas/epsMeas)^(1/2);       % Kolmogorov time scale [s]
muMeas = etaMEas/tauMeas;               % Kolmogorov velocity scale [m/s]

% Scaling ratios: The DNS data should be scaled using the factors below
ltScal = etaMEas/etaDNS;                % change in length scale
tmScl = tauMeas/tauDNS;                 % change in time scale
velScl = ltScal/tmScl;                  % change in velocity scale
epsScal = epsMeas/epsDNS;               % change in epsilon
TDRScal = 1;                            % change in TDR (needs to be formally computed; for now, set to 1)

% Calculate the DNS domain extent using measured Kolmogorov scale and scale separation
decSpac = log10(Nz/3) - log10(1);
Zscal = 10^(log10(etaMEas) + decSpac);	% Scaled Z DNS domain dimension [m]
Xscal = 10^(log10(etaMEas) + decSpac);	% Scaled X DNS domain dimension [m]
Yscal = 10^(log10(etaMEas) + decSpac);	% Scaled Y DNS domain dimension [m]
dxscal = Xscal/Nx;                      % Grid Resolution in X (scaled) [m]
dyscal = Yscal/Ny;                      % Grid Resolution in Y (scaled) [m]
dzscal = Zscal/Nz;                      % Grid Resolution in Z (scaled) [m]

%% calculate the turbulence Reynolds number -- [following eqn (6.7) from Stephen B Pope 2000]
Re = (Zscal/etaMEas)^(4/3);
int_L = DNSTKE^(3/2)/epsDNS;
Re_L = DNSTKE.^2/(epsDNS*nuDNS);
T_ms2 = sqrt(10)*(Re_L^(-0.5))*int_L;
var_ms = (1/sqrt(3))*sqrt((7.9508E-01)^2 + (1.4085E+00)^2 + (6.7310E-01)^2);
T_ms = sqrt(15*nuDNS/epsDNS)*var_ms;

Re_T1 = sqrt((20/3)*Re_L);
Re_T2 = (sqrt(Xl^2 + Yl^2 + Zl^2)/etaDNS)^(4/3);
Re_T3 = ((T_ms/etaDNS)*(1/sqrt(10)))^4;
Re_T4 = var_ms*T_ms2/nuDNS;

%% Load the dataset from file(s)
% for u
dirtry = [pwd '/'];
flPlnU = strcat(dirtry,'yz1_0182_U.txt');
temp = table2array(readtable(flPlnU));
PlnU = transpose(reshape(temp,[Nx+1,Ny+1]));
clear temp
% for v
flPlnV = strcat(dirtry,'yz1_0182_V.txt');
temp = table2array(readtable(flPlnV));
PlnV = transpose(reshape(temp,[Nx+1,Ny+1]));
clear temp
% for W
flPlnW = strcat(dirtry,'yz1_0182_W.txt');
temp = table2array(readtable(flPlnW));
PlnW = transpose(reshape(temp,[Nx+1,Ny+1]));
clear temp
% for T
flPlnT = strcat(dirtry,'yz1_0182_T.txt');
temp = table2array(readtable(flPlnT));
PlnT = transpose(reshape(temp,[Nx+1,Ny+1]));
clear temp
% for DNS TKE dissipation rate
flPlnE = strcat(dirtry,'yz1_0182_E.txt');
temp = table2array(readtable(flPlnE));
PlnE = transpose(reshape(temp,[Nx+1,Ny+1]));
clear temp
% for DNS temperature dissipation rate
flPlnTDR = strcat(dirtry,'yz1_0182_TDR.txt');
temp = table2array(readtable(flPlnTDR));
PlnTDR = transpose(reshape(temp,[Nx+1,Ny+1]));
clear temp
% Reshape the data arrays for the plane
for i = 1:1:length(PlnW(1,1:end-1))
    Uz(:,i) = PlnU(1:end-1,i).*velScl;
    Vz(:,i) = PlnV(1:end-1,i).*velScl;
	Wz(:,i) = PlnW(1:end-1,i).*velScl;
    Tz(:,i) = PlnT(1:end-1,i).*velScl;
	epsDNSz(:,i) = log10(PlnE(1:end-1,i).*epsScal);
    TDRDNSz(:,i) = log10(PlnTDR(1:end-1,i).*TDRScal);
end

%% Compute the number of trajectories and sample points per data record
tmpSz = size(Wz);
numTraj = tmpSz(2);
sampPts = tmpSz(1);
% Compute the number of sample points per data record
numInts = floor(Zscal/balRt);
numSmpls = Nx/numInts;                                          % Also equal to the sampling rate [Hz]

%% Compute frequency averaging and inputs for Spectral analysis
Nrec = Trec*numSmpls;                                           % Number of samples per record
f_low_avg = 1/Trec;                                             % low frequency limit for spectral averaging [Hz]
f_high_avg = numSmpls/2;                                        % high frequency limit for spectral averaging [Hz]
pts_in_dec = 1/3;                                               % [1/3 = 1/3rd decade averaging of the spectra] [set 1/10 for 1/10th decade averaging of the spectra]
%pit_NF = 10^(-4+log10(alpL)+(2/3)*log10(balRt));                % define the noise floor for pitot 
pit_NF = 10^(-20);                                              % define the noise floor for pitot 
min_fit_pts = 2;                                                % the minimum number of points to be used in the third pass of fitting
pass_2 = 1;                                                     % This switch turns on an additional level scan on the spectra -- provides conservative turbulence estimates
time_segment_center_inds = (Nrec/2):Nrec:length(Wz);
dt = 1/numSmpls;                                                % time resolution [sec]
f_low = f_low_avg;                                              % low frequency limit for spectral averaging [Hz]
f_high = f_high_avg;                                            % high frequency limit for spectral averaging [Hz]
freq = [1:Nrec/2]/(Nrec*dt);                                    % spectral frequency samples, up to Nyquist rate [Hz]
freqR = freq;                                                   % store the frequency for plotting
f_ind_low = find(freq >= f_low,1,'first');
f_ind_high = find(freq <= f_high,1,'last');
freqs_used = freq(f_ind_low:f_ind_high);
% Setup the frequency boundaries for fractional decade averaging 
decades = log10(f_high_avg)-log10(f_low_avg);
f_avg = logspace(log10(f_low_avg),log10(f_high_avg),floor(decades/pts_in_dec));
f_nth_dec = (f_avg(2:end)+f_avg(1:end-1))/2;                    % bin center frequencies
f_inds = zeros(1,length(f_avg));
for i = 1:length(f_avg)
    f_inds(i) = find(freq-f_avg(i) >= -1e-6,1,'First');         % index of first frequency in each bin
end
f_avgind_low = freq(f_inds(1:end-1));
f_avgind_high = freq(f_inds(2:end)-1);
freqs_used_avg = f_nth_dec;

%% Compute PSD and TKEDR for u' v' and w'
compSpecU;
compSpecV;
compSpecW;

% compute the average of [epsilon_spec_U epsilon_spec_V epsilon_spec_W]
for i = 1:1:numTraj
    for j = 1:1:(numInts/Trec)
        eps_avgUVW(j,i) = mean([log_Uepsilon_k(j,i) log_Vepsilon_k(j,i) log_Wepsilon_k(j,i)]);
    end
end

% also compute the standard deviations of the distributions of various estimates of TKEDR
epsDNSzSD = std(10.^(epsDNSz),0,'all','omitnan');
log_Wepsilon_kSD = std(10.^(log_Wepsilon_k),0,'all','omitnan');
mnepsDNSzSD = std(10.^(mnepsDNSz),0,'all','omitnan');

% Plot Histograms comparing spectrally derived TKEDR and averaged DNS pointwise extimates of TKEDR
figure(401)
clf
histogram(epsDNSz,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'LineWidth',2)
hold on
histogram(log_Wepsilon_k,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','b','LineWidth',2)
histogram(mnepsDNSz,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','r','LineWidth',2)
histogram(eps_avgUVW,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
plot([mean(mean(log_Wepsilon_k)) mean(mean(log_Wepsilon_k))],[0 0.25],'b','LineWidth',2)
plot([mean(mean(mnepsDNSz)) mean(mean(mnepsDNSz))],[0 0.25],'r','LineWidth',2)
plot([mean(mean(eps_avgUVW)) mean(mean(eps_avgUVW))],[0 0.25],'k','LineWidth',2)
xlim([edges(1)+1 edges(end)-1])
ylim([0 0.2])
xlabel('log10 \epsilon')
ylabel('probability')
legend('\epsilon_{DNS} - samples on YZ pln','\epsilon_{spec}','<\epsilon_{DNS}>','<\epsilon_{UVW}>','Location','NorthWest','FontSize',ftSz)
text(-2.5,0.175,['\sigma_{\epsilon_{DNS}} = ' num2str(epsDNSzSD,'%2.2e')],'FontSize',ftSz-2)
text(-2.5,0.160,['\sigma_{\epsilon_{spec}} = ' num2str(log_Wepsilon_kSD,'%2.2e')],'FontSize',ftSz-2)
text(-2.5,0.145,['\sigma_{<\epsilon_{DNS}>} = ' num2str(mnepsDNSzSD,'%2.2e')],'FontSize',ftSz-2)
text(-2.5,0.130,['\sigma_{\epsilon_{DNS}}/\sigma_{<\epsilon_{DNS}>} = ' num2str(epsDNSzSD/mnepsDNSzSD,'%2.2f')],'FontSize',ftSz-2)
grid on
set(gca,'FontSize',ftSz)
if saSp == 1
	savefig('prob_spec_DNS.fig')
    close all
end

figure(402)
clf
histogram(log_Uepsilon_k,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','r','LineWidth',2)
hold on
histogram(log_Vepsilon_k,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','g','LineWidth',2)
histogram(log_Wepsilon_k,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','b','LineWidth',2)
plot([mean(mean(log_Uepsilon_k)) mean(mean(log_Uepsilon_k))],[0 0.25],'r','LineWidth',2)
plot([mean(mean(log_Vepsilon_k)) mean(mean(log_Vepsilon_k))],[0 0.25],'g','LineWidth',2)
plot([mean(mean(log_Wepsilon_k)) mean(mean(log_Wepsilon_k))],[0 0.25],'b','LineWidth',2)
xlim([edges(1)+1 edges(end)-1])
ylim([0 0.2])
xlabel('log10 \epsilon')
ylabel('probability')
legend('\epsilon_{specU}','\epsilon_{specV}','\epsilon_{specW}','Location','NorthEast','FontSize',ftSz)
grid on
set(gca,'FontSize',ftSz)
if saSp == 1
	savefig('prob_spec_UVW.fig')
    close all
end

figure(403)
clf
histogram(log_Wepsilon_k,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','b','LineWidth',2)
hold on
histogram(mnepsDNSz,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','r','LineWidth',2)
histogram(log_Wepsilon_u,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'LineWidth',2)
histogram(log_Wepsilon_l,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'LineWidth',2)
plot([mean(mean(log_Wepsilon_k)) mean(mean(log_Wepsilon_k))],[0 0.25],'b','LineWidth',2)
plot([mean(mean(mnepsDNSz)) mean(mean(mnepsDNSz))],[0 0.25],'r','LineWidth',2)
xlim([edges(1)+1 edges(end)-1])
ylim([0 0.2])
xlabel('log10 \epsilon')
ylabel('probability')
legend('\epsilon_{spec}','\epsilon_{DNS}','\epsilon_{spec} - std. dev.','Location','NorthWest','FontSize',ftSz)
grid on
set(gca,'FontSize',ftSz)
if saSp == 1
	savefig('prob_specUNC_DNS.fig')
    close all
end

figure(4030)
clf
histogram(mnepsDNSz,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','r','LineWidth',2)
hold on
histogram(stdepsDNSzU,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'LineWidth',2)
histogram(stdepsDNSzL,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0.5 0.5 0.5],'LineWidth',2)
xlim([edges(1)+1 edges(end)-1])
ylim([0 0.2])
xlabel('log10 \epsilon')
ylabel('probability')
legend('\epsilon_{DNS}','range \epsilon_{DNS}','Location','NorthWest','FontSize',ftSz)
grid on
set(gca,'FontSize',ftSz)
if saSp == 1
	savefig('prob_DNS_rng.fig')
    close all
end

figure(404)
clf
histogram(log_Wepsilon_k,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','b','LineWidth',2)
hold on
histogram(mnepsDNSz,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','r','LineWidth',2)
histogram(log_Wepsilon_k2,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','m','LineWidth',2)
plot([mean(mean(log_Wepsilon_k)) mean(mean(log_Wepsilon_k))],[0 0.25],'b','LineWidth',2)
plot([mean(mean(log_Wepsilon_k2)) mean(mean(log_Wepsilon_k2))],[0 0.25],'m','LineWidth',2)
plot([mean(mean(mnepsDNSz)) mean(mean(mnepsDNSz))],[0 0.25],'r','LineWidth',2)
xlim([edges(1)+1 edges(end)-1])
ylim([0 0.2])
xlabel('log10 \epsilon')
ylabel('probability')
legend('\epsilon_{spec} - Theory','\epsilon_{DNS}','\epsilon_{spec} - Frehlich 2003','Location','NorthWest','FontSize',ftSz)
grid on
set(gca,'FontSize',ftSz)
if saSp == 1
	savefig('prob_specMODS_DNS.fig')
    close all
end

for pt = 1:1:(numInts/Trec)
    [rU,cU] = find(WTrlog_pwp_k_avg == max(WTrlog_pwp_k_avg(:,pt)));
    [rL,cL] = find(WTrlog_pwp_k_avg == min(WTrlog_pwp_k_avg(:,pt)));
    figure(405 + pt)
    clf
    semilogx(freq,log10(abs(WppdfPlt(:,rU,cU))),'r','LineWidth',0.5)
    hold on
    semilogx(freq,log10(abs(WppdfPlt(:,rL,cL))),'b','LineWidth',0.5)
    semilogx(f_nth_dec,avWlog_ppsd_avg(:,pt),'g','LineWidth',3)
    semilogx([freq(1) freq(end)],[log10(pit_NF) log10(pit_NF)],'m--','LineWidth',2)
    semilogx(f_nth_dec(avWk_inds_pit{pt}),avWlog_ppsd_avg(avWk_inds_pit{pt},pt),'ko','LineWidth',2)
    semilogx(freq,avWlog_pwp_k_avg(pt)+log10(freq.^(-5/3)),'k','LineWidth',2)
    semilogx(freq,(avWlog_pwp_k_avg(pt)+avWlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
    semilogx(freq,(avWlog_pwp_k_avg(pt)-avWlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
    axis([freq(1) freq(end) -15 2])
    xlabel('\kappa')
    ylabel('Wz PSD')
    grid on
    legend('raw spec. for \epsilon_{max}','raw spec. for \epsilon_{min}','avg. of all spec.',...
        'HW noise floor - 30 km [AGL]','qual. data pts.','f^{-5/3} fit','fit std. dev.',...
        'Location','SouthWest','FontSize',ftSz)
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig(['specBnds_Wz_' num2str(pt) '.fig'])
        close all
    end
end

[rU,cU] = find(WTrlog_pwp_k_avg == max(max(WTrlog_pwp_k_avg)));
[rL,cL] = find(WTrlog_pwp_k_avg == min(min(WTrlog_pwp_k_avg)));
figure(406)
clf
semilogx(freq,log10(abs(WppdfPlt(:,rU,cU))),'b','LineWidth',0.5)
hold on
semilogx(freq,log10(abs(WppdfPlt(:,rL,cL))),'b','LineWidth',0.5)
semilogx(f_nth_dec,mean(avWlog_ppsd_avg,2),'r','LineWidth',3)
axis([0.1 300 -15 2])
xlabel('\kappa')
ylabel('Wz PSD')
set(gca,'FontSize',ftSz)
if saSp == 1
    savefig(['specBnds_Wz_all.fig'])
    close all
end

% Plot averaged spectra from each data record height in the domain
if diag == 1 && lean == 1
    for pt = 1:1:(numInts/Trec)
        figure(500+pt)
        clf
        semilogx(freq,log10(abs(avWppsd(:,pt))),'b','LineWidth',2)
        hold on
        semilogx(freq,avWlog_pwp_k_avg(pt)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,avWlog_ppsd_avg(:,pt),'r*','LineWidth',4)
        semilogx(f_nth_dec(avWk_inds_pit{pt}),avWlog_ppsd_avg(avWk_inds_pit{pt},pt),'go','LineWidth',2)
        semilogx(freq,(avWlog_pwp_k_avg(pt)+avWlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(avWlog_pwp_k_avg(pt)-avWlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(f_nth_dec,avWlog_ppsd_avg(:,pt),'r','LineWidth',2)
        semilogx([freq(1) freq(end)],[log10(pit_NF) log10(pit_NF)],'m--','LineWidth',2)
        legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting',...
            'error bars','Location','NorthEast','FontSize',ftSz)
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa')
        ylabel('averaged Wz PSD')
        title('normalized wavenumber vs averaged Wz Spectra sampled on YZ plane (at X/2)')
        grid on
        set(gca,'FontSize',ftSz)
        if saSp == 1
            savefig(['spec_Wz_avg_',num2str(pt),'.fig'])
            close all
        end
        
        figure(600+pt)
        clf
        semilogx(freq,log10(abs(avUppsd(:,pt))),'b','LineWidth',2)
        hold on
        semilogx(freq,avUlog_pwp_k_avg(pt)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,avUlog_ppsd_avg(:,pt),'r*','LineWidth',4)
        semilogx(f_nth_dec(avUk_inds_pit{pt}),avUlog_ppsd_avg(avUk_inds_pit{pt},pt),'go','LineWidth',2)
        semilogx(freq,(avUlog_pwp_k_avg(pt)+avUlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(avUlog_pwp_k_avg(pt)-avUlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(f_nth_dec,avUlog_ppsd_avg(:,pt),'r','LineWidth',2)
        semilogx([freq(1) freq(end)],[log10(pit_NF) log10(pit_NF)],'m--','LineWidth',2)
        legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting',...
            'error bars','Location','NorthEast','FontSize',ftSz)
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa')
        ylabel('averaged Uz PSD')
        title('normalized wavenumber vs averaged Uz Spectra sampled on YZ plane (at X/2)')
        grid on
        set(gca,'FontSize',ftSz)
        if saSp == 1
            savefig(['spec_Uz_avg_',num2str(pt),'.fig'])
            close all
        end
        
        figure(700+pt)
        clf
        semilogx(freq,log10(abs(avVppsd(:,pt))),'b','LineWidth',2)
        hold on
        semilogx(freq,avVlog_pwp_k_avg(pt)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,avVlog_ppsd_avg(:,pt),'r*','LineWidth',4)
        semilogx(f_nth_dec(avVk_inds_pit{pt}),avVlog_ppsd_avg(avVk_inds_pit{pt},pt),'go','LineWidth',2)
        semilogx(freq,(avVlog_pwp_k_avg(pt)+avVlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(avVlog_pwp_k_avg(pt)-avVlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(f_nth_dec,avVlog_ppsd_avg(:,pt),'r','LineWidth',2)
        semilogx([freq(1) freq(end)],[log10(pit_NF) log10(pit_NF)],'m--','LineWidth',2)
        legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting',...
            'error bars','Location','NorthEast','FontSize',ftSz)
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa')
        ylabel('averaged Vz PSD')
        title('normalized wavenumber vs averaged Vz Spectra sampled on YZ plane (at X/2)')
        grid on
        set(gca,'FontSize',ftSz)
        if saSp == 1
            savefig(['spec_Vz_avg_',num2str(pt),'.fig'])
            close all
        end
        
        figure(800+pt)
        clf
        plot(mnepsDNSz(pt,:),'r')
        hold on
        plot(log_Wepsilon_k(pt,:),'b')
        plot(eps_avgUVW(pt,:),'k')
        plot([1 length(mnepsDNSz(pt,:))],[avTrepsDNSz(pt) avTrepsDNSz(pt)],'r','LineWidth',2)
        plot([1 length(log_Wepsilon_k(pt,:))],[log_avWepsilon_k(pt) log_avWepsilon_k(pt)],'b','LineWidth',2)
        legend('\epsilon_{DNS}','\epsilon_{specW}','<\epsilon_{UVW}>','\epsilon_{DNS} averaged','\epsilon_{DNS} averaged','Location','NorthEast','FontSize',ftSz)
        xlabel('Y')
        ylabel('\epsilon')
        grid on
        set(gca,'FontSize',ftSz)
        if saSp == 1
            savefig(['eps_DNS_Wspec_level_',num2str(pt),'.fig'])
            close all
        end
        
        figure(900+pt)
        clf
        plot(mnepsDNSz(pt,:),'r')
        hold on
        plot(log_Uepsilon_k(pt,:),'b')
        plot(eps_avgUVW(pt,:),'k')
        plot([1 length(mnepsDNSz(pt,:))],[avTrepsDNSz(pt) avTrepsDNSz(pt)],'r','LineWidth',2)
        plot([1 length(log_Uepsilon_k(pt,:))],[log_avUepsilon_k(pt) log_avUepsilon_k(pt)],'b','LineWidth',2)
        legend('\epsilon_{DNS}','\epsilon_{specU}','<\epsilon_{UVW}>','\epsilon_{DNS} averaged','\epsilon_{DNS} averaged','Location','NorthEast','FontSize',ftSz)
        xlabel('Y')
        ylabel('\epsilon')
        grid on
        set(gca,'FontSize',ftSz)
        if saSp == 1
            savefig(['eps_DNS_Uspec_level_',num2str(pt),'.fig'])
            close all
        end
        
        figure(1000+pt)
        clf
        plot(mnepsDNSz(pt,:),'r')
        hold on
        plot(log_Vepsilon_k(pt,:),'b')
        plot(eps_avgUVW(pt,:),'k')
        plot([1 length(mnepsDNSz(pt,:))],[avTrepsDNSz(pt) avTrepsDNSz(pt)],'r','LineWidth',2)
        plot([1 length(log_Vepsilon_k(pt,:))],[log_avVepsilon_k(pt) log_avVepsilon_k(pt)],'b','LineWidth',2)
        legend('\epsilon_{DNS}','\epsilon_{specV}','<\epsilon_{UVW}>','\epsilon_{DNS} averaged','\epsilon_{DNS} averaged','Location','NorthEast','FontSize',ftSz)
        xlabel('Y')
        ylabel('\epsilon')
        grid on
        set(gca,'FontSize',ftSz)
        if saSp == 1
            savefig(['eps_DNS_Vspec_level_',num2str(pt),'.fig'])
            close all
        end
    end
end

% Plot all individual specta
%{
if diag == 1 && lean == 0
    parfor i = 1:1:numTraj
        j = 1;
        figure(1100+i)
        clf
        subplot(1,2,1)
        plot(Wz(:,i))
        hold on
        plot(detrend(Wz(:,i)))
        plot(Wz(:,i)-mean(Wz(:,i)))
        grid on
        grid Minor
        ylabel('Wz')
        legend('Wz','Wz - detrended','Wz - mean removed','FontSize',ftSz)
        ylim([-3 3])
        subplot(1,2,2)
        semilogx(freq,log10(abs(WppdfPlt(:,j,i))),'b','LineWidth',2)
        hold on
        semilogx(freq,WTrlog_pwp_k_avg(j,i)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,WTrlog_ppsd_avg(:,j,i),'r*','LineWidth',4)
        semilogx(f_nth_dec(Wusd_inds{j,i}),WTrlog_ppsd_avg(Wusd_inds{j,i},j,i),'go','LineWidth',2)
        semilogx(freq,(WTrlog_pwp_k_avg(j,i)+WTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(WTrlog_pwp_k_avg(j,i)-WTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)        
        legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','NorthEast','FontSize',ftSz)
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa')
        ylabel('Wz PSD')
        title(['normalized wavenumber vs Wz Spectra, Trajectory = ',num2str(i)])
        grid on
        set(gca,'FontSize',ftSz)
        if saSp == 1
            savefig(['Wz_Tr',num2str(i),'.fig'])
            close all
        end
        
        figure(1200+i)
        clf
        subplot(1,2,1)
        plot(Uz(:,i))
        hold on
        plot(detrend(Uz(:,i)))
        plot(Uz(:,i)-mean(Uz(:,i)))
        grid on
        grid Minor
        ylabel('Uz')
        legend('Uz','Uz - detrended','Uz - mean removed','FontSize',ftSz)
        ylim([-3 3])
        subplot(1,2,2)
        semilogx(freq,log10(abs(UppdfPlt(:,j,i))),'b','LineWidth',2)
        hold on
        semilogx(freq,UTrlog_pwp_k_avg(j,i)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,UTrlog_ppsd_avg(:,j,i),'r*','LineWidth',4)
        semilogx(f_nth_dec(Uusd_inds{j,i}),UTrlog_ppsd_avg(Uusd_inds{j,i},j,i),'go','LineWidth',2)
        semilogx(freq,(UTrlog_pwp_k_avg(j,i)+UTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(UTrlog_pwp_k_avg(j,i)-UTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)        
        legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','NorthEast','FontSize',ftSz)
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa')
        ylabel('Uz PSD')
        title(['normalized wavenumber vs Uz Spectra, Trajectory = ',num2str(i)])
        grid on
        set(gca,'FontSize',ftSz)
        if saSp == 1
            savefig(['Uz_Tr',num2str(i),'.fig'])
            close all
        end
        
        figure(1300+i)
        clf
        subplot(1,2,1)
        plot(Vz(:,i))
        hold on
        plot(detrend(Vz(:,i)))
        plot(Vz(:,i)-mean(Vz(:,i)))
        grid on
        grid Minor
        ylabel('Vz')
        legend('Vz','Vz - detrended','Vz - mean removed','FontSize',ftSz)
        ylim([-3 3])
        subplot(1,2,2)
        semilogx(freq,log10(abs(VppdfPlt(:,j,i))),'b','LineWidth',2)
        hold on
        semilogx(freq,VTrlog_pwp_k_avg(j,i)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,VTrlog_ppsd_avg(:,j,i),'r*','LineWidth',4)
        semilogx(f_nth_dec(Vusd_inds{j,i}),VTrlog_ppsd_avg(Vusd_inds{j,i},j,i),'go','LineWidth',2)
        semilogx(freq,(VTrlog_pwp_k_avg(j,i)+VTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(VTrlog_pwp_k_avg(j,i)-VTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)        
        legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','NorthEast','FontSize',ftSz)
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa')
        ylabel('Vz PSD')
        title(['normalized wavenumber vs Vz Spectra, Trajectory = ',num2str(i)])
        grid on
        set(gca,'FontSize',ftSz)
        if saSp == 1
            savefig(['Vz_Tr',num2str(i),'.fig'])
            close all
        end
    end
end
%}

%% Compute Epsilon from Nx*Ny spectra of Wz
% Load the dataset from file(s)
load specUUz.mat
specUUz = specUUz.*(Zl/(2*pi));
load specVVz.mat
specVVz = specVVz.*(Zl/(2*pi));
load specWWz.mat
specWWz = specWWz.*(Zl/(2*pi));
load specEk.mat
specEk = 0.5*specEk(2:end-1).*(Zl/(2*pi));
planeSpecWWz;
planeSpecUUz;
planeSpecVVz;
planeSpecEk;

% plot the XZ all averaged PSD and fit
if diag == 1 && lean == 1
    figure(1105)
    clf
    loglog(freqE.*(etaDNS),volEkw.*(epsDNS)^(-2/3),'k','LineWidth',2)
    hold on
    loglog(freqU.*(etaDNS),volUzw.*(epsDNS)^(-2/3),'r','LineWidth',2)
    loglog(freqV.*(etaDNS),volVzw.*(epsDNS)^(-2/3),'m','LineWidth',2)
    loglog(freqW.*(etaDNS),volWzw.*(epsDNS)^(-2/3),'b','LineWidth',2)
    loglog([freqE(1).*(etaDNS) freqE(end).*(etaDNS)],[alp alp],'--k')
    loglog([freqE(1).*(etaDNS) freqE(end).*(etaDNS)],[alpT alpT],'--r')
    loglog([freqE(1).*(etaDNS) freqE(end).*(etaDNS)],[alpT alpT],'--m')
    loglog([freqE(1).*(etaDNS) freqE(end).*(etaDNS)],[alpL alpL],'--b')
    lgd = legend('{\epsilon}^{-2/3}\kappa^{5/3}E(\kappa)',...
        '{\epsilon}^{-2/3}{\kappa_{3}}^{5/3}E_{11}(\kappa_{3})',...
        '{\epsilon}^{-2/3}{\kappa_{3}}^{5/3}E_{22}(\kappa_{3})',...
        '{\epsilon}^{-2/3}{\kappa_{3}}^{5/3}E_{33}(\kappa_{3})',...
        'Location','SouthWest');
    lgd.FontSize = ftSz;
    xlabel('\kappa\eta')
    ylabel('{\epsilon}^{-2/3}\kappa^{5/3}E(\kappa)')
    title('Compensated spectra to verify universal constants')
    grid on
    grid Minor
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig('compenSpec_compare.fig')
        close all
    end
end

% Scale the SAM DNS averaged spectra derived TKEDR to compare with synthetic observation derived TKEDR
log_volepsilon_kW = log_volepsilon_kW+log10(epsScal);
log_volepsilon_kU = log_volepsilon_kU+log10(epsScal);
log_volepsilon_kV = log_volepsilon_kV+log10(epsScal);
log_volepsilon_kE = log_volepsilon_kE+log10(epsScal);
% Convert to log scale -- for plotting
logDNSEPSMN = log10(epsDNS.*epsScal);
PlnE_mean = mean(mean(epsDNSz));                                     % avg DNS on the sampled XZ plane

% Plot comparing individual sample Epsilon calculated using Wz profiles and
% the average epsilon for the plane
for pt = 1:1:(numInts/Trec)
    figure(2100+pt)
    clf
    plot(log_Wepsilon_k(pt,:),'g','LineWidth',2)        % spec - records
    hold on
    plot(mnepsDNSz(pt,:),'r','LineWidth',2)             % DNS - records
    plot([1 length(log_Wepsilon_k)],[log_avWepsilon_k(pt) log_avWepsilon_k(pt)],'g','LineWidth',2)  % avg spec pln
    plot([1 length(log_Wepsilon_k)],[PlnE_mean PlnE_mean],'r','LineWidth',2)                        % avg DNS pln
    plot([1 length(log_Wepsilon_k)],[logDNSEPSMN logDNSEPSMN],'k','LineWidth',2)                    % avg DNS full domain
    plot([1 length(log_Wepsilon_k)],[log_volepsilon_kW log_volepsilon_kW],'--b','LineWidth',2)      % avg spec all records for each XY -- WWz
    plot([1 length(log_Wepsilon_k)],[log_volepsilon_kU log_volepsilon_kU],'--r','LineWidth',2)      % avg spec all records for each XY -- UUz
    plot([1 length(log_Wepsilon_k)],[log_volepsilon_kV log_volepsilon_kV],'--m','LineWidth',2)      % avg spec all records for each XY -- VVz
    plot([1 length(log_Wepsilon_k)],[log_volepsilon_kE log_volepsilon_kE],'--k','LineWidth',2)      % avg spec all records for each XY -- Ek
    ylim([-6 1])
    grid on
    grid Minor
    xlabel('Profile number on the XZ plane')
    ylabel('log10 \epsilon_{Wz}')
    h = legend('$$\epsilon_{specW}$$ records','$$\epsilon_{DNS}$$ records','$$\overline{\epsilon}_{specW}$$ YZ pln avg','$$\overline{\epsilon}_{DNS}$$ YZ pln avg','$$<\epsilon> from all DNS grid points$$','$$<\epsilon> from E_{33}(\kappa_{3})$$','$$<\epsilon> from E_{11}(\kappa_{3})$$','$$<\epsilon> from E_{22}(\kappa_{3})$$','Location','NorthEast');
    title('Comparing \epsilon_{Wz} individual and average PSD for records on fixed YZ plane with domain DNS averaged \epsilon')
    set(h,'Interpreter','latex','FontSize',ftSz)
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig(['epsilon_Wz_pln_compare_',num2str(pt),'.fig'])
        close all
    end
    figure(2200+pt)
    clf
    plot(log_Uepsilon_k(pt,:),'g','LineWidth',2)        % spec - records
    hold on
    plot(mnepsDNSz(pt,:),'r','LineWidth',2)             % DNS - records
    plot([1 length(log_Uepsilon_k)],[log_avUepsilon_k(pt) log_avUepsilon_k(pt)],'g','LineWidth',2)  % avg spec pln
    plot([1 length(log_Uepsilon_k)],[PlnE_mean PlnE_mean],'r','LineWidth',2)                        % avg DNS pln
    plot([1 length(log_Uepsilon_k)],[logDNSEPSMN logDNSEPSMN],'k','LineWidth',2)                    % avg DNS full domain
    plot([1 length(log_Uepsilon_k)],[log_volepsilon_kW log_volepsilon_kW],'--b','LineWidth',2)      % avg spec all records for each XY -- WWz
    plot([1 length(log_Uepsilon_k)],[log_volepsilon_kU log_volepsilon_kU],'--r','LineWidth',2)      % avg spec all records for each XY -- UUz
    plot([1 length(log_Uepsilon_k)],[log_volepsilon_kV log_volepsilon_kV],'--m','LineWidth',2)      % avg spec all records for each XY -- VVz
    plot([1 length(log_Uepsilon_k)],[log_volepsilon_kE log_volepsilon_kE],'--k','LineWidth',2)      % avg spec all records for each XY -- Ek
    ylim([-6 1])
    grid on
    grid Minor
    xlabel('Profile number on the XZ plane')
    ylabel('log10 \epsilon_{Uz}')
    h = legend('$$\epsilon_{specU}$$ records','$$\epsilon_{DNS}$$ records','$$\overline{\epsilon}_{specU}$$ YZ pln avg','$$\overline{\epsilon}_{DNS}$$ YZ pln avg','$$<\epsilon> from all DNS grid points$$','$$<\epsilon> from E_{33}(\kappa_{3})$$','$$<\epsilon> from E_{11}(\kappa_{3})$$','$$<\epsilon> from E_{22}(\kappa_{3})$$','Location','NorthEast');
    title('Comparing \epsilon_{Uz} individual and average PSD for records on fixed YZ plane with domain DNS averaged \epsilon')
    set(h,'Interpreter','latex','FontSize',ftSz)
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig(['epsilon_Uz_pln_compare_',num2str(pt),'.fig'])
        close all
    end
    figure(2300+pt)
    clf
    plot(log_Vepsilon_k(pt,:),'g','LineWidth',2)        % spec - records
    hold on
    plot(mnepsDNSz(pt,:),'r','LineWidth',2)             % DNS - records
    plot([1 length(log_Vepsilon_k)],[log_avVepsilon_k(pt) log_avVepsilon_k(pt)],'g','LineWidth',2)  % avg spec pln
    plot([1 length(log_Vepsilon_k)],[PlnE_mean PlnE_mean],'r','LineWidth',2)                        % avg DNS pln
    plot([1 length(log_Vepsilon_k)],[logDNSEPSMN logDNSEPSMN],'k','LineWidth',2)                    % avg DNS full domain
    plot([1 length(log_Vepsilon_k)],[log_volepsilon_kW log_volepsilon_kW],'--b','LineWidth',2)      % avg spec all records for each XY -- WWz
    plot([1 length(log_Vepsilon_k)],[log_volepsilon_kU log_volepsilon_kU],'--r','LineWidth',2)      % avg spec all records for each XY -- UUz
    plot([1 length(log_Vepsilon_k)],[log_volepsilon_kV log_volepsilon_kV],'--m','LineWidth',2)      % avg spec all records for each XY -- VVz
    plot([1 length(log_Vepsilon_k)],[log_volepsilon_kE log_volepsilon_kE],'--k','LineWidth',2)      % avg spec all records for each XY -- Ek
    ylim([-6 1])
    grid on
    grid Minor
    xlabel('Profile number on the XZ plane')
    ylabel('log10 \epsilon_{Vz}')
    h = legend('$$\epsilon_{specV}$$ records','$$\epsilon_{DNS}$$ records','$$\overline{\epsilon}_{specV}$$ YZ pln avg','$$\overline{\epsilon}_{DNS}$$ YZ pln avg','$$<\epsilon> from all DNS grid points$$','$$<\epsilon> from E_{33}(\kappa_{3})$$','$$<\epsilon> from E_{11}(\kappa_{3})$$','$$<\epsilon> from E_{22}(\kappa_{3})$$','Location','NorthEast');
    title('Comparing \epsilon_{Vz} individual and average PSD for records on fixed YZ plane with domain DNS averaged \epsilon')
    set(h,'Interpreter','latex','FontSize',ftSz)
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig(['epsilon_Vz_pln_compare_',num2str(pt),'.fig'])
        close all
    end
end
% Plot comparing histogram of spectral and DNS derived TKEDR
for pt = 1:1:(numInts/Trec)
    figure(2400+pt)
    clf
    histogram(log_Wepsilon_k(pt,:),edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','b','LineWidth',2)
    hold on
    histogram(mnepsDNSz(pt,:),edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','r','LineWidth',2)
    plot([log_avWepsilon_k(pt) log_avWepsilon_k(pt)],[0 0.5],'g','LineWidth',2)
    plot([PlnE_mean PlnE_mean],[0 0.5],'r','LineWidth',2)
    plot([log_volepsilon_kW log_volepsilon_kW],[0 0.5],'m','LineWidth',2)
    plot([logDNSEPSMN logDNSEPSMN],[0 0.5],'k','LineWidth',2)
    xlim([edges(1)+1 edges(end)-1])
    ylim([0 0.2])
    grid on
    grid Minor
    xlabel('log10 \epsilon_{Wz}')
    ylabel('probability')
    h = legend('$$\epsilon_{specW}$$ records',...
        '$$\epsilon_{DNS}$$ averaged over records','$$\overline{\epsilon}_{specW}$$ YZ pln avg',...
        '$$\overline{\epsilon}_{DNS}$$ YZ pln avg','$$<\epsilon_{specW}>$$ all YZ records',...
        '$$<\epsilon_{DNS}>$$ full domain avg','Location','NorthWest');
    title('Comparing histograms of \epsilon_{specW} and \epsilon_{DNS} for all sampled records on the YZ plane')
    set(h,'Interpreter','latex','FontSize',ftSz)
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig(['hist_Wz_',num2str(pt),'.fig'])
        close all
    end
    figure(2500+pt)
    clf
    histogram(log_Uepsilon_k(pt,:),edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','b','LineWidth',2)
    hold on
    histogram(mnepsDNSz(pt,:),edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','r','LineWidth',2)
    plot([log_avUepsilon_k(pt) log_avUepsilon_k(pt)],[0 0.5],'g','LineWidth',2)
    plot([PlnE_mean PlnE_mean],[0 0.5],'r','LineWidth',2)
    plot([log_volepsilon_kU log_volepsilon_kU],[0 0.5],'m','LineWidth',2)
    plot([logDNSEPSMN logDNSEPSMN],[0 0.5],'k','LineWidth',2)
    xlim([edges(1)+1 edges(end)-1])
    ylim([0 0.2])
    grid on
    grid Minor
    xlabel('log10 \epsilon_{Uz}')
    ylabel('probability')
    h = legend('$$\epsilon_{specU}$$ records',...
        '$$\epsilon_{DNS}$$ averaged over records','$$\overline{\epsilon}_{specU}$$ YZ pln avg',...
        '$$\overline{\epsilon}_{DNS}$$ YZ pln avg','$$<\epsilon_{specU}>$$ all YZ records',...
        '$$<\epsilon_{DNS}>$$ full domain avg','Location','NorthWest');
    title('Comparing histograms of \epsilon_{specU} and \epsilon_{DNS} for all sampled records on the YZ plane')
    set(h,'Interpreter','latex','FontSize',ftSz)
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig(['hist_Uz_',num2str(pt),'.fig'])
        close all
    end
    figure(2600+pt)
    clf
    histogram(log_Vepsilon_k(pt,:),edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','b','LineWidth',2)
    hold on
    histogram(mnepsDNSz(pt,:),edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','r','LineWidth',2)
    plot([log_avVepsilon_k(pt) log_avVepsilon_k(pt)],[0 0.5],'g','LineWidth',2)
    plot([PlnE_mean PlnE_mean],[0 0.5],'r','LineWidth',2)
    plot([log_volepsilon_kV log_volepsilon_kV],[0 0.5],'m','LineWidth',2)
    plot([logDNSEPSMN logDNSEPSMN],[0 0.5],'k','LineWidth',2)
    xlim([edges(1)+1 edges(end)-1])
    ylim([0 0.2])
    grid on
    grid Minor
    xlabel('log10 \epsilon_{Vz}')
    ylabel('probability')
    h = legend('$$\epsilon_{specV}$$ records',...
        '$$\epsilon_{DNS}$$ averaged over records','$$\overline{\epsilon}_{specV}$$ YZ pln avg',...
        '$$\overline{\epsilon}_{DNS}$$ YZ pln avg','$$<\epsilon_{specW}>$$ all YZ records',...
        '$$<\epsilon_{DNS}>$$ full domain avg','Location','NorthWest');
    title('Comparing histograms of \epsilon_{specV} and \epsilon_{DNS} for all sampled records on the YZ plane')
    set(h,'Interpreter','latex','FontSize',ftSz)
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig(['hist_Vz_',num2str(pt),'.fig'])
        close all
    end
    
    figure(2700+pt)
    clf
    plot(mnepsDNSz(pt,:),log_Wepsilon_k(pt,:),'*k')
    hold on
    plot([edges(1) edges(end)],[edges(1) edges(end)],'k','LineWidth',2)
    xlim([edges(1) edges(end)])
    ylim([edges(1) edges(end)])
    grid on
    grid Minor
    xlabel('log10 \epsilon_{DNS}')
    ylabel('log10 \epsilon_{Wz}')
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig(['slope_DNS_Wspec_',num2str(pt),'.fig'])
        close all
    end
    
    figure(2800+pt)
    clf
    plot(mnepsDNSz(pt,:),log_Uepsilon_k(pt,:),'*k')
    hold on
    plot([edges(1) edges(end)],[edges(1) edges(end)],'k','LineWidth',2)
    xlim([edges(1) edges(end)])
    ylim([edges(1) edges(end)])
    grid on
    grid Minor
    xlabel('log10 \epsilon_{DNS}')
    ylabel('log10 \epsilon_{Uz}')
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig(['slope_DNS_Uspec_',num2str(pt),'.fig'])
        close all
    end
    
    figure(2900+pt)
    clf
    plot(mnepsDNSz(pt,:),log_Vepsilon_k(pt,:),'*k')
    hold on
    plot([edges(1) edges(end)],[edges(1) edges(end)],'k','LineWidth',2)
    xlim([edges(1) edges(end)])
    ylim([edges(1) edges(end)])
    grid on
    grid Minor
    xlabel('log10 \epsilon_{DNS}')
    ylabel('log10 \epsilon_{Vz}')
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig(['slope_DNS_Vspec_',num2str(pt),'.fig'])
        close all
    end
end

% plot and compare the fit slopes in the sampled spectra
figure(3000)
clf
subplot(1,3,1)
histogram(unConstSlpUTraj,edgesSlp,'Normalization','probability','DisplayStyle','stairs','EdgeColor','r','LineWidth',2)
hold on
histogram(unConstSlpVTraj,edgesSlp,'Normalization','probability','DisplayStyle','stairs','EdgeColor','g','LineWidth',2)
histogram(unConstSlpWTraj,edgesSlp,'Normalization','probability','DisplayStyle','stairs','EdgeColor','b','LineWidth',2)
plot([mean(mean(mean(unConstSlpUTraj))) mean(mean(mean(unConstSlpUTraj)))],[0 0.15],'r','LineWidth',2)
plot([mean(mean(mean(unConstSlpVTraj))) mean(mean(mean(unConstSlpVTraj)))],[0 0.15],'g','LineWidth',2)
plot([mean(mean(mean(unConstSlpWTraj))) mean(mean(mean(unConstSlpWTraj)))],[0 0.15],'b','LineWidth',2)
plot([-5/3 -5/3],[0 0.15],'--k','LineWidth',2)
xlim([edgesSlp(1) edgesSlp(end)])
ylim([0 0.02])
xlabel('fit slopes for all frequency bins')
ylabel('probability')
legend('slope u','slope v','slope w','mean slope u','mean slope v','mean slope w','ref. -5/3 slope','Location','NorthWest','FontSize',ftSz)
grid on
set(gca,'FontSize',ftSz)
subplot(1,3,2)
histogram(unConstSlpUFlSpTraj,edgesSlp,'Normalization','probability','DisplayStyle','stairs','EdgeColor','r','LineWidth',2)
hold on
histogram(unConstSlpVFlSpTraj,edgesSlp,'Normalization','probability','DisplayStyle','stairs','EdgeColor','g','LineWidth',2)
histogram(unConstSlpWFlSpTraj,edgesSlp,'Normalization','probability','DisplayStyle','stairs','EdgeColor','b','LineWidth',2)
plot([mean(mean(unConstSlpUFlSpTraj)) mean(mean(unConstSlpUFlSpTraj))],[0 0.15],'r','LineWidth',2)
plot([mean(mean(unConstSlpVFlSpTraj)) mean(mean(unConstSlpVFlSpTraj))],[0 0.15],'g','LineWidth',2)
plot([mean(mean(unConstSlpWFlSpTraj)) mean(mean(unConstSlpWFlSpTraj))],[0 0.15],'b','LineWidth',2)
plot([-5/3 -5/3],[0 0.15],'--k','LineWidth',2)
xlim([edgesSlp(1) edgesSlp(end)])
ylim([0 0.15])
xlabel('fit slopes')
ylabel('probability')
legend('slope u','slope v','slope w','mean slope u','mean slope v','mean slope w','ref. -5/3 slope','Location','NorthWest','FontSize',ftSz)
grid on
set(gca,'FontSize',ftSz)
subplot(1,3,3)
histogram(unConstSlpUusdTraj,edgesSlp,'Normalization','probability','DisplayStyle','stairs','EdgeColor','r','LineWidth',2)
hold on
histogram(unConstSlpVusdTraj,edgesSlp,'Normalization','probability','DisplayStyle','stairs','EdgeColor','g','LineWidth',2)
histogram(unConstSlpWusdTraj,edgesSlp,'Normalization','probability','DisplayStyle','stairs','EdgeColor','b','LineWidth',2)
plot([mean(mean(mean(unConstSlpUusdTraj,'omitnan'),'omitnan'),'omitnan') ...
    mean(mean(mean(unConstSlpUusdTraj,'omitnan'),'omitnan'),'omitnan')],...
    [0 0.15],'r','LineWidth',2)
plot([mean(mean(mean(unConstSlpVusdTraj,'omitnan'),'omitnan'),'omitnan') ...
    mean(mean(mean(unConstSlpVusdTraj,'omitnan'),'omitnan'),'omitnan')],...
    [0 0.15],'g','LineWidth',2)
plot([mean(mean(mean(unConstSlpWusdTraj,'omitnan'),'omitnan'),'omitnan') ...
    mean(mean(mean(unConstSlpWusdTraj,'omitnan'),'omitnan'),'omitnan')],...
    [0 0.15],'b','LineWidth',2)
plot([-5/3 -5/3],[0 0.15],'--k','LineWidth',2)
xlim([edgesSlp(1) edgesSlp(end)])
ylim([0 0.02])
xlabel('fit slopes')
ylabel('probability')
legend('slope u: usd. pts.','slope v: usd. pts.','slope w: usd. pts.',...
    'mean slope u: usd. pts.','mean slope v: usd. pts.','mean slope w: usd. pts.',...
    'ref. -5/3 slope','Location','NorthWest','FontSize',ftSz)
grid on
set(gca,'FontSize',ftSz)
if saSp == 1
	savefig('prob_fit_slopes.fig')
    close all
end

% Plot the contours of fit slopes and levels for all spectra and overlay the frequency bins used for qualified fits
smplN0 = 1:1:Nx;
pltdiffF = (f_nth_dec(1:end-1) + f_nth_dec(2:end)) / 2;
[Fx,Hz] = meshgrid(pltdiffF,smplN0);

% average the fit slopes to show a clear trend (if any)
for j = 1:1:(numInts/Trec)
    for k=1:length(pltdiffF)
        unConstSlpUTrajplt(k,j,:) = medfilt1(unConstSlpUTraj(k,j,:),3,'omitnan','truncate');
        unConstSlpVTrajplt(k,j,:) = medfilt1(unConstSlpVTraj(k,j,:),3,'omitnan','truncate');
        unConstSlpWTrajplt(k,j,:) = medfilt1(unConstSlpWTraj(k,j,:),3,'omitnan','truncate');
    end
end

slpArr = [-10 3*(-5/3) 2*(-5/3) 1.5*(-5/3) 1.4*(-5/3) 1.3*(-5/3) 1.2*(-5/3) 1.1*(-5/3) (-5/3) 0.9*(-5/3) 0.8*(-5/3) ...
    0.7*(-5/3) 0.6*(-5/3) 0.5*(-5/3) 0.1*(-5/3)];
tmpV = inferno;
cmap = tmpV(1:50:256,:);
% figure to plot the contours of fit slopes for HW data
for j = 1:1:(numInts/Trec)
    figure(3100+j)
    clf
    tlplt = tiledlayout(1,2,'TileSpacing','tight','Padding','compact');
    tl1 = nexttile;
    contourf(log10(Fx),Hz,reshape(unConstSlpUTrajplt(:,j,:),[length(pltdiffF) length(smplN0)])')
    clim([4*(-5/3) 0])
    col1 = colorbar('Ticks',[4*(-5/3) 3*(-5/3) 2*(-5/3) (-5/3) 0.5*(-5/3) 0],...
         'TickLabels',{'4*(-5/3)' '3*(-5/3)' '2*(-5/3)' '(-5/3)' '0.5*(-5/3)' '0'},'Location','southoutside');
    col1.Label.String = 'fit slope';
    colormap(cmap)
    xlim([min(log10(pltdiffF)) max(log10(pltdiffF))])
    ylim([min(smplN0) max(smplN0)])
    xlabel('log10(F)')
    ylabel('sample number')
    title('Contours of spectral slopes')
    tl2 = nexttile;
    contourf(log10(Fx),Hz,reshape(unConstSlpUTrajplt(:,j,:),[length(pltdiffF) length(smplN0)])',[-5/3-0.3 -5/3+0.3])
    colormap(cmap)
    xlim([min(log10(pltdiffF)) max(log10(pltdiffF))])
    ylim([min(smplN0) max(smplN0)])
    xlabel('log10(F) [Hz]')
    ylabel('sample number')
    title('Contours of slopes -5/3\pm0.3')
    % now superimpose the points used in fits
    hold on
    for i = 1:1:numTraj
        plot(log10(f_nth_dec(cell2mat(Uusd_inds(j,i)))),i.*ones([1 length(cell2mat(Uusd_inds(j,i)))]),'.r')
    end
    set(gcf, 'Position',  [1000, 50, 1000, 700])
    if saSp == 1
        savefig(['slp_cont_Uz_',num2str(j),'.fig'])
        close all
    end
    
    figure(3200+j)
    clf
    tlplt = tiledlayout(1,2,'TileSpacing','tight','Padding','compact');
    tl1 = nexttile;
    contourf(log10(Fx),Hz,reshape(unConstSlpVTrajplt(:,j,:),[length(pltdiffF) length(smplN0)])')
    clim([4*(-5/3) 0])
    col1 = colorbar('Ticks',[4*(-5/3) 3*(-5/3) 2*(-5/3) (-5/3) 0.5*(-5/3) 0],...
         'TickLabels',{'4*(-5/3)' '3*(-5/3)' '2*(-5/3)' '(-5/3)' '0.5*(-5/3)' '0'},'Location','southoutside');
    col1.Label.String = 'fit slope';
    colormap(cmap)
    xlim([min(log10(pltdiffF)) max(log10(pltdiffF))])
    ylim([min(smplN0) max(smplN0)])
    xlabel('log10(F)')
    ylabel('sample number')
    title('Contours of spectral slopes')
    tl2 = nexttile;
    contourf(log10(Fx),Hz,reshape(unConstSlpVTrajplt(:,j,:),[length(pltdiffF) length(smplN0)])',[-5/3-0.3 -5/3+0.3])
    colormap(cmap)
    xlim([min(log10(pltdiffF)) max(log10(pltdiffF))])
    ylim([min(smplN0) max(smplN0)])
    xlabel('log10(F) [Hz]')
    ylabel('sample number')
    title('Contours of slopes -5/3\pm0.3')
    % now superimpose the points used in fits
    hold on
    for i = 1:1:numTraj
        plot(log10(f_nth_dec(cell2mat(Vusd_inds(j,i)))),i.*ones([1 length(cell2mat(Vusd_inds(j,i)))]),'.r')
    end
    set(gcf, 'Position',  [1000, 50, 1000, 700])
    if saSp == 1
        savefig(['slp_cont_Vz_',num2str(j),'.fig'])
        close all
    end

    figure(3300+j)
    clf
    tlplt = tiledlayout(1,2,'TileSpacing','tight','Padding','compact');
    tl1 = nexttile;
    contourf(log10(Fx),Hz,reshape(unConstSlpWTrajplt(:,j,:),[length(pltdiffF) length(smplN0)])')
    clim([4*(-5/3) 0])
    col1 = colorbar('Ticks',[4*(-5/3) 3*(-5/3) 2*(-5/3) (-5/3) 0.5*(-5/3) 0],...
         'TickLabels',{'4*(-5/3)' '3*(-5/3)' '2*(-5/3)' '(-5/3)' '0.5*(-5/3)' '0'},'Location','southoutside');
    col1.Label.String = 'fit slope';
    colormap(cmap)
    xlim([min(log10(pltdiffF)) max(log10(pltdiffF))])
    ylim([min(smplN0) max(smplN0)])
    xlabel('log10(F)')
    ylabel('sample number')
    title('Contours of spectral slopes')
    tl2 = nexttile;
    contourf(log10(Fx),Hz,reshape(unConstSlpWTrajplt(:,j,:),[length(pltdiffF) length(smplN0)])',[-5/3-0.3 -5/3+0.3])
    colormap(cmap)
    xlim([min(log10(pltdiffF)) max(log10(pltdiffF))])
    ylim([min(smplN0) max(smplN0)])
    xlabel('log10(F) [Hz]')
    ylabel('sample number')
    title('Contours of slopes -5/3\pm0.3')
    % now superimpose the points used in fits
    hold on
    for i = 1:1:numTraj
        plot(log10(f_nth_dec(cell2mat(Wusd_inds(j,i)))),i.*ones([1 length(cell2mat(Wusd_inds(j,i)))]),'.r')
    end
    set(gcf, 'Position',  [1000, 50, 1000, 700])
    if saSp == 1
        savefig(['slp_cont_Wz_',num2str(j),'.fig'])
        close all
    end
end

%% plot the data on sampled plane
if diag == 1 && lean == 0
    figure(1)
    clf
    sUz = surface(Uz,'FaceColor','flat');
    sUz.EdgeColor = 'none';
    colorbar
    colormap('jet')
    xlim([1 Ny])
    ylim([1 Nz])
    xlabel('Y')
    ylabel('Z')
    grid on
    title('U sampled on YZ plane (at X/2)')
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig('YZU.fig')
        close all
    end
    figure(2)
    clf
    sVz = surface(Vz,'FaceColor','flat');
    sVz.EdgeColor = 'none';
    colorbar
    colormap('jet')
    xlim([1 Ny])
    ylim([1 Nz])
    xlabel('Y')
    ylabel('Z')
    grid on
    title('V sampled on YZ plane (at X/2)')
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig('YZV.fig')
        close all
    end
    figure(3)
    clf
    sWz = surface(Wz,'FaceColor','flat');
    sWz.EdgeColor = 'none';
    colorbar
    colormap('jet')
    xlim([1 Ny])
    ylim([1 Nz])
    xlabel('Y')
    ylabel('Z')
    grid on
    title('W sampled on YZ plane (at X/2)')
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig('YZW.fig')
        close all
    end
    figure(4)
    clf
    sTz = surface(Tz,'FaceColor','flat');
    sTz.EdgeColor = 'none';
    colorbar
    colormap('jet')
    xlim([1 Ny])
    ylim([1 Nz])
    xlabel('Y')
    ylabel('Z')
    grid on
    title('T sampled on YZ plane (at X/2)')
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig('YZT.fig')
        close all
    end
    figure(5)
    clf
    seps = surface(epsDNSz,'FaceColor','flat');
    seps.EdgeColor = 'none';
    colorbar
    colormap('jet')
    set(gca,'ColorScale','log')
    caxis([-5 -1.5])
    xlim([1 Ny])
    ylim([1 Nz])
    xlabel('Y')
    ylabel('Z')
    grid on
    title('log{\epsilon} sampled on YZ plane (at X/2)')
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig('YZE.fig')
        close all
    end
    figure(6)
    clf
    sTDR = surface(TDRDNSz,'FaceColor','flat');
    sTDR.EdgeColor = 'none';
    colorbar
    colormap('jet')
    set(gca,'ColorScale','log')
    caxis([-5 -1.5])
    xlim([1 Ny])
    ylim([1 Nz])
    xlabel('Y')
    ylabel('Z')
    grid on
    title('log{TDR} sampled on YZ plane (at X/2)')
    set(gca,'FontSize',ftSz)
    if saSp == 1
        savefig('YZTDR.fig')
        close all
    end
    % DIAGNOSTIC FIGURE TO CHECK BOX ORIENTATION! DO NOT DELETE
    %{
    figure(4)
    clf
    sWz = surface(Wz(:,1:Nx/2),'FaceColor','flat');
    sWz.EdgeColor = 'none';
    colorbar
    colormap('jet')
    caxis([-50*mean(mean(Wz)) 50*mean(mean(Wz))])
    xlim([1 Ny])
    ylim([1 Nz])
    set(gca,'FontSize',ftSz)
    %}
end
