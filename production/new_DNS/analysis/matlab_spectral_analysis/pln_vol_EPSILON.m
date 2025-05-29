%% script ot compute spectra, subsequently epsilon using the plane data
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

%% General Inputss
diag = 1;                               % switch to turn on program diagnostics
lean = 1;                               % unwanted figures are not plotted (saves time) 
Trec = 1;                               % Time record interval [s]
alp = 1.50;                             % Kolmogorov constant: 1.5 for 3D spectrum; 0.49-0.55 for longitudinal spectrum; 0.65 for transverse spectrum
alpL = 0.50;                            % Kolmogorov constant: 1.5 for 3D spectrum; 0.49-0.55 for longitudinal spectrum; 0.65 for transverse spectrum
alpT = 0.65;                            % Kolmogorov constant: 1.5 for 3D spectrum; 0.49-0.55 for longitudinal spectrum; 0.65 for transverse spectrum
winOn = 1;                              % switch to turn on or off the windowing
saSp = 1;                               % scitch to turn on or off saving figures and files
edges = -5:0.05:0;

% Simulation parameters set in sam.inp [ALL THESE QUANTITIES ARE SCALE NORMALIZED]
Xl = 01;                                % Xl = Domain X dimension
Yl = 01;                                % Yl = Domain Y dimension
Zl = 01;                                % Zl = Domain Z dimension
Nx = 3456;                              % Nx = Number of points in Domain X dimension
Ny = 3456;                              % Ny = Number of points in Domain Y dimension
Nz = 3456;                              % Nz = Number of points in Domain Z dimension
nuDNS = 1.5e-5;                         % The kinematic viscosity set in the DNS

% calculate DNS resolution parameters [ALL THESE QUANTITIES ARE SCALE NORMALIZED]
dx = Xl/Nx;                             % Grid Resolution in X
dy = Yl/Ny;                             % Grid Resolution in Y
dz = Zl/Nz;                             % Grid Resolution in Z

% Simulation outputs gathered from .out files [ALL THESE QUANTITIES ARE SCALE NORMALIZED]
resMet = mean([1.4602 1.4567 1.5124]);	% resMet = DNS resolution metric ( = grid spacing/ kolmogorov length scale)
DNSTKE = 1.4383;                        % DNS mean Energy at the sampled time step
epsDNS = 1.8641;                        % DNS mean TKEDR at the sampled time step
etaDNS = (nuDNS^3/epsDNS)^(1/4);        
tauDNS = (nuDNS/epsDNS)^(1/2);
muDNS = etaDNS/tauDNS;

% Scaling parameters applied from measurement [HYFLITS measurement at ~30 km AGL]
epsMeas = 2.225e-3;                     % epsMeas = TKE dissipation rate [m^3s^-2]
nuMeas = 4.0e-4;                        % nu = Kinematic viscosity [m^2s^-1]
balRt = 02.5;                           % balRt = HYFLITS balloon descent rate [m/s]
etaMEas = (nuMeas^3/epsMeas)^(1/4);     % Kolmogorov length scale [m]
tauMeas = (nuMeas/epsMeas)^(1/2);       % Kolmogorov time scale [s]
muMeas = etaMEas/tauMeas;               % Kolmogorov velocity scale [m/s]

% Scaling ratios: The DNS data should be scaled using the factors below
ltScal = etaMEas/etaDNS;                % change in length scale
tmScl = tauMeas/tauDNS;                 % change in time scale
velScl = ltScal/tmScl;                  % change in velocity scale
epsScal = epsMeas/epsDNS;               % change in epsilon

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

%% Load the dataset from file(s)
% for u
dirtry = [pwd '/'];
flPlnU = strcat(dirtry,'yz1_0140_U.txt');
temp = table2array(readtable(flPlnU));
PlnU = transpose(reshape(temp,[Nx+1,Ny+1]));
clear temp
% for v
flPlnV = strcat(dirtry,'yz1_0140_V.txt');
temp = table2array(readtable(flPlnV));
PlnV = transpose(reshape(temp,[Nx+1,Ny+1]));
clear temp
% for W
flPlnW = strcat(dirtry,'yz1_0140_W.txt');
temp = table2array(readtable(flPlnW));
PlnW = transpose(reshape(temp,[Nx+1,Ny+1]));
clear temp
% for DNS epsilon
flPlnE = strcat(dirtry,'yz1_0140_E.txt');
temp = table2array(readtable(flPlnE));
PlnE = transpose(reshape(temp,[Nx+1,Ny+1]));
clear temp
% Reshape the W array for the plane
for i = 1:1:length(PlnW(1,1:end-1))
    Uz(:,i) = PlnU(1:end-1,i).*velScl;
    Vz(:,i) = PlnV(1:end-1,i).*velScl;
	Wz(:,i) = PlnW(1:end-1,i).*velScl;
	epsDNSz(:,i) = PlnE(1:end-1,i).*epsScal;
end

%% plot diagnostic figures
if diag == 1 && lean == 1
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
    if saSp == 1
        savefig('YZW.fig')
        close all
    end
    figure(4)
    clf
    seps = surface(epsDNSz,'FaceColor','flat');
    seps.EdgeColor = 'none';
    colorbar
    colormap('jet')
    caxis([1e-5 10*mean(mean(epsDNSz))])
    xlim([1 Ny])
    ylim([1 Nz])
    xlabel('Y')
    ylabel('Z')
    grid on
    title('U\epsilon sampled on YZ plane (at X/2)')
    if saSp == 1
        savefig('YZE.fig')
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
    %}
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
pit_NF = 10^(-4+log10(alpL)+(2/3)*log10(balRt));                % define the noise floor for pitot 
min_fit_pts = 2;                                                % the minimum number of points to be used in the third pass of fitting
pass_3 = 1;                                                     % This switch turns on an additional level scan on the spectra -- provides conservative turbulence estimates
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

% Plot Histograms comparing spectrally derived TKEDR and averaged DNS pointwise extimates of TKEDR
figure(401)
clf
histogram(log_Wepsilon_k,edges,'Normalization','probability','FaceColor', 'b')
hold on
histogram(log10(mnepsDNSz),edges,'Normalization','probability','FaceColor', 'r')
plot([mean(mean(log_Wepsilon_k)) mean(mean(log_Wepsilon_k))],[0 0.25],'b','LineWidth',2)
plot([mean(mean(log10(mnepsDNSz))) mean(mean(log10(mnepsDNSz)))],[0 0.25],'r','LineWidth',2)
xlim([edges(1) edges(end)])
ylim([0 0.25])
xlabel('log10 \epsilon')
ylabel('probability')
legend('\epsilon_{spec}','\epsilon_{DNS} - averaged over data records','Location','NorthEast')
grid on
if saSp == 1
	savefig('prob_spec_DNS.fig')
    close all
end

figure(402)
clf
histogram(log_Uepsilon_k,edges,'Normalization','probability','FaceColor', 'r')
hold on
histogram(log_Vepsilon_k,edges,'Normalization','probability','FaceColor', 'g')
histogram(log_Wepsilon_k,edges,'Normalization','probability','FaceColor', 'b')
plot([mean(mean(log_Uepsilon_k)) mean(mean(log_Uepsilon_k))],[0 0.15],'r','LineWidth',2)
plot([mean(mean(log_Vepsilon_k)) mean(mean(log_Vepsilon_k))],[0 0.15],'g','LineWidth',2)
plot([mean(mean(log_Wepsilon_k)) mean(mean(log_Wepsilon_k))],[0 0.15],'b','LineWidth',2)
xlim([edges(1) edges(end)])
ylim([0 0.15])
xlabel('log10 \epsilon')
ylabel('probability')
legend('\epsilon_{specU}','\epsilon_{specV}','\epsilon_{specW}','Location','NorthEast')
grid on
if saSp == 1
	savefig('prob_spec_UVW.fig')
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
            'error bars','Location','NorthEast')
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa')
        ylabel('averaged Wz PSD')
        title('normalized wavenumber vs averaged Wz Spectra sampled on YZ plane (at X/2)')
        grid on
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
            'error bars','Location','NorthEast')
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa')
        ylabel('averaged Uz PSD')
        title('normalized wavenumber vs averaged Uz Spectra sampled on YZ plane (at X/2)')
        grid on
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
            'error bars','Location','NorthEast')
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa')
        ylabel('averaged Vz PSD')
        title('normalized wavenumber vs averaged Vz Spectra sampled on YZ plane (at X/2)')
        grid on
        if saSp == 1
            savefig(['spec_Vz_avg_',num2str(pt),'.fig'])
            close all
        end
        
        figure(800+pt)
        clf
        plot(log10(mnepsDNSz(pt,:)),'r')
        hold on
        plot(log_Wepsilon_k(pt,:),'b')
        plot([1 length(mnepsDNSz(pt,:))],[log10(avTrepsDNSz(pt)) log10(avTrepsDNSz(pt))],'r','LineWidth',2)
        plot([1 length(log_Wepsilon_k(pt,:))],[log_avWepsilon_k(pt) log_avWepsilon_k(pt)],'b','LineWidth',2)
        legend('\epsilon_{DNS}','\epsilon_{specW}','\epsilon_{DNS} averaged','\epsilon_{DNS} averaged','Location','NorthEast')
        xlabel('data index')
        ylabel('\epsilon')
        grid on
        if saSp == 1
            savefig(['eps_DNS_Wspec_level_',num2str(pt),'.fig'])
            close all
        end
        
        figure(900+pt)
        clf
        plot(log10(mnepsDNSz(pt,:)),'r')
        hold on
        plot(log_Uepsilon_k(pt,:),'b')
        plot([1 length(mnepsDNSz(pt,:))],[log10(avTrepsDNSz(pt)) log10(avTrepsDNSz(pt))],'r','LineWidth',2)
        plot([1 length(log_Uepsilon_k(pt,:))],[log_avUepsilon_k(pt) log_avUepsilon_k(pt)],'b','LineWidth',2)
        legend('\epsilon_{DNS}','\epsilon_{specU}','\epsilon_{DNS} averaged','\epsilon_{DNS} averaged','Location','NorthEast')
        xlabel('data index')
        ylabel('\epsilon')
        grid on
        if saSp == 1
            savefig(['eps_DNS_Uspec_level_',num2str(pt),'.fig'])
            close all
        end
        
        figure(1000+pt)
        clf
        plot(log10(mnepsDNSz(pt,:)),'r')
        hold on
        plot(log_Vepsilon_k(pt,:),'b')
        plot([1 length(mnepsDNSz(pt,:))],[log10(avTrepsDNSz(pt)) log10(avTrepsDNSz(pt))],'r','LineWidth',2)
        plot([1 length(log_Vepsilon_k(pt,:))],[log_avVepsilon_k(pt) log_avVepsilon_k(pt)],'b','LineWidth',2)
        legend('\epsilon_{DNS}','\epsilon_{specV}','\epsilon_{DNS} averaged','\epsilon_{DNS} averaged','Location','NorthEast')
        xlabel('data index')
        ylabel('\epsilon')
        grid on
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
        legend('Wz','Wz - detrended','Wz - mean removed')
        ylim([-3 3])
        subplot(1,2,2)
        semilogx(freq,log10(abs(WppdfPlt(:,j,i))),'b','LineWidth',2)
        hold on
        semilogx(freq,WTrlog_pwp_k_avg(j,i)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,WTrlog_ppsd_avg(:,j,i),'r*','LineWidth',4)
        semilogx(f_nth_dec(Wusd_inds{j,i}),WTrlog_ppsd_avg(Wusd_inds{j,i},j,i),'go','LineWidth',2)
        semilogx(freq,(WTrlog_pwp_k_avg(j,i)+WTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(WTrlog_pwp_k_avg(j,i)-WTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)        
        legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','NorthEast')
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa')
        ylabel('Wz PSD')
        title(['normalized wavenumber vs Wz Spectra, Trajectory = ',num2str(i)])
        grid on
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
        legend('Uz','Uz - detrended','Uz - mean removed')
        ylim([-3 3])
        subplot(1,2,2)
        semilogx(freq,log10(abs(UppdfPlt(:,j,i))),'b','LineWidth',2)
        hold on
        semilogx(freq,UTrlog_pwp_k_avg(j,i)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,UTrlog_ppsd_avg(:,j,i),'r*','LineWidth',4)
        semilogx(f_nth_dec(Uusd_inds{j,i}),UTrlog_ppsd_avg(Uusd_inds{j,i},j,i),'go','LineWidth',2)
        semilogx(freq,(UTrlog_pwp_k_avg(j,i)+UTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(UTrlog_pwp_k_avg(j,i)-UTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)        
        legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','NorthEast')
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa')
        ylabel('Uz PSD')
        title(['normalized wavenumber vs Uz Spectra, Trajectory = ',num2str(i)])
        grid on
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
        legend('Vz','Vz - detrended','Vz - mean removed')
        ylim([-3 3])
        subplot(1,2,2)
        semilogx(freq,log10(abs(VppdfPlt(:,j,i))),'b','LineWidth',2)
        hold on
        semilogx(freq,VTrlog_pwp_k_avg(j,i)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,VTrlog_ppsd_avg(:,j,i),'r*','LineWidth',4)
        semilogx(f_nth_dec(Vusd_inds{j,i}),VTrlog_ppsd_avg(Vusd_inds{j,i},j,i),'go','LineWidth',2)
        semilogx(freq,(VTrlog_pwp_k_avg(j,i)+VTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(VTrlog_pwp_k_avg(j,i)-VTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)        
        legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','NorthEast')
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa')
        ylabel('Vz PSD')
        title(['normalized wavenumber vs Vz Spectra, Trajectory = ',num2str(i)])
        grid on
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
    %for pt = 1:1:(numInts/Trec)
    %    semilogx(freqR.*(etaMEas),avWppsdw(:,pt).*(mean(10.^(log_avWepsilon_k)))^(-2/3),'g','LineWidth',2)
    %end
    loglog([freqE(1).*(etaDNS) freqE(end).*(etaDNS)],[alp alp],'--k')
    loglog([freqE(1).*(etaDNS) freqE(end).*(etaDNS)],[alpT alpT],'--r')
    loglog([freqE(1).*(etaDNS) freqE(end).*(etaDNS)],[alpT alpT],'--m')
    loglog([freqE(1).*(etaDNS) freqE(end).*(etaDNS)],[alpL alpL],'--b')
    lgd = legend('{\epsilon}^{-2/3}\kappa^{5/3}E(\kappa)',...
        '{\epsilon}^{-2/3}{\kappa_{3}}^{5/3}E_{11}(\kappa_{3})',...
        '{\epsilon}^{-2/3}{\kappa_{3}}^{5/3}E_{22}(\kappa_{3})',...
        '{\epsilon}^{-2/3}{\kappa_{3}}^{5/3}E_{33}(\kappa_{3})',...
        'Location','SouthWest');
    lgd.FontSize = 14;
    %axis([freqE(1).*(etaDNS) freqE(end).*(etaDNS) 0 2])
    xlabel('\kappa\eta','FontSize',24)
    ylabel('{\epsilon}^{-2/3}\kappa^{5/3}E(\kappa)','FontSize',24)
    title('Compensated spectra to verify universal constants','FontSize',24)
    grid on
    grid Minor
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
PlnE_mean = log10(mean(mean(epsDNSz)));                                     % avg DNS on the sampled XZ plane

% Plot comparing individual sample Epsilon calculated using Wz profiles and
% the average epsilon for the plane
for pt = 1:1:(numInts/Trec)
    figure(2100+pt)
    clf
    plot(log_Wepsilon_k(pt,:),'g','LineWidth',2)        % spec - records
    hold on
    plot(log10(mnepsDNSz(pt,:)),'r','LineWidth',2)             % DNS - records
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
    set(h,'Interpreter','latex','fontsize',16)
    if saSp == 1
        savefig(['epsilon_Wz_pln_compare_',num2str(pt),'.fig'])
        close all
    end
    figure(2200+pt)
    clf
    plot(log_Uepsilon_k(pt,:),'g','LineWidth',2)        % spec - records
    hold on
    plot(log10(mnepsDNSz(pt,:)),'r','LineWidth',2)             % DNS - records
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
    set(h,'Interpreter','latex','fontsize',16)
    if saSp == 1
        savefig(['epsilon_Uz_pln_compare_',num2str(pt),'.fig'])
        close all
    end
    figure(2300+pt)
    clf
    plot(log_Vepsilon_k(pt,:),'g','LineWidth',2)        % spec - records
    hold on
    plot(log10(mnepsDNSz(pt,:)),'r','LineWidth',2)             % DNS - records
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
    set(h,'Interpreter','latex','fontsize',16)
    if saSp == 1
        savefig(['epsilon_Vz_pln_compare_',num2str(pt),'.fig'])
        close all
    end
end
% Plot comparing histogram of spectral and DNS derived TKEDR
for pt = 1:1:(numInts/Trec)
    figure(2400+pt)
    clf
    histogram(log_Wepsilon_k(pt,:),edges,'Normalization','probability','FaceColor', 'b')
    hold on
    histogram(log10(mnepsDNSz(pt,:)),edges,'Normalization','probability','FaceColor', 'r')
    plot([log_avWepsilon_k(pt) log_avWepsilon_k(pt)],[0 0.5],'g','LineWidth',2)
    plot([PlnE_mean PlnE_mean],[0 0.5],'r','LineWidth',2)
    plot([log_volepsilon_kW log_volepsilon_kW],[0 0.5],'m','LineWidth',2)
    plot([logDNSEPSMN logDNSEPSMN],[0 0.5],'k','LineWidth',2)
    xlim([edges(1) edges(end)])
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
    set(h,'Interpreter','latex','fontsize',16)
    if saSp == 1
        savefig(['hist_Wz_',num2str(pt),'.fig'])
        close all
    end
    figure(2500+pt)
    clf
    histogram(log_Uepsilon_k(pt,:),edges,'Normalization','probability','FaceColor', 'b')
    hold on
    histogram(log10(mnepsDNSz(pt,:)),edges,'Normalization','probability','FaceColor', 'r')
    plot([log_avUepsilon_k(pt) log_avUepsilon_k(pt)],[0 0.5],'g','LineWidth',2)
    plot([PlnE_mean PlnE_mean],[0 0.5],'r','LineWidth',2)
    plot([log_volepsilon_kU log_volepsilon_kU],[0 0.5],'m','LineWidth',2)
    plot([logDNSEPSMN logDNSEPSMN],[0 0.5],'k','LineWidth',2)
    xlim([edges(1) edges(end)])
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
    set(h,'Interpreter','latex','fontsize',16)
    if saSp == 1
        savefig(['hist_Uz_',num2str(pt),'.fig'])
        close all
    end
    figure(2600+pt)
    clf
    histogram(log_Vepsilon_k(pt,:),edges,'Normalization','probability','FaceColor', 'b')
    hold on
    histogram(log10(mnepsDNSz(pt,:)),edges,'Normalization','probability','FaceColor', 'r')
    plot([log_avVepsilon_k(pt) log_avVepsilon_k(pt)],[0 0.5],'g','LineWidth',2)
    plot([PlnE_mean PlnE_mean],[0 0.5],'r','LineWidth',2)
    plot([log_volepsilon_kV log_volepsilon_kV],[0 0.5],'m','LineWidth',2)
    plot([logDNSEPSMN logDNSEPSMN],[0 0.5],'k','LineWidth',2)
    xlim([edges(1) edges(end)])
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
    set(h,'Interpreter','latex','fontsize',16)
    if saSp == 1
        savefig(['hist_Vz_',num2str(pt),'.fig'])
        close all
    end
    
    figure(2700+pt)
    clf
    plot(log10(mnepsDNSz(pt,:)),log_Wepsilon_k(pt,:),'*k')
    hold on
    plot([edges(1) edges(end)],[edges(1) edges(end)],'k','LineWidth',2)
    xlim([edges(1) edges(end)])
    ylim([edges(1) edges(end)])
    grid on
    grid Minor
    xlabel('log10 \epsilon_{DNS}')
    ylabel('log10 \epsilon_{Wz}')
    if saSp == 1
        savefig(['slope_DNS_Wspec_',num2str(pt),'.fig'])
        close all
    end
    
    figure(2800+pt)
    clf
    plot(log10(mnepsDNSz(pt,:)),log_Uepsilon_k(pt,:),'*k')
    hold on
    plot([edges(1) edges(end)],[edges(1) edges(end)],'k','LineWidth',2)
    xlim([edges(1) edges(end)])
    ylim([edges(1) edges(end)])
    grid on
    grid Minor
    xlabel('log10 \epsilon_{DNS}')
    ylabel('log10 \epsilon_{Uz}')
    if saSp == 1
        savefig(['slope_DNS_Uspec_',num2str(pt),'.fig'])
        close all
    end
    
    figure(2900+pt)
    clf
    plot(log10(mnepsDNSz(pt,:)),log_Vepsilon_k(pt,:),'*k')
    hold on
    plot([edges(1) edges(end)],[edges(1) edges(end)],'k','LineWidth',2)
    xlim([edges(1) edges(end)])
    ylim([edges(1) edges(end)])
    grid on
    grid Minor
    xlabel('log10 \epsilon_{DNS}')
    ylabel('log10 \epsilon_{Vz}')
    if saSp == 1
        savefig(['slope_DNS_Vspec_',num2str(pt),'.fig'])
        close all
    end
end