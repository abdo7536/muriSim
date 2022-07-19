%% Script to compare TKE dissipation Rate estimates from Synthetic Observations of DNS data
%% Author: Abhiram Doddi
%% Date Composed: 11th April 2022

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
exectry = pwd;
saSp = 1;
% Plot Inputs: All inputs are mandatory
ftsz = 22;

%% Setup calculations
% Calculate scale parameters for DNS data
eta = (nu^3/epsMeas)^(1/4);         % Kolmogorov length scale [m]
tau = (nu/epsMeas)^(1/2);           % Kolmogorov time scale [s]
mu = eta/tau;                       % Kolmogorov velocity scale [m/s]
epsStar = nu*(mu/eta)^2;            % Scaling to be applied for TKE dissipation rate [m^3/s^2]
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
epsScalDNS = ((dz/resMet)^(-4))*(nuDNS^(3));   % The DNS co-incidental scale with measurement EPSILON
% Compute the number of sample points per interval
numInts = floor(Zscal/balRt);
numSmpls = floor(balRt/dzscal);

%% Load Data
% For DNS Epsilon
flGxEx = strcat(dirtry,'GridXEx000100_004000.txt');
flGxEy = strcat(dirtry,'GridXEy000100_004000.txt');
flGxEz = strcat(dirtry,'GridXEz000100_004000.txt');
flGyEx = strcat(dirtry,'GridYEx000100_004000.txt');
flGyEy = strcat(dirtry,'GridYEy000100_004000.txt');
flGyEz = strcat(dirtry,'GridYEz000100_004000.txt');
flGzEx = strcat(dirtry,'GridZEx000100_004000.txt');
flGzEy = strcat(dirtry,'GridZEy000100_004000.txt');
flGzEz = strcat(dirtry,'GridZEz000100_004000.txt');
flEx = strcat(dirtry,'Ex000100_004000.txt');
flEy = strcat(dirtry,'Ey000100_004000.txt');
flEz = strcat(dirtry,'Ez000100_004000.txt');
GxEx = table2array(readtable(flGxEx));
GxEy = table2array(readtable(flGxEy));
GxEz = table2array(readtable(flGxEz));
GyEx = table2array(readtable(flGyEx));
GyEy = table2array(readtable(flGyEy));
GyEz = table2array(readtable(flGyEz));
GzEx = table2array(readtable(flGzEx));
GzEy = table2array(readtable(flGzEy));
GzEz = table2array(readtable(flGzEz));
Ex = table2array(readtable(flEx)).*epsStar;
Ey = table2array(readtable(flEy)).*epsStar;
Ez = table2array(readtable(flEz)).*epsStar;
% For U'
flGxUx = strcat(dirtry,'GridXUx000100_004000.txt');
flGyUx = strcat(dirtry,'GridYUx000100_004000.txt');
flGzUx = strcat(dirtry,'GridZUx000100_004000.txt');
flUx = strcat(dirtry,'Ux000100_004000.txt');
GxUx = table2array(readtable(flGxUx));
GyUx = table2array(readtable(flGyUx));
GzUx = table2array(readtable(flGzUx));
Ux = table2array(readtable(flUx)).*mu;
% For V'
flGxVy = strcat(dirtry,'GridXVy000100_004000.txt');
flGyVy = strcat(dirtry,'GridYVy000100_004000.txt');
flGzVy = strcat(dirtry,'GridZVy000100_004000.txt');
flVy = strcat(dirtry,'Vy000100_004000.txt');
GxVy = table2array(readtable(flGxVy));
GyVy = table2array(readtable(flGyVy));
GzVy = table2array(readtable(flGzVy));
Vy = table2array(readtable(flVy)).*mu;
% For W'
flGxWz = strcat(dirtry,'GridXWz000100_004000.txt');
flGyWz = strcat(dirtry,'GridYWz000100_004000.txt');
flGzWz = strcat(dirtry,'GridZWz000100_004000.txt');
flWz = strcat(dirtry,'Wz000100_004000.txt');
GxWz = table2array(readtable(flGxWz));
GyWz = table2array(readtable(flGyWz));
GzWz = table2array(readtable(flGzWz));
Wz = table2array(readtable(flWz)).*mu;

% Compute the number of trajectories
tmpSz = size(Wz);
numTraj = tmpSz(2);
sampPts = tmpSz(1);

%% Compute frequency averaging and inputs for Spectral analysis
Nrec = numSmpls;
f_low_avg = 2;          % low frequency limit for spectral averaging [Hz]
f_high_avg = Nrec/2;    % high frequency limit for spectral averaging [Hz]
pts_in_dec = 1/5;       % [1/3 = 1/3rd decade averaging of the spectra] [set 1/10 for 1/10th decade averaging of the spectra]
pit_NF = 10^-18;        % define the noise floor for pitot 
min_fit_pts = 2;        % the minimum number of points to be used in the third pass of fitting
pass_3 = 1;             % This switch turns on an additional level scan on the spectra -- provides conservative turbulence estimates
time_segment_center_inds = (Nrec/2):Nrec:length(GzEz);
dt = 1/Nrec;  % sampling period in [sec]
f_low = f_low_avg;              % low frequency limit for spectral averaging [Hz]
f_high = f_high_avg;            % high frequency limit for spectral averaging [Hz]
freq = 1:Nrec/2/(Nrec*dt);    % spectral frequency samples, up to Nyquist rate [Hz]
f_ind_low = find(freq >= f_low,1,'first');
f_ind_high = find(freq <= f_high,1,'last');
freqs_used = freq(f_ind_low:f_ind_high);
% Setup the frequency boundaries for fractional decade averaging 
decades = log10(f_high_avg)-log10(f_low_avg);
f_avg = logspace(log10(f_low_avg),log10(f_high_avg),floor(decades/pts_in_dec));
f_nth_dec = (f_avg(2:end)+f_avg(1:end-1))/2;           % bin center frequencies
f_inds = zeros(1,length(f_avg));
for i = 1:length(f_avg)
    f_inds(i) = find(freq-f_avg(i) >= -1e-6,1,'First') - 1; % index of first frequency in each bin
end

f_avgind_low = freq(f_inds(1:end-1));
f_avgind_high = freq(f_inds(2:end)-1);
freqs_used_avg = f_nth_dec;

%% Loop over U trajectories (includes averaging of TKE DNS)
for i = 1:1:numTraj
    Utrj = Ux(:,i);
    Extrj = Ex(:,i);
    GxUxtrj = GxUx(:,i);
    GyUxtrj = GyUx(:,i);
    GzUxtrj = GzUx(:,i);
    GxExtrj = GxEx(:,i);
    GyExtrj = GyEx(:,i);
    GzExtrj = GzEx(:,i);
    for j = 1:1:numInts
        start_index = time_segment_center_inds(j) - Nrec/2;
        stop_index =  time_segment_center_inds(j) + Nrec/2;
        
        Uprec(:,j) = abs(Utrj(start_index+1:stop_index));  % time records of relative wind in [m/s]
        % detrend and weight the data to reduce windowing artifacts, using a variance-preserving scale factor
        Uprecd(:,j) = detrend(Uprec(:,j)).*sqrt(2.66).*hanning(Nrec,'periodic');
    
        % fft the records to obtain the power spectral density
        Ups(:,j) = 2*fft(Uprecd(:,j))/Nrec;    % amplitude spectrum [m/s]
        Upps(:,j) = real(conj(Ups(:,j)).*Ups(:,j));   % power spectrum [m^2/s^2]
        Uppsd(:,j) = Upps(1:Nrec/2,j)*dt*Nrec/2;      % power spectral density up to Nyquist [m^2/s^2/Hz] = [m^2/s]
    
        % weight PSD by f^(5/3)   
        Uppsdw(:,j) =  Uppsd(:,j).*freq'.^(5/3);    % [m^2 s^(-8/3)] 
    
        % average the weighted PSD:
        for q = 1:length(f_inds)-1
            if q == 1
                if (f_inds(q+1) - f_inds(q)) == 0
                    Uppsdw_avg(q,j) = Uppsdw(f_inds(q+1),j);
                else
                    Uppsdw_avg(q,j) = mean( Uppsdw(f_inds(q):f_inds(q+1)-1,j));
                end
            else
                if (f_inds(q+1) - f_inds(q)) == 0
                    Uppsdw_avg(q,j) = Uppsdw(f_inds(q+1),j);
                else
                    Uppsdw_avg(q,j) = mean( Uppsdw(f_inds(q)+1:f_inds(q+1)-1,j));
                end
            end
        end
        % store the averaged spectral points in this array (for comparison)
        UPPSDW_AVG(:,j) = Uppsdw_avg(:,j);
   
        % average weighted power spectral density from f_low to f_high, [m^2 s^(-8/3)]
        Upwp(j) = sum(Uppsdw(f_ind_low:f_ind_high,j))/(f_ind_high-f_ind_low+1);
        Upwp_avg(j) = sum(Uppsdw_avg(:,j))/length(f_nth_dec);

        % remove data in the interval (f_ind_low:f_ind_high) at or below the sensor noise floor
        Upit_NFw = pit_NF*freq(f_ind_low:f_ind_high).^(5/3);
        Upit_NFw_avg = pit_NF*f_nth_dec.^(5/3);
        Upwp_est(:,j) = Uppsdw(f_ind_low:f_ind_high,j);
        Upwp_est_avg(:,j) = Uppsdw_avg(:,j);
        Uk_inds = find(Upwp_est(:,j) > Upit_NFw');   % indexes of data points to keep
        Uk_inds_avg = find(Upwp_est_avg(:,j) > Upit_NFw_avg');   % indexes of data points to keep
        Ulog_pwp_k(j) = mean(log10(abs(Upwp_est(Uk_inds,j))));      % average over log of kept indexes
        Ulog_pwp_ks(j) = std(log10(abs(Upwp_est(Uk_inds,j))));      % standard deviation of the log mean
        Ulog_pwp_k_avg(j) = mean(log10(abs(Upwp_est_avg(Uk_inds_avg,j)))); % average over log of kept indexes
        Ulog_pwp_ks_avg(j) = std(log10(abs(Upwp_est_avg(Uk_inds_avg,j)))); % standard deviation of the log mean
   
        % unweight these estimates for plotting and checking
        % for the average spectra
        Ulog_ppsd_est(:,j) = Ulog_pwp_k_avg(j) + log10(f_nth_dec.^(-5/3));
        Ulog_ppsd_avg(:,j) = log10(Uppsdw_avg(:,j))' + log10(f_nth_dec.^(-5/3));
        
        % second pass: include indexes where the fit function is above the noise floor
        Uk_inds_avg = find(Ulog_ppsd_est(:,j) > log10(Upit_NFw_avg.*f_nth_dec.^(-5/3))');   % indexes of data points to keep
        
        if pass_3 == 1
            % third pass: include the points that contribute towards reducing the
            % standard deviation
            if length(Uk_inds_avg) > min_fit_pts
                clear stand inds_to_use
                % use the first three points to begin with and then iterate over the remaining points 
                Uind_to_use = Uk_inds_avg(1:min_fit_pts);
                for p = 1:1:length(Uk_inds_avg)-(min_fit_pts)+1
                   Ustand(p) = std(log10(abs(Upwp_est_avg(Uind_to_use,j))));
                   if p>1 && p<length(Uk_inds_avg)-(min_fit_pts)+1
                       if Ustand(p) > 0 && Ustand(p) < 0.25
                          if Ustand(p) > 1.2*Ustand(p-1)
                             Uind_to_use = [Uk_inds_avg(Uind_to_use(1:end-1)); Uk_inds_avg(min_fit_pts+p)];
                             Ustand(p) = Ustand(p-1);
                          else
                             Uind_to_use = [Uk_inds_avg(Uind_to_use); Uk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Ustand(p) > 0.25 && Ustand(p) < 0.5
                          if Ustand(p) > 1.15*Ustand(p-1)
                             Uind_to_use = [Uk_inds_avg(Uind_to_use(1:end-1)); Uk_inds_avg(min_fit_pts+p)];
                             Ustand(p) = Ustand(p-1);
                          else
                             Uind_to_use = [Uk_inds_avg(Uind_to_use); Uk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Ustand(p) > 0.5 && Ustand(p) < 0.75
                          if Ustand(p) > 1.12*Ustand(p-1)
                             Uind_to_use = [Uk_inds_avg(Uind_to_use(1:end-1)); Uk_inds_avg(min_fit_pts+p)];
                             Ustand(p) = Ustand(p-1);
                          else
                             Uind_to_use = [Uk_inds_avg(Uind_to_use); Uk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Ustand(p) > 0.75 && Ustand(p) < 1
                          if Ustand(p) > 1.1*Ustand(p-1)
                             Uind_to_use = [Uk_inds_avg(Uind_to_use(1:end-1)); Uk_inds_avg(min_fit_pts+p)];
                             Ustand(p) = Ustand(p-1);
                          else
                             Uind_to_use = [Uk_inds_avg(Uind_to_use); Uk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Ustand(p) > 1 && Ustand(p) < 1.5
                          if Ustand(p) > 1.05*Ustand(p-1)
                             Uind_to_use = [Uk_inds_avg(Uind_to_use(1:end-1)); Uk_inds_avg(min_fit_pts+p)];
                             Ustand(p) = Ustand(p-1);
                          else
                             Uind_to_use = [Uk_inds_avg(Uind_to_use); Uk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Ustand(p) >= 1.5
                          if Ustand(p) > 1.02*Ustand(p-1)
                             Uind_to_use = [Uk_inds_avg(Uind_to_use(1:end-1)); Uk_inds_avg(min_fit_pts+p)];
                             Ustand(p) = Ustand(p-1);
                          else
                             Uind_to_use = [Uk_inds_avg(Uind_to_use); Uk_inds_avg(min_fit_pts+p)];
                          end
                       end
                   elseif p == 1
                       Uind_to_use = Uk_inds_avg(1:min_fit_pts+p);
                   elseif p == length(Uk_inds_avg)-(min_fit_pts)+1
                       if Ustand(p) > 0 && Ustand(p) < 0.25
                          if Ustand(p) > 1.20*Ustand(p-1)
                             Uind_to_use = Uk_inds_avg(Uind_to_use(1:end-1));
                             Ustand(p) = Ustand(p-1);
                          else
                             Uind_to_use = Uk_inds_avg(Uind_to_use);
                          end
                       end
                       if Ustand(p) > 0.25 && Ustand(p) < 0.5
                          if Ustand(p) > 1.15*Ustand(p-1)
                             Uind_to_use = Uk_inds_avg(Uind_to_use(1:end-1));
                             Ustand(p) = Ustand(p-1);
                          else
                             Uind_to_use = Uk_inds_avg(Uind_to_use);
                          end
                       end
                       if Ustand(p) > 0.5 && Ustand(p) < 0.75
                          if Ustand(p) > 1.12*Ustand(p-1)
                             Uind_to_use = Uk_inds_avg(Uind_to_use(1:end-1));
                             Ustand(p) = Ustand(p-1);
                          else
                             Uind_to_use = Uk_inds_avg(Uind_to_use);
                          end
                       end
                       if Ustand(p) > 0.75 && Ustand(p) < 1
                          if Ustand(p) > 1.1*Ustand(p-1)
                             Uind_to_use = Uk_inds_avg(Uind_to_use(1:end-1));
                             Ustand(p) = Ustand(p-1);
                          else
                             Uind_to_use = Uk_inds_avg(Uind_to_use);
                          end
                       end
                       if Ustand(p) > 1 && Ustand(p) < 1.5
                          if Ustand(p) > 1.05*Ustand(p-1)
                             Uind_to_use = Uk_inds_avg(Uind_to_use(1:end-1));
                             Ustand(p) = Ustand(p-1);
                          else
                             Uind_to_use = Uk_inds_avg(Uind_to_use);
                          end
                       end
                       if Ustand(p) >= 1.5
                          if Ustand(p) > 1.02*Ustand(p-1)
                             Uind_to_use = Uk_inds_avg(Uind_to_use(1:end-1));
                             Ustand(p) = Ustand(p-1);
                          else
                             Uind_to_use = Uk_inds_avg(Uind_to_use);
                          end
                       end
                       Uk_inds_avg = Uind_to_use;
                   end
                end
            end
        end
        
        if exist('stand','var')
            Ustand_pit{j} = Ustand;
        else
            Ustand_pit{j} = nan;
        end
        Uk_inds_pit{j} = Uk_inds_avg;            % store the indices used for each spectral fit in a structure (the array lengths may vary)
    
        % if the number of selected points are less than the min number of
        % poits needed (set by the user) to get a decent fit estimate, then
        % don't bother with such spectra [set to nan by default - generally indicates that all the points in the spectra are suspect]
        Ulog_pwp_k_avg(j) = mean(log10(abs(Upwp_est_avg(Uk_inds_avg,j)))); % average over log of kept indexes - in second pass
        Ulog_pwp_ks_avg(j) = std(log10(abs(Upwp_est_avg(Uk_inds_avg,j)))); % standard deviation of the log mean - in second pass

        % Also calculate the fit slope (using unweighted PSD) -- used for diagnostics
        Uact_fit_pit(:,j) = polyfit(log10(f_nth_dec(Uk_inds_pit{j})),log10(Uppsdw_avg((Uk_inds_pit{j}),j).*((f_nth_dec(Uk_inds_pit{j}))'.^(-5/3))),1);
    
        % Compute the mean parameters
        mUx(j) = sum(abs(Utrj(floor(start_index)+1:floor(stop_index))))/(Nrec);
        mEx(j) = sum(Extrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGxUx(j) = sum(GxUxtrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGyUx(j) = sum(GyUxtrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGzUx(j) = sum(GzUxtrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGxEx(j) = sum(GxExtrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGyEx(j) = sum(GyExtrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGzEx(j) = sum(GzExtrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
    end
    % Reassign the required variables before restarting analysis for a different trajectory
    mTrUx(:,i) = mUx;
    mTrEx(:,i) = mEx;
    mTrGxUx(:,i) = mGxUx;
    mTrGyUx(:,i) = mGyUx;
    mTrGzUx(:,i) = mGzUx;
    mTrGxEx(:,i) = mGxEx;
    mTrGyEx(:,i) = mGyEx;
    mTrGzEx(:,i) = mGzEx;
    UTrlog_pwp_k_avg(:,i) = Ulog_pwp_k_avg;
    UTrlog_pwp_ks_avg(:,i) = Ulog_pwp_ks_avg;
    
    % Compute Epsilons for each trajectory
    log_Uepsilon_k(:,i) = 3/2*UTrlog_pwp_k_avg(:,i) - 3/2*log10(.146169) - log10(mTrUx(:,i));
    % compute error bars for epsilon
    log_Uepsilon_d(:,i) = 3/2*UTrlog_pwp_ks_avg(:,i);
    log_Uepsilon_u(:,i) = log_Uepsilon_k(:,i) + log_Uepsilon_d(:,i);
    log_Uepsilon_l(:,i) = log_Uepsilon_k(:,i) - log_Uepsilon_d(:,i);
    
    % Variables to store for plotting spectra
    UppdfPlt(:,:,i) = Uppsd;
    UppsfPlt(:,:,i) = Upps;
    UTrlog_ppsd_avg(:,:,i) = Ulog_ppsd_avg;
    Uusd_inds(:,i) = Uk_inds_pit;
end

% Compute trajectory averaged spectra and epsilon for Ux
for pt = 1:1:numInts
    % calculate mean airspeed
    tmpU = mTrUx(pt,:);
    avmTrUx(pt) = mean(tmpU,2);
    % average the spectra for each altitude 
    tmppt = reshape(UppdfPlt(:,pt,:),[length(freq),numTraj]);
    avUppsd(:,pt) = mean(tmppt,2);
    % average epsilon Ex
    avgTrEx = mean(mTrEx,2); 
    % conduct the spectral fit
    % weight PSD by f^(5/3)
    avUppsdw(:,pt) =  avUppsd(:,pt).*freq'.^(5/3);    % [m^2 s^(-8/3)]

    % average the weighted PSD:
    for q = 1:length(f_inds)-1
        if q == 1
            if (f_inds(q+1) - f_inds(q)) == 0
                avUppsdw_avg(q,pt) = avUppsdw(f_inds(q+1),pt);
            else
                avUppsdw_avg(q,pt) = mean( avUppsdw(f_inds(q):f_inds(q+1)-1,pt));
            end
        else
            if (f_inds(q+1) - f_inds(q)) == 0
                avUppsdw_avg(q,pt) = avUppsdw(f_inds(q+1),pt);
            else
                avUppsdw_avg(q,pt) = mean( avUppsdw(f_inds(q)+1:f_inds(q+1)-1,pt));
            end
        end
    end
    % store the averaged spectral points in this array (for comparison)
    avUPPSDW_AVG(:,pt) = avUppsdw_avg(:,pt);

    % average weighted power spectral density from f_low to f_high, [m^2 s^(-8/3)]
    avUpwp(pt) = sum(avUppsdw(f_ind_low:f_ind_high,pt))/(f_ind_high-f_ind_low+1);
    avUpwp_avg(pt) = sum(avUppsdw_avg(:,pt))/length(f_nth_dec);

    % remove data in the interval (f_ind_low:f_ind_high) at or below the sensor noise floor
    avUpit_NFw = pit_NF*freq(f_ind_low:f_ind_high).^(5/3);
    avUpit_NFw_avg = pit_NF*f_nth_dec.^(5/3);
    avUpwp_est(:,pt) = avUppsdw(f_ind_low:f_ind_high,pt);
    avUpwp_est_avg(:,pt) = avUppsdw_avg(:,pt);
    avUk_inds = find(avUpwp_est(:,pt) > avUpit_NFw');   % indexes of data points to keep
    avUk_inds_avg = find(avUpwp_est_avg(:,pt) > avUpit_NFw_avg');   % indexes of data points to keep
    avUlog_pwp_k(pt) = mean(log10(abs(avUpwp_est(avUk_inds,pt))));      % average over log of kept indexes
    avUlog_pwp_ks(pt) = std(log10(abs(avUpwp_est(avUk_inds,pt))));      % standard deviation of the log mean
    avUlog_pwp_k_avg(pt) = mean(log10(abs(avUpwp_est_avg(avUk_inds_avg,pt)))); % average over log of kept indexes
    avUlog_pwp_ks_avg(pt) = std(log10(abs(avUpwp_est_avg(avUk_inds_avg,pt)))); % standard deviation of the log mean

    % unweight these estimates for plotting and checking
    % for the average spectra
    avUlog_ppsd_est(:,pt) = avUlog_pwp_k_avg(pt) + log10(f_nth_dec.^(-5/3));
    avUlog_ppsd_avg(:,pt) = log10(avUppsdw_avg(:,pt))' + log10(f_nth_dec.^(-5/3));

    % second pass: include indexes where the fit function is above the noise floor
    avUk_inds_avg = find(avUlog_ppsd_est(:,pt) > log10(avUpit_NFw_avg.*f_nth_dec.^(-5/3))');   % indexes of data points to keep

    if pass_3 == 1
        % third pass: include the points that contribute towards reducing the
        % standard deviation
        if length(avUk_inds_avg) > min_fit_pts
            clear stand inds_to_use
            % use the first three points to begin with and then iterate over the remaining points
            avUind_to_use = avUk_inds_avg(1:min_fit_pts);
            for p = 1:1:length(avUk_inds_avg)-(min_fit_pts)+1
                avUstand(p) = std(log10(abs(avUpwp_est_avg(avUind_to_use,pt))));
                if p>1 && p<length(avUk_inds_avg)-(min_fit_pts)+1
                    if avUstand(p) > 0 && avUstand(p) < 0.25
                        if avUstand(p) > 1.2*avUstand(p-1)
                            avUind_to_use = [avUk_inds_avg(avUind_to_use(1:end-1)); avUk_inds_avg(min_fit_pts+p)];
                            avUstand(p) = avUstand(p-1);
                        else
                            avUind_to_use = [avUk_inds_avg(avUind_to_use); avUk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avUstand(p) > 0.25 && avUstand(p) < 0.5
                        if avUstand(p) > 1.15*avUstand(p-1)
                            avUind_to_use = [avUk_inds_avg(avUind_to_use(1:end-1)); avUk_inds_avg(min_fit_pts+p)];
                            avUstand(p) = avUstand(p-1);
                        else
                            avUind_to_use = [avUk_inds_avg(avUind_to_use); avUk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avUstand(p) > 0.5 && avUstand(p) < 0.75
                        if avUstand(p) > 1.12*avUstand(p-1)
                            avUind_to_use = [avUk_inds_avg(avUind_to_use(1:end-1)); avUk_inds_avg(min_fit_pts+p)];
                            avUstand(p) = avUstand(p-1);
                        else
                            avUind_to_use = [avUk_inds_avg(avUind_to_use); avUk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avUstand(p) > 0.75 && avUstand(p) < 1
                        if avUstand(p) > 1.1*avUstand(p-1)
                            avUind_to_use = [avUk_inds_avg(avUind_to_use(1:end-1)); avUk_inds_avg(min_fit_pts+p)];
                            avUstand(p) = avUstand(p-1);
                        else
                            avUind_to_use = [avUk_inds_avg(avUind_to_use); avUk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avUstand(p) > 1 && avUstand(p) < 1.5
                        if avUstand(p) > 1.05*avUstand(p-1)
                            avUind_to_use = [avUk_inds_avg(avUind_to_use(1:end-1)); avUk_inds_avg(min_fit_pts+p)];
                            avUstand(p) = avUstand(p-1);
                        else
                            avUind_to_use = [avUk_inds_avg(avUind_to_use); avUk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avUstand(p) >= 1.5
                        if avUstand(p) > 1.02*avUstand(p-1)
                            avUind_to_use = [avUk_inds_avg(avUind_to_use(1:end-1)); avUk_inds_avg(min_fit_pts+p)];
                            avUstand(p) = avUstand(p-1);
                        else
                            avUind_to_use = [avUk_inds_avg(avUind_to_use); avUk_inds_avg(min_fit_pts+p)];
                        end
                    end
                elseif p == 1
                    avUind_to_use = avUk_inds_avg(1:min_fit_pts+p);
                elseif p == length(avUk_inds_avg)-(min_fit_pts)+1
                    if avUstand(p) > 0 && avUstand(p) < 0.25
                        if avUstand(p) > 1.20*avUstand(p-1)
                            avUind_to_use = avUk_inds_avg(avUind_to_use(1:end-1));
                            avUstand(p) = avUstand(p-1);
                        else
                            avUind_to_use = avUk_inds_avg(avUind_to_use);
                        end
                    end
                    if avUstand(p) > 0.25 && avUstand(p) < 0.5
                        if avUstand(p) > 1.15*avUstand(p-1)
                            avUind_to_use = avUk_inds_avg(avUind_to_use(1:end-1));
                            avUstand(p) = avUstand(p-1);
                        else
                            avUind_to_use = avUk_inds_avg(avUind_to_use);
                        end
                    end
                    if avUstand(p) > 0.5 && avUstand(p) < 0.75
                        if avUstand(p) > 1.12*avUstand(p-1)
                            avUind_to_use = avUk_inds_avg(avUind_to_use(1:end-1));
                            avUstand(p) = avUstand(p-1);
                        else
                            avUind_to_use = avUk_inds_avg(avUind_to_use);
                        end
                    end
                    if avUstand(p) > 0.75 && avUstand(p) < 1
                        if avUstand(p) > 1.1*avUstand(p-1)
                            avUind_to_use = avUk_inds_avg(avUind_to_use(1:end-1));
                            avUstand(p) = avUstand(p-1);
                        else
                            avUind_to_use = avUk_inds_avg(avUind_to_use);
                        end
                    end
                    if avUstand(p) > 1 && avUstand(p) < 1.5
                        if avUstand(p) > 1.05*avUstand(p-1)
                            avUind_to_use = avUk_inds_avg(avUind_to_use(1:end-1));
                            avUstand(p) = avUstand(p-1);
                        else
                            avUind_to_use = avUk_inds_avg(avUind_to_use);
                        end
                    end
                    if avUstand(p) >= 1.5
                        if avUstand(p) > 1.02*avUstand(p-1)
                            avUind_to_use = avUk_inds_avg(avUind_to_use(1:end-1));
                            avUstand(p) = avUstand(p-1);
                        else
                            avUind_to_use = avUk_inds_avg(avUind_to_use);
                        end
                    end
                    avUk_inds_avg = avUind_to_use;
                end
            end
        end
    end

    if exist('stand','var')
        avUstand_pit{pt} = avUstand;
    else
        avUstand_pit{pt} = nan;
    end
    avUk_inds_pit{pt} = avUk_inds_avg;            % store the indices used for each spectral fit in a structure (the array lengths may vary)

    % if the number of selected points are less than the min number of
    % poits needed (set by the user) to get a decent fit estimate, then
    % don't bother with such spectra [set to nan by default - generally indicates that all the points in the spectra are suspect]
    avUlog_pwp_k_avg(pt) = mean(log10(abs(avUpwp_est_avg(avUk_inds_avg,pt)))); % average over log of kept indexes - in second pass
    avUlog_pwp_ks_avg(pt) = std(log10(abs(avUpwp_est_avg(avUk_inds_avg,pt)))); % standard deviation of the log mean - in second pass

    % Also calculate the fit slope (using unweighted PSD) -- used for diagnostics
    avUact_fit_pit(:,pt) = polyfit(log10(f_nth_dec(avUk_inds_pit{pt})),log10(avUppsdw_avg((avUk_inds_pit{pt}),pt).*((f_nth_dec(avUk_inds_pit{pt}))'.^(-5/3))),1);
    
    % Compute Epsilons for each trajectory
    log_avUepsilon_k(pt) = 3/2*avUlog_pwp_k_avg(pt) - 3/2*log10(.146169) - log10(avmTrUx(pt));
    % compute error bars for epsilon
    log_avUepsilon_d(pt) = 3/2*avUlog_pwp_ks_avg(pt);
    log_avUepsilon_u(pt) = log_avUepsilon_k(pt) + log_avUepsilon_d(pt);
    log_avUepsilon_l(pt) = log_avUepsilon_k(pt) - log_avUepsilon_d(pt);
end

%% Loop over V trajectories (includes averaging of TKE DNS)
for i = 1:1:numTraj
    Vtrj = Vy(:,i);
    Eytrj = Ey(:,i);
    GxVytrj = GxVy(:,i);
    GyVytrj = GyVy(:,i);
    GzVytrj = GzVy(:,i);
    GxEytrj = GxEy(:,i);
    GyEytrj = GyEy(:,i);
    GzEytrj = GzEy(:,i);
    for j = 1:1:numInts
        start_index = time_segment_center_inds(j) - Nrec/2;
        stop_index =  time_segment_center_inds(j) + Nrec/2;
        
        Vprec(:,j) = abs(Vtrj(start_index+1:stop_index));  % time records of relative wind in [m/s]
        % detrend and weight the data to reduce windowing artifacts, using a variance-preserving scale factor
        Vprecd(:,j) = detrend(Vprec(:,j)).*sqrt(2.66).*hanning(Nrec,'periodic');
    
        % fft the records to obtain the power spectral density
        Vps(:,j) = 2*fft(Vprecd(:,j))/Nrec;    % amplitude spectrum [m/s]
        Vpps(:,j) = real(conj(Vps(:,j)).*Vps(:,j));   % power spectrum [m^2/s^2]
        Vppsd(:,j) = Vpps(1:Nrec/2,j)*dt*Nrec/2;      % power spectral density up to Nyquist [m^2/s^2/Hz] = [m^2/s]
    
        % weight PSD by f^(5/3)   
        Vppsdw(:,j) =  Vppsd(:,j).*freq'.^(5/3);    % [m^2 s^(-8/3)] 
    
        % average the weighted PSD:
        for q = 1:length(f_inds)-1
            if q == 1
                if (f_inds(q+1) - f_inds(q)) == 0
                    Vppsdw_avg(q,j) = Vppsdw(f_inds(q+1),j);
                else
                    Vppsdw_avg(q,j) = mean( Vppsdw(f_inds(q):f_inds(q+1)-1,j));
                end
            else
                if (f_inds(q+1) - f_inds(q)) == 0
                    Vppsdw_avg(q,j) = Vppsdw(f_inds(q+1),j);
                else
                    Vppsdw_avg(q,j) = mean( Vppsdw(f_inds(q)+1:f_inds(q+1)-1,j));
                end
            end
        end
        % store the averaged spectral points in this array (for comparison)
        VPPSDW_AVG(:,j) = Vppsdw_avg(:,j);
   
        % average weighted power spectral density from f_low to f_high, [m^2 s^(-8/3)]
        Vpwp(j) = sum(Vppsdw(f_ind_low:f_ind_high,j))/(f_ind_high-f_ind_low+1);
        Vpwp_avg(j) = sum(Vppsdw_avg(:,j))/length(f_nth_dec);

        % remove data in the interval (f_ind_low:f_ind_high) at or below the sensor noise floor
        Vpit_NFw = pit_NF*freq(f_ind_low:f_ind_high).^(5/3);
        Vpit_NFw_avg = pit_NF*f_nth_dec.^(5/3);
        Vpwp_est(:,j) = Vppsdw(f_ind_low:f_ind_high,j);
        Vpwp_est_avg(:,j) = Vppsdw_avg(:,j);
        Vk_inds = find(Vpwp_est(:,j) > Vpit_NFw');   % indexes of data points to keep
        Vk_inds_avg = find(Vpwp_est_avg(:,j) > Vpit_NFw_avg');   % indexes of data points to keep
        Vlog_pwp_k(j) = mean(log10(abs(Vpwp_est(Vk_inds,j))));      % average over log of kept indexes
        Vlog_pwp_ks(j) = std(log10(abs(Vpwp_est(Vk_inds,j))));      % standard deviation of the log mean
        Vlog_pwp_k_avg(j) = mean(log10(abs(Vpwp_est_avg(Vk_inds_avg,j)))); % average over log of kept indexes
        Vlog_pwp_ks_avg(j) = std(log10(abs(Vpwp_est_avg(Vk_inds_avg,j)))); % standard deviation of the log mean
   
        % unweight these estimates for plotting and checking
        % for the average spectra
        Vlog_ppsd_est(:,j) = Vlog_pwp_k_avg(j) + log10(f_nth_dec.^(-5/3));
        Vlog_ppsd_avg(:,j) = log10(Vppsdw_avg(:,j))' + log10(f_nth_dec.^(-5/3));
        
        % second pass: include indexes where the fit function is above the noise floor
        Vk_inds_avg = find(Vlog_ppsd_est(:,j) > log10(Vpit_NFw_avg.*f_nth_dec.^(-5/3))');   % indexes of data points to keep
        
        if pass_3 == 1
            % third pass: include the points that contribute towards reducing the
            % standard deviation
            if length(Vk_inds_avg) > min_fit_pts
                clear stand inds_to_use
                % use the first three points to begin with and then iterate over the remaining points 
                Vind_to_use = Vk_inds_avg(1:min_fit_pts);
                for p = 1:1:length(Vk_inds_avg)-(min_fit_pts)+1
                   Vstand(p) = std(log10(abs(Vpwp_est_avg(Vind_to_use,j))));
                   if p>1 && p<length(Vk_inds_avg)-(min_fit_pts)+1
                       if Vstand(p) > 0 && Vstand(p) < 0.25
                          if Vstand(p) > 1.2*Vstand(p-1)
                             Vind_to_use = [Vk_inds_avg(Vind_to_use(1:end-1)); Vk_inds_avg(min_fit_pts+p)];
                             Vstand(p) = Vstand(p-1);
                          else
                             Vind_to_use = [Vk_inds_avg(Vind_to_use); Vk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Vstand(p) > 0.25 && Vstand(p) < 0.5
                          if Vstand(p) > 1.15*Vstand(p-1)
                             Vind_to_use = [Vk_inds_avg(Vind_to_use(1:end-1)); Vk_inds_avg(min_fit_pts+p)];
                             Vstand(p) = Vstand(p-1);
                          else
                             Vind_to_use = [Vk_inds_avg(Vind_to_use); Vk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Vstand(p) > 0.5 && Vstand(p) < 0.75
                          if Vstand(p) > 1.12*Vstand(p-1)
                             Vind_to_use = [Vk_inds_avg(Vind_to_use(1:end-1)); Vk_inds_avg(min_fit_pts+p)];
                             Vstand(p) = Vstand(p-1);
                          else
                             Vind_to_use = [Vk_inds_avg(Vind_to_use); Vk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Vstand(p) > 0.75 && Vstand(p) < 1
                          if Vstand(p) > 1.1*Vstand(p-1)
                             Vind_to_use = [Vk_inds_avg(Vind_to_use(1:end-1)); Vk_inds_avg(min_fit_pts+p)];
                             Vstand(p) = Vstand(p-1);
                          else
                             Vind_to_use = [Vk_inds_avg(Vind_to_use); Vk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Vstand(p) > 1 && Vstand(p) < 1.5
                          if Vstand(p) > 1.05*Vstand(p-1)
                             Vind_to_use = [Vk_inds_avg(Vind_to_use(1:end-1)); Vk_inds_avg(min_fit_pts+p)];
                             Vstand(p) = Vstand(p-1);
                          else
                             Vind_to_use = [Vk_inds_avg(Vind_to_use); Vk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Vstand(p) >= 1.5
                          if Vstand(p) > 1.02*Vstand(p-1)
                             Vind_to_use = [Vk_inds_avg(Vind_to_use(1:end-1)); Vk_inds_avg(min_fit_pts+p)];
                             Vstand(p) = Vstand(p-1);
                          else
                             Vind_to_use = [Vk_inds_avg(Vind_to_use); Vk_inds_avg(min_fit_pts+p)];
                          end
                       end
                   elseif p == 1
                       Vind_to_use = Vk_inds_avg(1:min_fit_pts+p);
                   elseif p == length(Vk_inds_avg)-(min_fit_pts)+1
                       if Vstand(p) > 0 && Vstand(p) < 0.25
                          if Vstand(p) > 1.20*Vstand(p-1)
                             Vind_to_use = Vk_inds_avg(Vind_to_use(1:end-1));
                             Vstand(p) = Vstand(p-1);
                          else
                             Vind_to_use = Vk_inds_avg(Vind_to_use);
                          end
                       end
                       if Vstand(p) > 0.25 && Vstand(p) < 0.5
                          if Vstand(p) > 1.15*Vstand(p-1)
                             Vind_to_use = Vk_inds_avg(Vind_to_use(1:end-1));
                             Vstand(p) = Vstand(p-1);
                          else
                             Vind_to_use = Vk_inds_avg(Vind_to_use);
                          end
                       end
                       if Vstand(p) > 0.5 && Vstand(p) < 0.75
                          if Vstand(p) > 1.12*Vstand(p-1)
                             Vind_to_use = Vk_inds_avg(Vind_to_use(1:end-1));
                             Vstand(p) = Vstand(p-1);
                          else
                             Vind_to_use = Vk_inds_avg(Vind_to_use);
                          end
                       end
                       if Vstand(p) > 0.75 && Vstand(p) < 1
                          if Vstand(p) > 1.1*Vstand(p-1)
                             Vind_to_use = Vk_inds_avg(Vind_to_use(1:end-1));
                             Vstand(p) = Vstand(p-1);
                          else
                             Vind_to_use = Vk_inds_avg(Vind_to_use);
                          end
                       end
                       if Vstand(p) > 1 && Vstand(p) < 1.5
                          if Vstand(p) > 1.05*Vstand(p-1)
                             Vind_to_use = Vk_inds_avg(Vind_to_use(1:end-1));
                             Vstand(p) = Vstand(p-1);
                          else
                             Vind_to_use = Vk_inds_avg(Vind_to_use);
                          end
                       end
                       if Vstand(p) >= 1.5
                          if Vstand(p) > 1.02*Vstand(p-1)
                             Vind_to_use = Vk_inds_avg(Vind_to_use(1:end-1));
                             Vstand(p) = Vstand(p-1);
                          else
                             Vind_to_use = Vk_inds_avg(Vind_to_use);
                          end
                       end
                       Vk_inds_avg = Vind_to_use;
                   end
                end
            end
        end
        
        if exist('stand','var')
            Vstand_pit{j} = Vstand;
        else
            Vstand_pit{j} = nan;
        end
        Vk_inds_pit{j} = Vk_inds_avg;            % store the indices used for each spectral fit in a structure (the array lengths may vary)
    
        % if the number of selected points are less than the min number of
        % poits needed (set by the user) to get a decent fit estimate, then
        % don't bother with such spectra [set to nan by default - generally indicates that all the points in the spectra are suspect]
        Vlog_pwp_k_avg(j) = mean(log10(abs(Vpwp_est_avg(Vk_inds_avg,j)))); % average over log of kept indexes - in second pass
        Vlog_pwp_ks_avg(j) = std(log10(abs(Vpwp_est_avg(Vk_inds_avg,j)))); % standard deviation of the log mean - in second pass

        % Also calculate the fit slope (using unweighted PSD) -- used for diagnostics
        Vact_fit_pit(:,j) = polyfit(log10(f_nth_dec(Vk_inds_pit{j})),log10(Vppsdw_avg((Vk_inds_pit{j}),j).*((f_nth_dec(Vk_inds_pit{j}))'.^(-5/3))),1);
    
        % Compute the mean parameters
        mVy(j) = sum(abs(Vtrj(floor(start_index)+1:floor(stop_index))))/(Nrec);
        mEy(j) = sum(Eytrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGxVy(j) = sum(GxVytrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGyVy(j) = sum(GyVytrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGzVy(j) = sum(GzVytrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGxEy(j) = sum(GxEytrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGyEy(j) = sum(GyEytrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGzEy(j) = sum(GzEytrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
    end
    % Reassign the required variables before restarting analysis for a different trajectory
    mTrVy(:,i) = mVy;
    mTrEy(:,i) = mEy;
    mTrGxVy(:,i) = mGxVy;
    mTrGyVy(:,i) = mGyVy;
    mTrGzVy(:,i) = mGzVy;
    mTrGxEy(:,i) = mGxEy;
    mTrGyEy(:,i) = mGyEy;
    mTrGzEy(:,i) = mGzEy;
    VTrlog_pwp_k_avg(:,i) = Vlog_pwp_k_avg;
    VTrlog_pwp_ks_avg(:,i) = Vlog_pwp_ks_avg;
    
    % Compute Epsilons for each trajectory
    log_Vepsilon_k(:,i) = 3/2*VTrlog_pwp_k_avg(:,i) - 3/2*log10(.146169) - log10(mTrVy(:,i));
    % compute error bars for epsilon
    log_Vepsilon_d(:,i) = 3/2*VTrlog_pwp_ks_avg(:,i);
    log_Vepsilon_u(:,i) = log_Vepsilon_k(:,i) + log_Vepsilon_d(:,i);
    log_Vepsilon_l(:,i) = log_Vepsilon_k(:,i) - log_Vepsilon_d(:,i);
    
    % Variables to store for plotting spectra
    VppdfPlt(:,:,i) = Vppsd;
    VTrlog_ppsd_avg(:,:,i) = Vlog_ppsd_avg;
    Vusd_inds(:,i) = Vk_inds_pit;
end

% Compute trajectory averaged spectra and epsilon for Vy
for pt = 1:1:numInts
    % calculate mean airspeed
    tmpV = mTrVy(pt,:);
    avmTrVy(pt) = mean(tmpV,2);
    % average the spectra for each altitude 
    tmppt = reshape(VppdfPlt(:,pt,:),[length(freq),numTraj]);
    avVppsd(:,pt) = mean(tmppt,2);
    % average epsilon Ey
    avgTrEy = mean(mTrEy,2); 
    % conduct the spectral fit
    % weight PSD by f^(5/3)
    avVppsdw(:,pt) =  avVppsd(:,pt).*freq'.^(5/3);    % [m^2 s^(-8/3)]

    % average the weighted PSD:
    for q = 1:length(f_inds)-1
        if q == 1
            if (f_inds(q+1) - f_inds(q)) == 0
                avVppsdw_avg(q,pt) = avVppsdw(f_inds(q+1),pt);
            else
                avVppsdw_avg(q,pt) = mean( avVppsdw(f_inds(q):f_inds(q+1)-1,pt));
            end
        else
            if (f_inds(q+1) - f_inds(q)) == 0
                avVppsdw_avg(q,pt) = avVppsdw(f_inds(q+1),pt);
            else
                avVppsdw_avg(q,pt) = mean( avVppsdw(f_inds(q)+1:f_inds(q+1)-1,pt));
            end
        end
    end
    % store the averaged spectral points in this array (for comparison)
    avVPPSDW_AVG(:,pt) = avVppsdw_avg(:,pt);

    % average weighted power spectral density from f_low to f_high, [m^2 s^(-8/3)]
    avVpwp(pt) = sum(avVppsdw(f_ind_low:f_ind_high,pt))/(f_ind_high-f_ind_low+1);
    avVpwp_avg(pt) = sum(avVppsdw_avg(:,pt))/length(f_nth_dec);

    % remove data in the interval (f_ind_low:f_ind_high) at or below the sensor noise floor
    avVpit_NFw = pit_NF*freq(f_ind_low:f_ind_high).^(5/3);
    avVpit_NFw_avg = pit_NF*f_nth_dec.^(5/3);
    avVpwp_est(:,pt) = avVppsdw(f_ind_low:f_ind_high,pt);
    avVpwp_est_avg(:,pt) = avVppsdw_avg(:,pt);
    avVk_inds = find(avVpwp_est(:,pt) > avVpit_NFw');   % indexes of data points to keep
    avVk_inds_avg = find(avVpwp_est_avg(:,pt) > avVpit_NFw_avg');   % indexes of data points to keep
    avVlog_pwp_k(pt) = mean(log10(abs(avVpwp_est(avVk_inds,pt))));      % average over log of kept indexes
    avVlog_pwp_ks(pt) = std(log10(abs(avVpwp_est(avVk_inds,pt))));      % standard deviation of the log mean
    avVlog_pwp_k_avg(pt) = mean(log10(abs(avVpwp_est_avg(avVk_inds_avg,pt)))); % average over log of kept indexes
    avVlog_pwp_ks_avg(pt) = std(log10(abs(avVpwp_est_avg(avVk_inds_avg,pt)))); % standard deviation of the log mean

    % unweight these estimates for plotting and checking
    % for the average spectra
    avVlog_ppsd_est(:,pt) = avVlog_pwp_k_avg(pt) + log10(f_nth_dec.^(-5/3));
    avVlog_ppsd_avg(:,pt) = log10(avVppsdw_avg(:,pt))' + log10(f_nth_dec.^(-5/3));

    % second pass: include indexes where the fit function is above the noise floor
    avVk_inds_avg = find(avVlog_ppsd_est(:,pt) > log10(avVpit_NFw_avg.*f_nth_dec.^(-5/3))');   % indexes of data points to keep

    if pass_3 == 1
        % third pass: include the points that contribute towards reducing the
        % standard deviation
        if length(avVk_inds_avg) > min_fit_pts
            clear stand inds_to_use
            % use the first three points to begin with and then iterate over the remaining points
            avVind_to_use = avVk_inds_avg(1:min_fit_pts);
            for p = 1:1:length(avVk_inds_avg)-(min_fit_pts)+1
                avVstand(p) = std(log10(abs(avVpwp_est_avg(avVind_to_use,pt))));
                if p>1 && p<length(avVk_inds_avg)-(min_fit_pts)+1
                    if avVstand(p) > 0 && avVstand(p) < 0.25
                        if avVstand(p) > 1.2*avVstand(p-1)
                            avVind_to_use = [avVk_inds_avg(avVind_to_use(1:end-1)); avVk_inds_avg(min_fit_pts+p)];
                            avVstand(p) = avVstand(p-1);
                        else
                            avVind_to_use = [avVk_inds_avg(avVind_to_use); avVk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avVstand(p) > 0.25 && avVstand(p) < 0.5
                        if avVstand(p) > 1.15*avVstand(p-1)
                            avVind_to_use = [avVk_inds_avg(avVind_to_use(1:end-1)); avVk_inds_avg(min_fit_pts+p)];
                            avVstand(p) = avVstand(p-1);
                        else
                            avVind_to_use = [avVk_inds_avg(avVind_to_use); avVk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avVstand(p) > 0.5 && avVstand(p) < 0.75
                        if avVstand(p) > 1.12*avVstand(p-1)
                            avVind_to_use = [avVk_inds_avg(avVind_to_use(1:end-1)); avVk_inds_avg(min_fit_pts+p)];
                            avVstand(p) = avVstand(p-1);
                        else
                            avVind_to_use = [avVk_inds_avg(avVind_to_use); avVk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avVstand(p) > 0.75 && avVstand(p) < 1
                        if avVstand(p) > 1.1*avVstand(p-1)
                            avVind_to_use = [avVk_inds_avg(avVind_to_use(1:end-1)); avVk_inds_avg(min_fit_pts+p)];
                            avVstand(p) = avVstand(p-1);
                        else
                            avVind_to_use = [avVk_inds_avg(avVind_to_use); avVk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avVstand(p) > 1 && avVstand(p) < 1.5
                        if avVstand(p) > 1.05*avVstand(p-1)
                            avVind_to_use = [avVk_inds_avg(avVind_to_use(1:end-1)); avVk_inds_avg(min_fit_pts+p)];
                            avVstand(p) = avVstand(p-1);
                        else
                            avVind_to_use = [avVk_inds_avg(avVind_to_use); avVk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avVstand(p) >= 1.5
                        if avVstand(p) > 1.02*avVstand(p-1)
                            avVind_to_use = [avVk_inds_avg(avVind_to_use(1:end-1)); avVk_inds_avg(min_fit_pts+p)];
                            avVstand(p) = avVstand(p-1);
                        else
                            avVind_to_use = [avVk_inds_avg(avVind_to_use); avVk_inds_avg(min_fit_pts+p)];
                        end
                    end
                elseif p == 1
                    avVind_to_use = avVk_inds_avg(1:min_fit_pts+p);
                elseif p == length(avVk_inds_avg)-(min_fit_pts)+1
                    if avVstand(p) > 0 && avVstand(p) < 0.25
                        if avVstand(p) > 1.20*avVstand(p-1)
                            avVind_to_use = avVk_inds_avg(avVind_to_use(1:end-1));
                            avVstand(p) = avVstand(p-1);
                        else
                            avVind_to_use = avVk_inds_avg(avVind_to_use);
                        end
                    end
                    if avVstand(p) > 0.25 && avVstand(p) < 0.5
                        if avVstand(p) > 1.15*avVstand(p-1)
                            avVind_to_use = avVk_inds_avg(avVind_to_use(1:end-1));
                            avVstand(p) = avVstand(p-1);
                        else
                            avVind_to_use = avVk_inds_avg(avVind_to_use);
                        end
                    end
                    if avVstand(p) > 0.5 && avVstand(p) < 0.75
                        if avVstand(p) > 1.12*avVstand(p-1)
                            avVind_to_use = avVk_inds_avg(avVind_to_use(1:end-1));
                            avVstand(p) = avVstand(p-1);
                        else
                            avVind_to_use = avVk_inds_avg(avVind_to_use);
                        end
                    end
                    if avVstand(p) > 0.75 && avVstand(p) < 1
                        if avVstand(p) > 1.1*avVstand(p-1)
                            avVind_to_use = avVk_inds_avg(avVind_to_use(1:end-1));
                            avVstand(p) = avVstand(p-1);
                        else
                            avVind_to_use = avVk_inds_avg(avVind_to_use);
                        end
                    end
                    if avVstand(p) > 1 && avVstand(p) < 1.5
                        if avVstand(p) > 1.05*avVstand(p-1)
                            avVind_to_use = avVk_inds_avg(avVind_to_use(1:end-1));
                            avVstand(p) = avVstand(p-1);
                        else
                            avVind_to_use = avVk_inds_avg(avVind_to_use);
                        end
                    end
                    if avVstand(p) >= 1.5
                        if avVstand(p) > 1.02*avVstand(p-1)
                            avVind_to_use = avVk_inds_avg(avVind_to_use(1:end-1));
                            avVstand(p) = avVstand(p-1);
                        else
                            avVind_to_use = avVk_inds_avg(avVind_to_use);
                        end
                    end
                    avVk_inds_avg = avVind_to_use;
                end
            end
        end
    end

    if exist('stand','var')
        avVstand_pit{pt} = avVstand;
    else
        avVstand_pit{pt} = nan;
    end
    avVk_inds_pit{pt} = avVk_inds_avg;            % store the indices used for each spectral fit in a structure (the array lengths may vary)

    % if the number of selected points are less than the min number of
    % poits needed (set by the user) to get a decent fit estimate, then
    % don't bother with such spectra [set to nan by default - generally indicates that all the points in the spectra are suspect]
    avVlog_pwp_k_avg(pt) = mean(log10(abs(avVpwp_est_avg(avVk_inds_avg,pt)))); % average over log of kept indexes - in second pass
    avVlog_pwp_ks_avg(pt) = std(log10(abs(avVpwp_est_avg(avVk_inds_avg,pt)))); % standard deviation of the log mean - in second pass

    % Also calculate the fit slope (using unweighted PSD) -- used for diagnostics
    avVact_fit_pit(:,pt) = polyfit(log10(f_nth_dec(avVk_inds_pit{pt})),log10(avVppsdw_avg((avVk_inds_pit{pt}),pt).*((f_nth_dec(avVk_inds_pit{pt}))'.^(-5/3))),1);
    
    % Compute Epsilons for each trajectory
    log_avVepsilon_k(pt) = 3/2*avVlog_pwp_k_avg(pt) - 3/2*log10(.146169) - log10(avmTrVy(pt));
    % compute error bars for epsilon
    log_avVepsilon_d(pt) = 3/2*avVlog_pwp_ks_avg(pt);
    log_avVepsilon_u(pt) = log_avVepsilon_k(pt) + log_avVepsilon_d(pt);
    log_avVepsilon_l(pt) = log_avVepsilon_k(pt) - log_avVepsilon_d(pt);
end

%% Loop over W trajectories (includes averaging of TKE DNS)
for i = 1:1:numTraj
    Wtrj = Wz(:,i);
    Eztrj = Ez(:,i);
    GxWztrj = GxWz(:,i);
    GyWztrj = GyWz(:,i);
    GzWztrj = GzWz(:,i);
    GxEztrj = GxEz(:,i);
    GyEztrj = GyEz(:,i);
    GzEztrj = GzEz(:,i);
    for j = 1:1:numInts
        start_index = time_segment_center_inds(j) - Nrec/2;
        stop_index =  time_segment_center_inds(j) + Nrec/2;
        
        Wprec(:,j) = abs(Wtrj(start_index+1:stop_index));  % time records of relative wind in [m/s]
        % detrend and weight the data to reduce windowing artifacts, using a variance-preserving scale factor
        Wprecd(:,j) = detrend(Wprec(:,j)).*sqrt(2.66).*hanning(Nrec,'periodic');
    
        % fft the records to obtain the power spectral density
        Wps(:,j) = 2*fft(Wprecd(:,j))/Nrec;    % amplitude spectrum [m/s]
        Wpps(:,j) = real(conj(Wps(:,j)).*Wps(:,j));   % power spectrum [m^2/s^2]
        Wppsd(:,j) = Wpps(1:Nrec/2,j)*dt*Nrec/2;      % power spectral density up to Nyquist [m^2/s^2/Hz] = [m^2/s]
    
        % weight PSD by f^(5/3)   
        Wppsdw(:,j) =  Wppsd(:,j).*freq'.^(5/3);    % [m^2 s^(-8/3)] 
    
        % average the weighted PSD:
        for q = 1:length(f_inds)-1
            if q == 1
                if (f_inds(q+1) - f_inds(q)) == 0
                    Wppsdw_avg(q,j) = Wppsdw(f_inds(q+1),j);
                else
                    Wppsdw_avg(q,j) = mean( Wppsdw(f_inds(q):f_inds(q+1)-1,j));
                end
            else
                if (f_inds(q+1) - f_inds(q)) == 0
                    Wppsdw_avg(q,j) = Wppsdw(f_inds(q+1),j);
                else
                    Wppsdw_avg(q,j) = mean( Wppsdw(f_inds(q)+1:f_inds(q+1)-1,j));
                end
            end
        end
        % store the averaged spectral points in this array (for comparison)
        WPPSDW_AVG(:,j) = Wppsdw_avg(:,j);
   
        % average weighted power spectral density from f_low to f_high, [m^2 s^(-8/3)]
        Wpwp(j) = sum(Wppsdw(f_ind_low:f_ind_high,j))/(f_ind_high-f_ind_low+1);
        Wpwp_avg(j) = sum(Wppsdw_avg(:,j))/length(f_nth_dec);

        % remove data in the interval (f_ind_low:f_ind_high) at or below the sensor noise floor
        Wpit_NFw = pit_NF*freq(f_ind_low:f_ind_high).^(5/3);
        Wpit_NFw_avg = pit_NF*f_nth_dec.^(5/3);
        Wpwp_est(:,j) = Wppsdw(f_ind_low:f_ind_high,j);
        Wpwp_est_avg(:,j) = Wppsdw_avg(:,j);
        Wk_inds = find(Wpwp_est(:,j) > Wpit_NFw');   % indexes of data points to keep
        Wk_inds_avg = find(Wpwp_est_avg(:,j) > Wpit_NFw_avg');   % indexes of data points to keep
        Wlog_pwp_k(j) = mean(log10(abs(Wpwp_est(Wk_inds,j))));      % average over log of kept indexes
        Wlog_pwp_ks(j) = std(log10(abs(Wpwp_est(Wk_inds,j))));      % standard deviation of the log mean
        Wlog_pwp_k_avg(j) = mean(log10(abs(Wpwp_est_avg(Wk_inds_avg,j)))); % average over log of kept indexes
        Wlog_pwp_ks_avg(j) = std(log10(abs(Wpwp_est_avg(Wk_inds_avg,j)))); % standard deviation of the log mean
   
        % unweight these estimates for plotting and checking
        % for the average spectra
        Wlog_ppsd_est(:,j) = Wlog_pwp_k_avg(j) + log10(f_nth_dec.^(-5/3));
        Wlog_ppsd_avg(:,j) = log10(Wppsdw_avg(:,j))' + log10(f_nth_dec.^(-5/3));
        
        % second pass: include indexes where the fit function is above the noise floor
        Wk_inds_avg = find(Wlog_ppsd_est(:,j) > log10(Wpit_NFw_avg.*f_nth_dec.^(-5/3))');   % indexes of data points to keep
        
        if pass_3 == 1
            % third pass: include the points that contribute towards reducing the
            % standard deviation
            if length(Wk_inds_avg) > min_fit_pts
                clear stand inds_to_use
                % use the first three points to begin with and then iterate over the remaining points 
                Wind_to_use = Wk_inds_avg(1:min_fit_pts);
                for p = 1:1:length(Wk_inds_avg)-(min_fit_pts)+1
                   Wstand(p) = std(log10(abs(Wpwp_est_avg(Wind_to_use,j))));
                   if p>1 && p<length(Wk_inds_avg)-(min_fit_pts)+1
                       if Wstand(p) > 0 && Wstand(p) < 0.25
                          if Wstand(p) > 1.2*Wstand(p-1)
                             Wind_to_use = [Wk_inds_avg(Wind_to_use(1:end-1)); Wk_inds_avg(min_fit_pts+p)];
                             Wstand(p) = Wstand(p-1);
                          else
                             Wind_to_use = [Wk_inds_avg(Wind_to_use); Wk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Wstand(p) > 0.25 && Wstand(p) < 0.5
                          if Wstand(p) > 1.15*Wstand(p-1)
                             Wind_to_use = [Wk_inds_avg(Wind_to_use(1:end-1)); Wk_inds_avg(min_fit_pts+p)];
                             Wstand(p) = Wstand(p-1);
                          else
                             Wind_to_use = [Wk_inds_avg(Wind_to_use); Wk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Wstand(p) > 0.5 && Wstand(p) < 0.75
                          if Wstand(p) > 1.12*Wstand(p-1)
                             Wind_to_use = [Wk_inds_avg(Wind_to_use(1:end-1)); Wk_inds_avg(min_fit_pts+p)];
                             Wstand(p) = Wstand(p-1);
                          else
                             Wind_to_use = [Wk_inds_avg(Wind_to_use); Wk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Wstand(p) > 0.75 && Wstand(p) < 1
                          if Wstand(p) > 1.1*Wstand(p-1)
                             Wind_to_use = [Wk_inds_avg(Wind_to_use(1:end-1)); Wk_inds_avg(min_fit_pts+p)];
                             Wstand(p) = Wstand(p-1);
                          else
                             Wind_to_use = [Wk_inds_avg(Wind_to_use); Wk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Wstand(p) > 1 && Wstand(p) < 1.5
                          if Wstand(p) > 1.05*Wstand(p-1)
                             Wind_to_use = [Wk_inds_avg(Wind_to_use(1:end-1)); Wk_inds_avg(min_fit_pts+p)];
                             Wstand(p) = Wstand(p-1);
                          else
                             Wind_to_use = [Wk_inds_avg(Wind_to_use); Wk_inds_avg(min_fit_pts+p)];
                          end
                       end
                       if Wstand(p) >= 1.5
                          if Wstand(p) > 1.02*Wstand(p-1)
                             Wind_to_use = [Wk_inds_avg(Wind_to_use(1:end-1)); Wk_inds_avg(min_fit_pts+p)];
                             Wstand(p) = Wstand(p-1);
                          else
                             Wind_to_use = [Wk_inds_avg(Wind_to_use); Wk_inds_avg(min_fit_pts+p)];
                          end
                       end
                   elseif p == 1
                       Wind_to_use = Wk_inds_avg(1:min_fit_pts+p);
                   elseif p == length(Wk_inds_avg)-(min_fit_pts)+1
                       if Wstand(p) > 0 && Wstand(p) < 0.25
                          if Wstand(p) > 1.20*Wstand(p-1)
                             Wind_to_use = Wk_inds_avg(Wind_to_use(1:end-1));
                             Wstand(p) = Wstand(p-1);
                          else
                             Wind_to_use = Wk_inds_avg(Wind_to_use);
                          end
                       end
                       if Wstand(p) > 0.25 && Wstand(p) < 0.5
                          if Wstand(p) > 1.15*Wstand(p-1)
                             Wind_to_use = Wk_inds_avg(Wind_to_use(1:end-1));
                             Wstand(p) = Wstand(p-1);
                          else
                             Wind_to_use = Wk_inds_avg(Wind_to_use);
                          end
                       end
                       if Wstand(p) > 0.5 && Wstand(p) < 0.75
                          if Wstand(p) > 1.12*Wstand(p-1)
                             Wind_to_use = Wk_inds_avg(Wind_to_use(1:end-1));
                             Wstand(p) = Wstand(p-1);
                          else
                             Wind_to_use = Wk_inds_avg(Wind_to_use);
                          end
                       end
                       if Wstand(p) > 0.75 && Wstand(p) < 1
                          if Wstand(p) > 1.1*Wstand(p-1)
                             Wind_to_use = Wk_inds_avg(Wind_to_use(1:end-1));
                             Wstand(p) = Wstand(p-1);
                          else
                             Wind_to_use = Wk_inds_avg(Wind_to_use);
                          end
                       end
                       if Wstand(p) > 1 && Wstand(p) < 1.5
                          if Wstand(p) > 1.05*Wstand(p-1)
                             Wind_to_use = Wk_inds_avg(Wind_to_use(1:end-1));
                             Wstand(p) = Wstand(p-1);
                          else
                             Wind_to_use = Wk_inds_avg(Wind_to_use);
                          end
                       end
                       if Wstand(p) >= 1.5
                          if Wstand(p) > 1.02*Wstand(p-1)
                             Wind_to_use = Wk_inds_avg(Wind_to_use(1:end-1));
                             Wstand(p) = Wstand(p-1);
                          else
                             Wind_to_use = Wk_inds_avg(Wind_to_use);
                          end
                       end
                       Wk_inds_avg = Wind_to_use;
                   end
                end
            end
        end
        
        if exist('stand','var')
            Wstand_pit{j} = Wstand;
        else
            Wstand_pit{j} = nan;
        end
        Wk_inds_pit{j} = Wk_inds_avg;            % store the indices used for each spectral fit in a structure (the array lengths may vary)
    
        % if the number of selected points are less than the min number of
        % poits needed (set by the user) to get a decent fit estimate, then
        % don't bother with such spectra [set to nan by default - generally indicates that all the points in the spectra are suspect]
        Wlog_pwp_k_avg(j) = mean(log10(abs(Wpwp_est_avg(Wk_inds_avg,j)))); % average over log of kept indexes - in second pass
        Wlog_pwp_ks_avg(j) = std(log10(abs(Wpwp_est_avg(Wk_inds_avg,j)))); % standard deviation of the log mean - in second pass

        % Also calculate the fit slope (using unweighted PSD) -- used for diagnostics
        Vact_fit_pit(:,j) = polyfit(log10(f_nth_dec(Wk_inds_pit{j})),log10(Wppsdw_avg((Wk_inds_pit{j}),j).*((f_nth_dec(Wk_inds_pit{j}))'.^(-5/3))),1);
    
        % Compute the mean parameters
        mWz(j) = sum(abs(Wtrj(floor(start_index)+1:floor(stop_index))))/(Nrec);
        mEz(j) = sum(Eztrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGxWz(j) = sum(GxWztrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGyWz(j) = sum(GyWztrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGzWz(j) = sum(GzWztrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGxEz(j) = sum(GxEztrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGyEz(j) = sum(GyEztrj(floor(start_index)+1:floor(stop_index)))/(Nrec);
        mGzEz(j) = sum(GzEztrj(floor(start_index)+1:floor(stop_index)))/(Nrec);

    end
    % Reassign the required variables before restarting analysis for a different trajectory
    mTrWz(:,i) = mWz;
    mTrEz(:,i) = mEz;
    mTrGxWz(:,i) = mGxWz;
    mTrGyWz(:,i) = mGyWz;
    mTrGzWz(:,i) = mGzWz;
    mTrGxEz(:,i) = mGxEz;
    mTrGyEz(:,i) = mGyEz;
    mTrGzEz(:,i) = mGzEz;
    WTrlog_pwp_k_avg(:,i) = Wlog_pwp_k_avg;
    WTrlog_pwp_ks_avg(:,i) = Wlog_pwp_ks_avg;
    
    % Compute Epsilons for each trajectory
    log_Wepsilon_k(:,i) = 3/2*WTrlog_pwp_k_avg(:,i) - 3/2*log10(.146169) - log10(mTrWz(:,i));
    % compute error bars for epsilon
    log_Wepsilon_d(:,i) = 3/2*WTrlog_pwp_ks_avg(:,i);
    log_Wepsilon_u(:,i) = log_Wepsilon_k(:,i) + log_Wepsilon_d(:,i);
    log_Wepsilon_l(:,i) = log_Wepsilon_k(:,i) - log_Wepsilon_d(:,i);
    
    % Variables to store for plotting spectra
    WppdfPlt(:,:,i) = Wppsd;
    WTrlog_ppsd_avg(:,:,i) = Wlog_ppsd_avg;
    Wusd_inds(:,i) = Wk_inds_pit;
end

% Compute trajectory averaged spectra and epsilon for Wz
for pt = 1:1:numInts
    % calculate mean airspeed
    tmpW = mTrWz(pt,:);
    avmTrWz(pt) = mean(tmpW,2);
    % average the spectra for each altitude 
    tmppt = reshape(WppdfPlt(:,pt,:),[length(freq),numTraj]);
    avWppsd(:,pt) = mean(tmppt,2);
    % average epsilon Ez
    avgTrEz = mean(mTrEz,2);
    % conduct the spectral fit
    % weight PSD by f^(5/3)
    avWppsdw(:,pt) =  avWppsd(:,pt).*freq'.^(5/3);    % [m^2 s^(-8/3)]

    % average the weighted PSD:
    for q = 1:length(f_inds)-1
        if q == 1
            if (f_inds(q+1) - f_inds(q)) == 0
                avWppsdw_avg(q,pt) = avWppsdw(f_inds(q+1),pt);
            else
                avWppsdw_avg(q,pt) = mean( avWppsdw(f_inds(q):f_inds(q+1)-1,pt));
            end
        else
            if (f_inds(q+1) - f_inds(q)) == 0
                avWppsdw_avg(q,pt) = avWppsdw(f_inds(q+1),pt);
            else
                avWppsdw_avg(q,pt) = mean( avWppsdw(f_inds(q)+1:f_inds(q+1)-1,pt));
            end
        end
    end
    % store the averaged spectral points in this array (for comparison)
    avWPPSDW_AVG(:,pt) = avWppsdw_avg(:,pt);

    % average weighted power spectral density from f_low to f_high, [m^2 s^(-8/3)]
    avWpwp(pt) = sum(avWppsdw(f_ind_low:f_ind_high,pt))/(f_ind_high-f_ind_low+1);
    avWpwp_avg(pt) = sum(avWppsdw_avg(:,pt))/length(f_nth_dec);

    % remove data in the interval (f_ind_low:f_ind_high) at or below the sensor noise floor
    avWpit_NFw = pit_NF*freq(f_ind_low:f_ind_high).^(5/3);
    avWpit_NFw_avg = pit_NF*f_nth_dec.^(5/3);
    avWpwp_est(:,pt) = avWppsdw(f_ind_low:f_ind_high,pt);
    avWpwp_est_avg(:,pt) = avWppsdw_avg(:,pt);
    avWk_inds = find(avWpwp_est(:,pt) > avWpit_NFw');   % indexes of data points to keep
    avWk_inds_avg = find(avWpwp_est_avg(:,pt) > avWpit_NFw_avg');   % indexes of data points to keep
    avWlog_pwp_k(pt) = mean(log10(abs(avWpwp_est(avWk_inds,pt))));      % average over log of kept indexes
    avWlog_pwp_ks(pt) = std(log10(abs(avWpwp_est(avWk_inds,pt))));      % standard deviation of the log mean
    avWlog_pwp_k_avg(pt) = mean(log10(abs(avWpwp_est_avg(avWk_inds_avg,pt)))); % average over log of kept indexes
    avWlog_pwp_ks_avg(pt) = std(log10(abs(avWpwp_est_avg(avWk_inds_avg,pt)))); % standard deviation of the log mean

    % unweight these estimates for plotting and checking
    % for the average spectra
    avWlog_ppsd_est(:,pt) = avWlog_pwp_k_avg(pt) + log10(f_nth_dec.^(-5/3));
    avWlog_ppsd_avg(:,pt) = log10(avWppsdw_avg(:,pt))' + log10(f_nth_dec.^(-5/3));

    % second pass: include indexes where the fit function is above the noise floor
    avWk_inds_avg = find(avWlog_ppsd_est(:,pt) > log10(avWpit_NFw_avg.*f_nth_dec.^(-5/3))');   % indexes of data points to keep

    if pass_3 == 1
        % third pass: include the points that contribute towards reducing the
        % standard deviation
        if length(avWk_inds_avg) > min_fit_pts
            clear stand inds_to_use
            % use the first three points to begin with and then iterate over the remaining points
            avWind_to_use = avWk_inds_avg(1:min_fit_pts);
            for p = 1:1:length(avWk_inds_avg)-(min_fit_pts)+1
                avWstand(p) = std(log10(abs(avWpwp_est_avg(avWind_to_use,pt))));
                if p>1 && p<length(avWk_inds_avg)-(min_fit_pts)+1
                    if avWstand(p) > 0 && avWstand(p) < 0.25
                        if avWstand(p) > 1.2*avWstand(p-1)
                            avWind_to_use = [avWk_inds_avg(avWind_to_use(1:end-1)); avWk_inds_avg(min_fit_pts+p)];
                            avWstand(p) = avWstand(p-1);
                        else
                            avWind_to_use = [avWk_inds_avg(avWind_to_use); avWk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avWstand(p) > 0.25 && avWstand(p) < 0.5
                        if avWstand(p) > 1.15*avWstand(p-1)
                            avWind_to_use = [avWk_inds_avg(avWind_to_use(1:end-1)); avWk_inds_avg(min_fit_pts+p)];
                            avWstand(p) = avWstand(p-1);
                        else
                            avWind_to_use = [avWk_inds_avg(avWind_to_use); avWk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avWstand(p) > 0.5 && avWstand(p) < 0.75
                        if avWstand(p) > 1.12*avWstand(p-1)
                            avWind_to_use = [avWk_inds_avg(avWind_to_use(1:end-1)); avWk_inds_avg(min_fit_pts+p)];
                            avWstand(p) = avWstand(p-1);
                        else
                            avWind_to_use = [avWk_inds_avg(avWind_to_use); avWk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avWstand(p) > 0.75 && avWstand(p) < 1
                        if avWstand(p) > 1.1*avWstand(p-1)
                            avWind_to_use = [avWk_inds_avg(avWind_to_use(1:end-1)); avWk_inds_avg(min_fit_pts+p)];
                            avWstand(p) = avWstand(p-1);
                        else
                            avWind_to_use = [avWk_inds_avg(avWind_to_use); avWk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avWstand(p) > 1 && avWstand(p) < 1.5
                        if avWstand(p) > 1.05*avWstand(p-1)
                            avWind_to_use = [avWk_inds_avg(avWind_to_use(1:end-1)); avWk_inds_avg(min_fit_pts+p)];
                            avWstand(p) = avWstand(p-1);
                        else
                            avWind_to_use = [avWk_inds_avg(avWind_to_use); avWk_inds_avg(min_fit_pts+p)];
                        end
                    end
                    if avWstand(p) >= 1.5
                        if avWstand(p) > 1.02*avWstand(p-1)
                            avWind_to_use = [avWk_inds_avg(avWind_to_use(1:end-1)); avWk_inds_avg(min_fit_pts+p)];
                            avWstand(p) = avWstand(p-1);
                        else
                            avWind_to_use = [avWk_inds_avg(avWind_to_use); avWk_inds_avg(min_fit_pts+p)];
                        end
                    end
                elseif p == 1
                    avWind_to_use = avWk_inds_avg(1:min_fit_pts+p);
                elseif p == length(avWk_inds_avg)-(min_fit_pts)+1
                    if avWstand(p) > 0 && avWstand(p) < 0.25
                        if avWstand(p) > 1.20*avWstand(p-1)
                            avWind_to_use = avWk_inds_avg(avWind_to_use(1:end-1));
                            avWstand(p) = avWstand(p-1);
                        else
                            avWind_to_use = avWk_inds_avg(avWind_to_use);
                        end
                    end
                    if avWstand(p) > 0.25 && avWstand(p) < 0.5
                        if avWstand(p) > 1.15*avWstand(p-1)
                            avWind_to_use = avWk_inds_avg(avWind_to_use(1:end-1));
                            avWstand(p) = avWstand(p-1);
                        else
                            avWind_to_use = avWk_inds_avg(avWind_to_use);
                        end
                    end
                    if avWstand(p) > 0.5 && avWstand(p) < 0.75
                        if avWstand(p) > 1.12*avWstand(p-1)
                            avWind_to_use = avWk_inds_avg(avWind_to_use(1:end-1));
                            avWstand(p) = avWstand(p-1);
                        else
                            avWind_to_use = avWk_inds_avg(avWind_to_use);
                        end
                    end
                    if avWstand(p) > 0.75 && avWstand(p) < 1
                        if avWstand(p) > 1.1*avWstand(p-1)
                            avWind_to_use = avWk_inds_avg(avWind_to_use(1:end-1));
                            avWstand(p) = avWstand(p-1);
                        else
                            avWind_to_use = avWk_inds_avg(avWind_to_use);
                        end
                    end
                    if avWstand(p) > 1 && avWstand(p) < 1.5
                        if avWstand(p) > 1.05*avWstand(p-1)
                            avWind_to_use = avWk_inds_avg(avWind_to_use(1:end-1));
                            avWstand(p) = avWstand(p-1);
                        else
                            avWind_to_use = avWk_inds_avg(avWind_to_use);
                        end
                    end
                    if avWstand(p) >= 1.5
                        if avWstand(p) > 1.02*avWstand(p-1)
                            avWind_to_use = avWk_inds_avg(avWind_to_use(1:end-1));
                            avWstand(p) = avWstand(p-1);
                        else
                            avWind_to_use = avWk_inds_avg(avWind_to_use);
                        end
                    end
                    avWk_inds_avg = avWind_to_use;
                end
            end
        end
    end

    if exist('stand','var')
        avWstand_pit{pt} = avWstand;
    else
        avWstand_pit{pt} = nan;
    end
    avWk_inds_pit{pt} = avWk_inds_avg;            % store the indices used for each spectral fit in a structure (the array lengths may vary)

    % if the number of selected points are less than the min number of
    % poits needed (set by the user) to get a decent fit estimate, then
    % don't bother with such spectra [set to nan by default - generally indicates that all the points in the spectra are suspect]
    avWlog_pwp_k_avg(pt) = mean(log10(abs(avWpwp_est_avg(avWk_inds_avg,pt)))); % average over log of kept indexes - in second pass
    avWlog_pwp_ks_avg(pt) = std(log10(abs(avWpwp_est_avg(avWk_inds_avg,pt)))); % standard deviation of the log mean - in second pass

    % Also calculate the fit slope (using unweighted PSD) -- used for diagnostics
    avWact_fit_pit(:,pt) = polyfit(log10(f_nth_dec(avWk_inds_pit{pt})),log10(avWppsdw_avg((avWk_inds_pit{pt}),pt).*((f_nth_dec(avWk_inds_pit{pt}))'.^(-5/3))),1);
    
    % Compute Epsilons for each trajectory
    log_avWepsilon_k(pt) = 3/2*avWlog_pwp_k_avg(pt) - 3/2*log10(.146169) - log10(avmTrWz(pt));
    % compute error bars for epsilon
    log_avWepsilon_d(pt) = 3/2*avWlog_pwp_ks_avg(pt);
    log_avWepsilon_u(pt) = log_avWepsilon_k(pt) + log_avWepsilon_d(pt);
    log_avWepsilon_l(pt) = log_avWepsilon_k(pt) - log_avWepsilon_d(pt);
end

%-log10(epsStar)
%% Save data for plotting Histogram
% save data for epsilon Ux
sz_epsU = size(log_Uepsilon_k);
save_epsDNSU = reshape(log10(mTrEx),[1,sz_epsU(1)*sz_epsU(2)]);
save_eps_U = reshape(log_Uepsilon_k,[1,sz_epsU(1)*sz_epsU(2)]);
% save data for epsilon Vy
sz_epsV = size(log_Vepsilon_k);
save_epsDNSV = reshape(log10(mTrEy),[1,sz_epsU(1)*sz_epsU(2)]);
save_eps_V = reshape(log_Vepsilon_k,[1,sz_epsV(1)*sz_epsV(2)]);
% save data for epsilon Wz
sz_epsW = size(log_Wepsilon_k);
save_epsDNSW = reshape(log10(mTrEz),[1,sz_epsU(1)*sz_epsU(2)]);
save_eps_W = reshape(log_Wepsilon_k,[1,sz_epsW(1)*sz_epsW(2)]);
if saSp == 1
    cd(dirtry)
    save('Eps_hist_data','save_eps_U','save_eps_V','save_eps_W','save_epsDNSU',...
        'save_epsDNSV','save_epsDNSW','epsStar','epsScalDNS','Ex','Ey','Ez')
    cd(exectry)
end

%% Plot the profiles
% set limits for Epsilon plotting
Epslims = [-8 0];

% Overplot all spectra for u'
for j = 1:1:numInts
    figure(1000+j)
    clf
    for i = 1:1:numTraj
        figure(1000+j)
        semilogx(freq,log10(abs(UppdfPlt(:,j,i))),'LineWidth',0.5)
        hold on
        axis([freq(1) freq(end) -15 0])
        xlabel('F')
        ylabel('Ux PSD')
        a_1=title(['All Ux Spectra for Height = ', num2str(mTrGzUx(j,i))]);
        grid on
    end
    figure(1000+j)
    semilogx(freq,avUlog_pwp_k_avg(pt)+log10(freq.^(-5/3)),'k','LineWidth',2)
    semilogx(f_nth_dec,avUlog_ppsd_avg(:,pt),'r*','LineWidth',4)
    semilogx(f_nth_dec(avUk_inds_pit{pt}),avUlog_ppsd_avg(avUk_inds_pit{pt},pt),'go','LineWidth',2)
    semilogx(freq,(avUlog_pwp_k_avg(pt)+avUlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
    semilogx(freq,(avUlog_pwp_k_avg(pt)-avUlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
    if saSp == 1
        cd(dirtry) 
        savefig(['Ux_Tr_all_Ht',num2str(mTrGzUx(j,i)),'.fig'])
        cd(exectry)
    end
    close all
end
% Overplot all spectra for v'
for j = 1:1:numInts
    figure(2000+j)
    clf
    for i = 1:1:numTraj
        figure(2000+j)
        semilogx(freq,log10(abs(VppdfPlt(:,j,i))),'LineWidth',0.5)
        hold on
        axis([freq(1) freq(end) -15 0])
        xlabel('F')
        ylabel('Vy PSD')
        a_1=title(['All Vy Spectra for Height = ', num2str(mTrGzUx(j,i))]);
        grid on
    end
    figure(2000+j)
    semilogx(freq,VTrlog_pwp_k_avg(j,i)+log10(freq.^(-5/3)),'k','LineWidth',2)
	semilogx(f_nth_dec,VTrlog_ppsd_avg(:,j,i),'r*','LineWidth',4)
    semilogx(f_nth_dec(Vusd_inds{j,i}),VTrlog_ppsd_avg(Vusd_inds{j,i},j,i),'go','LineWidth',2)
    semilogx(freq,(VTrlog_pwp_k_avg(j,i)+VTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
    semilogx(freq,(VTrlog_pwp_k_avg(j,i)-VTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)        
    if saSp == 1
        cd(dirtry) 
        savefig(['Vy_Tr_all_Ht',num2str(mTrGzUx(j,i)),'.fig'])
        cd(exectry)
    end
    close all
end
% Overplot all spectra for w'
for j = 1:1:numInts
    figure(3000+j)
    clf
    for i = 1:1:numTraj
        figure(3000+j)
        semilogx(freq,log10(abs(WppdfPlt(:,j,i))),'LineWidth',0.5)
        hold on
        axis([freq(1) freq(end) -15 0])
        xlabel('F')
        ylabel('Wz PSD')
        a_1=title(['All Wz Spectra for Height = ', num2str(mTrGzUx(j,i))]);
        grid on
    end
    figure(3000+j)
    semilogx(freq,WTrlog_pwp_k_avg(j,i)+log10(freq.^(-5/3)),'k','LineWidth',2)
	semilogx(f_nth_dec,WTrlog_ppsd_avg(:,j,i),'r*','LineWidth',4)
    semilogx(f_nth_dec(Wusd_inds{j,i}),WTrlog_ppsd_avg(Wusd_inds{j,i},j,i),'go','LineWidth',2)
	semilogx(freq,(WTrlog_pwp_k_avg(j,i)+WTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
	semilogx(freq,(WTrlog_pwp_k_avg(j,i)-WTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)        
    if saSp == 1
        cd(dirtry) 
        savefig(['Wz_Tr_all_Ht',num2str(mTrGzUx(j,i)),'.fig'])
        cd(exectry)
    end
    close all
end

% Plot averaged trajectories
figure(2000)
clf
subplot(1,5,1)
plot(log10(avgTrEx),mTrGzEx(:,1),'k','LineWidth',2)
hold on
plot(log_avUepsilon_k,mTrGzUx(:,1),'r','LineWidth',2)
plot(log_avUepsilon_u,mTrGzUx(:,1),'--r')
plot(log_avUepsilon_l,mTrGzUx(:,1),'--r')
legend('Ex','Ux')
xlabel('averaged \epsilon')
ylabel('height')
grid on
grid Minor
xlim([Epslims(1) Epslims(2)])
ylim([0 Zscal])
a_1=title(['averaged Ux and Ex comparison']);
subplot(1,5,2)
plot(log10(avgTrEy),mTrGzEy(:,1),'k','LineWidth',2)
hold on
plot(log_avVepsilon_k,mTrGzVy(:,1),'r','LineWidth',2)
plot(log_avVepsilon_u,mTrGzVy(:,1),'--r')
plot(log_avVepsilon_l,mTrGzVy(:,1),'--r')
legend('Ey','Vy')
xlabel('averaged \epsilon')
ylabel('height')
grid on
grid Minor
xlim([Epslims(1) Epslims(2)])
ylim([0 Zscal])
a_1=title(['averaged Vy and Ey comparison']);
subplot(1,5,3)
plot(log10(avgTrEz),mTrGzEz(:,1),'k','LineWidth',2)
hold on
plot(log_avWepsilon_k,mTrGzWz(:,1),'r','LineWidth',2)
plot(log_avWepsilon_u,mTrGzWz(:,1),'--r')
plot(log_avWepsilon_l,mTrGzWz(:,1),'--r')
legend('Ez','Wz')
xlabel('averaged \epsilon')
ylabel('height')
grid on
grid Minor
xlim([Epslims(1) Epslims(2)])
ylim([0 Zscal])
a_1=title('averaged Wz and Ez comparison');
subplot(1,5,4)
plot(log10(avgTrEx),mTrGzEx(:,1),'r','LineWidth',2)
hold on
plot(log10(avgTrEy),mTrGzEy(:,1),'b','LineWidth',2)
plot(log10(avgTrEz),mTrGzEz(:,1),'m','LineWidth',2)
legend('Ex','Ey','Ez')
xlabel('averaged DNS \epsilon')
ylabel('height')
grid on
grid Minor
xlim([Epslims(1) Epslims(2)])
ylim([0 Zscal])
a_1=title('averaged Ex, Ey, Ez comparison');
subplot(1,5,5)
plot(log_avUepsilon_k,mTrGzUx(:,1),'r','LineWidth',2)
hold on
plot(log_avVepsilon_k,mTrGzVy(:,1),'b','LineWidth',2)
plot(log_avWepsilon_k,mTrGzWz(:,1),'m','LineWidth',2)
legend('Ux','Vy','Wz')
xlabel('averaged Spectral \epsilon')
ylabel('height')
grid on
grid Minor
xlim([Epslims(1) Epslims(2)])
ylim([0 Zscal])
a_1=title('averaged Ux, Vy, Wz comparison');
set(gcf, 'Position',[100, 100, 1800, 900])
if saSp == 1
    cd(dirtry)
    savefig(['avgTr_comp.fig'])
    cd(exectry)
end
close all

% Plot averaged spectra
for j = 1:1:numInts
    figure(100+j)
    clf
    semilogx(freq,log10(abs(avUppsd(:,j))),'b','LineWidth',2)
    hold on
    semilogx(freq,avUlog_pwp_k_avg(pt)+log10(freq.^(-5/3)),'k','LineWidth',2)
    semilogx(f_nth_dec,avUlog_ppsd_avg(:,pt),'r*','LineWidth',4)
    semilogx(f_nth_dec(avUk_inds_pit{pt}),avUlog_ppsd_avg(avUk_inds_pit{pt},pt),'go','LineWidth',2)
    semilogx(freq,(avUlog_pwp_k_avg(pt)+avUlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
    semilogx(freq,(avUlog_pwp_k_avg(pt)-avUlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
    semilogx(freq,(avUlog_pwp_k_avg(pt)+6)+log10(freq.^(-7)),'g','LineWidth',2)
    legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','SouthWest')
    axis([freq(1) freq(end) -15 0])
    xlabel('F')
    ylabel('averaged Ux PSD')
    a_1=title(['averaged Ux Spectra, Height = ', num2str(mTrGzUx(j,i))]);
    grid on
    if saSp == 1
        cd(dirtry)
        savefig(['Ux_avg',num2str(j),'_Ht',num2str(mTrGzUx(j,i)),'.fig'])
        cd(exectry)
    end
    close all
end
for j = 1:1:numInts
    figure(200+j)
    clf
    semilogx(freq,log10(abs(avVppsd(:,j))),'b','LineWidth',2)
    hold on
    semilogx(freq,avVlog_pwp_k_avg(pt)+log10(freq.^(-5/3)),'k','LineWidth',2)
    semilogx(f_nth_dec,avVlog_ppsd_avg(:,pt),'r*','LineWidth',4)
    semilogx(f_nth_dec(avVk_inds_pit{pt}),avVlog_ppsd_avg(avVk_inds_pit{pt},pt),'go','LineWidth',2)
    semilogx(freq,(avVlog_pwp_k_avg(pt)+avVlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
    semilogx(freq,(avVlog_pwp_k_avg(pt)-avVlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
    semilogx(freq,(avVlog_pwp_k_avg(pt)+6)+log10(freq.^(-7)),'g','LineWidth',2)
    legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','SouthWest')
    axis([freq(1) freq(end) -15 0])
    xlabel('F')
    ylabel('averaged Vy PSD')
    a_1=title(['averaged Vy Spectra, Height = ', num2str(mTrGzUx(j,i))]);
    grid on
    if saSp == 1
        cd(dirtry)
        savefig(['Vy_avg',num2str(j),'_Ht',num2str(mTrGzUx(j,i)),'.fig'])
        cd(exectry)
    end
    close all
end
for j = 1:1:numInts
    figure(300+j)
    clf
    semilogx(freq,log10(abs(avWppsd(:,j))),'b','LineWidth',2)
    hold on
    semilogx(freq,avWlog_pwp_k_avg(pt)+log10(freq.^(-5/3)),'k','LineWidth',2)
    semilogx(f_nth_dec,avWlog_ppsd_avg(:,pt),'r*','LineWidth',4)
    semilogx(f_nth_dec(avWk_inds_pit{pt}),avWlog_ppsd_avg(avWk_inds_pit{pt},pt),'go','LineWidth',2)
    semilogx(freq,(avWlog_pwp_k_avg(pt)+avWlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
    semilogx(freq,(avWlog_pwp_k_avg(pt)-avWlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
    semilogx(freq,(avWlog_pwp_k_avg(pt)+6)+log10(freq.^(-7)),'g','LineWidth',2)
    legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','SouthWest')
    axis([freq(1) freq(end) -15 0])
    xlabel('F')
    ylabel('averaged Wz PSD')
    a_1=title(['averaged Wz Spectra, Height = ', num2str(mTrGzUx(j,i))]);
    grid on
    if saSp == 1
        cd(dirtry)
        savefig(['Wz_avg',num2str(j),'_Ht',num2str(mTrGzUx(j,i)),'.fig'])
        cd(exectry)
    end
    close all
end

% Plot trajectories
for i = 1:1:numTraj
    pt = i;
    figure(i)
    clf
    subplot(1,5,1)
    plot(log10(mTrEx(:,pt)),mTrGzEx(:,pt),'k','LineWidth',2)
    hold on
    plot(log_Uepsilon_k(:,pt),mTrGzUx(:,pt),'r','LineWidth',2)
    scatter(log10(Ex(:,pt)),GzEx(:,pt),'ok','MarkerEdgeAlpha',0.1)
    plot(log_Uepsilon_u(:,pt),mTrGzUx(:,pt),'--r')
    plot(log_Uepsilon_l(:,pt),mTrGzUx(:,pt),'--r')
    legend('Ex','Ux')
    xlabel('\epsilon')
    ylabel('height')
    grid on
    grid Minor
    xlim([Epslims(1) Epslims(2)])
    ylim([0 Zscal])
    a_1=title(['Ux and Ex comparison, Trajectory = ',num2str(i)]);
    subplot(1,5,2)
    plot(log10(mTrEy(:,pt)),mTrGzEy(:,pt),'k','LineWidth',2)
    hold on
    plot(log_Vepsilon_k(:,pt),mTrGzVy(:,pt),'b','LineWidth',2)
    scatter(log10(Ey(:,pt)),GzEy(:,pt),'ok','MarkerEdgeAlpha',0.1)
    plot(log_Vepsilon_u(:,pt),mTrGzUx(:,pt),'--b')
    plot(log_Vepsilon_l(:,pt),mTrGzUx(:,pt),'--b')
    legend('Ey','Vy')
    xlabel('\epsilon')
    ylabel('height')
    grid on
    grid Minor
    xlim([Epslims(1) Epslims(2)])
    ylim([0 Zscal])
    a_1=title(['Vy and Ey comparison, Trajectory = ',num2str(i)]);
    subplot(1,5,3)
    plot(log10(mTrEz(:,pt)),mTrGzEz(:,pt),'k','LineWidth',2)
    hold on
    plot(log_Wepsilon_k(:,pt),mTrGzWz(:,pt),'m','LineWidth',2)
    scatter(log10(Ez(:,pt)),GzEz(:,pt),'ok','MarkerEdgeAlpha',0.1)
    plot(log_Wepsilon_u(:,pt),mTrGzUx(:,pt),'--m')
    plot(log_Wepsilon_l(:,pt),mTrGzUx(:,pt),'--m')
    legend('Ez','Wz')
    xlabel('\epsilon')
    ylabel('height')
    grid on
    grid Minor
    xlim([Epslims(1) Epslims(2)])
    ylim([0 Zscal])
    a_1=title(['Wz and Ez comparison, Trajectory = ',num2str(i)]);
    subplot(1,5,4)
    plot(log10(mTrEx(:,pt)),mTrGzEx(:,pt),'r','LineWidth',2)
    hold on
    plot(log10(mTrEy(:,pt)),mTrGzEy(:,pt),'b','LineWidth',2)
    plot(log10(mTrEz(:,pt)),mTrGzEz(:,pt),'m','LineWidth',2)
    legend('Ex','Ey','Ez')
    xlabel('DNS \epsilon')
    ylabel('height')
    grid on
    grid Minor
    xlim([Epslims(1) Epslims(2)])
    ylim([0 Zscal])
    a_1=title(['Ex, Ey, Ez comparison, Trajectory = ',num2str(i)]);
    subplot(1,5,5)
    plot(log_Uepsilon_k(:,pt),mTrGzUx(:,pt),'r','LineWidth',2)
    hold on
    plot(log_Vepsilon_k(:,pt),mTrGzVy(:,pt),'b','LineWidth',2)
    plot(log_Wepsilon_k(:,pt),mTrGzWz(:,pt),'m','LineWidth',2)
    legend('Ux','Vy','Wz')
    xlabel('Spectral \epsilon')
    ylabel('height')
    grid on
    grid Minor
    xlim([Epslims(1) Epslims(2)])
    ylim([0 Zscal])
    a_1=title(['Ux, Vy, Wz comparison, Trajectory = ',num2str(i)]);
    set(gcf, 'Position',[100, 100, 1800, 900])
    if saSp == 1
        cd(dirtry) 
        savefig(['Tr',num2str(i),'_comp.fig'])
        cd(exectry)
    end
    close all
end

% Plot spectra
for i = 1:1:numTraj
    for j = 1:1:numInts
        figure(i)
        clf
        semilogx(freq,log10(abs(UppdfPlt(:,j,i))),'b','LineWidth',2)
        hold on
        semilogx(freq,UTrlog_pwp_k_avg(j,i)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,UTrlog_ppsd_avg(:,j,i),'r*','LineWidth',4)
        semilogx(f_nth_dec(Uusd_inds{j,i}),UTrlog_ppsd_avg(Uusd_inds{j,i},j,i),'go','LineWidth',2)
        semilogx(freq,(UTrlog_pwp_k_avg(j,i)+UTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(UTrlog_pwp_k_avg(j,i)-UTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1) 
        legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','SouthWest')
        axis([freq(1) freq(end) -15 0])
        xlabel('F')
        ylabel('Ux PSD')
        a_1=title(['Ux Spectra, Trajectory = ',num2str(i), ', Height = ', num2str(mTrGzUx(j,i))]);
        grid on
        if saSp == 1
            cd(dirtry) 
            savefig(['Ux_Tr',num2str(i),'_Ht',num2str(mTrGzUx(j,i)),'.fig'])
            cd(exectry)
        end
        close all
    end
end
for i = 1:1:numTraj
    for j = 1:1:numInts
        figure(i)
        clf
        semilogx(freq,log10(abs(VppdfPlt(:,j,i))),'b','LineWidth',2)
        hold on
        semilogx(freq,VTrlog_pwp_k_avg(j,i)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,VTrlog_ppsd_avg(:,j,i),'r*','LineWidth',4)
        semilogx(f_nth_dec(Vusd_inds{j,i}),VTrlog_ppsd_avg(Vusd_inds{j,i},j,i),'go','LineWidth',2)
        semilogx(freq,(VTrlog_pwp_k_avg(j,i)+VTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(VTrlog_pwp_k_avg(j,i)-VTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)        
        legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','SouthWest')
        axis([freq(1) freq(end) -15 0])
        xlabel('F')
        ylabel('Vy PSD')
        a_1=title(['Vy Spectra, Trajectory = ',num2str(i), ', Height = ', num2str(mTrGzUx(j,i))]);
        grid on
        if saSp == 1
            cd(dirtry) 
            savefig(['Vy_Tr',num2str(i),'_Ht',num2str(mTrGzUx(j,i)),'.fig'])
            cd(exectry)
        end
        close all
    end
end
for i = 1:1:numTraj
    for j = 1:1:numInts
        figure(i)
        clf
        semilogx(freq,log10(abs(WppdfPlt(:,j,i))),'b','LineWidth',2)
        hold on
        semilogx(freq,WTrlog_pwp_k_avg(j,i)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,WTrlog_ppsd_avg(:,j,i),'r*','LineWidth',4)
        semilogx(f_nth_dec(Wusd_inds{j,i}),WTrlog_ppsd_avg(Wusd_inds{j,i},j,i),'go','LineWidth',2)
        semilogx(freq,(WTrlog_pwp_k_avg(j,i)+WTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(WTrlog_pwp_k_avg(j,i)-WTrlog_pwp_ks_avg(j,i))+ log10(freq.^(-5/3)),'k--','LineWidth',1)        
        legend('power spectrum','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','SouthWest')
        axis([freq(1) freq(end) -15 0])
        xlabel('F')
        ylabel('Wz PSD')
        a_1=title(['Wz Spectra, Trajectory = ',num2str(i), ', Height = ', num2str(mTrGzUx(j,i))]);
        grid on
        if saSp == 1
            cd(dirtry) 
            savefig(['Wz_Tr',num2str(i),'_Ht',num2str(mTrGzUx(j,i)),'.fig'])
            cd(exectry)
        end
        close all
    end
end