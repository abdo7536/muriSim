%$ Script to plot comparisons of Epsilon between DNS and UAS like profiles using spectra
%% House-cleaning
clear
close all
clc

%% inputs
u_max = 1;
v_max = 1;
w_max = 1;
ftsz = 22;

%% Create file name arrays
arr1 = dir('/Users/script_away/Projects/helix/subvol1helix*.txt');
arr2 = dir('/Users/script_away/Projects/helix/subvol1vert*.txt');
arr3 = dir('/Users/script_away/Projects/helix/subvol1xyavg*.txt');
for i = 1:numel(arr1)
   hel_names{i} = arr1(i).name;
   vert_names{i} = arr2(i).name;
   xyavg_names{i} = arr3(i).name;
end
clear arr1 arr2 arr3

%% Create altitude vectors
Z0 = 0;
Zl = 39;
Z_vert = linspace(Z0,Zl,1728);
Z_xyavg = linspace(Z0,Zl,1729);
Z_helix = linspace(Z0,Zl,99840);

% find the array points
st_ht = 16.5;
en_ht = 22.5;
z_hel_st = find(Z_helix>=st_ht,1,'First')-512;
z_hel_en = find(Z_helix<=en_ht,1,'Last')+512;
z_hel_tmp = Z_helix(z_hel_st:z_hel_en);

z_xyavg_st = find(Z_xyavg>=z_hel_tmp(1),1,'First')-1;
z_xyavg_en = find(Z_xyavg<=z_hel_tmp(end),1,'Last')+1;
z_xyavg_tmp = Z_xyavg(z_xyavg_st:z_xyavg_en);

z_vert_st = find(Z_vert>=z_hel_tmp(1),1,'First')-1;
z_vert_en = find(Z_vert<=z_hel_tmp(end),1,'Last')+1;
z_vert_tmp = Z_vert(z_vert_st:z_vert_en);

%% Some calcs for spectral stuff
Nrec = 2^7;
time_segment_center_inds = [(Nrec/2):Nrec:length(z_hel_tmp)];
% time_segment_center_inds = [160*launch:Nrec:160*landing-Nrec];
dk = (1/13)/576;  % sampling modes [1/m]
wave = dk+[0:Nrec/2-1]/(Nrec*dk); % spectral frequency samples, up to Nyquist rate [Hz]

%% Load data
for i = 1:1:length(vert_names)
    % Load true vertical profile
    fname_vert = vert_names{i};
    T = table2array(readtable(fname_vert));
    eps_vert(:,i) = T(:,7);
    clear T
    % load xy averaged array
    fname_xyavg = xyavg_names{i};
    T = table2array(readtable(fname_xyavg));
    eps_xyavg(:,i) = T;
    clear T
    % load helical profiles
    fname_helix = hel_names{i};
    T = table2array(readtable(fname_helix));
    eps_helix = T(:,7);
    u = T(:,1)*u_max;
    v = T(:,2)*v_max;
    w = T(:,3)*w_max;
    clear T fname_helix fname_vert fname_xyavg

    % trim the data arrays
    u_spec_all = u(z_hel_st:z_hel_en);
    v_spec_all = v(z_hel_st:z_hel_en);
    w_spec_all = w(z_hel_st:z_hel_en);
    airspeed = sqrt((u_spec_all.^2)+(v_spec_all.^2)+(w_spec_all.^2));
    eps_helix_spec = eps_helix(z_hel_st:z_hel_en);
    eps_xyavg_all(:,i) = eps_xyavg(z_xyavg_st:z_xyavg_en,i);
    eps_vert_all(:,i) = eps_vert(z_vert_st:z_vert_en,i);

    %% Setup the frequency boundaries for fractional decade averaging 
    for j = 1:length(time_segment_center_inds)-1
        start_index = time_segment_center_inds(j) - Nrec/2;
        stop_index =  time_segment_center_inds(j) + Nrec/2;  

        prec(:,j) = 12.5+airspeed(start_index+1:stop_index);  % time records of relative wind in [m/s]
        as_rec(:,j) = 12.5+airspeed(start_index+1:stop_index);
        ephel(:,j) = eps_helix_spec(start_index+1:stop_index);
        % detrend and weight the data to reduce windowing artifacts, using a variance-preserving scale factor
        precd(:,j) = detrend(prec(:,j)).*sqrt(2.66).*hanning(Nrec,'periodic');

        % fft the records to obtain the power spectral density
        ps(:,j) = 2*fft(precd(:,j))/Nrec;    % amplitude spectrum
        pps(:,j) = real(conj(ps(:,j)).*ps(:,j));   % power spectrum
        ppsd(:,j) = pps(1:Nrec/2,j)*dk*Nrec/2;      % power spectral density up to Nyquist

        % weight PSD by f^(5/3)   
        ppsdw(:,j) =  ppsd(:,j).*wave'.^(5/3);    

        % average weighted power spectral density from f_low to f_high,
        pwp(j) = sum(ppsdw(1:Nrec/2,j))/(wave(Nrec/2)-wave(1)+1);

        pwp_est(:,j) = ppsdw(1:Nrec/2,j);
        log_pwp_k(j) = mean(log10(abs(pwp_est(1:Nrec/2,j))));      % average over log of kept indexes
        log_pwp_ks(j) = std(log10(abs(pwp_est(1:Nrec/2,j))));      % standard deviation of the log mean

        % mean airspeed over the data record [m/s]
        mean_v(j) = sum(as_rec(:,j))/Nrec;

        % mean altitude over the data record [m AGL]
        if i == 1
            alt_E(j) = sum(z_hel_tmp(floor(start_index)+1:floor(stop_index)))/(Nrec);
        end
        
        % Also calculate the fit slope (using unweighted PSD) -- used for diagnostics
        act_fit_pit(:,j) = polyfit(log10(wave),log10(ppsdw(:,j).*((wave)'.^(-5/3))),1);

        % calculate mean epsilon
        eps_hel_mean(j) = sum(ephel(:,j))/Nrec;
    end

    %% calculate the dissipation rate from spectral estimates
    log_epsilon_k(:,i) = 3/2*log_pwp_k + 3/2*log10(mean_v) - 3/2*log10(0.1615); % log10  [m^(2*3/2-1) s^(-8/3*3/2+1)] or [m^2 s^(-3)]
    eps_k(:,i) = 10.^(log_epsilon_k(:,i));
    eps_spec_helix(:,i) = eps_hel_mean;
    alt_vec = alt_E;
    mean_velocity(:,i) = mean_v;
    if mod(i,200)==0
        break
    end
end

%% Calculate Means
mean_log_epsilon_k = mean(log_epsilon_k,2);
mean_eps_k = mean(eps_k,2);
mean_eps_spec_helix = mean(eps_spec_helix,2);
mean_eps_vert = mean(eps_vert_all,2);
mean_eps_xyavg = mean(eps_xyavg_all,2);

%% Plot figure
figure(1)
clf
plot(mean_log_epsilon_k,alt_vec,'k')
hold on 
plot(log10(mean_eps_spec_helix),alt_vec,'r')
plot(log10(mean_eps_vert),z_vert_tmp,'b')
plot(log10(mean_eps_xyavg),z_xyavg_tmp,'g')
grid on
grid Minor
xlabel('log_{10}(\epsilon)')
ylabel('Z')
legend('Mean Spectral Estimate - Helix', 'Mean of Average estimate - Helix','Mean of True vertical path','Mean of XY plane averaged')
set(gca,'FontSize',ftsz)
set(gcf, 'Position',[100, 100, 600, 900])
savefig('mean_eps_prof_200.fig')
saveas(gcf,'mean_eps_prof_200.png')