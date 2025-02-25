% Compute frequency averaging and inputs for Spectral analysis
Nrec = length(specEk);
f_low_avg = 2*pi/Zl;                                        % low frequency limit for spectral averaging [Hz]
f_high_avg = (2*pi/Zl)*(Nrec);                                     % high frequency limit for spectral averaging [Hz]
pts_in_dec = 1/4;                                           % [1/3 = 1/3rd decade averaging of the spectra] [set 1/10 for 1/10th decade averaging of the spectra]
pit_NF = 10^-20;                                            % define the noise floor for pitot 
min_fit_pts = 4;                                            % the minimum number of points to be used in the third pass of fitting
pass_3 = 1;                                                 % This switch turns on an additional level scan on the spectra -- provides conservative turbulence estimates
f_low = f_low_avg;                                          % low frequency limit for spectral averaging [Hz]
f_high = f_high_avg;                                        % high frequency limit for spectral averaging [Hz]
freqE = ((2*pi)/Zl).*([1:1:Nrec]);                          % spectral frequency samples, up to Nyquist rate [Hz]
f_ind_low = find(freqE >= f_low,1,'first');
f_ind_high = find(freqE <= f_high,1,'last');
freqs_used = freqE(f_ind_low:f_ind_high);
% Setup the frequency boundaries for fractional decade averaging 
decades = log10(f_high_avg)-log10(f_low_avg);
f_avg = logspace(log10(f_low_avg),log10(f_high_avg),floor(decades/pts_in_dec));
f_nth_decE = (f_avg(2:end)+f_avg(1:end-1))/2;                % bin center frequencies
f_inds = zeros(1,length(f_avg));
for i = 1:length(f_avg)
    f_inds(i) = find(freqE-f_avg(i) >= -1e-6,1,'First');     % index of first frequency in each bin
end
f_avgind_low = freqE(f_inds(1:end-1));
f_avgind_high = freqE(f_inds(2:end)-1);
freqs_used_avg = f_nth_decE;

% weight PSD by f^(5/3)   
volEkw =  specEk.*freqE'.^(5/3);

% average the weighted PSD
for q = 1:length(f_inds)-1
    if q == 1
        if (f_inds(q+1) - f_inds(q)) == 0
            volEkw_avg(q) = volEkw(f_inds(q+1));
        else
            volEkw_avg(q) = mean(volEkw(f_inds(q):f_inds(q+1)-1));
        end
    else
        if (f_inds(q+1) - f_inds(q)) == 0
            volEkw_avg(q) = volEkw(f_inds(q+1));
        else
            volEkw_avg(q) = mean(volEkw(f_inds(q)+1:f_inds(q+1)-1));
        end
    end
end

% remove data in the interval (f_ind_low:f_ind_high) at or below the sensor noise floor
volpit_NFE = pit_NF*freqE(f_ind_low:f_ind_high).^(5/3);
volpit_NFE_avg = pit_NF*f_nth_decE.^(5/3);
volpep_est = volEkw(f_ind_low:f_ind_high);
volpep_est_avg = volEkw_avg;
volk_inds = find(volpep_est > volpit_NFE);                                  % indexes of data points to keep
volk_inds_avg = find(volpep_est_avg > volpit_NFE_avg);                      % indexes of data points to keep
vollog_pwp_k_avgE = mean(log10(abs(volpep_est_avg(volk_inds_avg))));         % average over log of kept indexes
vollog_pwp_ks_avgE = std(log10(abs(volpep_est_avg(volk_inds_avg))));         % standard deviation of the log mean
   
% unweight these estimates for plotting and checking
% for the average spectra
vollog_ppsd_estE = vollog_pwp_k_avgE + log10(f_nth_decE.^(-5/3));
vollog_ppsd_avgE = log10(volEkw_avg) + log10(f_nth_decE.^(-5/3));

% second pass: include indexes where the fit function is above the noise floor
volk_inds_avg = find(vollog_ppsd_estE > log10(volpit_NFE_avg.*f_nth_decE.^(-5/3)));       % indexes of data points to keep

if pass_3 == 1
    % third pass: include the points that contribute towards reducing the
    % standard deviation
    if length(volk_inds_avg) > min_fit_pts
        clear stand inds_to_use
        % use the first three points to begin with and then iterate over the remaining points 
        volind_to_use = volk_inds_avg(1:min_fit_pts);
        for p = 1:1:length(volk_inds_avg)-(min_fit_pts)+1
           volstand(p) = std(log10(abs(volpep_est_avg(volind_to_use))));
           if p>1 && p<length(volk_inds_avg)-(min_fit_pts)+1
               if volstand(p) > 0 && volstand(p) < 0.25
                  if volstand(p) > 1.2*volstand(p-1)
                     volind_to_use = [volk_inds_avg(volind_to_use(1:end-1)) volk_inds_avg(min_fit_pts+p)];
                     volstand(p) = volstand(p-1);
                  else
                     volind_to_use = [volk_inds_avg(volind_to_use) volk_inds_avg(min_fit_pts+p)];
                  end
               end
               if volstand(p) > 0.25 && volstand(p) < 0.5
                  if volstand(p) > 1.15*volstand(p-1)
                     volind_to_use = [volk_inds_avg(volind_to_use(1:end-1)) volk_inds_avg(min_fit_pts+p)];
                     volstand(p) = volstand(p-1);
                  else
                     volind_to_use = [volk_inds_avg(volind_to_use) volk_inds_avg(min_fit_pts+p)];
                  end
               end
               if volstand(p) > 0.5 && volstand(p) < 0.75
                  if volstand(p) > 1.12*volstand(p-1)
                     volind_to_use = [volk_inds_avg(volind_to_use(1:end-1)) volk_inds_avg(min_fit_pts+p)];
                     volstand(p) = volstand(p-1);
                  else
                     volind_to_use = [volk_inds_avg(volind_to_use) volk_inds_avg(min_fit_pts+p)];
                  end
               end
               if volstand(p) > 0.75 && volstand(p) < 1
                  if volstand(p) > 1.1*volstand(p-1)
                     volind_to_use = [volk_inds_avg(volind_to_use(1:end-1)) volk_inds_avg(min_fit_pts+p)];
                     volstand(p) = volstand(p-1);
                  else
                     volind_to_use = [volk_inds_avg(volind_to_use) volk_inds_avg(min_fit_pts+p)];
                  end
               end
               if volstand(p) > 1 && volstand(p) < 1.5
                  if volstand(p) > 1.05*volstand(p-1)
                     volind_to_use = [volk_inds_avg(volind_to_use(1:end-1)) volk_inds_avg(min_fit_pts+p)];
                     volstand(p) = volstand(p-1);
                  else
                     volind_to_use = [volk_inds_avg(volind_to_use) volk_inds_avg(min_fit_pts+p)];
                  end
               end
               if volstand(p) >= 1.5
                  if volstand(p) > 1.02*volstand(p-1)
                     volind_to_use = [volk_inds_avg(volind_to_use(1:end-1)) volk_inds_avg(min_fit_pts+p)];
                     volstand(p) = volstand(p-1);
                  else
                     volind_to_use = [volk_inds_avg(volind_to_use) volk_inds_avg(min_fit_pts+p)];
                  end
               end
           elseif p == 1
               volind_to_use = volk_inds_avg(1:min_fit_pts+p);
           elseif p == length(volk_inds_avg)-(min_fit_pts)+1
               if volstand(p) > 0 && volstand(p) < 0.25
                  if volstand(p) > 1.20*volstand(p-1)
                     volind_to_use = volk_inds_avg(volind_to_use(1:end-1));
                     volstand(p) = volstand(p-1);
                  else
                     volind_to_use = volk_inds_avg(volind_to_use);
                  end
               end
               if volstand(p) > 0.25 && volstand(p) < 0.5
                  if volstand(p) > 1.15*volstand(p-1)
                     volind_to_use = volk_inds_avg(volind_to_use(1:end-1));
                     volstand(p) = volstand(p-1);
                  else
                     volind_to_use = volk_inds_avg(volind_to_use);
                  end
               end
               if volstand(p) > 0.5 && volstand(p) < 0.75
                  if volstand(p) > 1.12*volstand(p-1)
                     volind_to_use = volk_inds_avg(volind_to_use(1:end-1));
                     volstand(p) = volstand(p-1);
                  else
                     volind_to_use = volk_inds_avg(volind_to_use);
                  end
               end
               if volstand(p) > 0.75 && volstand(p) < 1
                  if volstand(p) > 1.1*volstand(p-1)
                     volind_to_use = volk_inds_avg(volind_to_use(1:end-1));
                     volstand(p) = volstand(p-1);
                  else
                     volind_to_use = volk_inds_avg(volind_to_use);
                  end
               end
               if volstand(p) > 1 && volstand(p) < 1.5
                  if volstand(p) > 1.05*volstand(p-1)
                     volind_to_use = volk_inds_avg(volind_to_use(1:end-1));
                     volstand(p) = volstand(p-1);
                  else
                     volind_to_use = volk_inds_avg(volind_to_use);
                  end
               end
               if volstand(p) >= 1.5
                  if volstand(p) > 1.02*volstand(p-1)
                     volind_to_use = volk_inds_avg(volind_to_use(1:end-1));
                     volstand(p) = volstand(p-1);
                  else
                     volind_to_use = volk_inds_avg(volind_to_use);
                  end
               end
               volk_inds_avg = volind_to_use;
           end
        end
    end
end

% if the number of selected points are less than the min number of
% poits needed (set by the user) to get a decent fit estimate, then
% don't bother with such spectra [set to nan by default - generally indicates that all the points in the spectra are suspect]
vollog_pwp_k_avgE = mean(log10(abs(volpep_est_avg(volk_inds_avg)))); % average over log of kept indexes - in second pass
vollog_pwp_ks_avgE = std(log10(abs(volpep_est_avg(volk_inds_avg)))); % standard deviation of the log mean - in second pass

% Compute Epsilons for each trajectory: TKE here has to be calculated using
% the theoretical model (use alpha for longitudinal spectrum)
log_volepsilon_kE = 3/2*vollog_pwp_k_avgE - 3/2*log10(alp);
% compute error bars for epsilon
log_volepsilon_dE = 3/2*vollog_pwp_ks_avgE;
log_volepsilon_uE = log_volepsilon_kE + log_volepsilon_dE;
log_volepsilon_lE = log_volepsilon_kE - log_volepsilon_dE;

if diag == 1 && lean == 1
    figure(1101)
    clf
    semilogx(freqW,log10(specWWz),'b','LineWidth',2)
    hold on
    semilogx(freqU,log10(specUUz),'r','LineWidth',2)
    semilogx(freqV,log10(specVVz),'m','LineWidth',2)
    %semilogx(linspace((1/sqrt(Xl^2 + Yl^2 + Zl^2)),freq(end),length(specEk)),log10(specEk(1:end)),'k','LineWidth',2)
    semilogx(freqE,log10(specEk),'k','LineWidth',2)
    semilogx(freqW,vollog_pwp_k_avgW+log10(freqW.^(-5/3)),'b','LineWidth',2)
    semilogx(f_nth_decW,vollog_ppsd_avgW,'*r','LineWidth',4)
    semilogx(f_nth_decW(volk_inds_avg),vollog_ppsd_avgW(volk_inds_avg),'go','LineWidth',4)
    semilogx(freqW,vollog_pwp_k_avgW+vollog_pwp_ks_avgW+log10(freqW.^(-5/3)),'b--','LineWidth',1)
    semilogx(freqW,vollog_pwp_k_avgW-vollog_pwp_ks_avgW+log10(freqW.^(-5/3)),'b--','LineWidth',1)
    legend('E_{33}(\kappa_{3})','E_{11}(\kappa_{3})','E_{22}(\kappa_{3})','E(\kappa)','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','SouthWest')
    %axis([(1/sqrt(Xl^2 + Yl^2 + Zl^2)) freq(end) -15 2])
    axis([freqW(1) freqW(end) -15 2])
    xlabel('\kappa')
    ylabel('E_{33}(\kappa_{3}), E_{11}(\kappa{3}), E_{22}(\kappa{3}), E(\kappa)')
    title('averaged Wz Spectra for all records on X,Y grid points')
    grid on
    grid Minor
    if saSp == 1
        savefig('EWWKz.fig')
        close all
    end
    figure(1102)
    clf
    semilogx(freqW(1:end),log10(specWWz(1:end)),'b','LineWidth',2)
    hold on
    semilogx(freqU(1:end),log10(specUUz(1:end)),'r','LineWidth',2)
    semilogx(freqV(1:end),log10(specVVz(1:end)),'m','LineWidth',2)
    %semilogx(linspace((1/sqrt(Xl^2 + Yl^2 + Zl^2)),freq(end),length(specEk)),log10(specEk(1:end)),'k','LineWidth',2)
    semilogx(freqE,log10(specEk),'k','LineWidth',2)
    semilogx(freqU,vollog_pwp_k_avgU+log10(freqU.^(-5/3)),'r','LineWidth',2)
    semilogx(f_nth_decU,vollog_ppsd_avgU,'*r','LineWidth',4)
    semilogx(f_nth_decU(volk_inds_avg),vollog_ppsd_avgU(volk_inds_avg),'go','LineWidth',4)
    semilogx(freqU,vollog_pwp_k_avgU+vollog_pwp_ks_avgU+log10(freqU.^(-5/3)),'r--','LineWidth',1)
    semilogx(freqU,vollog_pwp_k_avgU-vollog_pwp_ks_avgU+log10(freqU.^(-5/3)),'r--','LineWidth',1)
    legend('E_{33}(\kappa_{3})','E_{11}(\kappa_{3})','E_{22}(\kappa_{3})','E(\kappa)','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','SouthWest')
    %axis([(1/sqrt(Xl^2 + Yl^2 + Zl^2)) freq(end) -15 2])
    axis([freqU(1) freqU(end) -15 2])
    xlabel('\kappa')
    ylabel('E_{33}(\kappa_{3}), E_{11}(\kappa{3}), E_{22}(\kappa{3}), E(\kappa)')
    title('averaged Uz Spectra for all records on X,Y grid points')
    grid on
    grid Minor
    if saSp == 1
        savefig('EUUKz.fig')
        close all
    end
    figure(1103)
    clf
    semilogx(freqW(1:end),log10(specWWz(1:end)),'b','LineWidth',2)
    hold on
    semilogx(freqU(1:end),log10(specUUz(1:end)),'r','LineWidth',2)
    semilogx(freqV(1:end),log10(specVVz(1:end)),'m','LineWidth',2)
    %semilogx(linspace((1/sqrt(Xl^2 + Yl^2 + Zl^2)),freq(end),length(specEk)),log10(specEk(1:end)),'k','LineWidth',2)
    semilogx(freqE,log10(specEk),'k','LineWidth',2)
    semilogx(freqV,vollog_pwp_k_avgV+log10(freqV.^(-5/3)),'m','LineWidth',2)
    semilogx(f_nth_decV,vollog_ppsd_avgV,'*r','LineWidth',4)
    semilogx(f_nth_decV(volk_inds_avg),vollog_ppsd_avgV(volk_inds_avg),'go','LineWidth',4)
    semilogx(freqV,vollog_pwp_k_avgV+vollog_pwp_ks_avgV+log10(freqV.^(-5/3)),'m--','LineWidth',1)
    semilogx(freqV,vollog_pwp_k_avgV-vollog_pwp_ks_avgV+log10(freqV.^(-5/3)),'m--','LineWidth',1)
    legend('E_{33}(\kappa_{3})','E_{11}(\kappa_{3})','E_{22}(\kappa_{3})','E(\kappa)','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','SouthWest')
    %axis([(1/sqrt(Xl^2 + Yl^2 + Zl^2)) freq(end) -15 2])
    axis([freqV(1) freqV(end) -15 2])
    xlabel('\kappa')
    ylabel('E_{33}(\kappa_{3}), E_{11}(\kappa{3}), E_{22}(\kappa{3}), E(\kappa)')
    title('averaged Vz Spectra for all records on X,Y grid points')
    grid on
    grid Minor
    if saSp == 1
        savefig('EVVKz.fig')
        close all
    end
    figure(1104)
    clf
    semilogx(freqW,log10(specWWz),'b','LineWidth',2)
    hold on
    semilogx(freqU,log10(specUUz),'r','LineWidth',2)
    semilogx(freqV,log10(specVVz),'m','LineWidth',2)
    %semilogx(linspace((1/sqrt(Xl^2 + Yl^2 + Zl^2)),freq(end),length(specEk)),log10(specEk(1:end)),'k','LineWidth',2)
    semilogx(freqE,log10(specEk),'k','LineWidth',2)
    semilogx(freqE,vollog_pwp_k_avgE+log10(freqE.^(-5/3)),'k','LineWidth',2)
    semilogx(f_nth_decE,vollog_ppsd_avgE,'*r','LineWidth',4)
    semilogx(f_nth_decE(volk_inds_avg),vollog_ppsd_avgE(volk_inds_avg),'go','LineWidth',4)
    semilogx(freqE,vollog_pwp_k_avgE+vollog_pwp_ks_avgE+log10(freqE.^(-5/3)),'k--','LineWidth',1)
    semilogx(freqE,vollog_pwp_k_avgE-vollog_pwp_ks_avgE+log10(freqE.^(-5/3)),'k--','LineWidth',1)
    legend('E_{33}(\kappa_{3})','E_{11}(\kappa_{3})','E_{22}(\kappa_{3})','E(\kappa)','-5/3 fit line','points from bin. avg spectrum','points used in fitting', 'error bars','Location','SouthWest')
    %axis([(1/sqrt(Xl^2 + Yl^2 + Zl^2)) freq(end) -15 2])
    axis([freqE(1) freqE(end) -15 2])
    xlabel('\kappa')
    ylabel('E_{33}(\kappa_{3}), E_{11}(\kappa{3}), E_{22}(\kappa{3}), E(\kappa)')
    title('E(\kappa) Spectra for the DNS domain')
    grid on
    grid Minor
    if saSp == 1
        savefig('EK.fig')
        close all
    end
end