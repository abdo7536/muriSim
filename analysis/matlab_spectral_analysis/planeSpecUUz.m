% Compute frequency averaging and inputs for Spectral analysis
Nrec = length(specUUz);
f_low_avg = 2*pi;                                              % low frequency limit for spectral averaging [Hz]
f_high_avg = 2*pi*Nrec;                                          % high frequency limit for spectral averaging [Hz]
pts_in_dec = 1/4;                                           % [1/3 = 1/3rd decade averaging of the spectra] [set 1/10 for 1/10th decade averaging of the spectra]
pit_NF = 10^-20;                                            % define the noise floor for pitot 
min_fit_pts = 4;                                            % the minimum number of points to be used in the third pass of fitting
pass_3 = 1;                                                 % This switch turns on an additional level scan on the spectra -- provides conservative turbulence estimates
f_low = f_low_avg;                                          % low frequency limit for spectral averaging [Hz]
f_high = f_high_avg;                                        % high frequency limit for spectral averaging [Hz]
freqU = 2*pi*(1:1:Nrec);                                          % spectral frequency samples, up to Nyquist rate [Hz]
f_ind_low = find(freqU >= f_low,1,'first');
f_ind_high = find(freqU <= f_high,1,'last');
freqs_used = freqU(f_ind_low:f_ind_high);
% Setup the frequency boundaries for fractional decade averaging 
decades = log10(f_high_avg)-log10(f_low_avg);
f_avg = logspace(log10(f_low_avg),log10(f_high_avg),floor(decades/pts_in_dec));
f_nth_decU = (f_avg(2:end)+f_avg(1:end-1))/2;                % bin center frequencies
f_inds = zeros(1,length(f_avg));
for i = 1:length(f_avg)
    f_inds(i) = find(freqU-f_avg(i) >= -1e-6,1,'First');     % index of first frequency in each bin
end
f_avgind_low = freqU(f_inds(1:end-1));
f_avgind_high = freqU(f_inds(2:end)-1);
freqs_used_avg = f_nth_decU;

% weight PSD by f^(5/3)   
volUzw =  specUUz.*freqU'.^(5/3);

% average the weighted PSD
for q = 1:length(f_inds)-1
    if q == 1
        if (f_inds(q+1) - f_inds(q)) == 0
            volUzw_avg(q) = volUzw(f_inds(q+1));
        else
            volUzw_avg(q) = mean(volUzw(f_inds(q):f_inds(q+1)-1));
        end
    else
        if (f_inds(q+1) - f_inds(q)) == 0
            volUzw_avg(q) = volUzw(f_inds(q+1));
        else
            volUzw_avg(q) = mean(volUzw(f_inds(q)+1:f_inds(q+1)-1));
        end
    end
end

% remove data in the interval (f_ind_low:f_ind_high) at or below the sensor noise floor
volpit_NFu = pit_NF*freqU(f_ind_low:f_ind_high).^(5/3);
volpit_NFu_avg = pit_NF*f_nth_decU.^(5/3);
volpup_est = volUzw(f_ind_low:f_ind_high);
volpup_est_avg = volUzw_avg;
volk_inds = find(volpup_est > volpit_NFu);                                  % indexes of data points to keep
volk_inds_avg = find(volpup_est_avg > volpit_NFu_avg);                      % indexes of data points to keep
vollog_pwp_k_avgU = mean(log10(abs(volpup_est_avg(volk_inds_avg))));         % average over log of kept indexes
vollog_pwp_ks_avgU = std(log10(abs(volpup_est_avg(volk_inds_avg))));         % standard deviation of the log mean
   
% unweight these estimates for plotting and checking
% for the average spectra
vollog_ppsd_estU = vollog_pwp_k_avgU + log10(f_nth_decU.^(-5/3));
vollog_ppsd_avgU = log10(volUzw_avg) + log10(f_nth_decU.^(-5/3));

% second pass: include indexes where the fit function is above the noise floor
volk_inds_avg = find(vollog_ppsd_estU > log10(volpit_NFu_avg.*f_nth_decU.^(-5/3)));       % indexes of data points to keep

if pass_3 == 1
    % third pass: include the points that contribute towards reducing the
    % standard deviation
    if length(volk_inds_avg) > min_fit_pts
        clear stand inds_to_use
        % use the first three points to begin with and then iterate over the remaining points 
        volind_to_use = volk_inds_avg(1:min_fit_pts);
        for p = 1:1:length(volk_inds_avg)-(min_fit_pts)+1
           volstand(p) = std(log10(abs(volpup_est_avg(volind_to_use))));
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
vollog_pwp_k_avgU = mean(log10(abs(volpup_est_avg(volk_inds_avg)))); % average over log of kept indexes - in second pass
vollog_pwp_ks_avgU = std(log10(abs(volpup_est_avg(volk_inds_avg)))); % standard deviation of the log mean - in second pass

% Compute Epsilons for each trajectory: TKE here has to be calculated using
% the theoretical model (use alpha for longitudinal spectrum)
log_volepsilon_kU = 3/2*vollog_pwp_k_avgU - 3/2*log10(alpT);
% compute error bars for epsilon
log_volepsilon_dU = 3/2*vollog_pwp_ks_avgU;
log_volepsilon_uU = log_volepsilon_kU + log_volepsilon_dU;
log_volepsilon_lU = log_volepsilon_kU - log_volepsilon_dU;