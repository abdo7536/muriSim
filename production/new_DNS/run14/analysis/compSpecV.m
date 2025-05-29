%% Loop over V trajectories (includes averaging of TKE DNS)
for i = 1:1:numTraj
    Vtrj = Vz(:,i);
    for j = 1:1:(numInts/Trec)
        start_index = time_segment_center_inds(j) - Nrec/2 + 1;
        stop_index =  time_segment_center_inds(j) + Nrec/2;
        Vprec(:,j) = Vtrj(start_index:stop_index);                  % time records of relative wind in [m/s]
        % detrend and weight the data to reduce windowing artifacts, using a variance-preserving scale factor
        if winOn == 1
            VarbefV = 1/length(Vprec(:,j))*(Vprec(:,j)-mean(Vprec(:,j)))'*(Vprec(:,j)-mean(Vprec(:,j)));
            Vprecd(:,j) = Vprec(:,j).*hanning(Nrec,'periodic');
            VaraftV = 1/length(Vprecd(:,j))*(Vprecd(:,j)-mean(Vprecd(:,j)))'*(Vprecd(:,j)-mean(Vprecd(:,j)));
            Vprecd(:,j) = Vprecd(:,j)*(VarbefV/VaraftV)^(1/2);
            VarCorrV = 1/length(Vprecd(:,j))*(Vprecd(:,j)-mean(Vprecd(:,j)))'*(Vprecd(:,j)-mean(Vprecd(:,j)));
        else
            Vprecd(:,j) = Vprec(:,j);
        end

        % fft the records to obtain the power spectral density
        Vps(:,j) = 2*fft(Vprecd(:,j))/Nrec;             % amplitude spectrum [m/s]
        Vpps(:,j) = real(conj(Vps(:,j)).*Vps(:,j));     % power spectrum [m^2/s^2]
        Vppsd(:,j) = Vpps(1:Nrec/2,j)*dt*Nrec/2;        % power spectral density up to Nyquist [m^2/s^2/Hz] = [m^2/s]

        % weight PSD by f^(5/3)  
        Vppsdw(:,j) =  log10(Vppsd(:,j).*freq'.^(5/3));        % [m^2 s^(-8/3)] 

        % average the weighted PSD:
        for q = 1:length(f_inds)-1
            if q == 1
                if (f_inds(q+1) - f_inds(q)) == 0
                    Vppsdw_avg(q,j) = Vppsdw(f_inds(q+1),j);
                    Vppsd_avg(q,j) = log10(Vppsd(f_inds(q+1),j));
                else
                    Vppsdw_avg(q,j) = mean( Vppsdw(f_inds(q):f_inds(q+1)-1,j));
                    Vppsd_avg(q,j) = mean( log10(Vppsd(f_inds(q):f_inds(q+1)-1,j)));
                end
            else
                if (f_inds(q+1) - f_inds(q)) == 0
                    Vppsdw_avg(q,j) = Vppsdw(f_inds(q+1),j);
                    Vppsd_avg(q,j) = log10(Vppsd(f_inds(q+1),j));
                else
                    Vppsdw_avg(q,j) = mean( Vppsdw(f_inds(q)+1:f_inds(q+1)-1,j));
                    Vppsd_avg(q,j) = mean( log10(Vppsd(f_inds(q)+1:f_inds(q+1)-1,j)));
                end
            end
        end
        % store the averaged spectral points in this array (for comparison)
        VPPSDW_AVG(:,j) = Vppsdw_avg(:,j);
        VPPSD_AVG(:,j) = Vppsd_avg(:,j);
   
        % compute the fit slopes and levels (for non-weighted spectra)
        for q = 1:length(f_nth_dec)-1
            tmp1 = polyfit(log10(f_nth_dec(q:q+1)),VPPSD_AVG(q:q+1,j),1);
            unConstSlpV(j,q) = tmp1(:,1);
            unConstLvlV(j,q) = tmp1(:,2);
        end
        
        % average weighted power spectral density from f_low to f_high, [m^2 s^(-8/3)]
        Vpwp(j) = sum(Vppsdw(f_ind_low:f_ind_high,j))/(f_ind_high-f_ind_low+1);
        Vpwp_avg(j) = sum(Vppsdw_avg(:,j))/length(f_nth_dec);

        % remove data in the interval (f_ind_low:f_ind_high) at or below the sensor noise floor
        Vpit_NFw = log10(pit_NF*freq(f_ind_low:f_ind_high).^(5/3));
        Vpit_NFw_avg = log10(pit_NF*f_nth_dec.^(5/3));
        Vpwp_est(:,j) = Vppsdw(f_ind_low:f_ind_high,j);
        Vpwp_est_avg(:,j) = Vppsdw_avg(:,j);
        Vk_inds = find(Vpwp_est(:,j) > Vpit_NFw');   % indexes of data points to keep
        Vk_inds_avg = find(Vpwp_est_avg(:,j) > Vpit_NFw_avg');   % indexes of data points to keep
        Vlog_pwp_k(j) = mean((Vpwp_est(Vk_inds,j)));      % average over log of kept indexes
        Vlog_pwp_ks(j) = std((Vpwp_est(Vk_inds,j)));      % standard deviation of the log mean
        Vlog_pwp_k_avg(j) = mean((Vpwp_est_avg(Vk_inds_avg,j))); % average over log of kept indexes
        Vlog_pwp_ks_avg(j) = std((Vpwp_est_avg(Vk_inds_avg,j))); % standard deviation of the log mean
        %Vlog_pwp_k(j) = mean(log10(abs(Vpwp_est(Vk_inds,j))));      % average over log of kept indexes
        %Vlog_pwp_ks(j) = std(log10(abs(Vpwp_est(Vk_inds,j))));      % standard deviation of the log mean
        %Vlog_pwp_k_avg(j) = mean(log10(abs(Vpwp_est_avg(Vk_inds_avg,j)))); % average over log of kept indexes
        %Vlog_pwp_ks_avg(j) = std(log10(abs(Vpwp_est_avg(Vk_inds_avg,j)))); % standard deviation of the log mean
   
        % unweight these estimates for plotting and checking
        % for the average spectra
        Vlog_ppsd_est(:,j) = Vlog_pwp_k_avg(j) + log10(f_nth_dec.^(-5/3));
        Vlog_ppsd_avg(:,j) = Vppsdw_avg(:,j)' + log10(f_nth_dec.^(-5/3));
        
        % second pass: include indexes where the fit function is above the noise floor
        Vk_inds_avg = find(Vlog_ppsd_est(:,j) > (Vpit_NFw_avg+log10(f_nth_dec.^(-5/3)))');   % indexes of data points to keep
        
        if pass_3 == 1
            % third pass: include the points that contribute towards reducing the
            % standard deviation
            if length(Vk_inds_avg) > min_fit_pts
                clear stand inds_to_use
                % use the first three points to begin with and then iterate over the remaining points 
                Vind_to_use = Vk_inds_avg(1:min_fit_pts);
                for p = 1:1:length(Vk_inds_avg)-(min_fit_pts)+1
                   Vstand(p) = std((Vpwp_est_avg(Vind_to_use,j)));
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
        Vlog_pwp_k_avg(j) = mean((Vpwp_est_avg(Vk_inds_avg,j))); % average over log of kept indexes - in second pass
        Vlog_pwp_ks_avg(j) = std((Vpwp_est_avg(Vk_inds_avg,j))); % standard deviation of the log mean - in second pass
        
        % also compute the unconstrained fit slope and fit level for the
        % full frequency range of the spectra
        tmp2 = polyfit(log10(f_nth_dec(Vk_inds_avg)),VPPSD_AVG(Vk_inds_avg,j),1);
        unConstSlpVFlSp(j) = tmp2(:,1);
        unConstLvlVFlSp(j) = tmp2(:,2);

        % Compute the mean parameters
        mWz(j) = balRt;
    end
    % Reassign the required variables before restarting analysis for a different trajectory
    mTrWz(:,i) = mWz;
    VTrlog_pwp_k_avg(:,i) = Vlog_pwp_k_avg;
    VTrlog_pwp_ks_avg(:,i) = Vlog_pwp_ks_avg;
    
    % Compute Epsilons for each trajectory
    log_Vepsilon_k(:,i) = 3/2*VTrlog_pwp_k_avg(:,i) - 3/2*log10(alpT) - log10(mTrWz(:,i));
    % compute error bars for epsilon
    log_Vepsilon_d(:,i) = 3/2*VTrlog_pwp_ks_avg(:,i);
    log_Vepsilon_u(:,i) = log_Vepsilon_k(:,i) + log_Vepsilon_d(:,i);
    log_Vepsilon_l(:,i) = log_Vepsilon_k(:,i) - log_Vepsilon_d(:,i);
    
    % Variables to store for plotting spectra
    VppdfPlt(:,:,i) = Vppsd;
    VTrlog_ppsd_avg(:,:,i) = Vlog_ppsd_avg;
    Vusd_inds(:,i) = Vk_inds_pit;
    unConstSlpVTraj(:,:,i) = unConstSlpV;
	unConstLvlVTraj(:,:,i) = unConstLvlV;
    unConstSlpVFlSpTraj(:,i) = unConstSlpVFlSp;
	unConstLvlVFlSpTraj(:,i) = unConstLvlVFlSp;
end

% Compute trajectory averaged spectra and epsilon for Wz
for pt = 1:1:(numInts/Trec)
    % calculate mean airspeed
    tmpW = mTrWz(pt,:);
    avmTrWz(pt) = mean(tmpW,2);
    % average the spectra for each altitude 
    tmpptV = log10(reshape(VppdfPlt(:,pt,:),[length(freq),numTraj]));
    avVppsd(:,pt) = 10.^(mean(tmpptV,2));
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
                avVstand(p) = std(log10(avVpwp_est_avg(avVind_to_use,pt)));
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
    avVlog_pwp_k_avg(pt) = mean(log10(avVpwp_est_avg(avVk_inds_avg,pt))); % average over log of kept indexes - in second pass
    avVlog_pwp_ks_avg(pt) = std(log10(avVpwp_est_avg(avVk_inds_avg,pt))); % standard deviation of the log mean - in second pass

    % Compute Epsilons for each trajectory
    log_avVepsilon_k(pt) = 3/2*avVlog_pwp_k_avg(pt) - 3/2*log10(alpT) - log10(avmTrWz(pt));
    % compute error bars for epsilon
    log_avVepsilon_d(pt) = 3/2*avVlog_pwp_ks_avg(pt);
    log_avVepsilon_u(pt) = log_avVepsilon_k(pt) + log_avVepsilon_d(pt);
    log_avVepsilon_l(pt) = log_avVepsilon_k(pt) - log_avVepsilon_d(pt);
end