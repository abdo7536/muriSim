%% Loop over U trajectories (includes averaging of TKE DNS)
for i = 1:1:numTraj
    Utrj = Uz(:,i);
    for j = 1:1:(numInts/Trec)
        start_index = time_segment_center_inds(j) - Nrec/2 + 1;
        stop_index =  time_segment_center_inds(j) + Nrec/2;
        Uprec(:,j) = Utrj(start_index:stop_index);                  % time records of relative wind in [m/s]
        % detrend and weight the data to reduce windowing artifacts, using a variance-preserving scale factor
        if winOn == 1
            VarbefU = 1/length(Uprec(:,j))*(Uprec(:,j)-mean(Uprec(:,j)))'*(Uprec(:,j)-mean(Uprec(:,j)));
            Uprecd(:,j) = Uprec(:,j).*hanning(Nrec,'periodic');
            VaraftU = 1/length(Uprecd(:,j))*(Uprecd(:,j)-mean(Uprecd(:,j)))'*(Uprecd(:,j)-mean(Uprecd(:,j)));
            Uprecd(:,j) = Uprecd(:,j)*(VarbefU/VaraftU)^(1/2);
            VarCorrU = 1/length(Uprecd(:,j))*(Uprecd(:,j)-mean(Uprecd(:,j)))'*(Uprecd(:,j)-mean(Uprecd(:,j)));
        else
            Uprecd(:,j) = Uprec(:,j);
        end

        % fft the records to obtain the power spectral density
        Ups(:,j) = 2*fft(Uprecd(:,j))/Nrec;             % amplitude spectrum [m/s]
        Upps(:,j) = real(conj(Ups(:,j)).*Ups(:,j));     % power spectrum [m^2/s^2]
        Uppsd(:,j) = Upps(1:Nrec/2,j)*dt*Nrec/2;        % power spectral density up to Nyquist [m^2/s^2/Hz] = [m^2/s]

        % weight PSD by f^(5/3)  
        Uppsdw(:,j) =  Uppsd(:,j).*freq'.^(5/3);        % [m^2 s^(-8/3)] 

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
        Vact_fit_pitU(:,j) = polyfit(log10(f_nth_dec(Uk_inds_pit{j})),log10(Uppsdw_avg((Uk_inds_pit{j}),j).*((f_nth_dec(Uk_inds_pit{j}))'.^(-5/3))),1);
    
        % Compute the mean parameters
        mWz(j) = balRt;
    end
    % Reassign the required variables before restarting analysis for a different trajectory
    mTrWz(:,i) = mWz;
    UTrlog_pwp_k_avg(:,i) = Ulog_pwp_k_avg;
    UTrlog_pwp_ks_avg(:,i) = Ulog_pwp_ks_avg;
    
    % Compute Epsilons for each trajectory
    log_Uepsilon_k(:,i) = 3/2*UTrlog_pwp_k_avg(:,i) - 3/2*log10(alpT) - log10(mTrWz(:,i));
    % compute error bars for epsilon
    log_Uepsilon_d(:,i) = 3/2*UTrlog_pwp_ks_avg(:,i);
    log_Uepsilon_u(:,i) = log_Uepsilon_k(:,i) + log_Uepsilon_d(:,i);
    log_Uepsilon_l(:,i) = log_Uepsilon_k(:,i) - log_Uepsilon_d(:,i);
    
    % Variables to store for plotting spectra
    UppdfPlt(:,:,i) = Uppsd;
    UTrlog_ppsd_avg(:,:,i) = Ulog_ppsd_avg;
    Uusd_inds(:,i) = Uk_inds_pit;
end

% Compute trajectory averaged spectra and epsilon for Wz
for pt = 1:1:(numInts/Trec)
    % calculate mean airspeed
    tmpW = mTrWz(pt,:);
    avmTrWz(pt) = mean(tmpW,2);
    % average the spectra for each altitude 
    tmpptU = reshape(UppdfPlt(:,pt,:),[length(freq),numTraj]);
    avUppsd(:,pt) = mean(tmpptU,2);
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
    log_avUepsilon_k(pt) = 3/2*avUlog_pwp_k_avg(pt) - 3/2*log10(alpT) - log10(avmTrWz(pt));
    % compute error bars for epsilon
    log_avUepsilon_d(pt) = 3/2*avUlog_pwp_ks_avg(pt);
    log_avUepsilon_u(pt) = log_avUepsilon_k(pt) + log_avUepsilon_d(pt);
    log_avUepsilon_l(pt) = log_avUepsilon_k(pt) - log_avUepsilon_d(pt);
end