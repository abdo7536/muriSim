%% Loop over W trajectories (includes averaging of TKE DNS)
for i = 1:1:numTraj
    Wtrj = Wz(:,i);
    Etrj = epsDNSz(:,i);
    for j = 1:1:(numInts/Trec)
        start_index = time_segment_center_inds(j) - Nrec/2 + 1;
        stop_index =  time_segment_center_inds(j) + Nrec/2;
        Wprec(:,j) = Wtrj(start_index:stop_index);                  % time records of relative wind in [m/s]
        Eprec(:,j) = Etrj(start_index:stop_index);                  % time records of relative wind in [m/s]
        % detrend and weight the data to reduce windowing artifacts, using a variance-preserving scale factor
        if winOn == 1
            VarbefW = 1/length(Wprec(:,j))*(Wprec(:,j)-mean(Wprec(:,j)))'*(Wprec(:,j)-mean(Wprec(:,j)));
            Wprecd(:,j) = Wprec(:,j).*hanning(Nrec,'periodic');
            VaraftW = 1/length(Wprecd(:,j))*(Wprecd(:,j)-mean(Wprecd(:,j)))'*(Wprecd(:,j)-mean(Wprecd(:,j)));
            Wprecd(:,j) = Wprecd(:,j)*(VarbefW/VaraftW)^(1/2);
            VarCorrW = 1/length(Wprecd(:,j))*(Wprecd(:,j)-mean(Wprecd(:,j)))'*(Wprecd(:,j)-mean(Wprecd(:,j)));
        else
            Wprecd(:,j) = Wprec(:,j);
        end

        % fft the records to obtain the power spectral density
        Wps(:,j) = 2*fft(Wprecd(:,j))/Nrec;             % amplitude spectrum [m/s]
        Wpps(:,j) = real(conj(Wps(:,j)).*Wps(:,j));     % power spectrum [m^2/s^2]
        Wppsd(:,j) = Wpps(1:Nrec/2,j)*dt*Nrec/2;        % power spectral density up to Nyquist [m^2/s^2/Hz] = [m^2/s]

        % weight PSD by f^(5/3)  
        Wppsdw(:,j) =  log10(Wppsd(:,j).*freq'.^(5/3));        % [m^2 s^(-8/3)] 

        % average the weighted PSD:
        for q = 1:length(f_inds)-1
            if q == 1
                if (f_inds(q+1) - f_inds(q)) == 0
                    Wppsdw_avg(q,j) = Wppsdw(f_inds(q+1),j);
                    Wppsd_avg(q,j) = log10(Wppsd(f_inds(q+1),j));
                else
                    Wppsdw_avg(q,j) = mean( Wppsdw(f_inds(q):f_inds(q+1)-1,j));
                    Wppsd_avg(q,j) = mean( log10(Wppsd(f_inds(q):f_inds(q+1)-1,j)));
                end
            else
                if (f_inds(q+1) - f_inds(q)) == 0
                    Wppsdw_avg(q,j) = Wppsdw(f_inds(q+1),j);
                    Wppsd_avg(q,j) = log10(Wppsd(f_inds(q+1),j));
                else
                    Wppsdw_avg(q,j) = mean( Wppsdw(f_inds(q)+1:f_inds(q+1)-1,j));
                    Wppsd_avg(q,j) = mean( log10(Wppsd(f_inds(q)+1:f_inds(q+1)-1,j)));
                end
            end
        end
        % store the averaged spectral points in this array (for comparison)
        WPPSDW_AVG(:,j) = Wppsdw_avg(:,j);
        WPPSD_AVG(:,j) = Wppsd_avg(:,j);
        
        % compute the fit slopes and levels (for non-weighted spectra)
        for q = 1:length(f_nth_dec)-1
            tmp1 = polyfit(log10(f_nth_dec(q:q+1)),WPPSD_AVG(q:q+1,j),1);
            unConstSlpW(j,q) = tmp1(:,1);
            unConstLvlW(j,q) = tmp1(:,2);
        end
        
        % average weighted power spectral density from f_low to f_high, [m^2 s^(-8/3)]
        Wpwp(j) = sum(Wppsdw(f_ind_low:f_ind_high,j))/(f_ind_high-f_ind_low+1);
        Wpwp_avg(j) = sum(Wppsdw_avg(:,j))/length(f_nth_dec);

        % remove data in the interval (f_ind_low:f_ind_high) at or below the sensor noise floor
        Wpit_NFw = log10(pit_NF*freq(f_ind_low:f_ind_high).^(5/3));
        Wpit_NFw_avg = log10(pit_NF*f_nth_dec.^(5/3));
        Wpwp_est(:,j) = Wppsdw(f_ind_low:f_ind_high,j);
        Wpwp_est_avg(:,j) = Wppsdw_avg(:,j);
        Wk_inds = find(Wpwp_est(:,j) > Wpit_NFw');   % indexes of data points to keep
        Wk_inds_avg = find(Wpwp_est_avg(:,j) > Wpit_NFw_avg');   % indexes of data points to keep
        Wlog_pwp_k(j) = mean((Wpwp_est(Wk_inds,j)));      % average over log of kept indexes
        Wlog_pwp_ks(j) = std((Wpwp_est(Wk_inds,j)));      % standard deviation of the log mean
        Wlog_pwp_k_avg(j) = mean((Wpwp_est_avg(Wk_inds_avg,j))); % average over log of kept indexes
        Wlog_pwp_ks_avg(j) = std((Wpwp_est_avg(Wk_inds_avg,j))); % standard deviation of the log mean
        %Wlog_pwp_k(j) = mean(log10(abs(Wpwp_est(Wk_inds,j))));      % average over log of kept indexes
        %Wlog_pwp_ks(j) = std(log10(abs(Wpwp_est(Wk_inds,j))));      % standard deviation of the log mean
        %Wlog_pwp_k_avg(j) = mean(log10(abs(Wpwp_est_avg(Wk_inds_avg,j)))); % average over log of kept indexes
        %Wlog_pwp_ks_avg(j) = std(log10(abs(Wpwp_est_avg(Wk_inds_avg,j)))); % standard deviation of the log mean
        
        % unweight these estimates for plotting and checking
        % for the average spectra
        Wlog_ppsd_est(:,j) = Wlog_pwp_k_avg(j) + log10(f_nth_dec.^(-5/3));
        Wlog_ppsd_avg(:,j) = Wppsdw_avg(:,j)' + log10(f_nth_dec.^(-5/3));
        
        % second pass: include indexes where the fit function is above the noise floor
        Wk_inds_avg = find(Wlog_ppsd_est(:,j) > (Wpit_NFw_avg + log10(f_nth_dec.^(-5/3)))');   % indexes of data points to keep
        
        if pass_3 == 1
            % third pass: include the points that contribute towards reducing the
            % standard deviation
            if length(Wk_inds_avg) > min_fit_pts
                clear stand inds_to_use
                % use the first three points to begin with and then iterate over the remaining points 
                Wind_to_use = Wk_inds_avg(1:min_fit_pts);
                for p = 1:1:length(Wk_inds_avg)-(min_fit_pts)+1
                   Wstand(p) = std((Wpwp_est_avg(Wind_to_use,j)));
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
        Wlog_pwp_k_avg(j) = mean((Wpwp_est_avg(Wk_inds_avg,j))); % average over log of kept indexes - in second pass
        Wlog_pwp_ks_avg(j) = std((Wpwp_est_avg(Wk_inds_avg,j))); % standard deviation of the log mean - in second pass

        % also compute the unconstrained fit slope and fit level for the
        % full frequency range of the spectra
        tmp2 = polyfit(log10(f_nth_dec(Wk_inds_avg)),WPPSD_AVG(Wk_inds_avg,j),1);
        unConstSlpWFlSp(j) = tmp2(:,1);
        unConstLvlWFlSp(j) = tmp2(:,2);
        
        % Compute the mean parameters
        mWz(j) = balRt;
        mEzDNS(j) = mean(Eprec(:,j));
        sEzDNS(j) = std(Eprec(:,j));
        sEzDNSU(j) = max(Eprec(:,j));
        sEzDNSL(j) = min(Eprec(:,j));
    end
    % Reassign the required variables before restarting analysis for a different trajectory
    mTrWz(:,i) = mWz;
    WTrlog_pwp_k_avg(:,i) = Wlog_pwp_k_avg;
    WTrlog_pwp_ks_avg(:,i) = Wlog_pwp_ks_avg;
    
    % Compute Epsilons for each trajectory
    mnepsDNSz(:,i) = mEzDNS;
    stdepsDNSz(:,i) = sEzDNS;
    stdepsDNSzU(:,i) = sEzDNSU;
    stdepsDNSzL(:,i) = sEzDNSL;
    log_Wepsilon_k(:,i) = 3/2*WTrlog_pwp_k_avg(:,i) - 3/2*log10(alpL) - log10(mTrWz(:,i));
    log_Wepsilon_k2(:,i) = 3/2*WTrlog_pwp_k_avg(:,i) - 3/2*log10(cf_rad) - log10(mTrWz(:,i));
    % compute error bars for epsilon
    log_Wepsilon_d(:,i) = 3/2*WTrlog_pwp_ks_avg(:,i);
    log_Wepsilon_u(:,i) = log_Wepsilon_k(:,i) + log_Wepsilon_d(:,i);
    log_Wepsilon_l(:,i) = log_Wepsilon_k(:,i) - log_Wepsilon_d(:,i);
    
    % Variables to store for plotting spectra
    WppdfPlt(:,:,i) = Wppsd;
    WTrlog_ppsd_avg(:,:,i) = Wlog_ppsd_avg;
    Wusd_inds(:,i) = Wk_inds_pit;
    unConstSlpWTraj(:,:,i) = unConstSlpW;
	unConstLvlWTraj(:,:,i) = unConstLvlW;
    unConstSlpWFlSpTraj(:,i) = unConstSlpWFlSp;
	unConstLvlWFlSpTraj(:,i) = unConstLvlWFlSp;
end

% Compute trajectory averaged spectra and epsilon for Wz
for pt = 1:1:(numInts/Trec)
    % calculate mean airspeed
    tmpW = mTrWz(pt,:);
    avmTrWz(pt) = mean(tmpW,2);
    % average the spectra for each altitude 
    tmpptW = log10(reshape(WppdfPlt(:,pt,:),[length(freq),numTraj]));
    avWppsd(:,pt) = 10.^(mean(tmpptW,2));
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
                avWstand(p) = std(log10(avWpwp_est_avg(avWind_to_use,pt)));
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
    avWlog_pwp_k_avg(pt) = mean(log10(avWpwp_est_avg(avWk_inds_avg,pt))); % average over log of kept indexes - in second pass
    avWlog_pwp_ks_avg(pt) = std(log10(avWpwp_est_avg(avWk_inds_avg,pt))); % standard deviation of the log mean - in second pass

    % Compute Epsilons for each trajectory
    log_avWepsilon_k(pt) = 3/2*avWlog_pwp_k_avg(pt) - 3/2*log10(alpL) - log10(avmTrWz(pt));
    log_avWepsilon_k2(pt) = 3/2*avWlog_pwp_k_avg(pt) - 3/2*log10(.146169) - log10(avmTrWz(pt));
    % compute error bars for epsilon
    log_avWepsilon_d(pt) = 3/2*avWlog_pwp_ks_avg(pt);
    log_avWepsilon_u(pt) = log_avWepsilon_k(pt) + log_avWepsilon_d(pt);
    log_avWepsilon_l(pt) = log_avWepsilon_k(pt) - log_avWepsilon_d(pt);
    avTrepsDNSz(pt) = mean(mnepsDNSz(:,i));
end

% overplot all u' v' and w' spectra on the sampled plane
if diag == 1 && lean == 0
    for pt = 1:1:(numInts/Trec)
        figure(100+pt)
        clf
        for i = 1:1:numTraj
            figure(100+pt)
            semilogx(freq,log10(abs(UppdfPlt(:,j,i))),'LineWidth',0.5)
            hold on
        end
        figure(100+pt)
        semilogx(freq,avUlog_pwp_k_avg(pt)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,avUlog_ppsd_avg(:,pt),'r*','LineWidth',4)
        semilogx(f_nth_dec,avUlog_ppsd_avg(:,pt),'r','LineWidth',2)
        semilogx(f_nth_dec(avUk_inds_pit{pt}),avUlog_ppsd_avg(avUk_inds_pit{pt},pt),'go','LineWidth',2)
        semilogx(freq,(avUlog_pwp_k_avg(pt)+avUlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(avUlog_pwp_k_avg(pt)-avUlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx([freq(1) freq(end)],[log10(pit_NF) log10(pit_NF)],'m--','LineWidth',2)
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa','Fontsize',ftSz)
        ylabel('Uz PSD','Fontsize',ftSz)
        grid on
        title('normalized wavenumber vs all Uz Spectra on YZ plane (at X/2)','Fontsize',ftSz)
        set(gca,'FontSize',ftSz)
        if saSp == 1
            savefig(['spec_all_Uz_' num2str(pt) '.fig'])
            close all
        end
        figure(200+pt)
        clf
        for i = 1:1:numTraj
            figure(200+pt)
            semilogx(freq,log10(abs(VppdfPlt(:,j,i))),'LineWidth',0.5)
            hold on
        end
        figure(200+pt)
        semilogx(freq,avVlog_pwp_k_avg(pt)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,avVlog_ppsd_avg(:,pt),'r*','LineWidth',4)
        semilogx(f_nth_dec,avVlog_ppsd_avg(:,pt),'r','LineWidth',2)
        semilogx(f_nth_dec(avVk_inds_pit{pt}),avVlog_ppsd_avg(avVk_inds_pit{pt},pt),'go','LineWidth',2)
        semilogx(freq,(avVlog_pwp_k_avg(pt)+avVlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(avVlog_pwp_k_avg(pt)-avVlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx([freq(1) freq(end)],[log10(pit_NF) log10(pit_NF)],'m--','LineWidth',2)
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa','Fontsize',ftSz)
        ylabel('Vz PSD','Fontsize',ftSz)
        grid on
        title('normalized wavenumber vs all Vz Spectra on YZ plane (at X/2)','Fontsize',ftSz)
        set(gca,'FontSize',ftSz)
        if saSp == 1
            savefig(['spec_all_Vz_' num2str(pt) '.fig'])
            close all
        end
        figure(300+pt)
        clf
        for i = 1:1:numTraj
            figure(300+pt)
            semilogx(freq,log10(abs(WppdfPlt(:,j,i))),'LineWidth',0.5)
            hold on
        end
        figure(300+pt)
        semilogx(freq,avWlog_pwp_k_avg(pt)+log10(freq.^(-5/3)),'k','LineWidth',2)
        semilogx(f_nth_dec,avWlog_ppsd_avg(:,pt),'r*','LineWidth',4)
        semilogx(f_nth_dec(avWk_inds_pit{pt}),avWlog_ppsd_avg(avWk_inds_pit{pt},pt),'go','LineWidth',2)
        semilogx(f_nth_dec,avWlog_ppsd_avg(:,pt),'r','LineWidth',3)
        semilogx(freq,(avWlog_pwp_k_avg(pt)+avWlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx(freq,(avWlog_pwp_k_avg(pt)-avWlog_pwp_ks_avg(pt))+ log10(freq.^(-5/3)),'k--','LineWidth',1)
        semilogx([freq(1) freq(end)],[log10(pit_NF) log10(pit_NF)],'m--','LineWidth',2)
        axis([freq(1) freq(end) -15 2])
        xlabel('\kappa','Fontsize',ftSz)
        ylabel('Wz PSD','Fontsize',ftSz)
        grid on
        title('normalized wavenumber vs all Wz Spectra on YZ plane (at X/2)','Fontsize',ftSz)
        set(gca,'FontSize',ftSz)
        if saSp == 1
            savefig(['spec_all_Wz_' num2str(pt) '.fig'])
            close all
        end
    end
end