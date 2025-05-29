%% Loop over V trajectories (includes averaging of TKE DNS)
for i = 1:1:numTraj
    Vtrj = Vz(:,i);
    for j = 1:1:(numInts/Trec)
        % initialize arrays
        k_inds_avgV = nan(length(f_nth_dec),numInts/Trec);
        k_inds_pltV = nan(length(f_nth_dec),numInts/Trec);
        VSAUW_WSTDplt = nan(length(f_nth_dec),numInts/Trec);
        %unConstSlpV = nan(length(f_nth_dec),numInts/Trec);
        %unConstLvlV = nan(length(f_nth_dec),numInts/Trec);
        unConstSlpVusd = nan(length(f_nth_dec)-1,numInts/Trec);
        unConstLvlVusd = nan(length(f_nth_dec)-1,numInts/Trec);
        
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
            unConstSlpV(q,j) = tmp1(:,1);
            unConstLvlV(q,j) = mean(VPPSD_AVG(q:q+1,j),'omitnan');
        end
        
        % variables to convert bin-averaged spectra to/from weighted-non-weighted forms
        uw2w = log10(f_nth_dec.^(5/3))';         % add this to weight the spectra
        w2uw = log10(f_nth_dec.^(-5/3))';        % add this to unweight the spectra
        
        % remove data in the interval (f_ind_low:f_ind_high) at or below the sensor noise floor
        Vpit_NFw = log10(pit_NF*freq(f_ind_low:f_ind_high).^(5/3));
        Vpit_NFw_avg = log10(pit_NF) + uw2w;
        Vpwp_est(:,j) = Vppsdw(f_ind_low:f_ind_high,j);
        Vpwp_est_avg(:,j) = Vppsdw_avg(:,j);
        Vk_inds = find(Vpwp_est(:,j) > Vpit_NFw');                          % indexes of data points to keep
        Vk_inds_avg = find(Vpwp_est_avg(:,j) > Vpit_NFw_avg);               % indexes of data points to keep
        Vlog_pwp_k(j) = mean((Vpwp_est(Vk_inds,j)));                        % average over log of kept indexes
        Vlog_pwp_ks(j) = std((Vpwp_est(Vk_inds,j)));                        % standard deviation of the log mean
        Vlog_pwp_k_avg(j) = mean((Vpwp_est_avg(Vk_inds_avg,j)));            % average over log of kept indexes
        Vlog_pwp_ks_avg(j) = std((Vpwp_est_avg(Vk_inds_avg,j)));            % standard deviation of the log mean
   
        % unweight these estimates for plotting and checking
        % for the average spectra
        Vlog_ppsd_est(:,j) = Vlog_pwp_k_avg(j) + w2uw;
        Vlog_ppsd_avg(:,j) = Vppsdw_avg(:,j) + w2uw;
        
        % Spectral fitting procedure is carried out
        % Also use the custom function to pick the frequency bins to be used in fitting and the quality criterion
        [k_inds_avgV(:,j),V_selind(j),bd_slpV(j),fwerV(j),ftStdDevIterV{j},pasV(j),Vk_inds_pit{j}, ...
            VSAUW_WM(j),VSAUW_WS(j),VSAUW_WSTD(j),VSAUWfit(1:length(f_nth_dec),j),VSAUW_WSTDplt(:,j),k_inds_pltV(:,j)] = ...
            specfit2(f_nth_dec, Vppsdw_avg(:,j), Vppsdw_avg(:,j)+w2uw, Vpit_NFw(:,j)+w2uw, pass_2, min_fit_pts,k_inds_pltV(:,j),uw2w,VSAUW_WSTDplt(:,j));
    
        % if the number of selected points are less than the min number of
        % poits needed (set by the user) to get a decent fit estimate, then
        % don't bother with such spectra [set to nan by default - generally indicates that all the points in the spectra are suspect]
        Vlog_pwp_k_avg(j) = mean((Vppsdw_avg(Vk_inds_pit{j},j))); % average over log of kept indexes - in second pass
        Vlog_pwp_ks_avg(j) = std((Vppsdw_avg(Vk_inds_pit{j},j))); % standard deviation of the log mean - in second pass
        
        % compute the fit slopes and levels for non-weighted spectra (used points only)
        for q = 1:length(Vk_inds_pit{j})-1
            tmp1 = polyfit(log10(f_nth_dec(Vk_inds_pit{j}(q:q+1))),...
                VPPSD_AVG(Vk_inds_pit{j}(q:q+1),j),1);
            unConstSlpVusd(Vk_inds_pit{j}(q),j) = tmp1(:,1);           % successive fit slopes of bin-averaged spectra
            unConstLvlVusd(Vk_inds_pit{j}(q),j) = mean(VPPSD_AVG(Vk_inds_pit{j}(q:q+1),j),'omitnan');           % successive fit levels of bin-averaged spectra
        end

        % also compute the unconstrained fit slope and fit level for the
        % full frequency range of the spectra
        tmp2 = polyfit(log10(f_nth_dec(Vk_inds_avg)),VPPSD_AVG(Vk_inds_avg,j),1);
        unConstSlpVFlSp(j) = tmp2(:,1);
        unConstLvlVFlSp(j) = mean(VPPSD_AVG(Vk_inds_avg,j),'omitnan');

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
    unConstSlpVusdTraj(:,:,i) = unConstSlpVusd;
    unConstLvlVusdTraj(:,:,i) = unConstLvlVusd;
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
    avVlog_ppsd_est(:,pt) = avVlog_pwp_k_avg(pt) + w2uw;
    avVlog_ppsd_avg(:,pt) = log10(avVppsdw_avg(:,pt)) + w2uw;

    % initialize arrays
    k_inds_avgV = nan(length(f_nth_dec),numInts/Trec);
    k_inds_pltV_AVG = nan(length(f_nth_dec),numInts/Trec);
    VSAUW_WSTDplt_AVG = nan(length(f_nth_dec),numInts/Trec);
    unConstSlpV_AVG = nan(length(f_nth_dec),numInts/Trec);
    unConstLvlV_AVG = nan(length(f_nth_dec),numInts/Trec);
        
    % Spectral fitting for plane-averaged spectrum is carried out
    % Also use the custom function to pick the frequency bins to be used in fitting and the quality criterion
    [k_inds_avgV(:,pt),V_selind_AVG(j),bd_slpV_AVG(j),fwerV_AVG(j),ftStdDevIterV_AVG{j},pasV_AVG(j),avVk_inds_pit{pt}, ...
        VSAUW_WM_AVG(j),VSAUW_WS_AVG(j),VSAUW_WSTD_AVG(j),VSAUWfit_AVG(1:length(f_nth_dec),j),VSAUW_WSTDplt_AVG(:,j),k_inds_pltV_AVG(:,j)] = ...
        specfit2(f_nth_dec, Vppsdw_avg(:,j), Vppsdw_avg(:,j)+w2uw, Vpit_NFw(:,j)+w2uw, pass_2, min_fit_pts,k_inds_pltV_AVG(:,j),uw2w,VSAUW_WSTDplt_AVG(:,j));

    % if the number of selected points are less than the min number of
    % poits needed (set by the user) to get a decent fit estimate, then
    % don't bother with such spectra [set to nan by default - generally indicates that all the points in the spectra are suspect]
    avVlog_pwp_k_avg(pt) = mean(log10(avVppsdw_avg(avVk_inds_pit{pt},pt))); % average over log of kept indexes - in second pass
    avVlog_pwp_ks_avg(pt) = std(log10(avVppsdw_avg(avVk_inds_pit{pt},pt))); % standard deviation of the log mean - in second pass

    % Compute Epsilons for each trajectory
    log_avVepsilon_k(pt) = 3/2*avVlog_pwp_k_avg(pt) - 3/2*log10(alpT) - log10(avmTrWz(pt));
    % compute error bars for epsilon
    log_avVepsilon_d(pt) = 3/2*avVlog_pwp_ks_avg(pt);
    log_avVepsilon_u(pt) = log_avVepsilon_k(pt) + log_avVepsilon_d(pt);
    log_avVepsilon_l(pt) = log_avVepsilon_k(pt) - log_avVepsilon_d(pt);
end