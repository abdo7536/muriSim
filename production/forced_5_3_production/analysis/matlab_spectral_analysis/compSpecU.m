%% Loop over U trajectories (includes averaging of TKE DNS)
for i = 1:1:numTraj
    Utrj = Uz(:,i);
    for j = 1:1:(numInts/Trec)
        % initialize arrays
        k_inds_avgU = nan(length(f_nth_dec),numInts/Trec);
        k_inds_pltU = nan(length(f_nth_dec),numInts/Trec);
        USAUW_WSTDplt = nan(length(f_nth_dec),numInts/Trec);
        %unConstSlpU = nan(length(f_nth_dec)-1,numInts/Trec);
        %unConstLvlU = nan(length(f_nth_dec)-1,numInts/Trec);
        unConstSlpUusd = nan(length(f_nth_dec)-1,numInts/Trec);
        unConstLvlUusd = nan(length(f_nth_dec)-1,numInts/Trec);

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
        elseif winOn == 0
            Uprecd(:,j) = Uprec(:,j);
        end

        % fft the records to obtain the power spectral density
        Ups(:,j) = 2*fft(Uprecd(:,j))/Nrec;             % amplitude spectrum [m/s]
        Upps(:,j) = real(conj(Ups(:,j)).*Ups(:,j));     % power spectrum [m^2/s^2]
        Uppsd(:,j) = Upps(1:Nrec/2,j)*dt*Nrec/2;        % power spectral density up to Nyquist [m^2/s^2/Hz] = [m^2/s]

        % weight PSD by f^(5/3)  
        Uppsdw(:,j) =  log10(Uppsd(:,j).*freq'.^(5/3));        % [m^2 s^(-8/3)] 

        % average the weighted PSD:
        for q = 1:length(f_inds)-1
            if q == 1
                if (f_inds(q+1) - f_inds(q)) == 0
                    Uppsdw_avg(q,j) = Uppsdw(f_inds(q+1),j);
                    Uppsd_avg(q,j) = log10(Uppsd(f_inds(q+1),j));
                else
                    Uppsdw_avg(q,j) = mean( Uppsdw(f_inds(q):f_inds(q+1)-1,j));
                    Uppsd_avg(q,j) = mean(log10(Uppsd(f_inds(q):f_inds(q+1)-1,j)));
                end
            else
                if (f_inds(q+1) - f_inds(q)) == 0
                    Uppsdw_avg(q,j) = Uppsdw(f_inds(q+1),j);
                    Uppsd_avg(q,j) = log10(Uppsd(f_inds(q+1),j));
                else
                    Uppsdw_avg(q,j) = mean( Uppsdw(f_inds(q)+1:f_inds(q+1)-1,j));
                    Uppsd_avg(q,j) = mean(log10(Uppsd(f_inds(q)+1:f_inds(q+1)-1,j)));
                end
            end
        end
        % store the averaged spectral points in this array (for comparison)
        UPPSDW_AVG(:,j) = Uppsdw_avg(:,j);          % bin-averaged weighted spectra
        UPPSD_AVG(:,j) = Uppsd_avg(:,j);            % bin-averaged unweighted spectra
   
        % compute the fit slopes and levels for non-weighted spectra
        for q = 1:length(f_nth_dec)-1
            tmp1 = polyfit(log10(f_nth_dec(q:q+1)),UPPSD_AVG(q:q+1,j),1);
            unConstSlpU(q,j) = tmp1(:,1);           % successive fit slopes of bin-averaged spectra
            unConstLvlU(q,j) = mean(UPPSD_AVG(q:q+1,j),'omitnan');           % successive fit levels of bin-averaged spectra
        end

        % variables to convert bin-averaged spectra to/from weighted-non-weighted forms
        uw2w = log10(f_nth_dec.^(5/3))';         % add this to weight the spectra
        w2uw = log10(f_nth_dec.^(-5/3))';        % add this to unweight the spectra

        % remove data in the interval (f_ind_low:f_ind_high) at or below the sensor noise floor
        Upit_NFw = log10(pit_NF*freq(f_ind_low:f_ind_high).^(5/3));         % weighted noise-floor
        Upit_NFw_avg = log10(pit_NF) + uw2w;                                % bin-averaged weighted noise-floor
        Upwp_est(:,j) = Uppsdw(f_ind_low:f_ind_high,j);                     % all points in the user desired frequency range to be used for fitting
        Upwp_est_avg(:,j) = Uppsdw_avg(:,j);                                % all points in the user desired frequency range to be used for fitting (for bin-averaged spectra)
        Uk_inds = find(Upwp_est(:,j) > Upit_NFw');                          % indexes of data points to keep
        Uk_inds_avg = find(Upwp_est_avg(:,j) > Upit_NFw_avg);              % indexes of data points to keep for bin-averaged spectra
        Ulog_pwp_k(j) = mean((Upwp_est(Uk_inds,j)));                        % average over log of kept indexes (diagnostic: could be used for a crude epsilon estimate)
        Ulog_pwp_ks(j) = std((Upwp_est(Uk_inds,j)));                        % standard deviation of the log mean (diagnostic)
        Ulog_pwp_k_avg(j) = mean((Upwp_est_avg(Uk_inds_avg,j)));            % average over log of kept indexes for bin-averaged spectra
        Ulog_pwp_ks_avg(j) = std((Upwp_est_avg(Uk_inds_avg,j)));            % standard deviation of the log mean for bin-averaged spectra
   
        % unweight the average spectral level and standard deviations for plotting 
        % (creates a constant level spectra but in non-weighted space)
        Ulog_ppsd_est(:,j) = Ulog_pwp_k_avg(j) + w2uw;
        Ulog_ppsd_avg(:,j) = Uppsdw_avg(:,j) + w2uw;
        
        % Spectral fitting procedure is carried out
        % Also use the custom function to pick the frequency bins to be used in fitting and the quality criterion
        [k_inds_avgU(:,j),U_selind(j),bd_slpU(j),fwerU(j),ftStdDevIterU{j},pasU(j),Uk_inds_pit{j}, ...
            USAUW_WM(j),USAUW_WS(j),USAUW_WSTD(j),USAUWfit(1:length(f_nth_dec),j),USAUW_WSTDplt(:,j),k_inds_pltU(:,j)] = ...
            specfit2(f_nth_dec, Uppsdw_avg(:,j), Uppsdw_avg(:,j)+w2uw, Upit_NFw(:,j)+w2uw, pass_2, min_fit_pts,k_inds_pltU(:,j),uw2w,USAUW_WSTDplt(:,j));
    
        % if the number of selected points are less than the min number of
        % poits needed (set by the user) to get a decent fit estimate, then
        % don't bother with such spectra [set to nan by default - generally indicates that all the points in the spectra are suspect]
        Ulog_pwp_k_avg(j) = mean((Uppsdw_avg(Uk_inds_pit{j},j))); % average over log of kept indexes - in second pass
        Ulog_pwp_ks_avg(j) = std((Uppsdw_avg(Uk_inds_pit{j},j))); % standard deviation of the log mean - in second pass
        
        % compute the fit slopes and levels for non-weighted spectra (used points only)
        for q = 1:length(Uk_inds_pit{j})-1
            tmp1 = polyfit(log10(f_nth_dec(Uk_inds_pit{j}(q:q+1))),...
                UPPSD_AVG(Uk_inds_pit{j}(q:q+1),j),1);
            unConstSlpUusd(Uk_inds_pit{j}(q),j) = tmp1(:,1);           % successive fit slopes of bin-averaged spectra
            unConstLvlUusd(Uk_inds_pit{j}(q),j) = mean(UPPSD_AVG(Uk_inds_pit{j}(q:q+1),j),'omitnan');           % successive fit levels of bin-averaged spectra
        end

        % also compute the unconstrained fit slope and fit level for the
        % full frequency range of the spectra
        tmp2 = polyfit(log10(f_nth_dec(Uk_inds_avg)),UPPSD_AVG(Uk_inds_avg,j),1);
        unConstSlpUFlSp(j) = tmp2(:,1);
        unConstLvlUFlSp(j) = mean(UPPSD_AVG(Uk_inds_avg,j),'omitnan');
    
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
    unConstSlpUTraj(:,:,i) = unConstSlpU;
	unConstLvlUTraj(:,:,i) = unConstLvlU;
    unConstSlpUusdTraj(:,:,i) = unConstSlpUusd;
    unConstLvlUusdTraj(:,:,i) = unConstLvlUusd;
    unConstSlpUFlSpTraj(:,i) = unConstSlpUFlSp;
	unConstLvlUFlSpTraj(:,i) = unConstLvlUFlSp;
end

% Compute trajectory averaged spectra and epsilon for Wz
for pt = 1:1:(numInts/Trec)
    % calculate mean airspeed
    tmpW = mTrWz(pt,:);
    avmTrWz(pt) = mean(tmpW,2);
    % average the spectra for each altitude 
    tmpptU = log10(reshape(UppdfPlt(:,pt,:),[length(freq),numTraj]));
    avUppsd(:,pt) = 10.^(mean(tmpptU,2));
    % conduct the spectral fit
    % weight PSD by f^(5/3)
    avUppsdw(:,pt) =  avUppsd(:,pt).*(freq'.^(5/3));    % [m^2 s^(-8/3)]

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
    avUlog_ppsd_est(:,pt) = avUlog_pwp_k_avg(pt) + w2uw;
    avUlog_ppsd_avg(:,pt) = log10(avUppsdw_avg(:,pt)) + w2uw;

    % initialize arrays
    k_inds_avgU = nan(length(f_nth_dec),numInts/Trec);
    k_inds_pltU_AVG = nan(length(f_nth_dec),numInts/Trec);
    USAUW_WSTDplt_AVG = nan(length(f_nth_dec),numInts/Trec);
    unConstSlpU_AVG = nan(length(f_nth_dec)-1,numInts/Trec);
    unConstLvlU_AVG = nan(length(f_nth_dec)-1,numInts/Trec);
        
    % Spectral fitting for plane-averaged spectrum is carried out
    % Also use the custom function to pick the frequency bins to be used in fitting and the quality criterion
    [k_inds_avgU(:,pt),U_selind_AVG(j),bd_slpU_AVG(j),fwerU_AVG(j),ftStdDevIterU_AVG{j},pasU_AVG(j),avUk_inds_pit{pt}, ...
        USAUW_WM_AVG(j),USAUW_WS_AVG(j),USAUW_WSTD_AVG(j),USAUWfit_AVG(1:length(f_nth_dec),j),USAUW_WSTDplt_AVG(:,j),k_inds_pltU_AVG(:,j)] = ...
        specfit2(f_nth_dec, Uppsdw_avg(:,j), Uppsdw_avg(:,j)+w2uw, Upit_NFw(:,j)+w2uw, pass_2, min_fit_pts,k_inds_pltU_AVG(:,j),uw2w,USAUW_WSTDplt_AVG(:,j));
    
    % if the number of selected points are less than the min number of
    % poits needed (set by the user) to get a decent fit estimate, then
    % don't bother with such spectra [set to nan by default - generally indicates that all the points in the spectra are suspect]
    avUlog_pwp_k_avg(pt) = mean(log10(avUppsdw_avg(avUk_inds_pit{pt},pt))); % average over log of kept indexes - in second pass
    avUlog_pwp_ks_avg(pt) = std(log10(avUppsdw_avg(avUk_inds_pit{pt},pt))); % standard deviation of the log mean - in second pass

    % Compute Epsilons for each trajectory
    log_avUepsilon_k(pt) = 3/2*avUlog_pwp_k_avg(pt) - 3/2*log10(alpT) - log10(avmTrWz(pt));
    % compute error bars for epsilon
    log_avUepsilon_d(pt) = 3/2*avUlog_pwp_ks_avg(pt);
    log_avUepsilon_u(pt) = log_avUepsilon_k(pt) + log_avUepsilon_d(pt);
    log_avUepsilon_l(pt) = log_avUepsilon_k(pt) - log_avUepsilon_d(pt);
end