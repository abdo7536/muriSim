%% Loop over W trajectories (includes averaging of TKE DNS)
for i = 1:1:numTraj
    Wtrj = Wz(:,i);
    Etrj = epsDNSz(:,i);
    for j = 1:1:(numInts/Trec)
        % initialize arrays
        k_inds_avgW = nan(length(f_nth_dec),numInts/Trec);
        k_inds_pltW = nan(length(f_nth_dec),numInts/Trec);
        WSAUW_WSTDplt = nan(length(f_nth_dec),numInts/Trec);
        %unConstSlpW = nan(length(f_nth_dec),numInts/Trec);
        %unConstLvlW = nan(length(f_nth_dec),numInts/Trec);
        unConstSlpWusd = nan(length(f_nth_dec)-1,numInts/Trec);
        unConstLvlWusd = nan(length(f_nth_dec)-1,numInts/Trec);
        
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
            unConstSlpW(q,j) = tmp1(:,1);
            unConstLvlW(q,j) = mean(WPPSD_AVG(q:q+1,j),'omitnan');
        end
        
        % variables to convert bin-averaged spectra to/from weighted-non-weighted forms
        uw2w = log10(f_nth_dec.^(5/3))';         % add this to weight the spectra
        w2uw = log10(f_nth_dec.^(-5/3))';        % add this to unweight the spectra
        
        % remove data in the interval (f_ind_low:f_ind_high) at or below the sensor noise floor
        Wpit_NFw = log10(pit_NF*freq(f_ind_low:f_ind_high).^(5/3));
        Wpit_NFw_avg = log10(pit_NF) + uw2w;
        Wpwp_est(:,j) = Wppsdw(f_ind_low:f_ind_high,j);
        Wpwp_est_avg(:,j) = Wppsdw_avg(:,j);
        Wk_inds = find(Wpwp_est(:,j) > Wpit_NFw');   % indexes of data points to keep
        Wk_inds_avg = find(Wpwp_est_avg(:,j) > Wpit_NFw_avg);   % indexes of data points to keep
        Wlog_pwp_k(j) = mean((Wpwp_est(Wk_inds,j)));      % average over log of kept indexes
        Wlog_pwp_ks(j) = std((Wpwp_est(Wk_inds,j)));      % standard deviation of the log mean
        Wlog_pwp_k_avg(j) = mean((Wpwp_est_avg(Wk_inds_avg,j))); % average over log of kept indexes
        Wlog_pwp_ks_avg(j) = std((Wpwp_est_avg(Wk_inds_avg,j))); % standard deviation of the log mean
        
        % unweight these estimates for plotting and checking
        % for the average spectra
        Wlog_ppsd_est(:,j) = Wlog_pwp_k_avg(j) + w2uw;
        Wlog_ppsd_avg(:,j) = Wppsdw_avg(:,j) + w2uw;
        
        % Spectral fitting procedure is carried out
        % Also use the custom function to pick the frequency bins to be used in fitting and the quality criterion
        [k_inds_avgW(:,j),W_selind(j),bd_slpW(j),fwerW(j),ftStdDevIterW{j},pasW(j),Wk_inds_pit{j}, ...
            WSAUW_WM(j),WSAUW_WS(j),WSAUW_WSTD(j),WSAUWfit(1:length(f_nth_dec),j),WSAUW_WSTDplt(:,j),k_inds_pltW(:,j)] = ...
            specfit2(f_nth_dec, Wppsdw_avg(:,j), Wppsdw_avg(:,j)+w2uw, Wpit_NFw(:,j)+w2uw, pass_2, min_fit_pts,k_inds_pltW(:,j),uw2w,WSAUW_WSTDplt(:,j));
        
        % if the number of selected points are less than the min number of
        % poits needed (set by the user) to get a decent fit estimate, then
        % don't bother with such spectra [set to nan by default - generally indicates that all the points in the spectra are suspect]
        Wlog_pwp_k_avg(j) = mean((Wppsdw_avg(Wk_inds_pit{j},j))); % average over log of kept indexes - in second pass
        Wlog_pwp_ks_avg(j) = std((Wppsdw_avg(Wk_inds_pit{j},j))); % standard deviation of the log mean - in second pass

        % compute the fit slopes and levels for non-weighted spectra (used points only)
        for q = 1:length(Wk_inds_pit{j})-1
            tmp1 = polyfit(log10(f_nth_dec(Wk_inds_pit{j}(q:q+1))),...
                WPPSD_AVG(Wk_inds_pit{j}(q:q+1),j),1);
            unConstSlpWusd(Wk_inds_pit{j}(q),j) = tmp1(:,1);           % successive fit slopes of bin-averaged spectra
            unConstLvlWusd(Wk_inds_pit{j}(q),j) = mean(WPPSD_AVG(Wk_inds_pit{j}(q:q+1),j),'omitnan');           % successive fit levels of bin-averaged spectra
        end

        % also compute the unconstrained fit slope and fit level for the
        % full frequency range of the spectra
        tmp2 = polyfit(log10(f_nth_dec(Wk_inds_avg)),WPPSD_AVG(Wk_inds_avg,j),1);
        unConstSlpWFlSp(j) = tmp2(:,1);
        unConstLvlWFlSp(j) = mean(WPPSD_AVG(Wk_inds_avg,j),'omitnan');
        
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
    unConstSlpWusdTraj(:,:,i) = unConstSlpWusd;
    unConstLvlWusdTraj(:,:,i) = unConstLvlWusd;
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
    avWlog_ppsd_est(:,pt) = avWlog_pwp_k_avg(pt) + w2uw;
    avWlog_ppsd_avg(:,pt) = log10(avWppsdw_avg(:,pt)) + w2uw;

    % initialize arrays
    k_inds_avgW = nan(length(f_nth_dec),numInts/Trec);
    k_inds_pltW_AVG = nan(length(f_nth_dec),numInts/Trec);
    WSAUW_WSTDplt_AVG = nan(length(f_nth_dec),numInts/Trec);
    unConstSlpW_AVG = nan(length(f_nth_dec),numInts/Trec);
    unConstLvlW_AVG = nan(length(f_nth_dec),numInts/Trec);
        
    % Spectral fitting for plane-averaged spectrum is carried out
    % Also use the custom function to pick the frequency bins to be used in fitting and the quality criterion
    [k_inds_avgW(:,pt),W_selind_AVG(j),bd_slpW_AVG(j),fwerW_AVG(j),ftStdDevIterW_AVG{j},pasW_AVG(j),avWk_inds_pit{pt}, ...
        WSAUW_WM_AVG(j),WSAUW_WS_AVG(j),WSAUW_WSTD_AVG(j),WSAUWfit_AVG(1:length(f_nth_dec),j),WSAUW_WSTDplt_AVG(:,j),k_inds_pltW_AVG(:,j)] = ...
        specfit2(f_nth_dec, Wppsdw_avg(:,j), Wppsdw_avg(:,j)+w2uw, Wpit_NFw(:,j)+w2uw, pass_2, min_fit_pts,k_inds_pltW_AVG(:,j),uw2w,WSAUW_WSTDplt_AVG(:,j));
    
    % if the number of selected points are less than the min number of
    % poits needed (set by the user) to get a decent fit estimate, then
    % don't bother with such spectra [set to nan by default - generally indicates that all the points in the spectra are suspect]
    avWlog_pwp_k_avg(pt) = mean(log10(avWppsdw_avg(avWk_inds_pit{pt},pt))); % average over log of kept indexes - in second pass
    avWlog_pwp_ks_avg(pt) = std(log10(avWppsdw_avg(avWk_inds_pit{pt},pt))); % standard deviation of the log mean - in second pass

    % Compute Epsilons for each trajectory
    log_avWepsilon_k(pt) = 3/2*avWlog_pwp_k_avg(pt) - 3/2*log10(alpL) - log10(avmTrWz(pt));
    log_avWepsilon_k2(pt) = 3/2*avWlog_pwp_k_avg(pt) - 3/2*log10(cf_rad) - log10(avmTrWz(pt));
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