% function to fit f^5/3 slope line to spectra
%% variables needed
    % f_nth_dec -> the vector of sampled frequencies [Hz]
    % ppsdw_avg -> the PSD weighted by f^5/3 [m^2 s^(-8/3)]
    % ppsd_avg -> the unweighted PSD [m^2 s^-1]
    % pit_NFw -> the signal noise floor
    % pass_2 -> a switch turns on an additional level scan on the spectra
    % the second pass checks the PSD point-wise to include/omit points that
    % contribute towards reducing the fit std.
    % min_fit_pts -> the minimum number of points to be used in the third
    % pass of fitting (only used if pass_2 is set to 1)
    % also store only the frequency bins for qualified points to plot statistics

function [k_inds_avg,Q,bd_st_slp,fw_inds,iterStdDev,ps,tmpvar2,FSAUW_WM2,FSAUW_WS2,FSAUW_WSTD2,FSAUWfit2,FSAUW_WSTDplt2,k_inds_plt2] = ...
    specfit2(f_nth_dec, ppsdw_avg, ppsd_avg, pit_NFw, pass_2, min_fit_pts,k_inds_plt2,fbcw,FSAUW_WSTDplt2)
    %% assign a nan array for chosen indices
    k_inds_avg = nan(size(f_nth_dec));
    %% first pass: omit PSD points at/below the noise floor
    ps = 1;                     % a variable that keeps track of number of passes
    % indexes of data points to keep 
    % NOTE: care must be taken to ensure the use of unweighted PSD with unweighted NF to qualify the points in PSD above NF 
    k_inds_avg_tmp = find(ppsdw_avg > pit_NFw);
    % check if the indices above the noise floor are continuous. If not,
    % then only retain the frequency bins that are above the noise floor and continuous
    ind_chk = find(diff(k_inds_avg_tmp) ~= 1,1,'first');
    if ~isempty(ind_chk)
        k_inds_avg_tmp = k_inds_avg_tmp(1:ind_chk);
    end
    % average over log of kept indexes
    log_pwp_k_avg = mean(ppsdw_avg(k_inds_avg_tmp));
    % standard deviation of the log mean
    log_pwp_ks_avg = std(ppsdw_avg(k_inds_avg_tmp));
    % the spectral sample is qualified by default
    Q = 1;                      % spectral sample is qualified
    bd_st_slp = 0;              % flag if the best fit slope is unreasonable (ge 20% error)
    fw_inds = 0;                % flag if there are fewer than min_fit_pts available for fit
    
    %% second pass: include the PSD points that contribute towards reducing the standard deviation
    % compute the fit standard deviation for all successive frequency bins
    % by looping over the used frequency bins
    if length(k_inds_avg_tmp) >= min_fit_pts
        for i = 1:1:length(k_inds_avg_tmp)-1
            % constrained fit using f^5/3 weighted PSD
            tmp1(i,1) = -5/3;
            tmp1(i,2) = mean(ppsdw_avg(k_inds_avg_tmp(i:i+1)));        % fit level
            % also conduct a linear first order fit to the two points (unweighted) to
            % extimate the actual slope and fit level
            tmp2(i,1:2) = polyfit(log10(f_nth_dec(k_inds_avg_tmp(i:i+1))),ppsd_avg(k_inds_avg_tmp(i:i+1)),1);
        end
        % find the fit slope for all consecutive frequency bins
        %[bfFrBnIndMV,bfFrBnI] = min(abs(abs(tmp2(:,1)) - abs(tmp1(:,1))));
        [bfFrBnIndMV,bfFrBnI] = min(abs(tmp2(:,1) - tmp1(:,1)));
        % if the best fit slope is "unreasonably" poor, then reject the
        % spectral sample
        % Error Criterion: (greater than 20% error in fit slope compared to the expected ISR slope)
        % [(-5/3-0.33)<(-5/3)<(-5/3+0.33)] --> upto 20% error is tolerated
        if bfFrBnIndMV > (0.3*(5/3))            % tolerates 30% error in fitting the ISR
            Q = 0;                  % spectral sample is rejected
            k_inds_avg_tmp = [bfFrBnI bfFrBnI+1];   % also store the "best" fitting freuency bins for diagnosis
            bd_st_slp = 1;
        end
    else
        k_inds_avg_tmp = transpose(k_inds_avg_tmp);
        Q = 0;                      % spectral sample is rejected
        fw_inds = 1;
    end
    % only apply the second pass if the width of the usable frequency range is appreciable 
	% for ex: for 1/3rd decade bins, min_fit_pts = 3 will require atleast 0.9 decade for a reasonable fit
    if pass_2 == 1 && length(k_inds_avg_tmp) >= min_fit_pts && Q == 1
        ps = ps + 1;
        % declare arrays to vary tolerances for fit standard deviation if
        % necessary
        tol0 = [0.05 0.22];         % max of 0.22
        tol1 = [0.05 0.06 4];       % max of 0.24
        tol2 = [0.06 0.075 3];      % max of 0.225
        tol3 = [0.075 0.15 2];      % max of 0.3
        tol4 = [0.15 0.25 1.3];     % max of 0.325
        tol5 = [0.25 0.40 1.05];    % max of 0.42
        tol6 = [0.40 1.0];       
        
        clear stand inds_to_use         % clear recycled variables
        % loop over all the retained frequency bins such and evaluate their relative contributions to fit std. dev.
        % first, loop over the lower frequency end of the spectrum
        p = bfFrBnI;
        cntr = 1;
        ind_to_use = k_inds_avg_tmp(p:bfFrBnI+1);
        while p > k_inds_avg_tmp(1)-1
            % use the first three points to begin with and then iterate over the remaining points
            % NOTE: conduct all fit standard deviation checks on f^5/3 weighted PSD
            stand(cntr) = std(ppsdw_avg(ind_to_use));               % first compute the standard deviation using the "best" two frequency bin data points
            iterStdDevL(cntr) = stand(cntr);
            if p < bfFrBnI && p >= k_inds_avg_tmp(1)+1                  % for all frequency bin points that are below the lower frequency bin of the best fit
                if stand(cntr-1) < tol0(1)                          % if the fit standard deviation is 0<stand<0.25
                    if stand(cntr) > tol0(2)                        % and contributing to greater than 20% increase in standard deviation upon inclusion
                        ind_to_use = [k_inds_avg_tmp(p-1); k_inds_avg_tmp(ind_to_use(2:end))];  % exclude the data point in question
                        iterStdDevL(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = [k_inds_avg_tmp(p-1); k_inds_avg_tmp(ind_to_use)];         % include the point in question    
                    end    
                end
                if stand(cntr-1) > tol1(1) && stand(cntr-1) < tol1(2)   % if the fit standard deviation is 0<stand<0.25
                    if stand(cntr) > tol1(3)*stand(cntr-1)              % and contributing to greater than 20% increase in standard deviation upon inclusion
                        ind_to_use = [k_inds_avg_tmp(p-1); k_inds_avg_tmp(ind_to_use(2:end))];  % exclude the data point in question
                        iterStdDevL(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = [k_inds_avg_tmp(p-1); k_inds_avg_tmp(ind_to_use)];         % include the point in question    
                    end    
                end
                if stand(cntr-1) > tol2(1) && stand(cntr-1) < tol2(2)          % if the fit standard deviation is 0.25<stand<0.5
                    if stand(cntr) > tol2(3)*stand(cntr-1)             % and contributing to greater than 15% increase in standard deviation upon inclusion
                        ind_to_use = [k_inds_avg_tmp(p-1); k_inds_avg_tmp(ind_to_use(2:end))];  % exclude the data point in question
                        iterStdDevL(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = [k_inds_avg_tmp(p-1); k_inds_avg_tmp(ind_to_use)];         % include the point in question    
                    end    
                end
                if stand(cntr-1) > tol3(1) && stand(cntr-1) < tol3(2)          % if the fit standard deviation is 0.5<stand<0.75
                    if stand(cntr) > tol3(3)*stand(cntr-1)             % and contributing to greater than 12% increase in standard deviation upon inclusion
                        ind_to_use = [k_inds_avg_tmp(p-1); k_inds_avg_tmp(ind_to_use(2:end))];  % exclude the data point in question
                        iterStdDevL(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = [k_inds_avg_tmp(p-1); k_inds_avg_tmp(ind_to_use)];         % include the point in question    
                    end    
                end
                if stand(cntr-1) > tol4(1) && stand(cntr-1) < tol4(2)            % if the fit standard deviation is 0.5<stand<0.75
                    if stand(cntr) > tol4(3)*stand(cntr-1)             % and contributing to greater than 10% increase in standard deviation upon inclusion
                        ind_to_use = [k_inds_avg_tmp(p-1); k_inds_avg_tmp(ind_to_use(2:end))];  % exclude the data point in question
                        iterStdDevL(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = [k_inds_avg_tmp(p-1); k_inds_avg_tmp(ind_to_use)];         % include the point in question    
                    end    
                end
                if stand(cntr-1) > tol5(1) && stand(cntr-1) < tol5(2)             % if the fit standard deviation is 1<stand<1.5
                    if stand(cntr) > tol5(3)*stand(cntr-1)             % and contributing to greater than 5% increase in standard deviation upon inclusion
                        ind_to_use = [k_inds_avg_tmp(p-1); k_inds_avg_tmp(ind_to_use(2:end))];  % exclude the data point in question
                        iterStdDevL(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = [k_inds_avg_tmp(p-1); k_inds_avg_tmp(ind_to_use)];         % include the point in question    
                    end    
                end
                if stand(cntr-1) >= tol6(1)                               % if the fit standard deviation is 1.5<stand
                    if stand(cntr) > tol6(2)*stand(cntr-1)             % and contributing to greater than 2% increase in standard deviation upon inclusion
                        ind_to_use = [k_inds_avg_tmp(p-1); k_inds_avg_tmp(ind_to_use(2:end))];  % exclude the data point in question
                        iterStdDevL(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = [k_inds_avg_tmp(p-1); k_inds_avg_tmp(ind_to_use)];         % include the point in question
                    end    
                end
            elseif p == bfFrBnI && p == k_inds_avg_tmp(1)               % if the best fit happens to also be the first and second frequency bins
                ind_to_use = k_inds_avg_tmp(p:bfFrBnI+1);               % retain these points and continue iterating
                ptsL = transpose(ind_to_use);
                ptsL(find(ptsL==bfFrBnI)) = [];
                ptsL(find(ptsL==bfFrBnI+1)) = [];
                iterStdDevL(cntr) = stand(cntr);
            elseif p == bfFrBnI && p ~= k_inds_avg_tmp(1)               % for the two best fit frequency bin data points
                ind_to_use = k_inds_avg_tmp(p-1:bfFrBnI+1);             % retain these points and continue iterating
                iterStdDevL(cntr) = stand(cntr);
            elseif p == k_inds_avg_tmp(1)                               % for the lowest frequency bin data point
                if stand(cntr-1) < tol0(1)                          % if the fit standard deviation is 0<stand<0.25
                    if stand(cntr) > tol0(2)                        % and contributing to greater than 20% increase in standard deviation upon inclusion
                        ind_to_use = k_inds_avg_tmp(ind_to_use(2:end)); % exclude the data point in question
                        iterStdDevL(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = k_inds_avg_tmp(ind_to_use);        % include the point in question    
                    end    
                end
                if stand(cntr-1) > tol1(1) && stand(cntr-1) < tol1(2)            % if the fit standard deviation is 0<stand<0.25
                    if stand(cntr) > tol1(3)*stand(cntr-1)              % and contributing to greater than 20% increase in standard deviation upon inclusion
                        ind_to_use = k_inds_avg_tmp(ind_to_use(2:end)); % exclude the data point in question
                        iterStdDevL(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = k_inds_avg_tmp(ind_to_use);        % include the point in question    
                    end    
                end
                if stand(cntr-1) > tol2(1) && stand(cntr-1) < tol2(2)          % if the fit standard deviation is 0.25<stand<0.5
                    if stand(cntr) > tol2(3)*stand(cntr-1)             % and contributing to greater than 15% increase in standard deviation upon inclusion
                        ind_to_use = k_inds_avg_tmp(ind_to_use(2:end)); % exclude the data point in question
                        iterStdDevL(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = k_inds_avg_tmp(ind_to_use);        % include the point in question    
                    end    
                end
                if stand(cntr-1) > tol3(1) && stand(cntr-1) < tol3(2)          % if the fit standard deviation is 0.5<stand<0.75
                    if stand(cntr) > tol3(3)*stand(cntr-1)             % and contributing to greater than 12% increase in standard deviation upon inclusion
                        ind_to_use = k_inds_avg_tmp(ind_to_use(2:end)); % exclude the data point in question
                        iterStdDevL(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = k_inds_avg_tmp(ind_to_use);        % include the point in question
                    end
                end
                if stand(cntr-1) > tol4(1) && stand(cntr-1) < tol4(2)            % if the fit standard deviation is 0.5<stand<0.75
                    if stand(cntr) > tol4(3)*stand(cntr-1)             % and contributing to greater than 10% increase in standard deviation upon inclusion
                        ind_to_use = k_inds_avg_tmp(ind_to_use(2:end)); % exclude the data point in question
                        iterStdDevL(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = k_inds_avg_tmp(ind_to_use);        % include the point in question    
                    end    
                end
                if stand(cntr-1) > tol5(1) && stand(cntr-1) < tol5(2)             % if the fit standard deviation is 1<stand<1.5
                    if stand(cntr) > tol5(3)*stand(cntr-1)             % and contributing to greater than 5% increase in standard deviation upon inclusion
                        ind_to_use = k_inds_avg_tmp(ind_to_use(2:end)); % exclude the data point in question
                        iterStdDevL(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = k_inds_avg_tmp(ind_to_use);        % include the point in question    
                    end    
                end
                if stand(cntr-1) >= tol6(1)                               % if the fit standard deviation is 1.5<stand
                    if stand(cntr) > tol6(2)*stand(cntr-1)             % and contributing to greater than 2% increase in standard deviation upon inclusion
                        ind_to_use = k_inds_avg_tmp(ind_to_use(2:end)); % exclude the data point in question
                        iterStdDevL(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = k_inds_avg_tmp(ind_to_use);        % include the point in question
                    end    
                end
                ptsL = transpose(ind_to_use);                       % store the retained points
                ptsL(find(ptsL==bfFrBnI)) = [];
                ptsL(find(ptsL==bfFrBnI+1)) = [];
                % add more checks to mitigate false qualified spectra
                % 1) check if one or more data points are separated more than one intermediate data point
                if ~isempty(ptsL)
                    if ~isempty(find((bfFrBnI-ptsL) > 1))
                        ptsL(find((bfFrBnI-ptsL) > 1)) = [];
                    end
                end
            end
            p = p - 1;
            cntr = cntr + 1;
        end
        clear p cntr ind_to_use stand       % clear recycled variables
        % next, loop over the higher frequency end of the spectrum
        p = bfFrBnI+1;
        cntr = 1;
        ind_to_use = k_inds_avg_tmp(bfFrBnI:p);
        while p < k_inds_avg_tmp(end)+1
            % use the first three points to begin with and then iterate over the remaining points
            % NOTE: conduct all fit standard deviation checks on f^5/3 weighted PSD
            stand(cntr) = std(ppsdw_avg(ind_to_use));               % first compute the standard deviation using the "best" two frequency bin data points
            iterStdDevH(cntr) = stand(cntr);
            if p > bfFrBnI+1 && p <= k_inds_avg_tmp(end)-1              % for all frequency bin points that are below the lower frequency bin of the best fit
                if stand(cntr-1) < tol0(1)                          % if the fit standard deviation is 0<stand<0.25
                    if stand(cntr) > tol0(2)                        % and contributing to greater than 20% increase in standard deviation upon inclusion
                        ind_to_use = [k_inds_avg_tmp(ind_to_use(1:end-1)); k_inds_avg_tmp(p+1)];  % exclude the data point in question
                        iterStdDevH(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = [k_inds_avg_tmp(ind_to_use); k_inds_avg_tmp(p+1)];         % include the point in question
                    end    
                end
                if stand(cntr-1) > tol1(1) && stand(cntr-1) < tol1(2)            % if the fit standard deviation is 0<stand<0.25
                    if stand(cntr) > tol1(3)*stand(cntr-1)              % and contributing to greater than 20% increase in standard deviation upon inclusion
                        ind_to_use = [k_inds_avg_tmp(ind_to_use(1:end-1)); k_inds_avg_tmp(p+1)];  % exclude the data point in question
                        iterStdDevH(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = [k_inds_avg_tmp(ind_to_use); k_inds_avg_tmp(p+1)];         % include the point in question
                    end    
                end
                if stand(cntr-1) > tol2(1) && stand(cntr-1) < tol2(2)          % if the fit standard deviation is 0.25<stand<0.5
                    if stand(cntr) > tol2(3)*stand(cntr-1)             % and contributing to greater than 15% increase in standard deviation upon inclusion
                        ind_to_use = [k_inds_avg_tmp(ind_to_use(1:end-1)); k_inds_avg_tmp(p+1)];  % exclude the data point in question
                        iterStdDevH(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = [k_inds_avg_tmp(ind_to_use); k_inds_avg_tmp(p+1)];         % include the point in question
                    end    
                end
                if stand(cntr-1) > tol3(1) && stand(cntr-1) < tol3(2)          % if the fit standard deviation is 0.5<stand<0.75
                    if stand(cntr) > tol3(3)*stand(cntr-1)             % and contributing to greater than 12% increase in standard deviation upon inclusion
                        ind_to_use = [k_inds_avg_tmp(ind_to_use(1:end-1)); k_inds_avg_tmp(p+1)];  % exclude the data point in question
                        iterStdDevH(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = [k_inds_avg_tmp(ind_to_use); k_inds_avg_tmp(p+1)];         % include the point in question
                    end    
                end
                if stand(cntr-1) > tol4(1) && stand(cntr-1) < tol4(2)            % if the fit standard deviation is 0.5<stand<0.75
                    if stand(cntr) > tol4(3)*stand(cntr-1)             % and contributing to greater than 10% increase in standard deviation upon inclusion
                        ind_to_use = [k_inds_avg_tmp(ind_to_use(1:end-1)); k_inds_avg_tmp(p+1)];  % exclude the data point in question
                        iterStdDevH(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = [k_inds_avg_tmp(ind_to_use); k_inds_avg_tmp(p+1)];         % include the point in question
                    end    
                end
                if stand(cntr-1) > tol5(1) && stand(cntr-1) < tol5(2)             % if the fit standard deviation is 1<stand<1.5
                    if stand(cntr) > tol5(3)*stand(cntr-1)             % and contributing to greater than 5% increase in standard deviation upon inclusion
                        ind_to_use = [k_inds_avg_tmp(ind_to_use(1:end-1)); k_inds_avg_tmp(p+1)];  % exclude the data point in question
                        iterStdDevH(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = [k_inds_avg_tmp(ind_to_use); k_inds_avg_tmp(p+1)];         % include the point in question
                    end    
                end
                if stand(cntr-1) >= tol6(1)                               % if the fit standard deviation is 1.5<stand
                    if stand(cntr) > tol6(2)*stand(cntr-1)              % and contributing to greater than 2% increase in standard deviation upon inclusion
                        ind_to_use = [k_inds_avg_tmp(ind_to_use(1:end-1)); k_inds_avg_tmp(p+1)];  % exclude the data point in question
                        iterStdDevH(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = [k_inds_avg_tmp(ind_to_use); k_inds_avg_tmp(p+1)];         % include the point in question
                    end    
                end
            elseif p == bfFrBnI+1 && p == k_inds_avg_tmp(end)             % for the two best fit frequency bin data points
                ind_to_use = k_inds_avg_tmp(bfFrBnI:p);                 % retain these points and continue iterating
                ptsH = transpose(ind_to_use);
                ptsH(find(ptsH==bfFrBnI)) = [];
                ptsH(find(ptsH==bfFrBnI+1)) = [];
                iterStdDevH(cntr) = stand(cntr);
            elseif p == bfFrBnI+1 && p ~= k_inds_avg_tmp(end)             % for the two best fit frequency bin data points
                ind_to_use = k_inds_avg_tmp(bfFrBnI:p+1);               % retain these points and continue iterating
                iterStdDevH(cntr) = stand(cntr);
            elseif p == k_inds_avg_tmp(end)                             % for the lowest frequency bin data point
                if stand(cntr-1) < tol0(1)                          % if the fit standard deviation is 0<stand<0.25
                    if stand(cntr) > tol0(2)                        % and contributing to greater than 20% increase in standard deviation upon inclusion
                        ind_to_use = k_inds_avg_tmp(ind_to_use(1:end-1)); % exclude the data point in question
                        iterStdDevH(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = k_inds_avg_tmp(ind_to_use);        % include the point in question    
                    end    
                end
                if stand(cntr-1) > tol1(1) && stand(cntr-1) < tol1(2)            % if the fit standard deviation is 0<stand<0.25
                    if stand(cntr) > tol1(3)*stand(cntr-1)              % and contributing to greater than 20% increase in standard deviation upon inclusion
                        ind_to_use = k_inds_avg_tmp(ind_to_use(1:end-1)); % exclude the data point in question
                        iterStdDevH(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = k_inds_avg_tmp(ind_to_use);        % include the point in question    
                    end    
                end
                if stand(cntr-1) > tol2(1) && stand(cntr-1) < tol2(2)          % if the fit standard deviation is 0.25<stand<0.5
                    if stand(cntr) > tol2(3)*stand(cntr-1)             % and contributing to greater than 15% increase in standard deviation upon inclusion
                        ind_to_use = k_inds_avg_tmp(ind_to_use(1:end-1)); % exclude the data point in question
                        iterStdDevH(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = k_inds_avg_tmp(ind_to_use);        % include the point in question    
                    end    
                end
                if stand(cntr-1) > tol3(1) && stand(cntr-1) < tol3(2)          % if the fit standard deviation is 0.5<stand<0.75
                    if stand(cntr) > tol3(3)*stand(cntr-1)             % and contributing to greater than 12% increase in standard deviation upon inclusion
                        ind_to_use = k_inds_avg_tmp(ind_to_use(1:end-1)); % exclude the data point in question
                        iterStdDevH(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = k_inds_avg_tmp(ind_to_use);        % include the point in question    
                    end
                end
                if stand(cntr-1) > tol4(1) && stand(cntr-1) < tol4(2)            % if the fit standard deviation is 0.5<stand<0.75
                    if stand(cntr) > tol4(3)*stand(cntr-1)             % and contributing to greater than 10% increase in standard deviation upon inclusion
                        ind_to_use = k_inds_avg_tmp(ind_to_use(1:end-1)); % exclude the data point in question
                        iterStdDevH(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = k_inds_avg_tmp(ind_to_use);        % include the point in question    
                    end    
                end
                if stand(cntr-1) > tol5(1) && stand(cntr-1) < tol5(2)             % if the fit standard deviation is 1<stand<1.5
                    if stand(cntr) > tol5(3)*stand(cntr-1)             % and contributing to greater than 5% increase in standard deviation upon inclusion
                        ind_to_use = k_inds_avg_tmp(ind_to_use(1:end-1)); % exclude the data point in question
                        iterStdDevH(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = k_inds_avg_tmp(ind_to_use);        % include the point in question    
                    end    
                end
                if stand(cntr-1) >= tol6(1)                               % if the fit standard deviation is 1.5<stand
                    if stand(cntr) > tol6(2)*stand(cntr-1)             % and contributing to greater than 2% increase in standard deviation upon inclusion
                        ind_to_use = k_inds_avg_tmp(ind_to_use(1:end-1)); % exclude the data point in question
                        iterStdDevH(cntr) = stand(cntr);
                        stand(cntr) = stand(cntr-1);
                    else
                        ind_to_use = k_inds_avg_tmp(ind_to_use);        % include the point in question
                    end    
                end
                ptsH = transpose(ind_to_use);                       % store the retained points
                ptsH(find(ptsH==bfFrBnI)) = [];
                ptsH(find(ptsH==bfFrBnI+1)) = [];
                % add more checks to mitigate false qualified spectra
                % 1) check if one or more data points are separated more than one intermediate data point
                if ~isempty(ptsH)
                    if ~isempty(find((ptsH-(bfFrBnI+1)) > 1))
                        ptsH(find(((bfFrBnI+1)-ptsH) > 1)) = [];
                    end
                end
            end
            p = p + 1;
            cntr = cntr + 1;
        end
        % combine all the retained points into an array
        k_inds_avg_tmp = [ptsL bfFrBnI bfFrBnI+1 ptsH];
        % arrange the selected frequency bins in the order of their
        % appearence in the frequency scale
        k_inds_avg(k_inds_avg_tmp) = k_inds_avg_tmp;
        % discard the first data point because this is repeated from iterStdDevL(1)
        iterStdDevH = iterStdDevH(2:end);
        % combine the standard deviation arrays
        iterStdDev = [iterStdDevL iterStdDevH];
        % reject points that have fewer fit points than min_fit_pts
        if length(k_inds_avg_tmp) < min_fit_pts
            Q = 0;                      % spectral sample is qualified
            fw_inds = 1;
        end
    else
        % write out empty output arrays for the constrained and
        % unconstrained fit slopes and fit levels
        iterStdDev = [];
    end
    % also nan the rejected/non-qualified samples
    tmpvar2 = k_inds_avg;
    tmpvar2(isnan(k_inds_avg)) = [];
    % Fit a zero slope line to this weighted data
    FSAUW_WM2 = mean(ppsdw_avg(tmpvar2));
	FSAUW_WS2 = 0;
    FSAUW_WSTD2 = std(ppsdw_avg(tmpvar2));
    FSAUWfit2 = FSAUW_WM2.*ones(length(f_nth_dec),1) - fbcw;
    % also store only the frequency bins for qualified points to plot statistics
    if Q == 1
        FSAUW_WSTDplt2 = FSAUW_WSTD2;
        k_inds_plt2 = k_inds_avg;
    end
end