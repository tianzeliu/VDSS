%% Get phase shift and travel time for a SsPmp record
%
% History:
% Created.
% Tianze Liu, 01/07/2020

function [phs_best,tpmp_best,msft_best,phs_std,tpmp_std,crv,Msft,P_syn_best,T_syn_best] = GetPhaseShiftAndTravelTime(P_obs,S_obs,Src,b_obs,dt,Phs_rg,msft_bd)
    n_phs = length(Phs_rg);
    dphs = Phs_rg(2)-Phs_rg(1);
    npts_src = length(Src);
    npts_obs = length(S_obs);
    I_b_msft = zeros(n_phs,1);
    I_e_msft = zeros(n_phs,1);
    Msft = zeros(n_phs,1);
    PP_syn = zeros(npts_src,n_phs); % Vector to store synthetic P records
    Tpmp = zeros(n_phs,1); % Vector to store measured Tpmp
            
    for i = 1:n_phs
        % Apply phase shift to the source wavelet
        phs = Phs_rg(i);

        Src_hb = imag(hilbert(Src));
        P_syn = Src*cosd(phs)+Src_hb*sind(phs);
        P_syn = P_syn/max(P_syn);        
        PP_syn(:,i) = P_syn;

        % The time lag with maximum cross-correlation value between the modeled SsPmp and P component record. 
        Xc_p = xcorr(P_obs,P_syn);
        [~,i_p] = max(Xc_p);

        % The time lag with maximum cross-correlation value between the
        % source wavelet and the S component record
        S_syn = Src;
        Xc_s = xcorr(S_obs,S_syn);
        [~,i_s] = max(Xc_s);

        % Calculate the misfit
        i_b_msft = i_p-npts_obs+1;
        i_e_msft = i_b_msft+npts_src-1;
        P_obs_msft = P_obs(i_b_msft:i_e_msft);
        Msft(i) = sqrt(sum((P_obs_msft-P_syn).^2)/sum(P_obs_msft.^2)); % We calculate the normalized square error misfit
        I_b_msft(i) = i_b_msft;
        I_e_msft(i) = i_e_msft;

        % Calculate the travel time
        tpmp = (i_p-i_s)*dt;
        Tpmp(i) = tpmp;
    end

    % Find the phase shift with the smallest misfit 
    [msft_best,i_best] = min(Msft);
    phs_best = Phs_rg(i_best);
    tpmp_best = Tpmp(i_best);
    i_b_msft_best = I_b_msft(i_best);
    P_syn_best = PP_syn(:,i_best);

    % Compute the uncertainty for H and phi
    crv = (Msft(i_best+1)+Msft(i_best-1)-2*msft_best)/dphs^2;
    phs_var = 2*msft_bd/crv;
    phs_std = sqrt(phs_var);
    phs_up = phs_best+phs_std;
    phs_lo = phs_best-phs_std;
    tpmp_up = interp1(Phs_rg,Tpmp,phs_up);
    tpmp_lo = interp1(Phs_rg,Tpmp,phs_lo);
    tpmp_std = (abs(tpmp_up-tpmp_best)+abs(tpmp_best-tpmp_lo))/2;
    
    % Find the correct alignment between the synthetic and observed SsPmp
    b_syn_best = b_obs+(i_b_msft_best-1)*dt;
    e_syn_best = b_syn_best+(npts_src-1)*dt;
    T_syn_best = (b_syn_best:dt:e_syn_best)';
 end