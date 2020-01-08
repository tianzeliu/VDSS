%% Compute SsPmp phase shift as a function of lower-crustal and uppermost-mantle Vp
% Based on Aki and Richards, 1980, pp. 149-150
% 1 is lower crust. 2 is uppermost mantle.
% ray_eff is the effective ray parameter that takes into account Moho dip!
%
% History:
% Created.
% 01/11/2018, Tianze Liu
%
% The code is chanaged to a function with the velocities and ray parameter
% as inputs and phase shift as output.
% 05/03/2018, Tianze Liu

function phs = PhaseShiftMoho(rayp_eff,vp1,vp2,vs1,vs2,rho1,rho2)
    ita_p_1 = sqrt(1/vp1^2-rayp_eff^2);
    ita_p_2 = sqrt(1/vp2^2-rayp_eff^2);
    ita_s_1 = sqrt(1/vs1^2-rayp_eff^2);
    ita_s_2 = sqrt(1/vs2^2-rayp_eff^2);

    % First-level terms
    a = rho2*(1-2*vs2^2*rayp_eff^2)-rho1*(1-2*vs1^2*rayp_eff^2);
    b = rho2*(1-2*vs2^2*rayp_eff^2)+2*rho1*vs1^2*rayp_eff^2;
    c = rho1*(1-2*vs1^2*rayp_eff^2)+2*rho2*vs2^2*rayp_eff^2;
    d = 2*(rho2*vs2^2-rho1*vs1^2);

    % Second-level terms
    ee = b*ita_p_1+c*ita_p_2;
    ff = b*ita_s_1+c*ita_s_2;
    gg = a-d*ita_p_1*ita_s_2;
    hh = a-d*ita_p_2*ita_s_1;
    dd = ee*ff+gg*hh*rayp_eff^2;

    % The relfection and transmission coefficients
    p_to_p = ((b*ita_p_1-c*ita_p_2)*ff-(a+d*ita_p_1*ita_s_2)*hh*rayp_eff^2)/dd;
    p_to_p = -p_to_p; % This is to take into account the discrepancy in sign convention of incident S wave at the free surface

    % Measure the phase shift
    phs = rad2deg(angle(p_to_p));
end