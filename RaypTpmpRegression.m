%% Find crustal thickness and average vp simultaneously by linear regression of SsPmp move-out curve
%
% History:
% Created.
% Tianze Liu, 08/30/2018
%
% The code is changed to a function.
% Tianze Liu, 09/26/2018
% 
% Estimation of error is added.
% Tianze Liu, 01/16/2019
%
% The random samples of vp_best and h_best are computed
% Tianze Liu, 01/24/2019
%
% The observations are weighted by their standard deviation
% Tianze Liu, 05/09/2019

function [h_best,vp_best,h_std,vp_std,H_rnd_out,Vp_rnd_out] = RaypTpmpRegression(Rayp_obs,Tpmp_obs,Tpmp_obs_std,n_rnd)
    % Build the inverse problem
    n_obs = length(Rayp_obs);
    Rayp_obs = reshape(Rayp_obs,[],1);
    Tpmp_obs = reshape(Tpmp_obs,[],1);
    Tpmp_obs_std = reshape(Tpmp_obs_std,[],1);
    
    X_obs = Rayp_obs.^2;
    Y_obs = Tpmp_obs.^2/4;
    Y_obs_std = sqrt(Tpmp_obs.^2/4.*Tpmp_obs_std.^2);
    
    G = [ones(n_obs,1),X_obs];
    W = diag(1./Y_obs_std);
    G_w = W*G;
    Y_obs_w = W*Y_obs;
    
    % Perform linear regresion using the normal equation
    M = (G_w'*G_w)\G_w'*Y_obs_w;
    m0 = M(1);
    m1 = M(2);

    % Find the best-fit vp and h
    h_best = sqrt(-m1);
    vp_best = sqrt(h_best^2/m0);
    
    % Estimate uncertainty
    M_cov = inv(G_w'*G_w)*G_w'*eye(n_obs,n_obs)*G_w*inv(G_w'*G_w);
    
    m0_var = M_cov(1,1);
    m1_var = M_cov(2,2);
    m0m1_cov = M_cov(1,2);
    
    h_var = -m1_var/(4*m1);
    vp_var = -m1_var/(4*m1*m0)-m1*m0_var/(4*m0^3)+m0m1_cov/(2*m0^2);
    
    h_std = sqrt(h_var);
    vp_std = sqrt(vp_var);
    
    % Compute random samples of m0 and m1
    M_rnd = mvnrnd([m0,m1],[m0_var,m0m1_cov;m0m1_cov,m1_var],n_rnd);
    M0_rnd = M_rnd(:,1);
    M1_rnd = M_rnd(:,2);
    
    % Remove m1 samples above 0
    M0_rnd_out = M0_rnd(M1_rnd<0);
    M1_rnd_out = M1_rnd(M1_rnd<0);
    
    % Convert m0 and m1 samples to h and vp samples
    H_rnd_out = sqrt(-M1_rnd_out);
    Vp_rnd_out = sqrt(H_rnd_out.^2./M0_rnd_out);
end

