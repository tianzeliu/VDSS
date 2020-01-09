%% Trace a 2D PmP ray across a layered 2D model with a specific ray parameter
% At this point, the code only handles a single-layer crust and the reflection at the base of the crust
%
% The ray starts from location x at top of the model (free surface), with
% incident angle ang (in deg). ang > 0 stands for vertical to right and ang < 0
% stands for vertical to left. The structure array Interface contain interface
% information.
%
% History:
% Modified from RayTraceS_2Dlyr.
% Tianze Liu, 03/26/2018
%
% The output of effective ray parameters at the Moho is added.
% Tianze Liu, 05/03/2018

function [X_pt,Z_pt,pcr,rayp_eff] = RayTracePmP_2Dlyr(x_inc,ang_inc,Interface)    
    % Define the down-going ray at the top of the model
    interface_srf = Interface(end);
    Z_srf = interface_srf.Z;
    z_srf = mean(Z_srf);
    w_ray_in = -cotd(ang_inc);
    b_ray_in = z_srf-w_ray_in*x_inc;
    K_inc = [1;w_ray_in]/norm([1,w_ray_in]);
    
    X_pt = zeros(3,1);
    Z_pt = zeros(3,1);
    X_pt(1) = x_inc;
    Z_pt(1) = z_srf;
 
    % Loop over layers
    nlyr = length(Interface);
    interface_inc = Interface(nlyr-1);
    interface_tran = Interface(nlyr-2);
    vp_inc = interface_inc.vp;
    vp_tran = interface_tran.vp;
    X_bdr = interface_inc.X;
    Z_bdr = interface_inc.Z;
    n_bdr = length(Z_bdr);

    % Find the reflection point
    for j = 1:n_bdr-1
        x1_bdr = X_bdr(j);
        z1_bdr = Z_bdr(j);

        x2_bdr = X_bdr(j+1);
        z2_bdr = Z_bdr(j+1);           
        ind = (w_ray_in*x1_bdr-z1_bdr+b_ray_in)*(w_ray_in*x2_bdr-z2_bdr+b_ray_in);   

        % Find the crossing point between the ray and the interface
        if ind > 0
            continue
        else
            x_pt = ((z2_bdr-b_ray_in)*(x2_bdr-x1_bdr)-(z2_bdr-z1_bdr)*x2_bdr)/(w_ray_in*(x2_bdr-x1_bdr)-(z2_bdr-z1_bdr));
            z_pt = w_ray_in*x_pt+b_ray_in;
            c = [x_pt-x_inc,z_pt-z_srf]*K_inc;
            % Determine if the crossing point is at the correct side of
            % the initial point
            if c > 0
                X_pt(2) = x_pt; 
                Z_pt(2) = z_pt;
                break
            else
                continue
            end
        end

    end

    S_inc = K_inc/vp_inc;

    % Find the normal direction
    N = [-(z2_bdr-z1_bdr);x2_bdr-x1_bdr]/norm([-(z2_bdr-z1_bdr);x2_bdr-x1_bdr]);
    if N'*K_inc < 0
        N = -N;
    end

    % Calculate the reflection direction using Snell's law and determine if
    % the reflection is post-critical
    S_in_perp = (N'*S_inc)*N;
    S_in_para = S_inc-S_in_perp;
    S_ref_para = S_in_para;
    S_ref_perp = -S_in_perp;
    
    % The effective ray parameter at the Moho
    rayp_eff = norm(S_in_para);
    
    s_tran_perp = sqrt(1/vp_tran^2-norm(S_in_para)^2);
    
    if isreal(s_tran_perp)
        pcr = 0;
    else
        pcr = 1;
    end

    % Define the outgoing ray
    S_ref = S_ref_para+S_ref_perp;
    K_ref = S_ref/norm(S_ref);
    w_ray_ref = K_ref(2)/K_ref(1);
    b_ray_ref = (z_pt-K_ref(2)/K_ref(1)*x_pt);

    % Find the surfacing point
    x_end = (z_srf-b_ray_ref)/w_ray_ref;
    X_pt(3) = x_end;
    Z_pt(3) = z_srf;
end