%% Trace a 2D S ray across a layered 2D model with a specific incident angle at the bottom of the model
% The ray starts from location x at bottom of the model (z = 0), with
% incident angle ang (in deg). ang > 0 stands for vertical to right and ang < 0
% stands for vertical to left. The structure array Interface contain interface
% information.
%
%
% History:
% Created.
% Tianze Liu, 03/26/2018
%
% Treatment of post-critical reflection is added.
% Tianze Liu, 03/27/2018

function [X_pt,Z_pt] = RayTraceS_2Dlyr(x_inc,ang_inc,Interface)    
    % Define the incident ray at the bottom of the model
    w_ray_in = cotd(ang_inc);
    b_ray_in = -w_ray_in*x_inc;
    x0 = x_inc;
    z0 = 0;
    K_in = [1;w_ray_in]/norm([1,w_ray_in]);
    
    % Initialize the loop
    nlyr = length(Interface);
    X_pt = zeros(2*nlyr,1);
    Z_pt = zeros(2*nlyr,1);
    X_pt(1) = x0;
    Z_pt(1) = z0;
 
    
    % Loop over layers
    i = 2;
    n_pts = 1;
    while i <= nlyr
        interface_in = Interface(i-1);
        interface_out = Interface(i);
        vs_in = interface_in.vs;
        vs_out = interface_out.vs;
        X_bdr = interface_out.X;
        Z_bdr = interface_out.Z;
        n_bdr = length(Z_bdr);
        
        % Find the conversion point
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
                c = [x_pt-x0,z_pt-z0]*K_in;
                % Determine if the crossing point is at the correct side of
                % the initial point
                if c > 0
                    X_pt(n_pts+1) = x_pt; 
                    Z_pt(n_pts+1) = z_pt;
                    n_pts = n_pts+1;
                    break
                else
                    continue
                end
            end
    
        end
        
        if vs_out == 0
            X_pt = X_pt(1:n_pts);
            Z_pt = Z_pt(1:n_pts);
            return
        else

            S_in = K_in/vs_in;

            % Find the normal direction
            N = [-(z2_bdr-z1_bdr);x2_bdr-x1_bdr]/norm([-(z2_bdr-z1_bdr);x2_bdr-x1_bdr]);
            if N'*K_in < 0
                N = -N;
            end

            % Calculate the outgoing slowness using Snell's law, with
            % special
            % treatment of post-critical reflection
            S_in_perp = (N'*S_in)*N;
            S_in_para = S_in-S_in_perp;
            S_out_para = S_in_para;
            s_out_perp = sqrt(1/vs_out^2-norm(S_out_para)^2);
            if isreal(s_out_perp)
                S_out_perp = sqrt(1/vs_out^2-norm(S_out_para)^2)*N;
                i = i+1;
            else % In this case, post-critical reflection happens at the interface and the ray is trapped in the lower medium
                S_out_perp = -S_in_perp;
            end
            % Define the outgoing ray
            S_out = S_out_para+S_out_perp;
            K_out = S_out/norm(S_out);
            w_ray_out = K_out(2)/K_out(1);
            b_ray_out = (z_pt-K_out(2)/K_out(1)*x_pt);

            w_ray_in = w_ray_out;
            b_ray_in = b_ray_out;
            x0 = x_pt;
            z0 = z_pt;
            K_in = K_out;
        end
    end
end