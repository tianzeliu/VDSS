%% Apply 2D Ray-tracing of PmP wave with varying ray parameter
%
% History:
% Modified from AppRayTraceS_2Dlyr.m.
% Tianze Liu, 04/23/2018
% 
% Prediction of SsPmp phases is added.
% Tianze Liu, 05/03/2018
%
% Two set of Pmp rays are traced. One with the apparent S ray parameters
% (evenly distributed rays taking off from the free surface ) and the other
% with the theoretical S ray parameters (unevenly distributed rays taking off from the free surface)
%
% Pre-critical and post-critical reflection rays are ouput into different
% files.
% Tianze Liu, 11/06/2018
% 
% The estimated upper-mantle Vp is also added.
% Tianze Liu, 11/08/2018
%
% The predicted SsPmp travel times are added.
% Tianze Liu, 09/06/2019
%
% The Vp_um_test is reversed to descending order so that in the case of
% pre-critical reflection the highest velocity that enables pre-critical
% reflection is found.
% Tianze Liu, 09/16/2019
%
% The results are also written to GMT folder
% Tianze Liu, 09/18/2019

clear all
close all

%% Parameters
% The inputs
path_in_mod = '/Users/tianze/Research/WorkingDirectory/VDSSsynthetics/ModelByPart_2D/Archive/TXT/FlatLAB_HarmonicMoho400_LoLABchange_Ang35_Stn10/RayTrac_mod.txt';
path_in_rayp_app = '/Users/tianze/Research/WorkingDirectory/VDSSsynthetics/ModelByPart_2D/Archive/TXT/FlatLAB_HarmonicMoho400_LoLABchange_Ang35_Stn10/Rayp_app.txt';
path_in_rayp_theo = '/Users/tianze/Research/WorkingDirectory/VDSSsynthetics/ModelByPart_2D/Archive/TXT/FlatLAB_HarmonicMoho400_LoLABchange_Ang35_Stn10/Rayp_theo.txt';
path_in_res = '/Users/tianze/Research/WorkingDirectory/VDSSsynthetics/ModelByPart_2D/Archive/TXT/FlatLAB_HarmonicMoho400_LoLABchange_Ang35_Stn10/TravelTimeRes.txt';
path_in_phs = '/Users/tianze/Research/WorkingDirectory/VDSSsynthetics/ModelByPart_2D/Archive/TXT/FlatLAB_HarmonicMoho400_LoLABchange_Ang35_Stn10/PhaseShift_obs.txt';

% The output folder
dirnm_out_arc = '/Users/tianze/Research/WorkingDirectory/VDSSsynthetics/ModelByPart_2D/Archive/TXT/FlatLAB_HarmonicMoho400_LoLABchange_Ang35_Stn10';
dirnm_out_gmt = '/Users/tianze/Research/WorkingDirectory/VDSSsynthetics/ModelByPart_2D/GMT/TXTFile/FlatLAB_HarmonicMoho400_LoLABchange_Ang35_Stn10';

% The origin for output X distance
x0 = 1400;

% The height of the free surface
z0 = 300;

% Crustal velocity
vp_cr = 6.5;

% The initial ray parameter
rayp0 = 0.1275;

% Velocities across the Moho
vp_lc = 6.5;
vs_lc = 3.75;
rho_lc = 2.83;

vp_um = 8.30;
vs_um = 4.61;
rho_um = 3.4;

% Velocity bounds for searching for upper-mantle Vp
vp_um_min = 7.8;
vp_um_max = 8.8;
dvp_um = -0.005;
kappa_um = 1.8;
Vp_um_test = vp_um_max:dvp_um:vp_um_min; % Reverse order so that the highest velocity that makes pre-critical reflection happen is found
n_vp_um = length(Vp_um_test);

%% Read the ray parameters
% The apparent ray parameters
Input = load(path_in_rayp_app);
X_rayp_app = Input(:,1);
Rayp_app = Input(:,2);
X_rayp_app = X_rayp_app+x0; % Consider the shift of x origin

% The theoretcial ray parameters
Input = load(path_in_rayp_theo);
X_rayp_theo = Input(:,1);
Rayp_theo = Input(:,2);
X_rayp_theo = X_rayp_theo+x0; % Consider the shift of x origin

% The measured phase shift
Input = load(path_in_phs);
X_phs_obs = Input(:,1);
Phs_obs = Input(:,2);

% The travel time residues
% Read the travel time residues
fid = fopen(path_in_res);
Input = textscan(fid,'%s %f %f %f');
X_res = Input{2};
Res = Input{3};

%% Read and plot the model
fid = fopen(path_in_mod);
Input = textscan(fid,'%d',1,'CommentStyle','#');
nlyr = Input{1};
Interface = struct('vp',cell(nlyr,1),'vs',cell(nlyr,1),'X',cell(nlyr,1),'Z',cell(nlyr,1));

for i = 1:nlyr
    Input = textscan(fid,'%f %f',1,'CommentStyle','#');
    vp = Input{1};
    vs = Input{2};
    Input = textscan(fid,'%d',1,'CommentStyle','#');
    npts = Input{1};
    Input = textscan(fid,'%f %f',npts,'CommentStyle','#');
    X_bdr = Input{1};
    Z_bdr = Input{2};
    
    Interface(i).vp = vp;
    Interface(i).vs = vs;
    Interface(i).X = X_bdr;
    Interface(i).Z = Z_bdr;
end
fclose(fid);

%% Perform ray-tracing and write results to file. 
% Trace the rays with the apparent ray parameters and compute the phase
% shifts
nx_app = length(X_rayp_app);
Phs_syn = zeros(nx_app,1);
X_phs_syn = zeros(nx_app,1);
fid1 = fopen(fullfile(dirnm_out_arc,'RayPath_PmP_post_app.txt'),'w');
fid2 = fopen(fullfile(dirnm_out_arc,'RayPath_PmP_pre_app.txt'),'w');
fid3 = fopen(fullfile(dirnm_out_gmt,'RayPath_PmP_post_app.txt'),'w');
fid4 = fopen(fullfile(dirnm_out_gmt,'RayPath_PmP_pre_app.txt'),'w');

fig_path_app = figure;
title('Pmp rays traced with apparent ray parameter')
hold on
for i = 1:nlyr
    X_bdr = Interface(i).X;
    Z_bdr = Interface(i).Z;
    plot(X_bdr,Z_bdr,'b');
end
xlabel('Distance (km)')
zlabel('Height (km)')

Vp_um_best = zeros(nx_app,1);
X_vp_um_out = zeros(nx_app,1);
T_vdss_app = zeros(nx_app,1);
X_srf_app = zeros(nx_app,1);
k = 0;

for i = 1:nx_app
    x_inc = X_rayp_app(i);
    rayp_inc = Rayp_app(i);
    ang_inc = asind(rayp_inc*vp_cr);
    
    if ang_inc >= 90
        continue
    else
        % Cacluate the ray path
        [X_pt,Z_pt,pcr,rayp_eff] = RayTracePmP_2Dlyr(x_inc,ang_inc,Interface);
        
        % Plot the ray path
        plot(X_pt,Z_pt,'r');
        n_pts = length(Z_pt);
        X_pt_out = X_pt-x0;
        Z_pt_out = z0-Z_pt;
        Output = [X_pt_out,Z_pt_out];
        
        % Compute the Pmp travel time
        l = sum(sqrt(diff(X_pt).^2+diff(Z_pt).^2));
        t_pmp = l/vp_cr;
        
        % Compute the Ss travel time difference between the source and
        % surfacing point
        x_inc_out = X_pt_out(1);
        x_srf_out = X_pt_out(3);
        X_srf_app(i) = x_srf_out;
        
        % The predicted SsPmp travel time
        if x_srf_out < max(X_res)
            res_inc = interp1(X_res,Res,x_inc_out);
            res_srf = interp1(X_res,Res,x_srf_out);
            dt_ss = (res_inc-res_srf)+rayp0*(x_inc_out-x_srf_out);
            t_vdss = dt_ss+t_pmp;
            T_vdss_app(i) = t_vdss;
        end
             
        % Calculate SsPmp phase shift with the correct crustal and
        % upper-mantle velocities
        phs_syn = PhaseShiftMoho(rayp_eff,vp_lc,vp_um,vs_lc,vs_um,rho_lc,rho_um);
        Phs_syn(i) = phs_syn;
        X_phs_syn(i) = X_pt_out(end);
        
        % Determine the upper-mantle Vp
        
        if x_srf_out < max(X_phs_obs)
            phs_obs = interp1(X_phs_obs,Phs_obs,x_srf_out);
            k = k+1;
            Msft = zeros(n_vp_um,1);
            for j = 1:n_vp_um
                vp_um_test = Vp_um_test(j);
                vs_um_test = vp_um_test/kappa_um;
                rho_um_test = Vp2Rho(vp_um_test);
                phs_syn_test = PhaseShiftMoho(rayp_eff,vp_lc,vp_um_test,vs_lc,vs_um_test,rho_lc,rho_um_test);

                msft = abs(phs_syn_test-phs_obs);
                Msft(j) = msft;       
            end
            [~,j_m] = min(Msft);
            vp_um_best = Vp_um_test(j_m);
            Vp_um_best(i) = vp_um_best;
            x_ref = X_pt_out(2);
            X_vp_um_out(i) = x_ref;
            disp(phs_obs)
            disp(rayp_eff)
            disp(vp_um_best)
        end
        % Write the ray path to file
        if phs_syn < 180
            fprintf(fid1,'%f %f\n',Output');
            fprintf(fid1,'>\n');
            fprintf(fid3,'%f %f\n',Output');
            fprintf(fid3,'>\n');
        else
            fprintf(fid2,'%f %f\n',Output');
            fprintf(fid2,'>\n');
            fprintf(fid4,'%f %f\n',Output');
            fprintf(fid4,'>\n');            
        end
    end
end
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);

X_vp_um_out = X_vp_um_out(1:k);
Vp_um_best = Vp_um_best(1:k);
T_vdss_app = T_vdss_app(1:k);
X_srf_app = X_srf_app(1:k);

% Plot the estimated upper-mantle Vp as funciton of distance
figure;
plot(X_vp_um_out,Vp_um_best);
xlabel('Distance (km)');
ylabel('Upper-mantle V_p (km/s)')

fid = fopen(fullfile(dirnm_out_arc,'PhaseShift_mod_app.txt'),'w');
Output = [X_phs_syn,Phs_syn];
fprintf(fid,'%f %f\n',Output');
fclose(fid);

fid = fopen(fullfile(dirnm_out_gmt,'PhaseShift_mod_app.txt'),'w');
Output = [X_phs_syn,Phs_syn];
fprintf(fid,'%f %f\n',Output');
fclose(fid);

% Ray tracing with the theoretical ray parameters
nx_theo = length(X_rayp_theo);
Phs_theo = zeros(nx_theo,1);
X_phs_theo = zeros(nx_theo,1);
fid1 = fopen(fullfile(dirnm_out_arc,'RayPath_PmP_post_theo.txt'),'w');
fid2 = fopen(fullfile(dirnm_out_arc,'RayPath_PmP_pre_theo.txt'),'w');
fid3 = fopen(fullfile(dirnm_out_gmt,'RayPath_PmP_post_theo.txt'),'w');
fid4 = fopen(fullfile(dirnm_out_gmt,'RayPath_PmP_pre_theo.txt'),'w');

fig_path_theo = figure;
title('Pmp rays traced with theoretical ray parameter')
hold on
for i = 1:nlyr
    X_bdr = Interface(i).X;
    Z_bdr = Interface(i).Z;
    plot(X_bdr,Z_bdr,'b');
end
xlabel('Distance (km)')
zlabel('Height (km)')

for i = 1:nx_theo
    x_inc = X_rayp_theo(i);
    rayp_inc = Rayp_theo(i);
    ang_inc = asind(rayp_inc*vp_cr);
    
    if ang_inc >= 90
        continue
    else
        % Cacluate the ray path
        [X_pt,Z_pt,pcr,rayp_eff] = RayTracePmP_2Dlyr(x_inc,ang_inc,Interface);

        % Plot the ray path
        plot(X_pt,Z_pt,'r');
        n_pts = length(Z_pt);
        X_pt_out = X_pt-x0;
        Z_pt_out = z0-Z_pt;
        Output = [X_pt_out,Z_pt_out];

        % Calculate SsPmp phase shift
        phs = PhaseShiftMoho(rayp_eff,vp_lc,vp_um,vs_lc,vs_um,rho_lc,rho_um);
        Phs_syn(i) = phs;
        X_phs_syn(i) = X_pt_out(end);
                       
        % Write the ray path to file
        if phs < 180
            fprintf(fid1,'%f %f\n',Output');
            fprintf(fid1,'>\n');
            fprintf(fid3,'%f %f\n',Output');
            fprintf(fid3,'>\n');            
        else
            fprintf(fid2,'%f %f\n',Output');
            fprintf(fid2,'>\n');
            fprintf(fid4,'%f %f\n',Output');
            fprintf(fid4,'>\n');            
        end
    end
end
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);

fid = fopen(fullfile(dirnm_out_arc,'PhaseShift_mod_theo.txt'),'w');
Output = [X_phs_syn,Phs_syn];
fprintf(fid,'%f %f\n',Output');
fclose(fid);

% Plot the predicted T_vdss
fig_tvdss = figure;
plot(X_srf_app,T_vdss_app);

% Write the predicted T_vdss to files
fid = fopen(fullfile(dirnm_out_arc,'Tvdss_app.txt'),'w');
Out = [X_srf_app,T_vdss_app];
fprintf(fid,'%f %f\n',Out');
fclose(fid);

fid = fopen(fullfile(dirnm_out_gmt,'Tvdss_app.txt'),'w');
Out = [X_srf_app,T_vdss_app];
fprintf(fid,'%f %f\n',Out');
fclose(fid);

% Write the computed Vp_um to file
% Write results to file
fid = fopen(fullfile(dirnm_out_arc,'Vp_um_RayTrac.txt'),'w');
Output = [X_vp_um_out,Vp_um_best];
fprintf(fid,'%f %f\n',Output');
fclose(fid);

fid = fopen(fullfile(dirnm_out_gmt,'Vp_um_RayTrac.txt'),'w');
Output = [X_vp_um_out,Vp_um_best];
fprintf(fid,'%f %f\n',Output');
fclose(fid);
