%% Apply 2D Ray-tracing of plane S wave through a layered model
%
% History:
% Created.
% Tianze Liu, 03/26/2018
%
% The output ray paths are converted from height to depth.
% Tianze Liu, 04/11/2018

clear all
close all

%% Parameters
% The input model
path_in = '/Users/tianze/Research/WorkingDirectory/VDSSsynthetics/ModelByPart_2D/Archive/TXT/MohoMono+3deg_Ang35_Stn10/RayTrac_mod.txt';

% The output folder
dirnm_out = '/Users/tianze/Research/WorkingDirectory/VDSSsynthetics/ModelByPart_2D/Archive/TXT/MohoMono+3deg_Ang35_Stn10';

% The incident angle (in deg)
ang_inc = 35;

% The starting point of incident rays
X_inc = 1700:10:2800;

% The origin for output X distance
x0 = 1600;

% The height of the free surface
z0 = 400;


%% Read and plot the model
fid = fopen(path_in);
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

fig_path = figure;
hold on
for i = 1:nlyr
    X_bdr = Interface(i).X;
    Z_bdr = Interface(i).Z;
    plot(X_bdr,Z_bdr,'b');
end
xlabel('Distance (km)')
zlabel('Height (km)')

%% Perform ray-tracing and write results to file. Calculate ray parameters and also write it to file
nx = length(X_inc);
Rayp = zeros(nx,1);
X_srf = zeros(nx,1);

fid = fopen(fullfile(dirnm_out,'RayPath_S.txt'),'w');
for i = 1:nx
    x_inc = X_inc(i);
    
    % Cacluate the ray path
    [X_pt,Z_pt] = RayTraceS_2Dlyr(x_inc,ang_inc,Interface);
   
    plot(X_pt,Z_pt,'r');
    n_pts = length(Z_pt);
    ang_srf = acotd((Z_pt(n_pts)-Z_pt(n_pts-1))/(X_pt(n_pts)-X_pt(n_pts-1)));
    vs_srf = Interface(nlyr-1).vs;
    rayp = sind(ang_srf)/vs_srf;
    x_srf = X_pt(n_pts);
    Rayp(i) = rayp;
    X_srf(i) = x_srf;
    
    % Write the ray path to file
    X_pt_out = X_pt-x0;
    Z_pt_out = z0-Z_pt;
    Output = [X_pt_out,Z_pt_out];
    fprintf(fid,'%f %f\n',Output');
    fprintf(fid,'>\n');    
end
fclose(fid);

fig_rayp = figure;
%scatter(X_srf,Rayp,'filled');
plot(X_srf,Rayp);

% Write ray parameters to file
fid = fopen(fullfile(dirnm_out,'Rayp_theo.txt'),'w');
X_srf_out = X_srf-x0;
Rayp_out = Rayp;
Output = [X_srf_out,Rayp];
fprintf(fid,'%f %f\n',Output');



