%% Backproject 2D Observed SsPmp envelope function
%
% History:
% Modified from BackProject_gauss_Syn2D.m.
% Tianze Liu, 10/11/2018
%
% The energy on the observed traces are retrieved with linear
% interpolation as opposed to nearest-neighbor interpolation 
% Tianze Liu, 09/20/2019
%
% The phase-corrected traces are back-projected instead of the envelope
% functions.
% Tianze Liu, 11/1/2019

clear all
close all

%% Parameters
dirnm_in_sac = '/Volumes/TianzeResearch/OrdosChinArray/Archive/SAC/2017-080-23-10-25/PS_shift';
dirnm_in_txt = '/Volumes/TianzeResearch/OrdosChinArray/Archive/TXT/2017-080-23-10-25/Line1';

fnm_in_sac = 'SACinfo_p_raw_phs0_int.txt';
fnm_in_res = 'ResInt.txt';
fnm_in_rayp = 'Rayp_prj.txt';
fnm_in_rayp0 = 'Rayp0_prj.txt';

dirnm_out_grd = '/Volumes/TianzeResearch/OrdosChinArray/GMT/GRDFile/2017-080-23-10-25/Line1';
dirnm_out_txt = '/Volumes/TianzeResearch/OrdosChinArray/GMT/TXTFile/2017-080-23-10-25/Line1';

suf_phs = 'phs0_int';

d_src = 10;
vp = 6.0;

dx_dat = 1; % Grid spacing for data interpolation

x_src_min = 50; % Starting point of ray tracing 
x_src_max = 500;
dx_src = 10;

z_samp_min = 1;
z_samp_max = 70;
dz_samp = 0.5;

x_samp_min = 0;
x_samp_max = 600;
dx_samp = 1;

dx_ray = 1; % The X spacing between virtual source and receiver

dx_grd = 10;

% The range for Moho estimation
x_moho_min = 180;
x_moho_max = 430;

f_dm = 0.15; % The dominant frequency used to calculate Fresnel Zone width

pts = 1; % 1 means wave incident from left to right, -1 means right to left

c_max = 1; % The maximum in color bar

theta = 0; % The assumed interface dip

%% Read data
fid = fopen(fullfile(dirnm_in_txt,fnm_in_sac));
Input = textscan(fid,'%s %s %f %f %f %f %f');
Sacnm_in = Input{1};
X_sac = Input{5};
Sacst = SACST_fread(Sacnm_in,'prefix',[dirnm_in_sac,'/']);
dt = Sacst(1).delta;
tb = Sacst(1).b;
npts_in = Sacst(1).npts;
te = tb+dt*(npts_in-1);
n_sac = length(Sacst);

fid = fopen(fullfile(dirnm_in_txt,fnm_in_res));
Input = textscan(fid,'%f %f');
X_res = Input{1};
Res = Input{2};
Res0 = zeros(size(X_res));

fid = fopen(fullfile(dirnm_in_txt,fnm_in_rayp0));
Input = textscan(fid,'%f %f');
X_rayp0 = Input{1};
Rayp0 = Input{2};

fid = fopen(fullfile(dirnm_in_txt,fnm_in_rayp));
Input = textscan(fid,'%f %f');
X_rayp = Input{1};
Rayp = Input{2};

%% Bin data
x_grd_min = floor(min(X_sac)/dx_grd)*dx_grd;
x_grd_max = ceil(max(X_sac)/dx_grd)*dx_grd;

x_bin_min = x_grd_min+dx_grd/2;
x_bin_max = x_grd_max-dx_grd/2;

dx_bin = dx_grd;

X_grd = x_grd_min:dx_grd:x_grd_max;
X_bin = x_bin_min:dx_bin:x_bin_max;

nx_grd = length(X_grd);
nx_bin = nx_grd-1;

Data_bin = struct('data',{},'x',{});
for i = 1:nx_bin
    x_bin = X_bin(i);
    x_lo = X_grd(i);
    x_hi = X_grd(i+1);
    
    Sacst_bin = Sacst(X_sac < x_hi & X_sac > x_lo);
    if ~isempty(Sacst_bin)
        n_bin = length(Sacst_bin);
        P_bin = zeros(npts_in,1);
        for j = 1:n_bin
            sacst = Sacst_bin(j);
            P = sacst.data;
            P_bin = P_bin+P;
        end
        P_bin = P_bin/n_bin;
        data_bin.data = P_bin;
        data_bin.x = x_bin;
        Data_bin = [Data_bin,data_bin];
    end
end

n_dat_bin = length(Data_bin);
X_int = zeros(npts_in*n_dat_bin,1);
T_int = zeros(npts_in*n_dat_bin,1);
Data_int = zeros(npts_in*n_dat_bin,1);
k = 0;
for i = 1:n_dat_bin
    data_bin = Data_bin(i);
    x = data_bin.x;
    P_bin = data_bin.data;
       
    for j = 1:npts_in
        k = k+1;
        X_int(k) = x;
        T_int(k) = tb+(j-1)*dt;
        Data_int(k) = P_bin(j);
    end
end

% Interpolate to make the data matrix
T_bin = tb:dt:te;

[XX_bin,TT_bin] = meshgrid(X_bin,T_bin);
Img_dat = griddata(X_int,T_int,Data_int,XX_bin,TT_bin,'cubic');
Img_dat = Img_dat/max(max(Img_dat));

fig1 = figure;
colormap(parula)
image([x_bin_min,x_bin_max],[tb,te],Img_dat,'CDataMapping','scaled');
caxis([-1,1])
xlabel('Distance (km)')
ylabel('Time (s)')

%% Trace the rays for both the corrected and the original case
X_src = x_src_min:dx_src:x_src_max;

X_samp = x_samp_min:dx_samp:x_samp_max;
Z_samp = z_samp_min:dz_samp:z_samp_max;
n_src = length(X_src);

nx_samp = length(X_samp);
nz_samp = length(Z_samp);

% For the corrected case
disp('Back-projection for the uncorrected case...')
Img_mod = zeros(nz_samp,nx_samp);
for i = 1:n_src
    x_src = X_src(i);
    rayp_src = interp1(X_rayp,Rayp,x_src,'linear');
    disp(rayp_src);
    sin_i = vp*rayp_src;
    cos_i = sqrt(1-sin_i^2);
    
    if isnan(rayp_src)
        continue
    else
        for j = 1:nz_samp
            z_samp = Z_samp(j);
            
            l = z_samp*(1/cos_i+1/(cosd(2*theta)*cos_i-sind(2*theta)*sin_i));
            t_pmp = l/vp;
            d1 = z_samp*sin_i/cos_i;
            d2 = z_samp*(sind(2*theta)*cos_i+cosd(2*theta)*sin_i)/(cosd(2*theta)*cos_i-sind(2*theta)*sin_i);
            x_rec = x_src+pts*(d1+d2);
            x_mdpt = x_src+pts*d1;
            
            res_src = interp1(X_res,Res,x_src);
            
            if x_rec > max(X_res) || x_rec < min(X_res)
                continue
            else
                res_rec = interp1(X_res,Res,x_rec);
            end
            
            X_ray = (x_src:dx_ray:x_rec)';
            Rayp0_ray = interp1(X_rayp0,Rayp0,X_ray);
            dt_ss = (res_src-res_rec)-Rayp0_ray'*dx_ray*ones(size(X_ray));

            t_vdss = dt_ss+t_pmp;
            %disp(t_vdss)
%             it_amp = round((t_vdss-tb)/dt);
%             [~,ix_amp] = min(abs(X_sac-x_rec));
%             sacst_rec = Sacst(ix_amp);
%             P_rec = sacst_rec.data;
%             npts = length(P_rec);

            if t_vdss > te || x_rec > x_bin_max
                amp = 0;
            else
                amp = interp2(XX_bin,TT_bin,Img_dat,x_rec,t_vdss);
            end
%             disp(amp)

            %ix_samp = floor((x_mdpt-x_samp_min)/dx_samp+1); % Length of the ray path
            lambda = vp/f_dm;
            d_fn = sqrt(l*lambda)/2; % Half width of the Fresnel Zone.
            Img_mod(j,:) = Img_mod(j,:)+amp*exp(-(X_samp-x_mdpt).^2/2/d_fn^2);
        end
    end
end

Img_mod = Img_mod/max(max(Img_mod));

% For the uncorrected case
disp('Back-projection for the uncorrected case...')
Img0_mod = zeros(nz_samp,nx_samp);
for i = 1:n_src
    x_src = X_src(i);
    rayp_src = interp1(X_rayp0,Rayp0,x_src,'linear');
    disp(rayp_src);
    
    if isnan(rayp_src)
        continue
    else
        for j = 1:nz_samp
            z_samp = Z_samp(j);

            l = z_samp*(1/cos_i+1/(cosd(2*theta)*cos_i-sind(2*theta)*sin_i));
            t_pmp = l/vp;
            d1 = z_samp*sin_i/cos_i;
            d2 = z_samp*(sind(2*theta)*cos_i+cosd(2*theta)*sin_i)/(cosd(2*theta)*cos_i-sind(2*theta)*sin_i);
            x_rec = x_src+pts*(d1+d2);
            x_mdpt = x_src+pts*d1;
            res_src = interp1(X_res,Res0,x_src);
            
            if x_rec > max(X_res) || x_rec < min(X_res)
                continue
            else
                res_rec = interp1(X_res,Res0,x_rec);
            end
            
            X_ray = (x_src:dx_ray:x_rec)';
            Rayp0_ray = interp1(X_rayp0,Rayp0,X_ray);
            dt_ss = (res_src-res_rec)-Rayp0_ray'*dx_ray*ones(size(X_ray));

            t_vdss = dt_ss+t_pmp;
            it_amp = round((t_vdss-tb)/dt);
            
            if t_vdss > te || x_rec > x_bin_max
                amp = 0;
            else
                amp = interp2(XX_bin,TT_bin,Img_dat,x_rec,t_vdss);
            end
            
            lambda = vp/f_dm;
            d_fn = sqrt(l*lambda)/2; % Half width of the Fresnel Zone.
            Img0_mod(j,:) = Img0_mod(j,:)+amp*exp(-(X_samp-x_mdpt).^2/2/d_fn^2);
        end
    end
end

Img0_mod = Img0_mod/max(max(Img0_mod));

% Display results
x_out_min = x_samp_min+dx_samp/2;
x_out_max = x_samp_max+dx_samp/2;
dx_out = dx_samp;
X_out = x_out_min:dx_out:x_out_max;
nx_out = length(X_out);

z_out_min = z_samp_min;
z_out_max = z_samp_max;
dz_out = dz_samp;
Z_out = Z_samp;
nz_out = nz_samp;

%Find the Moho depth for the corrected case
Z_moho = zeros(nx_out,1);
for i = 1:nx_out
    [~,j_z] = max(Img_mod(:,i));
    z_moho = Z_out(j_z);
    Z_moho(i) = z_moho;
end
i_moho_min = floor((x_moho_min-x_samp_min)/dx_samp+1);
i_moho_max = floor((x_moho_max-x_samp_min)/dx_samp+1);
Z_moho_out = Z_moho(i_moho_min:i_moho_max);
X_moho_out = (x_moho_min:dx_samp:x_moho_max)';

% Find the Moho depth for the uncorrected case
Z_moho0 = zeros(nx_out,1);
for i = 1:nx_out
    [~,j_z] = max(Img0_mod(:,i));
    z_moho = Z_out(j_z);
    Z_moho0(i) = z_moho;
end
i_moho_min = floor((x_moho_min-x_samp_min)/dx_samp+1);
i_moho_max = floor((x_moho_max-x_samp_min)/dx_samp+1);
Z_moho0_out = Z_moho0(i_moho_min:i_moho_max);
X_moho0_out = (x_moho_min:dx_samp:x_moho_max)';

% The corrected image
fig2 = figure;
colormap(parula)
image([x_out_min,x_out_max],[z_out_min,z_out_max],Img_mod,'CDataMapping','scaled');
caxis([-1,1])
hold on
plot(X_moho_out,Z_moho_out,'red')

% The uncorrected image
fig3 = figure;
colormap(parula)
image([x_out_min,x_out_max],[z_out_min,z_out_max],Img0_mod,'CDataMapping','scaled');
caxis([-1,1])
hold on
plot(X_moho0_out,Z_moho0_out,'red')

%% Output images
disp('Writing files...')
fnm_out_grd = ['BackProjectImg_',suf_phs,'_gauss_',num2str(vp,'%.1f'),'.grd'];
fnm_out_grd0 = ['BackProjectImg0_',suf_phs,'_gauss_',num2str(vp,'%.1f'),'.grd'];

% Remove the original NetCDF files
unix(['rm ',fullfile(dirnm_out_grd,fnm_out_grd)]);
unix(['rm ',fullfile(dirnm_out_grd,fnm_out_grd0)]);

% Write NetCDF file
nccreate(fullfile(dirnm_out_grd,fnm_out_grd),'amplitude','Dimensions',{'x',nx_out,'z',nz_out});
nccreate(fullfile(dirnm_out_grd,fnm_out_grd),'x','Dimensions',{'x',nx_out});
nccreate(fullfile(dirnm_out_grd,fnm_out_grd),'z','Dimensions',{'z',nz_out});

ncwrite(fullfile(dirnm_out_grd,fnm_out_grd),'x',X_out)
ncwrite(fullfile(dirnm_out_grd,fnm_out_grd),'z',Z_out)
ncwrite(fullfile(dirnm_out_grd,fnm_out_grd),'amplitude',Img_mod')

nccreate(fullfile(dirnm_out_grd,fnm_out_grd0),'amplitude','Dimensions',{'x',nx_out,'z',nz_out});
nccreate(fullfile(dirnm_out_grd,fnm_out_grd0),'x','Dimensions',{'x',nx_out});
nccreate(fullfile(dirnm_out_grd,fnm_out_grd0),'z','Dimensions',{'z',nz_out});

ncwrite(fullfile(dirnm_out_grd,fnm_out_grd0),'x',X_out)
ncwrite(fullfile(dirnm_out_grd,fnm_out_grd0),'z',Z_out)
ncwrite(fullfile(dirnm_out_grd,fnm_out_grd0),'amplitude',Img0_mod')

%% Output Moho picks
fnm_out_txt = ['MohoDepth_',suf_phs,'_est_',num2str(vp,'%.1f'),'.txt'];
fnm_out_txt0 = ['MohoDepth0_',suf_phs,'_est_',num2str(vp,'%.1f'),'.txt'];

fid = fopen(fullfile(dirnm_out_txt,fnm_out_txt),'w');
Output = [X_moho_out,Z_moho_out];
fprintf(fid,'%f %f\n',Output');
fclose(fid);

fid = fopen(fullfile(dirnm_out_txt,fnm_out_txt0),'w');
Output = [X_moho0_out,Z_moho0_out];
fprintf(fid,'%f %f\n',Output');
fclose(fid);