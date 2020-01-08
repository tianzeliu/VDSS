%% Plot the trade-off curve between Vp/Vs ratio and crustal thickness derived from Ps arrival time assuming a fixed Vp
%
% History:
% Created.
% Tianze Liu, 01/09/2019
%
% The joint distribution of Vp and Vs is also output.
% Tianze Liu, 08/01/2019

close all
clear all

%% Inputs and outputs
% The input and output directory
dirnm_in = '/Users/tianze/Research/WorkingDirectory/POLARIS/RF/Tianze/Archive/TXT/EDZN';
dirnm_out_txt = '/Users/tianze/Research/WorkingDirectory/POLARIS/RF/Tianze/GMT/TXTFile';
dirnm_out_grd = '/Users/tianze/Research/WorkingDirectory/POLARIS/RF/Tianze/GMT/GRDFile';

% The input file containing the random simulaiton of H and Vp
fnm_in = 'HVp_rand5000_msft.txt';

% The output file containing the random simulation of k
fnm_out_txt_k_rnd = 'Kappa_rand5000.txt';
fnm_out_txt_k_best = 'Kappa_joint_mean.txt';
fnm_out_txt_k_hi = 'Kappa_joint_max.txt';
fnm_out_txt_k_lo = 'Kappa_joint_min.txt';
fnm_out_txt_h = 'Moho_VDSS.txt';

% The output file containing the joint distribution of random simulations of Vp and k
fnm_out_grd_k = 'VpKappa_rand5000_Pdf.grd';
fnm_out_grd_vs = 'VpVs_rand5000_Pdf.grd';

% The parameters 
t0 = 4.24; % The input Ps arrival time
vp = 6.67; % The input average crustal Vp
rayp = 0.0; % The ray parameter

h_min = 30;
h_max = 44;
dh = 0.1;

h_best = 38.2;

k_hist_min = 1.60;
k_hist_max = 2.00;
dk_hist = 0.01;

vp_hist_min = 5.5;
vp_hist_max = 7.6;
dvp_hist = 0.02;

vs_hist_min = 2.8;
vs_hist_max = 4.9;
dvs_hist = 0.02;

pr_thrsd = 0.68;

% The values used for the constant K curves
k1 = 1.7;
k2 = 1.75;
k3 = 1.8;

fnm_out_txt_k1 = ['ConstantK_',num2str(k1,'%.2f'),'.txt'];
fnm_out_txt_k2 = ['ConstantK_',num2str(k2,'%.2f'),'.txt'];
fnm_out_txt_k3 = ['ConstantK_',num2str(k3,'%.2f'),'.txt'];

%% Read in the random simulation of H and Vp
Input = load(fullfile(dirnm_in,fnm_in));
H_rnd = Input(:,1);
Vp_rnd = Input(:,2);

n_rnd = length(Vp_rnd);

%% Compute the H-k curve and find the k value for the best-fit h 
H = (h_min:0.1:h_max)';
K = sqrt((vp*t0./H+sqrt(1-vp^2*rayp^2)).^2+vp^2*rayp^2);
scatter(H,K,'filled')

k_best = interp1(H,K,h_best);

disp(['The best-fit k is ',num2str(k_best,'%.3f'),])

%% Compute the k and Vs values for the random simulations of H and Vp
K_rnd = sqrt((Vp_rnd*t0./H_rnd+sqrt(1-Vp_rnd.^2*rayp^2)).^2+Vp_rnd.^2*rayp^2);
Vs_rnd = Vp_rnd./K_rnd;

%% Compute the nominal standard error for Kappa
K_edge = k_hist_min:dk_hist:k_hist_max;
n_edge = length(K_edge);
[N_hist_k,~] = histcounts(K_rnd,K_edge);
Pr_hist_k = N_hist_k/n_rnd;

k_std = sqrt(var(K_rnd));
k_hi = k_best+k_std;
k_lo = k_best-k_std;
display(['The standard erorr of Kappa is ',num2str(k_std,'%.3f')]);

%% Plot the histogram of the random simulations of K
figure;
histogram(K_rnd,K_edge,'Normalization','probability');

%% Plot the Vp-K and Vp-Vs distribution 
Vp_edge = vp_hist_min:dvp_hist:vp_hist_max;
Vs_edge = vs_hist_min:dvs_hist:vs_hist_max;

[Pdf_hist_VpK,~] = histcounts2(K_rnd,Vp_rnd,K_edge,Vp_edge,'Normalization','pdf');
[Pdf_hist_VpVs,~] = histcounts2(Vs_rnd,Vp_rnd,Vs_edge,Vp_edge,'Normalization','pdf');

Vp_ax = Vp_edge(1:end-1)+diff(Vp_edge)/2;
vp_ax_max = max(Vp_ax);
vp_ax_min = min(Vp_ax);
nvp_ax = length(Vp_ax);

Vs_ax = Vs_edge(1:end-1)+diff(Vs_edge)/2;
vs_ax_max = max(Vs_ax);
vs_ax_min = min(Vs_ax);
nvs_ax = length(Vs_ax);

K_ax = K_edge(1:end-1)+diff(K_edge)/2;
k_ax_max = max(K_ax);
k_ax_min = min(K_ax);
nk_ax = length(K_ax);

figure;
colormap(hot);
image([vp_ax_min,vp_ax_max],[k_ax_min,k_ax_max],Pdf_hist_VpK,'CDataMapping','scaled');
colorbar

figure;
colormap(hot);
image([vp_ax_min,vp_ax_max],[vs_ax_min,vs_ax_max],Pdf_hist_VpVs,'CDataMapping','scaled');
colorbar

%% Output the results to TXT files
% The H-k curve
Output = [H,K];
fid = fopen(fullfile(dirnm_out_txt,'DepthKappaCurve.txt'),'w');
fprintf(fid,'%f %f\n',Output');
fclose(fid);

% The values for Vp/Vs ratios
fid = fopen(fullfile(dirnm_out_txt,fnm_out_txt_k_best),'w');
fprintf(fid,'0.0 %f\n',k_best);
fprintf(fid,'100.0 %f',k_best);
fclose(fid);

fid = fopen(fullfile(dirnm_out_txt,fnm_out_txt_k_hi),'w');
fprintf(fid,'0.0 %f\n',k_hi);
fprintf(fid,'100.0 %f',k_hi);
fclose(fid);

fid = fopen(fullfile(dirnm_out_txt,fnm_out_txt_k_lo),'w');
fprintf(fid,'0.0 %f\n',k_lo);
fprintf(fid,'100.0 %f',k_lo);
fclose(fid);

% The values for crustal thickness
fid = fopen(fullfile(dirnm_out_txt,fnm_out_txt_h),'w');
fprintf(fid,'%f 0.0\n',h_best);
fprintf(fid,'%f 10.0',h_best);
fclose(fid);

% The random simulation of k
fid = fopen(fullfile(dirnm_out_txt,fnm_out_txt_k_rnd),'w');
fprintf(fid,'%.3f\n',K_rnd);
fclose(fid);

%% Output the results to GRD file
% Output the Vp-K image
unix(['rm ',fullfile(dirnm_out_grd,fnm_out_grd_k)]);
Pdf_hist_VpK = Pdf_hist_VpK'; % When plotting with GMT, Kappa will be Y axis and Vp will be X axis

nccreate(fullfile(dirnm_out_grd,fnm_out_grd_k),'Probability','Dimensions',{'Vp',nvp_ax,'Kappa',nk_ax});
nccreate(fullfile(dirnm_out_grd,fnm_out_grd_k),'Kappa','Dimensions',{'Kappa',nk_ax});
nccreate(fullfile(dirnm_out_grd,fnm_out_grd_k),'Vp','Dimensions',{'Vp',nvp_ax});

ncwrite(fullfile(dirnm_out_grd,fnm_out_grd_k),'Probability',Pdf_hist_VpK)
ncwrite(fullfile(dirnm_out_grd,fnm_out_grd_k),'Kappa',K_ax)
ncwrite(fullfile(dirnm_out_grd,fnm_out_grd_k),'Vp',Vp_ax)

% Output the Vp-Vs image
unix(['rm ',fullfile(dirnm_out_grd,fnm_out_grd_vs)]);
Pdf_hist_VpVs = Pdf_hist_VpVs'; % When plotting with GMT, Vs will be Y axis and Vp will be X axis

nccreate(fullfile(dirnm_out_grd,fnm_out_grd_vs),'Probability','Dimensions',{'Vp',nvp_ax,'Vs',nvs_ax});
nccreate(fullfile(dirnm_out_grd,fnm_out_grd_vs),'Vs','Dimensions',{'Vs',nvs_ax});
nccreate(fullfile(dirnm_out_grd,fnm_out_grd_vs),'Vp','Dimensions',{'Vp',nvp_ax});

ncwrite(fullfile(dirnm_out_grd,fnm_out_grd_vs),'Probability',Pdf_hist_VpVs)
ncwrite(fullfile(dirnm_out_grd,fnm_out_grd_vs),'Vs',Vs_ax)
ncwrite(fullfile(dirnm_out_grd,fnm_out_grd_vs),'Vp',Vp_ax)

%% Output the constant K lines
Vs_k1 = Vp_ax/k1;
Vs_k2 = Vp_ax/k2;
Vs_k3 = Vp_ax/k3;

fid =fopen(fullfile(dirnm_out_txt,fnm_out_txt_k1),'w');
Out = [Vp_ax',Vs_k1'];
fprintf(fid,'%f %f\n',Out');
fclose(fid);

fid =fopen(fullfile(dirnm_out_txt,fnm_out_txt_k2),'w');
Out = [Vp_ax',Vs_k2'];
fprintf(fid,'%f %f\n',Out');
fclose(fid);

fid =fopen(fullfile(dirnm_out_txt,fnm_out_txt_k3),'w');
Out = [Vp_ax',Vs_k3'];
fprintf(fid,'%f %f\n',Out');
fclose(fid);
