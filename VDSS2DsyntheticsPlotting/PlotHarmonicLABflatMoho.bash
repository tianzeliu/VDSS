#!/bin/bash
gmt gmtset PS_MEDIA 10iX14i
gmt gmtset FONT_ANNOT_PRIMARY 8p,0,black # lat/lon, scale bar values, scale bar title
gmt gmtset FONT_LABEL 10p # scale bar label
gmt gmtset FONT_TITLE 12p # Titles

gmt gmtset MAP_FRAME_PEN 1p,black
gmt gmtset MAP_FRAME_TYPE plain
gmt gmtset MAP_VECTOR_SHAPE 1
gmt gmtset PS_CHAR_ENCODING Standard+

gmt gmtset PS_PAGE_ORIENTATION portrait

# Set the grid lines as dash lines
gmt gmtset MAP_GRID_PEN_PRIMARY 0.001i,gray,0.05_0.05:0i

# Input and output
dir_ps=../PSFile
dir_txt=../TXTFile/HarmonicLAB400_FlatMoho_LoLABchange_Ang35_Stn10
dir_grd=../GRDFile/HarmonicLAB400_FlatMoho_LoLABchange_Ang35_Stn10
dir_pdf=../PDFFile
dir_cpt=../CPTFile

ps=$dir_ps/VDSS_HarmonicLAB400_FlatMoho_LoLABchange_Ang35_phs90.ps

# Model blocks
txt_at=$dir_txt/Asthenosphere_mod.txt
txt_lt=$dir_txt/Lithosphere_mod.txt
txt_cr=$dir_txt/Crust_mod.txt

# Receiver coordinate
txt_rec=$dir_txt/Rec.txt

# Moho estimations
txt_moho0=$dir_txt/MohoDepth_true.txt
txt_moho1=$dir_txt/MohoDepth0_phs90_est.txt
txt_moho2=$dir_txt/MohoDepth_phs90_est.txt
txt_moho3=$dir_txt/MohoDepth_phs90_stk_est.txt
#txt_moho2=$dir_txt/MohoDepth_env_raw.txt
#txt_moho3=$dir_txt/MohoDepth_env_corr.txt
#txt_moho3=$dir_txt/MohoDepth_env_CorrPro.txt

# Data image
grd1_dat=$dir_grd/Array_wiggle_z_theo.grd
grd2_dat=$dir_grd/Array_wiggle_z_corr.grd

# Back projection image
grd=$dir_grd/BackProjectImg_phs90_gauss.grd
grd_stk=$dir_grd/BackProjectImg_phs90_gauss_stk.grd

# Vector orientation
txt_vct=$dir_txt/Vct_35.txt

# The wave incident angle
ang=35

# Travel time residues
txt_res=$dir_txt/TravelTimeRes.txt

# Phase shifts
txt_phs0=$dir_txt/PhaseShift0.txt
txt_phs1=$dir_txt/PhaseShift_mod_app.txt
txt_phs2=$dir_txt/PhaseShift_obs.txt

# Upper-mantle Vp
txt_vp0=$dir_txt/Vp_um_true.txt
txt_vp_est=$dir_txt/Vp_um_RayTrac.txt
txt_vp_hi=$dir_txt/Vp_um_hi.txt
txt_vp_lo=$dir_txt/Vp_um_lo.txt

# The hyperbolae
txt_hp_s=$dir_txt/Hyperbola_s.txt
txt_hp_p=$dir_txt/Hyperbola_p.txt

# X limits of the box
txt_bx=$dir_txt/BoxXlim.txt
x_min_bx=$(awk 'NR==1{print $1}' $txt_bx)
x_max_bx=$(awk 'NR==1{print $2}' $txt_bx)

# Stations to be plotted individually
txt_stn=$dir_txt/StnmSpec.txt
stnm1=$(awk 'NR==1{print $1}' $txt_stn)
stnm2=$(awk 'NR==2{print $1}' $txt_stn)
stnm3=$(awk 'NR==3{print $1}' $txt_stn)

# Wiggle files
txt_wg1=$dir_txt/Array_wiggle_z_theo.txt
txt_wg1_s1=$dir_txt/${stnm1}_wiggle_z_theo.txt
txt_wg1_s2=$dir_txt/${stnm2}_wiggle_z_theo.txt
txt_wg1_s3=$dir_txt/${stnm3}_wiggle_z_theo.txt
txt_wg1_s4=$dir_txt/${stnm4}_wiggle_z_theo.txt
txt_wg1_s5=$dir_txt/${stnm5}_wiggle_z_theo.txt
txt_wg2=$dir_txt/Array_wiggle_z_corr.txt
txt_wg2_s1=$dir_txt/${stnm1}_wiggle_z_corr.txt
txt_wg2_s2=$dir_txt/${stnm2}_wiggle_z_corr.txt
txt_wg2_s3=$dir_txt/${stnm3}_wiggle_z_corr.txt
txt_wg2_s4=$dir_txt/${stnm4}_wiggle_z_corr.txt
txt_wg2_s5=$dir_txt/${stnm5}_wiggle_z_corr.txt

# Waveform files
txt_wv1=$dir_txt/${stnm1}_wvfm_z_corr.txt
txt_wv2=$dir_txt/${stnm2}_wvfm_z_corr.txt
txt_wv3=$dir_txt/${stnm3}_wvfm_z_corr.txt

# Ray path files
txt_path_s=$dir_txt/RayPath_S.txt
txt_path_pre_theo=$dir_txt/RayPath_PmP_pre_theo.txt
txt_path_post_theo=$dir_txt/RayPath_PmP_post_theo.txt
txt_path_pre_app=$dir_txt/RayPath_PmP_pre_app.txt
txt_path_post_app=$dir_txt/RayPath_PmP_post_app.txt

# Ray parameter files
txt_rayp0=$dir_txt/Rayp0.txt
txt_rayp1=$dir_txt/Rayp_theo.txt
txt_rayp2=$dir_txt/Rayp_app.txt

# CPT file
cpt=${dir_cpt}/BackProjPol.cpt

# Plotting dimensions
it_wg=-0.06i
it_wv=0.14i
ia=0.3i # Waveform amplitude
ix_bg=0.002i
ix_sm=0.004i
iz_bg=-0.004i # Z scale for the velocity model
iz_bp_sm=-0.028i # Z scale for the velocity model
iz_mod_sm=-0.016i # Z scale for the velocity model
iz_bp=-0.015i # Z scale for the big BP image
iz_phs=0.006i # Z scale for phase shifts
ih=-0.05i # Moho depth scale
irayp=18i # Ray parameter scale
ivct=0.1i
ivp=1i # Vp scale
x_min_bx=850
x_max_bx=1550
z_min_bx_mod=-25
z_max_bx_mod=100
z_min_bx_bp=20
z_max_bx_bp=60

# Vertical exaggeration
VE=2.0;

# Axes attribute0.0s
B_t_wg=ya5f1
B_t_wv=xa5f1g5

# Coordinates of the phase labels
x_phs1=1
y_phs1=0.5

x_phs2=5
y_phs2=0.5

x_phs3=8.5
y_phs3=0.5

# Plotting ranges
R_mod_bg=0/2200/-50/300
R_mod_sm=${x_min_bx}/${x_max_bx}/${z_min_bx_mod}/${z_max_bx_mod}
R_wg=0/2200/-5/15
R_wv=-5/15/-2.1/2
R_bp_bg=0/2200/0/65
R_bp_sm=${x_min_bx}/${x_max_bx}/${z_min_bx_bp}/${z_max_bx_bp}
R_rayp=0/2200/0.11/0.18
R_vp=0/2200/7.8/8.8
R_phs=0/2200/0/180

# Linewidth
lw_wv=1.5p
lw_st_s=1p # Line width for stations in the zoom-out figures
lw_st_b=2p # Line width for stations in the zoom-out figures
lw_vct=1p
lw_wg=0.5p
lw_wg_s=1p
lw_moho=1.5p
lw_moho_bd=1.5p
lw_vp=1.5p # Line width for Vp
lw_path=0.25p # Line width for ray paths

# Marker sizes
a_st_s=0.1i # Station size for the zoom-out figures
a_st_b=0.2i # Station size for the zoom-in figures
a_wdg=0.5i
a_mk=0.1i
a_vct=0.1i
a_moho=0.075i # Moho marker siz_bges
a_rayp=0.025i # Ray parameter marker size
a_phs=0.025i # Phase shift marker size
a_res=0.025i # Travel time residue marker siz_bges

# Colors for different blocks
c_at=127/127/127
c_lt=0/0/0
c_cr=200/200/200

# Z scale of wiggles
zw=8

# Position of the legends
x_lgd=2i
y_lgd=0.7i
w_lgd=0.7i

# Fonts
ft_fig=12p
ft_ang=10p
ft_ve=10p
ft_phs=6p
ft_md=8p

# Position of figure numbers
x_mod=100
y_mod=325

x_wg=100
y_wg=-3

x_wv=-3
y_wv=0.8

# Position of angle siz_bge
x_ang=350
y_ang=270

# Position of V.E label
x_ve=2000
y_ve=270

# Position of phase label1
x_ss=300
y_ss=0

x_sspmp=250
y_sspmp=7.5

# Position of material properties
x_cr=2000
y_cr=20

x_lt=2000
y_lt=80

x_at=2000
y_at=220

# Spacing between plots

dx7=5i
dy7=8i

dy8=-2i

# Plot the velocity model
gmt psxy ${txt_at} -Jx$ix_bg/$iz_bg -R${R_mod_bg} -BWesN -Bxa500f100g0+l"Distance (km)" -Bya50f25g0+l"Depth (km)" -G${c_at} -X0.75i -Y11i -K -V > $ps
gmt psxy ${txt_lt} -J -R -G${c_lt} -O -K -V >> $ps
gmt psxy ${txt_cr} -J -R -G${c_cr} -O -K -V >> $ps
x_vct=$(awk '{print $1}' $txt_vct)
y_vct=$(awk '{print $2}' $txt_vct)
echo "50 280 ${x_vct} ${y_vct}" | gmt psxy -J -R -Sv${a_vct}+eA+z${ivct} -W${lw_vct},white -O -K -V >> $ps
echo "100 280 ${x_vct} ${y_vct}" | gmt psxy -J -R -Sv${a_vct}+eA+z${ivct} -W${lw_vct},white -O -K -V >> $ps
echo "150 280 ${x_vct} ${y_vct}" | gmt psxy -J -R -Sv${a_vct}+eA+z${ivct} -W${lw_vct},white -O -K -V >> $ps
#echo "$x_ang $y_ang $ang\312" | gmt pstext -J -R -F+f${ft_ang},1,white -O -K -V >> $ps
#awk '{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_s} -Gwhite -W${lw_st_s},127/127/127 -O -K -V >> $ps
awk -v stnm1=$stnm1 -v stnm2=$stnm2 -v stnm3=$stnm3 '$1!=stnm1 && $1!=stnm2 && $1!=stnm3{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_s} -Gwhite -W${lw_st_s},127/127/127 -O -K -V >> $ps
awk -v stnm1=$stnm1 '$1==stnm1{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_s} -Gwhite -W${lw_st_s},green -O -K -V >> $ps
awk -v stnm2=$stnm2 '$1==stnm2{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_s} -Gwhite -W${lw_st_s},cyan -O -K -V >> $ps
awk -v stnm3=$stnm3 '$1==stnm3{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_s} -Gwhite -W${lw_st_s},blue -O -K -V >> $ps
echo "$x_min_bx $z_min_bx_mod" > tmp.txt
echo "$x_max_bx $z_min_bx_mod" >> tmp.txt
echo "$x_max_bx $z_max_bx_mod" >> tmp.txt
echo "$x_min_bx $z_max_bx_mod" >> tmp.txt
gmt psxy tmp.txt -J -R -W1.5p,yellow,- -L -O -K -V >> $ps

#awk -v stnm1=$stnm1 '$1==stnm1{printf "%f 0",$4}' ${3} | gmt psxy -J -R -Sa${a_st_s} -Gwhite -W${lw_st_s},green -O -K -V >> $ps
#awk -v stnm2=$stnm2 '$1==stnm2{printf "%f 0",$4}' ${txt_moho3} | gmt psxy -J -R -Sa${a_st_s} -Gwhite -W${lw_st_s},cyan -O -K -V >> $ps
#awk -v stnm3=$stnm3 '$1==stnm3{printf "%f 0",$4}' ${txt_moho3} | gmt psxy -J -R -Sa${a_st_s} -Gwhite -W${lw_st_s},blue -O -K -V >> $ps

# Plot the data image
# gmt pswiggle ${txt_wg1} -Jx$ix_bg/${it_wg} -R${R_wg} -Z$zw -BWesn -Bxa500f100g0 -B${B_t_wg}+l"Time (s)" -A0 -gy10 -W${lw_wg},black -Y-1.5i -O -K -V >>$ps
# gmt pswiggle ${txt_wg1_s1} -J -R -Z$zw -A0 -gy10 -W${lw_wg_s},green -O -K -V >>$ps
# gmt pswiggle ${txt_wg1_s2} -J -R -Z$zw -A0 -gy10 -W${lw_wg_s},cyan -O -K -V >>$ps
# gmt pswiggle ${txt_wg1_s3} -J -R -Z$zw -A0 -gy10 -W${lw_wg_s},blue -O -K -V >>$ps
#gmt pswiggle ${txt_wg1_s4} -J -R -Z$zw -A0 -gy10 -W${lw_wg_s},cyan,3p_3p:0 -O -K -V >>$ps
#gmt pswiggle ${txt_wg1_s5} -J -R -Z$zw -A0 -gy10 -W${lw_wg_s},cyan,1p_1p:0 -O -K -V >>$ps
# gmt psxy $txt_mvt1 -J -R -W1.5p,red,2p_2p:0 -O -K -V >> $ps
# gmt psxy $txt_mvt2 -J -R -W1.5p,red,2p_2p:0 -O -K -V >> $ps
gmt grdimage ${grd1_dat} -Jx${ix_bg}/${it_wg} -R${R_wg} -BWsne -Bxa500f100g0 -Bya5f2.5g0+l"Time (s)" -C$cpt -Y-1.5i -O -K -V >> $ps
awk '{print $2,$3}' ${txt_res} | gmt psxy -J -R -Sc${a_res} -Gblack -O -K -V >> $ps
gmt psxy ${txt_hp_s} -J -R -W2p,black -O -K -V >> $ps
gmt psxy ${txt_hp_p} -J -R -W2p,black -O -K -V >> $ps

# gmt pswiggle ${txt_wg2} -Jx$ix_bg/${it_wg} -R${R_wg} -Z$zw -BWesn -Bxa500f100g0 -B${B_t_wg}+l"Time (s)" -A0 -gy10 -W${lw_wg},black -Y-1.5i -O -K -V >>$ps
# gmt pswiggle ${txt_wg2_s1} -J -R -Z$zw -A0 -gy10 -W${lw_wg_s},green -O -K -V >>$ps
# gmt pswiggle ${txt_wg2_s2} -J -R -Z$zw -A0 -gy10 -W${lw_wg_s},cyan -O -K -V >>$ps
# gmt pswiggle ${txt_wg2_s3} -J -R -Z$zw -A0 -gy10 -W${lw_wg_s},blue -O -K -V >>$ps
gmt grdimage ${grd2_dat} -Jx${ix_bg}/${it_wg} -R${R_wg} -BWsne -Bxa500f100g0 -Bya5f2.5g0+l"Time (s)" -C$cpt -Y-1.5i -O -K -V >> $ps
#gmt pswiggle ${txt_wg2_s4} -J -R -Z$zw -A0 -gy10 -W${lw_wg_s},cyan,3p_3p:0 -O -K -V >>$ps
#gmt pswiggle ${txt_wg2_s5} -J -R -Z$zw -A0 -gy10 -W${lw_wg_s},cyan,1p_1p:0 -O -K -V >>$ps

# Plot the ray parameters
gmt psxy ${txt_rayp0} -Jx$ix_bg/$irayp -R${R_rayp} -BWesn -Bxa500f100g0 -Bya0.02f0.01g0+l"Ray parameter (s/km)" -W1.5p,red -Y-1.7i -O -K -V >> $ps
gmt psxy ${txt_rayp1} -J -R -W3p,128/128/128 -O -K -V >> $ps
gmt psxy ${txt_rayp2} -J -R -W1.5p,black -O -K -V >> $ps

# Plot the observed and modeled phases
echo "Plotting Phase Shifts"
gmt psxy ${txt_phs0} -Jx$ix_bg/$iz_phs -R${R_phs} -BWsne -Bxa500f100g0+l"Distance (km)" -Bya60f30g0 -W1.5p,red -Y-1.35i -O -K -V >> $ps
gmt psxy ${txt_phs1}  -J -R -W3p,128/128/128 -O -K -V >> $ps
gmt psxy ${txt_phs2} -J -R -W1.5p,black -O -K -V >> $ps


#Plot the derived upper mantle Vp
# echo "Plotting estimated uppermost-mantle Vp"
# gmt psxy ${txt_vp_est} -Jx$ix_bg/$ivp -R${R_vp} -BWsen -Bxa500f100g0+l"Distance (km)" -Bya0.2f0.1g0 -W${lw_vp},black -Y-1.25i -O -K -V >> $ps
# gmt psxy ${txt_vp0} -J -R -BWSen -W1.5p,red -O -K -V >> $ps
# gmt psxy ${txt_vp_hi} -J -R -BWSen -W1.5p,red,- -O -K -V >> $ps
# gmt psxy ${txt_vp_lo} -J -R -BWSen -W1.5p,red,- -O -K -V >> $ps

# Plot the big back-projection image with one source
gmt grdimage ${grd} -Jx$ix_bg/$iz_bp -R${R_bp_bg} -BWsne -Bxa500f100g0 -Bya20f10g0+l"Depth (km)" -C$cpt -Y-1.5i -O -K -V >> $ps
gmt psxy ${txt_moho0} -J -R -W1.5p,128/128/128 -O -K -V >> $ps
echo "$x_min_bx $z_min_bx_bp" > tmp.txt
echo "$x_max_bx $z_min_bx_bp" >> tmp.txt
echo "$x_max_bx $z_max_bx_bp" >> tmp.txt
echo "$x_min_bx $z_max_bx_bp" >> tmp.txt
gmt psxy tmp.txt -J -R -W1.5p,black,- -L -O -K -V >> $ps

#Plot the big back-projection image with two sources from opposite sides
gmt grdimage ${grd_stk} -Jx$ix_bg/$iz_bp -R${R_bp_bg} -BWSne -Bxa500f100g0+l"Distance (km)" -Bya20f10g0+l"Depth (km)" -C$cpt -Y-1.35i -O -K -V >> $ps
gmt psxy ${txt_moho0} -J -R -W1.5p,128/128/128 -O -K -V >> $ps
echo "$x_min_bx $z_min_bx_bp" > tmp.txt
echo "$x_max_bx $z_min_bx_bp" >> tmp.txt
echo "$x_max_bx $z_max_bx_bp" >> tmp.txt
echo "$x_min_bx $z_max_bx_bp" >> tmp.txt
gmt psxy tmp.txt -J -R -W1.5p,black,- -L -O -K -V >> $ps

# Plot the waveforms of the 3 stations which show the interference of SsPmp by other phases
#gmt psxy ${txt_wv2} -Jx${it_wv}/$ia -R${R_wv} -B${B_t_wv}+l"Time (s)" -B${B_a} -BWeSn -W${lw_wv},cyan -Y-2i -O -K -V >> $ps
#gmt psxy ${txt_wv4} -J -R -W${lw_wv},cyan,3p_3p:0 -O -K -V >> $ps
#gmt psxy ${txt_wv5} -J -R -W${lw_wv},cyan,1p_1p:0 -O -K -V >> $ps

# Plot the ray paths of incident S waves
echo "Plotting the ray paths of incident S waves"
gmt psxy ${txt_at} -Jx$ix_sm/$iz_mod_sm -R${R_mod_sm} -BWesN -Bxa500f100g0+l"Distance (km)" -Bya50f25g0+l"Depth (km)" -G${c_at} -X5.5i -Y9i -O -K -V >> $ps
gmt psxy ${txt_lt} -J -R -G${c_lt} -O -K -V >> $ps
gmt psxy ${txt_cr} -J -R -G${c_cr} -O -K -V >> $ps
gmt psxy ${txt_path_s} -J -R -W${lw_path},red -O -K -V >> $ps
awk -v stnm1=$stnm1 -v stnm2=$stnm2 -v stnm3=$stnm3 '$1!=stnm1 && $1!=stnm2 && $1!=stnm3{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_b} -Gwhite -W${lw_st_b},127/127/127 -O -K -V >> $ps
awk -v stnm1=$stnm1 '$1==stnm1{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_b} -Gwhite -W${lw_st_b},green -O -K -V >> $ps
awk -v stnm2=$stnm2 '$1==stnm2{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_b} -Gwhite -W${lw_st_b},cyan -O -K -V >> $ps
awk -v stnm3=$stnm3 '$1==stnm3{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_b} -Gwhite -W${lw_st_b},blue -O -K -V >> $ps

# Plot the ray paths of reflection P waves (first theoretical then apparent)
echo "Plotting the ray paths of reflection P waves"
gmt psxy ${txt_at} -Jx$ix_sm/$iz_mod_sm -R${R_mod_sm} -BWesn -Bxa500f100g0 -Bya50f25g0+l"Depth (km)" -G${c_at} -Y-2.25i -O -K -V >> $ps
gmt psxy ${txt_lt} -J -R -G${c_lt} -O -K -V >> $ps
gmt psxy ${txt_cr} -J -R -G${c_cr} -O -K -V >> $ps
gmt psxy ${txt_path_post_theo} -J -R -W${lw_path},red -O -K -V >> $ps
gmt psxy ${txt_path_pre_theo} -J -R -W${lw_path},yellow -O -K -V >> $ps
awk -v stnm1=$stnm1 -v stnm2=$stnm2 -v stnm3=$stnm3 '$1!=stnm1 && $1!=stnm2 && $1!=stnm3{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_b} -Gwhite -W${lw_st_b},127/127/127 -O -K -V >> $ps
awk -v stnm1=$stnm1 '$1==stnm1{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_b} -Gwhite -W${lw_st_b},green -O -K -V >> $ps
awk -v stnm2=$stnm2 '$1==stnm2{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_b} -Gwhite -W${lw_st_b},cyan -O -K -V >> $ps
awk -v stnm3=$stnm3 '$1==stnm3{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_b} -Gwhite -W${lw_st_b},blue -O -K -V >> $ps

gmt psxy ${txt_at} -Jx$ix_sm/$iz_mod_sm -R${R_mod_sm} -BWesn -Bxa500f100g0 -Bya50f25g0+l"Depth (km)" -G${c_at} -Y-2.25i -O -K -V >> $ps
gmt psxy ${txt_lt} -J -R -G${c_lt} -O -K -V >> $ps
gmt psxy ${txt_cr} -J -R -G${c_cr} -O -K -V >> $ps
gmt psxy ${txt_path_post_app} -J -R -W${lw_path},red -O -K -V >> $ps
gmt psxy ${txt_path_pre_app} -J -R -W${lw_path},yellow -O -K -V >> $ps
awk -v stnm1=$stnm1 -v stnm2=$stnm2 -v stnm3=$stnm3 '$1!=stnm1 && $1!=stnm2 && $1!=stnm3{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_b} -Gwhite -W${lw_st_b},127/127/127 -O -K -V >> $ps
awk -v stnm1=$stnm1 '$1==stnm1{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_b} -Gwhite -W${lw_st_b},green -O -K -V >> $ps
awk -v stnm2=$stnm2 '$1==stnm2{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_b} -Gwhite -W${lw_st_b},cyan -O -K -V >> $ps
awk -v stnm3=$stnm3 '$1==stnm3{print $2,$3}' ${txt_rec} | gmt psxy -J -R -St${a_st_b} -Gwhite -W${lw_st_b},blue -O -K -V >> $ps

# Plot the small back-projection image for one source
echo "Plotting the small back-projection image for one source"
gmt grdimage ${grd} -Jx$ix_sm/$iz_bp_sm -R${R_bp_sm} -BWSne -Bxa500f100g0 -Bya20f10g0+l"Depth (km)" -C$cpt -Y-1.5i -O -K -V >> $ps
gmt psxy ${txt_moho0} -J -R -W1.5p,128/128/128 -O -K -V >> $ps
gmt psxy ${txt_moho1} -J -R -W1.5p,black,4p_4p:0 -O -K -V >> $ps
gmt psxy ${txt_moho2} -J -R -W1.5p,black -O -K -V >> $ps

# Plot the small back-projection image for two sources
echo "Plotting the small back-projection image for two sources from opposite sides"
gmt grdimage ${grd_stk} -Jx$ix_sm/$iz_bp_sm -R${R_bp_sm} -BWSne -Bxa500f100g0+l"Distance (km)" -Bya20f10g0+l"Depth (km)" -C$cpt -Y-1.5i -O -K -V >> $ps
gmt psxy ${txt_moho0} -J -R -W1.5p,128/128/128 -O -K -V >> $ps
gmt psxy ${txt_moho3} -J -R -W1.5p,black -O -K -V >> $ps

# Plot the waveforms of the 3 stations whose Moho depths are compared
echo "Plotting the individual waveforms"
gmt psxy ${txt_wv1} -Jx${it_wv}/$ia -R${R_wv} -Bxa5f1g0+l"Time (s)" -Bya1f1g0+l"Amplitude" -BWeSn -W${lw_wv},green -Y-2i -O -K -V >> $ps
gmt psxy ${txt_wv2} -J -R -W${lw_wv},cyan -O -K -V >> $ps
gmt psxy ${txt_wv3} -J -R -W${lw_wv},blue -O -V >> $ps

# Convert to PDF
gmt psconvert $ps -Tf -D$dir_pdf -V
