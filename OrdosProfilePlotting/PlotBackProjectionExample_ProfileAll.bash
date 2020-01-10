gmt gmtset PS_MEDIA  15iX30i
gmt gmtset FONT_ANNOT_PRIMARY 15 # lat/lon, scale bar values, scale bar title
gmt gmtset FONT_LABEL 15 # scale bar label
gmt gmtset FORMAT_GEO_MAP ddd

gmt gmtset MAP_FRAME_WIDTH 0.6i
gmt gmtset MAP_FRAME_TYPE plain

gmt gmtset COLOR_MODEL rgb

gmt gmtset PS_PAGE_ORIENTATION portrait

# Font info
fface=0 # 0=Helvetica; 1=Helv bold

# Event name
evt=2017-080-23-10-25

# Line name
lin=Line1

# Input directory
dir_txt=../TXTFile/${evt}/${lin}
dir_ps=../PSFile
dir_grd=../GRDFile/${evt}/${lin}
dir_cpt=../CPTFile
dir_pdf=../PDFFile

# The output file name
ps=$dir_ps/BPexample_2017-080-23-10-25_profile_raw_phs0_int.pdf

# The plot projections
J_trk=x0.015i/0.0012i
J_res=x0.015i/0.5i
J_rayp=x0.015i/55i
J_wg=x0.015i/-0.1i
J_bp=x0.015i/-0.045i
J_phs=x0.015i/0.015i

# The plot ranges
R_trk=0/550/0/1500
R_res=0/550/-2.5/2
R_rayp=0/550/0.11/0.15
R_wg=0/550/-10/20
R_bp=0/550/0/70
R_phs=0/550/-10/190

# The dimension for the wiggles
Z_wg=0.8

# The CPT file for the back-projection psimage
cpt_bp=$dir_cpt/BackProjPol.cpt

# The TXT files
txt_trk=$dir_txt/TopoTrack.txt
txt_res=$dir_txt/SACinfo_p_raw.txt
txt_rayp0=$dir_txt/Rayp0_prj.txt
txt_rayp=$dir_txt/Rayp_prj.txt
txt_fit=$dir_txt/ResInt.txt
txt1_phs=$dir_txt/PhaseShift_good.txt
txt2_phs=$dir_txt/PhaseShift_bad.txt
txt3_phs=$dir_txt/PhaseShift_interp.txt
txt_wg_p=$dir_txt/Wiggle_p_raw.txt
txt_wg_s=$dir_txt/Wiggle_s_raw.txt
txt0_moho=$dir_txt/MohoDepth0_phs0_int_est_6.2.txt
txt1_moho=$dir_txt/MohoDepth_phs0_int_est_6.2.txt
txt2_moho=$dir_txt/MohoDepth_phs0_int_est_6.0.txt
txt_moho_rf=$dir_txt/MohoDepthRF_FengDong.txt

# The GRD files
grd0_img=$dir_grd/BackProjectImg0_phs0_int_gauss_6.2.grd
grd1_img=$dir_grd/BackProjectImg_phs0_int_gauss_6.2.grd
grd2_img=$dir_grd/BackProjectImg_phs0_int_gauss_6.0.grd

# Plot the topography cross-section
gmt psxy $txt_trk -J${J_trk} -R${R_trk} -W2p,black -Bxa100f100+l"Distance (km)" -Bya500f500+l"Elevation (m)" -BWsNe -i2,3 -Y26i -K -V > $ps

# Plot the travel time residues
gmt psxy $txt_res -J${J_res} -R${R_res} -Bxa100f100 -Bya1f1+l"Traveltime residual (s)" -BWnse -Sc0.2i -Gcyan -W1p,black -i4,5 -Y-3i -O -K -V >> $ps
gmt psxy $txt_fit -J -R -W3p,red -O -K -V >> $ps

# Plot the ray parameters
gmt psxy $txt_rayp0 -J${J_rayp} -R${R_rayp} -Bxa100f100 -Bya0.01f0.01+l"Ray parameter (s/km)" -BWesn -W3p,black -Y-3i -O -K -V >> $ps
gmt psxy $txt_rayp -J${J_rayp} -R${R_rayp} -W3p,red -O -K -V >> $ps

# Plot the S wiggles
gmt pswiggle ${txt_wg_s} -J${J_wg} -R${R_wg} -Z${Z_wg} -Bxa100f100+l"Distance (km)" -Bya5f5+l"Time (s)" -BWesn -A0 -G+red -G-gray -gy10 -Y-3.75i -O -K -V >> $ps

# Plot the P wiggles
gmt pswiggle ${txt_wg_p} -J${J_wg} -R${R_wg} -Z${Z_wg} -Bxa100f100+l"Distance (km)" -Bya5f5+l"Time (s)" -BWesn -A0 -G+blue -G-gray -gy10 -Y-3.75i -O -K -V >> $ps

# Plot the measured phase shifts
gmt psxy $txt1_phs -J${J_phs} -R${R_phs} -Bxa100f100 -Bya30f30 -BWesn -Sd0.2i -Gcyan -W1p,black -O -K -V -Y-3.75i >> $ps
gmt psxy $txt2_phs -J -R -Sd0.2i -Ggray -W1p,black -O -K -V >> $ps
gmt psxy $txt3_phs -J -R -W3p,red -O -K -V >> $ps

# Plot the original backprojection image
gmt grdimage ${grd0_img} -J${J_bp} -R${R_bp} -C$cpt_bp -Bxa100f100+l"Distance (km)" -Bya10f10+l"Depth (km)" -BWesn -Y-3.75i -O -K -V >> $ps
gmt psxy ${txt0_moho} -J${J_bp} -R${R_bp} -Bxa100f100+l"Distance (km)" -Bya10f10+l"Depth (km)" -Bwesn -W3p,black,4p_4p:0 -O -K -V >> $ps
gmt psxy ${txt_moho_rf} -J -R -Sc0.15i -W1p,black -Gyellow -O -K -V >> $ps

# Plot the corrected backprojection image
gmt grdimage ${grd1_img} -J${J_bp} -R${R_bp} -C$cpt_bp -Bxa100f100+l"Distance (km)" -Bya10f10+l"Depth (km)" -BWeSn -Y-3.75i -O -K -V >> $ps
gmt psxy ${txt0_moho} -J -R -W3p,black,4p_4p:0 -O -K -V >> $ps
gmt psxy ${txt1_moho} -J -R -W3p,black -O -K -V >> $ps
gmt psxy ${txt_moho_rf} -J -R -Sc0.15i -W1p,black -Gyellow -O -V >> $ps

# Convert to PDF file
gmt psconvert $ps -D$dir_pdf -Tf -V
