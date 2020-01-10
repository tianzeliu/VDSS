# Plot the event distribution, waveforms, and tpmp and phase-shift as functions of ray parameters
#!/bin/bash

gmt gmtset PS_MEDIA 10iX11i
gmt gmtset FONT_ANNOT_PRIMARY 10p,0,black # lat/lon, scale bar values, scale bar title
gmt gmtset FONT_LABEL 8p # scale bar label
gmt gmtset PS_PAGE_ORIENTATION portrait

#gmt gmtset PS_PAGE_ORIENTATION portrait
gmt gmtset MAP_FRAME_PEN 1p,black
gmt gmtset MAP_FRAME_TYPE plain
gmt gmtset MAP_GRID_PEN_PRIMARY 0.01i,gray,0.05_0.05:0i
gmt gmtset FORMAT_GEO_MAP D
gmt gmtset MAP_ANNOT_OFFSET_PRIMARY 2p

# Input and output
dir_ps=../PSFile
dir_txt=../TXTFile
dir_pdf=../PDFFile
dir_grd=../GRDFile
dir_cpt=../CPTFile

ps=${dir_ps}/RockVpKappa_Tcorr250.ps

txt_rck_hacker=${dir_txt}/CrustRockVpk_Hacker_ValueOnly_Tcorr250.txt
txt_rck_chris=${dir_txt}/CrustRockVpk_Christensen_Tcorr250.txt

grd_pdf=${dir_grd}/VpKappa_rand5000_pdf.grd

# CPT file for H-k HkImage
cpt_pdf=${dir_cpt}/VpKappaPdf.cpt

# CPT file for rock silica content
cpt_rck=${dir_cpt}/Silica.cpt

# Line widths
lw_err=1.5p

# Plot the Vp-k PDF and values for rock samples
gmt grdimage ${grd_pdf} -Jx1.75i/10i -R5.8/7.5/1.62/1.92 -Bxa0.2f0.1g0 -Bya0.04f0.02g0 -BWesN -C$cpt_pdf -K -V > $ps
gmt psxy ${txt_rck_hacker} -J -R -Sc0.075i -C$cpt_rck -W1p,black -O -K -V >> $ps
gmt psxy ${txt_rck_chris} -J -R -Sd0.15i -Ggray -Exy+p$lw_err,gray -i1,2,3,4 -O -K -V >> $ps
gmt psscale -Dx3.25i/0i+w3i/0.1i -C$cpt_pdf -Ba20f10g0 -O -K -V >> $ps
gmt psscale -Dx0i/-0.25i+w3i/0.1i+h -C$cpt_rck -Ba10f10g0 -O -V >> $ps

# Convert to PDF file
gmt psconvert $ps -Tf -D$dir_pdf -V
