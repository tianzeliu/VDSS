gmt gmtset PS_MEDIA  63iX42i
gmt gmtset FONT_ANNOT_PRIMARY 25 # lat/lon, scale bar values, scale bar title
gmt gmtset FONT_LABEL 25 # scale bar label
gmt gmtset FORMAT_GEO_MAP ddd

gmt gmtset MAP_FRAME_WIDTH 2p
gmt gmtset MAP_FRAME_TYPE plain

gmt gmtset COLOR_MODEL rgb

gmt gmtset PS_PAGE_ORIENTATION portrait

# Font info:

fface=0 # 0=Helvetica; 1=Helv bold

# Event name
evt=2017-080-23-10-25

# Input directory
dir_txt=../TXTFile
dir_ps=../PSFile
dir_grd=../GRDFile
dir_cpt=../CPTFile
dir_pdf=../PDFFile

# The station to be plotted
stnm=61128

# The Moho depth and the trade-off curve
txt_moho=$dir_txt/Moho/$evt/MohoDepth_vdss_env.txt

# The station coordinates
txt_stn=$dir_txt/Stations/$evt/StnQual.txt

# GRD files
grd_topo=$dir_grd/NCC_extend_gtopo30.grd
grad_topo=$dir_grd/NCC_extend_gtopo30.grad
grd_vp=$dir_grd/TravelTimeResidual/Vp_30s_cubic_mean.grd

# Boundary file
txt_bdr=$dir_txt/Boundaries/OrdosBoundary_small.txt

# CPT files
cpt_vp=$dir_cpt/TurningVp.cpt

# The output file name
ps=$dir_ps/TurningVp_cubic_mean_${evt}.ps

# Plotting limits
west_op=106
east_op=115
south_op=33
north_op=42

center_x=$(echo "($west_op+$east_op)/2" | bc -l)
center_y=$(echo "($south_op+$north_op)/2" | bc -l)

# Plotting dimensions
lsc_op=1.2i

# Marker dimensions
a_stn=0.5i
a_mdpt=0.4i
a_src=0.5i

# Transverse Mercator projection for the Ordos map
center_x=$(echo "scale=1; (${west_op}+${east_op})/2" | bc -l)
center_y=$(echo "scale=1; (${south_op}+${north_op})/2" | bc -l)

# Topography map draped with turning velocity map
gmt grdview $grd_topo -Qi -R${west_op}/${east_op}/${south_op}/${north_op} -Jt${center_x}/${lsc_op} -Bxa2f1g0 -Bya2f1g0 -BWSne -G$grd_vp -C$cpt_vp -I$grad_topo -Y3i -K -V > $ps
gmt psxy $txt_bdr -J -R -W5p,black,- -O -K -V >> $ps
awk -v stnm=$stnm '$1==stnm{print $2,$3}' $txt_stn | gmt psxy -J -R -St$a_stn -W5p,black -O -K -V >> $ps
awk -v stnm=$stnm '$1==stnm{print $7,$8}' $txt_moho | gmt psxy -J -R -Sc$a_mdpt -W5p,black -O -K -V >> $ps
awk -v stnm=$stnm '$1==stnm{print $11,$12}' $txt_moho | gmt psxy -J -R -Sa$a_src -W5p,black -O -K -V >> $ps

# Plot the colorbar
gmt psscale -Dx10i/3i+w4i/0.3i -C$cpt_vp -Ba1f0.5g0 -O -V >> $ps

# Convert to PDF file
gmt psconvert $ps -D$dir_pdf -Tf -V
