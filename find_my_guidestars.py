#!/usr/bin/env python
'''Script to determine guide star options in a
given 2dF field of view.

17 December 2018 - Matt Taylor - alpha development
16 January  2019 - Matt Taylor - added UCAC for guide star catalogue
'''

#-----------------------------------------------------------------------
import numpy as np
from astropy.io import fits, ascii
from astropy.stats import sigma_clip
from astropy.coordinates import Angle, SkyCoord, match_coordinates_sky
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.visualization import (MinMaxInterval, LogStretch,
                                   ImageNormalize)
from matplotlib.patches import Circle
from astroquery.vizier import Vizier
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
import sys, os

import warnings
warnings.filterwarnings("ignore")
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------

version = '17 December 2018'
# coord_2dF   = SkyCoord(float(sys.argv[1]), float(sys.argv[2]), unit=(u.degree,u.degree))
if ":" in sys.argv[1]:
	coord_2dF   = SkyCoord(sys.argv[1], sys.argv[2], unit=(u.hourangle,u.degree))
else:
	coord_2dF   = SkyCoord(sys.argv[1], sys.argv[2], unit=(u.degree,u.degree))

def usage():
    print ''
    print 'NAME'
    print '       find_my_guidestars.py - Find guide star options for a given 2dF pointing\n'
    print 'SYNOPSIS'
    print '       find_my_guidestars.py RA Dec -> RA and Dec in either hh:mm:ss dd:mm:ss or decimal degree formats\n'
    print 'DESCRIPTION'
    print '       Find a list of guide star options from 2MASS based on SCABS catalogues.'
    print ' '
    print 'VERSION'
    print '       ', version
    print ''
    raise SystemExit

#-----------------------------------------------------------------------

def main():
	gs_coo_2mass, gs_vmag = get_gs_cat(coord_2dF) #find catalogue point sources within the 2dF field centred at the requested coordinates
	gv_gs_coords, gv_gs_mags = scabs_xmatch(gs_coo_2mass, gs_vmag) #match 2MASS sources to SCABS g-band catalogue
	save_results(gv_gs_coords, gv_gs_mags)
	plot_results(coord_2dF, gv_gs_coords, gv_gs_mags)
	make_stamps(gv_gs_coords)
	plt.show()


#----------Query 2MASS and Tycho-2 catalogues for stars in the given field-------------
def get_gs_cat(coo):
	Vizier.ROW_LIMIT = -1

	# print "querying gaia..."
	# result = Vizier.query_region(coo, radius=Angle(1.,"deg"), catalog="GSC2.3")
	# data = result["I/305/out"]
	# mask = np.logical_and(data["Vmag"] >= 10, data["Vmag"] <= 15)
	# data = data[mask]
	# print data
	# exit()
	
	#UCAC
	viz = Vizier(columns=["RAJ2000", "DEJ2000", "Vmag"])
	viz.ROW_LIMIT = -1
	result_ucac = viz.query_region(coo, radius=Angle(1.,"deg"), catalog="UCAC4")
	data_ucac = result_ucac['I/322A/out']

	mask       = np.logical_and(data_ucac["Vmag"] >= 10., data_ucac["Vmag"] <= 18.)
	data_ucac   = data_ucac[mask]
	gv_ucac_coo = SkyCoord(data_ucac["RAJ2000"], data_ucac["DEJ2000"], unit=(u.degree,u.degree))
	gv_ucac_mag = data_ucac["Vmag"]

	#2MASS point source coordinates
	viz = Vizier(columns=["RAJ2000", "DEJ2000", "Bmag", "Rmag"])
	viz.ROW_LIMIT = -1
	result_2mass = viz.query_region(coo, radius=Angle(1.,"deg"), catalog="2MASS-PSC")
	gv_2mass_coo = SkyCoord(result_2mass["II/246/out"]["RAJ2000"], result_2mass["II/246/out"]["DEJ2000"], unit=(u.degree,u.degree))

	#Tycho-2 coordinates and V-band mags	
	# result_tycho = Vizier.query_region(coo, radius=Angle(1.,"deg"), catalog="Tycho-2")
	# gv_tycho_coo = SkyCoord(result_tycho["I/259/tyc2"]["RA_ICRS_"], result_tycho["I/259/tyc2"]["DE_ICRS_"], unit=(u.degree,u.degree))
	# gv_tycho_mag = result_tycho["I/259/tyc2"]["VTmag"]

	#GSC2.3 coordinates and V mags
	# result_gsc = Vizier.query_region(coo, radius=Angle(1.,"deg"), catalog="GSC2.3")
	# data_gsc   = result_gsc["I/305/out"]
	# mask       = np.logical_and(data_gsc["Vmag"] >= 10., data_gsc["Vmag"] <= 25.)
	# data_gsc   = data_gsc[mask]
	# gv_gsc_coo = SkyCoord(data_gsc["RAJ2000"], data_gsc["DEJ2000"], unit=(u.degree,u.degree))
	# gv_gsc_mag = data_gsc["Vmag"]


	# gv_gs_inds, d2d, d3d   = match_coordinates_sky(gv_tycho_coo,gv_2mass_coo)
	# gv_gs_inds, d2d, d3d   = match_coordinates_sky(gv_gsc_coo,gv_2mass_coo)
	gv_gs_inds, d2d, d3d   = match_coordinates_sky(gv_ucac_coo,gv_2mass_coo)

	gv_2mass_coo = SkyCoord([gv_2mass_coo[ii].ra.deg for ii in gv_gs_inds], [gv_2mass_coo[ii].dec.deg for ii in gv_gs_inds], unit=(u.degree,u.degree))

	return gv_2mass_coo, gv_ucac_mag# gv_gsc_mag # gv_tycho_mag

#----------Cross-match to SCABS catalogues for sources within 0.1"-------------
def scabs_xmatch(gs_star_coo,gs_star_mag):
	#Load SCABS g-band source catalogue
	prefix = "/Users/mtaylor/Projects/SCABS/manuscripts/SCABS1-data/data"
	# prefix = "/Users/mtaylor/Dropbox/Projects/CenA_2dF_GCs/SCABS_Tile1-7_rband"
	scabs_data = fits.open(prefix+"/SCABS_Tiles1-7_sources_r-band.fits")[1].data
	scabs_coo = SkyCoord(scabs_data["ALPHA_J2000"], scabs_data["DELTA_J2000"], unit=(u.degree,u.degree))

	print "Begin matching..." #find a match for each catalogue point source
	gv_scabs_ind, gv_scabs_d2d, gv_scabs_d3d = match_coordinates_sky(gs_star_coo,scabs_coo)

	gv_sep = 0.3 #separation limit in arcsec
	mask = (gv_scabs_d2d.to(u.arcsec).value <= gv_sep)
	gv_scabs_ind = gv_scabs_ind[mask] ; gv_scabs_d2d = gv_scabs_d2d[mask] #keep only good matches
	print 'Limiting to matches within %.1f"....found %i matches.' % (gv_sep, len(gv_scabs_ind))
	gv_scabs_coo = SkyCoord([scabs_data["ALPHA_J2000"][ii] for ii in gv_scabs_ind], [scabs_data["DELTA_J2000"][ii] for ii in gv_scabs_ind], unit=(u.degree,u.degree))

	gv_inds, d2d, d3d = match_coordinates_sky(gv_scabs_coo,gs_star_coo)
	gv_star_mag = np.array([gs_star_mag[ii] for ii in gv_inds])

	#Limit sources to within given magnitude range
	mag_faint_lim  = 13.5 #faint limit in magnitude
	mag_bright_lim = 12.5 #bright limit

	mask = np.logical_and(gv_star_mag >= mag_bright_lim, gv_star_mag <= mag_faint_lim)

	gv_scabs_coo = gv_scabs_coo[mask]
	gv_star_mag = gv_star_mag[mask]
	print "Limiting to guide stars in %.1f <= V <= %.1f...found %i matches." % (mag_bright_lim, mag_faint_lim, len(gv_star_mag))

	# plt.figure()
	# plt.plot(gs_star_coo.ra.deg, gs_star_coo.dec.deg,'bo',alpha=0.5)
	# plt.plot(gv_scabs_coo.ra.deg, gv_scabs_coo.dec.deg,'ro',alpha=0.5)
	# plt.show()
	# exit()

	return gv_scabs_coo, gv_star_mag

#----------Save guidestar coordinates to file-------------
def save_results(gs_coo,gs_mag):
	fileout = open("2dF_guidestars_%s%s.txt" % (sys.argv[1], sys.argv[2]),"w")
	for ii in range(len(gs_coo)):
		print >> fileout, "%03d	%f 	%f 	%.2f" % (ii+1, gs_coo.ra.deg[ii], gs_coo.dec.deg[ii], gs_mag[ii])
	fileout.close()
	# exit()


#----------Plot guidestar options on 2dF field of view-------------
def plot_results(fov_coords, gs_coo, gs_mag):
	#Load the dss background image for wcs
	dss_hdu = fits.open("dss_image.fits")
	dss_im  = dss_hdu[0].data

	dss_wcs = WCS(dss_hdu[0].header)

	fig = plt.figure(figsize=(10,10))
	gs = gridspec.GridSpec(2,2, height_ratios=[0.03,1], width_ratios=[1,0.03])
	gs.update(left=0.05, right=0.95, bottom=0.08, top=0.93, wspace=0.02, hspace=0.03)

	ax = fig.add_subplot(gs[1,0], projection=dss_wcs)

	fov = SphericalCircle((fov_coords.ra.deg,fov_coords.dec.deg)*u.degree,1.*u.degree, edgecolor="None", facecolor="k",transform=ax.get_transform('fk5'), linestyle='solid', lw=2, alpha=0.1,zorder=1)
	ax.add_patch(fov)

	# ax.plot(gs_coo.ra.deg,gs_coo.dec.deg,'o',transform=ax.get_transform("fk5"))
	ax.plot([fov_coords.ra.deg], [fov_coords.dec.deg],'kX',ms=15,transform=ax.get_transform("fk5"))
	gs_ax = ax.scatter(gs_coo.ra.deg,gs_coo.dec.deg,c=gs_mag,
				marker='o', vmin=np.min(gs_mag), vmax=np.max(gs_mag),
				edgecolor="None",cmap="viridis",s=75,transform=ax.get_transform("fk5"))
	ax.grid()
	ax.set_xlabel(r"RA (J2000)",fontsize=18)
	ax.set_ylabel(r"Dec (J2000)",fontsize=18)

	#Label each GS with it's corresponding number to match them with postage stamps
	for ii in range(len(gs_coo)):
		ax.text(gs_coo[ii].ra.deg,gs_coo[ii].dec.deg, ii+1, color="k", fontsize=10, fontdict={"weight": "bold"}, transform=ax.get_transform('fk5'))#horizontalalignment="center", verticalalignment="center", transform=ax.get_transform('fk5'))

	cbax = fig.add_subplot(gs[1,1])
	cb = Colorbar(ax=cbax, mappable=gs_ax, orientation='vertical', ticklocation='right')
	cb.set_label(r"Guide Star $m_V$ (mag)",fontsize=18)#,labelpad=20)

	#Load GC catalogue
	gc_data = ascii.read("gcs_Tile1-7_coords_ugriz.txt")
	gc_coo  = SkyCoord(gc_data["RA"], gc_data["Dec"], unit=(u.hourangle,u.degree))
	#isolate GCs falling within 2dF fov
	seps    = fov_coords.separation(gc_coo)
	mask    = (seps.to(u.degree).value < 1.)
	gc_data = gc_data[mask] ; gc_coo = gc_coo[mask]
	#Calculate V-band mags and limit to 17 < V < 21.5
	gc_vmag = gc_data["g_mag"] - 0.58*(gc_data["g_mag"]-gc_data["r_mag"]) - 0.01
	mask = np.logical_and(gc_vmag >= 17., gc_vmag <= 21.5)
	gc_vmag = gc_vmag[mask] ; gc_data = gc_data[mask] ; gc_coo = gc_coo[mask]
	print "Number of GC targets available in the field = %i" % (len(gc_vmag))
	#Save the list of GC targets
	fileout = open("2dF_GCs_%s%s.txt" % (sys.argv[1], sys.argv[2]),"w")
	for ii in range(len(gc_coo)):
		print >> fileout, "%03d	%f 	%f 	%.2f" % (ii+1, gc_coo.ra.deg[ii], gc_coo.dec.deg[ii], gc_vmag[ii])
	fileout.close()

	#plot 'em
	gc_ax = ax.scatter(gc_coo.ra.deg,gc_coo.dec.deg,c=gc_vmag,
				marker='o', vmin=np.min(gc_vmag), vmax=np.max(gc_vmag),
				edgecolor="None",cmap="copper",s=20,transform=ax.get_transform("fk5"),alpha=0.75)

	cbax = fig.add_subplot(gs[0,0])
	cb = Colorbar(ax=cbax, mappable=gc_ax, orientation='horizontal', ticklocation='top')
	cb.set_label(r"GC Candidate $m_V$ (mag)",fontsize=18)#,labelpad=20)


	plt.savefig("2dF_guidestar_search_%s%s.pdf" % (sys.argv[1], sys.argv[2]),bbox_inches="tight",overwrite=True)
	# plt.show()

#----------Create r'-band postage stamp images for visual inspection-------------
def make_stamps(gs_coo):
	#Directory prefix for images
	prefix = "/Users/mtaylor/Projects/SCABS/stacks"
	# prefix = "/Users/mtaylor/Dropbox/Projects/CenA_2dF_GCs/SCABS_Tile1-7_rband"
	# hdu = fits.open(prefix+"/survey_tile1_r_short_ALIGNi.fits")
	# wcs = WCS(hdu[0].header)
	# data = hdu[0].data

	#Check which tile the cutout should come from
	tiles_coo = {"1": SkyCoord(201.365062792, -43.0191125833, unit=(u.degree,u.degree)),  #central tile coordinates
				 "2": SkyCoord(203.167159865, -41.8681243611, unit=(u.degree,u.degree)),
				 "3": SkyCoord(200.944693054, -41.2104168056, unit=(u.degree,u.degree)),
				 "4": SkyCoord(199.128469057, -42.3614050278, unit=(u.degree,u.degree)),
				 "5": SkyCoord(199.556484167, -44.1701008056, unit=(u.degree,u.degree)),
				 "6": SkyCoord(201.803072602, -44.8278083611, unit=(u.degree,u.degree)),
				 "7": SkyCoord(203.588227236, -43.6768201389, unit=(u.degree,u.degree))}

	#Check for the existence of the cutout directory and create if not there
	if not os.path.exists("cutouts_%s%s" % (sys.argv[1], sys.argv[2])): os.makedirs("cutouts_%s%s" % (sys.argv[1], sys.argv[2]))

	#Make cutouts
	size = u.Quantity((120.,150.), u.arcsec)  #Angular size of the cutout
	for ii in range(len(gs_coo)):
		#for each guide star, check tiles that include the guide star (e.g. is the GS coordinates within +/- 1 degree of tile centre?)
		temp_tiles = []
		for jj in tiles_coo:
			temp_coo = tiles_coo[str(jj)]
			if temp_coo.separation(gs_coo[ii]).to(u.degree).value <= 1.0:
				temp_tiles.append(str(jj))
		# print temp_tiles
		for tile_num in temp_tiles:
			#Open r-band image for the given tile
			# print tile_num
			hdu = fits.open(prefix+"/survey_tile%s_r_short_ALIGNi.fits" % tile_num)
			wcs = WCS(hdu[0].header)
			data = hdu[0].data

			cutout = Cutout2D(data, gs_coo[ii], size, wcs=wcs)

			norm = ImageNormalize(cutout.data, stretch=LogStretch())
			fig = plt.figure()
			ax1 = plt.subplot(projection=cutout.wcs)
			ax1.imshow(cutout.data,origin='lower',norm=norm,cmap='gray_r',vmin=-5,vmax=750,aspect='auto')
			ax1.grid(color='C0', ls='dashed', lw=1.5)

			#Add circle to indicate 2dF fibre extent
			fibre = SphericalCircle((gs_coo[ii].ra.deg,gs_coo[ii].dec.deg)*u.degree,1.*u.arcsec, edgecolor="green", facecolor="None",transform=ax1.get_transform('fk5'), linestyle='solid', lw=2, alpha=0.9,zorder=1)
			ax1.add_patch(fibre)

			lon = ax1.coords[0]
			lat = ax1.coords[1]
			lon.set_axislabel(r"$\alpha$ (J2000)",fontsize=22)
			lon.set_major_formatter("dd:mm")
			lon.set_ticks(size=6,width=2)
			lon.set_ticklabel(size=14)
			lon.set_ticks_position('b')

			lat.set_axislabel(r"$\delta$ (J2000)",fontsize=22)
			lat.set_major_formatter("dd:mm")
			lat.set_ticks(size=6,width=2)
			lat.set_ticklabel(size=14)
			lat.set_ticks_position('lr')
			plt.savefig("cutouts_%s%s/gs_%03d_tile%s.pdf" % (sys.argv[1], sys.argv[2], ii+1, tile_num),bbx_inches="tight",overwrite=True)
			hdu.close()
		plt.close("all")

if __name__ == "__main__":
    main()