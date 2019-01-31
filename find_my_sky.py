#!/usr/bin/env python
'''Script to find appropriate sky fibre placement 
options in a given 2dF field of view using Source
Extractor.

16 January 2018 - Matt Taylor - alpha development'''

#-----------------------------------------------------------------------
import numpy as np
from numpy import random
from astropy.io import fits, ascii
# from astropy.stats import sigma_clip
from astropy.coordinates import Angle, SkyCoord, match_coordinates_sky
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.visualization import (MinMaxInterval, LogStretch,
                                   ImageNormalize)
from regions import CircleSkyRegion
from matplotlib.patches import Circle
from astroquery.vizier import Vizier
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
import sys, os
import subprocess

import warnings
warnings.filterwarnings("ignore")
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------

version = '16 January 2019'
# coord_2dF   = SkyCoord(float(sys.argv[1]), float(sys.argv[2]), unit=(u.degree,u.degree))
if ":" in sys.argv[1]:
	coord_2dF   = SkyCoord(sys.argv[1], sys.argv[2], unit=(u.hourangle,u.degree))
else:
	coord_2dF   = SkyCoord(sys.argv[1], sys.argv[2], unit=(u.degree,u.degree))

# prefix = "/Users/mtaylor/Projects"
prefix = "/Users/mtaylor/Dropbox/Projects/CenA_2dF_GCs/SCABS_Tiles1-7_iband"


def usage():
    print ''
    print 'NAME'
    print '       find_my_sky.py - Find sky fibre placement options for a given 2dF pointing\n'
    print 'SYNOPSIS'
    print '       find_my_sky.py RA Dec -> RA and Dec in either hh:mm:ss dd:mm:ss or decimal degree formats\n'
    print 'DESCRIPTION'
    print '       Find a list of sky fibre placement options using prepared Source Extractor segmentation maps.'
    print ' '
    print 'VERSION'
    print '       ', version
    print ''
    raise SystemExit

#-----------------------------------------------------------------------
def main():
    #Check for the existence of the cutout directory and create if not there
    if not os.path.exists("sky_cutouts_%s%s" % (sys.argv[1], sys.argv[2])): os.makedirs("sky_cutouts_%s%s" % (sys.argv[1], sys.argv[2]))

    # make_segmentation(prefix+"/survey_tile7_i_short.fits")
    # exit()
    find_sky(coord_2dF)

#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
'''If desired, run Source Extractor on an image to create a new
segmentation map.'''
def make_segmentation(image):
    '''Run Source Extractor on the input image (using) the associated weight map
    to create a segmentation map to inform where sources are, and more importantly are not.'''
    se_command = "/Users/mtaylor/local/bin/./sex" #local call to source extractor
    command    = [se_command, image, "-c", "SE_config/ctio_decam.sex",
                "-PARAMETERS_NAME","SE_config/ctio_decam.param",
                "-FILTER_NAME","SE_config/gauss_3.5_7x7.conv",
                "-STARNNW_NAME","SE_config/default.nnw",
                "-CATALOG_NAME","deleteme.ldac",
                "-WEIGHT_IMAGE",image.replace(".fits",".WEIGHT.fits"),
                "-MAG_ZEROPOINT","30.",
                "-CHECKIMAGE_TYPE","SEGMENTATION",
                "-CHECKIMAGE_NAME",image.replace(".fits",".seg.fits")]
    print command
    subprocess.call(command)
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
'''Search the 2dF field of view centred on the given coordinates on a
segmentation map to find regions where there are no sources.'''
def find_sky(coo):
    '''Given a coordinate, generate random set of 50 coordinates uniformly spread across the 2dF field.
    For each coordinate, check against the segmentation image to see if a 3" radius surrounding it is empty
    of sources. If so, keep it. If not, [think about method of selecting new coordinate nearby.]'''
    #Directory prefix for images
    # prefix = "/Users/mtaylor/Projects/SCABS/stacks"
    prefix = "/Users/mtaylor/Dropbox/Projects/CenA_2dF_GCs/SCABS_Tiles1-7_iband"

    #Check which tile the cutout should come from
    tiles_coo = {"1": SkyCoord(201.365062792, -43.0191125833, unit=(u.degree,u.degree)),  #central tile coordinates
                 "2": SkyCoord(203.167159865, -41.8681243611, unit=(u.degree,u.degree)),
                 "3": SkyCoord(200.944693054, -41.2104168056, unit=(u.degree,u.degree)),
                 "4": SkyCoord(199.128469057, -42.3614050278, unit=(u.degree,u.degree)),
                 "5": SkyCoord(199.556484167, -44.1701008056, unit=(u.degree,u.degree)),
                 "6": SkyCoord(201.803072602, -44.8278083611, unit=(u.degree,u.degree)),
                 "7": SkyCoord(203.588227236, -43.6768201389, unit=(u.degree,u.degree))}

    #Open output file for writing
    fileout = open("sky_cutouts_%s%s/sky_fibres_%s%s.txt" % (sys.argv[1], sys.argv[2], sys.argv[1], sys.argv[2]),"w")
    print >> fileout, "#NUMBER  tile    RA              Dec"

    #Generate random set of 50 coordinates within the 2dF field of view
    n_coords = 50
    chk = 0
    while chk < n_coords:
        temp_ra = random.uniform(coo.ra.deg-1.,coo.ra.deg+1.,1) ; temp_dec = random.uniform(coo.dec.deg-1.,coo.dec.deg+1.,1)
        print temp_ra, temp_dec
        temp_sky_coo = SkyCoord(temp_ra[0], temp_dec[0], unit=(u.degree,u.degree))

        #Given this point, check the tile that it falls within
        for jj in tiles_coo:
                #Get image WCS
            hdu = fits.open(prefix+"/survey_tile%s_i_short.fits" % jj)
            wcs = WCS(hdu[0].header)
            hdu.close()

            temp_coo = tiles_coo[str(jj)]
            if temp_coo.separation(temp_sky_coo).to(u.degree).value <= 1.0:
                segmentation_map = prefix+"/survey_tile%s_i_short.seg.fits" % jj

                #Go find out if there is a source near the sky fibre coordinate. Shift around until a clear spot is found.

                sky_coo = find_clear_patch(temp_sky_coo,segmentation_map)
                #Save the segmentation map cutout
                plt.savefig("sky_cutouts_%s%s/sky%02d_tile%s_seg.pdf" % (sys.argv[1], sys.argv[2], chk+1, jj), bbox_inches="tight")
                plt.close()

                #Save the image cutout
                hdu = fits.open(prefix+"/survey_tile%s_i_short.fits" % jj)
                wcs = WCS(hdu[0].header)
                data = hdu[0].data
                size = u.Quantity((10.,10.), u.arcsec)

                cutout = Cutout2D(data, sky_coo, size, wcs=wcs)
                hdu.close()

                plt.figure(figsize=(10,10))
                norm = ImageNormalize(cutout.data, stretch=LogStretch())
                ax1 = plt.subplot(projection=cutout.wcs)
                ax1.imshow(cutout.data,origin='lower',norm=norm,cmap='gray_r',vmin=-5,vmax=750,aspect='auto')

                fibre = SphericalCircle((sky_coo.ra.deg,sky_coo.dec.deg)*u.degree,1.*u.arcsec, edgecolor="green", facecolor="None",transform=ax1.get_transform('fk5'), linestyle='solid', lw=2, alpha=0.9,zorder=1)
                ax1.add_patch(fibre)

                ax1.grid(color='C0', ls='dashed', lw=1.5)
                ax1.set_xlabel(r"$\alpha$ (J2000)")
                ax1.set_ylabel(r"$\delta$ (J2000)")
                plt.savefig("sky_cutouts_%s%s/sky%02d_tile%s_image.pdf" % (sys.argv[1], sys.argv[2], chk+1, jj), bbox_inches="tight")
                plt.close()

                print >> fileout, "S%02d         %s      %f  %f" % (chk+1, jj, sky_coo.ra.deg, sky_coo.dec.deg)
                temp_sky_coo = sky_coo

        chk += 1
    fileout.close()
    exit()
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
def find_clear_patch(temp_sky_coo,seg_map):
    #we found a tile that the point falls within, check the segmentation map
    #cut out a 2"x2" section of the segmentation map
    hdu = fits.open(seg_map)
    wcs = WCS(hdu[0].header)
    data = hdu[0].data
    size = u.Quantity((2.,2.), u.arcsec)

    #Cut out a 2"x2" section of segmentation map, check that cutout is free of sources.
    #If a source is included, shift randomly by within +/- 1" and try again until a clear patch is found
    clear_chk = 0
    while clear_chk == 0:

        cutout_fibre = Cutout2D(data, temp_sky_coo, size, wcs=wcs)

        #how many pixels are non-zero? (i.e. how many pixels are associated with a source?)
        n_sources = np.count_nonzero(cutout_fibre.data)
        print clear_chk, n_sources
        #If there is a source in the cutout, shift randomly by within +/- 1" and try again until a clear patch is found.
        if n_sources > 0:
            temp_ra = temp_sky_coo.ra.deg + (2.*random.random()-1.)/3600. ; temp_dec = temp_sky_coo.dec.deg + (2.*random.random()-1.)/3600.
            temp_sky_coo = SkyCoord(temp_ra, temp_dec, unit=(u.degree,u.degree))
        else:
            clear_chk += 1

    cutout_seg = Cutout2D(data, temp_sky_coo, u.Quantity((10.,10.), u.arcsec), wcs=wcs) 

    plt.figure(figsize=(10,10))
    ax1 = plt.subplot(projection=cutout_seg.wcs)
    ax1.imshow(cutout_seg.data,origin='lower',cmap='gray',vmin=-5,vmax=750,aspect='auto')

    fibre = SphericalCircle((temp_sky_coo.ra.deg,temp_sky_coo.dec.deg)*u.degree,1.*u.arcsec, edgecolor="green", facecolor="None",transform=ax1.get_transform('fk5'), linestyle='solid', lw=2, alpha=0.9,zorder=1)
    ax1.add_patch(fibre)

    ax1.set_xlabel(r"$\alpha$ (J2000)")
    ax1.set_ylabel(r"$\delta$ (J2000)")

    hdu.close()
    return temp_sky_coo
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
if __name__ == "__main__":
    main()