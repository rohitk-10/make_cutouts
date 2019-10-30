# Script to write the contours segments from pyplot to ds9 region file for around a given LOFAR source
# Commands to select directory based on hostname
from socket import gethostname

if gethostname() == 'colonsay':
    path_start = '/disk1/rohitk/ELN1_project/'
elif gethostname() == 'rohitk-elitebook':
    path_start = '/home/rohitk/Documents/PhD/Year1/ELN1_project/'

#################################################
# Add the path of useful functions at the start
import sys
sys.path.append(path_start+'basic_functions')
from useful_functions import return_hist_par, varstat, latest_dir, jytoabmag, field_filter
from plot_func import rc_def, make_fig, make_fig_multi
rc_def()
##################################################

import numpy as np
from matplotlib import pyplot as plt

from astropy.io import fits
from astropy.wcs import wcs, WCS
from os import remove
import time

import os
from astropy.table import Table
import glob

############################################


def cont_to_ds9reg(lofar_ra, lofar_dec, regfile_name, half_region_size, nlevels, im_data, im_wcs, rms_imdata):
    """
    Generate contours using pyplot and then parse them as ds9 polygon regions and write to file

    Inputs:
        lofar_ra:         RA (deg) of the LOFAR source
        lofar_dec:        DEC (deg) of the LOFAR source
        regfile_name:     Output filename of the region file
        half_region_size: Half the size of the region on which
        nlevels         : Number of contour levels

    Outputs:
        ds9 .reg file written to file 
    """
    # Convert world to pixel coordinates
    lof_x, lof_y = im_wcs.all_world2pix(lofar_ra, lofar_dec, 0)

    lof_x = int(lof_x)
    lof_y = int(lof_y)

    # Crop around this region
    cropped_region = im_data[lof_y - half_region_size:lof_y + half_region_size, lof_x - half_region_size: lof_x + half_region_size]

    # Generate a cropped region around the source in the RMS map and take the median to get the RMS value for contour levels
    rms_crop_region = rms_imdata[lof_y - 100:lof_y + 100, lof_x - 100: lof_x + 100]
    lofar_rms = np.nanmedian(rms_crop_region)

    # plt.imshow(cropped_region, cmap='rainbow', vmin=0.00008, vmax=0.0009)
    # Scale factor to generate contour levels from RMS
    sf = 3

    # Generate levels for this LOFAR source based on the RMS values
    cont_levels = np.sqrt(2)**(np.arange(nlevels))

    cont_levels = cont_levels * lofar_rms * sf
    cont_levels[0] = -1              # Set the first contour level to show "holes"

    # Generate the contour colours for each level
    cont_colours = ['green' for i in range(len(cont_levels))]
    cont_colours[0] = 'red'        # Green for the -1 level. The rest are coloured white

    # Get the size of the image to get x and y contour
    naxis1 = np.size(cropped_region, axis=0)
    naxis2 = np.size(cropped_region, axis=1)

    # Generate the contours
    cont = plt.contour(np.arange(naxis1), np.arange(naxis2), cropped_region, levels=cont_levels, linewidths=2.)

    # Get the contour segments
    cont_segs = cont.allsegs

    # Remove the file if it already exists
    try:
        remove(regfile_name)
    except OSError:
        pass

    # Generate the region file of the polygons for each level
    with open(regfile_name, 'w') as regout:
        regout.write("# Region file format: DS9 version 4.1\n")
        regout.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        regout.write("fk5\n")

        # Now write in a for loop, the polygons for each level
        for ii in range(len(cont_segs)):
            # Get the segments of the ii th level
            level_segs = cont_segs[ii]

            # For each polygon array in a given level
            for k in range(len(level_segs)):
                # Shift the pixel coordinates back to original positions before writing to file
                kth_polygon = level_segs[k]
                kth_polygon[:,0] = kth_polygon[:,0] + lof_x - half_region_size
                kth_polygon[:,1] = kth_polygon[:,1] + lof_y - half_region_size

                # Convert to world coordinates
                ra_cont,dec_cont = im_wcs.all_pix2world(kth_polygon[:,0],kth_polygon[:,1],0)

                regout.write("polygon(")

                # For each row in the polygon array
                # Write the first row out of the for loop
                regout.write(str(ra_cont[0]) + ',' + str(dec_cont[0]))

                for aa in range(1,len(ra_cont)):
                    regout.write(',' + str(ra_cont[aa]) + ',' + str(dec_cont[aa]))

                regout.write(") # color=" + cont_colours[ii] + '\n')

    # End of writing the region file
    return


# Make a function that reads in a FITS image and returns the wcs and image
def read_fits(fits_filename):
    """
    Function to read the FITS file and store the image array and WCS

    Inputs:
        fits_filename:        Name of the FITS filename
    Outputs:
        wcs_info:             WCS of the FITS image
        image_array:          Array of the FITS image
    """
    list_of_hdu = fits.open(fits_filename)

    # Get the image array
    image_array = fits.getdata(fits_filename)
    image_array = image_array[0][0]

    # Get the WCS of the image
    wcs_info = wcs.WCS(list_of_hdu[0].header, list_of_hdu).celestial

    return image_array, wcs_info


######################################################################################


LOFAR_IM_PATH = "/disk1/rohitk/ELN1_project/ELAIS-N1/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked-crop.fits"
LOFAR_RMS_PATH = "/disk1/rohitk/ELN1_project/ELAIS-N1/image_full_ampphase_di_m.NS_shift.blanked.scaled.rms.fits"

# Read in the image and the RMS map and store the data
lofar_image_data, lofar_im_wcs = read_fits(LOFAR_IM_PATH)
hdul_rms = fits.open(LOFAR_RMS_PATH)
rms_image_data = hdul_rms[0].data[0, 0, :, :]


# Read in the catalogue
cata = Table.read("workflow.txt", format='ascii')
sname_all = cata["Source_Name"]

final_fname = glob.glob("/disk3/rohitk/final_raido_catalogues/EN1/final-v*.fits")[-1]
mlfin_srl = Table.read(final_fname, character_as_bytes=False)

t1 = time.time()

BASE_OUT = "contour_regions"
if not os.path.exists(BASE_OUT):
    os.makedirs(BASE_OUT)

for ii, sname in enumerate(sname_all):
    print(sname)
    reg_fnames = "{0}/cont_{1}.reg".format(BASE_OUT, sname)
    full_cat_ind = np.isin(mlfin_srl["Source_Name"], sname)

    cont_to_ds9reg(mlfin_srl["RA"][full_cat_ind],mlfin_srl["DEC"][full_cat_ind], reg_fnames, 400, 9, lofar_image_data, lofar_im_wcs, rms_image_data)

print(time.time() - t1)

# plt.show()
