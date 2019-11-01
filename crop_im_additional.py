import numpy as np
from astropy.io import fits
from astropy.wcs import wcs
import glob

# For catalog matching
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import search_around_sky

from matplotlib import pyplot as plt
from astropy.table import Table

import time
import os
###################################################


def crop_write_fits(hdu_list, original_wcs, original_data, xcrop_array, ycrop_array, outfname):
    """
    Write HDUList to a new file with updated wcs and data
    """

    # Crop the data
    updated_data = original_data[yr, xr]

    # Crop the wcs information
    updated_wcs = original_wcs[yr, xr]

    # Print out the shape of the cropped data
    print("Cropped image has shape: " + str(np.shape(updated_data)))

    # Create a new fits PrimaryHDU class to store the data and header
    newf = fits.PrimaryHDU()
    # Assign the data to the FITS hdulist object
    newf.data = updated_data

    # Add the old header information
    newf.header = hdu_list[0].header

    # But, update the WCS to the cropped region
    newf.header.update(updated_wcs.to_header())

    # Write to FITS file
    newf.writeto(outfname, overwrite=True)
    return


img_dict = dict()
img_dict["i"] = "/disk3/rohitk/ELAIS_opt_swarped/iband_fits/EN1band_swarped/final/EL_EN1_iband.fits"
img_dict["sw2"] = "/disk3/rohitk/ELAIS_opt_swarped/sw2band_fits/EN1band_swarped/final/EL_EN1_sw2band.fits"
img_dict["radio"] = "/disk1/rohitk/ELN1_project/ELAIS-N1/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked-crop.fits"

# These are the optical and IR chi2 data
# img_dict["chi2o"] = "/disk3/rohitk/ELAIS_opt_swarped/chi2_ind_ugrizJK/EN1band_swarped/final/EL_EN1_chi2_ugrizJK.fits"
# img_dict["chi2s"] = "/disk3/rohitk/ELAIS_opt_swarped/chi2_ind_swse/EN1band_swarped/final/EL_EN1_chi2_swse.fits"


filts = list(img_dict.keys())

# Load in the txt file containing the source names and IDs
# mlfin_srl = Table.read("/disk1/rohitk/ELN1_project/eln1_workflow/iterated_endpoints/EN1_ML_RUN_fin_overlap_srl_workflow_fixed.fits")
final_fname = glob.glob("/disk3/rohitk/final_raido_catalogues/EN1/final-v*.fits")[-1]
mlfin_srl = Table.read(final_fname, character_as_bytes=False)


# prefilt_out = Table.read("../prefilt_output_test/workflow_rk_testing.txt", format='ascii')
# prefilt_out = prefilt_out[prefilt_out["col1"] == 6]

prefilt_out = Table.read("workflow.txt", format='ascii')
# prefilt_out = prefilt_out[(prefilt_out["Notes"] != -99) & (prefilt_out["Notes"] != -1)]

# print("No. of new imgs to crop: {0}".format(len(prefilt_out)))

ts = time.time()

BASE_OUT = "opt_mir_imgs"
if not os.path.exists(BASE_OUT):
    os.makedirs(BASE_OUT)

for ii in range(len(prefilt_out)):

    sname_indx = mlfin_srl["Source_Name"] == prefilt_out["Source_Name"][ii]
    ra = mlfin_srl["RA"][sname_indx]
    dec = mlfin_srl["DEC"][sname_indx]

    cent_coord = SkyCoord(ra, dec, unit='deg', frame='icrs')

    pa_pattern = [-135, 45]
    sep_pattern = [200, 200]
    corner_coord = cent_coord.directional_offset_by(pa_pattern*u.deg, sep_pattern*u.arcsec)

    for phot_band in filts:

        if phot_band == "radio":
            sep_pattern = [400, 400]
            corner_coord = cent_coord.directional_offset_by(pa_pattern*u.deg, sep_pattern*u.arcsec)
        else:
            sep_pattern = [300, 300]
            corner_coord = cent_coord.directional_offset_by(pa_pattern*u.deg, sep_pattern*u.arcsec)
            # Name the section based on Source_Name and phot_band
        crop_img_name = BASE_OUT + "/" + prefilt_out["Source_Name"][ii] + "_" + phot_band + ".fits"
        img_to_crop = img_dict[phot_band]

        if not os.path.exists(crop_img_name):
            hdul = fits.open(img_to_crop)
            img_wcs = wcs.WCS(hdul[0].header, hdul).celestial

            print(crop_img_name)

            xc, yc = img_wcs.all_world2pix(corner_coord.ra, corner_coord.dec, 0)
            xc = np.sort(xc)
            yc = np.sort(yc)

            xr = slice(int(np.round(xc[0])), int(np.floor(xc[1])))
            yr = slice(int(np.round(yc[0])), int(np.floor(yc[1])))

            if phot_band != "radio":
                crop_write_fits(hdul, img_wcs, hdul[0].data, xr, yr, crop_img_name)
            else:
                crop_write_fits(hdul, img_wcs, hdul[0].data[0, 0, :, :], xr, yr, crop_img_name)

print(time.time() - ts)
