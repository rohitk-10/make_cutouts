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

"""
Scipt to plot radio and optical/IR images onto ds9 and show:
	1. Radio Gaussians (in cyan ellipses)
	2. Radio contours (green regions)
	3. Radio source position (marked by red cross)
	4. Optical/IR catalogue detections (marked by small green circles)
	5. Optical host galaxy position found by LR or visual classification (cyan cross)
"""


def make_ds9_reg(ra_positions, dec_positions, output_regfile, marker_colours, lr_marker=None):
    first_regline = "# Region file format: DS9 version 4.1"
    global_reg_def = 'global color={0} dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'.format(marker_colours)

    with open(output_regfile, "w") as fout:
        fout.write(first_regline + "\n")
        fout.write(global_reg_def + "\n")
        fout.write("fk5\n")

        if lr_marker is "lrid":
            fout.write("point({0},{1}) # point=x 20 color={2} width=2\n".format(ra_positions, dec_positions, marker_colours))
        elif lr_marker is "radpos":
            fout.write("point({0},{1}) # point=cross 20 color={2} width=2\n".format(ra_positions, dec_positions, marker_colours))
        else:
            if marker_colours == "green":
                for ii, ra in enumerate(ra_positions):
                    fout.write('circle({0},{1},{2}")\n'.format(ra, dec_positions[ii], 0.3))

            else:
                for ii, ra in enumerate(ra_positions):
                    fout.write('box({0},{1},{2}",{3}",{4})\n'.format(ra, dec_positions[ii], 0.4, 0.4, 0.54026857))

    return


fname = "workflow.txt"
prefilt_out = Table.read(fname, format='ascii')
# prefilt_out = prefilt_out[prefilt_out["col1"] == 6]

all_ind = np.arange(len(prefilt_out))
inds_to_check = all_ind[prefilt_out["Notes"] == -1]

ind_to_plot = inds_to_check[0]

# mlfin_srl = Table.read("/disk1/rohitk/ELN1_project/eln1_workflow/iterated_endpoints/EN1_ML_RUN_fin_overlap_srl_workflow_fixed_fclean.fits")
final_fname = glob.glob("/disk3/rohitk/final_raido_catalogues/EN1/final-v*.fits")[-1]
mlfin_srl = Table.read(final_fname)

ra_s = mlfin_srl["RA"][mlfin_srl["Source_Name"] == prefilt_out["Source_Name"][ind_to_plot]]
dec_s = mlfin_srl["DEC"][mlfin_srl["Source_Name"] == prefilt_out["Source_Name"][ind_to_plot]]

cent_coord = SkyCoord(ra_s, dec_s, unit='deg', frame='icrs')

# Load in the optical and spitzer-chi2 catalogue
cata_opt = Table.read("/disk3/rohitk/ELAIS_opt_swarped/dual_analysis/combine/MASTER_catalogue/EN1_MASTER_opt_spitzer_merged_cedit_apcorr_adduncat_extra_run2_lite.fits")
# cata_mir = Table.read("/disk3/rohitk/ELAIS_opt_swarped/dual_sex/chi2salldet/18_12_2018_1/cat_chi2sall_g.fits", hdu=2)

opt_coords = SkyCoord(cata_opt["ALPHA_J2000"], cata_opt["DELTA_J2000"], unit='deg', frame='icrs')
# mir_coords = SkyCoord(cata_mir["ALPHA_J2000"], cata_mir["DELTA_J2000"], unit='deg', frame='icrs')

srad = 250.

opt_near_cent = cent_coord.separation(opt_coords).arcsec < srad
# mir_near_cent = cent_coord.separation(mir_coords).arcsec < srad

OUTDIR_POS = "cata_pos/"
if not os.path.exists(OUTDIR_POS):
    os.makedirs(OUTDIR_POS)

fin_fname = "{0}ellipse_{1}.reg".format(OUTDIR_POS, prefilt_out["Source_Name"][ind_to_plot])
if not os.path.exists(fin_fname):

    # Make a ds9 region file for both the optical and mir catalogue positions near this source
    make_ds9_reg(cata_opt["ALPHA_J2000"][opt_near_cent], cata_opt["DELTA_J2000"][opt_near_cent],
                 "{0}pos_{1}_opt_{2}.reg".format(OUTDIR_POS, ind_to_plot, prefilt_out["Source_Name"][ind_to_plot]), "green")
    # make_ds9_reg(cata_mir["ALPHA_J2000"][mir_near_cent], cata_mir["DELTA_J2000"][mir_near_cent],
    #              "{0}pos_{1}_mir_{2}.reg".format(OUTDIR_POS, ind_to_plot, prefilt_out["Source_Name"][ind_to_plot]), "red")

    # Get the position of the LR match too
    ra_lrid = mlfin_srl["ALPHA_J2000"][mlfin_srl["Source_Name"] == prefilt_out["Source_Name"][ind_to_plot]]
    dec_lrid = mlfin_srl["DELTA_J2000"][mlfin_srl["Source_Name"] == prefilt_out["Source_Name"][ind_to_plot]]
    make_ds9_reg(ra_lrid.tolist()[0], dec_lrid.tolist()[0], "{0}pos_{1}_LRID_{2}.reg".format(OUTDIR_POS, ind_to_plot, prefilt_out["Source_Name"][ind_to_plot]), "cyan", "lrid")

    # Get the region with the radio RA DEC position
    make_ds9_reg(ra_s.tolist()[0], dec_s.tolist()[0], "{0}pos_{1}_RADIO_{2}.reg".format(OUTDIR_POS, ind_to_plot, prefilt_out["Source_Name"][ind_to_plot]), "red", "radpos")

    print("Went through everything!")

# If instead all files exist, then skip the above and go straight to plotting

# Copy the latest positions files across to a different directory for ease of access
POS_LATEST = "cata_pos_latest/"
if os.path.exists(POS_LATEST):
    os.system("rm -rf " + POS_LATEST)
os.makedirs(POS_LATEST)

CONTOUR_PATH = "contour_regions"

os.system("cp {0}pos_{1}_opt_{2}.reg {3}".format(OUTDIR_POS, ind_to_plot, prefilt_out["Source_Name"][ind_to_plot], POS_LATEST))
os.system("cp {0}pos_{1}_LRID_{2}.reg {3}".format(OUTDIR_POS, ind_to_plot, prefilt_out["Source_Name"][ind_to_plot], POS_LATEST))
os.system("cp {0}pos_{1}_RADIO_{2}.reg {3}".format(OUTDIR_POS, ind_to_plot, prefilt_out["Source_Name"][ind_to_plot], POS_LATEST))
os.system("cp {0}/cont_{1}.reg {2}".format(CONTOUR_PATH, prefilt_out["Source_Name"][ind_to_plot], POS_LATEST))
os.system("cp {0}/ellipse_{1}.reg {2}".format(OUTDIR_POS, prefilt_out["Source_Name"][ind_to_plot], POS_LATEST))
# os.system("cp {0}pos_{1}_mir_{2}.reg {3}".format(OUTDIR_POS, ind_to_plot, prefilt_out["Source_Name"][ind_to_plot], POS_LATEST))

#################################################################


print("***** Plotting source: #{0}: {1} *****\n".format(ind_to_plot, prefilt_out["Source_Name"][ind_to_plot]))

# Show the ds9 plot of sources
os.system("ds9 {0}/{1}_radio.fits -zscale {0}/{1}_i.fits -zscale {0}/{1}_sw2.fits -zscale -match frame wcs -regions load all '{2}*.reg' &".format("opt_mir_imgs", prefilt_out["Source_Name"][ind_to_plot], POS_LATEST))

# Also the png plot
# os.system("display {0}_j.png & ".format(prefilt_out["Source_Name"][ind_to_plot]))

# Now ask for input on the flag for this source
print("***** Now enter some flag here which you can define later ***** \n",
      "***** This will be added to 'Notes' column 		 *****\n")

if "Sep" in prefilt_out.colnames:
    print("Separation between radio source and visual click: {0} arcsec".format(prefilt_out["Sep"][ind_to_plot]))
if "z1_median" in prefilt_out.colnames:
    print("Redshift of radio source (best possible radio source): {0}".format(prefilt_out["z1_median"][ind_to_plot]))

flag_update = int(float(input("Enter the flag for this source. If not sure, enter '-2' to move on: ")))

if (flag_update > 70) & (flag_update < 100):
    print("Unsuitable flag entered, setting {0} to -2".format(prefilt_out["Source_Name"][ind_to_plot]))
    prefilt_out["Notes"][ind_to_plot] = -2
else:
    prefilt_out["Notes"][ind_to_plot] = int(flag_update)

# Save a backup of the file before overwriting
os.system("cp {0} {0}.last".format(fname))

prefilt_out.write(fname, format='ascii', overwrite=True)
