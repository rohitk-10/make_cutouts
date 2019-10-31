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
import glob

# For catalog matching
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import search_around_sky

from astropy.table import Table
import os
##################################################

"""
Script to make radio source ellipses from the final catalogue based on a list of positions
"""


def make_ds9_reg(ra_positions, dec_positions, major_rad, minor_rad, pa, output_regfile, marker_colours):
    first_regline = "# Region file format: DS9 version 4.1"
    global_reg_def = 'global color={0} dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'.format(marker_colours)

    with open(output_regfile, "w") as fout:
        fout.write(first_regline + "\n")
        fout.write(global_reg_def + "\n")
        fout.write("fk5\n")

        for ii, ra in enumerate(ra_positions):
            fout.write('ellipse({0},{1},{2}",{3}",{4})\n'.format(ra, dec_positions[ii],
                                                                 major_rad[ii], minor_rad[ii], pa[ii]))

    return


# RA, DEC, MAJ rad, MIN RAD, angle(check this!)
# t['RA'],t['DEC'],t['Maj']*2/overlay_scale,t['Min']*2/overlay_scale, angle=90+t['PA']

# Read in radio catalogue
final_fname = glob.glob("/disk3/rohitk/final_raido_catalogues/EN1/final-v*.fits")[-1]
mlfin_srl = Table.read(final_fname, character_as_bytes=False)

# Workflow
cata = Table.read("workflow.txt", format='ascii')

wflow_coords = SkyCoord(cata["RA"], cata["DEC"], unit='deg', frame='icrs')
final_coords = SkyCoord(mlfin_srl["RA"], mlfin_srl["DEC"], unit='deg', frame='icrs')

srad = 300

OUTDIR_POS = "cata_pos/"
if not os.path.exists(OUTDIR_POS):
        os.makedirs(OUTDIR_POS)

for k in range(len(cata)):
    fin_fname = "{0}ellipse_{1}.reg".format(OUTDIR_POS, cata["Source_Name"][k])
    if not os.path.exists(fin_fname):
        print(fin_fname)
        cent_coord = SkyCoord(cata["RA"][k], cata["DEC"][k], unit='deg', frame='icrs')
        fin_near_cent = cent_coord.separation(final_coords).arcsec < srad
        # Now write ellipses for all sources in fin_near_cent
        make_ds9_reg(mlfin_srl["RA"][fin_near_cent].tolist(), mlfin_srl["DEC"][fin_near_cent].tolist(), (mlfin_srl["Maj"]*3600).tolist(),
                     (mlfin_srl["Min"]*3600).tolist(), (mlfin_srl["PA"]+90).tolist(), fin_fname, "cyan")
