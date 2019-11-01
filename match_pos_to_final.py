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
##################################################

"""
Script to match a list of positions to the radio catalogue
	- And, outptus a workflow.txt file with Source_Name for visual checks
	- Also, creates a subset of radio catalogue for these sources
"""

# Read in radio catalogue
final_fname = glob.glob("/disk3/rohitk/final_raido_catalogues/EN1/final-v*.fits")[-1]
mlfin_srl = Table.read(final_fname, character_as_bytes=False)

# Read in file with a list of positions
tab_to_match = Table.read("EN1_radio_cool_sources_pos.txt", format='ascii')

radio_coords = SkyCoord(mlfin_srl["RA"], mlfin_srl["DEC"], unit='deg', frame='icrs')
match_coords = SkyCoord(tab_to_match["RA"], tab_to_match["DEC"], unit='deg', frame='icrs')

# ind_r, ind_m, sep2d, _ = search_around_sky(radio_coords, match_coords, seplimit=5*u.arcsec)
ind_r, sep2d, _ = match_coordinates_sky(match_coords, radio_coords, nthneighbor=1)

mlfin_srl["Notes"] = -1
mlfin_srl["Sep"] = np.nan
mlfin_srl["Sep"][ind_r] = sep2d.arcsec
mlfin_srl["Source_Name", "RA", "DEC", "Sep", "z1_median", "Notes"][ind_r].write("workflow.txt", format='ascii.commented_header')
