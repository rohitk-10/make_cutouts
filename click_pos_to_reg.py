# coding: utf-8
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
                    fout.write('point({0},{1}) # point=diamond 20 color=magenta width=2\n'.format(ra, dec_positions[ii]))

    return
    
cata = Table.read("EN1_radio_cool_sources_pos.txt", format='ascii')
make_ds9_reg(cata["RA"].tolist(), cata["DEC"].tolist(), "EN1_radio_cool_sources_pos.reg", "magenta", None)

