# Setup codes
ipython click_pos_to_reg.py  	# Take a list of radio source click positions and write to ds9 region
ipython match_pos_to_final.py	# Take output of previous code and match to NN radio sources and write to workflow.txt
ipython make_radio_ellipses.py	# Make ellipses of radio sources in workflow.txt and write to ds9 region file
ipython cont_to_reg.py  	# Make contours of sources in workflow.txt
ipython crop_im_additional.py	# Make cutouts of radio and N optical/IR images
# Now that everything is set-up, run plot_ds9.py to view each source
# ipython plot_ds9.py
