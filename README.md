# Codes to make and display cutouts of FITS images

## Codes
First generate a list of sources of interest and put them in 'workflow.txt' using two codes below
- click_pos_to_reg.py  	# Take a list of radio source click positions and write to ds9 region
- match_pos_to_final.py	# Take output of previous code and match to NN radio sources and write to workflow.txt
Now use the workflow.txt for the rest of the codes:
- cont_to_reg.py : Read in source positions from 'workflow.txt' and use radio image + rms map to generate contours and write to ds9 region file
- crop_im_additional.py : At source positions in 'workflow.txt', make cutouts of radio, optical and IR (and any other supplied) FITS data
- make_radio_ellipses.py : Generate ellipses based on radio source position and sizes, etc. to overlay 
- plot_ds9.py : Convenience code to use outputs from previous codes and plot cutouts + contours to ds9 - run multiple times for each new source

### For Rachel:
- plot_ds9_rachel.py : Convenience code to use outputs from previous codes and plot cutouts + contours to ds9 - has small changes for your needs

## Running plot_ds9.py
This will show the radio and optical/IR cutouts of a single radio source. Additional overlays: Nearby radio source contours (green), nearby radio source ellipses (cyan ellipses), radio source of interest (marked by red plus), optical-ID of radio source of interest (cyan cross), all nearby optical/IR catalogue detections (small green circles). And, click positions (red diamond) - only in plot_ds9.py

- Also allows you to flag each source and update to "Notes" in workflow.txt. Use any flag value (except -1 and -2). -1 reserved for sources that haven't been inspected yet. Use -2 if you're not sure what to flag and want to move on to next source.

Code can be run in two modes:
1. Specify source to inspect:
	- "python plot_ds9.py -s source_name" : Here, 'source_name' is the Source_Name of radio source
2. Work through each source in workflow.txt
	- "python plot_ds9.py" : Will work through each source with "Notes" = -1 in workflow.txt. Need to run code again to move on to next source.
