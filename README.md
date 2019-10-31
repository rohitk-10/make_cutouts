# Codes to make and display cutouts of FITS images

## Running the code
First generate a list of sources of interest and put them in 'workflow.txt'
- cont_to_reg.py : Read in source positions from 'workflow.txt' and use radio image + rms map to generate contours and write to ds9 region file
- crop_im_additional.py : At source positions in 'workflow.txt', make cutouts of radio, optical and IR (and any other supplied) FITS data
- make_radio_ellipses.py : Generate ellipses based on radio source position and sizes, etc. to overlay 
- plot_ds9.py : Convenience code to use outputs from previous codes and plot cutouts + contours to ds9
