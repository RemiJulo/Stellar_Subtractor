###################
# PARAMETERS FILE #
###################

class argmax: pass

###################################################################################################
###################################################################################################
# List of paths to .fits files (or to folders full of .fits files) to be post-processed

input_paths = ["data/pre_processed/"]

###################################################################################################
# Path to the (new) saving folder from the input file locations
# Any needed series of folders will be created in the operation

output_path = "../post_processed/"

###################################################################################################
###################################################################################################
# Origin coordinates relocation ('int' for pixel values or 'float' for physical values)
# Class 'argmax' is also usable to recover the coordinates of the brightest pixel
# Origin coordinates are printed in the terminal if 'print_verbose' = 1
# Built as a dictionary associating axis names with origin coordinates
# It is also possible to choose units with the same axis names as keys
# (the list of compatible units is printed in case of bad unit values)
# These units are also the ones used in the terminal

origin = {

    'RA' : argmax,
    'DEC' : argmax,

    'AWAV' : 0.0

}

origin_units = {

    'AWAV' : 'Å'

}

###################################################################################################
# Zooms along the different axes ('int' for pixel values or 'float' for physical values)
# Built as a dictionary associating axis names with pairs of lower and upper zoom bounds
# For instance, the entry " 'RA' : (-10,0.01) " allows to zoom in on the right ascension
# 10 pixels east of the origin of 'RA' to 0.01 deg west of the origin of 'RA'
# Axis names and units are printed in the terminal if 'print_verbose' = 1
# It is also possible to choose units with the same axis names as keys
# (the list of compatible units is printed in case of bad unit values)
# Both the upper and lower bounds are included in the selection

zooms = {
    
    'RA' : (-0.6,0.6),
    'DEC' : (-0.6,0.6),

    'AWAV' : (6462.8,6662.8)

}

zooms_units = {

    'RA' : 'arcsec',
    'DEC' : 'arcsec',

    'AWAV' : 'Å'

}

###################################################################################################
# Choice of the grids units for all the following parameters and shows
# (the list of compatible units is printed in case of bad unit values)

grids_units = {

    'RA' : 'arcsec',
    'DEC' : 'arcsec',

    'AWAV' : 'Å'
    
}

###################################################################################################
###################################################################################################
# Height of the stellar spread compared to what is consequently defined as the 1-sigma noise level
# The stellar component of the pixels with a flux below this threshold is assumed to be negligible
# The corresponding pixels are not used for estimating the stellar spectrum by spatial integration
# The unused pixels are masked in the image of the stellar spectrum estimate if 'show_verbose' = 1
# This helps to limit electronic noise contamination, especially with large selected field of view
# Too large values of 'stellar_sigma' lead to a selection of too many too noisy pixels in practice
# (thus could lead to strange results in the following because of various numerical instabilities)

stellar_sigma = 500

###################################################################################################
# Misleading a priori in the form of a dictionary associating axis names with lists of doublets
# between which pixels are masked (and disregarded) during the stellar spectrum estimation step

stellar_mask = {



}

###################################################################################################
###################################################################################################
# Instrument-dependent hyperparameter, dealing with the smoothness of its point spread function
# Singleton to set the degree of a Legendre polynomial to fit the spectral behaviour of the PSF
# Empty set to estimate directly the PSF by division of the data by a stellar spectrum estimate

psf_degrees = [4]

###################################################################################################
# Misleading a priori in the form of a dictionary associating axis names with lists of doublets
# between which pixels are masked (and disregarded) during the point spread function estimation

psf_mask = {

    'AWAV' : [(6552.8,6572.8)],

}

###################################################################################################
###################################################################################################
# 'save_verbose' for backup feedback saved in 'saving_folders' ('0' / '1' for no / full display)
# 'print_verbose' for textual feedback printed in the terminal ('0' / '1' for no / full display)
# 'show_verbose' for visual feedback showed in .pyplot windows ('0' / '1' for no / full display)

save_verbose = 1

print_verbose = 1

show_verbose = 1

###################################################################################################
# Data filters for visual feedbacks (as products from Gaussians with height 1)
# Built as a dictionary associating axis names with lists of (location, scale)
# Gaussians parameters doublets with which to filter images and spectra
# ('int' for pixel values or 'float' for physical values)

spectra_filter = {

}

images_filter = {

    'AWAV' : [(6562.8,1.25)],

}

###################################################################################################
# Data blocker for visual feedbacks (as products from 1 - Gaussians with height 1)
# Built as a dictionary associating axis names with lists of (location, scale)
# Gaussians parameters doublets with which to block images and spectra
# ('int' for pixel values or 'float' for physical values)

spectra_blocker = {

    'RA' : [(0.0,0.1)],
    'DEC' : [(0.0,0.1)],

}

images_blocker = {

    'RA' : [(0.0,0.1)],
    'DEC' : [(0.0,0.1)],

}

###################################################################################################
###################################################################################################
