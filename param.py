###################
# PARAMETERS FILE #
###################

class argmax: pass

###################################################################################################
###################################################################################################
# List of paths to .fits files (or to folders full of .fits files) to be post-processed

input_paths = ["data/"]

###################################################################################################
# Path to the (new) saving folder from the input file locations
# Any needed series of folders will be created in the operation

output_path = "post_processing"

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
    
    'AWAV' : (6062.8,7062.8)

}

zooms_units = {

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

stellar_sigma = 100

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
# 'save_verbose' for backup feedback saved in output directory ('0' / '1' for no / full display)
# 'print_verbose' for textual feedback printed in the terminal ('0' / '1' for no / full display)
# 'show_verbose' for visual feedback showed in .pyplot windows ('0' / '1' for no / full display)

save_verbose = 1

print_verbose = 1

show_verbose = 1

###################################################################################################
# Saturation of the feedback images as quantiles of their pixel values (if show_verbose = 1)
# Built as the doublet of the 'vmin' then 'vmax' saturation quantile

vminmax_quantiles = (0.01,0.99)

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

    'RA' : [(0.0,0.25)],
    'DEC' : [(0.0,0.25)],

}

images_blocker = {

    'RA' : [(0.0,0.25)],
    'DEC' : [(0.0,0.25)],

}

################################################################################################### PARAMETERS CHECKING
################################################################################################### -- DO NOT MODIFY --

if not type(input_paths) is list:
    raise ValueError("\x1B[31m" + "'input_paths' should be a list of 'str'" + "\x1B[0m")
for input_path in input_paths:
    if not type(input_path) is str:
        raise ValueError("\x1B[31m" + "'input_paths' should be a list of 'str'" + "\x1B[0m")

if not output_path or not type(output_path) is str:
    raise ValueError("\x1B[31m" + "'output_path' should be a un-empty 'str'" + "\x1B[0m")

if not type(origin) is dict:
    raise ValueError("\x1B[31m" + "'origin' should be a 'dict'" + "\x1B[0m")
for key, value in origin.items():
    if not key or not type(key) is str:
        raise ValueError("\x1B[31m" + "'origin' keys should be un-empty 'str'" + "\x1B[0m")
    if not value is argmax and not type(value) in (int, float):
        raise ValueError("\x1B[31m" + "'origin' values should be 'int/float/arg'" + "\x1B[0m")
if not type(origin_units) is dict:
    raise ValueError("\x1B[31m" + "'origin_units' should be a 'dict'" + "\x1B[0m")
for key, value in origin_units.items():
    if not key or not type(key) is str:
        raise ValueError("\x1B[31m" + "'origin_units' keys should be un-empty 'str'" + "\x1B[0m")
    if not type(value) is str:
        raise ValueError("\x1B[31m" + "'origin_units' values should be 'str'" + "\x1B[0m")

if not type(zooms) is dict:
    raise ValueError("\x1B[31m" + "'zooms' should be a 'dict'" + "\x1B[0m")
for key, value in zooms.items():
    if not key or not type(key) is str:
        raise ValueError("\x1B[31m" + "'zooms' keys should be un-empty 'str'" + "\x1B[0m")
    if not type(value) is tuple or not len(value) == 2:
        raise ValueError("\x1B[31m" + "'zooms' values should be 'int/float'-doublets" + "\x1B[0m")
    if not type(value[0]) in (int, float) or not type(value[1]) in (int, float):
        raise ValueError("\x1B[31m" + "'zooms' values should be 'int/float'-doublets" + "\x1B[0m")
if not type(zooms_units) is dict:
    raise ValueError("\x1B[31m" + "'zooms_units' should be a 'dict'" + "\x1B[0m")
for key, value in zooms_units.items():
    if not key or not type(key) is str:
        raise ValueError("\x1B[31m" + "'zooms_units' keys should be un-empty 'str'" + "\x1B[0m")
    if not type(value) is str:
        raise ValueError("\x1B[31m" + "'zooms_units' values should be 'str'" + "\x1B[0m")

if not type(stellar_sigma) in (int, float):
    raise ValueError("\x1B[31m" + "'stellar_sigma' should be 'int/float'" + "\x1B[0m")

if not type(stellar_mask) is dict:
    raise ValueError("\x1B[31m" + "'stellar_mask' should be a 'dict'" + "\x1B[0m")
for key, value in stellar_mask.items():
    if not key or not type(key) is str:
        raise ValueError("\x1B[31m" + "'stellar_mask' keys should be un-empty 'str'" + "\x1B[0m")
    if not type(value) is list:
        raise ValueError("\x1B[31m" + "'stellar_mask' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not type(entry) is tuple or not len(entry) == 2:
            error_message = "'stellar_mask' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not type(entry[0]) in (int, float) or not type(entry[1]) in (int, float):
            error_message = "'stellar_mask' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

if not type(psf_degrees) is list:
    raise ValueError("\x1B[31m" + "'psf_degrees' should be a list" + "\x1B[0m")
for entry in psf_degrees:
    if not type(entry) is int or not entry > 0:
        raise ValueError("\x1B[31m" + "'psf_degrees' entries should be a 'int' > 0" + "\x1B[0m")

if not type(psf_mask) is dict:
    raise ValueError("\x1B[31m" + "'psf_mask' should be a 'dict'" + "\x1B[0m")
for key, value in psf_mask.items():
    if not key or not type(key) is str:
        raise ValueError("\x1B[31m" + "'psf_mask' keys should be un-empty 'str'" + "\x1B[0m")
    if not type(value) is list:
        raise ValueError("\x1B[31m" + "'psf_mask' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not type(entry) is tuple or not len(entry) == 2:
            error_message = "'psf_mask' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not type(entry[0]) in (int, float) or not type(entry[1]) in (int, float):
            error_message = "'psf_mask' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

if not save_verbose in [0,1]:
    raise ValueError("\x1B[31m" + "'save_verbose' should be in {'0','1'}" + "\x1B[0m")
if not print_verbose in [0,1]:
    raise ValueError("\x1B[31m" + "'print_verbose' should be in {'0','1'}" + "\x1B[0m")
if not show_verbose in [0,1]:
    raise ValueError("\x1B[31m" + "'show_verbose' should be in {'0','1'}" + "\x1B[0m")

if not type(vminmax_quantiles) is tuple and not len(vminmax_quantiles) == 2:
    raise ValueError("\x1B[31m" + "'vminmax_quantiles' should be a doublet" + "\x1B[0m")
if not type(vminmax_quantiles[0]) in (int,float) and not type(vminmax_quantiles[1]) in (int,float):
    raise ValueError("\x1B[31m" + "'vminmax_quantiles' entries should be 'int/float'" + "\x1B[0m")

if not type(spectra_filter) is dict:
    raise ValueError("\x1B[31m" + "'spectra_filter' should be a 'dict'" + "\x1B[0m")
for key, value in spectra_filter.items():
    if not key or not type(key) is str:
        raise ValueError("\x1B[31m" + "'spectra_filter' keys should be un-empty 'str'" + "\x1B[0m")
    if not type(value) is list:
        raise ValueError("\x1B[31m" + "'spectra_filter' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not type(entry) is tuple or not len(entry) == 2:
            error_message = "'spectra_filter' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not type(entry[0]) in (int, float) or not type(entry[1]) in (int, float):
            error_message = "'spectra_filter' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
if not type(images_filter) is dict:
    raise ValueError("\x1B[31m" + "'images_filter' should be a 'dict'" + "\x1B[0m")
for key, value in images_filter.items():
    if not key or not type(key) is str:
        raise ValueError("\x1B[31m" + "'images_filter' keys should be un-empty 'str'" + "\x1B[0m")
    if not type(value) is list:
        raise ValueError("\x1B[31m" + "'images_filter' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not type(entry) is tuple or not len(entry) == 2:
            error_message = "'images_filter' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not type(entry[0]) in (int, float) or not type(entry[1]) in (int, float):
            error_message = "'images_filter' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

if not type(spectra_blocker) is dict:
    raise ValueError("\x1B[31m" + "'spectra_blocker' should be a 'dict'" + "\x1B[0m")
for key, value in spectra_blocker.items():
    if not key or not type(key) is str:
        raise ValueError("\x1B[31m" + "'spectra_blocker' keys should be un-empty 'str'" + "\x1B[0m")
    if not type(value) is list:
        raise ValueError("\x1B[31m" + "'spectra_blocker' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not type(entry) is tuple or not len(entry) == 2:
            error_message = "'spectra_blocker' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not type(entry[0]) in (int, float) or not type(entry[1]) in (int, float):
            error_message = "'spectra_blocker' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
if not type(images_blocker) is dict:
    raise ValueError("\x1B[31m" + "'images_blocker' should be a 'dict'" + "\x1B[0m")
for key, value in images_blocker.items():
    if not key or not type(key) is str:
        raise ValueError("\x1B[31m" + "'images_blocker' keys should be un-empty 'str'" + "\x1B[0m")
    if not type(value) is list:
        raise ValueError("\x1B[31m" + "'images_blocker' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not type(entry) is tuple or not len(entry) == 2:
            error_message = "'images_blocker' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not type(entry[0]) in (int, float) or not type(entry[1]) in (int, float):
            error_message = "'images_blocker' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
