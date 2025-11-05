###################
# PARAMETERS FILE #
###################

class argmax: pass

###################################################################################################
###################################################################################################
# List of paths to .fits files (or to folders full of .fits files) to be post-processed

inputs_paths = ["inputs/"]

###################################################################################################
# Path to the saving folder from the location of the code files
# Any needed series of folders will be created in the operation
# (not ending with '/' or '\' allows to add a prefix)

outputs_path = "outputs/"

###################################################################################################
###################################################################################################
# Origin coordinates relocation ('int' for pixel values or 'float' for physical values)
# Herewith, 'argmax()' is also usable to recover the coordinates of the brightest pixel
# For instance, the entry " 'AWAV' : 0.0 " allows to replace the air wavelengths origin
# from the air wavelength of the first spectral channel to the one of ~ 0.0 m
# Built as a dictionary associating axis names with origin coordinates
# It is also possible to change the units used for the physical values
# (the list of compatible units is printed in case of bad unit values)

origin = {

    'AWAV' : 0.0

}

origin_units = {

    'AWAV' : 'Å'

}

###################################################################################################
# Zooms along the different axes ('int' for pixel values or 'float' for physical values)
# Built as a dictionary associating axis names with lower and upper zoom bounds doublets
# For instance, the entry " 'RA' : (0.01,-10) " allows to zoom in on the right ascension
# 0.01 deg east of the origin of 'RA' to 10 pixels east of the right 'RA' border
# It is also possible to change the units used for the physical values
# (the list of compatible units is printed in case of bad unit values)
# These units are used too in all the following dictionaries

zooms = {
    
    'RA' : (0.6,-0.6),
    'DEC' : (-0.6,0.6),

    'AWAV' : (6462.8,6762.8)

}

zooms_units = {

    'RA' : 'arcsec',
    'DEC' : 'arcsec',

    'AWAV' : 'Å'

}

###################################################################################################
###################################################################################################
# Detection of the outliers for their masking during the full duration of the post-processings
# Built as a dictionary associating axis names with numbers of standard deviations above means

outliers_sigmas = {



}

###################################################################################################
# Misleading a priori in the form of a dictionary associating axis names with lists of doublets
# between which pixels are masked (and disregarded) during the outliers detection and selection
# The final global mask is the intersection of the sub total masks defined for each of the axes

outliers_mask = {



}

###################################################################################################
###################################################################################################
# Height of the stellar level compared to the one of the noise level (setting noise/stellar areas)
# The stellar component of the pixels with a flux below this threshold is assumed to be negligible
# The corresponding pixels are not used for estimating the stellar spectrum by spatial integration
# The unused pixels are masked in the image of the stellar spectrum estimate if 'show_verbose' = 1
# This helps to limit electronic noise contamination, especially with large selected field of view
# Too large values of 'stellar_sigma' lead to a selection of too many too noisy pixels in practice
# (thus could lead to uncanny results in the following because of various numerical instabilities)

stellar_level = 50

###################################################################################################
# Misleading a priori in the form of a dictionary associating axis names with lists of doublets
# between which pixels are masked (and disregarded) during the stellar spectrum estimation step
# The final global mask is the intersection of the sub total masks defined for each of the axes

stellar_mask = {



}

###################################################################################################
###################################################################################################
# Instrument-dependent hyperparameters dealing with the smoothness of its point spread function
# Empty set to estimate directly the PSF by division of the data by a stellar spectrum estimate
# Singleton to set the degree of a Legendre polynomial to fit the spectral behaviour of the PSF
# Doublet to set the degrees of Legendre polynomials to fit the behaviour of the PSF as a whole
# (this last option being only realistic and thus relevant in the next version of 'process.py')

psf_degrees = [4]

###################################################################################################
# Misleading a priori in the form of a dictionary associating axis names with lists of doublets
# between which pixels are masked (and disregarded) during the point spread function estimation
# The final global mask is the intersection of the sub total masks defined for each of the axes

psf_mask = {

    'AWAV' : [(6552.8,6572.8)],

}

###################################################################################################
###################################################################################################
# Number of PSF estimation loopback (taking into account an increasingly better noise estimate)

noise_loopback = 0

###################################################################################################
# Misleading a priori in the form of a dictionary associating axis names with lists of doublets
# between which pixels are masked (and disregarded) during the 2 covariance matrices estimation
# The final global mask is the intersection of the sub total masks defined for each of the axes

noise_mask = {



}

###################################################################################################
###################################################################################################
# Detection of the anomalies for their highlight in the post-processing report in the terminal
# Built as a dictionary associating axis names with numbers of standard deviations above means

anomalies_sigmas = {



}

###################################################################################################
# Misleading a priori in the form of a dictionary associating axis names with lists of doublets
# between which pixels are masked (and disregarded) during the anomalies detection and highight
# The final global mask is the intersection of the sub total masks defined for each of the axes

anomalies_mask = {



}

###################################################################################################
# 'save_verbose' for backup feedback in 'outputs_path' ('0' / '1' for no / full display)
# 'print_verbose' for textual feedback in the terminal ('0' / '1' for no / full display)
# 'show_verbose' for visual feedback in python windows ('0' / '1' for no / full display)

save_verbose = 1

print_verbose = 1

show_verbose = 1

###################################################################################################
# Data filters for visual feedbacks (as products with unit-height Gaussians)
# Built as dictionaries associating axis names with (location, scale) lists
# of Gaussians parameters doublets with which to filter images and spectra

spectra_filter = {



}

images_filter = {



}

###################################################################################################
# Data maskers for visual feedbacks (as NaNing of the intersection of the selected windows)
# Built as dictionaries associating axis names with lower and upper bounds doublets lists
# between which pixels are masked just before the spectra and images means calculations

spectra_masker = {



}

images_masker = {



}

###################################################################################################
# Data hinders for visual feedbacks (as products from unit-height Gaussians)
# Built as dictionaries associating axis names with (location, scale) lists
# of Gaussians parameters doublets with which to hinder images and spectra

spectra_hinder = {



}

images_hinder = {



}

################################################################################################### PARAMETERS CHECKING
################################################################################################### -- DO NOT MODIFY --

if not isinstance(inputs_paths, list):
    raise ValueError("\x1B[31m" + "'inputs_paths' should be a list of 'str'" + "\x1B[0m")
for input_path in inputs_paths:
    if not isinstance(input_path, str):
        raise ValueError("\x1B[31m" + "'inputs_paths' should be a list of 'str'" + "\x1B[0m")

if not isinstance(outputs_path, str):
    raise ValueError("\x1B[31m" + "'outputs_path' should be a 'str'" + "\x1B[0m")

###################################################################################################

if not isinstance(origin, dict):
    raise ValueError("\x1B[31m" + "'origin' should be a 'dict'" + "\x1B[0m")
for key, value in origin.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'origin' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, (int, float, argmax)):
        raise ValueError("\x1B[31m" + "'origin' values should be 'int/float/argmax'" + "\x1B[0m")
if not isinstance(origin_units, dict):
    raise ValueError("\x1B[31m" + "'origin_units' should be a 'dict'" + "\x1B[0m")
for key, value in origin_units.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'origin_units' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, str):
        raise ValueError("\x1B[31m" + "'origin_units' values should be 'str'" + "\x1B[0m")

if not isinstance(zooms, dict):
    raise ValueError("\x1B[31m" + "'zooms' should be a 'dict'" + "\x1B[0m")
for key, value in zooms.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'zooms' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, tuple) or not len(value) == 2:
        raise ValueError("\x1B[31m" + "'zooms' values should be 'int/float'-doublets" + "\x1B[0m")
    if not isinstance(value[0], (int, float)) or not isinstance(value[1], (int, float)):
        raise ValueError("\x1B[31m" + "'zooms' values should be 'int/float'-doublets" + "\x1B[0m")
if not isinstance(zooms_units, dict):
    raise ValueError("\x1B[31m" + "'zooms_units' should be a 'dict'" + "\x1B[0m")
for key, value in zooms_units.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'zooms_units' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, str):
        raise ValueError("\x1B[31m" + "'zooms_units' values should be 'str'" + "\x1B[0m")
    
###################################################################################################

if not isinstance(outliers_sigmas, dict):
    raise ValueError("\x1B[31m" + "'outliers_sigmas' should be a 'dict'" + "\x1B[0m")
for key, value in outliers_sigmas.items():
    if not key or not isinstance(key, str):
        error_message = "'outliers_sigmas' keys should be un-empty 'str'"
        raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
    if not isinstance(value, (int, float)) or not value >= 0:
        error_message = "'outliers_sigmas' values should be positive 'int/float'"
        raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

if not isinstance(outliers_mask, dict):
    raise ValueError("\x1B[31m" + "'outliers_mask' should be a 'dict'" + "\x1B[0m")
for key, value in outliers_mask.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'outliers_mask' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, list):
        raise ValueError("\x1B[31m" + "'outliers_mask' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not isinstance(entry, tuple) or not len(entry) == 2:
            error_message = "'outliers_mask' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(entry[0], (int, float)) or not isinstance(entry[1], (int, float)):
            error_message = "'outliers_mask' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

###################################################################################################

if not isinstance(stellar_level, (int, float)) and not stellar_level >= 0:
    raise ValueError("\x1B[31m" + "'stellar_level' should be a postive 'int/float'" + "\x1B[0m")

if not isinstance(stellar_mask, dict):
    raise ValueError("\x1B[31m" + "'stellar_mask' should be a 'dict'" + "\x1B[0m")
for key, value in stellar_mask.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'stellar_mask' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, list):
        raise ValueError("\x1B[31m" + "'stellar_mask' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not isinstance(entry, tuple) or not len(entry) == 2:
            error_message = "'stellar_mask' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(entry[0], (int, float)) or not isinstance(entry[1], (int, float)):
            error_message = "'stellar_mask' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        
###################################################################################################

if not isinstance(psf_degrees, list):
    raise ValueError("\x1B[31m" + "'psf_degrees' should be a list" + "\x1B[0m")
for entry in psf_degrees:
    if not isinstance(entry, int) or not entry >= 0:
        raise ValueError("\x1B[31m" + "'psf_degrees' entries should be positive 'int'" + "\x1B[0m")

if not isinstance(psf_mask, dict):
    raise ValueError("\x1B[31m" + "'psf_mask' should be a 'dict'" + "\x1B[0m")
for key, value in psf_mask.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'psf_mask' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, list):
        raise ValueError("\x1B[31m" + "'psf_mask' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not isinstance(entry, tuple) or not len(entry) == 2:
            error_message = "'psf_mask' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(entry[0], (int, float)) or not isinstance(entry[1], (int, float)):
            error_message = "'psf_mask' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

###################################################################################################

if not isinstance(noise_loopback, int) or not noise_loopback >= 0:
    raise ValueError("\x1B[31m" + "'noise_loopback' should be a positive 'int'" + "\x1B[0m")

if not isinstance(noise_mask, dict):
    raise ValueError("\x1B[31m" + "'noise_mask' should be a 'dict'" + "\x1B[0m")
for key, value in noise_mask.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'noise_mask' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, list):
        raise ValueError("\x1B[31m" + "'noise_mask' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not isinstance(entry, tuple) or not len(entry) == 2:
            error_message = "'noise_mask' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(entry[0], (int, float)) or not isinstance(entry[1], (int, float)):
            error_message = "'noise_mask' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

###################################################################################################

if not isinstance(anomalies_sigmas, dict):
    raise ValueError("\x1B[31m" + "'anomalies_sigmas' should be a 'dict'" + "\x1B[0m")
for key, value in anomalies_sigmas.items():
    if not key or not isinstance(key, str):
        error_message = "'anomalies_sigmas' keys should be un-empty 'str'"
        raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
    if not isinstance(value, (int, float)) or not value >= 0:
        error_message = "'anomalies_sigmas' values should be positive 'int/float'"
        raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

if not isinstance(anomalies_mask, dict):
    raise ValueError("\x1B[31m" + "'anomalies_mask' should be a 'dict'" + "\x1B[0m")
for key, value in anomalies_mask.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'anomalies_mask' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, list):
        raise ValueError("\x1B[31m" + "'anomalies_mask' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not isinstance(entry, tuple) or not len(entry) == 2:
            error_message = "'anomalies_mask' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(entry[0], (int, float)) or not isinstance(entry[1], (int, float)):
            error_message = "'anomalies_mask' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

###################################################################################################

if not save_verbose in [0,1]:
    raise ValueError("\x1B[31m" + "'save_verbose' should be in {'0','1'}" + "\x1B[0m")
if not print_verbose in [0,1]:
    raise ValueError("\x1B[31m" + "'print_verbose' should be in {'0','1'}" + "\x1B[0m")
if not show_verbose in [0,1]:
    raise ValueError("\x1B[31m" + "'show_verbose' should be in {'0','1'}" + "\x1B[0m")

if not isinstance(spectra_filter, dict):
    raise ValueError("\x1B[31m" + "'spectra_filter' should be a 'dict'" + "\x1B[0m")
for key, value in spectra_filter.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'spectra_filter' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, list):
        raise ValueError("\x1B[31m" + "'spectra_filter' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not isinstance(entry, tuple) or not len(entry) == 2:
            error_message = "'spectra_filter' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(entry[0], (int, float)) or not isinstance(entry[1], (int, float)):
            error_message = "'spectra_filter' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
if not isinstance(images_filter, dict):
    raise ValueError("\x1B[31m" + "'images_filter' should be a 'dict'" + "\x1B[0m")
for key, value in images_filter.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'images_filter' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, list):
        raise ValueError("\x1B[31m" + "'images_filter' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not isinstance(entry, tuple) or not len(entry) == 2:
            error_message = "'images_filter' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(entry[0], (int, float)) or not isinstance(entry[1], (int, float)):
            error_message = "'images_filter' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

if not isinstance(spectra_masker, dict):
    raise ValueError("\x1B[31m" + "'spectra_masker' should be a 'dict'" + "\x1B[0m")
for key, value in spectra_masker.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'spectra_masker' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, list):
        raise ValueError("\x1B[31m" + "'spectra_masker' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not isinstance(entry, tuple) or not len(entry) == 2:
            error_message = "'spectra_masker' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(entry[0], (int, float)) or not isinstance(entry[1], (int, float)):
            error_message = "'spectra_masker' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
if not isinstance(images_masker, dict):
    raise ValueError("\x1B[31m" + "'images_masker' should be a 'dict'" + "\x1B[0m")
for key, value in images_masker.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'images_masker' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, list):
        raise ValueError("\x1B[31m" + "'images_masker' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not isinstance(entry, tuple) or not len(entry) == 2:
            error_message = "'images_masker' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(entry[0], (int, float)) or not isinstance(entry[1], (int, float)):
            error_message = "'images_masker' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

if not isinstance(spectra_hinder, dict):
    raise ValueError("\x1B[31m" + "'spectra_hinder' should be a 'dict'" + "\x1B[0m")
for key, value in spectra_hinder.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'spectra_hinder' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, list):
        raise ValueError("\x1B[31m" + "'spectra_hinder' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not isinstance(entry, tuple) or not len(entry) == 2:
            error_message = "'spectra_hinder' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(entry[0], (int, float)) or not isinstance(entry[1], (int, float)):
            error_message = "'spectra_hinder' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
if not isinstance(images_hinder, dict):
    raise ValueError("\x1B[31m" + "'images_hinder' should be a 'dict'" + "\x1B[0m")
for key, value in images_hinder.items():
    if not key or not isinstance(key, str):
        raise ValueError("\x1B[31m" + "'images_hinder' keys should be un-empty 'str'" + "\x1B[0m")
    if not isinstance(value, list):
        raise ValueError("\x1B[31m" + "'images_hinder' values should be 'list" + "\x1B[0m")
    for entry in value:
        if not isinstance(entry, tuple) or not len(entry) == 2:
            error_message = "'images_hinder' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(entry[0], (int, float)) or not isinstance(entry[1], (int, float)):
            error_message = "'images_hinder' values entries should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
