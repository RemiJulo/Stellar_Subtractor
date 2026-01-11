###################
# PARAMETERS FILE #
###################

class argmax: pass

################################################################################################### ###################
################################################################################################### IN/OUTPUTS FEATURES
################################################################################################### ###################
# List of paths to .fits files (or to folders full of .fits files) to be post-processed

inputs_paths = ["inputs/"]

###################################################################################################
# Path to the saving folder from the location of the code files
# Any needed series of folders will be created in the operation
# (not ending with '/' or '\' allows to add a prefix)

outputs_path = "outputs/"

################################################################################################### ###################
################################################################################################### AXIS GRIDS FEATURES
################################################################################################### ###################
# Origin coordinates relocation ('int' for pixel values or 'float' for physical values)
# Herewith, 'argmax()' is also usable to recover the coordinates of the brightest pixel
# For instance, the entry " 'AWAV' : 0.0 " allows to replace the air wavelengths origin
# from the air wavelength of the first spectral channel to the one of ~ 0.0 m
# Built as a dictionary associating axis names with origin coordinates
# It is also possible to change the units used for the physical values
# (the list of compatible units is printed in case of bad unit values)

origin = {



}

origin_units = {



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



}

zooms_units = {



}

################################################################################################### ###################
################################################################################################### ESTIMATION FEATURES
################################################################################################### ###################
# Number of standard deviations above the mean beyond which a value is considered to be an outlier

outliers_sigma = 5

###################################################################################################
# Misleading a priori in the form of a list of dictionaries associating axis names with doublets
# between which pixels are masked (thus disregarded) during the outliers detection and selection
# Each dictionnary defines a mask as an intersection of intervals along the different given axes

outliers_mask = [


    
]

###################################################################################################
# Height of the stellar level compared to the one of the noise level (setting noise/stellar areas)
# The stellar component of the pixels with a flux below this threshold is assumed to be negligible
# The corresponding pixels are not used for estimating the stellar spectrum by spatial integration
# The unused pixels are masked in the image of the stellar spectrum estimate if 'show_verbose' = 1
# This helps to limit electronic noise contamination, especially with large selected field of view
# Too large values of 'stellar_sigma' lead to a selection of too many too noisy pixels in practice
# (thus could lead to uncanny results in the following because of various numerical instabilities)

star_fluxlevel = 50

###################################################################################################
# Misleading a priori in the form of a list of dictionaries associating axis names with doublets
# between which pixels are masked (thus disregarded) during the stellar spectrum estimation step
# Each dictionnary defines a mask as an intersection of intervals along the different given axes

star_mask = [



]

###################################################################################################
# Instrument-dependent hyperparameter dealing with the smoothness of its point spread function
# (knowing that this smoothness is relative to the spectral window sizes ultimately selected)
# It is the degree of truncation of the Legendre polynomial basis used for the decomposition
# of the chromatic PSF (i.e. the spectrum deformations given by the PSF spectral evolution)
# This can also be seen as the degree of polynomials used to fit each of the chromatic PSF
# The appendix file 'app.py' allows to validate the relevance of this choice afterward
# by verifying the decrease in the energy of the coefficients down to the noise levels
# as well as the shape of the patterns actually estimated within the stellar component

psf_polydegree = 5

###################################################################################################
# Misleading a priori in the form of a list of dictionaries associating axis names with doublets
# between which pixels are masked (thus disregarded) during the point spread function estimation
# Each dictionnary defines a mask as an intersection of intervals along the different given axes

psf_mask = [



]

###################################################################################################
# Number of PSF estimation loopback (taking into account an increasingly better noise estimate)

noise_loopback = 0

###################################################################################################
# Misleading a priori in the form of a list of dictionaries associating axis names with doublets
# between which pixels are masked (thus disregarded) during the 2 covariance matrices estimation
# Each dictionnary defines a mask as an intersection of intervals along the different given axes

noise_mask = [



]

###################################################################################################
# Number of standard deviations above the mean beyond which a value is considered to be an anomaly

anomalies_sigma = 0

###################################################################################################
# Misleading a priori in the form of a list of dictionaries associating axis names with doublets
# between which pixels are masked (thus disregarded) during the anomalies detection and highight
# Each dictionnary defines a mask as an intersection of intervals along the different given axes

anomalies_mask = [



]

################################################################################################### ###################
################################################################################################### MONITORING FEATURES
################################################################################################### ###################
# 'save_verbose' for backup feedback in 'outputs_path' ('0' / '1' for no / full display)
# 'print_verbose' for textual feedback in the terminal ('0' / '1' for no / full display)
# 'show_verbose' for visual feedback in python windows ('0' / '1' for no / full display)

save_verbose = 1

print_verbose = 1

show_verbose = 1

###################################################################################################
# Data filters for visual feedbacks (as a product from unit-height Gaussians)
# Built as lists of dictionaries associating axis names with (location, scale)
# Gaussians parameters with which to multiply images and spectra to filter them

spectra_filter = [



]

images_filter = [



]

###################################################################################################
# Data maskers for visual feedbacks (as NaNing of the intersection of the selected intervals)
# Built as lists of dictionaries associating axis names with lower and upper bounds doublets
# whose intersection of the corresponding intervals is masked for each of the dictionnaries

spectra_masker = [



]

images_masker = [



]

###################################################################################################
# Data hinders for visual feedbacks (as a product from unit-height Gaussians)
# Built as lists of dictionaries associating axis names with (location, scale)
# Gaussians parameters with which to multiply images and spectra to hinder them

spectra_hinder = [



]

images_hinder = [



]

################################################################################################### ###################
################################################################################################### PARAMETERS CHECKING
################################################################################################### -- DO NOT MODIFY --
################################################################################################### ###################

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

if not isinstance(outliers_sigma, (int, float)) or not outliers_sigma >= 0:
    error_message = "'outliers_sigma' should be a positive 'int/float'"
    raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

if not isinstance(outliers_mask, list):
    raise ValueError("\x1B[31m" + "'outliers_mask' should be a 'list'" + "\x1B[0m")
for entry in outliers_mask:
    if not isinstance(entry, dict):
        raise ValueError("\x1B[31m" + "'outliers_mask' entries should be 'dict'" + "\x1B[0m")
    for key, value in entry.items():
        if not key or not isinstance(key, str):
            error_message = "'outliers_mask' entries keys should be un-empty 'str'"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value, tuple) or not len(value) == 2:
            error_message = "'outliers_mask' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value[0], (int, float)) or not isinstance(value[1], (int, float)):
            error_message = "'outliers_mask' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

if not isinstance(star_fluxlevel, (int, float)) and not star_fluxlevel >= 0:
    raise ValueError("\x1B[31m" + "'star_fluxlevel' should be a postive 'int/float'" + "\x1B[0m")

if not isinstance(star_mask, list):
    raise ValueError("\x1B[31m" + "'star_mask' should be a 'list'" + "\x1B[0m")
for entry in star_mask:
    if not isinstance(entry, dict):
        raise ValueError("\x1B[31m" + "'star_mask' entries should be 'dict'" + "\x1B[0m")
    for key, value in entry.items():
        if not key or not isinstance(key, str):
            error_message = "'star_mask' entries keys should be un-empty 'str'"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value, tuple) or not len(value) == 2:
            error_message = "'star_mask' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value[0], (int, float)) or not isinstance(value[1], (int, float)):
            error_message = "'star_mask' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

if not isinstance(psf_polydegree, int) or not psf_polydegree >= 0:
    raise ValueError("\x1B[31m" + "'psf_polydegree' should be a positive 'int'" + "\x1B[0m")

if not isinstance(psf_mask, list):
    raise ValueError("\x1B[31m" + "'psf_mask' should be a 'list'" + "\x1B[0m")
for entry in psf_mask:
    if not isinstance(entry, dict):
        raise ValueError("\x1B[31m" + "'psf_mask' entries should be 'dict'" + "\x1B[0m")
    for key, value in entry.items():
        if not key or not isinstance(key, str):
            error_message = "'psf_mask' entries keys should be un-empty 'str'"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value, tuple) or not len(value) == 2:
            error_message = "'psf_mask' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value[0], (int, float)) or not isinstance(value[1], (int, float)):
            error_message = "'psf_mask' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        
if not isinstance(noise_loopback, int) or not noise_loopback >= 0:
    raise ValueError("\x1B[31m" + "'noise_loopback' should be a positive 'int'" + "\x1B[0m")

if not isinstance(noise_mask, list):
    raise ValueError("\x1B[31m" + "'noise_mask' should be a 'list'" + "\x1B[0m")
for entry in noise_mask:
    if not isinstance(entry, dict):
        raise ValueError("\x1B[31m" + "'noise_mask' entries should be 'dict'" + "\x1B[0m")
    for key, value in entry.items():
        if not key or not isinstance(key, str):
            error_message = "'noise_mask' entries keys should be un-empty 'str'"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value, tuple) or not len(value) == 2:
            error_message = "'noise_mask' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value[0], (int, float)) or not isinstance(value[1], (int, float)):
            error_message = "'noise_mask' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

if not isinstance(anomalies_sigma, (int, float)) or not anomalies_sigma >= 0:
    error_message = "'anomalies_sigma' should be a positive 'int/float'"
    raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

if not isinstance(anomalies_mask, list):
    raise ValueError("\x1B[31m" + "'anomalies_mask' should be a 'list'" + "\x1B[0m")
for entry in anomalies_mask:
    if not isinstance(entry, dict):
        raise ValueError("\x1B[31m" + "'anomalies_mask' entries should be 'dict'" + "\x1B[0m")
    for key, value in entry.items():
        if not key or not isinstance(key, str):
            error_message = "'anomalies_mask' entries keys should be un-empty 'str'"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value, tuple) or not len(value) == 2:
            error_message = "'anomalies_mask' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value[0], (int, float)) or not isinstance(value[1], (int, float)):
            error_message = "'anomalies_mask' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        
###################################################################################################

if not save_verbose in [0,1]:
    raise ValueError("\x1B[31m" + "'save_verbose' should be in {'0','1'}" + "\x1B[0m")
if not print_verbose in [0,1]:
    raise ValueError("\x1B[31m" + "'print_verbose' should be in {'0','1'}" + "\x1B[0m")
if not show_verbose in [0,1]:
    raise ValueError("\x1B[31m" + "'show_verbose' should be in {'0','1'}" + "\x1B[0m")

if not isinstance(spectra_filter, list):
    raise ValueError("\x1B[31m" + "'spectra_filter' should be a 'list'" + "\x1B[0m")
for entry in spectra_filter:
    if not isinstance(entry, dict):
        raise ValueError("\x1B[31m" + "'spectra_filter' entries should be 'dict'" + "\x1B[0m")
    for key, value in entry.items():
        if not key or not isinstance(key, str):
            error_message = "'spectra_filter' entries keys should be un-empty 'str'"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value, tuple) or not len(value) == 2:
            error_message = "'spectra_filter' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value[0], (int, float)) or not isinstance(value[1], (int, float)):
            error_message = "'spectra_filter' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
if not isinstance(images_filter, list):
    raise ValueError("\x1B[31m" + "'images_filter' should be a 'list'" + "\x1B[0m")
for entry in images_filter:
    if not isinstance(entry, dict):
        raise ValueError("\x1B[31m" + "'images_filter' entries should be 'dict'" + "\x1B[0m")
    for key, value in entry.items():
        if not key or not isinstance(key, str):
            error_message = "'images_filter' entries keys should be un-empty 'str'"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value, tuple) or not len(value) == 2:
            error_message = "'images_filter' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value[0], (int, float)) or not isinstance(value[1], (int, float)):
            error_message = "'images_filter' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        
if not isinstance(spectra_masker, list):
    raise ValueError("\x1B[31m" + "'spectra_masker' should be a 'list'" + "\x1B[0m")
for entry in spectra_masker:
    if not isinstance(entry, dict):
        raise ValueError("\x1B[31m" + "'spectra_masker' entries should be 'dict'" + "\x1B[0m")
    for key, value in entry.items():
        if not key or not isinstance(key, str):
            error_message = "'spectra_masker' entries keys should be un-empty 'str'"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value, tuple) or not len(value) == 2:
            error_message = "'spectra_masker' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value[0], (int, float)) or not isinstance(value[1], (int, float)):
            error_message = "'spectra_masker' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
if not isinstance(images_masker, list):
    raise ValueError("\x1B[31m" + "'images_masker' should be a 'list'" + "\x1B[0m")
for entry in images_masker:
    if not isinstance(entry, dict):
        raise ValueError("\x1B[31m" + "'images_masker' entries should be 'dict'" + "\x1B[0m")
    for key, value in entry.items():
        if not key or not isinstance(key, str):
            error_message = "'images_masker' entries keys should be un-empty 'str'"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value, tuple) or not len(value) == 2:
            error_message = "'images_masker' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value[0], (int, float)) or not isinstance(value[1], (int, float)):
            error_message = "'images_masker' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        
if not isinstance(spectra_hinder, list):
    raise ValueError("\x1B[31m" + "'spectra_hinder' should be a 'list'" + "\x1B[0m")
for entry in spectra_hinder:
    if not isinstance(entry, dict):
        raise ValueError("\x1B[31m" + "'spectra_hinder' entries should be 'dict'" + "\x1B[0m")
    for key, value in entry.items():
        if not key or not isinstance(key, str):
            error_message = "'spectra_hinder' entries keys should be un-empty 'str'"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value, tuple) or not len(value) == 2:
            error_message = "'spectra_hinder' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value[0], (int, float)) or not isinstance(value[1], (int, float)):
            error_message = "'spectra_hinder' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
if not isinstance(images_hinder, list):
    raise ValueError("\x1B[31m" + "'images_hinder' should be a 'list'" + "\x1B[0m")
for entry in images_hinder:
    if not isinstance(entry, dict):
        raise ValueError("\x1B[31m" + "'images_hinder' entries should be 'dict'" + "\x1B[0m")
    for key, value in entry.items():
        if not key or not isinstance(key, str):
            error_message = "'images_hinder' entries keys should be un-empty 'str'"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value, tuple) or not len(value) == 2:
            error_message = "'images_hinder' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
        if not isinstance(value[0], (int, float)) or not isinstance(value[1], (int, float)):
            error_message = "'images_hinder' entries values should be 'int/float'-doublets"
            raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
