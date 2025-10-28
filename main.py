import os
import time

import warnings

import numpy as np

from astropy.io import fits
from astropy.wcs import WCS

from astropy import units as u

import param
import process

###################################################################################################

from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category = AstropyWarning)

sep = os.path.sep # '/' for Linux/Mac or '\' for Windows

################################################################################################### Time starting

start = time.time()

################################################################################################### Paths expanding

expanded_input_paths = [] # To expand paths to folders in list of paths to files

if not type(param.input_paths) is list:
    raise ValueError("\x1B[31m" + "'input_paths' should be a list of 'str'" + "\x1B[0m")

for input_path in param.input_paths:

    if not type(input_path) is str:
        raise ValueError("\x1B[31m" + "'input_paths' should be a list of 'str'" + "\x1B[0m")

    input_path = input_path.replace('/', sep).replace('\\', sep)

    if not os.path.exists(input_path):
        raise FileNotFoundError("\x1B[31m" + input_path + "\x1B[0m")

    if not '.' in input_path:
        for input_file in os.listdir(input_path):
            input_file_path = input_path + sep + input_file
            expanded_input_paths.append(input_file_path)
    else : expanded_input_paths.append(input_path)

################################################################################################### Parameters checking

if not param.output_path or not type(param.output_path) is str:
    raise ValueError("\x1B[31m" + "'output_path' should be a un-empty 'str'" + "\x1B[0m")

if not type(param.origin) is dict:
    raise ValueError("\x1B[31m" + "'origin' should be a 'dict'" + "\x1B[0m")
for key, value in param.origin.items():
    if not key or not type(key) is str:
        raise ValueError("\x1B[31m" + "'origin' keys should be un-empty 'str'" + "\x1B[0m")
    if not value is param.argmax and not type(value) in (int, float):
        raise ValueError("\x1B[31m" + "'origin' values should be 'int/float/arg'" + "\x1B[0m")
if not type(param.origin_units) is dict:
    raise ValueError("\x1B[31m" + "'origin_units' should be a 'dict'" + "\x1B[0m")
for key, value in param.origin_units.items():
    if not key or not type(key) is str:
        raise ValueError("\x1B[31m" + "'origin_units' keys should be un-empty 'str'" + "\x1B[0m")
    if not type(value) is str:
        raise ValueError("\x1B[31m" + "'origin_units' values should be 'str'" + "\x1B[0m")

if not type(param.zooms) is dict:
    raise ValueError("\x1B[31m" + "'zooms' should be a 'dict'" + "\x1B[0m")
for key, value in param.zooms.items():
    if not key or not type(key) is str:
        raise ValueError("\x1B[31m" + "'zooms' keys should be un-empty 'str'" + "\x1B[0m")
    if not type(value) is tuple or not len(value) == 2:
        raise ValueError("\x1B[31m" + "'zooms' values should be 'int/float'-doublets" + "\x1B[0m")
    if not type(value[0]) in (int, float) or not type(value[1]) in (int, float):
        raise ValueError("\x1B[31m" + "'zooms' values should be 'int/float'-doublets" + "\x1B[0m")
if not type(param.zooms_units) is dict:
    raise ValueError("\x1B[31m" + "'zooms_units' should be a 'dict'" + "\x1B[0m")
for key, value in param.zooms_units.items():
    if not key or not type(key) is str:
        raise ValueError("\x1B[31m" + "'zooms_units' keys should be un-empty 'str'" + "\x1B[0m")
    if not type(value) is str:
        raise ValueError("\x1B[31m" + "'zooms_units' values should be 'str'" + "\x1B[0m")

if not param.save_verbose in [0,1]:
    raise ValueError("\x1B[31m" + "'save_verbose' should be in {'0','1'}" + "\x1B[0m")
if not param.print_verbose in [0,1]:
    raise ValueError("\x1B[31m" + "'print_verbose' should be in {'0','1'}" + "\x1B[0m")
if not param.show_verbose in [0,1]:
    raise ValueError("\x1B[31m" + "'show_verbose' should be in {'0','1'}" + "\x1B[0m")

if not type(param.spectra_filter) is dict:
    raise ValueError("\x1B[31m" + "'spectra_filter' should be a 'dict'" + "\x1B[0m")
for key, value in param.spectra_filter.items():
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
if not type(param.images_filter) is dict:
    raise ValueError("\x1B[31m" + "'images_filter' should be a 'dict'" + "\x1B[0m")
for key, value in param.images_filter.items():
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

if not type(param.spectra_blocker) is dict:
    raise ValueError("\x1B[31m" + "'spectra_blocker' should be a 'dict'" + "\x1B[0m")
for key, value in param.spectra_blocker.items():
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
if not type(param.images_blocker) is dict:
    raise ValueError("\x1B[31m" + "'images_blocker' should be a 'dict'" + "\x1B[0m")
for key, value in param.images_blocker.items():
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

if not type(param.stellar_sigma) in (int, float):
    raise ValueError("\x1B[31m" + "'stellar_sigma' should be 'int/float'" + "\x1B[0m")

if not type(param.stellar_mask) is dict:
    raise ValueError("\x1B[31m" + "'stellar_mask' should be a 'dict'" + "\x1B[0m")
for key, value in param.stellar_mask.items():
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

if not type(param.psf_degrees) is list:
    raise ValueError("\x1B[31m" + "'psf_degrees' should be a list" + "\x1B[0m")
for entry in param.psf_degrees:
    if not type(entry) is int or not entry > 0:
        raise ValueError("\x1B[31m" + "'psf_degrees' entries should be a 'int' > 0" + "\x1B[0m")

if not type(param.psf_mask) is dict:
    raise ValueError("\x1B[31m" + "'psf_mask' should be a 'dict'" + "\x1B[0m")
for key, value in param.psf_mask.items():
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

################################################################################################### Data files opening

for input_path in expanded_input_paths:
    if not input_path.endswith(".fits"):
        continue # .fits files only

    input_path = input_path.replace('/', sep).replace('\\', sep)

    print("\x1B[1;34m" + f"\nFile opening: {input_path.split(sep)[-1]}\n" + "\x1B[0m")

    hdul = fits.open(input_path)

    for hdu in hdul:
        data, header = hdu.data, hdu.header

################################################################################################### Input preparing

        # WCS stands for 'World Coordinate System'
        if type(data) is np.ndarray:
            wcs = WCS(header)

            # Axis coordinate types (among 'spectral', 'celestial', 'temporal', ...)
            coordinate_types = [axis_type["coordinate_type"] for axis_type in wcs.get_axis_types()]

            # Existence and uniqueness of 1 spectral axis
            if coordinate_types.count('spectral') != 1: break
            spectral_axis_index = coordinate_types.index('spectral')

            # Axes names and units
            axis_names = wcs.axis_type_names
            axis_units = wcs.world_axis_units

################################################################################################### Origin replacement

            # Default origin from the header
            array_origin = np.rint(wcs.wcs.crpix-1).astype(int)
            world_origin = np.array(wcs.wcs.crval).astype(float)
            
            # 'argmax' indexes save
            argmax_indexes_list = []

            # Origin coordinates replacement
            # (ignoring the non-existent axes)
            for key, value in param.origin.items():
                if key not in axis_names: continue

                # Replacement of the associated origin coordinate with 'int' values
                # (and update of all other coordinates to deal with correlations between axes)
                if type(value) is int:
                    array_origin[axis_names.index(key)] = value
                    [world_origin] = np.array(wcs.array_index_to_world_values(*[[array_origin]]))
                    if key in param.origin_units.keys():
                        irreducible_unit = u.Unit(axis_units[axis_names.index(key)])
                        try: # Unit conversion if possible
                            target_unit = u.Unit(param.origin_units[key])
                            world_origin = (world_origin * irreducible_unit).to(target_unit).value
                        except (ValueError, u.UnitConversionError):
                            error_message = f"Axis {key} unit should be in:\n"
                            error_message += f"{irreducible_unit.find_equivalent_units()}"
                            raise ValueError("\x1B[31m" + error_message + "\x1B[0m") from None
                
                # Replacement of the associated origin coordinate with 'float' values
                # (and update of all other coordinates to deal with correlations between axes)
                if type(value) is float:
                    if key in param.origin_units.keys():
                        irreducible_unit = u.Unit(axis_units[axis_names.index(key)])
                        try: # Unit conversion if possible
                            target_unit = u.Unit(param.origin_units[key])
                            value = (value * target_unit).to(irreducible_unit).value
                        except (ValueError, u.UnitConversionError):
                            error_message = f"Axis {key} unit should be in:\n"
                            error_message += f"{irreducible_unit.find_equivalent_units()}"
                            raise ValueError("\x1B[31m" + error_message + "\x1B[0m") from None
                    world_origin[axis_names.index(key)] = value
                    [array_origin] = np.array(wcs.world_to_array_index_values(*[[world_origin]]))

                # Subsquent assignation of the 'argmax' values
                # to ensure that the final values are indeed 'argmax'
                # (despite the correlations with the potential following axes)
                if value is param.argmax:
                    argmax_indexes_list.append(axis_names.index(key))

            # Preliminary search for the minimum and maximum flux pixels, if necessary only
            # (to avoid the redundancy of the same costly calculation with multiple 'argmax' uses)
            if argmax_indexes_list:
                if param.print_verbose:
                    print("    Search for the maximum flux pixel", end = '\r')

                # Search for the maximum flux pixel
                argmax_array = np.unravel_index(np.nanargmax(data), data.shape)[::-1]
            
                # Origin coordinates replacement
                for index in argmax_indexes_list:
                    array_origin[index] = argmax_array[index]

                # Update of all other coordinates to deal with correlations between axes
                [world_origin] = np.array(wcs.array_index_to_world_values(*[[array_origin]]))

################################################################################################### Data zooming

            # Slice, world and array zooms initialization
            array_zooms = np.array([[0]*len(data.shape[::-1]), list(data.shape[::-1])])
            world_zooms = np.array(wcs.array_index_to_world_values(*[array_zooms]))
            slice_zooms = [slice(0, axis_len) for axis_len in wcs.pixel_shape]

            # Data zooming from the origin
            # (ignoring the non-existent axes)
            for key, value in param.zooms.items():
                if key not in axis_names: continue

                # Zoom lower and upper bounds
                lower_bound, upper_bound = value

                # Replacement of the associated zoom lower and upper bounds with 'int' values
                if type(lower_bound) is int:
                    array_origin_key = array_origin[axis_names.index(key)]
                    array_zooms[0][axis_names.index(key)] = array_origin_key + lower_bound
                if type(upper_bound) is int:
                    origin_key = array_origin[axis_names.index(key)]
                    array_zooms[1][axis_names.index(key)] = array_origin_key + upper_bound

                # Replacement of the associated zoom lower and upper bounds with 'float' values
                if type(lower_bound) is float:
                    if key in param.zooms_units.keys():
                        irreducible_unit = u.Unit(axis_units[axis_names.index(key)])
                        try: # Unit conversion if possible
                            target_unit = u.Unit(param.zooms_units[key])
                            lower_bound = (lower_bound * target_unit).to(irreducible_unit).value
                        except (ValueError, u.UnitConversionError):
                            error_message = f"Axis {key} unit should be in:\n"
                            error_message += f"{irreducible_unit.find_equivalent_units()}"
                            raise ValueError("\x1B[31m" + error_message + "\x1B[0m") from None
                    world_origin_key = world_origin[axis_names.index(key)]
                    world_zooms[0][axis_names.index(key)] = world_origin_key + lower_bound
                if type(upper_bound) is float:
                    if key in param.zooms_units.keys():
                        irreducible_unit = u.Unit(axis_units[axis_names.index(key)])
                        try: # Unit conversion if possible
                            target_unit = u.Unit(param.zooms_units[key])
                            upper_bound = (upper_bound * target_unit).to(irreducible_unit).value
                        except (ValueError, u.UnitConversionError):
                            error_message = f"Axis {key} unit should be in:\n"
                            error_message += f"{irreducible_unit.find_equivalent_units()}"
                            raise ValueError("\x1B[31m" + error_message + "\x1B[0m") from None
                    world_origin_key = world_origin[axis_names.index(key)]
                    world_zooms[1][axis_names.index(key)] = world_origin_key + upper_bound

            # Add of the world zooms to the array zooms
            # The data shape is also subtracted to the array zooms
            # (to deal with the duplicate caused by double array-world initialization)
            world_zooms_to_array_zooms = wcs.world_to_array_index_values(*[world_zooms])
            array_zooms += np.array(world_zooms_to_array_zooms)[::-1]
            array_zooms[1] -= np.array(data.shape[::-1])
            array_bounds = list(zip(*array_zooms))

            # Zoom validity checking (to avoid errors caused by flat axes)
            for axis, (lower_bound, upper_bound) in enumerate(array_bounds):
                if lower_bound > upper_bound: # Decreasing axis cases
                    lower_bound, upper_bound = upper_bound-1, lower_bound+1
                    array_bounds[axis] = (lower_bound, upper_bound)
                if lower_bound > wcs.pixel_shape[axis]: # Overflow beyond the number of pixels
                    error_message = f"Lower bound too large for the '{axis_names[axis]}' axis"
                    raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
                if upper_bound < 0 : # Overflow beyond the index zero origin of the array axes
                    error_message = f"Upper bound too small for the '{axis_names[axis]}' axis"
                    raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

            # Zoom arrays to slices conversion
            for key, value in param.zooms.items():
                if key not in axis_names: continue
                list_array_zooms_key = list(array_bounds[axis_names.index(key)])
                list_array_zooms_key[0] = max(0, list_array_zooms_key[0]) # For positive indexes
                list_array_zooms_key[1] += 1 # To include the upper bounds in the zoom selection
                slice_zooms[axis_names.index(key)] = slice(*list_array_zooms_key)

            # Data and WCS zooming
            tuple_slice_zooms = tuple(slice_zooms)[::-1]
            zoomed_wcs = wcs.slice(tuple_slice_zooms)
            zoomed_data = data[tuple_slice_zooms]

################################################################################################### Grids sum up

            if param.print_verbose:

                world_min_max = wcs.array_index_to_world_values(*[array_zooms])
                [world_origin] = wcs.array_index_to_world_values(*[[array_origin]])

                axes_information = [axis_names, axis_units]
                axes_information += [array_origin, world_origin]
                axes_information += list(array_zooms) + list(world_min_max)

                i = zip(*tuple(axes_information)) # Information list for all axes

                for name, unit, array_0, world_0, array_min, array_max, world_min, world_max in i:

                    if name in param.origin_units.keys():
                        irreducible_unit = u.Unit(axis_units[axis_names.index(name)])
                        try: # Unit conversion if possible
                            unit = u.Unit(param.origin_units[name])
                            world_0 = (world_0 * irreducible_unit).to(unit).value
                            world_min = (world_min * irreducible_unit).to(unit).value
                            world_max = (world_max * irreducible_unit).to(unit).value
                        except (ValueError, u.UnitConversionError):
                            error_message = f"Axis {key} unit should be in:\n"
                            error_message += f"{irreducible_unit.find_equivalent_units()}"
                            raise ValueError("\x1B[31m" + error_message + "\x1B[0m") from None

                    axis_information = " -> " + "\x1B[1;34m" + f"{name}" + "\x1B[0;34m"
                    axis_information += f" from {world_min} {unit} (pixel {array_min})"
                    axis_information += f" to {world_max} {unit} (pixel {array_max})"

                    axis_bounds = " with " + "\x1B[1;34m" +  "origin" +  "\x1B[0;34m"
                    axis_bounds += f" at {world_0} {unit} (pixel {array_0})"

                    print("\x1B[34m" + axis_information + axis_bounds + "\x1B[0m")

################################################################################################### World grids build

            # Correlation matrix between axes
            # to avoid unnecessary massive calculations
            correlation_matrix = zoomed_wcs.axis_correlation_matrix

            # Correlation groups for world grids build from needed axes only
            # (with shapes order inversion to deal with wcs-array reversed order)
            needed_axes = [np.where(correlation)[0] for correlation in correlation_matrix]
            needed_shapes = [np.array(zoomed_data.shape[::-1])[axis][::-1] for axis in needed_axes]

            # World grids completion
            # correlated-axes-by-correlated-axes
            world_grids = [[[], None]]*len(zoomed_data.shape[::-1])
            for axis, (corr_axes, corr_shapes) in enumerate(zip(needed_axes, needed_shapes)):

                # Optimization by avoiding the calculation of grids
                # already determined for axes with the same correlations
                if list(world_grids[axis][0]) == list(corr_axes): continue

                # World grids determination from 'zoomed_wcs' and array grids of the correlated axes
                needed_wcs = zoomed_wcs.sub(tuple(int(corr_axis) + 1 for corr_axis in corr_axes))
                needed_grids = needed_wcs.array_index_to_world_values(*np.indices(corr_shapes))

                # Save of the grid for the given axis
                # as well as of the grids for the correlated axes
                # (generally without additional correlations with other axes)
                for corr_axis, needed_grid in zip(corr_axes, np.atleast_2d(needed_grids)):
                    needed_grid -= world_origin[corr_axis] # Shift to the origin
                    if axis_names[corr_axis] in param.grids_units.keys():
                        irreducible_unit = u.Unit(axis_units[corr_axis])
                        try: # Unit conversion if possible
                            target_unit = u.Unit(param.grids_units[axis_names[corr_axis]])
                            needed_grid = (needed_grid * irreducible_unit).to(target_unit).value
                        except (ValueError, u.UnitConversionError):
                            error_message = f"Axis {key} unit should be in:\n"
                            error_message += f"{irreducible_unit.find_equivalent_units()}"
                            raise ValueError("\x1B[31m" + error_message + "\x1B[0m") from None
                    world_grids[corr_axis] = (list(corr_axes), needed_grid)

################################################################################################### Data matricization
            
            # Move of the spectral axis to the last position
            # Optimization by putting the observations along the rows
            # (with a transpose '.T' against the wcs-array reversed order)
            moved_data = np.moveaxis(zoomed_data.T, spectral_axis_index, -1)
            
            # Move of the spectral axis to the last position for grids too
            # (and associated correction of the grid axis indexes where needed)
            for axis_index in range(len(world_grids)):
                for corr_axis_index in range(len(world_grids[axis_index][0])):
                    if world_grids[axis_index][0][corr_axis_index] == spectral_axis_index:
                        world_grids[axis_index][0][corr_axis_index] = len(world_grids) - 1
                    elif world_grids[axis_index][0][corr_axis_index] > spectral_axis_index:
                        world_grids[axis_index][0][corr_axis_index] -= 1
            world_grids.append(world_grids[spectral_axis_index])
            world_grids.pop(spectral_axis_index)

            # Move of the spectral axis
            # to the last position for axis names
            axis_names.append(axis_names[spectral_axis_index])
            axis_names.pop(spectral_axis_index)

            # Move of the spectral axis
            # to the last position for axis units
            axis_units.append(axis_units[spectral_axis_index])
            axis_units.pop(spectral_axis_index)

            # Pixels axis
            # name and unit
            if param.show_verbose:
                try:
                    axis_names += [header['BTYPE']]
                except KeyError:
                    axis_names += ["Flux"]
                try:
                    axis_units += [header['BUNIT']]
                except KeyError:
                    axis_units += ["< u >"]
            
            # Vectorization of the first (non-spectral) axes
            # Optimization by putting the features along the columns
            matricized_data = moved_data.reshape(-1, moved_data.shape[-1])

################################################################################################### Main processing

            # Main processing from the observation matrix 'matricized_data'
            processed_data = process.process(matricized_data, world_grids, axis_names, axis_units)

################################################################################################### Data tensorization

            # Devectorization of the first (non-spectral) axes
            tensorized_data = processed_data.reshape(tuple(moved_data.shape))

            # Move of the spectral axis to its original pre-move position
            # (with a transpose '.T' against the wcs-array reversed order)
            data = np.moveaxis(tensorized_data.T, -1, spectral_axis_index)

################################################################################################### Output saving

        hdu.data, hdu.header = data, header

    output_folders_sep = param.output_path.replace('/', sep).replace('\\', sep)

    file_name = input_path.split(sep)[-1]
    output_path = input_path[:-len(file_name)]
    for output_folder in output_folders_sep.split(sep):
        output_path = os.path.join(output_path, output_folder)
        if not os.path.exists(output_path):
            os.mkdir(output_path)

    hdul.writeto(output_path + file_name, overwrite = True)

################################################################################################### Data files closing

    hdul.close()

################################################################################################### Time ending

print("\x1B[32m" + f"\nExecution time: " + f"{time.time() - start:.2f} seconds\n" + "\x1B[0m")

###################################################################################################
