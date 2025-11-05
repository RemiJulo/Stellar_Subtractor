import os
import time
import warnings

import numpy as np

from astropy import units
from astropy.io import fits
from astropy.wcs import WCS

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

for input_path in param.inputs_paths:

    input_path = input_path.replace('/', sep).replace('\\', sep)

    if not os.path.exists(input_path):
        raise FileNotFoundError("\x1B[31m" + input_path + "\x1B[0m")

    if os.path.isdir(input_path):
        for input_file in os.listdir(input_path):
            input_file_path = input_path + sep + input_file
            expanded_input_paths.append(input_file_path)
    else : expanded_input_paths.append(input_path)

outputs_path = "" # To create any needed series of outputs folders

for outputs_folder in param.outputs_path.replace('/', sep).replace('\\', sep).split(sep):
    outputs_path = os.path.join(outputs_path, outputs_folder)
    if outputs_path and not os.path.exists(outputs_path):
        os.mkdir(outputs_path)

################################################################################################### Data files opening

for input_path in expanded_input_paths:
    if not input_path.endswith(".fits"):
        continue # .fits files only

    input_path = input_path.replace('/', sep).replace('\\', sep)

    print("\x1B[1;34m" + f"\nFile opening: {input_path.split(sep)[-1]}\n" + "\x1B[0m")

    hdul = fits.open(input_path)

    for hdu in hdul:
        data = hdu.data
        header = hdu.header

################################################################################################### Labels extraction

        # WCS stands for 'World Coordinate System'
        if isinstance(data, np.ndarray):
            wcs = WCS(header)

            # Axis coordinate types (among 'spectral', 'celestial', 'temporal', ...)
            coordinate_types = [axis_type["coordinate_type"] for axis_type in wcs.get_axis_types()]

            # Existence and uniqueness of 1 spectral axis
            if coordinate_types.count('spectral') != 1: break
            spectral_axis = coordinate_types.index('spectral')

            # Axes names and units
            axis_names = wcs.axis_type_names
            axis_units = wcs.world_axis_units

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
            
################################################################################################### Origin relocation

            # Default origin from the header
            array_origin = np.rint(wcs.wcs.crpix-1).astype(int)
            world_origin = np.array(wcs.wcs.crval).astype(float)
            
            # 'argmax' axes save
            argmax_axes_list = []

            # Origin coordinates replacement
            # (ignoring the non-existent axes)
            for axis_name, origin in param.origin.items():
                if axis_name not in axis_names: continue
                axis_index = axis_names.index(axis_name)

                # Replacement of the associated origin coordinate with 'int' values
                # (and update of all other coordinates to deal with correlations between axes)
                if isinstance(origin, int):
                    array_origin[axis_index] = origin
                    [world_origin] = np.array(wcs.array_index_to_world_values(*[[array_origin]]))
                    if axis_name in param.origin_units.keys():
                        irreducible_unit = units.Unit(axis_units[axis_index])
                        try: # Unit conversion if possible
                            target_unit = units.Unit(param.origin_units[axis_name])
                            world_origin = (world_origin * irreducible_unit).to(target_unit).value
                        except (ValueError, units.UnitConversionError):
                            error_message = f"Axis {axis_name} unit should be in:\n"
                            error_message += f"{irreducible_unit.find_equivalent_units()}"
                            raise ValueError("\x1B[31m" + error_message + "\x1B[0m") from None
                
                # Replacement of the associated origin coordinate with 'float' values
                # (and update of all other coordinates to deal with correlations between axes)
                if isinstance(origin, float):
                    if axis_name in param.origin_units.keys():
                        irreducible_unit = units.Unit(axis_units[axis_index])
                        try: # Unit conversion if possible
                            target_unit = units.Unit(param.origin_units[axis_name])
                            origin = (origin * target_unit).to(irreducible_unit).value
                        except (ValueError, units.UnitConversionError):
                            error_message = f"Axis {axis_name} unit should be in:\n"
                            error_message += f"{irreducible_unit.find_equivalent_units()}"
                            raise ValueError("\x1B[31m" + error_message + "\x1B[0m") from None
                    world_origin[axis_index] = origin
                    [array_origin] = np.array(wcs.world_to_array_index_values(*[[world_origin]]))

                # Subsquent assignation of the 'argmax' values
                # to ensure that the final values are indeed 'argmax'
                # (despite the correlations with the potential following axes)
                if isinstance(origin, param.argmax):
                    argmax_axes_list.append(axis_index)

            # Preliminary search for the minimum and maximum flux pixels, if necessary only
            # (to avoid the redundancy of the same costly calculation with multiple 'argmax' uses)
            if argmax_axes_list:
                if param.print_verbose:
                    print("    Search for the maximum flux pixel", end = '\r')

                # Search for the brightest pixel and origin coordinates replacement
                argmax_array = np.unravel_index(np.nanargmax(data), data.shape)[::-1]
                for axis in argmax_axes_list: array_origin[axis] = argmax_array[axis]

                # Update of all other coordinates to deal with correlations between axes
                [world_origin] = np.array(wcs.array_index_to_world_values(*[[array_origin]]))

################################################################################################### Zoom selection

            # Slice, world and array zooms initialization
            array_zooms = np.array([[0]*len(data.shape[::-1]), list(data.shape[::-1])])
            world_zooms = np.array(wcs.array_index_to_world_values(*[array_zooms]))
            slice_zooms = [slice(0, axis_len) for axis_len in wcs.pixel_shape]

            # Data zooming from the origin
            # (ignoring the non-existent axes)
            for axis_name, bounds in param.zooms.items():
                if axis_name not in axis_names: continue
                axis_index = axis_names.index(axis_name)

                # Zoom lower and upper bounds
                lower_bound, upper_bound = bounds

                # Replacement of the associated zoom lower and upper bounds with 'int' values
                if isinstance(lower_bound, int):
                    if lower_bound >= 0: array_zooms[0][axis_index] = lower_bound
                    else : array_zooms[0][axis_index] = array_zooms[1][axis_index] + lower_bound
                if isinstance(upper_bound, int):
                    if upper_bound >= 0: array_zooms[1][axis_index] = upper_bound
                    else : array_zooms[1][axis_index] = array_zooms[1][axis_index] + upper_bound

                # Replacement of the associated zoom lower and upper bounds with 'float' values
                if isinstance(lower_bound, float):
                    if axis_name in param.zooms_units.keys():
                        irreducible_unit = units.Unit(axis_units[axis_index])
                        try: # Unit conversion if possible
                            target_unit = units.Unit(param.zooms_units[axis_name])
                            lower_bound = (lower_bound * target_unit).to(irreducible_unit).value
                        except (ValueError, units.UnitConversionError):
                            error_message = f"Axis {axis_name} unit should be in:\n"
                            error_message += f"{irreducible_unit.find_equivalent_units()}"
                            raise ValueError("\x1B[31m" + error_message + "\x1B[0m") from None
                    world_zooms[0][axis_index] = world_origin[axis_index] + lower_bound
                if isinstance(upper_bound, float):
                    if axis_name in param.zooms_units.keys():
                        irreducible_unit = units.Unit(axis_units[axis_index])
                        try: # Unit conversion if possible
                            target_unit = units.Unit(param.zooms_units[axis_name])
                            upper_bound = (upper_bound * target_unit).to(irreducible_unit).value
                        except (ValueError, units.UnitConversionError):
                            error_message = f"Axis {axis_name} unit should be in:\n"
                            error_message += f"{irreducible_unit.find_equivalent_units()}"
                            raise ValueError("\x1B[31m" + error_message + "\x1B[0m") from None
                    world_zooms[1][axis_index] = world_origin[axis_index] + upper_bound

            # Add of the world zooms to the array zooms
            # The data shape is also subtracted to the array zooms
            # (to deal with the duplicate caused by double array-world initialization)
            world_zooms_to_array_zooms = wcs.world_to_array_index_values(*[world_zooms])
            array_zooms += np.array(world_zooms_to_array_zooms)[::-1]
            array_zooms[1] -= np.array(data.shape[::-1])
            array_bounds = list(zip(*array_zooms))

            # Zoom validity checking
            # (to avoid errors caused by flat axes)
            for axis_index, bounds in enumerate(array_bounds):
                axis_name, (lower_bound, upper_bound) = axis_names[axis_index], bounds
                if lower_bound > upper_bound: # Decreasing order axis cases
                    lower_bound, upper_bound = upper_bound, lower_bound
                    array_bounds[axis_index] = (lower_bound, upper_bound)
                if lower_bound > wcs.pixel_shape[axis_index]: # Overflow beyond the pixels number
                    error_message = f"Lower bound too large for the '{axis_name}' axis"
                    raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
                if upper_bound < 0 : # Overflow beyond the index zero origin of the array axes
                    error_message = f"Upper bound too small for the '{axis_name}' axis"
                    raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

            # Zoom arrays to slices conversion
            for axis_name, bounds in param.zooms.items():
                if axis_name not in axis_names: continue
                axis_index = axis_names.index(axis_name)
                list_array_zooms_key = list(array_bounds[axis_index])
                list_array_zooms_key[0] = max(0, list_array_zooms_key[0]) # For positive indexes
                list_array_zooms_key[1] += 1 # To include the upper bounds in the zoom selection
                slice_zooms[axis_index] = slice(*list_array_zooms_key)

            # WCS, data and header zooming
            tuple_slice_zooms = tuple(slice_zooms)[::-1]
            zoomed_wcs = wcs.slice(tuple_slice_zooms)
            zoomed_data = data[tuple_slice_zooms]
            header = zoomed_wcs.to_header()
            
################################################################################################### Data sum up

            if param.print_verbose:

                world_min_max = wcs.array_index_to_world_values(*[array_zooms])
                [world_origin] = wcs.array_index_to_world_values(*[[array_origin]])

                axes_information = [axis_names, axis_units]
                axes_information += [array_origin, world_origin]
                axes_information += array_zooms.tolist() + list(world_min_max)

                i = zip(*tuple(axes_information)) # Information list for all axes

                for name, unit, array_0, world_0, array_min, array_max, world_min, world_max in i:

                    if name in param.origin_units.keys():
                        irreducible_unit = units.Unit(axis_units[axis_names.index(name)])
                        try: # Unit conversion if possible
                            unit = units.Unit(param.origin_units[name])
                            world_0 = (world_0 * irreducible_unit).to(unit).value
                            world_min = (world_min * irreducible_unit).to(unit).value
                            world_max = (world_max * irreducible_unit).to(unit).value
                        except (ValueError, units.UnitConversionError):
                            error_message = f"Axis {name} unit should be in:\n"
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

                # Optimization by avoiding the calculation of the grids
                # already determined for the axes with the same correlations
                if world_grids[axis][0] == corr_axes.tolist(): continue

                # World grids determination from 'zoomed_wcs' and the correlated axes array grids
                needed_wcs = zoomed_wcs.sub(tuple(int(corr_axis) + 1 for corr_axis in corr_axes))
                needed_grids = needed_wcs.array_index_to_world_values(*np.indices(corr_shapes))

                # Save of the grid for the given axis
                # as well as of the grids for the correlated axes
                # (generally without additional correlations with other axes)
                for corr_axis, needed_grid in zip(corr_axes, np.atleast_2d(needed_grids)):
                    needed_grid -= world_origin[corr_axis] # Shift to the origin coordinate
                    if axis_names[corr_axis] in param.zooms_units.keys():
                        irreducible_unit = units.Unit(axis_units[corr_axis])
                        try: # Unit conversion if possible
                            target_unit = units.Unit(param.zooms_units[axis_names[corr_axis]])
                            needed_grid = (needed_grid * irreducible_unit).to(target_unit).value
                        except (ValueError, units.UnitConversionError):
                            error_message = f"Axis {axis_names[corr_axis]} unit should be in:\n"
                            error_message += f"{irreducible_unit.find_equivalent_units()}"
                            raise ValueError("\x1B[31m" + error_message + "\x1B[0m") from None
                    world_grids[corr_axis] = (corr_axes.tolist(), needed_grid)

################################################################################################### Main processing

            data = process.process(zoomed_data, spectral_axis, world_grids, axis_names, axis_units)

################################################################################################### Outputs saving

            hdu.header.update(header)
            hdu.data = data
    
    if param.save_verbose:
        hdul.writeto(outputs_path + input_path.split(sep)[-1], overwrite = True)

################################################################################################### Data files closing

    hdul.close()

################################################################################################### Time ending

print("\x1B[32m" + f"\nExecution time: " + f"{time.time() - start:.2f} seconds\n" + "\x1B[0m")

###################################################################################################
