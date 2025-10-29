import matplotlib.pyplot as plt
import numpy as np

from astropy import units as u

import param

################################################################################################### Algorithm

def process(input_data: np.ndarray, world_grids: list,
            spectral_axis: int, axis_names: list, axis_units: list) -> np.ndarray:
    """ Stellar Subtraction methods - shared architecture
    ---
    """

    # DATA MATRICIZATION

    if param.show_verbose:
        pixels_axis_name = axis_names.pop()
        pixels_axis_unit = axis_units.pop()

    matricized = matricization(input_data, world_grids, spectral_axis, axis_names, axis_units)
    data, pre_matricization_shape, world_grids, axis_names, axis_units = matricized

    if param.show_verbose:
        show_dematricized(data, pre_matricization_shape,
                          pixels_axis_name, pixels_axis_unit,
                          world_grids, axis_names, axis_units,
                          "Pre-processed data\n", 'plasma', 'b')

    # STELLAR SPECTRUM ESTIMATION

    if param.print_verbose:
        print("\x1B[1;35m" + "\nStellar spectrum estimation" + "\x1B[0m")

    D = misleaders_masking(data.copy(), pre_matricization_shape,
                           param.stellar_mask, axis_names, world_grids)

    s_est, stellar_selection = stellar_spectrum_estimate(D)

    if param.show_verbose:

        additional_stellar_mask = ~ stellar_selection
        D[additional_stellar_mask] = np.full(D.shape[1], np.nan)

        show_dematricized(D, pre_matricization_shape,
                          pixels_axis_name, pixels_axis_unit,
                          world_grids, axis_names, axis_units,
                          "Stellar spectrum estimation\n", 'magma', 'm')

    # POINT SPREAD FUNCTION ESTIMATION

    if param.print_verbose:
        print("\x1B[1;35m" + "\nPoint Spread Function estimation" + "\x1B[0m")

    D = misleaders_masking(data.copy(), pre_matricization_shape,
                           param.psf_mask, axis_names, world_grids)

    match len(param.psf_degrees):
        case 0: psf_est = point_spread_function_estimate_S3(D, s_est, world_grids)
        case 1: psf_est = point_spread_function_estimate_S4(D, s_est, world_grids)
        case 2: psf_est = point_spread_function_estimate_S5(D, s_est, world_grids)

    if param.show_verbose:
        show_dematricized(psf_est, pre_matricization_shape,
                          pixels_axis_name + "ratio", "no unit",
                          world_grids, axis_names, axis_units,
                          "Point Spread Function estimation\n", 'cividis', 'y')

    # STELLAR COMPONENT ESTIMATION

    if param.print_verbose:
        print("\x1B[1;35m" + "\nStellar component estimation" + "\x1B[0m")

    S_est = stellar_component_estimate(psf_est, s_est)

    if param.show_verbose:
        show_dematricized(S_est, pre_matricization_shape,
                          pixels_axis_name, pixels_axis_unit,
                          world_grids, axis_names, axis_units,
                          "Stellar component estimation\n", 'inferno', 'r')

    # NOISE COMPONENT ESTIMATION

    E_est = data - S_est

    if param.show_verbose:
        show_dematricized(E_est, pre_matricization_shape,
                          pixels_axis_name, pixels_axis_unit,
                          world_grids, axis_names, axis_units,
                          "Noise component estimation\n", 'viridis', 'c')

    # DATA TENSORIZATION

    output_data = tensorization(E_est, spectral_axis, pre_matricization_shape)

    return output_data

def matricization(data_tensor: np.ndarray, world_grids: list,
                  spectral_axis: int, axis_names: list, axis_units: list) -> np.ndarray:
    """ Data matricization with spectral observations along the rows and features along the columns
    ---
    """

    # Move of the spectral axis to the last position
    # Optimization by putting the observations along the rows
    # (with a transpose '.T' against the wcs-array reversed order)
    moved_data_tensor = np.moveaxis(data_tensor.T, spectral_axis, -1)
    
    # Move of the spectral axis to the last position for grids too
    # (and associated correction of the grid axis indexes where needed)
    for axis_index in range(len(world_grids)):
        for corr_axis_index in range(len(world_grids[axis_index][0])):
            if world_grids[axis_index][0][corr_axis_index] == spectral_axis:
                world_grids[axis_index][0][corr_axis_index] = len(world_grids) - 1
            elif world_grids[axis_index][0][corr_axis_index] > spectral_axis:
                world_grids[axis_index][0][corr_axis_index] -= 1
    world_grids.append(world_grids[spectral_axis])
    world_grids.pop(spectral_axis)

    # Move of the spectral axis
    # to the last position for axis names
    axis_names.append(axis_names[spectral_axis])
    axis_names.pop(spectral_axis)

    # Move of the spectral axis
    # to the last position for axis units
    axis_units.append(axis_units[spectral_axis])
    axis_units.pop(spectral_axis)

    # Pre-matricization shape save for tensorization
    pre_matricization_shape = moved_data_tensor.shape

    # Vectorization of the first (non-spectral) axes
    # Optimization by putting the features along the columns
    matricized_data = moved_data_tensor.reshape(-1, pre_matricization_shape[-1])

    return matricized_data, pre_matricization_shape, world_grids, axis_names, axis_units

def tensorization(data_matrix: np.ndarray, spectral_axis: int,
                  pre_matricization_shape: tuple[int]) -> np.ndarray:
    """ Data tensorization
    ---
    """

    # Devectorization of the first (non-spectral) axes
    tensorized_data = data_matrix.reshape(pre_matricization_shape)

    # Move of the spectral axis to its original pre-move position
    # (with a transpose '.T' against the wcs-array reversed order)
    tensorized_data = np.moveaxis(tensorized_data.T, -1, spectral_axis)

    return tensorized_data

################################################################################################### Preparators

def misleaders_masking(matricized_data: np.ndarray,
                       pre_matricization_shape: tuple[int],
                       mask_dict: dict[str, list[tuple]],
                       axis_names: list, world_grids: list) -> np.ndarray:
    """ Misleading pixels masking
    ---
    """

    # Dematricization of the matricized data to its pre-matricization shape
    dematricized_data = matricized_data.reshape(pre_matricization_shape)

    pixels_to_be_masked = np.zeros_like(dematricized_data)

    for axis_name, masks_bounds_list in mask_dict.items():
        if axis_name not in axis_names: continue
        axis_index = axis_names.index(axis_name)

        non_axis_grid_axes = list(range(len(world_grids[axis_index][0])))
        del non_axis_grid_axes[world_grids[axis_index][0].index(axis_index)]
        axis_grid = np.mean(world_grids[axis_index][1].T, axis = tuple(non_axis_grid_axes))

        for lower_bound, upper_bound in masks_bounds_list:

            if type(lower_bound) is float:
                if axis_grid[0] < axis_grid[-1]: # In/decreasing axes cases
                    lower_bound = np.searchsorted(axis_grid, lower_bound, side = 'right') - 1
                else : lower_bound = np.searchsorted(-axis_grid, -lower_bound, side = 'left')
            if type(upper_bound) is float:
                if axis_grid[0] < axis_grid[-1]: # In/decreasing axes cases
                    upper_bound = np.searchsorted(axis_grid, upper_bound, side = 'left')
                else : upper_bound = np.searchsorted(-axis_grid, -upper_bound, side = 'right') - 1

            if lower_bound > upper_bound: lower_bound, upper_bound = upper_bound, lower_bound

            slicer = pixels_to_be_masked.ndim * [slice(None)]
            slicer[axis_index] = slice(lower_bound, upper_bound)
            pixels_to_be_masked[tuple(slicer)] += 1

    masks_intersection = (pixels_to_be_masked == max(1, np.max(pixels_to_be_masked)))

    dematricized_data[masks_intersection] = np.nan

    matricized_data = dematricized_data.reshape(matricized_data.shape)

    return matricized_data

################################################################################################### Estimators

def stellar_spectrum_estimate(D: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """ Stellar spectrum estimation
    ---
    """

    # Stellar positions selection
    data_spectral_integration = np.nansum(D, axis = 1)
    noise_level = np.nanmax(data_spectral_integration) / param.stellar_sigma
    stellar_selection = noise_level <= data_spectral_integration

    # Stellar spectrum estimation
    s_est = np.nansum(D[stellar_selection], axis = 0)

    return s_est, stellar_selection

def point_spread_function_estimate_S3(D: np.ndarray, s_est: np.ndarray,
                                      world_grids: list) -> np.ndarray:
    """ PSF estimation - version Stellar Spread Subtraction (S3)
    ---
    """

    psf_est = D / s_est

    return psf_est

def point_spread_function_estimate_S4(D: np.ndarray, s_est: np.ndarray,
                                      world_grids: list) -> np.ndarray:
    """ PSF estimation - version Stellar Spread Spectral Subtraction (S4)
    ---
    """

    # Spectrum grid construction
    # (averaging along correlated non-z axes)
    non_z_grid_axes = list(range(len(world_grids[-1][0])))
    del non_z_grid_axes[world_grids[-1][0].index(len(world_grids)-1)]
    z_grid = np.nanmean(world_grids[-1][1].T, axis = tuple(non_z_grid_axes))

    # Shift of the 'z_grid' in the Legendre [-1,1] range
    legendre_grid = 2 * (z_grid - min(z_grid)) / (max(z_grid) - min(z_grid)) - 1

    # Matrix of the Legendre polynomials
    # by following Bonnet's recursion formula
    # -> (n+1)Pn+1(x) = (2n+1)xPn(x) - nPn-1(x)
    lower_spectral_dimension = param.psf_degrees[0] + 1
    M = np.tile(legendre_grid, (lower_spectral_dimension, 1))
    M[0] = 1 # Degree 0 (total flux) Legendre polynomial
    for i in range(2, lower_spectral_dimension):
        M[i] *= M[i-1] * (2 * i - 1) / i
        M[i] -= M[i-2] * (i - 1) / i

    # Matrix of the Legendre polynomial modulation
    # of the stellar spectrum estimate 's_est'
    M_s_est = M * s_est

    # 'NaN' and '0' pixels of 's_est' masking
    # Data 'D' matrix columns full of 'NaN' masking
    s_est_mask = ~ (np.isnan(s_est) | (s_est == 0))
    all_nan_D_sheets = ~ np.all(np.isnan(D), axis = 0)
    cols_masking = all_nan_D_sheets & s_est_mask

    # Scaling of the first row (regressor) to unit norm
    # so that it weights as the others in the regression
    # (it is 's_est' as the 1st legendre polynomial is 1)
    regressor_0_norm = np.linalg.norm(s_est[cols_masking])
    M_s_est[0] /= regressor_0_norm
    M[0] /= regressor_0_norm
    # Orthogonalization of variation rows with the (first) total flux row
    # (inevitably fully correlated because of the multiplication with 's_est')
    # so that it comes '<M_s_est_ortho[i],M_s_est[0]> = M_s_est_ortho[i].T.s_est
    # is equal to '(M_s_est[i]-(<M_s_est[i],s_est>/<s_est,s_est>) s_est).T.s_est'
    # and 'M_s_est[i].T.s_est-(M_s_est[i].T.s_est/s_est.T.s_est) s_est.T.s_est'
    # and 'M_s_est[i].T.s_est-M_s_est[i].T.s_est = 0' hence the orthogonality
    # 'M' is subtracted in a way it keeps the same relation with 'M_s_est'
    colinearity = M_s_est[1:, cols_masking] @ s_est[cols_masking].T
    normalization_s_est = s_est[cols_masking] @ s_est[cols_masking].T
    normalized_colinearity = colinearity / normalization_s_est
    M_s_est[1:] -= normalized_colinearity[:, None] * s_est
    M[1:] -= normalized_colinearity[:, None]
    # Scaling of all other rows (regressors) to unit norms
    # so that each one of them weight equally in the regression
    regressors_norms = np.linalg.norm(M_s_est[1:, cols_masking], axis = 1)
    M_s_est[1:] /= regressors_norms[:, None]
    M[1:] /= regressors_norms[:, None]

    # Matrix of orthogonal projection onto the column space of 'M_s_est'
    # 'M_s_est' is well conditionned as based on the Legendre polynomials
    # Orthogonality could also ease constraining based on physical priors
    M_s_est_Gram_matrix = M_s_est[:, cols_masking] @ M_s_est[:, cols_masking].T
    Pi_M_s_est = M_s_est[:, cols_masking].T @ np.linalg.inv(M_s_est_Gram_matrix)

    # Search for data spaxels full of 'NaN'
    # (for 'NaN' reseting after 'NaN' zeroing)
    all_nan_spaxels = np.all(np.isnan(D), axis = 1)
    # 'NaN' pixels zeroing (for the matrix multiplication)
    D[np.isnan(D)] = 0
    # Estimation of the polynomial coefficients
    coeffs_est = D[:, cols_masking] @ Pi_M_s_est
    # Reset to 'NaN' of the data spaxels full of 'NaN'
    coeffs_est[all_nan_spaxels] = np.nan

    # Estimation of the PSF
    psf_est = coeffs_est @ M

    return psf_est

def point_spread_function_estimate_S5(D: np.ndarray, s_est: np.ndarray,
                                      world_grids: list) -> np.ndarray:
    """ PSF estimation - version Stellar Spread S< private > Spectral Subtraction (S5)
    ---
    """

    raise NotImplementedError("\x1B[31m" + "\nPrivate S5 method at the moment" + "\x1B[0m")

def stellar_component_estimate(psf_est: np.ndarray, s_est: np.ndarray) -> np.ndarray:
    """ Stellar component estimation
    ---
    """

    S_est = psf_est * s_est

    return S_est

################################################################################################### Shows

def show_dematricized(matricized_data: np.ndarray,
                      pre_matricization_shape: tuple[int],
                      pixels_axis_name: str, pixels_axis_unit: str,
                      world_grids: list, axis_names: list, axis_units: list,
                      title: str = "", cmap: str = 'binary', color: str = 'k') -> None:
    """ Dematricized data - image and spectrum integrations shows
    ---
    """
    
    # Dematricization of the matricized data to its pre-matricization shape
    dematricized_data = matricized_data.reshape(pre_matricization_shape)

    # Alert message for too large number of dimension
    if param.print_verbose and dematricized_data.ndim > 3:
        alert_message = "Too large number of dimension for image display:"
        alert_message += "Truncation at the first two axes (excluding the spectral axis)"
        print("\x1B[33m" + '\n' + alert_message + '\n' + "\x1B[0m")

    # Spectrum grid construction
    # (averaging along correlated non-z axes)
    non_z_grid_axes = list(range(len(world_grids[-1][0])))
    del non_z_grid_axes[world_grids[-1][0].index(len(world_grids)-1)]
    z_grid = np.mean(world_grids[-1][1].T, axis = tuple(non_z_grid_axes))
    # Image grid construction
    # (averaging along correlated non-x and non-y axes)
    non_x_grid_axes = list(range(len(world_grids[0][0])))
    del non_x_grid_axes[world_grids[0][0].index(0)]
    x_grid = np.mean(world_grids[0][1].T, axis = tuple(non_x_grid_axes))
    non_y_grid_axes = list(range(len(world_grids[1][0])))
    del non_y_grid_axes[world_grids[1][0].index(1)]
    y_grid = np.mean(world_grids[1][1].T, axis = tuple(non_y_grid_axes))

    # Filtering and blocking for spectra plots
    dematricized_data_to_plot = filtering_and_blocking(dematricized_data.copy(),
                                                       param.spectra_filter,
                                                       param.spectra_blocker,
                                                       axis_names, world_grids)
    # Filtering and blocking for images shows
    dematricized_data_to_show = filtering_and_blocking(dematricized_data.copy(),
                                                       param.images_filter,
                                                       param.images_blocker,
                                                       axis_names, world_grids)

    # Spectra construction
    non_z_axes = tuple(axis for axis in range(dematricized_data.ndim -1))
    z_nan = np.all(np.isnan(dematricized_data_to_plot), axis = non_z_axes)
    dematricized_data_to_plot[..., z_nan] = 0 # To avoid 'nanmean' warnings
    plot_spectrum = np.nanmean(dematricized_data_to_plot, axis = non_z_axes)
    plot_spectrum[z_nan] = np.nan # Erasure of fully 'NaN' sheets
    # Images construction
    non_xy_axes = tuple(axis for axis in range(2, dematricized_data.ndim))
    xy_nan = np.all(np.isnan(dematricized_data_to_show), axis = non_xy_axes)
    dematricized_data_to_show[xy_nan, ...] = 0 # To avoid 'nanmean' warnings
    pcolor_image = np.nanmean(dematricized_data_to_show, axis = non_xy_axes)
    pcolor_image[xy_nan] = np.nan # Erasure of fully 'NaN' spaxels

    # Figure configuration
    plt.figure("Dematricized show")
    plt.suptitle(title, fontsize = 24)
    # Fullscreen (depending on the graphics engine - backend - used by matplotlib)
    manager, backend = plt.get_current_fig_manager(), plt.get_backend().lower()
    if 'tkagg' in backend: manager.window.state('zoomed')
    elif 'qt' in backend: manager.window.showMaximized()
    elif "wx" in backend: manager.frame.Maximize(True)

    # World axes labels determination
    plot_xlabel_name = axis_names[-1]
    plot_xlabel_unit = axis_units[-1]
    if plot_xlabel_name in param.grids_units.keys():
        plot_xlabel_unit = param.grids_units[plot_xlabel_name]
    show_xlabel_name = axis_names[0]
    show_xlabel_unit = axis_units[0]
    if show_xlabel_name in param.grids_units.keys():
        show_xlabel_unit = param.grids_units[show_xlabel_name]
    show_ylabel_name = axis_names[1]
    show_ylabel_unit = axis_units[1]
    if show_ylabel_name in param.grids_units.keys():
        show_ylabel_unit = param.grids_units[show_ylabel_name]
    
    # Pixels axis labels determination
    try: pixels_axis_unit = u.Unit(pixels_axis_unit).to_string(format = "unicode")
    except ValueError: pass # An unrecognized unit (as "without unit") is left as is

    # Spectrum to plot
    plt.subplot(1, 2, 1)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    if z_grid[0] > z_grid[-1]: plt.gca().invert_xaxis()
    plt.plot(z_grid, plot_spectrum, linewidth = 3, color = color)
    plt.xlabel(f"{plot_xlabel_name} ({plot_xlabel_unit})", fontsize = 18)
    plt.ylabel(f"{pixels_axis_name} ({pixels_axis_unit})", fontsize = 18)
    # Image to show
    plt.subplot(1, 2, 2)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    if x_grid[0] > x_grid[-1]: plt.gca().invert_xaxis()
    if y_grid[0] > y_grid[-1]: plt.gca().invert_yaxis()
    vmin = np.nanquantile(pcolor_image, param.vminmax_quantiles[0])
    vmax = np.nanquantile(pcolor_image, param.vminmax_quantiles[1])
    plt.xlabel(f"{show_xlabel_name} ({show_xlabel_unit})", fontsize = 18)
    plt.ylabel(f"{show_ylabel_name} ({show_ylabel_unit})", fontsize = 18)
    plt.pcolor(x_grid, y_grid, pcolor_image.T, vmin = vmin, vmax = vmax, cmap = cmap)
    clrbar = plt.colorbar(fraction = 0.025, aspect = 50)
    clrbar_label = f"{pixels_axis_name} ({pixels_axis_unit})"
    clrbar.set_label(clrbar_label, fontsize = 18, rotation = 270, labelpad = 36)

    # Figure show
    plt.show()

def filtering_and_blocking(dematricized_data: np.ndarray,
                           filters: dict[str,list[tuple]],
                           blockers: dict[str,list[tuple]],
                           axis_names: list, world_grids: list) -> np.ndarray:
    """ Filtering from Gaussians with height 1 and blocking from 1 - Gaussians with height 1
    ---
    """

    # Gaussian with height 1

    gaussian = lambda grid, loc, scale: np.exp(-1/2 * ((grid-loc) / scale)**2)
 
    # Filters

    filter = np.ones_like(dematricized_data)

    for axis_name, filters_list in filters.items():
        if axis_name not in axis_names: continue
        axis_index = axis_names.index(axis_name)

        non_axis_grid_axes = list(range(len(world_grids[axis_index][0])))
        del non_axis_grid_axes[world_grids[axis_index][0].index(axis_index)]
        axis_grid = np.mean(world_grids[axis_index][1].T, axis = tuple(non_axis_grid_axes))

        axis_filter = np.zeros_like(axis_grid)

        for loc_filter, scale_filter in filters_list:

            if type(loc_filter) is int:
                loc_filter = axis_grid[loc_filter]
            if type(scale_filter) is int:
                if axis_grid[0] < axis_grid[-1]: # In/decreasing axes cases
                    loc_filter_index = np.searchsorted(axis_grid, loc_filter)
                else: loc_filter_index = np.searchsorted(axis_grid, loc_filter)
                scale_filter_right = axis_grid[loc_filter_index + scale_filter]
                scale_filter_left = axis_grid[loc_filter_index - scale_filter]
                scale_filter = np.abs(scale_filter_right - scale_filter_left) / 2

            axis_filter += gaussian(axis_grid, loc_filter, scale_filter)

        shape = [1] * dematricized_data.ndim
        shape[axis_index] = dematricized_data.shape[axis_index]

        if axis_filter.any():
            filter *= axis_filter.reshape(shape)

    dematricized_data *= filter

    # Blockers

    reversed_blocker = np.ones_like(dematricized_data)

    for axis_name, blockers_list in blockers.items():
        if axis_name not in axis_names: continue
        axis_index = axis_names.index(axis_name)

        non_axis_grid_axes = list(range(len(world_grids[axis_index][0])))
        del non_axis_grid_axes[world_grids[axis_index][0].index(axis_index)]
        axis_grid = np.mean(world_grids[axis_index][1].T, axis = tuple(non_axis_grid_axes))

        axis_blocker = np.zeros_like(axis_grid)

        for loc_blocker, scale_blocker in blockers_list:

            if type(loc_blocker) is int:
                loc_blocker = axis_filter[loc_blocker]
            if type(scale_blocker) is int:
                if axis_grid[0] < axis_grid[-1]: # In/decreasing axes cases
                    loc_blocker_index = np.searchsorted(axis_grid, loc_blocker)
                else: loc_blocker_index = np.searchsorted(axis_grid, loc_blocker)
                scale_blocker_right = axis_blocker[loc_blocker_index + scale_blocker]
                scale_blocker_left = axis_blocker[loc_blocker_index - scale_blocker]
                scale_blocker = np.abs(scale_blocker_right - scale_blocker_left) / 2

            axis_blocker += gaussian(axis_grid, loc_blocker, scale_blocker)

        shape = [1] * dematricized_data.ndim
        shape[axis_index] = dematricized_data.shape[axis_index]

        if axis_blocker.any():
            reversed_blocker *= axis_blocker.reshape(shape)

    blocker = np.ones_like(reversed_blocker) - reversed_blocker

    if blocker.any():
        dematricized_data *= blocker

    return dematricized_data
