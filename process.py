import matplotlib.pyplot as plt
import numpy as np

import param

################################################################################################### Algorithm

def process(data, world_grids, axis_labels):
    """ Stellar Subtraction methods - shared architecture
    ---
    """

    # Data recall

    if param.show_verbose:
        
        (pixels_axis_name, pixels_axis_unit) = axis_labels[0]
        pixels_axis_unit = pixels_axis_unit.replace('**', '^')
        pixels_axis_unit = pixels_axis_unit.replace('(','{')
        pixels_axis_unit = pixels_axis_unit.replace(')','}')
        pixels_axis_unit = f"${pixels_axis_unit}$"

        if pixels_axis_name == "< name >":
            pixels_axis_name = "Flux"

        show_dematricized(data, world_grids, axis_labels,
                          "Flux", pixels_axis_unit,
                          "Input data\n", 'plasma', 'b')

    # Stellar spectrum estimation

    if param.print_verbose:
        print("\x1B[1;35m" + "\nStellar spectrum estimation" + "\x1B[0m")

    s_est, stellar_selection = stellar_spectrum_estimate(data)

    if param.show_verbose:

        stellar_heart = data.copy()
        nan_spaxels = np.full(data.shape[1], np.nan)
        stellar_heart[~stellar_selection] = nan_spaxels

        show_dematricized(stellar_heart, world_grids, axis_labels,
                          "Flux", pixels_axis_unit,
                          "Stellar spectrum estimates\n", 'magma', 'm')

    # Point Spread Function estimation

    if param.print_verbose:
        print("\x1B[1;35m" + "\nPoint Spread Function estimation" + "\x1B[0m")

    match len(param.hyperparameter):
        case 0: psf_est = point_spread_function_estimate_S3(data, s_est, world_grids)
        case 1: psf_est = point_spread_function_estimate_S4(data, s_est, world_grids)
        case 2: psf_est = point_spread_function_estimate_S5(data, s_est, world_grids)

    if param.show_verbose:
        show_dematricized(psf_est, world_grids, axis_labels,
                          "Flux" + " ratio", "$cm^{-2}.angstrom^{-1}$",
                          "PSF estimates\n", 'cividis', 'y')

    # Stellar component estimation

    if param.print_verbose:
        print("\x1B[1;35m" + "\nStellar component estimation" + "\x1B[0m")

    S_est = stellar_component_estimate(psf_est, s_est)

    if param.show_verbose:
        show_dematricized(S_est, world_grids, axis_labels,
                          "Flux", pixels_axis_unit,
                          "Stellar estimates\n", 'inferno', 'r')

    # Residuals estimation

    residual_component = data - S_est

    if param.show_verbose:
        show_dematricized(residual_component, world_grids, axis_labels,
                          "Flux", pixels_axis_unit,
                          "Output residuals\n", 'viridis', 'c')

    return residual_component

################################################################################################### Estimators

def stellar_spectrum_estimate(data):
    """ Stellar spectrum estimation
    ---
    """

    # Stellar positions selection
    data_spectral_integration = np.nansum(data, axis = 1)
    noise_level = np.nanmax(data_spectral_integration) / param.stellar_sigma
    stellar_selection = noise_level <= data_spectral_integration

    # Stellar spectrum estimation
    s_est = np.nansum(data[stellar_selection], axis = 0)

    return s_est, stellar_selection

def point_spread_function_estimate_S3(data, s_est, world_grids):
    """ PSF estimation - version Stellar Spread Subtraction (S3)
    ---
    """

    psf_est = data / s_est

    return psf_est

def point_spread_function_estimate_S4(data, s_est, world_grids):
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
    lower_dimension = param.hyperparameter[0] + 1
    M = np.tile(legendre_grid, (lower_dimension, 1))
    M[0] = 1 # Degree 0 Legendre polynomial
    for i in range(2, lower_dimension):
        M[i] *= M[i-1] * (2 * i - 1) / i
        M[i] -= M[i-2] * (i - 1) / i

    # Matrix of the Legendre polynomial modulation
    # of the stellar spectrum estimate 's_est'
    M_s_est = M * s_est

    # 'NaN' and '0' pixels of 's_est' masking
    s_est_mask = ~ (np.isnan(s_est) | (s_est == 0))

    # Scaling of the first row (regressor) to unit norm
    # so that it weights as the others in the regression
    # (it is 's_est' as the 1st legendre polynomial is 1)
    regressor_0_norm = np.linalg.norm(s_est[s_est_mask])
    M_s_est[0] /= regressor_0_norm
    M[0] /= regressor_0_norm
    # Orthogonalization of variation rows with the (first) total flux row
    # (inevitably fully correlated because of the multiplication with 's_est')
    # so that it comes '<M_s_est_ortho[i],M_s_est[0]> = M_s_est_ortho[i].T.s_est
    # is equal to '(M_s_est[i]-(<M_s_est[i],s_est>/<s_est,s_est>) s_est).T.s_est'
    # and 'M_s_est[i].T.s_est-(M_s_est[i].T.s_est/s_est.T.s_est) s_est.T.s_est'
    # and 'M_s_est[i].T.s_est-M_s_est[i].T.s_est = 0' hence the orthogonality
    # 'M' is subtracted in a way it keeps the same relation with 'M_s_est'
    colinearity = M_s_est[1:, s_est_mask] @ s_est[s_est_mask].T
    normalization_s_est = s_est[s_est_mask] @ s_est[s_est_mask].T
    normalized_colinearity = colinearity / normalization_s_est
    M_s_est[1:] -= normalized_colinearity[:, None] * s_est
    M[1:] -= normalized_colinearity[:, None]
    # Scaling of all other rows (regressors) to unit norms
    # so that each one of them weight equally in the regression
    regressors_norms = np.linalg.norm(M_s_est[1:, s_est_mask], axis = 1)
    M_s_est[1:] /= regressors_norms[:, None]
    M[1:] /= regressors_norms[:, None]

    # Matrix of orthogonal projection onto the column space of 'M_s_est'
    # 'M_s_est' is well conditionned as based on the Legendre polynomials
    # Orthogonality could also ease constraining based on physical priors
    M_s_est_Gram_matrix = M_s_est[:, s_est_mask] @ M_s_est[:, s_est_mask].T
    Pi_M_s_est = M_s_est[:, s_est_mask].T @ np.linalg.inv(M_s_est_Gram_matrix)

    # Search for data spaxels full of 'NaN'
    # (for 'NaN' reseting after 'NaN' zeroing)
    all_nan_data = np.all(np.isnan(data), axis = 1)
    # 'NaN' data pixels zeroing
    # (for the matrix product)
    data[np.isnan(data)] = 0
    # Estimation of the polynomial coefficients
    coeffs_est = data[:, s_est_mask] @ Pi_M_s_est
    # Reset to 'NaN' of the data spaxels full of 'NaN'
    coeffs_est[all_nan_data] = np.nan

    # Estimation of the PSF
    psf_est = coeffs_est @ M

    return psf_est

def point_spread_function_estimate_S5(data, s_est, world_grids):
    """ PSF estimation - version Stellar Spread S< private > Spectral Subtraction (S5)
    ---
    """

    raise NotImplementedError("\x1B[31m" + "\nPrivate S5 method at the moment" + "\x1B[0m")

def stellar_component_estimate(psf_est, s_est):
    """ Stellar component estimation
    ---
    """

    S_est = psf_est * s_est

    return S_est

################################################################################################### Shows

def show_dematricized(matricized_data, world_grids, axis_labels,
                      pixels_axis_name = "< pixels axis name >",
                      pixels_axis_unit = "< pixels axis unit >",
                      title = "", cmap = 'binary', color = 'k'):
    """ Dematricized data - image and spectrum integrations shows
    ---
    """

    # Pre-matricization shape recovery from world_grids
    pre_matricization_shape = np.full(len(world_grids), -1)
    for axes_grid, grid in world_grids:
        for indice_axe, axe in enumerate(axes_grid):
            pre_matricization_shape[axe] = grid.T.shape[indice_axe]
    # Pre-matricization reconstruction from pre-matricization shape
    dematricized_data = matricized_data.reshape(pre_matricization_shape)

    # Alert message for too large number of dimension
    if param.print_verbose and len(pre_matricization_shape) > 3:
        alert_message = "Too large number of dimension for image display:"
        alert_message += "Truncation at the first two axes (excluding the spectral axis)"
        print("\x1B[33m" + '\n' + alert_message + '\n' + "\x1B[0m")

    # Spectrum grid construction
    # (averaging along correlated non-z axes)
    non_z_grid_axes = list(range(len(world_grids[-1][0])))
    del non_z_grid_axes[world_grids[-1][0].index(len(world_grids)-1)]
    z_grid = np.nanmean(world_grids[-1][1].T, axis = tuple(non_z_grid_axes))

    # Image grid construction
    # (averaging along correlated non-x and non-y axes)
    non_x_grid_axes = list(range(len(world_grids[0][0])))
    del non_x_grid_axes[world_grids[0][0].index(0)]
    non_y_grid_axes = list(range(len(world_grids[1][0])))
    del non_y_grid_axes[world_grids[1][0].index(1)]
    x_grid = np.nanmean(world_grids[0][1].T, axis = tuple(non_x_grid_axes))
    y_grid = np.nanmean(world_grids[1][1].T, axis = tuple(non_y_grid_axes))

    # Spectrum construction
    non_z_axes = tuple(axis for axis in range(len(dematricized_data.shape)-1))
    z_nan = np.all(np.isnan(dematricized_data), axis = tuple(non_z_axes))
    spectrum_to_plot = np.nansum(dematricized_data, axis = non_z_axes)
    spectrum_to_plot[z_nan] = np.nan # Erasure of fully 'NaN' sheets
    
    # Image construction
    non_xy_axes = [axis for axis in range(2, len(dematricized_data.shape))]
    xy_nan = np.all(np.isnan(dematricized_data), axis = tuple(non_xy_axes))
    image_to_show = np.nansum(dematricized_data, axis = tuple(non_xy_axes))
    image_to_show[xy_nan] = np.nan # Erasure of fully 'NaN' spaxels

    # Figure configuration
    plt.figure("Dematricized show")
    plt.suptitle(title, fontsize = 24)
    # Fullscreen (depending on the graphics engine - backend - used by matplotlib)
    manager, backend = plt.get_current_fig_manager(), plt.get_backend().lower()
    if 'tkagg' in backend: manager.window.state('zoomed')
    elif 'qt' in backend: manager.window.showMaximized()
    elif "wx" in backend: manager.frame.Maximize(True)

    # Labels determination
    (plot_xlabel_name, plot_xlabel_unit) = axis_labels[-1]
    if plot_xlabel_name in param.grids_units.keys():
        plot_xlabel_unit = param.grids_units[plot_xlabel_name]
    (show_xlabel_name, show_xlabel_unit) = axis_labels[1]
    if show_xlabel_name in param.grids_units.keys():
        show_xlabel_unit = param.grids_units[show_xlabel_name]
    (show_ylabel_name, show_ylabel_unit) = axis_labels[2]
    if show_ylabel_name in param.grids_units.keys():
        show_ylabel_unit = param.grids_units[show_ylabel_name]
    if ".angstrom^{-1}" in pixels_axis_unit:
        pixels_axis_unit_spectrum = pixels_axis_unit.replace(".angstrom^{-1}", "")
    else : pixels_axis_unit_spectrum = "< unit >"
    if "cm^{-2}." in pixels_axis_unit:
        pixels_axis_unit_image = pixels_axis_unit.replace("cm^{-2}.", "")
    else : pixels_axis_unit_image = "< unit >"

    # Spectrum to plot
    plt.subplot(1, 2, 1)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    if z_grid[0] > z_grid[-1] : plt.gca().invert_xaxis()
    plt.plot(z_grid, spectrum_to_plot, linewidth = 3, color = color)
    plt.xlabel(f"{plot_xlabel_name} ({plot_xlabel_unit})", fontsize = 18)
    plt.ylabel(f"{pixels_axis_name} ({pixels_axis_unit_spectrum})", fontsize = 18)
    # Image to show
    plt.subplot(1, 2, 2)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    if x_grid[0] > x_grid[-1] : plt.gca().invert_xaxis()
    if y_grid[0] > y_grid[-1] : plt.gca().invert_yaxis()
    plt.pcolor(x_grid, y_grid, image_to_show.T, cmap = cmap)
    plt.xlabel(f"{show_xlabel_name} ({show_xlabel_unit})", fontsize = 18)
    plt.ylabel(f"{show_ylabel_name} ({show_ylabel_unit})", fontsize = 18)
    clrbar = plt.colorbar(fraction = 0.025, aspect = 50)
    clrbar_label = f"{pixels_axis_name} ({pixels_axis_unit_image})"
    clrbar.set_label(clrbar_label, fontsize = 18, rotation = 270, labelpad = 36)

    # Figure show
    plt.show()
