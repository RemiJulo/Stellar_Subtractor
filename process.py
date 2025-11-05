import matplotlib.pyplot as plt
import numpy as np

import os

from astropy.units import Unit
from astropy.stats import sigma_clip

import param

###################################################################################################

type float_arr_type = np.typing.NDArray[np.float64]
type mask_type = dict[str, list[tuple[int | float, int | float]]]
type tuple_float_arr_type = tuple[float_arr_type, float_arr_type]

################################################################################################### General algorithm

def process(data_tensor: float_arr_type, spectral_axis: int,
            world_grids: list[tuple[list[int], float_arr_type]],
            axis_names: list[str], axis_units: list[str]) -> float_arr_type:
    """ Stellar Subtraction (S2) methods - General algorithm
    ---

    Shared structure for the different stellar subtraction methods:

    - Outliers detection (for their disregard in the following processing)
    - Data matricization (with observations along the rows and features along the columns)
    - Start of the estimation loop (with increasingly better noise estimates)
        - Stellar spectrum estimation (by averaging a selection of the least noisy pixels)
        - Point Spread Function estimation (with the S3, S4 or S5 method)
        - Stellar component estimation (from the stellar spectrum and PSF estimates)
        - Noise component estimation (with covariances matrices estimations)
    - End of the estimation loop (in practice, after two or even just one iteration)
    - Data tensorization (to the original shape of the input tensor data)
    - Anomalies detection (for their highlight in the post-subtraction report)
    
    Args:
        data_tensor (np.typing.NDArray[np.float64]):
            Data tensor from which to subtract the stellar component
        spectral_axis (int):
            Index of the spectral axis on which both the S4 and the S5 methods are based
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the grid of the corresponding axes, in the same access order
        axis_names (list[str]):
            Names of the different axes as extracted in the data header
        axis_units (list[str]):
            Units of the different axes as extracted in the data header \\
            or as choosen in `param.zooms_units`

    Returns:
        data_tensor (np.typing.NDArray[np.float64]):
            Data tensor from which an estimate of the stellar component has been subtracted
    """

    # Outliers detection

    data_tensor[outliers_pixels_estimate(data_tensor, axis_names, world_grids)] = np.nan

    # Data matricization

    tensor_args = locals() # Tensor data args for transformation into matrix data args

    if param.show_verbose:
        pixels_axis_name = axis_names.pop() # The pixels axis name and unit is added just
        pixels_axis_unit = axis_units.pop() # at the end of the axis names and units lists

    data, data_tensor_shape = matricization(**tensor_args)

    if param.show_verbose:
        show_dematricized(data.reshape(data_tensor_shape),
                          pixels_axis_name, pixels_axis_unit,
                          axis_names, axis_units, world_grids,
                          "Pre-processed data\n", 'b', 'plasma')

    # Estimation loop

    covariances = [np.nan, np.nan] # Rows and columns covariances

    for noise_loop in range(1 + param.noise_loopback):

        if param.print_verbose and param.noise_loopback:
            print("\x1B[1;35m" + f"\nEstimation loop: {noise_loop + 1}" + "\x1B[0m")

        # Stellar spectrum estimation

        if param.print_verbose:
            print("\x1B[1;35m" + "\n    Stellar spectrum estimation" + "\x1B[0m")

        D = masking_tool(data.reshape(data_tensor_shape).copy(), axis_names,
                         world_grids, param.stellar_mask).reshape(data.shape)

        s_est = stellar_spectrum_estimate(D)

        if param.show_verbose:
            show_dematricized(D.reshape(data_tensor_shape),
                              pixels_axis_name, pixels_axis_unit,
                              axis_names, axis_units, world_grids,
                              "Stellar spectrum estimation\n", 'm', 'magma')

        # Point Spread Function estimation

        if param.print_verbose:
            print("\x1B[1;35m" + "\n    Point Spread Function estimation" + "\x1B[0m")

        D = masking_tool(data.reshape(data_tensor_shape).copy(), axis_names,
                         world_grids, param.psf_mask).reshape(data.shape)

        match len(param.psf_degrees):
            case 0: psf_est = point_spread_function_estimate_S3(D, s_est, world_grids, covariances)
            case 1: psf_est = point_spread_function_estimate_S4(D, s_est, world_grids, covariances)
            case 2: psf_est = point_spread_function_estimate_S5(D, s_est, world_grids, covariances)

        if param.show_verbose:
            show_dematricized(psf_est.reshape(data_tensor_shape),
                              pixels_axis_name + "ratio", "no unit",
                              axis_names, axis_units, world_grids,
                              "Point Spread Function estimation\n", 'y', 'cividis')

        # Stellar component estimation

        if param.print_verbose:
            print("\x1B[1;35m" + "\n    Stellar component estimation" + "\x1B[0m")

        S_est = stellar_component_estimate(psf_est, s_est)

        if param.show_verbose:
            show_dematricized(S_est.reshape(data_tensor_shape),
                              pixels_axis_name, pixels_axis_unit,
                              axis_names, axis_units, world_grids,
                              "Stellar component estimation\n", 'r', 'inferno')

        # Noise component estimation

        E_est = data - S_est

        if param.show_verbose:
            show_dematricized(E_est.reshape(data_tensor_shape),
                              pixels_axis_name, pixels_axis_unit,
                              axis_names, axis_units, world_grids,
                              "Noise component estimation\n", 'c', 'viridis')

        if param.noise_loopback:
            covariances = noise_covariance_estimate(E_est)

    # Anomalies detection

    if param.show_verbose:

        anomalies_selection = E_est.reshape(data_tensor_shape).copy()

        anomalies = outliers_pixels_estimate(anomalies_selection, axis_names, world_grids,
                                             param.anomalies_sigmas, param.anomalies_mask)
        
        anomalies_selection[ ~ anomalies] = np.nan

        show_dematricized(anomalies_selection,
                          pixels_axis_name, pixels_axis_unit,
                          axis_names, axis_units, world_grids,
                          "Anomaly detection\n", 'g', 'viridis')

    # Data tensorization

    data_tensor = tensorization(E_est, spectral_axis, data_tensor_shape)

    return data_tensor

def matricization(data_tensor: float_arr_type, spectral_axis: int,
                  world_grids: list[tuple[list[int], float_arr_type]],
                  axis_names: list[str], axis_units: list[str]) -> tuple[float_arr_type,
                                                                         tuple[int]]:
    """ Data matricization
    ---

    In this `process.py`, the features are the spectral channels \\
    while the observations are the various spaxels collected following any other type of axis

    **Note:** The observations are set along the rows and the features are set along the columns \\
    for optimization purposes (the observations/spaxels being used as a single block usually, \\
    and therefore preferably contiguously - the exception being when using spectral masks)

    Args:
        data_tensor (np.typing.NDArray[np.float64]):
            Data tensor from which to subtract the stellar component \\
            and transform into a data matrix in this function
        spectral_axis (int):
            Index of the spectral axis on which both the S4 and the S5 methods are based
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the grid of the corresponding axes, in the same access order
        axis_names (list[str]):
            Names of the different axes as extracted in the data header
        axis_units (list[str]):
            Units of the different axes as extracted in the data header \\
            or as choosen in `param.zooms_units`

    Returns:
        matricix_args (tuple[np.typing.NDArray[np.float64]]], tuple[int]]):
        Matricized data with observations along the rows and features along the columns, \\
        as well as pre-matricization (but post spectral axis move) shape of the tensor data \\
        (importantly, this also modifies `world_grids`, `axis_names` and `axis_units`)
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

    # Data tensor shape save for tensorization
    data_tensor_shape = moved_data_tensor.shape

    # Vectorization of the first (non-spectral) axes
    # Optimization by putting the features along the columns
    matricized_data = moved_data_tensor.reshape(-1, data_tensor_shape[-1])

    # Matricized data and pre-matricization shape
    matricix_args = matricized_data, data_tensor_shape

    return matricix_args

def tensorization(data_matrix: float_arr_type, spectral_axis: int,
                  pre_matricization_shape: tuple[int]) -> float_arr_type:
    """ Data tensorization
    ---

    Reshape of the matricized data to its original tensor shape \\
    and with reposition of the spectral axis to its original position
    
    Args:
        data_matrix (np.typing.NDArray[np.float64]):
            Matricized data with observations along the rows and features along the columns \\
            (following the processing previously done in the `matricization` function)
        spectral_axis (int):
            Original index of the spectral axis (where to move the last - features/spectral - axis)
        pre_matricization_shape (tuple[int]):
            Original shape of the input tensor data (how to reshape the matriced data)

    Returns:
        tensorized_data (np.typing.NDArray[np.float64]):
            Tensorized data matrix to the original shape of the input tensor data \\
            (with reposition of the spectral axis to its original position)
    """

    # Devectorization of the first (non-spectral) axes
    pre_move_tensorized_data = data_matrix.reshape(pre_matricization_shape)

    # Move of the spectral axis to its original position
    # (with a transpose '.T' against the wcs-array reversed order)
    tensorized_data = np.moveaxis(pre_move_tensorized_data.T, -1, spectral_axis)

    return tensorized_data

################################################################################################### Estimators

def outliers_pixels_estimate(data_tensor: float_arr_type, axis_names: list[str],
                             world_grids: list[tuple[list[int], float_arr_type]],
                             standard_deviations: dict[str, int | float] = param.outliers_sigmas,
                             mask: mask_type = param.outliers_mask) -> np.typing.NDArray[np.bool]:
    """ Outliers pixels detection
    ---

    Extraction of the outliers pixels, whether for their removal or their highlight:
    
    - A pixel is an outlier when its pixel value deviates from the mean of the pixel values \\
    by as many standard deviations of the pixel values as set in `standard_deviations`, \\
    each time following the corresponding axis (for each entry of `standard_deviations`)

    **Note:** The order of the `standard_deviations` entries can theoretically \\
    impact the outliers actually detected, since these are identified axis by axis \\
    (though, with a sufficient number of pixels, this is not supposed to happen)

    Args:
        data_tensor (np.typing.NDArray[np.float64]):
            Data tensor from which to identify the outlier pixels
        axis_names (list[str]):
            Names of the different axes as extracted in the data header
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the grid of the corresponding axes, in the same access order
        standard_deviations (dict[str, int | float], optional):
            Numbers of standard deviations above the means above which the values of the pixels \\
            classifie them as outliers (and lead to returning a `True` entry at their position), \\
            for each given axis (as long as its name is found in the header via `axis_names`) \\
            Defaults to `param.outliers_sigmas` (since used only for outliers by default)
        mask (dict[str, list[tuple[int | float, int | float]]], optional):
            Lists of intervals between which pixels are ignored when searching for outliers, \\
            for each given axis (as long as its name is found in the header via `axis_names`) \\
            Defaults to `param.outliers_mask` (since used only for outliers by default)

    Returns:
        outliers_pixels (np.typing.NDArray[np.bool]):
            Mask identifying the outlier pixels
    """

    outliers_pixels = np.full_like(data_tensor, False, dtype = bool)

    data_tensor = masking_tool(data_tensor, axis_names, world_grids, mask)

    for axis_name, sigma in standard_deviations.items():
        if axis_name not in axis_names: continue
        axis_index = axis_names.index(axis_name)

        data_tensor = sigma_clip(data_tensor, sigma, cenfunc = 'mean', axis = axis_index)

        outliers_pixels = data_tensor.mask

    if param.print_verbose:
        outliers_number = np.count_nonzero(outliers_pixels)
        if outliers_number:
            print("\x1B[33m" + f"\n-- {outliers_number} outliers identified --" + "\x1B[0m")

    return outliers_pixels

def stellar_spectrum_estimate(D: float_arr_type) -> float_arr_type:
    """ Stellar spectrum estimation
    ---
    
    Estimation of the spectrum of the main (stellar) component of the data, \\
    by mean of the pixels whose spectrally integrated flux is greater than the threshold \\
    given by `param.stellar_level` (by comparing each of these fluxes with the maximum of them)

    Args:
        D (np.typing.NDArray[np.float64]):
            Data matrix from which to estimate the stellar spectrum

    Returns:
        s_est (np.typing.NDArray[np.float64]):
            Stellar spectrum estimate by mean of the selected pixels
    """

    # Stellar positions selection
    data_spectral_integration = np.nansum(D, axis = 1)
    noise_level = np.nanmax(data_spectral_integration) / param.stellar_level
    stellar_selection = noise_level <= data_spectral_integration

    # Stellar spectrum estimation
    s_est = np.nansum(D[stellar_selection], axis = 0)

    # Masking for shows
    if param.show_verbose:
        D[~stellar_selection] = np.full(D.shape[1], np.nan)

    return s_est

def point_spread_function_estimate_S3(D: float_arr_type, s_est: float_arr_type,
                                      world_grids: list[tuple[list[int], float_arr_type]],
                                      covariances: tuple_float_arr_type) -> float_arr_type:
    """ PSF estimation - Stellar Spread Subtraction (S3) version
    ---

    Estimation of the Point Spread Function matrix  using the S3 method, \\
    that is to say directly by division of the data with the stellar spectrum estimate

    - As there is no regularization with this method, the subtraction of the stellar component \\
    in the following steps should lead to the null matrix following the mathematical model: \\
    in practice, this reveals the numerical instabilities inherent to any machine calculation
    - Yet, this approach is an usefull tool to get an idea of the smoothness of the PSF

    Args:
        D (np.typing.NDArray[np.float64]):
            Data matrix from which to estimate the Point Spread Function
        s_est (np.typing.NDArray[np.float64]):
            Stellar spectrum estimate by mean of the selected pixels \\
            (with the `stellar_spectrum_estimate` function)
        covariances (tuple[np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]):
            Rows and columns noise covariances matrices estimates \\
            (with the `noise_covariance_estimate` function) \\
            or `np.nan` to indicate unestimated covariances
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the grid of the corresponding axes, in the same access order

    Returns:
        psf_est (np.typing.NDArray[np.float64]):
            Point Spread Function estimate with the S3 method
    """

    out = np.full_like(D, np.nan) # Default 'psf_est' value
    where = (s_est != 0) & ~ np.isnan(s_est) # Divisibility
    psf_est = np.divide(D, s_est, where = where, out = out)

    # Save of the PSF
    if param.save_verbose:
        # Dictionnary associating the spectral channels with their flux values
        psf_est_dict = {f"Channel {j}" : psf_est[:, j] for j in range(len(psf_est[0]))}
        # Name of the .npz file in which are saved the spectral channels flux values
        save_name = f"_PSF_{int(len(os.listdir(param.outputs_path))/2)}_.npz"
        # Save of the Legendre coefficients representing the reduced PSF
        np.savez(param.outputs_path + save_name, **psf_est_dict)

    return psf_est

def point_spread_function_estimate_S4(D: float_arr_type, s_est: float_arr_type,
                                      world_grids: list[tuple[list[int], float_arr_type]],
                                      covariances: tuple_float_arr_type) -> float_arr_type:
    """ PSF estimation - Stellar Spread Subtraction (S4) version
    ---

    Estimation of the Point Spread Function matrix using the S4 method, \\
    that is to say by factorizing the PSF matrix with one Legendre polynomial matrix, \\
    on the features / spectral side (using the assumption of PSF smoothness in this direction)

    - The regularization that this factorization causes allows a much robust (yet biaised) \\
    estimation of the PSF with regard to the noise (or even to any other source in the data)
    - This comes to estimate the stellar component of the data with a polynomial modulation \\
    of the stellar spectrum estimate (in a robust way thanks to the Legendre Basis)

    With this factorization and the assumption of a Gaussian noise, \\
    the maximum likehood estimate of the Point Spread Function matrix becomes:

        D C_w Ms (Ms C_w Ms^T)^-1 Ms

    with `Ms` the matrix of Legendre polynomial modulation of the stellar spectrum estimate \\
    and `C_w` the covariance matrix of the noise in the rows direction

    Args:
        D (np.typing.NDArray[np.float64]):
            Data matrix from which to estimate the Point Spread Function
        s_est (np.typing.NDArray[np.float64]):
            Stellar spectrum estimate by mean of the selected pixels \\
            (with the `stellar_spectrum_estimate` function)
        covariances (tuple[np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]):
            Rows and columns noise covariances matrices estimates \\
            (with the `noise_covariance_estimate` function) \\
            or `np.nan` to indicate unestimated covariances
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the grid of the corresponding axes, in the same access order

    Returns:
        psf_est (np.typing.NDArray[np.float64]):
            Point Spread Function estimate with the S4 method
    """

    # Spectrum grid construction
    # (averaging along correlated non-rows axes)
    non_rows_grid_axes = list(range(len(world_grids[-1][0])))
    del non_rows_grid_axes[world_grids[-1][0].index(len(world_grids)-1)]
    rows_grid = np.nanmean(world_grids[-1][1].T, axis = tuple(non_rows_grid_axes))

    # Shift of the 'rows_grid' in the Legendre [-1,1] range
    legendre_grid = 2 * (rows_grid - min(rows_grid)) / (max(rows_grid) - min(rows_grid)) - 1

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
    s_est_mask = ~ (np.isnan(s_est) | (s_est == 0))
    # Data 'D' matrix columns full of 'NaN' masking
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

    # Data 'D' matrix 'NaN' values zeroing
    all_nan_D_spaxels = np.all(np.isnan(D), axis = 1)
    D[np.isnan(D)] = 0 # For the matrix multiplications
    
    # Precision matrix weighting
    if covariances[0] is np.nan: covariances[0] = np.identity(D.shape[1])
    cols_masking = cols_masking & ~ np.all(np.isnan(covariances[0]), axis = 0)
    masked_rows_covariance = covariances[0][cols_masking, :][:, cols_masking]
    rows_precision_matrix = np.linalg.inv(masked_rows_covariance)
    M_s_est = M_s_est[:, cols_masking]

    # Matrix of orthogonal projection onto the column space of 'M_s_est'
    # 'M_s_est' is well conditionned as based on the Legendre polynomials
    # Orthogonality could also ease constraining based on physical priors
    M_s_est_Gram_matrix = M_s_est @ rows_precision_matrix @ M_s_est.T
    Pi_M_s_est = M_s_est.T @ np.linalg.inv(M_s_est_Gram_matrix)

    # Estimation of the Legendre polynomial coefficients of the data decomposition
    coeffs_est = D[:, cols_masking] @ rows_precision_matrix @ Pi_M_s_est
    # Reset to 'NaN' of the data spaxels full of 'NaN'
    coeffs_est[all_nan_D_spaxels] = np.nan

    # Estimation of the PSF
    psf_est = coeffs_est @ M

    # Save of the Legendre coefficients
    # being a reduced estimate of the PSF
    if param.save_verbose:
        # Coefficients normalization (for homogeneous comparisons)
        coeffs_normalization = np.linalg.norm(D[:, cols_masking], axis = 1)[:, None]
        coeffs_normalization[np.isnan(coeffs_normalization) | (coeffs_normalization == 0)] = 1
        coeffs_est /= coeffs_normalization
        # Add of the previous coefficients if they exist
        # (to avoid overloading 'outputs' with too many files)
        if os.path.exists(param.outputs_path + "_PSF_rows_.npy"):
            coeffs_est += np.load(param.outputs_path + "_PSF_rows_.npy")
        # Save of these new coefficients representing the reduced PSF
        np.save(param.outputs_path + "_PSF_rows_.npy", coeffs_est)

    return psf_est

def point_spread_function_estimate_S5(D: float_arr_type, s_est: float_arr_type,
                                      world_grids: list[tuple[list[int], float_arr_type]],
                                      covariances: tuple_float_arr_type) -> float_arr_type:
    """ PSF estimation - Stellar Spread Subtraction (S5) version
    ---

    Estimation of the Point Spread Function matrix using the S5 method, \\
    that is to say by factorizing the PSF matrix with two Legendre polynomial matrix, \\
    on the features / spectral side (using the assumption of PSF smoothness in this direction), \\
    but also on the observations side (e.g. assuming the PSF smoothness in the radial direction)

    - The regularization that this factorization causes allows a much robust (yet biaised) \\
    estimation of the PSF with regard to the noise (or even to any other source in the data)
    - This comes to estimate the stellar component of the data with a polynomial modulation \\
    of the stellar spectrum estimate (in a robust way thanks to the Legendre Basis), \\
    but also with a polynomial modulation of the PSF in the observations axes direction \\
    (e.g. in the radial direction if the axes are the relative right ascensions / declinaisons)

    With this factorization and the assumption of a Gaussian noise, \\
    the maximum likehood estimate of the Point Spread Function matrix becomes:

        N (N^T C_v N)^-1 N C_v D C_w Ms (Ms C_w Ms^T)^-1 Ms

    with `Ms` the matrix of Legendre polynomial modulation of the stellar spectrum estimate, \\
    and `C_w` the covariance matrix of the noise in the rows direction, \\
    as well as `N` the design matrix of the Legendre polynomials, \\
    and `C_v` the covariance matrix of the noise in the cols direction

    Args:
        D (np.typing.NDArray[np.float64]):
            Data matrix from which to estimate the Point Spread Function
        s_est (np.typing.NDArray[np.float64]):
            Stellar spectrum estimate by mean of the selected pixels \\
            (with the `stellar_spectrum_estimate` function)
        covariances (tuple[np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]):
            Rows and columns noise covariances matrices estimates \\
            (with the `noise_covariance_estimate` function) \\
            or `np.nan` to indicate unestimated covariances
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the grid of the corresponding axes, in the same access order

    Returns:
        psf_est (np.typing.NDArray[np.float64]):
            Point Spread Function estimate with the S5 method
    """
    
    # Legendre polynomial fit of the data rows
    psf_est = point_spread_function_estimate_S4(**locals())

    # Non-cols grids construction
    # (averaging along correlated non-cols axes)
    cols_grids_axes = [] # Liste of the axis grids
    for axis_index in range(len(world_grids) - 1):
        non_cols_grid_axes = list(range(len(world_grids[axis_index][0])))
        del non_cols_grid_axes[world_grids[axis_index][0].index(axis_index)]
        cols_grid_axis = np.mean(world_grids[axis_index][1].T, axis = tuple(non_cols_grid_axes))
        cols_grids_axes.append(cols_grid_axis)
    # Euclidean merge of the grids of the different cols axes
    merged_grids = np.meshgrid(*cols_grids_axes, indexing = 'ij')
    euclidian_merged_grids = np.sqrt(sum(grid**2 for grid in merged_grids))
    cols_grid = euclidian_merged_grids.reshape(-1, 1) # Vectorization as the columns axes

    # Shift of the 'cols_grid' in the Legendre [-1,1] range
    legendre_grid = 2 * (cols_grid - min(cols_grid)) / (max(cols_grid) - min(cols_grid)) - 1

    # Matrix of the Legendre polynomials
    # by following Bonnet's recursion formula
    # -> (n+1)Pn+1(x) = (2n+1)xPn(x) - nPn-1(x)
    lower_columns_dimension = param.psf_degrees[1] + 1
    N = np.tile(legendre_grid, (1, lower_columns_dimension))
    N[:, 0] = 1 # Degree 0 (total flux) Legendre polynomial
    for j in range(2, lower_columns_dimension):
        N[:, j] *= N[:, j-1] * (2 * j - 1) / j
        N[:, j] -= N[:, j-2] * (j - 1) / j

    # Data 'D' matrix cols full of 'NaN' masking
    rows_masking = ~ np.all(np.isnan(D), axis = 1)

    # Scaling of all other cols (regressors) to unit norms
    # so that each one of them weight equally in the regression
    regressors_norms = np.linalg.norm(N[rows_masking, 1:], axis = 0)
    N[:, 1:] /= regressors_norms

    # Data 'D' matrix 'NaN' values zeroing
    all_nan_D_images = np.all(np.isnan(D), axis = 0)
    D[np.isnan(D)] = 0 # For the matrix multiplications
    
    # Precision matrix weighting
    if covariances[1] is np.nan: covariances[1] = np.identity(D.shape[0])
    rows_masking = rows_masking & ~ np.all(np.isnan(covariances[1]), axis = 0)
    masked_cols_covariance = covariances[1][rows_masking, :][:, rows_masking]
    rows_precision_matrix = np.linalg.inv(masked_cols_covariance)
    N = N[rows_masking]

    # Matrix of orthogonal projection onto the column space of 'N'
    # 'N' is well conditionned as based on the Legendre polynomials
    # Orthogonality could also ease constraining based on physical priors
    N_Gram_matrix = N.T @ rows_precision_matrix @ N
    Pi_N = np.linalg.inv(N_Gram_matrix) @ N.T

    # Estimation of the PSF
    all_nan_psf_spaxels = np.all(np.isnan(psf_est), axis = 1)
    psf_est = N @ Pi_N[:, ~ all_nan_psf_spaxels] @ psf_est[ ~ all_nan_psf_spaxels]
    psf_est[all_nan_psf_spaxels, :] = np.nan
    psf_est[:, all_nan_D_images] = np.nan

    # Save of the Legendre coefficients
    # being a reduced estimate of the PSF
    if param.save_verbose:
        # Estimation of the Legendre polynomial coefficients of the data decomposition
        coeffs_est = Pi_N @ rows_precision_matrix @ D[rows_masking]
        # Reset to 'NaN' of the data spaxels full of 'NaN'
        coeffs_est[:, all_nan_D_images] = np.nan
        # Coefficients normalization (for homogeneous comparisons)
        coeffs_normalization = np.linalg.norm(D[rows_masking], axis = 0)
        coeffs_normalization[np.isnan(coeffs_normalization) | (coeffs_normalization == 0)] = 1
        coeffs_est /= coeffs_normalization
        # Add of the previous coefficients if they exist
        # (to avoid overloading 'outputs' with too many files)
        if os.path.exists(param.outputs_path + "_PSF_cols_.npy"):
            coeffs_est += np.load(param.outputs_path + "_PSF_cols_.npy")
        # Save of these new coefficients representing the reduced PSF
        np.save(param.outputs_path + "_PSF_cols_.npy", coeffs_est)

    return psf_est

def stellar_component_estimate(psf_est: float_arr_type, s_est: float_arr_type) -> float_arr_type:
    """ Stellar component estimation
    ---

    Estimation of the stellar component of the data following the model,\\
    by multiplication of the PSF estimate with the stellar spectrum estimate

    Args:
        psf_est (np.typing.NDArray[np.float64]):
            Point Spread Function matrix estimate with the S3, S4 or S5 method \\
            (i.e. with one of the functions `point_spread_function_estimate_S3`, \\
            `point_spread_function_estimate_S4` or `point_spread_function_estimate_S5` \\)
        s_est (np.typing.NDArray[np.float64]):
            Stellar spectrum estimate by mean of the selected pixels \\
            (with the `stellar_spectrum_estimate` function)

    Returns:
        S_est (np.typing.NDArray[np.float64]):
            Stellar component estimate
    """

    S_est = psf_est * s_est

    return S_est

def noise_covariance_estimate(E_est: float_arr_type) -> tuple[float_arr_type, float_arr_type]:
    """ Rows and columns noise covariances matrices estimation
    ---

    Estimation of the covariance matrices `C` from `np.cov`, \\
    then with optimal shrinkage towards the corresponding diagonal matrix with:

        tr(C²) + tr²(C) - 2 tr(C) ☉ tr(C)
        ---------------------------------- (Flasseur et al. 2018a)
          | C | (tr(C²) - tr(C) ☉ tr(C))

    with `tr` the trace operator, `☉` the entry-wise Hadamard product  \\
    and `| . |` the number of samples used for the estimation of C

    Args:
        E_est (np.typing.NDArray[np.float64]):
            Noise component matrix estimate by subtraction of the stellar component matrix estimate

    Raises:
        ValueError:
            Impossible shrinkage if the rows covariance matrix is only 0 outside of the diagonal \\
            (and in any case, not only does weighting become unnecessary, \\
            but it is also very likely a sign that something has gone wrong)
        ValueError:
            Impossible shrinkage if the cols covariance matrix is only 0 outside of the diagonal \\
            (and in any case, not only does weighting become unnecessary, \\
            but it is also very likely a sign that something has gone wrong)

    Returns:
        noise_covariances (tuple[np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]]):
            Rows and columns noise covariances matrices estimates
    """

    E_est[np.isnan(E_est)] = 0

    n_rows, n_cols = E_est.shape

    # Rows covariance matrix

    C_rows = np.cov(E_est, rowvar = False)

    squared_trace_C_rows = np.trace(C_rows) **2
    trace_squared_C_rows = np.trace(C_rows @ C_rows)
    sum_diag_squared_C_rows = np.sum(np.diag(C_rows)**2)

    if trace_squared_C_rows - sum_diag_squared_C_rows == 0:
        error_message = "Zero covariances: shrinkage impossible"
        error_message += " (and noise statistics weighting thus useless anyway)"
        error_message += " > Something wrong probably happened during this estimation loop"
        raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

    rows_shrinkage = squared_trace_C_rows + trace_squared_C_rows - 2*sum_diag_squared_C_rows
    rows_shrinkage /= (n_cols + 1) * (trace_squared_C_rows - sum_diag_squared_C_rows)

    C_rows = (1 - rows_shrinkage) * C_rows + rows_shrinkage * np.diag(np.diag(C_rows))

    # Columns covariance matrix

    C_cols = np.cov(E_est, rowvar = True)

    squared_trace_C_cols = np.trace(C_rows) **2
    trace_squared_C_cols = np.trace(C_rows @ C_rows)
    sum_diag_squared_C_cols = np.sum(np.diag(C_rows)**2)

    if trace_squared_C_cols - sum_diag_squared_C_cols == 0:
        error_message = "Zero covariances: shrinkage impossible"
        error_message += " (and noise statistics weighting thus useless anyway)"
        error_message += " > Something wrong probably happened during this estimation loop"
        raise ValueError("\x1B[31m" + error_message + "\x1B[0m")

    cols_shrinkage = squared_trace_C_cols + trace_squared_C_cols - 2*sum_diag_squared_C_cols
    cols_shrinkage /= (n_rows + 1) * (trace_squared_C_cols - sum_diag_squared_C_cols)

    C_cols = (1 - cols_shrinkage) * C_cols + cols_shrinkage * np.diag(np.diag(C_cols))

    # Covariance matrices grouping

    noise_covariances = (C_rows, C_cols)

    return noise_covariances

################################################################################################### Tools

def masking_tool(data_tensor: float_arr_type, axis_names: list[str],
                 world_grids: list[tuple[list[int], float_arr_type]],
                 intervals_to_be_masked: mask_type) -> float_arr_type:
    """ Masking tool - replacements with NaN values
    ---

    Replacement of pixel values of a data tensor with the NaN value (to mask these pixels), \\
    on the intersection of different given intervals along different given axes

    **Note:** If no interval is specified, `data_tensor` is returned as is \\
    thanks to `max(1, np.max(pixels_to_be_masked)))` \\
    instead of `np.max(pixels_to_be_masked)`

    Args:
        data_tensor (np.typing.NDArray[np.float64]):
            Data tensor for which the given intervals have to masked
        axis_names (list[str]):
            Names of the different axes as extracted in the data header
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the grid of the corresponding axes, in the same access order
        intervals_to_be_masked (dict[str, list[tuple[int | float, int | float]]]):
            List of intervals between which pixels are replaced by NaN values after intersection \\
            for each given axis (as long as its name is found in the header via `axis_names`)

    Returns:
        data_tensor (np.typing.NDArray[np.float64]):
            Masked data tensor, i.e. with NaN at the intersection of the given intervals
    """

    pixels_to_be_masked = np.zeros_like(data_tensor) # To identify the intersection

    for axis_name, intervals_bounds in intervals_to_be_masked.items():
        if axis_name not in axis_names: continue
        axis_index = axis_names.index(axis_name)

        non_axis_grid_axes = list(range(len(world_grids[axis_index][0])))
        del non_axis_grid_axes[world_grids[axis_index][0].index(axis_index)]
        axis_grid = np.mean(world_grids[axis_index][1].T, axis = tuple(non_axis_grid_axes))

        for lower_bound, upper_bound in intervals_bounds:

            if isinstance(lower_bound, float):
                if axis_grid[0] < axis_grid[-1]:
                    lower_bound = np.searchsorted(axis_grid, lower_bound, side = 'right') - 1
                else : lower_bound = np.searchsorted(-axis_grid, -lower_bound, side = 'left')
            if isinstance(upper_bound, float):
                if axis_grid[0] < axis_grid[-1]:
                    upper_bound = np.searchsorted(axis_grid, upper_bound, side = 'left')
                else : upper_bound = np.searchsorted(-axis_grid, -upper_bound, side = 'right') - 1

            if lower_bound > upper_bound: lower_bound, upper_bound = upper_bound, lower_bound

            slicer = pixels_to_be_masked.ndim * [slice(None)]
            slicer[axis_index] = slice(lower_bound, upper_bound)
            pixels_to_be_masked[tuple(slicer)] += 1

    masks_intersection = (pixels_to_be_masked == max(1, np.max(pixels_to_be_masked)))

    data_tensor[masks_intersection] = np.nan

    return data_tensor

def filtering_tool(data_tensor: float_arr_type, axis_names: list[str],
                   world_grids: list[tuple[list[int], float_arr_type]],
                   loc_scale_filter_lists: mask_type) -> float_arr_type:
    """ Filtering tool - multiplication with unit-height Gaussians
    ---

    Multiplication of the pixel values of a data tensor with a filter built as \\
    the multiplication along the different axes of the sums of the Gaussians \\
    given by `loc_scale_filter_lists` through their locations and scales

    **Note:** The Gaussians used being unit-height Gaussians, \\
    it is possible to filter some areas more or less than others \\
    by adding the same Gaussian location-scale multiple times

    Args:
        data_tensor (np.typing.NDArray[np.float64]):
            Data tensor for which the pixels values have to be multiplied \\
            with Gaussians at the given locations and scales parameters
        axis_names (list[str]):
            Names of the different axes as extracted in the data header
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the grid of the corresponding axes, in the same access order
        loc_scale_filter_lists (mask_type):
            List of locations and scales of Gaussians to sum to build a global filter \\
            for each given axis (as long as its name is found in the header via `axis_names`)

    Returns:
        data_tensor (np.typing.NDArray[np.float64]):
            Filtered data tensor, i.e. mutliplied with the products of the sums\\
            of unitary-height Gaussians at the given locations and scales
    """

    gaussian = lambda grid, loc, scale: np.exp(-1/2 * ((grid-loc) / scale)**2)

    global_filter = np.ones_like(data_tensor)

    for axis_name, loc_scale_filter_list in loc_scale_filter_lists.items():
        if axis_name not in axis_names: continue
        axis_index = axis_names.index(axis_name)

        non_axis_grid_axes = list(range(len(world_grids[axis_index][0])))
        del non_axis_grid_axes[world_grids[axis_index][0].index(axis_index)]
        axis_grid = np.mean(world_grids[axis_index][1].T, axis = tuple(non_axis_grid_axes))

        axis_filter = np.zeros_like(axis_grid)

        for loc_filter, scale_filter in loc_scale_filter_list:

            if isinstance(loc_filter, int):
                loc_filter = axis_grid[loc_filter]
            if isinstance(scale_filter, int):
                if axis_grid[0] < axis_grid[-1]:
                    loc_filter_index = np.searchsorted(axis_grid, loc_filter)
                else: loc_filter_index = np.searchsorted(axis_grid, loc_filter)
                scale_filter_right = axis_grid[loc_filter_index + scale_filter]
                scale_filter_left = axis_grid[loc_filter_index - scale_filter]
                scale_filter = np.abs(scale_filter_right - scale_filter_left) / 2

            axis_filter += gaussian(axis_grid, loc_filter, scale_filter)

        shape = [1] * data_tensor.ndim
        shape[axis_index] = data_tensor.shape[axis_index]

        if axis_filter.any():
            global_filter *= axis_filter.reshape(shape)

    data_tensor *= global_filter

    return data_tensor

def hindering_tool(data_tensor: float_arr_type, axis_names: list[str],
                   world_grids: list[tuple[list[int], float_arr_type]],
                   loc_scale_hinder_lists: mask_type) -> float_arr_type:
    """ Hindering tool - multiplication with reversed unit-height Gaussians
    ---

    Multiplication of the pixel values of a data tensor with a hinder built from unit-arrays \\
    minus the multiplication along the different axes of the sums of the Gaussians \\
    given by `loc_scale_hinder_lists` through their locations and scales

    **Note:** The Gaussians used being unit-height Gaussians, the add of a unit-array \\
    multiplied by the maximum height reached by the global reversed hinder \\
    allows to avoid any negative values in the global hinder obtained by \\
    opposite of the global reversed hinder obtained by sums products \\
    (this can be used to hinder some areas more or less than others \\
    by adding the same Gaussian location-scale multiple times)

    Args:
        data_tensor (np.typing.NDArray[np.float64]):
            Data tensor for which the pixels values have to be multiplied \\
            with Gaussians at the given locations and scales parameters
        axis_names (list[str]):
            Names of the different axes as extracted in the data header
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the grid of the corresponding axes, in the same access order
        loc_scale_filter_lists (mask_type):
            List of locations and scales of Gaussians to sum to build a global reversed hinder \\
            for each given axis (as long as its name is found in the header via `axis_names`), \\
            itself subtracted from the unit-arrays to build the global hinder ultimately used

    Returns:
        data_tensor (np.typing.NDArray[np.float64]):
            Hindered data tensor, i.e. mutliplied with from unit-arrays minus \\
            the products of the sums of unitary-height Gaussians at the given locations and scales
    """

    gaussian = lambda grid, loc, scale: np.exp(-1/2 * ((grid-loc) / scale)**2)

    global_reversed_hinder = np.ones_like(data_tensor)

    for axis_name, loc_scale_hinder_list in loc_scale_hinder_lists.items():
        if axis_name not in axis_names: continue
        axis_index = axis_names.index(axis_name)

        non_axis_grid_axes = list(range(len(world_grids[axis_index][0])))
        del non_axis_grid_axes[world_grids[axis_index][0].index(axis_index)]
        axis_grid = np.mean(world_grids[axis_index][1].T, axis = tuple(non_axis_grid_axes))

        axis_reversed_hinder = np.zeros_like(axis_grid)

        for loc_hinder, scale_hinder in loc_scale_hinder_list:

            if isinstance(loc_hinder, int):
                loc_hinder = axis_grid[loc_hinder]
            if isinstance(scale_hinder, int):
                if axis_grid[0] < axis_grid[-1]:
                    loc_hinder_index = np.searchsorted(axis_grid, loc_hinder)
                else: loc_hinder_index = np.searchsorted(axis_grid, loc_hinder)
                scale_hinder_right = axis_grid[loc_hinder_index + scale_hinder]
                scale_hinder_left = axis_grid[loc_hinder_index - scale_hinder]
                scale_hinder = np.abs(scale_hinder_right - scale_hinder_left) / 2

            axis_reversed_hinder += gaussian(axis_grid, loc_hinder, scale_hinder)

        shape = [1] * data_tensor.ndim
        shape[axis_index] = data_tensor.shape[axis_index]

        if axis_reversed_hinder.any():
            global_reversed_hinder *= axis_reversed_hinder.reshape(shape)

    global_hinder = np.max(global_reversed_hinder)*np.ones_like(global_reversed_hinder)
    global_hinder -= global_reversed_hinder

    if global_hinder.any(): data_tensor *= global_hinder

    return data_tensor

################################################################################################### Show

def show_dematricized(dematricized_data: float_arr_type,
                      pixels_axis_name: str, pixels_axis_unit: str,
                      axis_names: list[str], axis_units: list[str],
                      world_grids: list[tuple[list[int], float_arr_type]],
                      suptitle: str = "", color: str = 'k', cmap: str = 'binary') -> None:
    """ Dematricized data shows - images and spectra means
    ---

    Show of a mean spectrum (along the last axis) and of a mean image (along the two first axes) \\
    of the dematricized data, for monitoring purposes, with labels/units extracted from the header

    Args:
        dematricized_data (np.typing.NDArray[np.float64]):
            Dematricized data matrix for which a mean spectrum and a mean image must be shown
        pixels_axis_name (str):
            Name of the pixel axis values as extracted in the data header \\
            (or specifically modified depending on the type of data displayed)
        pixels_axis_unit (str):
            Unit of the pixel axis values as extracted in the data header \\
            (or specifically modified depending on the type of data displayed)
        axis_names (list[str]):
            Names of the different axes as extracted in the data header
        axis_units (list[str]):
            Units of the different axes as extracted in the data header \\
            or as choosen in `param.zooms_units`
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the grid of the corresponding axes, in the same access order
        suptitle (str, optional):
            Title surmounting the spectrum and the image display \\
            Defaults to ""
        color (str, optional):
            Color of the spectrum line \\
            Defaults to 'k' (for black)
        cmap (str, optional):
            Color map of the image cells \\
            Defaults to 'binary' (for white and black)
    """
    
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

    # Copy of `dematricized_data` to avoid overwriting it
    to_plot, to_show = dematricized_data.copy(), dematricized_data.copy()
    # Masking, filtering and hindering for spectra plots
    to_plot = masking_tool(to_plot, axis_names, world_grids, param.spectra_masker)
    to_plot = filtering_tool(to_plot, axis_names, world_grids, param.spectra_filter)
    to_plot = hindering_tool(to_plot, axis_names, world_grids, param.spectra_hinder)
    # Masking, filtering and hindering for images shows
    to_show = masking_tool(to_show, axis_names, world_grids, param.images_masker)
    to_show = filtering_tool(to_show, axis_names, world_grids, param.images_filter)
    to_show = hindering_tool(to_show, axis_names, world_grids, param.images_hinder)
    # Separate processing for spectra to plot and images to show
    dematricized_data_to_plot, dematricized_data_to_show = to_plot, to_show

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
    plt.suptitle(suptitle, fontsize = 24)
    # Fullscreen (depending on the graphics engine - backend - used by matplotlib)
    manager, backend = plt.get_current_fig_manager(), plt.get_backend().lower()
    if 'tkagg' in backend:
        try: manager.window.state('zoomed')
        except: manager.window.wm_state('zoomed')
    elif 'qt' in backend: manager.window.showMaximized()
    elif "wx" in backend: manager.frame.Maximize(True)
    elif 'gtk' in backend: manager.window.maximize()
    elif backend == 'macosx': pass
    else: pass

    # World axes labels determination
    plot_xlabel_name = axis_names[-1]
    plot_xlabel_unit = axis_units[-1]
    if plot_xlabel_name in param.zooms_units.keys():
        plot_xlabel_unit = param.zooms_units[plot_xlabel_name]
    show_xlabel_name = axis_names[0]
    show_xlabel_unit = axis_units[0]
    if show_xlabel_name in param.zooms_units.keys():
        show_xlabel_unit = param.zooms_units[show_xlabel_name]
    show_ylabel_name = axis_names[1]
    show_ylabel_unit = axis_units[1]
    if show_ylabel_name in param.zooms_units.keys():
        show_ylabel_unit = param.zooms_units[show_ylabel_name]
    
    # Pixels axis labels determination
    try: pixels_axis_unit = Unit(pixels_axis_unit).to_string(format = "unicode")
    except ValueError: pass # An unrecognized unit (as "without unit") is left as is

    # Colorbar saturation
    if not np.all(np.isnan(pcolor_image)):
        xy_anomalies_sigmas = [0,0]
        if show_xlabel_name in param.anomalies_sigmas.keys():
            xy_anomalies_sigmas[0] = param.anomalies_sigmas[show_xlabel_name]
        if show_ylabel_name in param.anomalies_sigmas.keys():
            xy_anomalies_sigmas[1] = param.anomalies_sigmas[show_ylabel_name]
        if xy_anomalies_sigmas[0] or xy_anomalies_sigmas[1]:
            vmax = np.nanmean(pcolor_image) + max(xy_anomalies_sigmas)*np.nanstd(pcolor_image)
        else: vmax = np.nanmax(pcolor_image)
    else: vmax = None

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
    plt.xlabel(f"{show_xlabel_name} ({show_xlabel_unit})", fontsize = 18)
    plt.ylabel(f"{show_ylabel_name} ({show_ylabel_unit})", fontsize = 18)
    plt.pcolor(x_grid, y_grid, pcolor_image.T, vmin = 0, vmax = vmax, cmap = cmap)
    clrbar = plt.colorbar(fraction = 0.025, aspect = 50)
    clrbar_label = f"{pixels_axis_name} ({pixels_axis_unit})"
    clrbar.set_label(clrbar_label, fontsize = 18, rotation = 270, labelpad = 36)

    # Figure show
    plt.show()

###################################################################################################
