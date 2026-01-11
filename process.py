import matplotlib.pyplot as plt
import numpy as np

import warnings
import os

import param

###################################################################################################

type tensor_type = np.typing.NDArray[np.float64]
type w_grids_type = list[tuple[list[int], tensor_type]]
type boolean_type = np.bool | np.typing.NDArray[np.bool]
type moment_type = list[int | float, tensor_type, tensor_type]
type mask_type = list[dict[str, tuple[int | float, int | float]]]

################################################################################################### Algorithm

def process(data_tensor: tensor_type, world_grids: w_grids_type,
            axis_names: list[str], axis_units: list[str]) -> tensor_type:
    """ Stellar Subtraction - General algorithm (SLOPES)
    ---

    The **S.L.O.P.E.S.** acronym stands for \\
    **S**tellar **L**egendre **O**rthogonal **P**olynomials **E**stimation **S**ubtractor \\
    and symbolizes both the *slopes* of the chromatic PSFs \\
    and the algorithm's alpine origin

    - Outliers detection (for their disregard in the post-processing)
    - Start of the estimation loop (with increasingly better noise estimates)
        - Stellar spectrum estimation (by averaging the selected least noisy pixels)
        - Point Spread Function estimation (through decomposition on Legendre basis)
        - Stellar component estimation (from the stellar spectrum and PSF estimates)
    - End of the estimation loop (in practice after two or even just one iteration)
    - Anomalies detection (for their highlight in the post-subtraction report)
    
    Args:
        data_tensor (np.typing.NDArray[np.float64]):
            Data from which to subtract the stellar component \\
            The last axis of this tensor must be its only spectral axis
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values along each of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the joint grid of the corresponding axes, in the same axis access order
        axis_names (list[str]):
            Names of the different axes as extracted in the data header \\
            The spectral axis entry must be the panultimate one \\
            (brightness axis entry being the ultimate one)
        axis_units (list[str]):
            Units of the different axes as extracted in the data header \\
            The spectral axis entry must be the panultimate one \\
            (brightness axis entry being the ultimate one)
    Returns:
        subtracted_tensor (np.typing.NDArray[np.float64]):
            Data tensor from which a stellar component estimate has been subtracted
    """

    b_axis_name = axis_names.pop() # The brightness axis name and unit are at
    b_axis_unit = axis_units.pop() # the end of the axis names and units lists

    if param.show_verbose:
        
        if data_tensor.ndim > 3:
            warning_message = "Too large number of dimension for image display:"
            warning_message += "Truncation of the tensor at the first two axes"
            print("\x1B[33m" + '\n\n' + warning_message + '\n\n' + "\x1B[0m")

        images_grid = joint_grid([0, 1], world_grids)
    spectra_grid = joint_grid([-1], world_grids)

    if param.show_verbose:

        data_show(data_tensor, world_grids, spectra_grid, images_grid,
                  axis_names, axis_units, b_axis_name, b_axis_unit,
                  "\nRaw pre-processed data\n", 'b', 'plasma')

    standard = rank1_standardization(data_tensor, *rank1_moments(data_tensor, True))
    standard = masking(standard, axis_names, world_grids, param.outliers_mask)
    outliers = abnormal_standardization(standard, param.outliers_sigma)

    if param.show_verbose:

        data_show(standard, world_grids, spectra_grid, images_grid,
                  axis_names, axis_units, b_axis_name, b_axis_unit,
                  "\nOutliers in pre-processed data\n", 'k', 'hot')
    
    noise_est = np.full_like(data_tensor, np.nan)
    for loopback in range(1 + param.noise_loopback):

        if param.print_verbose and param.noise_loopback:
            print("\x1B[1;35m" + f"\nLoop {loopback + 1}" + "\x1B[0m")
        
        noise_est = estimation(data_tensor, world_grids, spectra_grid, images_grid,
                               axis_names, axis_units, b_axis_name, b_axis_unit,
                               outliers, param.star_mask, param.psf_mask,
                               *rank1_moments(noise_est, loopback))

    if param.show_verbose:
        
        std_noise = rank1_standardization(noise_est, *rank1_moments(noise_est, True))
        std_noise = masking(std_noise, axis_names, world_grids, param.anomalies_mask)
        noise_est[~ abnormal_standardization(std_noise, param.anomalies_sigma)] = 0

        data_show(std_noise, world_grids, spectra_grid, images_grid,
                  axis_names, axis_units, b_axis_name, b_axis_unit,
                  "\nAnomalies in subtraction\n", 'g', 'viridis')

    return noise_est

def estimation(data_tensor: tensor_type, world_grids: w_grids_type,
               spectra_grid: tensor_type, images_grid: tensor_type,
               axis_names: list[str], axis_units: list[str],
               b_axis_name: str, b_axis_unit: str,
               mask: boolean_type, star_mask: mask_type, psf_mask: mask_type,
               means: moment_type, stds: moment_type) -> tensor_type:
    """ Stellar Subtraction - Estimation loop
    ---

    Estimation of the noise component from the estimation of the stellar one

    - Stellar spectrum estimation (by averaging the selected least noisy pixels)
    - Point Spread Function estimation (through decomposition on Legendre basis)
    - Stellar component estimation (from the stellar spectrum and PSF estimates)
    
    Args:
        data_tensor (np.typing.NDArray[np.float64]):
            Data from which to subtract the stellar component \\
            The last axis of this tensor must be its only spectral axis
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values along each of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the joint grid of the corresponding axes, in the same axis access order
        axis_names (list[str]):
            Names of the different axes as extracted in the data header \\
            The spectral axis entry must be the panultimate one \\
            (brightness axis entry being the ultimate one)
        axis_units (list[str]):
            Units of the different axes as extracted in the data header \\
            The spectral axis entry must be the panultimate one \\
            (brightness axis entry being the ultimate one)
        b_axis_name (str):
            Name of the brightness axis as extracted in the data header
        b_axis_unit (str):
            Unit of the brightness axis as extracted in the data header
        mask (bool | np.typing.NDArray[np.bool]):
            Array of pixels to be masked (replacing them with 'NaN' values), \\
            for stellar spectrum and Point Spread Function estimates, \\
            or 'False' to avoid masking anything pixel-by-pixel
        star_mask (list[dict[str, tuple[int | float, int | float]]]):
            Dictionnaries of intervals along multiples axes for which to mask \\
            the pixels of the intersections (replacing them with 'NaN' values) \\
            for stellar spectrum (but not Point Spread Function) estimate, \\
            or '[]' to avoid masking anything interval-by-interval
        psf_mask (list[dict[str, tuple[int | float, int | float]]]):
            Dictionnaries of intervals along multiples axes for which to mask \\
            the pixels of the intersections (replacing them with 'NaN' values) \\
            for Point Spread Function (but not stellar spectrum) estimate, \\
            or '[]' to avoid masking anything interval-by-interval
        means (list[int | float, np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]]):
            Estimation of the order 1 moments from rank 1 assumption \\
            in the form of a list containing a common factor, the vectorized non-spectral part, \\
            and lastly the spectral part (to multiply tensorially to reconstruct the rank 1 moment)
        stds (list[int | float, np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]]):
            Estimation of the square roots of the centered order 2 moments \\
            in the form of a list containing a common factor, the vectorized non-spectral part, \\
            and lastly the spectral part (to multiply tensorially to reconstruct the rank 1 moment)

    Returns:
        E_est (np.typing.NDArray[np.float64]):
            Estimate of the noise component of the given data tensor
    """
    
    D = masking(data_tensor, world_grids, axis_names, star_mask, mask)

    s_est = stellar_spectrum_estimate(D)

    if param.show_verbose:

        data_show(D, world_grids, spectra_grid, images_grid,
                  axis_names, axis_units, b_axis_name, b_axis_unit,
                  "\nStellar spectrum estimation\n", 'm', 'inferno')

    D = masking(data_tensor, world_grids, axis_names, psf_mask, mask)
    
    if param.show_verbose:

        out = np.full_like(D, np.nan)
        where = (s_est != 0) & ~ np.isnan(s_est)
        ratio = np.divide(D, s_est, where = where, out = out)

        data_show(ratio, world_grids, spectra_grid, images_grid,
                  axis_names, axis_units, b_axis_name, b_axis_unit,
                  "\nUnregularized PSF estimation\n", 'y', 'cividis')

    psf_est = psf_estimate(D, s_est, spectra_grid, means, stds)

    if param.show_verbose:

        data_show(psf_est, world_grids, spectra_grid, images_grid,
                  axis_names, axis_units, b_axis_name, b_axis_unit,
                  "\nRegularized PSF estimation\n", 'y', 'cividis')

    S_est = stellar_component_estimate(s_est, psf_est)

    if param.show_verbose:
        
        data_show(S_est, world_grids, spectra_grid, images_grid,
                  axis_names, axis_units, b_axis_name, b_axis_unit,
                  "\nStellar component estimation\n", 'r', 'magma')

    E_est = data_tensor - S_est

    if param.show_verbose:

        data_show(E_est, world_grids, spectra_grid, images_grid,
                  axis_names, axis_units, b_axis_name, b_axis_unit,
                  "\nPost stellar subtraction data\n", 'c', 'cool')

    return E_est

################################################################################################### Estimation

def stellar_spectrum_estimate(D: tensor_type) -> tensor_type:
    """ Stellar spectrum estimation
    ---
    
    Estimation of the spectrum of the main (stellar) component of the data, \\
    by median of the pixels whose spectral median flux is greater than the threshold \\
    given by `param.stellar_level` (by comparing each of these fluxes with the maximum of them)

    **Note:** Despite the stellar spectrum normally being given by total flux integration, \\
    the median of the brightest pixels is preferred in this function for greater robustness, \\
    knowing that there are multiplicative degenerations with the chromatic PSFs in any case

    Args:
        D (np.typing.NDArray[np.float64]):
            Data tensor from which to estimate the stellar spectrum

    Returns:
        s_est (np.typing.NDArray[np.float64]):
            Stellar spectrum estimate by median of the selected pixels
    """

    if param.print_verbose:
        print("\x1B[1;35m" + "\n    Stellar spectrum estimation" + "\x1B[0m")

    # Stellar positions selection
    star_flux_max = np.nanmax(np.nanmedian(D, axis = -1))
    noise_flux_level = star_flux_max / param.star_fluxlevel
    selection = noise_flux_level <= np.nanmedian(D, axis = -1)

    # Stellar spectrum estimation
    s_est = np.nanmedian(D[selection], axis = 0)

    # Masking for shows
    if param.show_verbose: D[~ selection] = np.full(D.shape[-1], np.nan)

    return s_est

def psf_estimate(D: tensor_type, s_est: tensor_type, spectra_grid: tensor_type,
                 means: moment_type, stds: moment_type) -> tensor_type:
    """ Point Spread Function estimation
    ---

    Estimation of the Point Spread Function matrix, \\
    by factorization of the PSF matrix with a Legendre polynomial matrix, \\
    on the features / spectral side (using the assumption of PSF smoothness in this direction)

    - The regularization that this factorization causes allows a much robust (yet biaised) \\
    estimation of the PSF with regard to the noise (or even to any other source in the data)
    - This comes to estimate the stellar component of the data with a polynomial modulation \\
    of the stellar spectrum estimate (in a robust way thanks to the Legendre Basis)

    With this factorization and the assumption of a Gaussian noise, \\
    the maximum likehood estimate of the Point Spread Function matrix becomes:

        (D - M1) M2^-1 Ls (Ls M2^-1 Ls^T)^-1 L

    with `Ls` the Legendre Orthogonal Polynomials of modulation of the stellar spectrum estimate \\
    `M1` the mean of the noise and `M2` the covariance of the noise in the spectral direction

    Args:
        D (np.typing.NDArray[np.float64]):
            Data matrix from which to estimate the Point Spread Function
        s_est (np.typing.NDArray[np.float64]):
            Stellar spectrum estimate by median of the selected pixels \\
            (with the `stellar_spectrum_estimate` function)
        spectra_grid (np.typing.NDArray[np.float64]):
            Grid of the world (physical) values along the last (spectral) axis
        means (list[int | float, np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]]):
            Estimation of the order 1 moments from rank 1 assumption \\
            in the form of a list containing a common factor, the vectorized non-spectral part, \\
            and lastly the spectral part (to multiply tensorially to reconstruct the rank 1 moment)
        stds (list[int | float, np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]]):
            Estimation of the square roots of the centered order 2 moments \\
            in the form of a list containing a common factor, the vectorized non-spectral part, \\
            and lastly the spectral part (to multiply tensorially to reconstruct the rank 1 moment)

    Returns:
        psf_est (np.typing.NDArray[np.float64]):
            Point Spread Function estimate
    """

    if param.print_verbose:
        print("\x1B[1;35m" + "\n    Point Spread Function estimation" + "\x1B[0m")

    # Masking of the 'NaN' and '0' pixels of 's_est'
    spectra_mask = ~ (np.isnan(s_est) | (s_est == 0))
    # Masking of the 'NaN' pixels of the spectral std
    spectra_mask &= ~ np.all(np.isnan(stds[2]), axis = 0)
    # Masking of the spectral channels of 'D' full of 'NaN'
    non_spectra_axes = tuple(axis for axis in range(D.ndim - 1))
    spectra_mask &= ~ np.all(np.isnan(D), axis = non_spectra_axes)

    L, L_s_est = legendre_modulation_matrix(spectra_grid[0], s_est, spectra_mask)

    # Data 'D' matrix 'NaN' values zeroing
    all_nan_D_spaxels = np.all(np.isnan(D), axis = -1)
    D[np.isnan(D)] = 0 # For the matrix multiplications
    
    # Data 'D' matrix 'NaN' values zeroing
    rank1_mean = means[0] * np.outer(means[1], means[2])
    spectral_cov_inv = np.linalg.inv(np.diag(stds[2]))

    # Masking and standardization of the data 'D' and Legendre matrix 'L_s_est'
    D = (D - rank1_mean.reshape(D.shape))[..., spectra_mask] @ spectral_cov_inv[spectra_mask]
    L_s_est_tilde = L_s_est[:, spectra_mask] @ spectral_cov_inv[spectra_mask, :][:, spectra_mask]

    # Generalized inverse 'L_s_est_plus' of the Legendre matrix 'L_s_est'
    # 'L_s_est' is well conditionned as based on Legendre Orthogonal Polynomials
    L_s_est_plus = L_s_est_tilde.T @ np.linalg.inv(L_s_est_tilde @ L_s_est_tilde.T)

    # Estimation of the PSF
    psf_est = np.tensordot(D[..., spectra_mask], L_s_est_plus @ L , axes = ([-1], [0]))
    # Reset to 'NaN' of the data spaxels full of 'NaN'
    psf_est[all_nan_D_spaxels] = np.nan

    # Save of the Legendre coefficients
    # being a reduced estimate of the PSF
    if param.save_verbose:

        # Estimation of the Legendre Orthogonal Polynomial coefficients of the data decomposition
        coeffs_est = np.tensordot(D[..., spectra_mask], L_s_est_plus, axes = ([-1], [0]))
        # Reset to 'NaN' of the data spaxels full of 'NaN'
        coeffs_est[all_nan_D_spaxels] = np.nan

        # Coefficients normalization (for homogeneous comparisons)
        coeffs_normalization = np.linalg.norm(D[..., spectra_mask], axis = -1)[..., None]
        coeffs_normalization[np.isnan(coeffs_normalization) | (coeffs_normalization == 0)] = 1
        # Save of the normalized coefficients representing the reduced PSF
        np.save("PSF.npy", coeffs_est / coeffs_normalization)

    return psf_est

def stellar_component_estimate(s_est: tensor_type, psf_est: tensor_type) -> tensor_type:
    """ Stellar component estimation
    ---

    Estimation of the stellar component of the data following the model, \\
    by multiplication of the PSF estimate with the stellar spectrum estimate

    Args:
        s_est (np.typing.NDArray[np.float64]):
            Stellar spectrum estimate by median of the selected pixels \\
            (with the `stellar_spectrum_estimate` function)
        psf_est (np.typing.NDArray[np.float64]):
            PSF estimate by polynomial modulation of the stellar spectrum estimate \\
            (with the `psf_estimate` and `legendre_modulation_matrix` functions)

    Returns:
        S_est (np.typing.NDArray[np.float64]):
            Stellar component estimate
    """

    if param.print_verbose:
        print("\x1B[1;35m" + "\n    Stellar component estimation" + "\x1B[0m")

    S_est = psf_est * s_est

    return S_est

################################################################################################### Factorization

def legendre_modulation_matrix(world_grid: tensor_type, s_est: tensor_type,
                               mask: boolean_type) -> tuple[tensor_type, tensor_type]:
    """ Matrix of the Legendre Orthogonal Polynomials
    ---

    Matrix containing the successive Legendre polynomials on its different rows \\
    and its row-by-row multiplication by the stellar spectrum estimate, \\
    in order to modulate this latter robustly thanks to \\
    the orthogonality of the built basis

    Args:
        world_grid (np.typing.NDArray[np.float64]):
            Grid shifted in the [-1,1] range for the Legendre polynomials sampling
        s_est (np.typing.NDArray[np.float64]):
            Stellar spectrum estimate by median of the selected pixels \\
            (with the `stellar_spectrum_estimate` function)
        mask (np.bool | np.typing.NDArray[np.bool]):
            Spectral channels of the `s_est` spectrum to disregard \\
            during the matrices scaling and orthogonalization

    Returns:
        L (np.typing.NDArray[np.float64]):
            Matrix of the Legendre Orthogonal Polynomials
        L_s_est (np.typing.NDArray[np.float64]):
            Matrix of the Legendre Orthogonal Polynomials \\
            of modulation of the stellar spectrum estimate `s_est`
    """
    
    # Shift of the 'world_grid' in the Legendre [-1,1] range
    legendre_grid = 2 * (world_grid - min(world_grid)) / (max(world_grid) - min(world_grid)) - 1

    # Matrix of the Legendre polynomials
    # by following Bonnet's recursion formula
    # -> (n+1)Pn+1(x) = (2n+1)xPn(x) - nPn-1(x)
    lower_spectral_dimension = param.psf_polydegree + 1
    L = np.tile(legendre_grid, (lower_spectral_dimension, 1))
    L[0] = 1 # Degree 0 (total flux) Legendre polynomial
    for i in range(2, param.psf_polydegree + 1):
        L[i] *= L[i-1] * (2 * i - 1) / i
        L[i] -= L[i-2] * (i - 1) / i

    # Matrix of the Legendre polynomial modulation
    # of the stellar spectrum estimate 's_est'
    L_s_est = L * s_est

    # Scaling of the first row (regressor) to unit norm
    # so that it weights as the others in the regression
    # (it is 's_est' as the 1st legendre polynomial is 1)
    regressor_0_norm = np.linalg.norm(s_est[mask])
    L_s_est[0] /= regressor_0_norm
    L[0] /= regressor_0_norm
    # Orthogonalization of variation rows with the (first) total flux row
    # (inevitably fully correlated because of the multiplication with 's_est')
    # so that it comes '<L_s_est_ortho[i],L_s_est[0]> = L_s_est_ortho[i].T.s_est
    # is equal to '(L_s_est[i]-(<L_s_est[i],s_est>/<s_est,s_est>) s_est).T.s_est'
    # and 'L_s_est[i].T.s_est-(L_s_est[i].T.s_est/s_est.T.s_est) s_est.T.s_est'
    # and 'L_s_est[i].T.s_est-L_s_est[i].T.s_est = 0' hence the orthogonality
    # 'L' is subtracted in a way it keeps the same relation with 'L_s_est'
    colinearity = L_s_est[1:, mask] @ s_est[mask].T
    normalization_s_est = s_est[mask] @ s_est[mask].T
    normalized_colinearity = colinearity / normalization_s_est
    L_s_est[1:] -= normalized_colinearity[:, None] * s_est
    L[1:] -= normalized_colinearity[:, None]
    # Scaling of all other rows (regressors) to unit norms
    # so that each one of them weight equally in the regression
    regressors_norms = np.linalg.norm(L_s_est[1:, mask], axis = 1)
    L_s_est[1:] /= regressors_norms[:, None]
    L[1:] /= regressors_norms[:, None]

    return L, L_s_est

################################################################################################### Statistics

def rank1_moments(data_tensor: tensor_type, svd: bool | int) -> tuple[moment_type, moment_type]:
    """ Rank1 means and standard deviations estimation
    ---

    Estimation, with rank 1 assumptions, of the means (order 1 moments) \\
    and standard deviations (square roots of the centered order 2 moments)

    **Note:** To avoid the extremely costly calculation of all data pixel covariances \\
    (followed by rank 1 "non-spectral x spectral" decomposition by SVD as for the mean), \\
    only pixel variances are considered (inter-pixel covariances being close to zero anyway), \\
    approximated empirically axis by axis (despite the inaccuracy by a factor that this means), \\
    and via their square roots directly (being the ones used in the standardization calculations)

    Args:
        data_tensor (np.typing.NDArray[np.float64]):
            Data from which to estimate the rank 1 moments \\
            The last axis of this tensor must be its only spectral axis
        svd (bool | int):
            Activation of rank 1 moments calculation from SVD, \\
            on data vectorized along the non-spectral axes \\
            (default zero and identity matrices else)

    Returns:
        means (list[int | float, np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]]):
            Estimation of the order 1 moments from rank 1 assumption \\
            in the form of a list containing a common factor, the vectorized non-spectral part, \\
            and lastly the spectral part (to multiply tensorially to reconstruct the rank 1 moment)
        stds (list[int | float, np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]]):
            Estimation of the square roots of the centered order 2 moments \\
            in the form of a list containing a common factor, the vectorized non-spectral part, \\
            and lastly the spectral part (to multiply tensorially to reconstruct the rank 1 moment)
    """

    data_tensor_shape = data_tensor.shape

    if svd:
        data_tensor[np.isnan(data_tensor)] = 0

        vectors = data_tensor.reshape(-1, data_tensor_shape[-1])
        left_sv, sv, right_sv_T = np.linalg.svd(vectors, full_matrices = False)

        means = [sv[0], left_sv[:, 0], right_sv_T[0, :]]
        vectors -= means[0] * np.outer(means[1], means[2])

        stds = [1, np.std(vectors, axis = 1), np.std(vectors, axis = 0)]

    else:
        means = [1, np.zeros(np.prod(data_tensor_shape[:-1])), np.zeros(data_tensor_shape[-1])]
        stds = [1, np.ones(np.prod(data_tensor_shape[:-1])), np.ones(data_tensor_shape[-1])]

    return means, stds

def rank1_standardization(data_tensor: tensor_type,
                          means: moment_type, stds: moment_type) -> tensor_type:
    """ Data standardization with rank 1 mean and rank 1 standard deviation
    ---

    Subtraction of a rank 1 mean estimate and division by a rank 1 standard deviation estimate \\
    (the rank 1 separation being between the non-spectral and spectral axes)

    Args:
        data_tensor (np.typing.NDArray[np.float64]):
            Data tensor from which to standardized the spaxels \\
            The last axis of this tensor must be its (only) spectral axis
        means (list[int | float, np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]]):
            Estimation of the order 1 moments from rank 1 assumption \\
            in the form of a list containing a common factor, the vectorized non-spectral part, \\
            and lastly the spectral part (to multiply tensorially to reconstruct the rank 1 moment)
        stds (list[int | float, np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]]):
            Estimation of the square roots of the centered order 2 moments \\
            in the form of a list containing a common factor, the vectorized non-spectral part, \\
            and lastly the spectral part (to multiply tensorially to reconstruct the rank 1 moment)

    Returns:
        standardized_data (np.typing.NDArray[np.float64]):
            Data tensor from which the spaxels have been standardized with rank 1 moments
    """

    rank1_mean = means[0] * np.outer(means[1], means[2]).reshape(data_tensor.shape)
    rank1_std = stds[0] * np.outer(stds[1], stds[2]).reshape(data_tensor.shape)
    
    out = np.full_like(data_tensor, np.nan)
    where = (rank1_std != 0) & ~ np.isnan(rank1_std)

    standardized_data = np.divide(data_tensor - rank1_mean, rank1_std, out = out, where = where)

    return standardized_data

def abnormal_standardization(data_tensor: tensor_type,
                             sigmas: int | float) -> boolean_type:
    """ Abnormal pixels extraction
    ---

    Extraction of abnormal pixels, whether along the spectral or non-spectral axes, \\
    abnormalily being when a value deviates from the mean of the pixel values \\
    by more standard deviations of the pixel values than fixed by `sigmas`

    **Note:** The non-abnormal pixels are directly filled with 'NaN' values

    Args:
        data_tensor (np.typing.NDArray[np.float64]):
            Data tensor from which to identify the abnormal pixels
        sigmas (int | float):
            Numbers of standard deviations above the mean above which the values of the pixels \\
            classifie them as abnormal (leading to a `True` entry at their position)
            
    Returns:
        abnormal_pixels (bool | np.typing.NDArray[np.bool]):
            Mask identifying the abnormal pixels
    """

    non_spectral = tuple(axis for axis in range(data_tensor.ndim - 1))

    data_tensor[..., np.all(np.isnan(data_tensor), axis = non_spectral)] = 0
    non_spectral_deviations = np.abs(data_tensor - np.nanmean(data_tensor, axis = non_spectral))
    abnormal_pixels = non_spectral_deviations > sigmas*np.nanstd(data_tensor, axis = non_spectral)
    data_tensor[..., np.all(data_tensor == 0, axis = non_spectral)] = np.nan

    data_tensor[np.all(np.isnan(data_tensor), axis = -1)] = 0
    spectral_deviations = np.abs(data_tensor - np.nanmean(data_tensor, axis = -1)[..., None])
    abnormal_pixels &= spectral_deviations > sigmas*np.nanstd(data_tensor, axis = -1)[..., None]
    data_tensor[np.all(data_tensor == 0, axis = -1)] = np.nan

    data_tensor[~ abnormal_pixels] = np.nan

    return abnormal_pixels

################################################################################################### Tools

def masking(data_tensor: tensor_type, world_grids: w_grids_type, axis_names: list[str],
            intervals: mask_type = [], pixels: boolean_type = False) -> tensor_type:
    """ Masking tool - replacements with 'NaN' values
    ---

    Replacement of pixel values of a data tensor with the 'NaN' value, \\
    on the intersections of different given intervals along different given axes

    Args:
        data_tensor (np.typing.NDArray[np.float64]):
            Data tensor for which the given intervals have to masked
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values along each of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the joint grid of the corresponding axes, in the same axis access order
        axis_names (list[str]):
            Names of the different axes as extracted in the data header
        intervals (list[dict[str, tuple[int | float, int | float]]], optionnal):
            Dictionnaries of intervals along multiples axes for which to mask \\
            the pixels of the intersections (replacing them with 'NaN' values), \\
            or '[]' by default to avoid masking anything interval-by-interval
        pixels (bool | np.typing.NDArray[np.bool], optionnal):
            Array of pixels to be masked (replacing them with 'NaN' values), \\
            or 'False' by default to avoid masking anything pixel-by-pixel

    Returns:
        masked_data_tensor (np.typing.NDArray[np.float64]):
            Masked data tensor, i.e. with 'NaN' at the given pixels and intervals
    """

    masked_data_tensor = data_tensor.copy()

    for interval in intervals:

        intersections = np.zeros_like(data_tensor)

        for axis_name, (lower_bound, upper_bound) in interval.items():

            if axis_name not in axis_names: continue
            axis_index = axis_names.index(axis_name)

            axis_grid = joint_grid([axis_index], world_grids)[0]
           
            if isinstance(lower_bound, float):
                if axis_grid[0] < axis_grid[-1]:
                    lower_bound = np.searchsorted(axis_grid, lower_bound, side = 'right') - 1
                else : lower_bound = np.searchsorted(-axis_grid, -lower_bound, side = 'left')
            if isinstance(upper_bound, float):
                if axis_grid[0] < axis_grid[-1]:
                    upper_bound = np.searchsorted(axis_grid, upper_bound, side = 'left')
                else : upper_bound = np.searchsorted(-axis_grid, -upper_bound, side = 'right') - 1

            if lower_bound > upper_bound: lower_bound, upper_bound = upper_bound, lower_bound

            slicer = intersections.ndim * [slice(None)]
            slicer[axis_index] = slice(lower_bound, upper_bound)
            intersections[tuple(slicer)] += 1 # To identify the intersection

        mask_intersections = (intersections == max(1, np.max(intersections)))

        masked_data_tensor[mask_intersections] = np.nan

    masked_data_tensor[pixels] = np.nan

    return masked_data_tensor

def filtering(data_tensor: tensor_type, world_grids: w_grids_type,
              axis_names: list[str], filters: mask_type) -> tensor_type:
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
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values along each of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the joint grid of the corresponding axes, in the same axis access order
        axis_names (list[str]):
            Names of the different axes as extracted in the data header
        filters (mask_type):
            List of locations and scales of Gaussians to sum to build a global filter \\
            for each given axis (as long as its name is found in the header via `axis_names`)

    Returns:
        filtered_data_tensor (np.typing.NDArray[np.float64]):
            Filtered data tensor, i.e. mutliplied with the products of the sums\\
            of unitary-height Gaussians at the given locations and scales
    """

    gaussian = lambda grid, loc, scale: np.exp(-1/2 * ((grid-loc) / scale)**2)

    filtered_data_tensor = data_tensor.copy()

    for filter in filters:

        global_filter = np.ones_like(data_tensor)

        for axis_name, (loc_filter, scale_filter) in filter.items():

            if axis_name not in axis_names: continue
            axis_index = axis_names.index(axis_name)

            axis_grid = joint_grid([axis_index], world_grids)[0]
           
            if isinstance(loc_filter, int):
                loc_filter = axis_grid[loc_filter]
            if isinstance(scale_filter, int):
                if axis_grid[0] < axis_grid[-1]:
                    loc_filter_index = np.searchsorted(axis_grid, loc_filter)
                else: loc_filter_index = np.searchsorted(axis_grid, loc_filter)
                scale_filter_right = axis_grid[loc_filter_index + scale_filter]
                scale_filter_left = axis_grid[loc_filter_index - scale_filter]
                scale_filter = np.abs(scale_filter_right - scale_filter_left) / 2

            shape = [1] * data_tensor.ndim
            shape[axis_index] = data_tensor.shape[axis_index]

            global_filter *= gaussian(axis_grid, loc_filter, scale_filter).reshape(shape)

        filtered_data_tensor *= global_filter

    return filtered_data_tensor

def hindering(data_tensor: tensor_type, world_grids: w_grids_type,
                   axis_names: list[str], hinders: mask_type) -> tensor_type:
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
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values along each of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the joint grid of the corresponding axes, in the same axis access order
        axis_names (list[str]):
            Names of the different axes as extracted in the data header
        hinders (mask_type):
            List of locations and scales of Gaussians to sum to build a global reversed hinder \\
            for each given axis (as long as its name is found in the header via `axis_names`), \\
            itself subtracted from the unit-arrays to build the global hinder ultimately used

    Returns:
        hindered_data_tensor (np.typing.NDArray[np.float64]):
            Hindered data tensor, i.e. mutliplied with from unit-arrays minus \\
            the products of the sums of unitary-height Gaussians at the given locations and scales
    """

    gaussian = lambda grid, loc, scale: np.exp(-1/2 * ((grid-loc) / scale)**2)

    hindered_data_tensor = data_tensor.copy()

    for hinder in hinders:

        global_reversed_hinder = np.ones_like(data_tensor)

        for axis_name, (loc_hinder, scale_hinder) in hinder.items():

            if axis_name not in axis_names: continue
            axis_index = axis_names.index(axis_name)

            axis_grid = joint_grid([axis_index], world_grids)[0]
           
            if isinstance(loc_hinder, int):
                loc_hinder = axis_grid[loc_hinder]
            if isinstance(scale_hinder, int):
                if axis_grid[0] < axis_grid[-1]:
                    loc_hinder_index = np.searchsorted(axis_grid, loc_hinder)
                else: loc_hinder_index = np.searchsorted(axis_grid, loc_hinder)
                scale_hinder_right = axis_grid[loc_hinder_index + scale_hinder]
                scale_hinder_left = axis_grid[loc_hinder_index - scale_hinder]
                scale_hinder = np.abs(scale_hinder_right - scale_hinder_left) / 2

            shape = [1] * data_tensor.ndim
            shape[axis_index] = data_tensor.shape[axis_index]

            global_reversed_hinder *= gaussian(axis_grid, loc_hinder, scale_hinder).reshape(shape)

        global_hinder = np.max(global_reversed_hinder)*np.ones_like(global_reversed_hinder)
        global_hinder -= global_reversed_hinder

        if global_hinder.any(): hindered_data_tensor *= global_hinder

    return hindered_data_tensor

################################################################################################### Figures

def data_show(data_tensor: tensor_type, world_grids: w_grids_type,
              spectra_grid: tensor_type, images_grid: tensor_type,
              axis_names: list[str], axis_units: list[str],
              b_axis_name: str, b_axis_unit: str,
              title: str, color: str, cmap: str):
    """ Data show - images and spectra means
    ---
 
    Mean spectrum (along the last axis) and mean image (along the two first axes)

    Args:
        data_tensor (np.typing.NDArray[np.float64]):
            Data tensor for which to show a mean spectrum and a mean image
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values along each of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the joint grid of the corresponding axes, in the same axis access order
        spectra_grid (np.typing.NDArray[np.float64]):
            Grid of the world (physical) values along the last (spectral) axis
        images_grid (np.typing.NDArray[np.float64]):
            Grid of the world (physical) values along the two first axis
        axis_names (list[str]):
            Names of the different axes as extracted in the data header \\
            The spectral axis entry must be the ultimate one
        axis_units (list[str]):
            Units of the different axes as extracted in the data header \\
            The spectral axis entry must be the ultimate one
        b_axis_name (str):
            Name of the brightness axis as extracted in the data header
        b_axis_unit (str):
            Unit of the brightness axis as extracted in the data header
        title (str):
            Title surmounting the spectrum and the image \\
            It also names the matplotlib window
        color (str):
            Color of the spectrum line
        cmap (str):
            Color map of the image
    """
    
    fullscreen(title)

    plt.subplot(121)

    to_plot = data_tensor

    to_plot = masking(to_plot, world_grids, axis_names, param.spectra_masker)
    to_plot = filtering(to_plot, world_grids, axis_names, param.spectra_filter)
    to_plot = hindering(to_plot, world_grids, axis_names, param.spectra_hinder)

    non_spectra_axes = tuple(axis for axis in range(to_plot.ndim - 1))
    spectra_nan = np.all(np.isnan(to_plot), axis = non_spectra_axes)
    to_plot[..., spectra_nan] = 0 # To avoid any 'nanmean' warnings
    mean_spectrum = np.nanmean(to_plot, axis = non_spectra_axes)
    mean_spectrum[spectra_nan] = np.nan # 'NaN' images erasure

    plt.plot(spectra_grid[0], mean_spectrum, color = color, linewidth = 3)
    
    if spectra_grid[0][0] > spectra_grid[0][-1]: plt.gca().invert_xaxis()

    spectrum_show(axis_names, axis_units, b_axis_name, b_axis_unit)

    plt.subplot(122)

    to_pcolormesh = data_tensor

    to_pcolormesh = masking(to_pcolormesh, world_grids, axis_names, param.images_masker)
    to_pcolormesh = filtering(to_pcolormesh, world_grids, axis_names, param.images_filter)
    to_pcolormesh = hindering(to_pcolormesh, world_grids, axis_names, param.images_hinder)

    non_images_axes = tuple(axis for axis in range(2, to_pcolormesh.ndim))
    images_nan = np.all(np.isnan(to_pcolormesh), axis = non_images_axes)
    to_pcolormesh[images_nan, ...] = 0 # To avoid any nanmean warnings
    mean_image = np.nanmean(to_pcolormesh, axis = non_images_axes)
    mean_image[images_nan] = np.nan # 'NaN' spectra erasure 
    
    warnings.filterwarnings("ignore", # For grids nearly constant along some axes
    message = "The input coordinates to pcolormesh are interpreted as cell centers")

    plt.pcolormesh(images_grid[0], images_grid[1], mean_image, cmap = cmap, vmin = 0)

    if images_grid[0][0][0] > images_grid[0][0][-1]: plt.gca().invert_xaxis()
    if images_grid[1][0][0] > images_grid[1][0][-1]: plt.gca().invert_yaxis()
    
    image_show(axis_names, axis_units, b_axis_name, b_axis_unit)

    plt.show()

def spectrum_show(axis_names: list[str], axis_units: list[str],
                  b_axis_name: str, b_axis_unit: str):
    """ Spectrum labeling with names and units
    ---

    Label a spectrum with names and units as extracted in the data header \\
    (brightness unit is converted to 'unicode' format)

    Args:
        axis_names (list[str]):
            Names of the different axes as extracted in the data header \\
            The spectral axis entry must be the ultimate one
        axis_units (list[str]):
            Units of the different axes as extracted in the data header \\
            The spectral axis entry must be the ultimate one
        b_axis_name (str):
            Name of the brightness axis as extracted in the data header
        b_axis_unit (str):
            Unit of the brightness axis as extracted in the data header
    """

    x_axis_name = axis_names[-1]
    if x_axis_name in param.zooms_units.keys():
        x_axis_unit = param.zooms_units[x_axis_name]
    else: x_axis_unit = axis_units[-1]
    
    plt.xlabel(f"{x_axis_name} ({x_axis_unit})", fontsize = 18)
    plt.ylabel(f"{b_axis_name} ({b_axis_unit})", fontsize = 18)

    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)

def image_show(axis_names: list[str], axis_units: list[str],
                  b_axis_name: str, b_axis_unit: str):
    """ Image labeling with names and units
    ---

    Label an image with names and units as extracted in the data header \\
    (brightness unit is converted to 'unicode' format)

    **Note:** This also adds a small colorbar labeling the brightness

    Args:
        axis_names (list[str]):
            Names of the different axes as extracted in the data header \\
            The spectral axis entry must be the ultimate one
        axis_units (list[str]):
            Units of the different axes as extracted in the data header \\
            The spectral axis entry must be the ultimate one
        b_axis_name (str):
            Name of the brightness axis as extracted in the data header
        b_axis_unit (str):
            Unit of the brightness axis as extracted in the data header
    """

    x_axis_name = axis_names[0]
    if x_axis_name in param.zooms_units.keys():
        x_axis_unit = param.zooms_units[x_axis_name]
    else: x_axis_unit = axis_units[0]
    y_axis_name = axis_names[1]
    if y_axis_name in param.zooms_units.keys():
        y_axis_unit = param.zooms_units[y_axis_name]
    else: y_axis_unit = axis_units[1]
    
    plt.xlabel(f"{x_axis_name} ({x_axis_unit})", fontsize = 18)
    plt.ylabel(f"{y_axis_name} ({y_axis_unit})", fontsize = 18)

    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)

    cbar = plt.colorbar(fraction = 0.025, aspect = 50) # Small colorbar labeling the brightness
    cbar.set_label(f"{b_axis_name} ({b_axis_unit})", fontsize = 18, rotation = 270, labelpad = 36)

def fullscreen(title: str):
    """ Fullscreen figure initialization
    ---

    Fullscreen depending on the graphics engine / backend used by matplotlib

    **Note:** If fullscreen is impossible, the figure is just initialized in the usual way

    Args:
        title (str):
            Title of the figure
    """
    
    plt.figure(title)
    plt.suptitle(title, fontsize = 24)
    
    manager, backend = plt.get_current_fig_manager(), plt.get_backend().lower()

    try:

        if 'gtk' in backend: manager.window.maximize()
        elif 'wx' in backend: manager.frame.Maximize(True)
        elif 'qt' in backend: manager.window.showMaximized()
        elif 'tk' in backend: # Windows ('nt') or Linux/macOS
            if os.name == 'nt': manager.window.state('zoomed')
            match os.uname().sysname: # Linux or macOS ('Darwin')
                case 'Linux': manager.window.attributes('-zoomed', True)
                case 'Darwin': manager.window.attributes('-fullscreen', True)

    except Exception: pass

###################################################################################################

def joint_grid(grid_axes: list[int], world_grids: w_grids_type) -> tensor_type:
    """ Builds the joint grid of given axes
    ---

    Build the joint grid of the world (physical) values along the different given axes \\
    (the dimension of the returned grid is equal to the number of element of `grid_axes`)

    For example, `world_grids = [([0,1], x_grid_xy), ([0,1], y_grid_xy), ([2], z_grid_z)]` \\
    leads to `[x_grid_xy, y_grid_xy]` with `grid_axes = [0,1]` \\
    and to `[z_grid_z]` with `grid_axes = [2]`

    Args:
        grid_axes (list[int]):
            List of axes from which to construct the joint grid \\
            (-1 corresponds to the ultimate, -2 to the panultimate, ...)
        world_grids (list[tuple[list[int], np.typing.NDArray[np.float64]]]):
            Grids of the world (physical) values along each of the different axis, \\
            grouped by correlated axes so that each tuple contains a list of axis indexes \\
            followed by the joint grid of the corresponding axes, in the same axis access order

    Returns:
        joint_grid (np.typing.NDArray[np.float64]):
            Joint grid of the world (physical) values along each of the given axes
    """
    
    grid_axes = [axis if axis >= 0 else len(world_grids) + axis for axis in grid_axes]

    joint_grid = []

    for grid_axis in grid_axes:

        corr_axes, corr_grid = world_grids[grid_axis]

        non_grid_axes = [index for index, axis in enumerate(corr_axes) if axis not in grid_axes]
        axis_length = lambda axis: world_grids[axis][1].shape[world_grids[axis][0].index(axis)]
        expansion = tuple(1 if axis in corr_axes else axis_length(axis) for axis in grid_axes)
        
        axis_grid = np.mean(corr_grid.T, axis = tuple(non_grid_axes))
        expanded_axis_grid = np.tile(axis_grid, expansion)

        joint_grid.append(expanded_axis_grid)

    return joint_grid
