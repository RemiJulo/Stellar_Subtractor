import matplotlib.pyplot as plt
import numpy as np
import warnings
import os

from os import path
from astropy.io import fits
from astropy.wcs import WCS

import process
import param

###################################################################################################

from astropy.wcs import FITSFixedWarning
warnings.simplefilter('ignore', FITSFixedWarning)

################################################################################################### Merging

if not path.exists("add"): os.mkdir("add")

outputs = param.outputs_path

hdu_data = None
hdu_header = None

print("\nMerging outputs:\n")

for fits_file in [f for f in os.listdir(outputs) if f.endswith(".fits")]:
    hdul = fits.open(path.join(outputs, fits_file))
    print(f"  {fits_file}")

    for hdu in hdul:
        data = hdu.data
        header = hdu.header

        if isinstance(data, np.ndarray):

            if hdu_header is None:
                hdu_header = header.copy()
                hdu_data = data.copy()

            else:
                if data.shape != hdu_data.shape:
                    error_message = f"{data.shape} new shape"
                    error_message += f" (instead of {hdu_data.shape})"
                    raise ValueError("\x1B[31m" + error_message + "\x1B[0m")
                
                hdu_header.update(header)
                hdu_data += data

    hdul.close()

if hdu_header is not None:
    wcs = WCS(header)

    hdu = fits.PrimaryHDU(hdu_data, hdu_header)
    hdu.writeto(path.join("add", "Estimation.fits"), overwrite = True)

################################################################################################### Whitening

print("\n> Whitening\n")

coordinate_types = [axis_type["coordinate_type"] for axis_type in wcs.get_axis_types()]
if coordinate_types.count('spectral') != 1: raise ValueError("No spectral axis")
hdu_data = np.moveaxis(hdu_data.T, coordinate_types.index('spectral'), -1)

hdu_data = process.rank1_standardization(hdu_data, *process.rank1_moments(hdu_data, True))
hdu_data = np.moveaxis(hdu_data, -1, coordinate_types.index('spectral')).T

hdu = fits.PrimaryHDU(hdu_data, hdu_header)

hdu.writeto(path.join("add", "Detection.fits"), overwrite = True)

################################################################################################### LOP coeffs maps

psf_list = [np.load(path.join(outputs, d, "PSF.npy")) for d in os.listdir(outputs)
            if path.isdir(path.join(outputs, d)) and path.isfile(path.join(outputs, d, "PSF.npy"))]

if len(psf_list):
    psf = np.median(psf_list, axis = 0)
    np.save(path.join("add", "PSF.npy"), psf)

    print("---\n\nLOP decomposition\n")

    plt.figure("Legendre Orthogonal Polynomials coefficients maps", figsize = (16,9))
    plt.suptitle("Legendre Orthogonal Polynomials coefficients maps\n", fontsize = 20)

    legendre_coeffs_energies = np.zeros(psf.shape[-1])
    for degree, coeffs_map in enumerate(np.moveaxis(psf, -1, 0)):
        legendre_coeffs_energies[degree] = np.quantile(coeffs_map**2, 0.75)
        
        if coeffs_map.ndim > 2:
            gt2 = tuple(axis for axis in range(coeffs_map.ndim) if axis > 2)
            coeffs_map = np.nansum(coeffs_map, axis = gt2)
        if coeffs_map.ndim < 2: coeffs_map = coeffs_map[:, None]

        vmin = np.quantile(coeffs_map, 0.25)
        vmax = np.quantile(coeffs_map, 0.75)

        if degree != 0 :
            vmax = np.max((np.max(-vmin), np.max(vmax)))
            vmin = -np.max((np.max(-vmin), np.max(vmax)))
        
        plt.subplot(2, (psf.shape[-1] + 1) // 2, degree + 1)

        plt.imshow(coeffs_map.T, origin = 'lower', cmap = 'bwr', vmin = vmin, vmax = vmax)
        
        plt.colorbar(orientation = 'horizontal')
        plt.title(f"Degree {degree}", fontsize = 15)
        plt.xlabel(wcs.axis_type_names[1], fontsize = 10)
        plt.ylabel(wcs.axis_type_names[0], fontsize = 10)
        plt.xticks([])
        plt.yticks([])

    print("    LOP_maps.pdf", end = '\r')
    plt.savefig(path.join("add", "LOP_maps.pdf"), dpi = 150)

################################################################################################### LOP coeffs energies

    plt.figure("Legendre Orthogonal Polynomials coefficients energies", figsize = (16,9))
    plt.title("Legendre Orthogonal Polynomial coefficients energies\n", fontsize = 20)

    degrees, energies = np.arange(1, len(legendre_coeffs_energies)), legendre_coeffs_energies[1:]
    markerline, stemlines, baseline = plt.stem(degrees, energies, linefmt = 'b', markerfmt = 'r')

    markerline.set_markersize(15)
    stemlines.set_linewidth(5)
    baseline.set_visible(0)

    elbow = np.argmax(np.diff(np.diff(np.log(energies)))) + 1 + 1 + 1
    plt.scatter(degrees[:elbow], energies[:elbow], 500, linewidths = 2, fc = 'none', ec = 'red')
    
    plt.xlabel("Legendre Orthogonal Polynomial degree\n", fontsize = 15)
    plt.ylabel("Legendre Orthogonal Polynomial energy\n", fontsize = 15)
    plt.grid(True, which = "major", linestyle = '-')
    plt.grid(True, which = "minor", linestyle = ':')
    plt.xticks(degrees, fontsize = 10)
    plt.yticks(fontsize = 10)
    plt.yscale('log')

    print("    LOP_energies.pdf", end = '\r')
    plt.savefig(path.join("add", "LOP_energies.pdf"), dpi = 150)

###################################################################################################
