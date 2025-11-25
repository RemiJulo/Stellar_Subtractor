import os
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

import param

sep = os.path.sep # '/' for Linux/macOS or '\' for Windows

###################################################################################################

outputs_folder = param.outputs_path
hdu_data_merge, hdu_header_merge = None, None

fits_files = [f for f in os.listdir(outputs_folder) if f.endswith(".fits")]

if not fits_files:
    exit(f"No .fits found in the {param.outputs_path} folder")

for output_file in fits_files:
    file_path = os.path.join(outputs_folder, output_file)
    try:
        with fits.open(file_path) as hdul:
            for hdu_index, hdu in enumerate(hdul):
                data, header = hdu.data, hdu.header
                if isinstance(data, np.ndarray):
                    if hdu_data_merge is None:
                        hdu_data_merge = data.copy()
                        hdu_header_merge = header.copy()
                    else:
                        if data.shape != hdu_data_merge.shape:
                            raise ValueError(f"The data haven't all the same shape")
                        hdu_data_merge += data
                    for key in header:
                        if key not in hdu_header_merge:
                            hdu_header_merge[key] = header[key]
    except Exception as e:
        print(f"Error while reading the file {output_file}: {e}")
        continue

if hdu_data_merge is not None:
    hdu = fits.PrimaryHDU(hdu_data_merge, header = hdu_header_merge)
    output_path = os.path.join(outputs_folder, "Fusion.fits")
    hdu.writeto(output_path, overwrite = True)
