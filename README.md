# STELLAR SUBTRACTOR (S2)

## > Overview

### Prerequisites

The following libraries must be installed:
- `numpy`
- `astropy`
- `matplotlib`

This project is currently developed with Python 3.12.11.

### List of files

- `main.py` is the only file to execute
- `param.py` is the only file to modify
- `process.py` is the only process file

### Documentation

Analytical and on-sky demonstration:  
https://arxiv.org/pdf/2509.09878

### Contact

Feel free to contact me with any question at `remi.julo@univ-grenoble-alpes.fr`

### Much more to come!

## > First use

1. Modify `param.py`

    [x] Assign `input_paths` a sub-list of paths to .fits files (or to folders full of .fits files) to be processed  
    [x] Make sure that `save_verbose`, `print_verbose` and `show_verbose` are assigned to `1` to monitor execution  
    [x] Choose the folder names of the backup path `saving_folders` in which all results will be saved

2. Execute `main.py`

    [x] Either open a terminal and type `<path_to_your_Python_version> main.py`  
    [x] Or open `main.py` with your favorite IDE to execute it

3. Refine `param.py`

    [x] Adapt the `zooms` and `origin` parameters based on textual, visual, and backup feedback  
    [x] Change any other relevant parameters (see **Advanced parameters** section for wise choices)

4. Reexecute `main.py` on a sub-list of paths to .fits files (or to folders full of .fits files) to be processed

    [x] Same as 2.

5. Clean `param.py`

    [x] If the returns are still not satisfactory, go back to step 3  
    [x] Else, put `print_verbose` and `show_verbose` to `0` to speed up processing

6. Reexecute `main.py` on the full-list of paths to .fits files (or to folders full of .fits files) to be processed

    [x] Same as 2. and 4.  
    [x] Open `output_path`  
    [x] Enjoy your post-processing 😎

## > Advanced parameters

Upcoming...
