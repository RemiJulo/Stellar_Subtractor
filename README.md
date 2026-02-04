# Stellar Legendre Orthogonal Polynomials Estimation Subtractor (SLOPES)

## > Overview

### Prerequisites

This project is currently being  
developed with `Python` `3.12.10`,  
using only the 3 following libraries:

- `numpy` `2.3.3`
- `astropy` `7.1.0`
- `matplotlib` `3.10.6`

If you encounter a compatibility issue,  
please update to these verified versions.

### List of files

Each repository of this project only consists in 4 files:

- `process.py` is a processing file, it contains algorithms, estimators, tools and display functions.
- `param.py` is the file to modify, it contains inputs, outputs, zoom, estimation, and monitoring options.
- `main.py` is the file to execute, it reads the contents of the inputs files and prepare their further treatment.
- `add.py` is an additional file, it reads the contents of the outputs folders and show a few more advanced results.

The `process.py` file is designed to be specific to the *method*.  
The `param.py` file is designed to be specific to the *account*.  
The `main.py` file is designed to be specific to the *project*.  
The `add.py` file is designed to be specific to the *results*.

### Documentation

Analytical and on-sky demonstration:  
https://arxiv.org/pdf/2509.09878  
(publication in A&A).

### Contact

Feel free to contact me at  
`remi.julo@univ-grenoble-alpes.fr`  
with any question, request or even simple remark.

## > First use

* **Download the .py code files of the project**

If *Git* is installed on your computer, open a terminal in your working directory and run:

> git clone https://github.com/RemiJulo/Stellar_Subtractor

Else, follow https://github.com/RemiJulo/Stellar_Subtractor > `<> Code` > `Download ZIP`  
and unzip the obtained *Stellar_Subtractor-main.zip* archive at the location of your choice.

* **Verify that you comply with the prerequisites**

If you haven't already done so, you will need to install `Python`  
as well as the 3 different libraries indicated in ***Prerequisites***  
(directly on your computer or via a virtual environment).

Regarding the versions, you can always give it a try  
if yours are not too far from those specified in ***Prerequisites***   
(these are the verified development versions, not the minimal ones),  
but if you encounter any problems, first update both your libraries and `Python`  
(and second contact me as indicated in ***Contact*** so that I can correct the problem for everyone).

* **Choose .fits data to be post-processed**

These must contain at least one Header Data Unit with:  
*.data* whose type is *np.ndarray* (regardless of the dimensions)  
*.header* whose number of axis of type *'spectral'* is one (and only one)

If you do not have any available, you can consult the ***Step-by-step guides***  
for which demonstration data are available in each of the different cases (below).

* **Set the *inputs_paths* parameter**

Open `param.py` and set *inputs_paths* with the list of paths to .fits files to be post-processed  
or with the (list of) path(s) to the folder(s) in which are stored the .fits files to be post-processed.

Alternatively, create an "input/" folder at the same location as your .py code files  
and move each of your .fits data to be post-processed inside of it  
(keeping the default value *inputs_paths = ["inputs/"]*).

* **Set the *outputs_paths* parameter**

Open `param.py` and set *outputs_paths* with the name of the folder in which to save the results  
or with the path to the folder in which to save the results from the .py code files location.  
Else, all your post-processed data should appear in a folder named *"outputs/"*  
(as long as your keep the default value *outputs_path = "outputs/"*).

* **Enable all monitoring options**

Open `param.py` and make sure that:

> save_verbose = 1

> print_verbose = 1

> show_verbose = 1

This will allow you to monitor the execution progress to verify that everything is going well  
and retrieve the information needed for subsequent more detailed post-processing.

On particular, *show_verbose = 1* will show spectra and images in .py windows,  
while *print_verbose = 1* will print essential information in your terminal,  
and *save_verbose = 1* will save additional results in *outputs_path*.

* **Finally execute the code**

If your IDE allows it, open `main.py` and run it.  
Else, equivalently, open a terminal and run:

> < path_to_python > < path_to_main.py >

If your `Python` executable is in the `PATH` environment variable,  
and your terminal is running from the .py code files location, this should be:

> python3 main.py

To understand your results and refine your post-processing with additional parameters,  
refer to ***Step-by-step guides*** for concrete examples explained from start to finish.

## > Step-by-step guides

### The detection of PDS 70 b and c with the MUSE instrument

In this subsection, we propose an emblematic case study, i.e. the detection of the two protoplanets PDS 70 b and c  
in the VLT/MUSE observations of June 20, 2018, thanks to this new method and as presented in ***Documentation***.

This is also a good way to discover in more practical terms the range of parameters available with this code,  
as well as the variety of problems encounterable, in this case but also in each of your own applications.

In order to do so, your first task will be the download of the pre-processed data made available at:  
https://github.com/RemiJulo/Stellar_Subtractor/releases/tag/MUSE/  
(this consists in 6 cubes of ~ 1.9 Go each).  

---

#### First try: understanding the monitoring displays

---

Thereupon, if you have kept all the default settings in `param.py` (especially with all empty lists and dictionaries),  
then correctly followed the ***First use*** instructions (in particular by activating all monitoring options),  
you should therefore see the following messages in your terminal:

> **File opening: V1032_Cen_2018-06-20T00_02_15.fits**

During execution, this message will appear each time a new cube from the list to be post-processed actually begins to be post-processed,  
each time with the character string ***File opening:*** followed by the name of the cube actually post-processed.

> -> **RA** from 212.0442047603703 deg (pixel 0) to 212.0412477952566 deg (pixel 317) with **origin** at 212.04274024518807 deg (pixel 157)

> -> **DEC** from -41.39912414572433 deg (pixel 0) to -41.39696344529658 deg (pixel 306) with **origin** at -41.3980579271964 deg (pixel 151)

> -> **AWAV** from 6.0509604492188e-07 m (pixel 0) to 9.350960449218801e-07 m (pixel 2640) with **origin** at 4.7497104492188e-07 m (pixel -1041)

In this series of messages, each line corresponds to an axis (e.g. *3* here, since it is hyperspectral cubes that are being post-processed).  
It is very important to understand these ones and how to modify them in order to use all of the parameters that follow.  
In particular, each first keyword in bold will identify the corresponding axis for the following parameters  
(**RA** for Right ascension, **DEC** for DEClination, and **AWAV** for Air WAVelength here).

On each of these lines, 3 *floats* are displayed, each followed by an *int* in parentheses.  
The first 2 correspond to the lower and upper boundaries of the given cube along the corresponding axis,  
while the last correspond to the origin chosen for the corresponding grid axis (for all the parameters that follow),  
in absolute physical values, and in the unit indicated immediately after (e.g. *deg* for **RA** and **DEC** and *m* for **AWAV** here).  
This is the *zooms* and *origin* (as well as *zooms_units* and *origin_units*) parameters that will allow the modification of these values.  
The values in parentheses indicate the indexes of the corresponding pixels in the original cube (before the subsequent zooms and origin relocation).  

---

After that (potentially after a few seconds given the size of a MUSE cube),  
a .py window should appear to display a visualization of the data cube,  
with the mean spaxel on the left and the mean image on the right  
(the whites are the spaxels and images full of NaN values).

---

Then, the stellar spectrum estimation starts with the following message:

> **Stellar spectrum estimation**

This is the first of the two quantities to be estimated in order to construct the estimate of the stellar component of the data.  
This latter having been estimated, a new .py window opens with the stellar spectrum estimate on the left (to the number of spaxels used)  
and the selection of the given spaxels on the right (the estimate being the integration of the values of these pixels, wavelength by wavelength).

---

Then, the Point Spread Function estimation starts with the following message:

> **Point Spread Function estimation**

This is the second of the two quantities to be estimated to construct the estimate of the stellar component of the data.  
Even before its estimation begins, a first .py window opens to show its unregularized estimation  
(through the entrywise division of the data spaxels by the stellar spectrum estimate).  
The .py window showing its regularized estimation appears just after.

The unregularized estimation gives an idea of the *slopes* of the chromatic PSF but includes all data noise  
(as well as any other component resulting from another source in the field).

The regularized estimation, though, keeps only the lowest pseudo-frequencies of these *slopes*  
by decomposing them on the basis of the *Legendre orthogonal polynomials*  
and retaining only the degrees less than or equal to *psf_degree*.  
This can be seen too as a "robust" polynomial fit.

At this point, a .npy file in a folder named after the post-processed .fits file appears (at the location specified by *outputs_path*).  
It contains, for each of the *slopes*, the coefficients estimated for the retained axes of the basis of the *Legendre orthogonal polynomials*.  
These results, a little more technical, are formatted by `app.py` and explained in ***Advanced results*** (they enable a wise choice of *psf_degree*).

---

Then, the stellar component estimation starts with the following message:

> **Stellar component estimation**

This is with the two previous quantities that this component is estimated  
(through the entrywise multiplication of the *slopes* by the stellar spectrum estimate).  
After its estimation, another .py window display the stellar estimate mean spaxel and mean image.

---

Once all this has been done, all that remains is to subtract the stellar component estimate from the data  
to obtain the noise component estimate (a.k.a. the post-subtraction residuals) according to the crafted model.  
This is the mean spaxel and the mean image of this noise component estimate that is displayed in the last .py window.

By hypothesis, it should only be a centered Gaussian white noise: any anomaly can thus be interpreted as the presence of another source.  
Yet, for now, neither PDS 70 b nor c is emerging from the noise (brighter sources should already be doing so).  
A closer look to the parameters available  in `param.py` is therefore necessary.

---

#### Fine tuning: understanding the parameters features

---

To begin with, it is wise to first choose a more comfortable framework  
by reducing the size of the fied of view under consideration,  
in order to reduce the execution and display times.  
To do so, two steps can be taken there.

---

On the one hand, if the cubes are not correctly centered on the star (as it is the case for the PDS 70 star here),  
it is possible to change the origin of the grids to facilitate the definition of the zooms  
(as well as all the other definitions of parameters in the following).

By visual inspection, the position of the star seems to be approximatley at coordinates (212.04273 deg, -41.39805 deg).  
This first parameter then allows to enter these new origin values for each of the considered axes,  
by associating the keywords found above in the terminal with the corresponding *floats*.

> origin = {
> 
>     'RA': 212.04273,
>     'DEC': -41.39805
> 
> }

This second parameter allows to specify the units in which the quantities above have been entered.  
(incidentally, specifying the units here also change them in the terminal, even if *origin* is blank).  

> origin_units = {
> 
>     'RA': 'deg',
>     'DEC': 'deg'
> 
> }

Note that leaving *origin_units* blank would amounts to the same here since 'deg' is already the default unit in both cases.

Alternatively, although this has not been primarily designed for this purpose,  
it is also possible to enter *int* values to indicate pixel indexes:

> origin = {
> 
>     'RA': 158,
>     'DEC': 152
> 
> }

This avoids the potential rounding errors by specifying the exact coordinates of the brightest pixel.

> origin_units = {
> 
> 
> 
> }

The units are then no longer considered.

Finally, for this parameter (and only for this parameter), it is also possible to use the sentinel class *argmax()*  
to indicate the desire to extract the coordinates of the brightest pixel and to set them there.

> origin = {
> 
>     'RA': argmax(),
>     'DEC': argmax()
> 
> }

It is an undeniable gain in comfort, but it is still important to keep in mind that using  
pre-centered data is always a better choice (especially for data that are more disturbed),  
both to save a few seconds here and to use a potentially more sophisticated pre-centering.

> origin_units = {
> 
> 
> 
> }

Once again, the units are there no longer considered.

---

On the other hand, now that an origin has been clearly established,  
it is possible to set a zoom level of his choice, from this origin.

For example, in the case of the MUSE instrument,  
the correction radius of the adaptive optics being around 0.5",  
it is not expected to detect any faint companion beyond this separation.  
Only a small square around the position of the star, defined just above as the new origin, can thus be selected.

> zooms = {
> 
>     'RA': (0.6,-0.6),
>     'DEC': (-0.6,0.6),
> 
> }

In the **RA** case, the upper physical bound is smaller than the lower physical bound  
because this axis increases from right to left (in the opposite direction to the pixel indexes).  
That being said, the bounds are automatically reversed if identified to have been defined the wrong way.

> zooms_units = {
> 
>     'RA': 'arcsec',
>     'DEC': 'arcsec'
> 
> }

The unit change is important this time so that the values above are indeed considered in arcseconds and not in degrees.  
In addition, all values associated with **RA** and **DEC** in further dictionnaries will also be in arcseconds from now on.

Moreover, once again it is possible to select the zooms based on pixel indices rather than physical values.

> zooms = {
> 
>     'RA': (-20,20),
>     'DEC': (-20,20),
> 
> }

Consequently, note that if you want to select a field of view of 1 arcsecond on either side of the center,  
then you must enter *(-1.0,1.0)* and not *(-1,1)*, in which case only a field of 1 pixel  
on either side of the center would be selected.

> zooms_units = {
> 
>     'RA': 'arcsec',
>     'DEC': 'arcsec'
> 
> }

Importantly, it is necessary to leave the custom units for all the parameters that follow this time  
(if the entered unit is not supported, then the list of available units will automatically display).

---

It is best not to use pixels outside the adaptive optics correction radius, for example.
The choice of the *psf_degree* parameter is a trade-off between the complexity of the chromatic PSF and the correlation distance of the noise.

---

### The characterization of YSES 1 b with the UVES instrument

Upcoming...

## > Advanced results

Upcoming...
