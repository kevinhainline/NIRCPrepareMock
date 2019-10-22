# NIRCPrepareMock
Code for taking the JAGUAR Mock Catalog and preparing it for photo-z codes

### `JAGUAR_to_Subsample.py`
This script takes in the mock catalogs from JAGUAR, and then spits out a file
that has a subsample of the full objects, such that the subsample tries to 
equally sample redshift and mass space (instead of being weighted towards 
low redshift, low mass objects). It's run with a variety of flags.
		
```
usage: JAGUAR_to_Subsample.py [-h] -in INPUT_FOLDER [-iID INPUTIDS]
                              [-sfo SF_OUTPUT_FILENAME]
                              [-qo Q_OUTPUT_FILENAME] [-iIDf INPUTIDFILENAME]
                              [-rid RANDOMIZE_IDS]
                              [-nz NUMBER_PER_REDSHIFT_BIN] [-fi FILTERS_FILE]
                              [-co COMBINE_FILENAME] [-mf]

optional arguments:
  -h, --help            show this help message and exit
  -in INPUT_FOLDER, --inputfolder INPUT_FOLDER
                        Input Folder With Mock Catalog Files
  -iID INPUTIDS, --inputIDs INPUTIDS
                        Input ID Filename
  -sfo SF_OUTPUT_FILENAME, --sfoutputfilename SF_OUTPUT_FILENAME
                        SF Output Filename
  -qo Q_OUTPUT_FILENAME, --qoutputfilename Q_OUTPUT_FILENAME
                        Q Output Filename
  -iIDf INPUTIDFILENAME, --inputIDfilename INPUTIDFILENAME
                        Filename when using an ID numbers file?
  -rid RANDOMIZE_IDS, --randomizeids RANDOMIZE_IDS
                        Randomize the ID Numbers? (1 = Yes, 0 = No)
  -nz NUMBER_PER_REDSHIFT_BIN, --nobjectsperz NUMBER_PER_REDSHIFT_BIN
                        Number of Objects Per Redshift Bin)
  -fi FILTERS_FILE, --filters_file FILTERS_FILE
                        Filters File Name
  -co COMBINE_FILENAME, --coutputfilename COMBINE_FILENAME
                        Filename for combined file?
  -mf, --make_fits      Make Fits File?

```

Here is an example of the input:
		
`% python JAGUAR_to_Subsample.py -in /Path/To/Your/Mock_Catalog_Files/ -rid 0 -sfo sf_output.dat -qo q_output.dat -co all_output.dat -fi filters.dat -mf`

In this example, the path to the mock catalog file is specified, the ID numbers
are not randomized, and then the output SF and Q files are specified (which
triggers the creation of both files). I've also specified that a final combined
file should be created. A filters file is specified, which is one
column including the filters that are desired. Finally, the make fits flag is
set, so two fits files will also be created. 

It is also possible to supply a list of IDs to create a subsample, using the 
`-iIDs` flag:

`% python JAGUAR_to_Subsample.py -in /Path/To/Your/Mock_Catalog_Files/ -rid 0 -iID your_list_of_IDs.dat -fi filters.dat -mf`

This will produce the file `ID_output_list.dat` (and a fits version if the `-mf` flag is
set).

### `JAGUAR_to_Region_Sample.py`
This (optional) script takes in the mock catalogs from JAGUAR, along with a width and a
height, and then produces a list of ID numbers for objects within some random region of
that size within the 11'x11' JAGUAR input file. This is useful if you want to retain
the number density of objects but want to focus on a smaller subregion.
		
```
usage: JAGUAR_to_Region_Sample.py [-h] -in INPUT_FOLDER -width WIDTH -height
                                  HEIGHT

optional arguments:
  -h, --help            show this help message and exit
  -in INPUT_FOLDER, --inputfolder INPUT_FOLDER
                        Input Folder With Mock Catalog Files
  -width WIDTH          Width of subsample (in arcminutes, max 11')
  -height HEIGHT        Height of subsample (in arcminutes, max 11')
```

Here is an example of the input:
		
`% python JAGUAR_to_Region_Sample.py -in /Path/To/Your/Mock_Catalog_Files/ -width 4.0 -height 4.0 > All_4.0_by_4.0_Objects.IDs.dat`

In this example, the path to the mock catalog file is specified, and the code will
then just pick a random region within the full JAGUAR catalog of the specified size,
which will print to screen, so it's helpful to pipe this to a file, here called 
`All_4.0_by_4.0_Objects.IDs.dat`

You can then use these IDs to create the subsample with `JAGUAR_to_Subsample.py`, using the 
`-iIDs` flag:

`% python JAGUAR_to_Subsample.py -in /Path/To/Your/Mock_Catalog_Files/ -rid 0 -iID All_4.0_by_4.0_Objects.IDs.dat -fi filters.dat -mf`



### `Subsample_to_NoisySubsample.py`
This script takes in the output files from `JAGUAR_to_Subsample.py` and then
adds noise to the fluxes, producing noisy fluxes and uncertainties depending
on the HST and/or NIRCam depths you specify. Because these depths are driven
by the CANDELs survey (Grogin et al) and the JADES survey, only those filters
will be produced. Currently, the HST and NIRCam uncertainties are derived
based on calculating the flux and errors through an aperture based on the size and 
sersic index, and stacking individual images.
```
usage: Subsample_to_NoisySubsample.py [-h] -in INPUT_FILE -out OUTPUT_FILE -ni
                                      NIRCAM_DEPTH -fi FILTERS_FILE [-mf]
                                      [-nircamuserdepths NIRCAMUSERDEPTHS]

optional arguments:
  -h, --help            show this help message and exit
  -in INPUT_FILE, --inputfile INPUT_FILE
                        Input file with subample fluxes
  -out OUTPUT_FILE, --outputfile OUTPUT_FILE
                        Output file with noisy fluxes
  -ni NIRCAM_DEPTH, --nircam_depth NIRCAM_DEPTH
                        NIRCam survey: deep or medium
  -fi FILTERS_FILE, --filters_file FILTERS_FILE
                        Filters File Name
  -mf, --make_fits      Make Fits File?
  -nircamuserdepths NIRCAMUSERDEPTHS, --nircamuserdepths NIRCAMUSERDEPTHS
                        NIRCam User Defined Depths

```

`% python Subsample_to_NoisySubsample.py -in sf_output.fits -o sf_noisy.dat -ni deep -fi filters.dat -mf`

In this example, the name of the input file is specified (text or .fits), and 
the NIRCam depths are set. Next, the filters file is specified (any filter
that is not in the CANDELs or JADES surveys are not used), and finally, the make 
fits flag is set, so a fits file will also be created. You can also specify NIRCam depths
in this way:

`% python Subsample_to_NoisySubsample.py -in sf_output.fits -o sf_noisy.dat -ni user -nircamuserdepths '90, 90, 90, 90, 90, 90, 90, 90, 90, 90' -fi filters.dat -mf`

In this way, you set the `-ni` flag to `user`, and then specify the user depths in
each filter separately. In this example, you're stacking 90 images in each filter
with a base exposure time set within the code in the `base_exposure_time_nircam_user` 
array, which is currently set to 1385 seconds.

### `NoisySubsample_to_PhotoZInput.py`
This script takes in the output files from `Subsample_to_NoisySubsample.py` and
prepares it for EAZY, BPZ, BEAGLE, Le Phare, or ZEBRA photometric redshift codes.
		
```
usage: NoisySubsample_to_PhotoZInput.py [-h] -in INPUT_FILE [-out OUTPUT_NAME]
                                        [-beagle] [-eazy] [-bpz] [-lep] [-zeb]

optional arguments:
  -h, --help            show this help message and exit
  -in INPUT_FILE, --inputfile INPUT_FILE
                        Input file with noisy fluxes
  -out OUTPUT_NAME, --output_name OUTPUT_NAME
                        Output Filename Extension?
  -beagle, --beagle_input_file
                        Make BEAGLE file?
  -eazy, --eazy_input_file
                        Make EAZY file?
  -bpz, --bpz_input_file
                        Make BPZ file?
  -lep, --lephare_input_file, --lephare
                        Make Le Phare file?
  -zeb, --zebra_input_file, --zebra
                        Make Zebra file?

```


`% python NoisySubsample_to_PhotoZInput.py -in input_filename.fits -eazy -beagle -bpz -lep -out 5_1_18`

In this example, the name of the input file is specified (text or .fits), and 
then flags for eazy, beagle, bpz, and lephare to create those input files. Finally,
we provide a string to put into the output files (`5_1_18`)

This also requires another file, `NoisySubsample_to_PhotoZInput_filters.dat`, which specifies
various filter files for the individual photometric redshift codes. This should match up
against the filters for the input file.

```
HST_F435W F435W HST-ACS-wfc_F435W_resampled.res HST_ACS_F435W.res
HST_F606W F606W HST-ACS-wfc_F606W_resampled.res HST_ACS_F606W.res
HST_F775W F775W HST-ACS-wfc_F775W_resampled.res HST_ACS_F775W.res
HST_F814W F814W HST-ACS-wfc_F814W_resampled.res HST_ACS_F814W.res
HST_F850LP F850LP HST-ACS-wfc_F850LP_resampled.res HST_ACS_F850LP.res
NRC_F090W F090W F090W_NRConly_ModAB_mean_resampled.res NIRCam_F090W.res
NRC_F115W F115W F115W_NRConly_ModAB_mean_resampled.res NIRCam_F115W.res
NRC_F150W F150W F150W_NRConly_ModAB_mean_resampled.res NIRCam_F150W.res
NRC_F200W F200W F200W_NRConly_ModAB_mean_resampled.res NIRCam_F200W.res
NRC_F277W F277W F277W_NRConly_ModAB_mean_resampled.res NIRCam_F277W.res
NRC_F335M F335M F335M_NRConly_ModAB_mean_resampled.res NIRCam_F335M.res
NRC_F356W F356W F356W_NRConly_ModAB_mean_resampled.res NIRCam_F356W.res
NRC_F410M F410M F410M_NRConly_ModAB_mean_resampled.res NIRCam_F410M.res
NRC_F444W F444W F444W_NRConly_ModAB_mean_resampled.res NIRCam_F444W.res
```

### `ComparePhotoZ_to_SpecZ.py`
This script takes in the output files from EAZY and
calculates statistics and produces plots that show the photometric redshifts compared to the 
spectroscopic redshifts as a function of both SNR, and a second parameter from the JAGUAR 
catalog, if you specify. The user can specify whether they want to make an SNR cut (and
which filter to use), or a `prob_z` cut, or a `q_z` cut, for the analysis. The user
can also specify a `q_z` cut over different redshift ranges, if they want to make a more
liberal cut at higher redshift. Finally, the program
produces a summary file of the important statistics from the EAZY output file, and the
user can set the multiple flag to print out other chi-square minima redshifts (and their
associated chi-square values), and they can set the kde flag to run a kde smoothing 
(sigma = 5) on the chi-square surface to help against printing out a huge amount of 
multiple solutions. 
```
usage: ComparePhotoZ_to_SpecZ.py [-h] -input INPUT_PHOTOMETRY -nrcf
                                 NIRCAM_FILTER -eazy EAZY_OUTPUT_FILE
                                 [-snrl SNR_LIMIT] [-eazyprob EAZYPROB_LIMIT]
                                 [-eazyqz EAZYQZ_LIMIT]
                                 [-eazyqz_zrange EAZYQZ_ZRANGE]
                                 [-outf OPT_OUTPUT_FOLDER] [-minz MINIMUM_Z]
                                 [-maxz MAXIMUM_Z] [-magmin MAG_MIN]
                                 [-magmax MAG_MAX] [-mp] [-outliers]
                                 [-pointdensity] [-multiple] [-kde]
                                 [-jaguar JAGUAR_PATH] [-jparam JAGUAR_PARAM]

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_PHOTOMETRY, --input_photometry INPUT_PHOTOMETRY
                        Input Photometry for analysis
  -nrcf NIRCAM_FILTER, --nircam_filter NIRCAM_FILTER
                        NIRCam filter for SNR analysis?
  -eazy EAZY_OUTPUT_FILE, --eazy_output_file EAZY_OUTPUT_FILE
                        EAZY output file path?
  -snrl SNR_LIMIT, --snr_limit SNR_LIMIT
                        SNR Limit for analysis?
  -eazyprob EAZYPROB_LIMIT, --eazyprob_limit EAZYPROB_LIMIT
                        SNR Limit for analysis?
  -eazyqz EAZYQZ_LIMIT, --eazyqz_limit EAZYQZ_LIMIT
                        EAZY Q_z limit?
  -eazyqz_zrange EAZYQZ_ZRANGE, --eazyqz_zrange EAZYQZ_ZRANGE
                        EAZY Q_z limit (as a function of redshift)?
  -outf OPT_OUTPUT_FOLDER, --opt_output_folder OPT_OUTPUT_FOLDER
                        Optional Output Folder?
  -minz MINIMUM_Z, --minimum_z MINIMUM_Z
                        Minimum Redshift for Analysis
  -maxz MAXIMUM_Z, --maximum_z MAXIMUM_Z
                        Maximum Redshift for Analysis
  -magmin MAG_MIN, --magmin MAG_MIN
                        Minimum Mag for Error Histogram
  -magmax MAG_MAX, --magmax MAG_MAX
                        Maximum Mag for Error Histogram
  -mp, --make_plots     Make Plots?
  -outliers             Calculate Catastrophic Outliers?
  -pointdensity         Make Point Density Plots?
  -multiple             Calculate Multiple Photo-z Solutions?
  -kde                  Use a kde on the chi-square values to calculate
                        multiple solutions?
  -jaguar JAGUAR_PATH, --jaguar_path JAGUAR_PATH
                        Path to JAGUAR Catalogs?
  -jparam JAGUAR_PARAM, --jaguar_param JAGUAR_PARAM
                        JAGUAR Parameter for Coloring Plots?
```

`% python ComparePhotoZ_to_SpecZ.py -input /Path/To/Input/Noisy/Photometry/all_fluxes_5_1_18_noisy.dat -nrcf NRC_F200W -snrl 5 -eazy /Path/To/EAZY/output/photz.zout -jaguar /Path/To/Your/Mock_Catalog_Files/ -jparam sSFR -mp -outliers -outf /Path/To/Optional/Output/Folder/`
`% python ComparePhotoZ_to_SpecZ.py -input /Path/To/Input/Noisy/Photometry/all_fluxes_5_1_18_noisy.dat -nrcf NRC_F200W -eazyqz 3.0 -eazyprob 0.0 -eazy /Path/To/EAZY/output/photz.zout -multiple -kde -jaguar /Path/To/Your/Mock_Catalog_Files/ -jparam sSFR -mp -outliers -outf /Path/To/Optional/Output/Folder/`

In this example, we point to the full set of noisy photometry from `Subsample_to_NoisySubsample.py`, and
then I specify a NIRCam filter for printing out the SNR information. Then I specify that I want
to focus on objects with `q_z < 3`, and `prob_z > 0.0`. I could also put in `-eazyqz_zrange '6.0, 3.0, 15.0, 7.5'`
which would do `q_z < 3` for (zphot = 0 - 6) and `q_z < 7.5` for (zphot = 6 - 15.0). 

Next, we point to the EAZY output file, and set the flags to find multiple chi-square
solutions, using a kde method. I also provide the JAGUAR file, and the parameter
for making spec-z vs. photo-z plots colored by that parameter. Next, I specify
that we'd like to make the plots instead of just producing statistics, and I then specify the
optional output folder for putting all of the files.  

The statistics that are produced are:

```
------------------------------------
ALL OBJECTS (N = 117981)
 bias = -0.017 +/- 1.069
 sigma_68 = 0.124
 NMAD = 0.074
 fraction (> 0.15) = 0.283
------------------------------------
prob_z > 0.0, q_z < 3.0 (N = 55423)
 bias = 0.005 +/- 0.122
 sigma_68 = 0.034
 NMAD = 0.029
 fraction (> 0.15) = 0.043
------------------------------------
```

The first is the bias, which is the average and standard deviation of 
`delta_z = (z_spec - z_phot) / (1 + z_spec)`. Next is `sigma_68`, the value of `delta_z` 
that encompasses 68% of the residuals around 0. NMAD is Normalized Median Absolute Deviation
of the residuals, which is defined as `NMAD(delta_z) = 1.48 * Median(abs(delta_z))`. Then
there's the fraction of outliers with `abs(delta_z) > 0.15`. This is done for all objects, 
and then those with various photo-z code flags.

### `Explore_Outliers.py`
This script allows you to look at the galaxy properties of individual selected outliers.
The script takes in the summary file from `ComparePhotoZ_to_SpecZ` and shows a photometric
redshift vs. spectroscopic redshift plot for those objects with a SNR above a specified
limit, and a probability above a given limit, allowing you to select subsamples via the lasso
tool, after which it pops up a subplot showing a histogram of one specified parameter or
a scatter plot comparing two parameters for the subsample as compared to those objects
from the full set of high SNR, high prob objects in the same redshift range as the
selected objects. You can use this to easily see if outlier galaxies share a specific
galaxy property. 
```
% python Explore_Outliers.py -h
usage: Explore_Outliers.py [-h] -snrl SNR_LIMIT -jaguar JAGUAR_PATH -gpar
                           GALAXY_PARAM [-eazy EAZY_OUTPUT_FILE]
                           [-bpz BPZ_OUTPUT_FILE] [-lep LEPHARE_OUTPUT_FILE]
                           [-sgpar SECOND_GALAXY_PARAM] [-minz MINIMUM_Z]
                           [-maxz MAXIMUM_Z]

optional arguments:
  -h, --help            show this help message and exit
  -snrl SNR_LIMIT, --snr_limit SNR_LIMIT
                        SNR Limit for analysis?
  -jaguar JAGUAR_PATH, --jaguar_path JAGUAR_PATH
                        Path to JAGUAR Catalogs?
  -gpar GALAXY_PARAM, --galaxy_param GALAXY_PARAM
                        Galaxy Parameter For Subplot
  -eazy EAZY_OUTPUT_FILE, --eazy_output_file EAZY_OUTPUT_FILE
                        Analyze EAZY output?
  -bpz BPZ_OUTPUT_FILE, --bpz_output_file BPZ_OUTPUT_FILE
                        Analyze BPZ output?
  -lep LEPHARE_OUTPUT_FILE, --lephare_output_file LEPHARE_OUTPUT_FILE
                        Analyze Le Phare output?
  -sgpar SECOND_GALAXY_PARAM, --second_galaxy_param SECOND_GALAXY_PARAM
                        Galaxy Parameter For Subplot
  -minz MINIMUM_Z, --minimum_z MINIMUM_Z
                        Minimum Redshift for Analysis
  -maxz MAXIMUM_Z, --maximum_z MAXIMUM_Z
                        Maximum Redshift for Analysis
```

`% python Explore_Outliers.py -bpz /Path/to/BPZ_results_summary.dat -snrl 5 -jaguar /Path/To/Your/Mock_Catalog_Files/ -gpar mStar -sgpar max_stellar_age`

In this example, the program will plot the BPZ results for those objects with a SNR (in
the band from the run of `ComparePhotoZ_to_SpecZ.py`) above the `snrl` limit, and
then allow you to select objects, wherein it will plot a subplot showing maximum stellar
age to galaxy stellar mass. Right now, the full list of galaxy properties that can be
compared are: `mStar`, `sSFR`, `tau`, and `max_stellar_age`, but more will be supported
soon.  



### `NoisySubsample_to_ColorPlots.py`
This script looks at color-color plots from the Noisy Subsample as a way of exploring
how to select high-redshift objects by color. The user specifies the input catalog,
the filters, and the redshift limit under review, and the code will then find all
objects with a SNR above 3.0 (default, or the user specified value) in the reddest
two filters, at which point it makes a color-color plot, a color-redshift plot, and
a plot showing two curves. The first is the fraction of objects above a given selection 
limit that are above the specified redshift limit as a function of selection limit, and 
the second is the fraction of objects selected this way out of the total number of
objects above the redshift limit. In this way, the user can see where the contamination
is minimized and completeness is maximized, which is calculated and provided. 
```
usage: NoisySubsample_to_ColorPlots.py [-h] -in INPUT_FILE -yf1 YFILTER1 -yf2
                                       YFILTER2 -xf1 XFILTER1 -xf2 XFILTER2
                                       -zlim ZLIMIT [-idlist ID_NUMBER_LIST]
                                       [-outf OPT_OUTPUT_FOLDER]
                                       [-xlim XLIMIT] [-xplim XPLIMIT]
                                       [-xnlim XNLIMIT] [-yslope YSLOPE]
                                       [-yint YINT_VALUE] [-snr SNR] [-ps]
                                       [-verb]

optional arguments:
  -h, --help            show this help message and exit
  -in INPUT_FILE, --inputfile INPUT_FILE
                        Input file with noisy fluxes
  -yf1 YFILTER1, --yfilter1 YFILTER1
                        Y-Axis Filter 1
  -yf2 YFILTER2, --yfilter2 YFILTER2
                        Y-Axis Filter 2
  -xf1 XFILTER1, --xfilter1 XFILTER1
                        X-Axis Filter 1
  -xf2 XFILTER2, --xfilter2 XFILTER2
                        X-Axis Filter 2
  -zlim ZLIMIT, --zlimit ZLIMIT
                        Redshift Limit?
  -idlist ID_NUMBER_LIST, --id_number_list ID_NUMBER_LIST
                        List of ID Numbers?
  -outf OPT_OUTPUT_FOLDER, --opt_output_folder OPT_OUTPUT_FOLDER
                        Optional Output Folder?
  -xlim XLIMIT, --xlimit XLIMIT
                        Limit on X-Axis Color?
  -xplim XPLIMIT, --xplimit XPLIMIT
                        Positive Limit on X-Axis Color?
  -xnlim XNLIMIT, --xnlimit XNLIMIT
                        Negative Limit on X-Axis Color?
  -yslope YSLOPE, --yslope YSLOPE
                        Slope on the Selection Line?
  -yint YINT_VALUE, --yint_value YINT_VALUE
                        Y-Intercept on the Selection Line?
  -snr SNR, --snr SNR   Filter SNR (default = 3.0)
  -ps, --plot_to_screen
                        Display Plot on Screen (Don't Save)?
  -verb, --verbose      Verbose Mode?

```

`% python ignore NoisySubsample_to_ColorPlots.py -in /Path/To/Noisy/Output/File.fits -yf1 NRC_F090W -yf2 NRC_F115W -xf1 NRC_F115W -xf2 NRC_F150W -zlim 7.0 -yslope 0.4 -snr 5.0 -ps`

In this example, the program will look at a noisy output file and plot, on the
y-axis, `NRC_F090W - NRC_F115W`, and on the x-axis, `NRC_F115W - NRC_F150W`. It will
look at how these colors can be used to select galaxies at `z > 7.0` by exploring
color cuts with a slope of 0.2, for those objects with SNR > 5.0 in `NRC_F115W` and `NRC_F150W`
filters. Finally, it plots to the screen instead of saving a plot to disk. 


### `NoisySubsample_to_ColorCutAnalysis.py`
This script allows you to explore how objects detected at or above a given SNR in two 
red filters (f2 and f3), can select for dropout candidate galaxies above a given redshift
(often the redshift where the Lyman break enters the f1 filter). The spectroscopic redshift
distribution is shown for galaxies selected with these color cuts, and the user can specify
a distribution in a given galaxy parameter (like stellar mass) to see what types of
interloper objects are being selected with a given red color. The user can also specify
a rejection filter and SNR such that the object should not be detected above that significance
for it to count as a dropout candidate.
```
usage: NoisySubsample_to_ColorCutAnalysis.py [-h] -in INPUT_FILE -f1 FILTER1
                                             -f2 FILTER2 -f3 FILTER3 -clim
                                             CLIMIT [-snr SNR]
                                             [-maglim MAGLIM] [-bf BFILTER]
                                             [-bsnr BSNR]
                                             [-outf OPT_OUTPUT_FOLDER] -zlim
                                             ZLIMIT [-ps] [-verb] [-prout]
                                             [-hp HISTPARAM] [-xclim XCLIMIT]

optional arguments:
  -h, --help            show this help message and exit
  -in INPUT_FILE, --inputfile INPUT_FILE
                        Input file with noisy fluxes
  -f1 FILTER1, --filter1 FILTER1
                        Filter 1
  -f2 FILTER2, --filter2 FILTER2
                        Filter 2
  -f3 FILTER3, --filter3 FILTER3
                        Filter 3
  -clim CLIMIT, --climit CLIMIT
                        Color Limit?
  -snr SNR, --snr SNR   Filter SNR (default = 0.0)
  -maglim MAGLIM, --maglim MAGLIM
                        Magnitude Limit (will supersede SNR limits)
  -bf BFILTER, --bfilter BFILTER
                        X-Axis Filter 1
  -bsnr BSNR, --bsnr BSNR
                        Blue Filter SNR (default = 0.0)
  -outf OPT_OUTPUT_FOLDER, --opt_output_folder OPT_OUTPUT_FOLDER
                        Optional Output Folder?
  -zlim ZLIMIT, --zlimit ZLIMIT
                        Redshift Limit?
  -ps, --plot_to_screen
                        Display Plot on Screen (Don't Save)?
  -verb, --verbose      Verbose Mode?
  -prout, --proutput    Print the output objects?
  -hp HISTPARAM, --hparam HISTPARAM
                        Galaxy Parameter for Histogram?
  -xclim XCLIMIT, --xclimit XCLIMIT
                        Upper limit on filter2 - filter3 color?
```

`% python -W ignore NoisySubsample_to_ColorCutAnalysis.py -in /Path/To/Noisy/Output/File.fits -bf HST_F775W -f1 NRC_F090W -f2 NRC_F115W -f3 NRC_F150W -clim 1.9 -snr 5.0 -bsnr 2.0 -zlim 5.6 -hp mStar -xclim 1.0`

In this example, the Noisy catalog is specified with `-in`, and then the script looks
at objects with F115W and F150W fluxes above a SNR of 5.0, and with a F775W SNR less
then 2.0. It then finds all objects with F090W - F115W > 5.0, and F115W - F150W < 1.0, 
and explores how many are at z > 5.6 vs z < 5.6 (the wavelength where the Lyman break 
enters the F090W filter), and plots the spectroscopic redshift histogram of these two samples. 
It also plots a mass histogram showing the mass distribution of these two samples. You could
also print the output objects to the screen if the sample is small enough. 


### `JAGUAR_plot_SEDs.py`
This script, which is found in the `More_In_Prep_Scripts` directory, allows you to plot
SF or Q (or both) SEDs from the JAGUAR catalog.
```
usage: JAGUAR_plot_SEDs.py [-h] -in INPUT_FOLDER -sedin SED_INPUT_FOLDER
                           [-SFID SF_ID] [-QID Q_ID] [-filt] [-noLyA] [-abmag]
                           [-sp]

optional arguments:
  -h, --help            show this help message and exit
  -in INPUT_FOLDER, --inputfolder INPUT_FOLDER
                        Input Folder With Mock Catalog Files
  -sedin SED_INPUT_FOLDER, --sedinputfolder SED_INPUT_FOLDER
                        Input Folder With Mock Catalog SED Files
  -SFID SF_ID, --SFID SF_ID
                        SF Object ID
  -QID Q_ID, --QID Q_ID
                        Q Object ID
  -filt, --plot_filters
                        Plot the Filters on Top of the Plot?
  -noLyA, --noLyA       Use the non-Lyman alpha photometry?
  -abmag, --abmag       Plot in AB Magnitude Units?
  -sp, --save_plot      Save plot to png?
```

It would be perhaps run in this way:

`% python JAGUAR_plot_SEDs.py -in /Path/To/Your/Mock_Catalog_Files/ -sedin /Path/To/Your/Mock_Catalog_SEDs/ -SFID 269865 -QID 309979 -filt -sp`

This would plot the SF object 269865 and the Q object 309979, and save these to a png file. 
You can specify one SF object or one Q object, or both, as in the example. Also, the `-filt` 
flag is set, which plots the HST and NIRCam filters on top of the SED. By default, the code 
plots in units of nJy, but you can specify AB magnitudes with the `-abmag` flag.


### `NIRCam_Mock_Image_Cutouts.py`
This program generates cutouts of objects or regions from the mock images, following the
RGB color image method from Robert Lupton (http://adsabs.harvard.edu/abs/2004PASP..116..133L). 
The user can specify whether they want a cross designating a specific object, or a scalebar,
or whether or not they want to change the RGB color scaling on the images. This requires 
astropy, as it uses the `make_lupton_rgb` function. 

```
usage: NIRCam_Mock_Image_Cutouts.py [-h] -in INPUT_FOLDER -image IMAGE_FOLDER
                                    [-objid OBJID] [-ra CENTER_RA]
                                    [-dec CENTER_DEC] [-rasize RASIZE]
                                    [-decsize DECSIZE] [-rgb RGBLIST] [-cross]
                                    [-info] [-filters]

optional arguments:
  -h, --help            show this help message and exit
  -in INPUT_FOLDER, --inputfolder INPUT_FOLDER
                        Input Folder With Mock Catalog Files (and images)
  -image IMAGE_FOLDER, --imagefolder IMAGE_FOLDER
                        (Optional) Alternate Image Folder
  -objid OBJID, --objid OBJID
                        Object ID
  -ra CENTER_RA, --ra CENTER_RA
                        Center RA of Box
  -dec CENTER_DEC, --dec CENTER_DEC
                        Center DEC of Box
  -rasize RASIZE, --rasize RASIZE
                        Width of Cutout in Arcseconds (Default = 5)
  -decsize DECSIZE, --decsize DECSIZE
                        Height of Cutout in Arcseconds (Default = 5)
  -rgb RGBLIST, --rgblist RGBLIST
                        RGB Color Scaling List, e.g. 1,1,1
  -cross, --cross       Show crosshairs?
  -info, --info         Show info on image?
  -filters, --filters   Show filter names on image?
```

It would be perhaps run in this way:

`% python NIRCam_Mock_Image_Cutouts.py -in /Path/To/Your/Mock_Catalog_Files/ -image /Path/To/Your/Mock_Catalog_Images/ -objid 280879 -rasize 10 -decsize 10 -rgb 0.7,1,1 -cross`

Here, after pointing to where the input SF and Q mock catalog files are, and to where the
images live on your machine, you can specify an Object ID (or an RA and a DEC), and a 
size of the region you want the image to span. Here, I've also set it to put a cross 
around the object, and I'm scaling the r-band image by 0.7. Sometimes,
because the mock images do not span the entire JAGUAR catalog, the object you specify will
not be found in the image, and the program will quit with an error.

The images are hardwired in the code, if you want to change which filters are being used,
you'll want to change lines 19-21. Right now, they're set to:

```
r_name = 'goods_s_F200W_2018_08_29.fits'
g_name = 'goods_s_F115W_2018_08_29.fits'
b_name = 'goods_s_F090W_2018_08_29.fits'
```
