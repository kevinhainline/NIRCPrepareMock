# NIRCPrepareMock
Code for taking the JAGUAR Mock Catalog and preparing it for photo-z codes

### `JAGUAR_to_Subsample.py`
This script takes in the mock catalogs from JAGUAR, and then spits out a file
that has a subsample of the full objects, such that the subsample tries to 
equally sample redshift and mass space (instead of being weighted towards 
low redshift, low mass objects). It's run with a variety of flags.
		
```
usage: JAGUAR_to_Subsample.py [-h] -in INPUT_FOLDER [-sfo SF_OUTPUT_FILENAME]
                       [-qo Q_OUTPUT_FILENAME] [-rid RANDOMIZE_IDS]
                       [-nz NUMBER_PER_REDSHIFT_BIN] [-fi FILTERS_FILE]
                       [-mf]

optional arguments:
  -h, --help            show this help message and exit
  -in INPUT_FOLDER, --inputfolder INPUT_FOLDER
                        Input Folder With Mock Catalog Files
  -sfo SF_OUTPUT_FILENAME, --sfoutputfilename SF_OUTPUT_FILENAME
                        SF Output Filename
  -qo Q_OUTPUT_FILENAME, --qoutputfilename Q_OUTPUT_FILENAME
                        Q Output Filename
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

### `Subsample_to_NoisySubsample.py`
This script takes in the output files from `JAGUAR_to_Subsample.py` and then
adds noise to the fluxes, producing noisy fluxes and uncertainties depending
on the HST and/or NIRCam depths you specify. Because these depths are driven
by the CANDELs survey (Grogin et al) and the JADES survey, only those filters
will be produced. Currently, the HST noise addition is very simple, while the
NIRCam addition will take into account the size, estimating an aperture
correction. 
		
```
usage: Subsample_to_NoisySubsample.py [-h] -in INPUT_FILE -out OUTPUT_FILE
		                                      -hst HST_DEPTH -ni NIRCAM_DEPTH -fi
		                                      FILTERS_FILE [-mf]

optional arguments:
  -h, --help            show this help message and exit
  -in INPUT_FILE, --inputfile INPUT_FILE
                        Input file with subample fluxes
  -out OUTPUT_FILE, --outputfile OUTPUT_FILE
                        Output file with noisy fluxes
  -hst HST_DEPTH, --hst_depth HST_DEPTH
                        HST survey: deeppointsource, deepextended, or
                        flankingpointsource
  -ni NIRCAM_DEPTH, --nircam_depth NIRCAM_DEPTH
                        NIRCam survey: deep or medium
  -fi FILTERS_FILE, --filters_file FILTERS_FILE
                        Filters File Name
  -mf, --make_fits      Make Fits File?
```


`% python Subsample_to_NoisySubsample.py -in sf_output.fits -o sf_noisy.dat -hst flankingpointsource -ni deep -fi filters.dat -mf`

In this example, the name of the input file is specified (text or .fits), and 
the HST and NIRCam depths are set. Next, the filters file is specified (any filter
that is not in the CANDELs or JADES surveys are not used), and finally, the make 
fits flag is set, so a fits file will also be created. 


### `NoisySubsample_to_PhotoZInput.py`
This script takes in the output files from `Subsample_to_NoisySubsample.py` and
prepares it for EAZY, BPZ, BEAGLE, Le Phare, or ZEBRA photometric redshift codes.
		
```
usage: NoisySubsample_to_PhotoZInput.py [-h] -in INPUT_FILE -fi FILTERS_FILE
                                        [-out OUTPUT_NAME] [-beagle] [-eazy]
                                        [-bpz] [-lep] [-zeb]

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
This script takes in the output files from EAZY, BPZ, and Le Phare (at the moment) and
calculates statistics and produces plots that show the photometric redshifts compared to the 
spectroscopic redshifts as a function of both SNR, and a second parameter from the JAGUAR 
catalog, if you specify.


```
usage: ComparePhotoZ_to_SpecZ.py [-h] -input INPUT_PHOTOMETRY -nrcf
                                 NIRCAM_FILTER -snrl SNR_LIMIT
                                 [-eazy EAZY_OUTPUT_FILE]
                                 [-bpz BPZ_OUTPUT_FILE]
                                 [-lep LEPHARE_OUTPUT_FILE]
                                 [-zeb ZEBRA_OUTPUT_FILE] [-minz MINIMUM_Z]
                                 [-maxz MAXIMUM_Z] [-mp] [-jaguar JAGUAR_PATH]
                                 [-jparam JAGUAR_PARAM]

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_PHOTOMETRY, --input_photometry INPUT_PHOTOMETRY
                        Input Photometry for analysis
  -nrcf NIRCAM_FILTER, --nircam_filter NIRCAM_FILTER
                        NIRCam filter for SNR analysis?
  -snrl SNR_LIMIT, --snr_limit SNR_LIMIT
                        SNR Limit for analysis?
  -outf OPT_OUTPUT_FOLDER, --opt_output_folder OPT_OUTPUT_FOLDER
                        Optional Output Folder?
  -eazy EAZY_OUTPUT_FILE, --eazy_output_file EAZY_OUTPUT_FILE
                        Analyze EAZY output?
  -bpz BPZ_OUTPUT_FILE, --bpz_output_file BPZ_OUTPUT_FILE
                        Analyze BPZ output?
  -lep LEPHARE_OUTPUT_FILE, --lephare_output_file LEPHARE_OUTPUT_FILE
                        Analyze Le Phare output?
  -zeb ZEBRA_OUTPUT_FILE, --zebra_output_file ZEBRA_OUTPUT_FILE
                        Analyze Zebra output?
  -minz MINIMUM_Z, --minimum_z MINIMUM_Z
                        Minimum Redshift for Analysis
  -maxz MAXIMUM_Z, --maximum_z MAXIMUM_Z
                        Maximum Redshift for Analysis
  -mp, --make_plots     Make Plots?
  -outliers             Calculate Catastrophic Outliers?
  -jaguar JAGUAR_PATH, --jaguar_path JAGUAR_PATH
                        Path to JAGUAR Catalogs?
  -jparam JAGUAR_PARAM, --jaguar_param JAGUAR_PARAM
                        JAGUAR Parameter for Coloring Plots?

```

`% python ComparePhotoZ_to_SpecZ.py -input /Path/To/Input/Noisy/Photometry/all_fluxes_5_1_18_noisy.dat -nrcf NRC_F200W -snrl 5 -eazy /Path/To/EAZY/output/photz.zout -jaguar /Path/To/Your/Mock_Catalog_Files/ -jparam sSFR -mp -outliers  -outf /Path/To/Optional/Output/Folder/`

In this example, we point to the full set of noisy photometry from `Subsample_to_NoisySubsample.py`, and
then I specify the SNR filter, and the SNR level, for files that cut down on noisy, low
SNR objects. Then I specify that I want to do an EAZY analysis by pointing to the EAZY
output file. I then point to the JAGUAR mock file and specify a galaxy parameter
for making spec-z vs. photo-z plots colored by that parameter. Next, I specify
that we'd like to make the plots instead of just producing statistics, and I say that I
want to produce outlier files of the objects with `delta_z > 0.15`. I then specify the
optional output folder for putting all of the files.  

The statistics that are produced are:

```
------------------
ALL OBJECTS
 bias = 0.121 +/- 1.643
 sigma_68 = 0.066
 NMAD = 0.038
 fraction (> 0.15) = 0.191
------------------
NRC_F200W_SNR > 5
 bias = 0.035 +/- 0.183
 sigma_68 = 0.034
 NMAD = 0.025
 fraction (> 0.15) = 0.067
------------------
```

The first is the bias, which is the average and standard deviation of 
`delta_z = (z_spec - z_phot) / (1 + z_spec)`. Next is `sigma_68`, the value of `delta_z` 
that encompasses 68% of the residuals around 0. NMAD is Normalized Median Absolute Deviation
of the residuals, which is defined as `NMAD(delta_z) = 1.48 * Median(abs(delta_z))`. Then
there's the fraction of outliers with `abs(delta_z) > 0.15`.