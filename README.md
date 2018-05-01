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
  -mf, --make_fits      Make Fits File?
```

Here is an example of the input:
		
`% python JAGUAR_to_Subsample.py -in /Path/To/Your/Mock_Catalog_Files/ -rid 0 -sfo sf_output.dat -qo q_output.dat -fi filters.dat -mf`

In this example, the path to the mock catalog file is specified, the ID numbers
are not randomized, and then the output SF and Q files are specified (which
triggers the creation of both files). A filters file is specified, which is one
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
we provide a string to put into the output files.  

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
