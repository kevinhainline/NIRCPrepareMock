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


`% python Subsample_to_NoisySubsample.py -in sf_fluxes_4_30_18.fits -o test_noisy.dat -hst flankingpointsource -ni deep -fi filters.dat -mf`

In this example, the name of the input file is specified (text or .fits), and 
the HST and NIRCam depths are set. Next, the filters file is specified (any filter
that is not in the CANDELs or JADES surveys are not used), and finally, the make 
fits flag is set, so a fits file will also be created. 


### `NoisySubsample_to_PhotoZInput.py`
This script takes in the output files from `Subsample_to_NoisySubsample.py` and
prepares it for EAZY, BPZ, BEAGLE, Le Phare, or ZEBRA photometric redshift codes.
		
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


`% python Subsample_to_NoisySubsample.py -in sf_fluxes_4_30_18.fits -o test_noisy.dat -hst flankingpointsource -ni deep -fi filters.dat -mf`

In this example, the name of the input file is specified (text or .fits), and 
the HST and NIRCam depths are set. Next, the filters file is specified (any filter
that is not in the CANDELs or JADES surveys are not used), and finally, the make 
fits flag is set, so a fits file will also be created. 
