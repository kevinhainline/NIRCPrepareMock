import os
import sys
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from PlottingFunctions import *

filter_file_name = 'NoisySubsample_to_PhotoZInput_filters.dat'

# Parser for figuring out the various input parameters
parser = argparse.ArgumentParser()

######################
# Required Arguments #
######################

# Input Photometry
parser.add_argument(
  '-input','--input_photometry',
  help="Input Photometry for analysis",
  action="store",
  type=str,
  dest="input_photometry",
  required=True
)

# NIRCam SNR Filter
parser.add_argument(
  '-nrcf','--nircam_filter',
  help="NIRCam filter for SNR analysis?",
  action="store",
  type=str,
  dest="nircam_filter",
  required=True
)

# NIRCam SNR Filter
parser.add_argument(
  '-snrl','--snr_limit',
  help="SNR Limit for analysis?",
  action="store",
  type=float,
  dest="snr_limit",
  required=True
)

######################
# Optional Arguments #
######################

# Analyze EAZY output file
parser.add_argument(
  '-eazy','--eazy_output_file',
  help="Analyze EAZY output?",
  action="store",
  type=str,
  dest="eazy_output_file",
  required=False
)

# Generate BPZ input file
parser.add_argument(
  '-bpz','--bpz_output_file',
  help="Analyze BPZ output?",
  action="store",
  type=str,
  dest="bpz_output_file",
  required=False
)

# Generate Le Phare input file
parser.add_argument(
  '-lep','--lephare_output_file',
  help="Analyze Le Phare output?",
  action="store",
  type=str,
  dest="lephare_output_file",
  required=False
)

# Generate Zebra input file
parser.add_argument(
  '-zeb','--zebra_output_file',
  help="Analyze Zebra output?",
  action="store",
  type=str,
  dest="zebra_output_file",
  required=False
)

# Minimum Redshift
parser.add_argument(
  '-minz','--minimum_z',
  help="Minimum Redshift for Analysis",
  action="store",
  type=float,
  dest="minimum_z",
  required=False
)

# Maximum Redshift
parser.add_argument(
  '-maxz','--maximum_z',
  help="Maximum Redshift for Analysis",
  action="store",
  type=float,
  dest="maximum_z",
  required=False
)


args=parser.parse_args()

# Set the minimum and maximum redshifts
if (args.minimum_z):
	min_redshift = args.minimum_z
else:
	min_redshift = 0
	
if (args.maximum_z):
	max_redshift = args.maximum_z
else:
	max_redshift = 15
	
nircam_SNR_filter = args.nircam_filter
# SNR_limit_flag = 1 # 1 for yes, 0 for no
SNR_limit = args.snr_limit

# Let's open up the input file

# First, we have to open the filters file
# OPEN THE FILTERS FILE
filters_full = np.loadtxt(filter_file_name, dtype='str')
header_filters = filters_full[:,0]
number_filters = len(header_filters)

filter_index_array = np.where(header_filters == args.nircam_filter)[0]
if not filter_index_array:
  sys.exit("Specified filter: "+args.nircam_filter+", is not in the list")
else:
	filter_index = filter_index_array[0]

# args.input_photometry
if (args.input_photometry.endswith('.fits')):
	fitsinput = fits.open(args.input_file)
	ID_values = fitsinput[1].data['ID']
	redshifts = fitsinput[1].data['redshift']
	number_objects = ID_values.size

	NIRc_flux = fitsinput[1].data[args.nircam_filter]
	NIRc_err = fitsinput[1].data[args.nircam_filter+'_err']

else:
	full_input_file = np.loadtxt(args.input_photometry)
	
	ID_values = full_input_file[:,0]
	redshifts = full_input_file[:,1]
	NIRc_flux = full_input_file[:,(filter_index+1)*2]
	NIRc_err = full_input_file[:,((filter_index+1)*2)+1]
	number_objects = len(ID_values)


# Now, let's calculate the SNR ratio for our filter of choice:
NIRc_SNR = NIRc_flux/NIRc_err
n_objects = len(NIRc_SNR)

logNIRc_SNR = np.zeros(number_objects)
n_highSNR_objects = 0
for x in range(0, number_objects):

	if (NIRc_SNR[x] >= args.snr_limit):
		n_highSNR_objects = n_highSNR_objects+1

	if (NIRc_SNR[x] <= 0):
		NIRc_SNR[x] = 1e-5

	logNIRc_SNR[x] = math.log10(NIRc_SNR[x])

# Now, let's open up the output files!

# EAZY
if (args.eazy_output_file):
	output_folder = 'EAZY_analysis/'
	if not os.path.exists(output_folder):
	    os.makedirs(output_folder)

	output_zs = np.loadtxt(args.eazy_output_file) 
	ids = output_zs[:,0]
	z_spec = output_zs[:,1]
	z_phot = output_zs[:,13]

	#nircam_depth = depth_of_observation+'.'+hst_depth_val #'rieke.28.5'
	#file_nircam_depth = depth_of_observation #'rieke'
	title_for_plot = 'EAZY Results'


zspec_highSNR = np.zeros(n_highSNR_objects)
zphot_highSNR = np.zeros(n_highSNR_objects)
logNIRc_highSNR = np.zeros(n_highSNR_objects)
x_highSNR_objects = 0
for x in range(0, n_objects):
	
	if (NIRc_SNR[x] >= SNR_limit):
		zspec_highSNR[x_highSNR_objects] = z_spec[x]
		zphot_highSNR[x_highSNR_objects] = z_phot[x]

		logNIRc_highSNR[x_highSNR_objects] = logNIRc_SNR[x]
		x_highSNR_objects = x_highSNR_objects + 1

b = np.arange(0,max_redshift+1)

# And now let's do some plotting. 	
name = output_folder+'/z_phot_vs_z_spec_SNR_'+str(SNR_limit)+'.png'
plotty(zspec_highSNR, zphot_highSNR, logNIRc_highSNR, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift)

name = output_folder+'/z_phot_vs_z_spec_SNR_all.png'
plotty(z_spec, z_phot, logNIRc_SNR, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift)

#xmin = -2
#xmax = 2
#bins = 10

name = output_folder+'/z_phot_vs_z_spec_deltazvsz.png'
plot_deltazvsz(z_spec, z_phot, logNIRc_SNR, name, args.nircam_filter)

name = output_folder+'/SNR_histograms.png'
plot_histograms(z_spec, z_phot, logNIRc_SNR, name, args.nircam_filter, title_for_plot)

