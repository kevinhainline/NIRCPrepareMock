import os
import sys
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import astropy.stats
import scipy.stats
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from PlottingFunctions import *

filter_file_name = 'NoisySubsample_to_PhotoZInput_filters.dat'

# Which version of the mock are you using
JAGUAR_version = 'r1_v1.1'

def FluxtoABMag(flux):
	return (-5.0 / 2.0) * np.log10(flux) - 48.60

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

# Make Plots?
parser.add_argument(
  '-mp','--make_plots',
  help="Make Plots?",
  action="store_true",
  dest="make_plots",
  required=False
)

# Original JAGUAR Catalog Path
parser.add_argument(
  '-jaguar','--jaguar_path',
  help="Path to JAGUAR Catalogs?",
  action="store",
  type=str,
  dest="jaguar_path",
  required=False
)

# JAGUAR Parameter For Coloring Plots
parser.add_argument(
  '-jparam','--jaguar_param',
  help="JAGUAR Parameter for Coloring Plots?",
  action="store",
  type=str,
  dest="jaguar_param",
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

NIRc_mag = np.zeros(number_objects)
NIRc_SNR = np.zeros(number_objects)
for x in range(0, number_objects):
	if (NIRc_flux[x] > 0):
		NIRc_mag[x] = FluxtoABMag(NIRc_flux[x]*1e-23*1e-9)
		NIRc_SNR[x] = NIRc_flux[x]/NIRc_err[x]
	else:
		NIRc_mag[x] = 99.99
		NIRc_SNR[x] = 0.00
#NIRc_mag = FluxtoABMag(NIRc_flux*1e-23*1e-9)


# Now, let's calculate the SNR ratio for our filter of choice:
#NIRc_SNR = NIRc_flux/NIRc_err
n_objects = len(NIRc_SNR)

logNIRc_SNR = np.zeros(number_objects)
n_highSNR_objects = 0
for x in range(0, number_objects):

	if (NIRc_SNR[x] >= args.snr_limit):
		n_highSNR_objects = n_highSNR_objects+1

	if (NIRc_SNR[x] <= 0):
		NIRc_SNR[x] = 1e-5

	logNIRc_SNR[x] = math.log10(NIRc_SNR[x])

# Let's open up the catalogs, if they're desired

if (args.jaguar_path):
	path_to_jaguar_cat_given = 1
else:
	path_to_jaguar_cat_given = 0
if (args.jaguar_param):
	jaguar_param_given = 1
else:
	jaguar_param_given = 0

if ((path_to_jaguar_cat_given == 0) and (jaguar_param_given == 1)):
	sys.exit("You can't specify a JAGUAR parameter but not provide a path to the JAGUAR catalogs!")
if ((path_to_jaguar_cat_given == 1) and (jaguar_param_given == 0)):
	print "No JAGUAR parameter specified, so no additional plot will be made."

if ((path_to_jaguar_cat_given == 1) and (jaguar_param_given == 1)):
	# Full star-forming catalogue, fits file
	full_sf_mock_file = args.jaguar_path+'JADES_SF_mock_'+JAGUAR_version+'.fits'
	sftestfits = fits.open(full_sf_mock_file)

	# Get redshifts and galaxy properties
	catalog_sf_IDs = sftestfits[1].data['ID']
	catalog_sf_param = sftestfits[1].data[args.jaguar_param]
	n_sf_objects = len(catalog_sf_IDs)
	
	# Full quiescent catalogue, fits file
	full_q_mock_file = args.jaguar_path+'JADES_Q_mock_'+JAGUAR_version+'.fits'
	qtestfits = fits.open(full_q_mock_file)
	
	# Get redshifts and galaxy properties
	catalog_q_IDs = qtestfits[1].data['ID']
	catalog_q_param = qtestfits[1].data[args.jaguar_param]
	n_q_objects = len(catalog_q_IDs)
	
	catalog_all_IDs = np.append(catalog_sf_IDs, catalog_q_IDs)
	catalog_all_param = np.append(catalog_sf_param, catalog_q_param)
	
	photoz_param_indices = np.zeros(number_objects, dtype = int)
	for x in range(0, number_objects):
		photoz_param_indices[x] = np.where(catalog_all_IDs == ID_values[x])[0][0]

	catalog_photoz_param = catalog_all_param[photoz_param_indices]

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

	title_for_plot = 'EAZY Results'

# Le Phare
if (args.lephare_output_file):
	output_folder = 'LePhare_analysis/'
	if not os.path.exists(output_folder):
	    os.makedirs(output_folder)

	output_zs = np.loadtxt(args.lephare_output_file) 
	ids = output_zs[:,0]
	z_spec = redshifts
	z_phot = output_zs[:,1]

	title_for_plot = 'Le Phare Results'

# BPZ
if (args.bpz_output_file):
	output_folder = 'BPZ_analysis/'
	if not os.path.exists(output_folder):
	    os.makedirs(output_folder)

	output_zs = np.loadtxt(args.bpz_output_file) 
	ids = output_zs[:,0]
	z_spec = redshifts
	z_phot = output_zs[:,1]

	title_for_plot = 'BPZ Results'


zspec_highSNR = np.zeros(n_highSNR_objects)
zphot_highSNR = np.zeros(n_highSNR_objects)
NIRc_mag_highSNR = np.zeros(n_highSNR_objects)
logNIRc_highSNR = np.zeros(n_highSNR_objects)
if ((path_to_jaguar_cat_given == 1) and (jaguar_param_given == 1)):
	catalog_photoz_param_highSNR = np.zeros(n_highSNR_objects)
x_highSNR_objects = 0
for x in range(0, n_objects):
	
	if (NIRc_SNR[x] >= SNR_limit):
		zspec_highSNR[x_highSNR_objects] = z_spec[x]
		zphot_highSNR[x_highSNR_objects] = z_phot[x]
		NIRc_mag_highSNR[x_highSNR_objects] = NIRc_mag[x]

		logNIRc_highSNR[x_highSNR_objects] = logNIRc_SNR[x]
		if ((path_to_jaguar_cat_given == 1) and (jaguar_param_given == 1)):
			catalog_photoz_param_highSNR[x_highSNR_objects] = catalog_photoz_param[x]

		x_highSNR_objects = x_highSNR_objects + 1

b = np.arange(0,max_redshift+1)

# # # # # # # # # # # # #
#  Calculate Statisics  #
# # # # # # # # # # # # #
print "------------------"
# ALL OBJECTS
residuals = (z_spec - z_phot) / (1 + z_spec)

residual_mean = np.mean(residuals)
residual_std = np.std(residuals)
residual_68_value = residuals_68(residuals)
NMAD = 1.48 * astropy.stats.median_absolute_deviation(residuals)
fraction_gt_15 = len(np.where(abs(residuals) > 0.15)[0])*1.0/len(abs(residuals))*1.0
#skewness = scipy.stats.skew(residuals)

print "ALL OBJECTS"
print " bias = %.3f +/- %.3f" % (residual_mean, residual_std)
print " sigma_68 = %.3f" % (residual_68_value)
print " NMAD = %.3f" % (NMAD)
print " fraction (> 0.15) = %.3f" % (fraction_gt_15)
#print " skewness = %.3f" % (skewness)
print "------------------"
# HIGH SNR OBJECTS

residuals_SNR_5 = (zspec_highSNR - zphot_highSNR) / (1 + zspec_highSNR)

residual_mean_SNR_5 = np.mean(residuals_SNR_5)
residual_std_SNR_5 = np.std(residuals_SNR_5)
residual_68_value_SNR_5 = residuals_68(residuals_SNR_5)
NMAD_SNR_5 = 1.48 * astropy.stats.median_absolute_deviation(residuals_SNR_5)
fraction_gt_15_SNR_5 = len(np.where(abs(residuals_SNR_5) > 0.15)[0])*1.0/len(abs(residuals_SNR_5))*1.0
#skewness_SNR_5 = scipy.stats.skew(residuals_SNR_5)

print ""+args.nircam_filter+"_SNR > 5"
print " bias = %.3f +/- %.3f" % (residual_mean_SNR_5, residual_std_SNR_5)
print " sigma_68 = %.3f" % (residual_68_value_SNR_5)
print " NMAD = %.3f" % (NMAD_SNR_5)
print " fraction (> 0.15) = %.3f" % (fraction_gt_15_SNR_5)
#print " skewness = %.3f" % (skewness)
print "------------------"



# # # # # # # # #
#  Make  Plots  #
# # # # # # # # #

if (args.make_plots):

	n_bins = 200
	b = plt.hist(residuals[np.where((residuals > -1.0) & (residuals < 1.0))[0]], bins = n_bins)
	max_value = np.max(b[0])+np.max(b[0])*0.1
	plt.plot([residual_mean, residual_mean],[-1, max_value], label = 'bias = '+str(round(residual_mean,4)))
	plt.plot([NMAD, NMAD],[-1, max_value], label = 'NMAD = '+str(round(NMAD,4)))
	plt.plot([residual_68_value, residual_68_value],[-1, max_value], color = 'grey', alpha = 0.3)
	plt.plot([-1.0*residual_68_value, -1.0*residual_68_value],[-1, max_value], label = 'sigma_68 = '+str(round(residual_68_value,4)), color = 'grey', alpha = 0.3)
	plt.xlabel('(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$)')
	plt.ylabel('N')
	plt.title(title_for_plot)
	plt.axis([-0.5,1.0,0.0,max_value])
	plt.legend()
	plt.savefig(output_folder+'/residual_histogram.png', dpi=300)

	plt.clf()
	n_bins = 200
	b = plt.hist(residuals_SNR_5[np.where((residuals_SNR_5 > -1.0) & (residuals_SNR_5 < 1.0))[0]], bins = n_bins)
	max_value = np.max(b[0])+np.max(b[0])*0.1
	plt.plot([residual_mean_SNR_5, residual_mean_SNR_5],[-1, max_value], label = 'bias = '+str(round(residual_mean_SNR_5,4)))
	plt.plot([NMAD_SNR_5, NMAD_SNR_5],[-1, max_value], label = 'NMAD = '+str(round(NMAD_SNR_5,4)))
	plt.plot([residual_68_value_SNR_5, residual_68_value_SNR_5],[-1, max_value], color = 'grey', alpha = 0.3)
	plt.plot([-1.0*residual_68_value_SNR_5, -1.0*residual_68_value_SNR_5],[-1, max_value], label = 'sigma_68 = '+str(round(residual_68_value_SNR_5,4)), color = 'grey', alpha = 0.3)
	plt.xlabel('(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$)')
	plt.ylabel('N')
	plt.title(title_for_plot+', SNR > '+str(SNR_limit))
	plt.axis([-0.5,1.0,0.0,max_value])
	plt.legend()
	plt.savefig(output_folder+'/residual_histogram_SNR_'+str(SNR_limit)+'.png', dpi=300)


	colorlabel = 'log(SNR$_{'+args.nircam_filter+'}$)'
	name = output_folder+'/z_phot_vs_z_spec.png'
	photz_specz_offset(z_spec, z_phot, logNIRc_SNR, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
	name = output_folder+'/z_phot_vs_z_spec_SNR_'+str(SNR_limit)+'.png'
	photz_specz_offset(zspec_highSNR, zphot_highSNR, logNIRc_highSNR, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
	#name = output_folder+'/z_phot_vs_z_spec_SNR_'+str(SNR_limit)+'.png'
	#plotty(zspec_highSNR, zphot_highSNR, logNIRc_highSNR, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
	#name = output_folder+'/z_phot_vs_z_spec_SNR_all.png'
	#plotty(z_spec, z_phot, logNIRc_SNR, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
	
	if ((path_to_jaguar_cat_given == 1) and (jaguar_param_given == 1)):
		colorlabel = args.jaguar_param
		name = output_folder+'/z_phot_vs_z_spec_'+args.jaguar_param+'_SNR_'+str(SNR_limit)+'.png'
		plotty(zspec_highSNR, zphot_highSNR, catalog_photoz_param_highSNR, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
		name = output_folder+'/z_phot_vs_z_spec_'+args.jaguar_param+'_all.png'
		plotty(z_spec, z_phot, catalog_photoz_param, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
	
	name = output_folder+'/z_phot_vs_z_spec_deltazvsz.png'
	plot_deltazvsz(z_spec, z_phot, logNIRc_SNR, name, args.nircam_filter)
	
	name = output_folder+'/SNR_histograms.png'
	plot_histograms(z_spec, z_phot, logNIRc_SNR, name, args.nircam_filter, title_for_plot)
	
	#name = output_folder+'/offset_vs_redshift.png'
	#offset_vs_specz(z_spec, z_phot, title_for_plot, name, min_redshift, max_redshift)
	
	#name = output_folder+'/offset_vs_redshift_SNR_'+str(SNR_limit)+'.png'
	#offset_vs_specz(zspec_highSNR, zphot_highSNR, title_for_plot, name, min_redshift, max_redshift)
	
	name = output_folder+'/outlier_fraction_mag_'+args.nircam_filter+'.png'
	outlier_fraction_vs_mag(z_spec, z_phot, NIRc_mag, name, args.nircam_filter, title_for_plot)
	
	#name = output_folder+'/outlier_fraction_mag_'+args.nircam_filter+'_SNR_'+str(SNR_limit)+'.png'
	#outlier_fraction_vs_mag(zspec_highSNR, zphot_highSNR, NIRc_mag_highSNR, name, args.nircam_filter, title_for_plot)
	
	name = output_folder+'/outlier_fraction_z.png'
	outlier_fraction_vs_redshift(z_spec, z_phot, name, title_for_plot)
	
	name = output_folder+'/outlier_fraction_z_SNR_'+str(SNR_limit)+'.png'
	outlier_fraction_vs_redshift(zspec_highSNR, zphot_highSNR, name, title_for_plot)

