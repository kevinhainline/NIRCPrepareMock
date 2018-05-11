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
noisy_jades_filters = ['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 'HST_F850LP', 'NRC_F070W', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W']
number_jades_filters = len(noisy_jades_filters)

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

# Optional Output Folder
parser.add_argument(
  '-outf','--opt_output_folder',
  help="Optional Output Folder?",
  action="store",
  type=str,
  dest="opt_output_folder",
  required=False
)

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

# Explore Outliers?
parser.add_argument(
  '-outliers',
  help="Calculate Catastrophic Outliers?",
  action="store_true",
  dest="outliers",
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

if (args.opt_output_folder):
	optional_output_folder = args.opt_output_folder
	if not optional_output_folder.endswith("/"):
		optional_output_folder = optional_output_folder + "/"
	if not os.path.exists(optional_output_folder):
	    os.makedirs(optional_output_folder)
else:
	optional_output_folder = ""
	

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

# args.input_photometry
if (args.input_photometry.endswith('.fits')):
	fitsinput = fits.open(args.input_photometry)
	ID_values = fitsinput[1].data['ID']
	redshifts = fitsinput[1].data['redshift']
	number_objects = ID_values.size

	NIRc_flux = fitsinput[1].data[args.nircam_filter]
	NIRc_err = fitsinput[1].data[args.nircam_filter+'_err']

else:
	filter_index_array = -99
	for x in range(0, number_jades_filters):
		if (noisy_jades_filters[x] == args.nircam_filter):
			filter_index_array = x
	if (filter_index_array < 0):
		sys.exit("Specified filter: "+args.nircam_filter+", is not in the list")
	else:
		filter_index = filter_index_array

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
	version = 'EAZY'
	output_folder = optional_output_folder+'EAZY_analysis/'
	if not os.path.exists(output_folder):
	    os.makedirs(output_folder)

	output_zs = np.loadtxt(args.eazy_output_file) 
	ids = output_zs[:,0]
	z_spec = output_zs[:,1]
	z_phot = output_zs[:,13]

	title_for_plot = 'EAZY Results'
	

# Le Phare
if (args.lephare_output_file):
	version = 'LePhare'
	output_folder = optional_output_folder+'LePhare_analysis/'
	if not os.path.exists(output_folder):
	    os.makedirs(output_folder)

	output_zs = np.loadtxt(args.lephare_output_file) 
	ids = output_zs[:,0]
	z_spec = redshifts
	z_phot = output_zs[:,1]

	title_for_plot = 'Le Phare Results'

# BPZ
if (args.bpz_output_file):
	version = 'BPZ'
	output_folder = optional_output_folder+'BPZ_analysis/'
	if not os.path.exists(output_folder):
	    os.makedirs(output_folder)

	output_zs = np.loadtxt(args.bpz_output_file) 
	ids = output_zs[:,0]
	z_spec = redshifts
	z_phot = output_zs[:,1]

	title_for_plot = 'BPZ Results'

IDs_highSNR = np.zeros(n_highSNR_objects)
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
		IDs_highSNR[x_highSNR_objects] = ID_values[x]
		
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

print ""+args.nircam_filter+"_SNR > "+str(args.snr_limit)
print " bias = %.3f +/- %.3f" % (residual_mean_SNR_5, residual_std_SNR_5)
print " sigma_68 = %.3f" % (residual_68_value_SNR_5)
print " NMAD = %.3f" % (NMAD_SNR_5)
print " fraction (> 0.15) = %.3f" % (fraction_gt_15_SNR_5)
#print " skewness = %.3f" % (skewness)
print "------------------"


# # # # # # # # # # # # # #
#  Catastrophic Outliers  #
# # # # # # # # # # # # # #

if (args.outliers):
	catastrophic_outlier_IDs = find_catastrophic_outliers(ID_values, z_spec, z_phot)
	n_catastrophic_outliers = len(catastrophic_outlier_IDs)
	catastrophic_outlier_IDs_highSNR = find_catastrophic_outliers(IDs_highSNR, zspec_highSNR, zphot_highSNR)
	n_catastrophic_outliers_highSNR = len(catastrophic_outlier_IDs_highSNR)

	f = open(output_folder+version+'_catastrophic_outliers.dat', 'a')
	g = open(output_folder+version+'_catastrophic_outliers_SNR_gt_'+str(args.snr_limit)+'.dat', 'a')

	# Now I need to output these IDs to a noisy input file, for creating Photo-Z Input Models. 
	if (args.input_photometry.endswith('.fits')):
				
		for x in range(0, n_catastrophic_outliers):
			ID_index_cat = np.where(ID_values == catastrophic_outlier_IDs[x])[0][0]
			ID_value_cat = ID_values[ID_index_cat]
			redshift_cat = redshifts[ID_index_cat]
			f.write(str(ID_value_cat)+' '+str(redshift_cat)+'  ')
			for y in range(0, number_jades_filters):
				flux_value_cat = fitsinput[1].data[noisy_jades_filters[y]][ID_index_cat]
				flux_value_err_cat = fitsinput[1].data[noisy_jades_filters[y]+'_err'][ID_index_cat]
#			for y in range(0, number_filters):
#				flux_value_cat = fitsinput[1].data[header_filters[y]][ID_index_cat]
#				flux_value_err_cat = fitsinput[1].data[header_filters[y]+'_err'][ID_index_cat]
				f.write(str(flux_value_cat)+' '+str(flux_value_err_cat)+'  ')
			f.write('\n')
		
		for x in range(0, n_catastrophic_outliers_highSNR):
			ID_index_cat_highSNR = np.where(ID_values == catastrophic_outlier_IDs_highSNR[x])[0][0]
			ID_value_cat_highSNR = ID_values[ID_index_cat_highSNR]
			redshift_cat_highSNR = redshifts[ID_index_cat_highSNR]
			g.write(str(ID_value_cat_highSNR)+' '+str(redshift_cat_highSNR)+'  ')
			for y in range(0, number_jades_filters):
				flux_value_cat_highSNR = fitsinput[1].data[noisy_jades_filters[y]][ID_index_cat_highSNR]
				flux_value_err_cat_highSNR = fitsinput[1].data[noisy_jades_filters[y]+'_err'][ID_index_cat_highSNR]
#			for y in range(0, number_filters):
#				flux_value_cat_highSNR = fitsinput[1].data[header_filters[y]][ID_index_cat_highSNR]
#				flux_value_err_cat_highSNR = fitsinput[1].data[header_filters[y]+'_err'][ID_index_cat_highSNR]
				g.write(str(flux_value_cat_highSNR)+' '+str(flux_value_err_cat_highSNR)+'  ')
			g.write('\n')
	
	else:

		for x in range(0, n_catastrophic_outliers):
			ID_index_cat = np.where(ID_values == catastrophic_outlier_IDs[x])[0][0]
			ID_value_cat = ID_values[ID_index_cat]
			redshift_cat = redshifts[ID_index_cat]
			f.write(str(ID_value_cat)+' '+str(redshift_cat)+'  ')
			for y in range(0, number_jades_filters):
				flux_value_cat = full_input_file[ID_index_cat,(y+1)*2]
				flux_value_err_cat = full_input_file[ID_index_cat,((y+1)*2)+1]
#			for y in range(0, number_filters):
#				for k in range(0, number_jades_filters):
#					if (header_filters[y] == noisy_jades_filters[k]):
#						flux_value_cat = full_input_file[ID_index_cat,(k+1)*2]
#						flux_value_err_cat = full_input_file[ID_index_cat,((k+1)*2)+1]
				f.write(str(flux_value_cat)+' '+str(flux_value_err_cat)+'  ')
			f.write('\n')
			
		for x in range(0, n_catastrophic_outliers_highSNR):
			ID_index_cat_highSNR = np.where(ID_values == catastrophic_outlier_IDs_highSNR[x])[0][0]
			ID_value_cat_highSNR = ID_values[ID_index_cat_highSNR]
			redshift_cat_highSNR = redshifts[ID_index_cat_highSNR]
			g.write(str(ID_value_cat_highSNR)+' '+str(redshift_cat_highSNR)+'  ')
			for y in range(0, number_jades_filters):
				flux_value_cat_highSNR = full_input_file[ID_index_cat_highSNR,(y+1)*2]
				flux_value_err_cat_highSNR = full_input_file[ID_index_cat_highSNR,((y+1)*2)+1]
#			for y in range(0, number_filters):
#				for k in range(0, number_jades_filters):
#					if (header_filters[y] == noisy_jades_filters[k]):
#						flux_value_cat_highSNR = full_input_file[ID_index_cat_highSNR,(k+1)*2]
#						flux_value_err_cat_highSNR = full_input_file[ID_index_cat_highSNR,((k+1)*2)+1]
				g.write(str(flux_value_cat_highSNR)+' '+str(flux_value_err_cat_highSNR)+'  ')
			g.write('\n')
			
		#full_input_file = np.loadtxt(args.input_photometry)
		
		#ID_values = full_input_file[:,0]
		#redshifts = full_input_file[:,1]
		#NIRc_flux = full_input_file[:,(filter_index+1)*2]
		#NIRc_err = full_input_file[:,((filter_index+1)*2)+1]
		#number_objects = len(ID_values)

	f.close()
	g.close()

# # # # # # # # #
#  Make  Plots  #
# # # # # # # # #

if (args.make_plots):

	# Residual Histograms
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

	# Plot photo-z vs. spec-z and delta_z vs. spec-z 
	colorlabel = 'log(SNR$_{'+args.nircam_filter+'}$)'
	name = output_folder+'/z_phot_vs_z_spec.png'
	photz_specz_offset(z_spec, z_phot, logNIRc_SNR, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
	name = output_folder+'/z_phot_vs_z_spec_SNR_'+str(SNR_limit)+'.png'
	photz_specz_offset(zspec_highSNR, zphot_highSNR, logNIRc_highSNR, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
	
	# Plot photo-z vs. spec-z and delta_z vs. spec-z, colored by a JAGUAR parameter
	if ((path_to_jaguar_cat_given == 1) and (jaguar_param_given == 1)):
		colorlabel = args.jaguar_param
		name = output_folder+'/z_phot_vs_z_spec_'+args.jaguar_param+'_all.png'
		photz_specz_offset(z_spec, z_phot, catalog_photoz_param, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
		name = output_folder+'/z_phot_vs_z_spec_'+args.jaguar_param+'_SNR_'+str(SNR_limit)+'.png'
		photz_specz_offset(zspec_highSNR, zphot_highSNR, catalog_photoz_param_highSNR, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
	
	# Plot delta_z against zspec as a function of SNR
	name = output_folder+'/z_phot_vs_z_spec_deltazvsz.png'
	plot_deltazvsz(z_spec, z_phot, logNIRc_SNR, name, args.nircam_filter)
	
	# Plot the SNR histograms
	name = output_folder+'/SNR_histograms.png'
	plot_histograms(z_spec, z_phot, logNIRc_SNR, name, args.nircam_filter, title_for_plot)
	
	# Plot the Outlier Fraction as a function of magnitude
	name = output_folder+'/outlier_fraction_mag_'+args.nircam_filter+'.png'
	outlier_fraction_vs_mag(z_spec, z_phot, NIRc_mag, name, args.nircam_filter, title_for_plot)
	
	# Plot the Outlier Fraction as a function of redshift
	name = output_folder+'/outlier_fraction_z.png'
	outlier_fraction_vs_redshift(z_spec, z_phot, name, title_for_plot)
	
	# Plot the Outlier Fraction as a function of redshift for the high SNR objects
	name = output_folder+'/outlier_fraction_z_SNR_'+str(SNR_limit)+'.png'
	outlier_fraction_vs_redshift(zspec_highSNR, zphot_highSNR, name, title_for_plot)

# Unused
	#name = output_folder+'/z_phot_vs_z_spec_SNR_'+str(SNR_limit)+'.png'
	#plotty(zspec_highSNR, zphot_highSNR, logNIRc_highSNR, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
	#name = output_folder+'/z_phot_vs_z_spec_SNR_all.png'
	#plotty(z_spec, z_phot, logNIRc_SNR, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
#		colorlabel = args.jaguar_param
#		name = output_folder+'/z_phot_vs_z_spec_'+args.jaguar_param+'_SNR_'+str(SNR_limit)+'.png'
#		plotty(zspec_highSNR, zphot_highSNR, catalog_photoz_param_highSNR, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
#		name = output_folder+'/z_phot_vs_z_spec_'+args.jaguar_param+'_all.png'
#		plotty(z_spec, z_phot, catalog_photoz_param, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
	#name = output_folder+'/offset_vs_redshift.png'
	#offset_vs_specz(z_spec, z_phot, title_for_plot, name, min_redshift, max_redshift)
	
	#name = output_folder+'/offset_vs_redshift_SNR_'+str(SNR_limit)+'.png'
	#offset_vs_specz(zspec_highSNR, zphot_highSNR, title_for_plot, name, min_redshift, max_redshift)
	#name = output_folder+'/outlier_fraction_mag_'+args.nircam_filter+'_SNR_'+str(SNR_limit)+'.png'
	#outlier_fraction_vs_mag(zspec_highSNR, zphot_highSNR, NIRc_mag_highSNR, name, args.nircam_filter, title_for_plot)
