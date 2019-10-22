import os
import sys
import math
import astropy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import argparse
from astropy.io import fits
from astropy.table import Table

noisy_jades_filters = ['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 'HST_F850LP', 'NRC_F070W', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W']
number_jades_filters = len(noisy_jades_filters)

def FluxtoABMag(flux):
	return (-5.0 / 2.0) * np.log10(flux) - 48.60

# Which version of the mock are you using
JAGUAR_version = 'r1_v1.2'

######################
# Required Arguments #
######################

parser = argparse.ArgumentParser()

# Input folder
parser.add_argument(
  '-in','--inputfolder',
  help="Input Folder With Mock Catalog Files (and images)",
  action="store",
  type=str,
  dest="input_folder",
  required=True
)

# Input file
parser.add_argument(
  '-no','--noisyfile',
  help="Input file with noisy fluxes",
  action="store",
  type=str,
  dest="noisy_file",
  required=True
)

# Filter 1
parser.add_argument(
  '-f1','--filter1',
  help="Filter 1",
  action="store",
  type=str,
  dest="filter1",
  required=True
)

# Filter 2
parser.add_argument(
  '-f2','--filter2',
  help="Filter 2",
  action="store",
  type=str,
  dest="filter2",
  required=True
)

# Color Selection
parser.add_argument(
  '-colorcut','--colorcut',
  help="Color Cut?",
  action="store",
  type=float,
  dest="colorcut",
  required=True
)

# Comp Parameter
parser.add_argument(
  '-comp_param','--comp_param',
  help="Comparison Parameter?",
  action="store",
  type=str,
  dest="comp_param",
  required=True
)


######################
# Optional Arguments #
######################

# Redshift Minimum 
parser.add_argument(
  '-zmin','--zminimum',
  help="Redshift Minimum? ",
  action="store",
  type=float,
  dest="zminimum",
  required=False
)

# Redshift Minimum 
parser.add_argument(
  '-zmax','--zmaximum',
  help="Redshift Maximum? ",
  action="store",
  type=float,
  dest="zmaximum",
  required=False
)

# SNR Limit 
parser.add_argument(
  '-snr','--snr',
  help="Filter1 AND 2 SNR (default = 0.0)",
  action="store",
  type=float,
  dest="snrfilter",
  required=False
)

# SNR Limit, Filter One? 
parser.add_argument(
  '-snrf1','--snrf1',
  help="Filter1 SNR (default = 0.0)",
  action="store",
  type=float,
  dest="snrfilter1",
  required=False
)

# SNR Limit, Filter Two? 
parser.add_argument(
  '-snrf2','--snrf2',
  help="Filter2 SNR (default = 0.0)",
  action="store",
  type=float,
  dest="snrfilter2",
  required=False
)

# Plot to Screen? 
parser.add_argument(
  '-ps','--plot_to_screen',
  help="Display Plot on Screen (Don't Save)?",
  action="store_true",
  dest="plot_to_screen",
  required=False
)

# Photometric Catalog? 
parser.add_argument(
  '-pz','--pz',
  help="Photometric Catalog To Use? (Column 1: ID, Column 2: z_phot)",
  action="store",
  type=str,
  dest="photozcat",
  required=False
)

# Set up the parsing from above.
args=parser.parse_args()

# Set the JAGUAR Path
JAGUAR_Path = args.input_folder#'/Users/knh/Desktop/NIRCam/mock_catalog_for_public/'

# Set the filters to explore
filter1 = args.filter1
filter2 = args.filter2

# Set the redshift minimum and maximum
if (args.zminimum):
	redshift_min = args.zminimum
else:
	redshift_min = 0.0
	
if (args.zmaximum):
	redshift_max = args.zmaximum
else:
	redshift_max = 15.0

if (args.snrfilter):
	snr_limit_one = args.snrfilter
	snr_limit_two = args.snrfilter
else:
	snr_limit_one = 0.0
	snr_limit_two = 0.0

# Set the SNR limits on the filters
if (args.snrfilter1):
	snr_limit_one = args.snrfilter1
	
if (args.snrfilter2):
	snr_limit_two = args.snrfilter2

# Input Noisy Photometry
input_photometry = args.noisy_file

# Now, let's open the noisy photometry
fitsinput = fits.open(input_photometry)
ID_values = fitsinput[1].data['ID']
redshifts = fitsinput[1].data['redshift']
number_objects = ID_values.size

filter_flux_one = fitsinput[1].data[filter1]
filter_flux_one_err = fitsinput[1].data[filter1+'_err']
filter_flux_two = fitsinput[1].data[filter2]
filter_flux_two_err = fitsinput[1].data[filter2+'_err']

# Let's do a SNR check. 
filter_flux_one_snr = filter_flux_one / filter_flux_one_err
filter_flux_two_snr = filter_flux_two / filter_flux_two_err

# Get the photometric redshifts if the catalog has been specified.
if (args.photozcat):

	# Set up raw arrays
	ID_values_raw = ID_values
	redshifts_raw = redshifts
	number_objects_raw = number_objects
	filter_flux_one_raw = filter_flux_one
	filter_flux_one_err_raw = filter_flux_one_err
	filter_flux_two_raw = filter_flux_two
	filter_flux_two_err_raw = filter_flux_two_err
	filter_flux_one_snr_raw = filter_flux_one_snr
	filter_flux_two_snr_raw = filter_flux_two_snr
	
	photoz_file = np.loadtxt(args.photozcat)
	photoz_IDs = photoz_file[:,0]
	photometric_redshift = photoz_file[:,1]
	number_photoz = len(photoz_IDs)
	# Match the photo-z catalog to the noisy catalog
	photoz_object_indices = np.zeros(number_photoz, dtype = int)
	for x in range(0, number_photoz):
		#print ID_values[0], photoz_IDs[x]
		#print np.where(ID_values == photoz_IDs[x])
		photoz_object_indices[x] = np.where(ID_values == photoz_IDs[x])[0][0]

	ID_values = ID_values[photoz_object_indices]
	spec_redshifts = redshifts[photoz_object_indices]
	redshifts = photometric_redshift
	filter_flux_one = filter_flux_one[photoz_object_indices]
	filter_flux_one_err = filter_flux_one_err[photoz_object_indices]
	filter_flux_two = filter_flux_two[photoz_object_indices]
	filter_flux_two_err = filter_flux_two_err[photoz_object_indices]
	filter_flux_one_snr = filter_flux_one_snr[photoz_object_indices]
	filter_flux_two_snr = filter_flux_two_snr[photoz_object_indices]
	number_objects = number_photoz

SNR_check_objects = np.where((filter_flux_two_snr > snr_limit_two) & (filter_flux_one_snr > snr_limit_one))[0]

filter_one = FluxtoABMag(filter_flux_one[SNR_check_objects]*1e-23*1e-9)
filter_one[np.where(filter_flux_one[SNR_check_objects] < 0)[0]] = 99.00
filter_two = FluxtoABMag(filter_flux_two[SNR_check_objects]*1e-23*1e-9)
filter_two[np.where(filter_flux_two[SNR_check_objects] < 0)[0]] = 99.00

ID_values_SNR = ID_values[SNR_check_objects]	

correct_redshifts = redshifts[SNR_check_objects]
redshift_range = np.where((correct_redshifts > redshift_min) & (correct_redshifts < redshift_max))[0]

filter_color = filter_one - filter_two


# Let's open up the full SF catalog
full_sf_mock_file = JAGUAR_Path+'JADES_SF_mock_'+JAGUAR_version+'.fits'
sftestfits = fits.open(full_sf_mock_file)

# Get redshifts and galaxy properties
catalog_sf_IDs = sftestfits[1].data['ID']
catalog_sf_redshifts = sftestfits[1].data['redshift']
catalog_sf_logmass = sftestfits[1].data['mStar']
catalog_sf_ssfr = sftestfits[1].data['sSFR']
catalog_sf_sfr_100 = sftestfits[1].data['SFR_100']
catalog_sf_max_stellar_age = sftestfits[1].data['max_stellar_age']
catalog_sf_OIII_EW = sftestfits[1].data['O3_5007_EW']

sfhdr = sftestfits[1].header
parameter_found = 0
for z in range(0, len(sfhdr)):
	if (sfhdr[z] == args.comp_param):
		parameter_found = 1
	
if (parameter_found == 1):
	catalog_sf_comparison_parameter = sftestfits[1].data[args.comp_param]
else:
	print "Parameter "+str(args.comp_param)+" not found in Star-forming Catalog."
	catalog_sf_comparison_parameter = np.zeros(len(catalog_sf_IDs))
	
n_sf_objects = len(catalog_sf_IDs)

sf_clean_filter_flux_one = sftestfits[1].data[filter1+'_fnu']
sf_clean_filter_flux_two = sftestfits[1].data[filter2+'_fnu']
sf_clean_filter_one = FluxtoABMag(sf_clean_filter_flux_one*1e-23*1e-9)
sf_clean_filter_one[np.where(sf_clean_filter_flux_one < 0)[0]] = 99.00
sf_clean_filter_two = FluxtoABMag(sf_clean_filter_flux_two*1e-23*1e-9)
sf_clean_filter_two[np.where(sf_clean_filter_flux_two < 0)[0]] = 99.00


# Let's open up the full Q catalog
full_q_mock_file = JAGUAR_Path+'JADES_Q_mock_'+JAGUAR_version+'.fits'
qtestfits = fits.open(full_q_mock_file)

# Get redshifts and galaxy properties
catalog_q_IDs = qtestfits[1].data['ID']
catalog_q_redshifts = qtestfits[1].data['redshift']
catalog_q_logmass = qtestfits[1].data['mStar']
catalog_q_sfr_10 = np.zeros(len(catalog_q_IDs))-300.0
catalog_q_ssfr = qtestfits[1].data['sSFR']
catalog_q_sfr_100 = catalog_q_ssfr + catalog_q_logmass
catalog_q_max_stellar_age = qtestfits[1].data['max_stellar_age']

qhdr = qtestfits[1].header
parameter_found = 0
for z in range(0, len(qhdr)):
	if (qhdr[z] == args.comp_param):
		parameter_found = 1
	
if (parameter_found == 1):
	catalog_q_comparison_parameter = qtestfits[1].data[args.comp_param]
else:
	print "Parameter "+str(args.comp_param)+" not found in Quiescent Catalog."
	catalog_q_comparison_parameter = np.zeros(len(catalog_q_IDs))


n_q_objects = len(catalog_q_IDs)

q_clean_filter_flux_one = qtestfits[1].data[filter1+'_fnu']
q_clean_filter_flux_two = qtestfits[1].data[filter2+'_fnu']
q_clean_filter_one = FluxtoABMag(q_clean_filter_flux_one*1e-23*1e-9)
q_clean_filter_one[np.where(q_clean_filter_flux_one < 0)[0]] = 99.00
q_clean_filter_two = FluxtoABMag(q_clean_filter_flux_two*1e-23*1e-9)
q_clean_filter_two[np.where(q_clean_filter_flux_two < 0)[0]] = 99.00

# Combine the IDs, redshifts, masses, etc, from the two catalogs.
catalog_all_IDs = np.append(catalog_sf_IDs, catalog_q_IDs)
catalog_all_redshifts = np.append(catalog_sf_redshifts, catalog_q_redshifts)
catalog_all_logmass = np.append(catalog_sf_logmass, catalog_q_logmass)
catalog_all_sfr_100 = np.append(catalog_sf_sfr_100, catalog_q_sfr_100)
catalog_all_comparison_parameter = np.append(catalog_sf_comparison_parameter, catalog_q_comparison_parameter)
catalog_all_ssfr = np.append(catalog_sf_ssfr, catalog_q_ssfr)
catalog_all_max_stellar_age = np.append(catalog_sf_max_stellar_age, catalog_q_max_stellar_age)

all_clean_filter_one = np.append(sf_clean_filter_one, q_clean_filter_one)
all_clean_filter_two = np.append(sf_clean_filter_two, q_clean_filter_two)

n_all_objects = len(catalog_all_IDs)

# What is the clean catalog color (for comparison)
clean_filter_color = all_clean_filter_one - all_clean_filter_two

# Match the noisy catalog to the full catalog
all_object_indices = np.zeros(number_objects, dtype = int)
for x in range(0, number_objects):
	all_object_indices[x] = np.where(catalog_all_IDs == ID_values[x])[0][0]


color_selection = args.colorcut#0.4#0.55

#plt.clf()
plt.figure(figsize=(20, 6))
plt.subplot(131)
plt.scatter(correct_redshifts, filter_color, color = 'black', s = 3, label = 'All Objects')
plt.scatter(correct_redshifts[redshift_range], filter_color[redshift_range], color = 'red', s = 2, label = 'In Redshift Range')
if (color_selection >= 0):
	above_color = np.where(filter_color[redshift_range] >= color_selection)[0]
	plt.scatter(correct_redshifts[redshift_range][above_color], filter_color[redshift_range][above_color], color = 'orange', s = 1, label = 'In Redshift Range, Above '+filter1+' - '+filter2+' > '+str(color_selection))
if (color_selection < 0):
	below_color = np.where(filter_color[redshift_range] < color_selection)[0]
	plt.scatter(correct_redshifts[redshift_range][below_color], filter_color[redshift_range][below_color], color = 'orange', s = 1, label = 'In Redshift Range, Below '+filter1+' - '+filter2+' < '+str(color_selection))
plt.axvline(redshift_min, color = 'red')
plt.axvline(redshift_max, color = 'red')
if (args.photozcat):
	plt.xlabel('Photometric Redshift')
else:
	plt.xlabel('Redshift')
plt.ylabel(filter1+' - '+filter2)
plt.legend()
#plt.show()

# Make the plot! 
comparison_parameter = catalog_all_comparison_parameter
parameter_name = args.comp_param#'[OIII] 5007 EW'#'sSFR'
x_limits = [-1.5, 1.5]
#y_limits = [-9.0, -7.8]
y_limits = [np.min(comparison_parameter[all_object_indices][SNR_check_objects]), np.max(comparison_parameter[all_object_indices][SNR_check_objects])]

#plt.clf()
plt.subplot(132)
plt.scatter(filter_color, comparison_parameter[all_object_indices][SNR_check_objects], s=2, c = 'black', alpha = 0.1, label = 'All Objects')
#plt.scatter(clean_filter_color[all_object_indices][SNR_check_objects], comparison_parameter[all_object_indices][SNR_check_objects], s=2, color = 'blue', alpha = 0.1)
plt.scatter(filter_color[redshift_range], comparison_parameter[all_object_indices][SNR_check_objects][redshift_range], s=8, c = 'red', label = 'In Redshift Range')
plt.axvline(color_selection, color = 'orange')
if (color_selection >= 0):
	plt.scatter(filter_color[redshift_range][above_color], comparison_parameter[all_object_indices][SNR_check_objects][redshift_range][above_color], s=2, c = 'orange', label = 'In Redshift Range, Above '+filter1+' - '+filter2+' > '+str(color_selection))#'red')
if (color_selection < 0):
	plt.scatter(filter_color[redshift_range][below_color], comparison_parameter[all_object_indices][SNR_check_objects][redshift_range][below_color], s=2, c = 'orange', label = 'In Redshift Range, Below '+filter1+' - '+filter2+' < '+str(color_selection))#'red')
plt.xlabel(filter1+' - '+filter2)
plt.ylabel(parameter_name)
plt.axis([x_limits[0],x_limits[1],y_limits[0], y_limits[1]])
#plt.legend()
#plt.show()

#plt.clf()
plt.subplot(133)
parameter_bins = np.arange(y_limits[0],y_limits[1],(y_limits[1] - y_limits[0])/100.0)#,0.1)
above_color = np.where(filter_color[redshift_range] > color_selection)[0]
plt.hist(comparison_parameter[all_object_indices][SNR_check_objects][redshift_range], bins = parameter_bins, color = 'red', label = 'In Redshift Range')
if (color_selection >= 0):
	plt.hist(comparison_parameter[all_object_indices][SNR_check_objects][redshift_range][above_color], bins = parameter_bins, color = 'orange', label = 'In Redshift Range, Above '+filter1+' - '+filter2+' > '+str(color_selection))
if (color_selection < 0):
	plt.hist(comparison_parameter[all_object_indices][SNR_check_objects][redshift_range][below_color], bins = parameter_bins, color = 'orange', label = 'In Redshift Range, Below '+filter1+' - '+filter2+' < '+str(color_selection))
plt.xlabel(parameter_name)
plt.ylabel('N per redshift bin')
#plt.legend()
plt.tight_layout()

if (args.plot_to_screen):
	plt.show()
else:
	# -f1 NRC_F335M -snrf1 10.0 -f2 NRC_F356W -snrf2 10.0 -zmin 3.82 -zmax 4.4 -color -0.4 -comp_param HBaA_6563_EW
	output_filename_root = 'color_vs_param_'+args.filter1+'-'+args.filter2
	if (color_selection >= 0):
		output_filename_root = output_filename_root+'_gt_'+str(color_selection)
	else:
		output_filename_root = output_filename_root+'_lt_'+str(color_selection)
	output_filename_root = output_filename_root+'_vs_'+args.comp_param
	if (args.snrfilter):
		output_filename_root = output_filename_root+'_SNR_'+str(args.snrfilter)
	if (args.snrfilter1):
		output_filename_root = output_filename_root+'_'+args.filter1+'_SNR_'+str(args.snrfilter1)
	if (args.snrfilter2):
		output_filename_root = output_filename_root+'_'+args.filter2+'_SNR_'+str(args.snrfilter2)
	output_filename = output_filename_root+'.png'
	plt.savefig(output_filename, dpi = 200)

print "The IDs of the objects selected with this color cut are:"
if (color_selection >= 0):
	print ID_values_SNR[redshift_range][above_color]
if (color_selection < 0):
	print ID_values_SNR[redshift_range][below_color]