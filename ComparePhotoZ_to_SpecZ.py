import os
import ast
import sys
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import astropy.stats
import scipy.stats
from scipy import ndimage
from scipy.signal import argrelextrema
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from PlottingFunctions import *

filter_file_name = 'NoisySubsample_to_PhotoZInput_filters.dat'
#filter_file_name = 'NoisySubsample_to_PhotoZInput_medium_filters.dat'
#noisy_jades_filters = ['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 'HST_F850LP', 'NRC_F070W', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W']
noisy_jades_filters = ['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 'HST_F850LP', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W']
number_jades_filters = len(noisy_jades_filters)

# Which version of the mock are you using
JAGUAR_version = 'r1_v1.2'

def FluxtoABMag(flux):
	return (-5.0 / 2.0) * np.log10(flux) - 48.60

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


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

# Analyze EAZY output file
parser.add_argument(
  '-eazy','--eazy_output_file',
  help="EAZY output file path?",
  action="store",
  type=str,
  dest="eazy_output_file",
  required=True
)

######################
# Optional Arguments #
######################

# NIRCam SNR Limit
parser.add_argument(
  '-snrl','--snr_limit',
  help="SNR Limit for analysis?",
  action="store",
  type=float,
  dest="snr_limit",
  required=False
)

# EAZY prob Limit
parser.add_argument(
  '-eazyprob','--eazyprob_limit',
  help="SNR Limit for analysis?",
  action="store",
  type=float,
  dest="eazyprob_limit",
  required=False
)

# EAZY qz Limit
parser.add_argument(
  '-eazyqz','--eazyqz_limit',
  help="EAZY Q_z limit?",
  action="store",
  type=float,
  dest="eazyqz_limit",
  required=False
)

# EAZY qz Limit as a function of redshift
parser.add_argument(
  '-eazyqz_zrange','--eazyqz_zrange',
  help="EAZY Q_z limit (as a function of redshift)?",
  action="store",
  type=str,
  dest="eazyqz_zrange",
  required=False
)

# Optional Output Folder
parser.add_argument(
  '-outf','--opt_output_folder',
  help="Optional Output Folder?",
  action="store",
  type=str,
  dest="opt_output_folder",
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

# Minimum Mag For Error Histogram
parser.add_argument(
  '-magmin','--magmin',
  help="Minimum Mag for Error Histogram",
  action="store",
  type=float,
  dest="mag_min",
  required=False
)

# Maximum Mag For Error Histogram
parser.add_argument(
  '-magmax','--magmax',
  help="Maximum Mag for Error Histogram",
  action="store",
  type=float,
  dest="mag_max",
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

# Calculate Multiple Solutions?
parser.add_argument(
  '-multiple',
  help="Calculate Multiple Photo-z Solutions?",
  action="store_true",
  dest="multiple",
  required=False
)

# Calculate Multiple Solutions?
parser.add_argument(
  '-kde',
  help="Use a kde on the chi-square values to calculate multiple solutions?",
  action="store_true",
  dest="kde",
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

print "Setting up the inputs. "
if (args.opt_output_folder):
	optional_output_folder = args.opt_output_folder
	if not optional_output_folder.endswith("/"):
		optional_output_folder = optional_output_folder + "/"
	if not os.path.exists(optional_output_folder):
	    os.makedirs(optional_output_folder)
else:
	optional_output_folder = ""
	

# Set the minimum and maximum redshifts
if (args.minimum_z is not None):
	min_redshift = args.minimum_z
else:
	min_redshift = 0
	
if (args.maximum_z is not None):
	max_redshift = args.maximum_z
else:
	max_redshift = 15

nircam_SNR_filter = args.nircam_filter
# SNR_limit_flag = 1 # 1 for yes, 0 for no
if (args.snr_limit is not None):
	print "Setting SNR limit"
	SNR_limit_set = True
	SNR_limit = args.snr_limit
else:
	print "No SNR limit"
	SNR_limit_set = False
	
# Let's open up the input file

print "Opening up the filters file. "
# First, we have to open the filters file
# OPEN THE FILTERS FILE
filters_full = np.loadtxt(filter_file_name, dtype='str')
header_filters = filters_full[:,0]
number_filters = len(header_filters)

print "Opening up the input noisy photometry file. "
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
NIRc_SNR_raw = np.zeros(number_objects)
n_highSNR_objects = 0
for x in range(0, number_objects):
	NIRc_SNR_raw[x] = NIRc_SNR[x]
	if (args.snr_limit):
		if (NIRc_SNR[x] >= args.snr_limit):
			n_highSNR_objects = n_highSNR_objects+1

	if (NIRc_SNR[x] <= 0):
		NIRc_SNR[x] = 1e-5

	logNIRc_SNR[x] = math.log10(NIRc_SNR[x])

if (args.mag_min is not None):
	mag_min = np.float(args.mag_min)
if (args.mag_max is not None):
	mag_max = np.float(args.mag_max)


# Let's open up the catalogs, if they're desired

if (args.jaguar_path is not None):
	path_to_jaguar_cat_given = 1
else:
	path_to_jaguar_cat_given = 0
if (args.jaguar_param is not None):
	jaguar_param_given = 1
else:
	jaguar_param_given = 0

if ((path_to_jaguar_cat_given == 0) and (jaguar_param_given == 1)):
	sys.exit("You can't specify a JAGUAR parameter but not provide a path to the JAGUAR catalogs!")
if ((path_to_jaguar_cat_given == 1) and (jaguar_param_given == 0)):
	print "No JAGUAR parameter specified, so no additional plot will be made."

if ((path_to_jaguar_cat_given == 1) and (jaguar_param_given == 1)):
	print "Opening up the full JAGUAR SF and Q mock files for parameter comparison. "
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
	if ((args.jaguar_param == 'beta') or (args.jaguar_param == 'SFR_10') or (args.jaguar_param == 'SFR_100') or (args.jaguar_param == 'tauV_eff') or (args.jaguar_param == 'A1500')):
		# I need to replace this with something better than just zeroing these
		# values.
		catalog_q_param = qtestfits[1].data['ID']*0
	else:
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
version = 'EAZY'
output_folder = optional_output_folder+'EAZY_analysis/'
if not os.path.exists(output_folder):
	os.makedirs(output_folder)
	
print "Opening up the EAZY output file. "
# id z_spec z_a z_m1 chi_a l68 u68  l95 u95  l99 u99  nfilt q_z z_peak peak_prob z_mc
output_zs = np.loadtxt(args.eazy_output_file) 
ids = output_zs[:,0]
z_spec = output_zs[:,1]
z_a = output_zs[:,2]
chi_a = output_zs[:,4]
z_m1 = output_zs[:,3]
z_peak = output_zs[:,13]
# We use z_peak as the "correct" redshift, but we also will report z_a. 
z_phot = z_peak
if (args.eazyprob_limit is not None):
	prob_limit_set = True
	prob_limit = args.eazyprob_limit
else:
	prob_limit_set = False
	prob_limit = 0.0
z_prob = output_zs[:,14]
if (args.eazyqz_limit is not None):
	q_z_limit_set = True
	q_z_limit = args.eazyqz_limit
else:
	q_z_limit_set = False
	q_z_limit = 0.0
if (args.eazyqz_zrange is not None):
	q_z_range_set = True
	eazyqz_zrange = np.array(ast.literal_eval(args.eazyqz_zrange))
	if (len(eazyqz_zrange) % 2 != 0):
		sys.exit("I don't understand your eazyqz_range input? It should be \'zmax, qz, zmax, qz\', etc")
	number_q_z_values = int(len(eazyqz_zrange)/2)
else:
	q_z_range_set = False

q_z = output_zs[:,12]

print "Starting to create the EAZY summary file. "
w = open(output_folder+version+'_results_summary.dat', 'a')
if (args.multiple):
	from Spectrum_Compare_Functions import *
	photo_z_path = args.eazy_output_file.replace("photz.zout","")
	tempfilt, z_grid, obs_sed, templam, temp_sed, pz, chi2fit = generate_sed_arrays(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY=photo_z_path, CACHE_FILE='Same')

	w.write('#ID  z_spec  z_phot  z_prob  q_z  '+args.nircam_filter+'_SNR  z_a  chisq         z_a_other  chisq_other  \n')
else:
	w.write('#ID  z_spec  z_phot  z_prob  q_z  '+args.nircam_filter+'_SNR  z_a  chisq  \n')

if (args.kde):
	kde_version = 1
else:
	kde_version = 0

if (args.multiple):
	print "Finding Multiple Chi-square minima..."
	if (kde_version == 0):
		print "   (finding *all* chi-quare minima...)"
	else:
		print "   (using a kde to find the chi-quare minima...)"
		
for x in range(0, n_objects):
	if (args.multiple):
		idx_value = np.where(ids == ids[x])[0]
		idx = idx_value[0]
		
		# First, do the version where I just look at all objects within range_frac of the best chi-square
		if (kde_version == 0):
			range_frac = 0.10
	
			minima_indices = []
			minima = []
			close_fits = []
			# Find all of the minima
			minima_indices = argrelextrema(chi2fit[:,x], np.less)[0]
			minima = tempfilt['zgrid'][argrelextrema(chi2fit[:,x], np.less)]
			if len(chi2fit[minima_indices,x] > 0):
				best_fit_chisq = np.min(chi2fit[minima_indices,x])
				best_fit_range = best_fit_chisq + range_frac * best_fit_chisq
				close_fits = np.where(chi2fit[minima_indices,x] < best_fit_range)[0]
				
				w.write(str(np.int(ids[x]))+'  '+str(z_spec[x])+'  '+str(z_phot[x])+'  '+str(z_prob[x])+'  '+str(q_z[x])+'  '+str(round(NIRc_SNR_raw[x],4))+'  '+str(z_a[x])+'  '+str(chi_a[x])+'  ')
				if (len(close_fits) == 1):
					w.write('\n')
				else:
					sorted_indices = np.argsort(chi2fit[minima_indices[close_fits][:],x])
					for z in range(1, len(close_fits)):
						w.write('     '+str(round(minima[close_fits][sorted_indices][z],4))+'  '+str(round(chi2fit[minima_indices[close_fits][sorted_indices][z],x],5))+'  ')
						if (z == (len(close_fits)-1)):
							w.write('\n')
			else:
				# these are objects without EAZY fits. 
				w.write(str(np.int(ids[x]))+'  '+str(z_spec[x])+'  '+str(z_phot[x])+'  '+str(z_prob[x])+'  '+str(q_z[x])+'  '+str(round(NIRc_SNR_raw[x],4))+'  '+str(z_a[x])+'  '+str(chi_a[x])+'\n')	
		# Or do the version where you use a smoothed kde to find the chi-square minima
		else:
			smooth_sigma = 5.0
			minima_indices = []
			minima = []
			close_fits = []
			# Find all of the minima
			minima_indices = argrelextrema(chi2fit[:,x], np.less)[0]
			minima = tempfilt['zgrid'][argrelextrema(chi2fit[:,x], np.less)]
			
			if len(chi2fit[minima_indices,x] > 0):
			
				# Make a smooth version
				chi2fitsmooth = ndimage.gaussian_filter1d(chi2fit[:,x], smooth_sigma)
				# Find the smooth minima redshifts
				minima_smooth = tempfilt['zgrid'][argrelextrema(chi2fitsmooth, np.less)]
				# And their corresponding chi-square values
				
				nearest_redshifts = np.zeros(len(minima_smooth))
				nearest_redshift_chisqs = np.zeros(len(minima_smooth))
				for z in range(0, len(minima_smooth)):
					nearest_redshifts[z] = find_nearest(minima, minima_smooth[z])
					nearest_redshift_chisqs[z] = chi2fit[minima_indices,idx][np.where(minima == nearest_redshifts[z])[0]]
					#nearest_value_indices[z] = np.where(minima == minima_smooth[z])
				
				w.write(str(np.int(ids[x]))+'  '+str(z_spec[x])+'  '+str(z_phot[x])+'  '+str(z_prob[x])+'  '+str(q_z[x])+'  '+str(round(NIRc_SNR_raw[x],4))+'  '+str(z_a[x])+'  '+str(chi_a[x])+'  ')
				sorted_indices = np.argsort(nearest_redshift_chisqs)
				for z in range(0, len(nearest_redshifts)):
					if (abs(round(nearest_redshifts[sorted_indices][z],3) - z_a[x]) > 0.001):
						w.write('     '+str(round(nearest_redshifts[sorted_indices][z],4))+'  '+str(round(nearest_redshift_chisqs[sorted_indices][z],5))+'  ')
					if (z == (len(nearest_redshifts)-1)):
						w.write('\n')
			else:
				# these are objects without EAZY fits. 
				w.write(str(np.int(ids[x]))+'  '+str(z_spec[x])+'  '+str(z_phot[x])+'  '+str(z_prob[x])+'  '+str(q_z[x])+'  '+str(round(NIRc_SNR_raw[x],4))+'  '+str(z_a[x])+'  '+str(chi_a[x])+'\n')	
			
	else:
		w.write(str(np.int(ids[x]))+'  '+str(z_spec[x])+'  '+str(z_phot[x])+'  '+str(z_prob[x])+'  '+str(q_z[x])+'  '+str(round(NIRc_SNR_raw[x],4))+'  '+str(z_a[x])+'  '+str(chi_a[x])+' \n')
w.close()


title_for_plot = 'EAZY Results'

# Start by setting up the selected objects. 
IDs_selected = ids
zspec_selected = z_spec
zphot_selected = z_phot
zprob_selected = z_prob
NIRc_mag_selected = NIRc_mag
logNIRc_SNR_selected = logNIRc_SNR
q_z_selected = q_z
if ((path_to_jaguar_cat_given == 1) and (jaguar_param_given == 1)):
	catalog_photoz_param_selected = catalog_photoz_param


# Now, let's look at SNR. 
if (SNR_limit_set == True):
	print "Cutting by SNR."
	SNR_selected_indices = np.where(NIRc_SNR >= SNR_limit)[0]

	# Now transfer these over to the "selected" arrays
	IDs_selected = IDs_selected[SNR_selected_indices]
	zspec_selected = zspec_selected[SNR_selected_indices]
	zphot_selected = zphot_selected[SNR_selected_indices]
	zprob_selected = zprob_selected[SNR_selected_indices]
	NIRc_mag_selected = NIRc_mag_selected[SNR_selected_indices]
	logNIRc_SNR_selected = logNIRc_SNR_selected[SNR_selected_indices]
	q_z_selected = q_z_selected[SNR_selected_indices]
	if ((path_to_jaguar_cat_given == 1) and (jaguar_param_given == 1)):
		catalog_photoz_param_selected = catalog_photoz_param_selected[SNR_selected_indices]
	
#b = np.arange(0,max_redshift+1)

print "Making cuts based on..."
# And now, let's look at the probability and q_z limits.
print 
if ((prob_limit_set == True) & (q_z_limit_set == False) & (q_z_range_set == False)):
	print "   prob_z only."
	selected_indices = np.where(zprob_selected > prob_limit)[0]
if ((prob_limit_set == True) & (q_z_limit_set == True) & (q_z_range_set == False)):
	print "   prob_z and q_z."
	selected_indices = np.where((zprob_selected > prob_limit) & (q_z_selected <= q_z_limit))[0]
if ((prob_limit_set == True) & (q_z_range_set == True)):
	print "   prob_z and a q_z range."
	for x in range (0, number_q_z_values):
		if (x == 0):
			selected_indices = np.where((zprob_selected > prob_limit) & (zphot_selected <= eazyqz_zrange[0]) & (q_z_selected <= eazyqz_zrange[1]))[0]
		else:
			selected_indices = np.append(selected_indices, np.where((zprob_selected > prob_limit) & (zphot_selected > eazyqz_zrange[(2*x)-2]) & (zphot_selected <= eazyqz_zrange[2*x]) & (q_z_selected <= eazyqz_zrange[(2*x)+1]))[0])
if ((prob_limit_set == False) & (q_z_limit_set == True) & (q_z_range_set == False)):
	print "   q_z only."
	selected_indices = np.where(q_z_selected <= q_z_limit)[0]
if ((prob_limit_set == False) & (q_z_range_set == True)):
	print "   a q_z range only."
	for x in range (0, number_q_z_values):
		if (x == 0):
			selected_indices = np.where((zphot_selected <= eazyqz_zrange[0]) & (q_z_selected <= eazyqz_zrange[1]))[0]
		else:
			selected_indices = np.append(selected_indices, np.where((zphot_selected > eazyqz_zrange[(2*x)-2]) & (zphot_selected <= eazyqz_zrange[2*x]) & (q_z_selected <= eazyqz_zrange[(2*x)+1]))[0])



IDs_selected = IDs_selected[selected_indices]
zspec_selected = zspec_selected[selected_indices]
zphot_selected = zphot_selected[selected_indices]
zprob_selected = zprob_selected[selected_indices]
NIRc_mag_selected = NIRc_mag_selected[selected_indices]
logNIRc_SNR_selected = logNIRc_SNR_selected[selected_indices]
q_z_selected = q_z_selected[selected_indices]
if ((path_to_jaguar_cat_given == 1) and (jaguar_param_given == 1)):
	catalog_photoz_param_selected = catalog_photoz_param_selected[selected_indices]


# # # # # # # # # # # # #
#  Calculate Statisics  #
# # # # # # # # # # # # #
print "Calculating statistics. "
print "------------------------------------"
# ALL OBJECTS
residuals = (z_spec - z_phot) / (1 + z_spec)

residual_mean = np.mean(residuals)
residual_std = np.std(residuals)
residual_68_value = residuals_68(residuals)
NMAD = 1.48 * astropy.stats.median_absolute_deviation(residuals)
fraction_gt_15 = len(np.where(abs(residuals) > 0.15)[0])*1.0/len(abs(residuals))*1.0
#skewness = scipy.stats.skew(residuals)

print "ALL OBJECTS (N = "+str(len(residuals))+")"
print " bias = %.3f +/- %.3f" % (residual_mean, residual_std)
print " sigma_68 = %.3f" % (residual_68_value)
print " NMAD = %.3f" % (NMAD)
print " fraction (> 0.15) = %.3f" % (fraction_gt_15)
#print " skewness = %.3f" % (skewness)
print "------------------------------------"

if (SNR_limit_set == True):
	if ((prob_limit_set == True) & (q_z_limit_set == False) & (q_z_range_set == False)):
		print ""+args.nircam_filter+"_SNR > "+str(args.snr_limit)+", prob_z > "+str(args.eazyprob_limit)+" (N = "+str(len(IDs_selected))+")"
	if ((prob_limit_set == True) & (q_z_limit_set == True) & (q_z_range_set == False)):
		print ""+args.nircam_filter+"_SNR > "+str(args.snr_limit)+", prob_z > "+str(args.eazyprob_limit)+", q_z < "+str(args.eazyqz_limit)+" (N = "+str(len(IDs_selected))+")"
	if ((prob_limit_set == True) & (q_z_range_set == True)):
		title_all = ""+args.nircam_filter+"_SNR > "+str(args.snr_limit)+", prob_z > "+str(args.eazyprob_limit)+", "
		for x in range (0, number_q_z_values):
			if (x == 0):
				title_all = title_all+'q_z < '+str(eazyqz_zrange[1])+' (0 < z < '+str(eazyqz_zrange[0])+')'
			else:
				title_all = title_all+', q_z < '+str(eazyqz_zrange[(2*x)+1])+' ('+str(eazyqz_zrange[(2*x)-2])+' < z < '+str(eazyqz_zrange[2*x])+') '				
		title_all = title_all+",  (N = "+str(len(IDs_selected))+")"
		print title_all
	if ((prob_limit_set == False) & (q_z_limit_set == True) & (q_z_range_set == False)):
		print ""+args.nircam_filter+"_SNR > "+str(args.snr_limit)+", q_z < "+str(args.eazyqz_limit)+" (N = "+str(len(IDs_selected))+")"
	if ((prob_limit_set == False) & (q_z_range_set == True)):
		title_all = ""+args.nircam_filter+"_SNR > "+str(args.snr_limit)+", "
		for x in range (0, number_q_z_values):
			if (x == 0):
				title_all = title_all+'q_z < '+str(eazyqz_zrange[1])+' (0 < z < '+str(eazyqz_zrange[0])+')'
			else:
				title_all = title_all+', q_z < '+str(eazyqz_zrange[(2*x)+1])+' ('+str(eazyqz_zrange[(2*x)-2])+' < z < '+str(eazyqz_zrange[2*x])+') '				
		title_all = title_all+",  (N = "+str(len(IDs_selected))+")"
		print title_all
else:
	if ((prob_limit_set == True) & (q_z_limit_set == False) & (q_z_range_set == False)):
		print "prob_z > "+str(args.eazyprob_limit)+" (N = "+str(len(IDs_selected))+")"
	if ((prob_limit_set == True) & (q_z_limit_set == True) & (q_z_range_set == False)):
		print "prob_z > "+str(args.eazyprob_limit)+", q_z < "+str(args.eazyqz_limit)+" (N = "+str(len(IDs_selected))+")"
	if ((prob_limit_set == True) & (q_z_range_set == True)):
		title_all = "prob_z > "+str(args.eazyprob_limit)+", "
		for x in range (0, number_q_z_values):
			if (x == 0):
				title_all = title_all+'q_z < '+str(eazyqz_zrange[1])+' (0 < z < '+str(eazyqz_zrange[0])+')'
			else:
				title_all = title_all+', q_z < '+str(eazyqz_zrange[(2*x)+1])+' ('+str(eazyqz_zrange[(2*x)-2])+' < z < '+str(eazyqz_zrange[2*x])+') '				
		title_all = title_all+",  (N = "+str(len(IDs_selected))+")"
		print title_all
	if ((prob_limit_set == False) & (q_z_limit_set == True) & (q_z_range_set == False)):
		print "q_z < "+str(args.eazyqz_limit)+" (N = "+str(len(IDs_selected))+")"
	if ((prob_limit_set == False) & (q_z_range_set == True)):
		title_all = ""
		for x in range (0, number_q_z_values):
			if (x == 0):
				title_all = title_all+'q_z < '+str(eazyqz_zrange[1])+' (0 < z < '+str(eazyqz_zrange[0])+')'
			else:
				title_all = title_all+', q_z < '+str(eazyqz_zrange[(2*x)+1])+' ('+str(eazyqz_zrange[(2*x)-2])+' < z < '+str(eazyqz_zrange[2*x])+') '				
		title_all = title_all+",  (N = "+str(len(IDs_selected))+")"
		print title_all


residuals_selected = (zspec_selected - zphot_selected) / (1 + zspec_selected)

residual_mean_selected = np.mean(residuals_selected)
residual_std_selected = np.std(residuals_selected)
residual_68_value_selected = residuals_68(residuals_selected)
NMAD_selected = 1.48 * astropy.stats.median_absolute_deviation(residuals_selected)
fraction_gt_15_selected = len(np.where(abs(residuals_selected) > 0.15)[0])*1.0/len(abs(residuals_selected))*1.0
#skewness_selected = scipy.stats.skew(residuals_selected)

print " bias = %.3f +/- %.3f" % (residual_mean_selected, residual_std_selected)
print " sigma_68 = %.3f" % (residual_68_value_selected)
print " NMAD = %.3f" % (NMAD_selected)
print " fraction (> 0.15) = %.3f" % (fraction_gt_15_selected)
#print " skewness = %.3f" % (skewness)
print "------------------------------------"




# # # # # # # # # # # # # #
#  Catastrophic Outliers  #
# # # # # # # # # # # # # #

if (args.outliers):
	outlier_faction_value = 0.15
	catastrophic_outlier_IDs, catastrophic_outlier_z_spec, catastrophic_outlier_z_phot = find_catastrophic_outliers(ID_values, z_spec, z_phot, outlier_faction_value)
	n_catastrophic_outliers = len(catastrophic_outlier_IDs)
	
	catastrophic_folder = 'catastrophic_outliers/'
	if not os.path.exists(output_folder+catastrophic_folder):
	    os.makedirs(output_folder+catastrophic_folder)

	# Files with spec_z's and photo_z's
	q = open(output_folder+catastrophic_folder+version+'_catastrophic_outliers_summary.dat', 'a')

	if (args.eazy_output_file):
		q.write('#ID  z_spec  z_phot  z_prob  q_z '+args.nircam_filter+'_SNR \n')

	# Files with the individual objects
	f = open(output_folder+catastrophic_folder+version+'_catastrophic_outliers.dat', 'a')

	# Now I need to output these IDs to a noisy input file, for creating Photo-Z Input Models. 
	if (args.input_photometry.endswith('.fits')):
		
		# All of the outliers		
		for x in range(0, n_catastrophic_outliers):
			ID_index_cat = np.where(ID_values == catastrophic_outlier_IDs[x])[0][0]
			ID_value_cat = ID_values[ID_index_cat]
			redshift_cat = redshifts[ID_index_cat]
			ID_index_photz_cat = np.where(ids == catastrophic_outlier_IDs[x])[0][0]
			if (args.eazy_output_file):
				q.write(str(np.int(ids[ID_index_photz_cat]))+'  '+str(z_spec[ID_index_photz_cat])+'  '+str(z_phot[ID_index_photz_cat])+'  '+str(z_prob[ID_index_photz_cat])+'  '+str(q_z[ID_index_photz_cat])+'  '+str(round(NIRc_SNR_raw[ID_index_photz_cat],4))+'   \n')
			
			#q.write(str(catastrophic_outlier_IDs[x])+' '+str(catastrophic_outlier_z_spec[x])+'  '+str(catastrophic_outlier_z_phot[x]))
			f.write(str(ID_value_cat)+' '+str(redshift_cat)+'  ')
			for y in range(0, number_jades_filters):
				flux_value_cat = fitsinput[1].data[noisy_jades_filters[y]][ID_index_cat]*1e-23*1e-9
				flux_value_err_cat = fitsinput[1].data[noisy_jades_filters[y]+'_err'][ID_index_cat]*1e-23*1e-9
				f.write(str(flux_value_cat)+' '+str(flux_value_err_cat)+'  ')
			#q.write('\n')
			f.write('\n')
		

	else:

		# All of the outliers		
		for x in range(0, n_catastrophic_outliers):
			ID_index_cat = np.where(ID_values == catastrophic_outlier_IDs[x])[0][0]
			ID_value_cat = ID_values[ID_index_cat]
			redshift_cat = redshifts[ID_index_cat]
			ID_index_photz_cat = np.where(ids == catastrophic_outlier_IDs[x])[0][0]
			if (args.eazy_output_file):
				q.write(str(np.int(ids[ID_index_photz_cat]))+'  '+str(z_spec[ID_index_photz_cat])+'  '+str(z_phot[ID_index_photz_cat])+'  '+str(z_prob[ID_index_photz_cat])+'  '+str(q_z[ID_index_photz_cat])+'  '+str(round(NIRc_SNR_raw[ID_index_photz_cat],4))+'   \n')

			#q.write(str(catastrophic_outlier_IDs[x])+' '+str(catastrophic_outlier_z_spec[x])+'  '+str(catastrophic_outlier_z_phot[x]))
			f.write(str(ID_value_cat)+' '+str(redshift_cat)+'  ')
			for y in range(0, number_jades_filters):
				flux_value_cat = full_input_file[ID_index_cat,(y+1)*2]*1e-23*1e-9
				flux_value_err_cat = full_input_file[ID_index_cat,((y+1)*2)+1]*1e-23*1e-9
				f.write(str(flux_value_cat)+' '+str(flux_value_err_cat)+'  ')
			q.write('\n')
			f.write('\n')
			
		
	q.close()
	f.close()
	
	name = output_folder+catastrophic_folder+'/catastrophic_outlier_zspec_vs_zphot.png'
	catastrophic_outlier_plots(z_spec, z_phot, catastrophic_outlier_z_spec, catastrophic_outlier_z_phot, name, title_for_plot, min_redshift, max_redshift)

# # # # # # # # #
#  Make  Plots  #
# # # # # # # # #

if (args.make_plots):

	# Residual Histograms
	plt.clf()
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
	b = plt.hist(residuals_selected[np.where((residuals_selected > -1.0) & (residuals_selected < 1.0))[0]], bins = n_bins)
	max_value = np.max(b[0])+np.max(b[0])*0.1
	plt.plot([residual_mean_selected, residual_mean_selected],[-1, max_value], label = 'bias = '+str(round(residual_mean_selected,4)))
	plt.plot([NMAD_selected, NMAD_selected],[-1, max_value], label = 'NMAD = '+str(round(NMAD_selected,4)))
	plt.plot([residual_68_value_selected, residual_68_value_selected],[-1, max_value], color = 'grey', alpha = 0.3)
	plt.plot([-1.0*residual_68_value_selected, -1.0*residual_68_value_selected],[-1, max_value], label = 'sigma_68 = '+str(round(residual_68_value_selected,4)), color = 'grey', alpha = 0.3)
	plt.xlabel('(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$)')
	plt.ylabel('N')
	if (SNR_limit_set == True):
		if ((prob_limit_set == True) & (q_z_limit_set == False) & (q_z_range_set == False)):
			plt.title(title_for_plot+", "+args.nircam_filter+"_SNR > "+str(args.snr_limit)+", prob_z > "+str(args.eazyprob_limit))
		if ((prob_limit_set == True) & (q_z_limit_set == True) & (q_z_range_set == False)):
			plt.title(title_for_plot+", "+args.nircam_filter+"_SNR > "+str(args.snr_limit)+", prob_z > "+str(args.eazyprob_limit)+", q_z < "+str(args.eazyqz_limit))
		if ((prob_limit_set == True) & (q_z_range_set == True)):
			title_all = title_for_plot+", "+args.nircam_filter+"_SNR > "+str(args.snr_limit)+", prob_z > "+str(args.eazyprob_limit)+", "
			for x in range (0, number_q_z_values):
				if (x == 0):
					title_all = title_all+'q_z < '+str(eazyqz_zrange[1])+' (0 < z < '+str(eazyqz_zrange[0])+')'
				else:
					title_all = title_all+', q_z < '+str(eazyqz_zrange[(2*x)+1])+' ('+str(eazyqz_zrange[(2*x)-2])+' < z < '+str(eazyqz_zrange[2*x])+') '				
			plt.title(title_all)
		if ((prob_limit_set == False) & (q_z_limit_set == True) & (q_z_range_set == False)):
			plt.title(title_for_plot+", "+args.nircam_filter+"_SNR > "+str(args.snr_limit)+", q_z < "+str(args.eazyqz_limit))
		if ((prob_limit_set == False) & (q_z_range_set == False)):
			title_all = title_for_plot+", "+args.nircam_filter+"_SNR > "+str(args.snr_limit)+", "
			for x in range (0, number_q_z_values):
				if (x == 0):
					title_all = title_all+'q_z < '+str(eazyqz_zrange[1])+' (0 < z < '+str(eazyqz_zrange[0])+')'
				else:
					title_all = title_all+', q_z < '+str(eazyqz_zrange[(2*x)+1])+' ('+str(eazyqz_zrange[(2*x)-2])+' < z < '+str(eazyqz_zrange[2*x])+') '				
			plt.title(title_all, fontsize=5)
	else:
		if ((prob_limit_set == True) & (q_z_limit_set == False) & (q_z_range_set == False)):
			plt.title(title_for_plot+", prob_z > "+str(args.eazyprob_limit))
		if ((prob_limit_set == True) & (q_z_limit_set == True) & (q_z_range_set == False)):
			plt.title(title_for_plot+", prob_z > "+str(args.eazyprob_limit)+", q_z < "+str(args.eazyqz_limit))
		if ((prob_limit_set== True) & (q_z_range_set == True)):
			title_all = title_for_plot+", prob_z > "+str(args.eazyprob_limit)+", "
			for x in range (0, number_q_z_values):
				if (x == 0):
					title_all = title_all+'q_z < '+str(eazyqz_zrange[1])+' (0 < z < '+str(eazyqz_zrange[0])+')'
				else:
					title_all = title_all+', q_z < '+str(eazyqz_zrange[(2*x)+1])+' ('+str(eazyqz_zrange[(2*x)-2])+' < z < '+str(eazyqz_zrange[2*x])+') '				
			plt.title(title_all)
		if ((prob_limit_set == False) & (q_z_limit_set == True) & (q_z_range_set == False)):
			plt.title(title_for_plot+", q_z < "+str(args.eazyqz_limit))
		if ((prob_limit_set == False) & (q_z_range_set == True)):
			title_all = title_for_plot+", "
			for x in range (0, number_q_z_values):
				if (x == 0):
					title_all = title_all+'q_z < '+str(eazyqz_zrange[1])+' (0 < z < '+str(eazyqz_zrange[0])+')'
				else:
					title_all = title_all+', q_z < '+str(eazyqz_zrange[(2*x)+1])+' ('+str(eazyqz_zrange[(2*x)-2])+' < z < '+str(eazyqz_zrange[2*x])+') '				
			plt.title(title_all, fontsize=5)

	plt.axis([-0.5,1.0,0.0,max_value])
	plt.legend()
	plt.savefig(output_folder+'/residual_histogram_selected.png', dpi=300)

	# Plot photo-z vs. spec-z and delta_z vs. spec-z 
	colorlabel = 'log(SNR$_{'+args.nircam_filter+'}$)'
	name = output_folder+'/z_phot_vs_z_spec.png'
	photz_specz_offset(z_spec, z_phot, logNIRc_SNR, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
	name = output_folder+'/z_phot_vs_z_spec_selected.png'
	photz_specz_offset(zspec_selected, zphot_selected, logNIRc_SNR_selected, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
	
	# Plot photo-z vs. spec-z and delta_z vs. spec-z, colored by a JAGUAR parameter
	if ((path_to_jaguar_cat_given == 1) and (jaguar_param_given == 1)):
		colorlabel = args.jaguar_param
		name = output_folder+'/z_phot_vs_z_spec_'+args.jaguar_param+'_all.png'
		photz_specz_offset(z_spec, z_phot, catalog_photoz_param, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
		name = output_folder+'/z_phot_vs_z_spec_'+args.jaguar_param+'_selected.png'
		photz_specz_offset(zspec_selected, zphot_selected, catalog_photoz_param_selected, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
		
	# Plot the Outlier Fraction as a function of magnitude
	name = output_folder+'/outlier_fraction_mag_'+args.nircam_filter+'.png'
	outlier_fraction_vs_mag(z_spec, z_phot, NIRc_mag, name, args.nircam_filter, title_for_plot)

	z_gt_6_obj = np.where(z_spec >= 6.0)[0]
	# Plot the Outlier Fraction as a function of magnitude, z > 6
	name = output_folder+'/outlier_fraction_mag_'+args.nircam_filter+'_z_gt_6.png'
	outlier_fraction_vs_mag(z_spec[z_gt_6_obj], z_phot[z_gt_6_obj], NIRc_mag[z_gt_6_obj], name, args.nircam_filter, title_for_plot+', z_spec > 6')

	z_gt_8_obj = np.where(z_spec >= 8.0)[0]
	# Plot the Outlier Fraction as a function of magnitude, z > 8
	name = output_folder+'/outlier_fraction_mag_'+args.nircam_filter+'_z_gt_8.png'
	outlier_fraction_vs_mag(z_spec[z_gt_8_obj], z_phot[z_gt_8_obj], NIRc_mag[z_gt_8_obj], name, args.nircam_filter, title_for_plot+', z_spec > 8')

	z_gt_10_obj = np.where(z_spec >= 10.0)[0]
	# Plot the Outlier Fraction as a function of magnitude, z > 10
	name = output_folder+'/outlier_fraction_mag_'+args.nircam_filter+'_z_gt_10.png'
	outlier_fraction_vs_mag(z_spec[z_gt_10_obj], z_phot[z_gt_10_obj], NIRc_mag[z_gt_10_obj], name, args.nircam_filter, title_for_plot+', z_spec > 10')

	# Plot the Outlier Fraction as a function of redshift
	name = output_folder+'/outlier_fraction_z.png'
	outlier_fraction_vs_redshift(z_spec, z_phot, name, title_for_plot)
	
	# Plot the Outlier Fraction as a function of redshift for the high SNR objects
	name = output_folder+'/outlier_fraction_z_selected.png'
	outlier_fraction_vs_redshift(zspec_selected, zphot_selected, name, title_for_plot)

	# Plot a histogram of Outlier Fraction for objects with a specific magnitude range
	if (args.mag_min):
		if (args.mag_max):
			mag_objects = np.where((NIRc_mag > mag_min) & (NIRc_mag < mag_max))[0]
			name = output_folder+'/error_distribution_'+str(mag_min)+'-'+str(mag_max)+'.png'
			mag_error_histogram(zspec[mag_objects], zphot[mag_objects], mag_min, mag_max, name)
	
	# Plot Point Density Plots
	colorlabel = 'Number of Points'
	name = output_folder+'/z_phot_vs_z_spec_density.png'
	photz_specz_density(z_spec, z_phot, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
	name = output_folder+'/z_phot_vs_z_spec_selected_density.png'
	photz_specz_density(zspec_selected, zphot_selected, name, args.nircam_filter, title_for_plot, min_redshift, max_redshift, colorlabel)
	
	# Plot q_z analysis figures
	if ((q_z_limit_set == True) & (q_z_range_set == False)):
		name_of_file = output_folder+'/Q_z_analysis_limit_'+str(args.eazyqz_limit)+'.png'
		q_z_analysis_plots(z_spec, z_phot, q_z, args.eazyqz_limit, name_of_file, title_for_plot)	
	if (q_z_range_set == True):
		for x in range (0, number_q_z_values):
			if (x == 0):
				name_of_file = output_folder+'/Q_z_analysis_limit_'+str(eazyqz_zrange[1])+'_z_lt_'+str(eazyqz_zrange[0])+'.png'
				good_objects = np.where(z_phot < eazyqz_zrange[1])[0]
				plot_information = 'z < '+str(eazyqz_zrange[0])+', Q_z < '+str(eazyqz_zrange[1])
				q_z_analysis_plots(z_spec[good_objects], z_phot[good_objects], q_z[good_objects], eazyqz_zrange[1], name_of_file, plot_information)
			else:
				name_of_file = output_folder+'/Q_z_analysis_limit_'+str(eazyqz_zrange[(2*x)+1])+'_z_'+str(eazyqz_zrange[(2*x)-2])+'-'+str(eazyqz_zrange[2*x])+'.png'
				good_objects = np.where((z_phot > eazyqz_zrange[(2*x)-2]) & (z_phot <= eazyqz_zrange[2*x]))[0]
				plot_information = str(eazyqz_zrange[(2*x)-2])+' < z < '+str(eazyqz_zrange[2*x])+', Q_z < '+str(eazyqz_zrange[(2*x)+1])
				q_z_analysis_plots(z_spec[good_objects], z_phot[good_objects], q_z[good_objects], eazyqz_zrange[(2*x)+1], name_of_file, plot_information)
