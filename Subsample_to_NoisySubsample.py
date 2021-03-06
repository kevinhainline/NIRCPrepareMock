import os
import ast
import sys
import math
import argparse
import numpy as np
import scipy.special as sc
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table

# Calculate the sersic flux given an effective radius and aperture size. 
def sersic(sersic_nval, r_effval, r_aperval):

	if (sersic_nval <= 0.36):
	    b = (0.01945 - 0.8902*sersic_nval + 10.95*sersic_nval*sersic_nval 
	    	- 19.67*sersic_nval*sersic_nval*sersic_nval 
	    	+ 13.43*sersic_nval*sersic_nval*sersic_nval*sersic_nval)
	else:
		b = (2.0 * sersic_nval - 0.3333333333 + 4.0/(405.0 * sersic_nval) 
			+ 46.0/(25515.0*sersic_nval*sersic_nval) + 131.0/(1148175*sersic_nval*sersic_nval*sersic_nval) 
			- 2194697.0/(3069071775 * sersic_nval * sersic_nval * sersic_nval * sersic_nval))
		
	del_r = r_aperval / 100.0;	
	pb = np.power(b, 2.0*sersic_nval)
	temp = 2.0 * 3.14159 * sc.gamma(2.0*sersic_nval) * r_effval * r_effval
	I_tot = pb / temp
		
	r = np.zeros(100)
	sigma = np.zeros(100)
	ap_flux = np.zeros(100)
	for i in range(0, 100):
		r[i] = (1.0)*(i+1) * del_r
		eta = r[i] / r_effval
		sigma[i] = pb * np.exp(-1.0*b*np.power(eta, 1.0/sersic_nval)) / (temp * sersic_nval)

		if (i == 0):
			ap_flux[i] = sigma[i] * 3.14159 * r[i] * r[i]
		else:
			ap_flux[i] = ap_flux[i-1] + 3.14159 *(r[i]*r[i] - r[i-1]*r[i-1]) * ( sigma[i] + sigma[i-1])/2.0
		
		flux = ap_flux[i]
		
	return flux

# Calculate the sersic flux given an effective radius and aperture size. 
def updated_sersic(sersic_nval, r_effval, r_aperval):

	if (sersic_nval <= 0.36):
	    b = (0.01945 - 0.8902*sersic_nval + 10.95*sersic_nval*sersic_nval 
	    	- 19.67*sersic_nval*sersic_nval*sersic_nval 
	    	+ 13.43*sersic_nval*sersic_nval*sersic_nval*sersic_nval)
	else:
		b = (2.0 * sersic_nval - 0.3333333333 + 4.0/(405.0 * sersic_nval) 
			+ 46.0/(25515.0*sersic_nval*sersic_nval) + 131.0/(1148175*sersic_nval*sersic_nval*sersic_nval) 
			- 2194697.0/(3069071775 * sersic_nval * sersic_nval * sersic_nval * sersic_nval))
	
	a = 2.0*sersic_nval
	ratio = r_aperval / r_effval
	x = b * np.power(ratio, 1.0/sersic_nval)
	
	# print r_effval, sersic_nval, b, a, ratio, x
	ap_corr = sc.gammainc(a, x)

	return ap_corr


# NIRCam photometric values
nfilters = 15 # 
n_aper = 4   # number_of_apertures__set_to_0_for_adaptive
filt_name = ['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 'HST_F850LP', 'NRC_F070W', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W']
pixel_size = np.array([0.049, 0.049, 0.049, 0.049, 0.049, 0.032, 0.032, 0.032, 0.032, 0.032, 0.064, 0.064, 0.064, 0.064, 0.064])
bkgnd_nJy = np.array([4.778, 9.237, 12.240, 12.350, 12.678, 6.39, 6.39, 5.92, 4.97, 3.79, 18.76, 13.82, 14.81, 21.72, 31.60])     # /* backgrounds_for_GOODS_S_in_nJy/pix */
e_sec_nJy = np.array([0.005119048, 0.010952381, 0.005, 0.006506024, 0.002435897, 0.0189, 0.0268, 0.0297, 0.0392, 0.0452, 0.0375, 0.0164, 0.0433, 0.0177, 0.0379])   # /* e/sec/nJy  */
read_noise = np.array([5.6, 5.6, 5.6, 5.6, 5.6, 6.0, 6.0, 6.0, 6.0, 6.0, 9.0, 9.0, 9.0, 9.0, 9.0]) # /* read_noise */

aperture = np.array([0.16, 0.24, 0.32, 0.64])   # photom_apertures_in_arcsec  */

# HST DEPTH
# HST with XDF depths, Table 3, Illingworth et al (2013)
n_exposures_xdf = np.array([164, 286, 460, 362, 700])
base_exposure_time_xdf = np.array([929., 609., 821., 140., 602.])

# HST with  CANDELs depths, Table 6, Grogin et al (2011)
n_exposures_candels = np.array([18, 27, 22, 80, 72])
base_exposure_time_candels = np.array([400, 400, 400, 400, 400])

# HST/ACS with CANDELs depths, GOODS-S, Flanking, Table 6, Grogin et al (2011)
n_exposures_candels_flanking = np.array([18, 22, 18, 19, 75])
base_exposure_time_candels_flanking = np.array([400, 400, 400, 400, 400])

# NIRCam DEPTH
# Medium (Average) Survey Depth
n_exposures_medium_nircam = np.array([ 5, 8, 8, 7, 5, 7, 5, 5, 8, 8])
# Deep (Average) Survey Depth
n_exposures_deep_nircam = np.array([ 0, 44, 58, 43, 28, 35, 22, 28, 44, 44])
# User Survey Depth
user_value = 65
n_exposures_deep_user = np.array([ user_value, user_value, user_value, user_value, user_value, user_value, user_value, user_value, user_value, user_value])
base_exposure_time_nircam_user = np.array([1385., 1385., 1385., 1385., 1385., 1385., 1385., 1385., 1385., 1385.])

# If using DEEP8 with N_groups = 7
base_exposure_time_nircam = np.array([1385., 1385., 1385., 1385., 1385., 1385., 1385., 1385., 1385., 1385.])

# If using DEEP8 with N_groups = 9
#base_exposure_time_nircam = np.array([1804., 1804., 1804., 1804., 1804., 1804., 1804., 1804., 1804., 1804.])

# DEEP
n_exposures_deep = np.append(n_exposures_xdf, n_exposures_deep_nircam)
#n_exposures_deep = np.append(n_exposures_candels, n_exposures_deep_nircam)
base_exposure_time_deep = np.append(base_exposure_time_xdf, base_exposure_time_nircam)
#base_exposure_time_deep = np.append(base_exposure_time_candels, base_exposure_time_nircam)

# MEDIUM CANDELS
#n_exposures_medium = np.append(n_exposures_candels, n_exposures_medium_nircam)
#base_exposure_time_medium = np.append(base_exposure_time_candels, base_exposure_time_nircam)

# MEDIUM CANDELS flanking
n_exposures_medium = np.append(n_exposures_candels_flanking, n_exposures_medium_nircam)
base_exposure_time_medium = np.append(base_exposure_time_candels_flanking, base_exposure_time_nircam)

# USER 
n_exposures_user = np.append(n_exposures_xdf, n_exposures_deep_user)
base_exposure_time_user = np.append(base_exposure_time_xdf, base_exposure_time_nircam_user)


######################
# Required Arguments #
######################

parser = argparse.ArgumentParser()

# Input file
parser.add_argument(
  '-in','--inputfile',
  help="Input file with subample fluxes",
  action="store",
  type=str,
  dest="input_file",
  required=True
)

# Output file
parser.add_argument(
  '-out','--outputfile',
  help="Output file with noisy fluxes",
  action="store",
  type=str,
  dest="output_file",
  required=True
)

# NIRCam Deep or Medium
parser.add_argument(
  '-ni','--nircam_depth',
  help="NIRCam survey: deep or medium",
  action="store",
  type=str,
  dest="nircam_depth",
  required=True
)

# Filters File
parser.add_argument(
  '-fi','--filters_file',
  help="Filters File Name",
  action="store",
  type=str,
  dest="filters_file",
  required=True
)

######################
# Optional Arguments #
######################

# Generate Fits File
parser.add_argument(
  '-mf','--make_fits',
  help="Make Fits File?",
  action="store_true",
  dest="make_fits",
  required=False
)

# User Depths 
parser.add_argument(
  '-nircamuserdepths','--nircamuserdepths',
  help="NIRCam User Defined Depths",
  action="store",
  type=str,
  dest="nircamuserdepths",
  required=False
)

args=parser.parse_args()

if ((args.output_file.endswith('.txt')) or (args.output_file.endswith('.dat'))):
	filenameroot = ('.').join(args.output_file.split('.')[:-1])
else:
	filenameroot = args.output_file


# Read in filters file
filter_file_name = args.filters_file
filters_full = np.loadtxt(filter_file_name, dtype='str')
filters = filters_full
number_filters = len(filters)


# Parse NIRCam depth
if (args.nircam_depth == 'deep'):
	#time = deep_time
	time = base_exposure_time_deep
	n_exposures = n_exposures_deep
elif (args.nircam_depth == 'medium'):
	#time = medium_time
	time = base_exposure_time_medium
	n_exposures = n_exposures_medium
elif (args.nircam_depth == 'user'):
	#time = medium_time
	if (args.nircamuserdepths):
		n_exposures_deep_user = np.array(ast.literal_eval(args.nircamuserdepths))
		n_exposures_user = np.append(n_exposures_xdf, n_exposures_deep_user)	
	time = base_exposure_time_user
	n_exposures = n_exposures_user
else:
	sys.exit("NIRCam depth should be deep or medium or user")

output_file_name = args.output_file

if args.input_file.endswith('.fits'):
	#sys.exit("this is a fits file!")
	fitsinput = fits.open(args.input_file)
	ID_numbers = fitsinput[1].data['ID']
	redshifts = fitsinput[1].data['redshift']
	re_major = fitsinput[1].data['Re_maj']
	sersic_n = fitsinput[1].data['sersic_n']
	axis_ratio = fitsinput[1].data['axis_ratio']
	
	n_objects = ID_numbers.size

	apparent_flux = np.zeros([number_filters, n_objects])
	for j in range(0, number_filters):
		apparent_flux[:][j] = fitsinput[1].data[filters[j]]
	
else:
	# Open up the subsample catalogue file
	catalogue_file = args.input_file
	cat_file_full = np.loadtxt(catalogue_file)
	ID_numbers = cat_file_full[:,0]
	redshifts = cat_file_full[:,1]
	re_major = cat_file_full[:,3]
	sersic_n = cat_file_full[:,4]
	axis_ratio = cat_file_full[:,5]
	n_objects = ID_numbers.size

	apparent_flux = np.zeros([number_filters, n_objects])
	for j in range(0, number_filters):
		apparent_flux[:][j] = cat_file_full[:,6+j]
	
f = open(filenameroot+'.NOFIRSTLINE.dat', 'a')

average_sig_value = np.zeros(number_filters)

# Calculate the circularized half-light radius 
#re_circ = re_major * np.sqrt(axis_ratio)
re_circ = re_major

# Set up aperture correction and aperture size arrays. 
ap_corr = np.zeros(n_objects)
ap_radius = np.zeros(n_objects)


# Go through all of the galaxies. 
for x in range(0, n_objects):

	flux_value = 0
	flux_value_noisy = 0
	
	f.write(str(ID_numbers[x])+' '+str(redshifts[x])+'   ')
	# Determine, based on the re_maj value, what the aperture that we should use to 
	# contain the flux in the object.	
#	if (re_major[x] <= aperture[0]):
#		i_ap = 0
	
#	for j in range(1, n_aper):
#		if ((aperture[j-1] < re_major[x]) & (re_major[x] <= aperture[j])):
#			i_ap = j
	
#	if ((aperture[n_aper-1] < re_major[x])):
#		i_ap = n_aper - 1
			
#	ap_radius = aperture[i_ap]

#	ap_corr = updated_sersic(sersic_n[x], re_major[x], ap_radius)

	# Determine, based on the re_circ value, what the aperture that we should use to 
	# contain the flux in the object.	
	if (re_circ[x] <= aperture[0]):
		i_ap = 0
	
	for j in range(1, n_aper):
		if ((aperture[j-1] < re_circ[x]) & (re_circ[x] <= aperture[j])):
			i_ap = j
	
	if ((aperture[n_aper-1] < re_circ[x])):
		i_ap = n_aper - 1
			
	ap_radius[x] = aperture[i_ap]

	ap_corr[x] = updated_sersic(sersic_n[x], re_circ[x], ap_radius[x])

	f.write(str(ap_radius[x])+' '+str(ap_corr[x])+'   ')
	
	if (ap_corr[x] > 1.0):
		ap_corr[x] = 1.0
	
	final_number_filters = 0
	final_filters = np.empty([0], dtype='S20')
	# Go through the various filters
	for j in range(0, number_filters):
		filter_found = 0
		
		# HST FILTERS AND FLUX ERROR ESTIMATES
		if (filters[j] == 'HST_F435W'):
			filter_index = 0
			if (n_exposures[filter_index] > 0):
				final_number_filters = final_number_filters + 1
				final_filters = np.append(final_filters, ['HST_F435W'])
				filter_found = 1
				flux_value = apparent_flux[j][x]
				if (flux_value < 0):
					flux_value = 0
				pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
				filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
				flux = ap_corr[x] * flux_value + filt_bkgnd
				ftemp = flux * e_sec_nJy[filter_index]
				ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
				ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
				
				sum_flux = 0.0
				sum_flux2 = 0.0
				for k in range(0, n_exposures[filter_index]):
					ftemp = ferror * np.random.normal()
					flux  = ap_corr[x] * flux_value + ftemp
					sum_flux = sum_flux + flux
					sum_flux2 = sum_flux2 + (flux * flux)
	
				filter_flux_value = sum_flux / n_exposures[filter_index]
				filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
				filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
				average_sig_value[j] = average_sig_value[j] + filter_sig_value
			else:
				filter_flux_value = -9999
				filter_sig_value = -9999
				average_sig_value[j] = -9999

		if (filters[j] == 'HST_F606W'):
			filter_index = 1
			if (n_exposures[filter_index] > 0):
				final_number_filters = final_number_filters + 1
				final_filters = np.append(final_filters, ['HST_F606W'])
				filter_found = 1
				flux_value = apparent_flux[j][x]
				if (flux_value < 0):
					flux_value = 0
				pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
				filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
				flux = ap_corr[x] * flux_value + filt_bkgnd
				ftemp = flux * e_sec_nJy[filter_index]
				ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
				ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
				
				sum_flux = 0.0
				sum_flux2 = 0.0
				for k in range(0, n_exposures[filter_index]):
					ftemp = ferror * np.random.normal()
					flux  = ap_corr[x] * flux_value + ftemp
					sum_flux = sum_flux + flux
					sum_flux2 = sum_flux2 + (flux * flux)
		
				filter_flux_value = sum_flux / n_exposures[filter_index]
				filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
				filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
				average_sig_value[j] = average_sig_value[j] + filter_sig_value
			else:
				filter_flux_value = -9999
				filter_sig_value = -9999
				average_sig_value[j] = -9999

		if (filters[j] == 'HST_F775W'):
			filter_index = 2
			if (n_exposures[filter_index] > 0):
				final_number_filters = final_number_filters + 1
				final_filters = np.append(final_filters, ['HST_F775W'])
				filter_found = 1
				flux_value = apparent_flux[j][x]
				if (flux_value < 0):
					flux_value = 0
				pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
				filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
				flux = ap_corr[x] * flux_value + filt_bkgnd
				ftemp = flux * e_sec_nJy[filter_index]
				ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
				ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
					
				sum_flux = 0.0
				sum_flux2 = 0.0
				for k in range(0, n_exposures[filter_index]):
					ftemp = ferror * np.random.normal()
					flux  = ap_corr[x] * flux_value + ftemp
					sum_flux = sum_flux + flux
					sum_flux2 = sum_flux2 + (flux * flux)
		
				filter_flux_value = sum_flux / n_exposures[filter_index]
				filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
				filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
				average_sig_value[j] = average_sig_value[j] + filter_sig_value
			else:
				filter_flux_value = -9999
				filter_sig_value = -9999
				average_sig_value[j] = -9999

		if (filters[j] == 'HST_F814W'):
			filter_index = 3
			if (n_exposures[filter_index] > 0):
				final_number_filters = final_number_filters + 1
				final_filters = np.append(final_filters, ['HST_F814W'])
				filter_found = 1
				flux_value = apparent_flux[j][x]
				if (flux_value < 0):
					flux_value = 0
				pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
				filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
				flux = ap_corr[x] * flux_value + filt_bkgnd
				ftemp = flux * e_sec_nJy[filter_index]
				ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
				ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
					
				sum_flux = 0.0
				sum_flux2 = 0.0
				for k in range(0, n_exposures[filter_index]):
					ftemp = ferror * np.random.normal()
					flux  = ap_corr[x] * flux_value + ftemp
					sum_flux = sum_flux + flux
					sum_flux2 = sum_flux2 + (flux * flux)
		
				filter_flux_value = sum_flux / n_exposures[filter_index]
				filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
				filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
				average_sig_value[j] = average_sig_value[j] + filter_sig_value
			else:
				filter_flux_value = -9999
				filter_sig_value = -9999
				average_sig_value[j] = -9999

		if (filters[j] == 'HST_F850LP'):
			filter_index = 4
			if (n_exposures[filter_index] > 0):
				final_number_filters = final_number_filters + 1
				final_filters = np.append(final_filters, ['HST_F850LP'])
				filter_found = 1
				flux_value = apparent_flux[j][x]
				if (flux_value < 0):
					flux_value = 0
				pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
				filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
				flux = ap_corr[x] * flux_value + filt_bkgnd
				ftemp = flux * e_sec_nJy[filter_index]
				ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
				ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
					
				sum_flux = 0.0
				sum_flux2 = 0.0
				for k in range(0, n_exposures[filter_index]):
					ftemp = ferror * np.random.normal()
					flux  = ap_corr[x] * flux_value + ftemp
					sum_flux = sum_flux + flux
					sum_flux2 = sum_flux2 + (flux * flux)
		
				filter_flux_value = sum_flux / n_exposures[filter_index]
				filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
				filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
				average_sig_value[j] = average_sig_value[j] + filter_sig_value
			else:
				filter_flux_value = -9999
				filter_sig_value = -9999
				average_sig_value[j] = -9999

		# NIRCam FILTERS AND FLUX ERROR ESTIMATES
		if (filters[j] == 'NRC_F070W'):
			if ((args.nircam_depth == 'medium') or (args.nircam_depth == 'user')):
				filter_index = 5
				final_number_filters = final_number_filters + 1
				final_filters = np.append(final_filters, ['NRC_F070W'])
				filter_found = 1
				flux_value = apparent_flux[j][x]
				if (flux_value < 0):
					flux_value = 0
				pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
				filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
				flux = ap_corr[x] * flux_value + filt_bkgnd
				ftemp = flux * e_sec_nJy[filter_index]
				ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
				ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
					
				sum_flux = 0.0
				sum_flux2 = 0.0
				for k in range(0, n_exposures[filter_index]):
					ftemp = ferror * np.random.normal()
					flux  = ap_corr[x] * flux_value + ftemp
					sum_flux = sum_flux + flux
					sum_flux2 = sum_flux2 + (flux * flux)
		
				filter_flux_value = sum_flux / n_exposures[filter_index]
				filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
				filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
				average_sig_value[j] = average_sig_value[j] + filter_sig_value
			else:
				filter_flux_value = -9999
				filter_sig_value = -9999
				average_sig_value[j] = -9999

		if (filters[j] == 'NRC_F090W'):
			filter_index = 6
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F090W'])
			filter_found = 1
			flux_value = apparent_flux[j][x]
			if (flux_value < 0):
				flux_value = 0
			pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
			filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
			flux = ap_corr[x] * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[filter_index]
			ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
			ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
			
			sum_flux = 0.0
			sum_flux2 = 0.0
			for k in range(0, n_exposures[filter_index]):
				ftemp = ferror * np.random.normal()
				flux  = ap_corr[x] * flux_value + ftemp
				sum_flux = sum_flux + flux
				sum_flux2 = sum_flux2 + (flux * flux)
	
			filter_flux_value = sum_flux / n_exposures[filter_index]
			filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
			filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
			average_sig_value[j] = average_sig_value[j] + filter_sig_value


		if (filters[j] == 'NRC_F115W'):
			filter_index = 7
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F115W'])
			filter_found = 1
			flux_value = apparent_flux[j][x]
			if (flux_value < 0):
				flux_value = 0
			pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
			filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
			flux = ap_corr[x] * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[filter_index]
			ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
			ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
			
			sum_flux = 0.0
			sum_flux2 = 0.0
			for k in range(0, n_exposures[filter_index]):
				ftemp = ferror * np.random.normal()
				flux  = ap_corr[x] * flux_value + ftemp
				sum_flux = sum_flux + flux
				sum_flux2 = sum_flux2 + (flux * flux)
	
			filter_flux_value = sum_flux / n_exposures[filter_index]
			filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
			filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
			average_sig_value[j] = average_sig_value[j] + filter_sig_value

		if (filters[j] == 'NRC_F150W'):
			filter_index = 8
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F150W'])
			filter_found = 1
			flux_value = apparent_flux[j][x]
			if (flux_value < 0):
				flux_value = 0
			pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
			filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
			flux = ap_corr[x] * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[filter_index]
			ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
			ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
			
			sum_flux = 0.0
			sum_flux2 = 0.0
			for k in range(0, n_exposures[filter_index]):
				ftemp = ferror * np.random.normal()
				flux  = ap_corr[x] * flux_value + ftemp
				sum_flux = sum_flux + flux
				sum_flux2 = sum_flux2 + (flux * flux)
	
			filter_flux_value = sum_flux / n_exposures[filter_index]
			filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
			filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
			average_sig_value[j] = average_sig_value[j] + filter_sig_value

		if (filters[j] == 'NRC_F200W'):
			filter_index = 9
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F200W'])
			filter_found = 1
			flux_value = apparent_flux[j][x]
			if (flux_value < 0):
				flux_value = 0
			pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
			filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
			flux = ap_corr[x] * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[filter_index]
			ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
			ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
				
			sum_flux = 0.0
			sum_flux2 = 0.0
			for k in range(0, n_exposures[filter_index]):
				ftemp = ferror * np.random.normal()
				flux  = ap_corr[x] * flux_value + ftemp
				sum_flux = sum_flux + flux
				sum_flux2 = sum_flux2 + (flux * flux)
	
			filter_flux_value = sum_flux / n_exposures[filter_index]
			filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
			filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
			average_sig_value[j] = average_sig_value[j] + filter_sig_value

		if (filters[j] == 'NRC_F277W'):
			filter_index = 10
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F277W'])
			filter_found = 1
			flux_value = apparent_flux[j][x]
			if (flux_value < 0):
				flux_value = 0
			pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
			filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
			flux = ap_corr[x] * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[filter_index]
			ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
			ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
				
			sum_flux = 0.0
			sum_flux2 = 0.0
			for k in range(0, n_exposures[filter_index]):
				ftemp = ferror * np.random.normal()
				flux  = ap_corr[x] * flux_value + ftemp
				sum_flux = sum_flux + flux
				sum_flux2 = sum_flux2 + (flux * flux)
	
			filter_flux_value = sum_flux / n_exposures[filter_index]
			filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
			filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
			average_sig_value[j] = average_sig_value[j] + filter_sig_value

		if (filters[j] == 'NRC_F335M'):
			filter_index = 11
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F335M'])
			filter_found = 1
			flux_value = apparent_flux[j][x]
			if (flux_value < 0):
				flux_value = 0
			pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
			filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
			flux = ap_corr[x] * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[filter_index]
			ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
			ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
				
			sum_flux = 0.0
			sum_flux2 = 0.0
			for k in range(0, n_exposures[filter_index]):
				ftemp = ferror * np.random.normal()
				flux  = ap_corr[x] * flux_value + ftemp
				sum_flux = sum_flux + flux
				sum_flux2 = sum_flux2 + (flux * flux)
	
			filter_flux_value = sum_flux / n_exposures[filter_index]
			filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
			filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
			average_sig_value[j] = average_sig_value[j] + filter_sig_value

		if (filters[j] == 'NRC_F356W'):
			filter_index = 12
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F356W'])
			filter_found = 1
			flux_value = apparent_flux[j][x]
			if (flux_value < 0):
				flux_value = 0
			pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
			filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
			flux = ap_corr[x] * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[filter_index]
			ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
			ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
			
			sum_flux = 0.0
			sum_flux2 = 0.0
			for k in range(0, n_exposures[filter_index]):
				ftemp = ferror * np.random.normal()
				flux  = ap_corr[x] * flux_value + ftemp
				sum_flux = sum_flux + flux
				sum_flux2 = sum_flux2 + (flux * flux)

			filter_flux_value = sum_flux / n_exposures[filter_index]
			filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
			filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
			average_sig_value[j] = average_sig_value[j] + filter_sig_value

		if (filters[j] == 'NRC_F410M'):
			filter_index = 13
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F410M'])
			filter_found = 1
			flux_value = apparent_flux[j][x]
			if (flux_value < 0):
				flux_value = 0
			pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
			filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
			flux = ap_corr[x] * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[filter_index]
			ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
			ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
				
			sum_flux = 0.0
			sum_flux2 = 0.0
			for k in range(0, n_exposures[filter_index]):
				ftemp = ferror * np.random.normal()
				flux  = ap_corr[x] * flux_value + ftemp
				sum_flux = sum_flux + flux
				sum_flux2 = sum_flux2 + (flux * flux)
	
			filter_flux_value = sum_flux / n_exposures[filter_index]
			filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
			filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
			average_sig_value[j] = average_sig_value[j] + filter_sig_value

		if (filters[j] == 'NRC_F444W'):
			filter_index = 14
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F444W'])
			filter_found = 1
			flux_value = apparent_flux[j][x]
			if (flux_value < 0):
				flux_value = 0
			pixel_area = 3.14159265 * ap_radius[x] * ap_radius[x] / (pixel_size[filter_index] * pixel_size[filter_index])
			filt_bkgnd = pixel_area * bkgnd_nJy[filter_index]
			flux = ap_corr[x] * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[filter_index]
			ferror = np.sqrt(ftemp * time[filter_index] + (read_noise[filter_index]*pixel_area)*(read_noise[filter_index]))							
			ferror = ferror / (e_sec_nJy[filter_index] * time[filter_index])
				
			sum_flux = 0.0
			sum_flux2 = 0.0
			for k in range(0, n_exposures[filter_index]):
				ftemp = ferror * np.random.normal()
				flux  = ap_corr[x] * flux_value + ftemp
				sum_flux = sum_flux + flux
				sum_flux2 = sum_flux2 + (flux * flux)
	
			filter_flux_value = sum_flux / n_exposures[filter_index]
			filter_sig_value = (sum_flux2 - (sum_flux * sum_flux)/n_exposures[filter_index])/(n_exposures[filter_index]-1)
			filter_sig_value = np.sqrt(filter_sig_value / n_exposures[filter_index])
			average_sig_value[j] = average_sig_value[j] + filter_sig_value

		if(filter_found == 1):
			#print filters[j],str(flux_value),str(flux_value_noisy),str(flux_error)
			f.write( '%.3f %.3f  ' % (filter_flux_value, filter_sig_value))

	f.write(' \n')
			
f.close()

for j in range(0, number_filters):
	average_sig_value[j] =  average_sig_value[j]/(n_objects)

f = open(filenameroot+'.FIRSTLINE.dat', 'a')
f.write('# ID   redshift   aperture_radius   aperture_corr   ')
for j in range(0, final_number_filters):
	f.write(final_filters[j]+' '+final_filters[j]+'_err   ')
f.write(' \n')
f.close()

os.system('cat '+filenameroot+'.FIRSTLINE.dat '+filenameroot+'.NOFIRSTLINE.dat > '+filenameroot+'.dat')
os.system('rm '+filenameroot+'.NOFIRSTLINE.dat')
os.system('rm '+filenameroot+'.FIRSTLINE.dat')


if (args.make_fits):
	# First, let's make the  dtype and colnames arrays
	colnames = np.zeros(4+(final_number_filters*2), dtype ='S20')
	dtype = np.zeros(4+(final_number_filters*2), dtype ='str')
	colnames[0] = 'ID'
	colnames[1] = 'redshift'
	colnames[2] = 'ap_radius'
	colnames[3] = 'ap_corr'

	dtype[0] = 'I'
	dtype[1] = 'f'
	dtype[2] = 'f'
	dtype[3] = 'f'
	for j in range(0, final_number_filters):
		colnames[(j+2)*2] = final_filters[j]
		colnames[((j+2)*2)+1] = final_filters[j]+'_err'
		dtype[(j+2)*2] = 'f8'
		dtype[((j+2)*2)+1] = 'f8'

	catalogue_file = output_file_name
	cat_file_full = np.loadtxt(catalogue_file)
	ID_numbers = cat_file_full[:,0]
	redshifts = cat_file_full[:,1]
	aperture_radius = cat_file_full[:,2]
	aperture_correction = cat_file_full[:,3]
	n_objects = ID_numbers.size
	
	apparent_flux = np.zeros([final_number_filters, n_objects])
	apparent_flux_err = np.zeros([final_number_filters, n_objects])
	for j in range(0, final_number_filters):
		apparent_flux[:][j] = cat_file_full[:,(j+2)*2]
		apparent_flux_err[:][j] = cat_file_full[:,((j+2)*2)+1]
		
	# And now let's assemble the data array
	output_data = np.zeros([n_objects, 4+(final_number_filters*2)])
	output_data[:,0] = ID_numbers
	output_data[:,1] = redshifts
	output_data[:,2] = aperture_radius
	output_data[:,3] = aperture_correction
	for j in range(0, final_number_filters):
		output_data[:,(j+2)*2] = apparent_flux[:][j]
		output_data[:,((j+2)*2)+1] = apparent_flux_err[:][j]
		
	# And finally, let's write out the output file.
	outtab = Table(output_data, names=colnames, dtype=dtype)
	outtab.write(filenameroot+'.fits')



