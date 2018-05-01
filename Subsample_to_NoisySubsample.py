import os
import sys
import math
import argparse
import numpy as np
import scipy.special as sc
import matplotlib.pyplot as plt
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

# NIRCam photometric values
nfilters = 9 # 
n_aper = 4   # number_of_apertures__set_to_0_for_adaptive
filt_name = ['NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W']
bkgnd_nJy = np.array([6.39, 5.92, 4.97, 3.79, 18.76, 13.82, 14.81, 21.72, 31.60])     # backgrounds_for_GOODS_S_in_nJy/pix */
e_sec_nJy = np.array([0.0268, 0.0297, 0.0392, 0.0452, 0.0375, 0.0164, 0.0433, 0.0177, 0.0371])   # e/sec/nJy_from_as-built_spreadsheet */ 
medium_time = np.array([11700., 11700., 9100., 7500., 9200., 6700., 7500., 11700., 11700.]) # Deep_survey_exp_times_per_filter */
deep_time = np.array([60500., 80800., 59300., 38800., 49100., 30900., 38800., 60500., 60500.]) # Deep_survey_exp_times_per_filter */
aperture = np.array([0.16, 0.24, 0.32, 0.64])   # photom_apertures_in_arcsec  */

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

# HST Depth
parser.add_argument(
  '-hst','--hst_depth',
  help="HST survey: deeppointsource, deepextended, or flankingpointsource",
  action="store",
  type=str,
  dest="hst_depth",
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

# Go through and see which NIRCam filters are the ones in our survey
#number_matching_survey = 0
#for j in range(0, number_filters):
#	for k in range(0, nfilters):
#		if (filters[j] == filt_name[k]):
#			# The number of objects in the input subsample file that are in our survey
#			number_matching_survey = number_matching_survey + 1


# Parse HST depth
if (args.hst_depth == 'deeppointsource'):
	hst_depth_header = 'deeppointsource'
elif (args.hst_depth == 'deepextended'):
	hst_depth_header = 'deepextended'
elif (args.hst_depth == 'flankingpointsource'):
	hst_depth_header = 'flankingpointsource'
else:
	sys.exit("HST depth should be deeppointsource, deepextended, or flankingpointsource")

# Parse NIRCam depth
if (args.nircam_depth == 'deep'):
	time = deep_time
elif (args.nircam_depth == 'medium'):
	time = medium_time
else:
	sys.exit("NIRCam depth should be deep or medium")

output_file_name = args.output_file

if args.input_file.endswith('.fits'):
	#sys.exit("this is a fits file!")
	fitsinput = fits.open(args.input_file)
	ID_numbers = fitsinput[1].data['ID']
	redshifts = fitsinput[1].data['redshift']
	logmasses = fitsinput[1].data['logmass']
	re_major = fitsinput[1].data['re_major']
	sersic_n = fitsinput[1].data['sersic_n']
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
	logmasses = cat_file_full[:,2]
	re_major = cat_file_full[:,3]
	sersic_n = cat_file_full[:,4]
	n_objects = ID_numbers.size

	# Get the star-forming fluxes. 
	#print number_filters
	#print n_objects
	apparent_flux = np.zeros([number_filters, n_objects])
	for j in range(0, number_filters):
		apparent_flux[:][j] = cat_file_full[:,5+j]
	
#filt_flux = np.zeros([nfilters, n_objects])
#filt_sig = np.zeros([nfilters, n_objects])

# Open up the HST noise file
noise_file = '/Users/knh/Desktop/NIRCam/photometric_redshifts/prep_code/old_codes/grogin_5_sigma_sensitivities.txt'
tab = ascii.read(noise_file, names=['filter_names','deeppointsource','deepextended','flankingpointsource','ps_test1','ps_test2','ps_test3'], data_start=0)
filter_names = tab['filter_names'] 

hst_noise_five_sigma = tab[hst_depth_header] 

hst_noise_five_sigma_flux = 10**((-2.0/5.0) * (hst_noise_five_sigma + 48.60))
hst_noise_one_sigma_flux = hst_noise_five_sigma_flux/1e-23 / 1e-9 / 5.0

f = open(filenameroot+'.dat', 'a')
f.write('# ID   redshift   ')
for j in range(0, number_filters):
	f.write(filters[j]+' err   ')
f.write(' \n')

#pix_area = np.zeros(number_matching_survey)
#ap_corr = np.zeros(n_objects)
#filt_bkgnd = np.zeros([number_matching_survey, n_objects_flux])

# Go through all of the galaxies. 
for x in range(0, n_objects):

	flux_value = 0
	flux_value_noisy = 0
	
	f.write(str(ID_numbers[x])+' '+str(redshifts[x])+'   ')
	# Determine, based on the re_maj value, what the aperture that we should use to 
	# contain the flux in the object.	
	if (re_major[x] <= aperture[0]):
		i_ap = 0
	
	for j in range(1, n_aper):
		if ((aperture[j-1] < re_major[x]) & (re_major[x] <= aperture[j])):
			i_ap = j
	
	if ((aperture[n_aper-1] < re_major[x])):
		i_ap = n_aper - 1
			
	ap_radius = aperture[i_ap]

	ap_corr = sersic(sersic_n[x], re_major[x], ap_radius)
	
	if (ap_corr > 1.0):
		ap_corr = 1.0
	
	final_number_filters = 0
	final_filters = np.empty([0], dtype='S20')
	# Go through the various filters
	for j in range(0, number_filters):
		filter_found = 0
		
		# HST FILTERS AND FLUX ERROR ESTIMATES
		if (filters[j] == 'HST_F435W'):
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['HST_F435W'])
			filter_found = 1
			flux_value = apparent_flux[j][x]
			flux_error = hst_noise_one_sigma_flux[0]
			if (flux_value < 0):
				flux_value_noisy = np.random.normal(0.0, flux_error, 1)[0]
			else: 
				flux_value_noisy = np.random.normal(flux_value, flux_error, 1)[0]

		if (filters[j] == 'HST_F606W'):
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['HST_F606W'])
			filter_found = 1
			flux_value = apparent_flux[j][x]
			flux_error = hst_noise_one_sigma_flux[1]
			if (flux_value < 0):
				flux_value_noisy = np.random.normal(0.0, flux_error, 1)[0]
			else: 
				flux_value_noisy = np.random.normal(flux_value, flux_error, 1)[0]

		if (filters[j] == 'HST_F775W'):
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['HST_F775W'])
			filter_found = 1
			flux_value = apparent_flux[j][x]
			flux_error = hst_noise_one_sigma_flux[2]
			if (flux_value < 0):
				flux_value_noisy = np.random.normal(0.0, flux_error, 1)[0]
			else: 
				flux_value_noisy = np.random.normal(flux_value, flux_error, 1)[0]

		if (filters[j] == 'HST_F814W'):
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['HST_F814W'])
			filter_found = 1
			flux_value = apparent_flux[j][x]
			flux_error = hst_noise_one_sigma_flux[3]
			if (flux_value < 0):
				flux_value_noisy = np.random.normal(0.0, flux_error, 1)[0]
			else: 
				flux_value_noisy = np.random.normal(flux_value, flux_error, 1)[0]

		if (filters[j] == 'HST_F850LP'):
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['HST_F850LP'])
			filter_found = 1
			flux_value = apparent_flux[j][x]
			flux_error = hst_noise_one_sigma_flux[4]
			if (flux_value < 0):
				flux_value_noisy = np.random.normal(0.0, flux_error, 1)[0]
			else: 
				flux_value_noisy = np.random.normal(flux_value, flux_error, 1)[0]

		pix_area_SW = 3.14159265 * ap_radius * ap_radius / (.032 * .032 )
		pix_area_LW = 3.14159265 * ap_radius * ap_radius / (.064 * .064 )
		# NIRCam FILTERS AND FLUX ERROR ESTIMATES
		if (filters[j] == 'NRC_F090W'):
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F090W'])
			filter_found = 1
			nrc_filter_index = 0
			flux_value = apparent_flux[j][x]
			filt_bkgnd = pix_area_SW * bkgnd_nJy[nrc_filter_index]
			flux = ap_corr * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[nrc_filter_index]
			ferror = np.sqrt(ftemp * time[nrc_filter_index] + (6.0*pix_area_SW)*(6.0*pix_area_SW))							
			ferror = ferror / (e_sec_nJy[nrc_filter_index] * time[nrc_filter_index])
			ftemp = ferror * np.random.normal()
			flux_value_noisy = ap_corr * flux_value + ftemp
			flux_error = ferror

		if (filters[j] == 'NRC_F115W'):
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F115W'])
			filter_found = 1
			nrc_filter_index = 1
			flux_value = apparent_flux[j][x]
			filt_bkgnd = pix_area_SW * bkgnd_nJy[nrc_filter_index]
			flux = ap_corr * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[nrc_filter_index]
			ferror = np.sqrt(ftemp * time[nrc_filter_index] + (6.0*pix_area_SW)*(6.0*pix_area_SW))							
			ferror = ferror / (e_sec_nJy[nrc_filter_index] * time[nrc_filter_index])
			ftemp = ferror * np.random.normal()
			flux_value_noisy = ap_corr * flux_value + ftemp
			flux_error = ferror

		if (filters[j] == 'NRC_F150W'):
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F150W'])
			filter_found = 1
			nrc_filter_index = 2
			flux_value = apparent_flux[j][x]
			filt_bkgnd = pix_area_SW * bkgnd_nJy[nrc_filter_index]
			flux = ap_corr * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[nrc_filter_index]
			ferror = np.sqrt(ftemp * time[nrc_filter_index] + (6.0*pix_area_SW)*(6.0*pix_area_SW))							
			ferror = ferror / (e_sec_nJy[nrc_filter_index] * time[nrc_filter_index])
			ftemp = ferror * np.random.normal()
			flux_value_noisy = ap_corr * flux_value + ftemp
			flux_error = ferror

		if (filters[j] == 'NRC_F200W'):
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F200W'])
			filter_found = 1
			nrc_filter_index = 3
			flux_value = apparent_flux[j][x]
			filt_bkgnd = pix_area_SW * bkgnd_nJy[nrc_filter_index]
			flux = ap_corr * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[nrc_filter_index]
			ferror = np.sqrt(ftemp * time[nrc_filter_index] + (6.0*pix_area_SW)*(6.0*pix_area_SW))							
			ferror = ferror / (e_sec_nJy[nrc_filter_index] * time[nrc_filter_index])
			ftemp = ferror * np.random.normal()
			flux_value_noisy = ap_corr * flux_value + ftemp
			flux_error = ferror

		if (filters[j] == 'NRC_F277W'):
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F277W'])
			filter_found = 1
			nrc_filter_index = 4
			flux_value = apparent_flux[j][x]
			filt_bkgnd = pix_area_LW * bkgnd_nJy[nrc_filter_index]
			flux = ap_corr * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[nrc_filter_index]
			ferror = np.sqrt(ftemp * time[nrc_filter_index] + (9.0*pix_area_LW)*(9.0*pix_area_LW))							
			ferror = ferror / (e_sec_nJy[nrc_filter_index] * time[nrc_filter_index])
			ftemp = ferror * np.random.normal()
			flux_value_noisy = ap_corr * flux_value + ftemp
			flux_error = ferror

		if (filters[j] == 'NRC_F335M'):
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F335M'])
			filter_found = 1
			nrc_filter_index = 5
			flux_value = apparent_flux[j][x]
			filt_bkgnd = pix_area_LW * bkgnd_nJy[nrc_filter_index]
			flux = ap_corr * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[nrc_filter_index]
			ferror = np.sqrt(ftemp * time[nrc_filter_index] + (9.0*pix_area_LW)*(9.0*pix_area_LW))						
			ferror = ferror / (e_sec_nJy[nrc_filter_index] * time[nrc_filter_index])
			ftemp = ferror * np.random.normal()
			flux_value_noisy = ap_corr * flux_value + ftemp
			flux_error = ferror

		if (filters[j] == 'NRC_F356W'):
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F356W'])
			filter_found = 1
			nrc_filter_index = 6
			flux_value = apparent_flux[j][x]
			filt_bkgnd = pix_area_LW * bkgnd_nJy[nrc_filter_index]
			flux = ap_corr * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[nrc_filter_index]
			ferror = np.sqrt(ftemp * time[nrc_filter_index] + (9.0*pix_area_LW)*(9.0*pix_area_LW))						
			ferror = ferror / (e_sec_nJy[nrc_filter_index] * time[nrc_filter_index])
			ftemp = ferror * np.random.normal()
			flux_value_noisy = ap_corr * flux_value + ftemp
			flux_error = ferror

		if (filters[j] == 'NRC_F410M'):
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F410M'])
			filter_found = 1
			nrc_filter_index = 7
			flux_value = apparent_flux[j][x]
			filt_bkgnd = pix_area_LW * bkgnd_nJy[nrc_filter_index]
			flux = ap_corr * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[nrc_filter_index]
			ferror = np.sqrt(ftemp * time[nrc_filter_index] + (9.0*pix_area_LW)*(9.0*pix_area_LW))					
			ferror = ferror / (e_sec_nJy[nrc_filter_index] * time[nrc_filter_index])
			ftemp = ferror * np.random.normal()
			flux_value_noisy = ap_corr * flux_value + ftemp
			flux_error = ferror

		if (filters[j] == 'NRC_F444W'):
			final_number_filters = final_number_filters + 1
			final_filters = np.append(final_filters, ['NRC_F444W'])
			filter_found = 1
			nrc_filter_index = 8
			flux_value = apparent_flux[j][x]
			filt_bkgnd = pix_area_LW * bkgnd_nJy[nrc_filter_index]
			flux = ap_corr * flux_value + filt_bkgnd
			ftemp = flux * e_sec_nJy[nrc_filter_index]
			ferror = np.sqrt(ftemp * time[nrc_filter_index] + (9.0*pix_area_LW)*(9.0*pix_area_LW))							
			ferror = ferror / (e_sec_nJy[nrc_filter_index] * time[nrc_filter_index])
			ftemp = ferror * np.random.normal()
			flux_value_noisy = ap_corr * flux_value + ftemp
			flux_error = ferror

		if(filter_found == 1):
			#print filters[j],str(flux_value),str(flux_value_noisy),str(flux_error)
			f.write( '%.3f %.3f  ' % (flux_value_noisy, flux_error))

	f.write(' \n')
			
f.close()

print final_number_filters
print final_filters

if (args.make_fits):
	# First, let's make the  dtype and colnames arrays
	colnames = np.zeros(2+(final_number_filters*2), dtype ='S20')
	dtype = np.zeros(2+(final_number_filters*2), dtype ='str')
	colnames[0] = 'ID'
	colnames[1] = 'redshift'

	dtype[0] = 'I'
	dtype[1] = 'd'
	for j in range(0, final_number_filters):
		colnames[(j+1)*2] = final_filters[j]
		colnames[((j+1)*2)+1] = final_filters[j]+'_err'
		dtype[(j+1)*2] = 'd'
		dtype[((j+1)*2)+1] = 'd'

	catalogue_file = output_file_name
	cat_file_full = np.loadtxt(catalogue_file)
	ID_numbers = cat_file_full[:,0]
	redshifts = cat_file_full[:,1]
	n_objects = ID_numbers.size
	
	apparent_flux = np.zeros([final_number_filters, n_objects])
	apparent_flux_err = np.zeros([final_number_filters, n_objects])
	for j in range(0, number_filters):
		apparent_flux[:][j] = cat_file_full[:,(j+1)*2]
		apparent_flux_err[:][j] = cat_file_full[:,((j+1)*2)+1]
		
	# And now let's assemble the data array
	output_data = np.zeros([n_objects, 2+(final_number_filters*2)])
	output_data[:,0] = ID_numbers
	output_data[:,1] = redshifts
	for j in range(0, final_number_filters):
		output_data[:,(j+1)*2] = apparent_flux[:][j]
		output_data[:,((j+1)*2)+1] = apparent_flux_err[:][j]
		
	# And finally, let's write out the output file.
	outtab = Table(output_data, names=colnames, dtype=dtype)
	outtab.write(filenameroot+'.fits')



