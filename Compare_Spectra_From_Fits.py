import os
import sys
import math
import time
import argparse
import numpy as np
import matplotlib.pyplot as plt
import astropy.stats
import scipy.stats
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from Spectrum_Compare_Functions import *


# Which version of the mock are you using
JAGUAR_version = 'r1_v1.1'

def FluxtoABMag(flux):
	return (-5.0 / 2.0) * np.log10(flux) - 48.60

# speed of light in angstroms per second
c = 2.99792458e18

noisy_jades_filters = ['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 'HST_F850LP', 'NRC_F070W', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W']
noisy_jades_wave = np.array([0.4297, 0.5907, 0.7764, 0.8057, 0.9033, 0.704, 0.902, 1.154, 1.501, 1.989, 2.762, 3.362, 3.568, 4.082, 4.408])
number_jades_filters = len(noisy_jades_filters)


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
  '-eazy','--eazy_output_folder',
  help="Analyze EAZY output?",
  action="store",
  type=str,
  dest="eazy_output_folder",
  required=False
)

# Generate BPZ input file
parser.add_argument(
  '-bpz','--bpz_output_folder',
  help="Analyze BPZ output?",
  action="store",
  type=str,
  dest="bpz_output_folder",
  required=False
)

# Generate Le Phare input file
parser.add_argument(
  '-lep','--lephare_output_folder',
  help="Analyze Le Phare output?",
  action="store",
  type=str,
  dest="lephare_output_folder",
  required=False
)

# ID number
parser.add_argument(
  '-id','--id_number',
  help="ID Number?",
  action="store",
  type=int,
  dest="id_number",
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

if (args.id_number):
	ID_number = args.id_number
else:
	sys.exit("Need to specify an ID number! ")


print "Opening up Input Photometry"
# Open up the input photometry 
# args.input_photometry
if (args.input_photometry.endswith('.fits')):
	fitsinput = fits.open(args.input_photometry)
	ID_values = fitsinput[1].data['ID']
	redshifts = fitsinput[1].data['redshift']
	number_objects = ID_values.size
	flux_value_cat = np.zeros([number_objects, number_jades_filters])
	flux_value_err_cat = np.zeros([number_objects, number_jades_filters])

	for j in range(0, number_objects):
		for y in range(0, number_jades_filters):
			flux_value_cat[j,y] = fitsinput[1].data[noisy_jades_filters[y]][j]
			flux_value_err_cat[j,y] = fitsinput[1].data[noisy_jades_filters[y]+'_err'][j]

else:
	full_input_file = np.loadtxt(args.input_photometry)
	ID_values = full_input_file[:,0]
	redshifts = full_input_file[:,1]
	number_objects = len(ID_values)
	flux_value_cat = np.zeros([number_objects, number_jades_filters])
	flux_value_err_cat = np.zeros([number_objects, number_jades_filters])

	for j in range(0, number_objects):
		for y in range(0, number_jades_filters):
			flux_value_cat[j,y] = full_input_file[j,(y+1)*2]
			flux_value_err_cat[j,y] = full_input_file[j,((y+1)*2)+1]

# Look at the output from the Photo-Z programs. 

# EAZY
if (args.eazy_output_folder):
	print "EAZY Analysis on object: "+str(ID_number)
	version = 'EAZY'
	output_folder = optional_output_folder+'EAZY_analysis/'
	#if not os.path.exists(output_folder):
	#    os.makedirs(output_folder)

	output_zs = np.loadtxt(args.eazy_output_folder+'photz.zout') 
	ids = output_zs[:,0]
	z_spec = output_zs[:,1]
	z_phot = output_zs[:,13]

	title_for_plot = 'EAZY Results'

	#output_directory = '/Users/knh/Desktop/NIRCam/photometric_redshifts/eazy-photoz/early_tests/mock_v1.1_tests/mock_v1.1_5_8_18_outliers/OUTPUT/'

	tempfilt, z_grid, obs_sed, templam, temp_sed = generate_sed_arrays(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY=args.eazy_output_folder, CACHE_FILE='Same')
        
	#idx = 4864
	#ID_number = 76901

	idx_value = np.where(ids == ID_number)[0]
	if not idx_value:
	  sys.exit("Need to specify a valid ID number! ")
	else:
		idx = idx_value[0]
	plt.errorbar(tempfilt['lc'], tempfilt['fnu'][:,idx]/1e-23/1e-9,
		tempfilt['efnu'][:,idx]/1e-23/1e-9, color='green', markersize='60', alpha = 0.8, 
		linestyle='None', zorder=1000)
#	plt.errorbar(tempfilt['lc'], tempfilt['fnu'][:,idx]/1e-23/1e-9,
#		tempfilt['efnu'][:,idx]/1e-23/1e-9, marker='.', color='green', alpha = 0.8, 
#		linestyle='None', zorder=1000, label='Data')

	# GETTING THE TRUE SPECTRUM
	start_time = time.time()
	print "Starting to get true spectrum. "
	print int(ID_values[idx])
	true_wavelength, spectrum_obj = return_true_spectrum(int(ID_values[idx]))
	end_time = time.time()
	delta_time = end_time - start_time
	print "Done! It took: "+str(delta_time)

	                  	
	#
	# START PLOTTING
	# 
	plt.scatter(tempfilt['lc'], obs_sed[:,idx]/1e-23/1e-9, color='r', alpha=0.8, s = 60, label='Template Photometry')
	plt.plot(templam*(1+z_grid[idx]), temp_sed[:,idx]/1e-23/1e-9, color='b', alpha=0.8, label='Template Spectrum')
	plt.scatter(noisy_jades_wave*10000, flux_value_cat[idx,:], color = 'green', s = 60, label='Noisy Input Flux', alpha = 0.8)
	
	# PLOT THE TRUE SPECTRUM	
	sf_spec_array_fnu = (spectrum_obj/(1.0+z_spec[idx]) * ((1.0+z_spec[idx]) * true_wavelength)**2 / c) 
	sf_wave_z = np.array((1.0+z_spec[idx]) * true_wavelength)
	sf_flux_nJy = np.array(sf_spec_array_fnu/1e-23/1e-9)
	plt.plot(sf_wave_z, sf_flux_nJy, color = 'orange', label = 'Input Spectrum', alpha = 0.6)
	
	# GET THE MIN/MAX VALUES
	min_flux = np.min(flux_value_cat[idx,:])
	if (min_flux < 1e-2):
		min_flux = 1e-2
	max_flux = np.max(flux_value_cat[idx,:])
	y_axis_difference = max_flux - min_flux
	
	plt.text(3200, min_flux + 1.8*y_axis_difference, 'ID = '+str(int(ID_values[idx])))
	plt.text(3200, min_flux + 1.2*y_axis_difference, 'z$_{spec}$ = '+str(z_spec[idx]))
	plt.text(3200, min_flux + 0.8*y_axis_difference, 'z$_{phot}$ = '+str(z_phot[idx]))
	plt.loglog()
	plt.xlim(3000,6.e4)
	plt.ylim(0.5*min_flux, max_flux+2*max_flux)
	plt.legend(loc='lower right')
	plt.xlabel(r'Observed wavelength, $\mathrm{\AA}$')
	plt.ylabel(r'Observed flux, $f_\nu$ @ PRIOR_ABZP')
	plt.tight_layout()
	plt.show()

# BPZ 
if (args.bpz_output_folder):

	bpz_template_path = '/Users/knh/Desktop/NIRCam/photometric_redshifts/BPZ/bpz-1.99.3/SED/'
	bpz_template_file = 'new_templates.list'
	bpz_filters = ['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 'HST_F850LP', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W']
	bpz_wave = np.array([0.4297, 0.5907, 0.7764, 0.8057, 0.9033, 0.902, 1.154, 1.501, 1.989, 2.762, 3.362, 3.568, 4.082, 4.408])
	n_bpz_filters = len(bpz_wave)
	
	idx = 39

	# Open up the templates 
	full_list_of_templates = np.loadtxt(bpz_template_path+bpz_template_file, dtype = str)
	template_name = full_list_of_templates[:]
	number_templates = len(template_name)

	output_zs = np.loadtxt(args.bpz_output_folder+'.bpz') 
	ids = output_zs[:,0]
	z_phot = output_zs[:,1]
	template_number = output_zs[:,4]
	z_spec = output_zs[:,9]
	
	ID_index = np.where(ids == idx)[0][0]
	
	# Open up flux comparison. 
	# BPZ_input_5_8_18.flux_comparison

	flux_comparison = np.loadtxt(args.bpz_output_folder+'.flux_comparison') 
	id_fluxcomp = int(flux_comparison[ID_index,0])
	mult_factor = flux_comparison[ID_index,4]
	template_fluxes = np.zeros(n_bpz_filters)
	for j in range(0, n_bpz_filters):
		template_fluxes[j] = flux_comparison[ID_index,j+5]

	# I don't exactly know what the units are here. 
	multiplied_fluxes = mult_factor * template_fluxes

	# Maybe we should also multiply by a factor of 1.00e-22 ? 
	test_factor = 1.00e-22 
	test_multiplied_fluxes = multiplied_fluxes * test_factor

	# Get the correct template	
	if ((int(template_number[ID_index]) - template_number[ID_index]) == 0):
		template_index = int(template_number[ID_index])
		template_file = np.loadtxt(bpz_template_path+template_name[template_index-1]) 
		temp_wave_final = template_file[:,0]
		temp_sed_final = template_file[:,1]		
	else:
		template_one_index = int(template_number[ID_index])
		template_two_index = int(template_number[ID_index]) + 1
		template_one_scale = template_two_index*1.0 - template_number[ID_index]
		template_two_scale = 1 - template_one_scale
		template_one_file = np.loadtxt(bpz_template_path+template_name[template_one_index-1]) 
		temp_one_wave = template_one_file[:,0]
		temp_one_sed = template_one_file[:,1]
		template_two_file = np.loadtxt(bpz_template_path+template_name[template_two_index-1]) 
		temp_two_wave = template_two_file[:,0]
		temp_two_sed = template_two_file[:,1]

		temp_one_sed_interp = np.interp(temp_two_wave, temp_one_wave, temp_one_sed)
		temp_wave_final = temp_one_wave * (1+z_phot[ID_index])
		temp_sed_final = (temp_one_sed_interp * template_one_scale) + (temp_two_sed * template_two_scale)
		temp_sed_final_normalized = temp_sed_final * mult_factor * test_factor * 1e-8
	
	#plt.plot(temp_wave_final, temp_sed_final_normalized)
    #plt.show()
	  
