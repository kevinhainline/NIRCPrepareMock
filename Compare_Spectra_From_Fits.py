import os
import sys
import math
import time
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy.stats
import scipy.stats
from scipy.interpolate import CubicSpline
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from Spectrum_Compare_Functions import *


# Which version of the mock are you using
JAGUAR_version = 'r1_v1.1'

def flambda_to_fnu(waveang,flambda):
	return (waveang**2.0) / 2.99792458e+18 * flambda

def FluxtoABMag(flux):
	return (-5.0 / 2.0) * np.log10(flux) - 48.60

def ABMagtoFlux(abmag):
	return np.power(10, (abmag + 48.60)/-2.5)


# speed of light in angstroms per second
c = 2.99792458e18

#noisy_jades_filters = np.array(['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 'HST_F850LP', 'NRC_F070W', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W'], dtype = 'str')
#noisy_jades_wave = np.array([0.4297, 0.5907, 0.7764, 0.8057, 0.9033, 0.704, 0.902, 1.154, 1.501, 1.989, 2.762, 3.362, 3.568, 4.082, 4.408])
noisy_jades_filters = np.array(['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 'HST_F850LP', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W'], dtype = 'str')
noisy_jades_wave = np.array([0.4297, 0.5907, 0.7764, 0.8057, 0.9033, 0.902, 1.154, 1.501, 1.989, 2.762, 3.362, 3.568, 4.082, 4.408])
number_jades_filters = len(noisy_jades_filters)

# Medium Filters
#filters_to_use =  np.array(['NRC_F070W', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W'], dtype = 'str')
# Deep Filters
filters_to_use =  np.array(['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 'HST_F850LP', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W'], dtype = 'str')
filters_to_use_indices = np.zeros(len(filters_to_use), dtype = 'int')
for z in range(0, len(filters_to_use)):
	filters_to_use_indices[z] = np.where(noisy_jades_filters == filters_to_use[z])[0][0]


normalization_filter = 'NRC_F200W'
for j in range(0, number_jades_filters):
	if (noisy_jades_filters[j] == normalization_filter):
		normalization_index = j

# Parser for figuring out the various input parameters
parser = argparse.ArgumentParser()

######################
# Required Arguments #
######################

# Input folder
parser.add_argument(
  '-in','--inputfolder',
  help="Input Folder With Mock Catalog Files",
  action="store",
  type=str,
  dest="input_folder",
  required=True
)

# Noisy Photometry
parser.add_argument(
  '-ni','--noisy_input',
  help="Input Photometry for analysis",
  action="store",
  type=str,
  dest="noisy_input",
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

# ID list
parser.add_argument(
  '-idlist','--id_number_list',
  help="List of ID Numbers?",
  action="store",
  type=str,
  dest="id_number_list",
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
	ID_numbers = np.zeros(1)
	ID_numbers[0] = args.id_number
	number_input_objects = len(ID_numbers)
	if (args.id_number_list):
		"You can't specify an individual ID and a list, ignoring the list."

if not (args.id_number):
	if not (args.id_number_list):
		sys.exit("Need to specify an ID number or a list of ID numbers!")

if (args.id_number_list):
	ID_input_file = np.loadtxt(args.id_number_list)
	ID_numbers = ID_input_file[:,0]
	list_z_spec = ID_input_file[:,1]
	list_z_phot = ID_input_file[:,2]
	number_input_objects = len(ID_numbers)
	avg_z_spec = np.mean(list_z_spec)
	avg_z_phot = np.mean(list_z_phot)

print "Opening up Catalog Photometry"
mock_filters = ['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 
	'HST_F850LP', 'NRC_F070W', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W',
	'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 
	'NRC_F410M', 'NRC_F444W']

mock_filter_wavelength = np.array([0.4297, 0.5907, 0.7764, 0.8057, 0.9033, 
	0.704, 0.902, 1.154, 1.501, 1.989, 2.762, 3.362, 3.568, 4.082, 4.408])


number_mock_filters = len(mock_filters)

# Full star-forming catalogue, fits file
full_sf_mock_file = args.input_folder+'JADES_SF_mock_'+JAGUAR_version+'.fits'
sftestfits = fits.open(full_sf_mock_file)

# Get redshifts and galaxy properties
full_sf_IDs = sftestfits[1].data['ID']
full_sf_redshifts = sftestfits[1].data['redshift']
full_sf_logmass = sftestfits[1].data['mStar']
full_sf_re_maj = sftestfits[1].data['Re_maj']
full_sf_sersic_n = sftestfits[1].data['sersic_n']
n_sf_objects = len(full_sf_IDs)

# Get the star-forming fluxes. 
sf_filt_flux = np.zeros([number_mock_filters, n_sf_objects])
for j in range(0, number_mock_filters):
	sf_filt_flux[:][j] = sftestfits[1].data[mock_filters[j]+'_fnu']#*1e-23*1e-9

# Full quiescent catalogue, fits file
full_q_mock_file = args.input_folder+'JADES_Q_mock_'+JAGUAR_version+'.fits'
qtestfits = fits.open(full_q_mock_file)

# Get redshifts and galaxy properties
full_q_IDs = qtestfits[1].data['ID']
full_q_redshifts = qtestfits[1].data['redshift']
full_q_logmass = qtestfits[1].data['mStar']
full_q_re_maj = qtestfits[1].data['Re_maj']
full_q_sersic_n = qtestfits[1].data['sersic_n']
n_q_objects = len(full_q_IDs)

# Get the quiescent fluxes. 
q_filt_flux = np.zeros([number_mock_filters, n_q_objects])
for j in range(0, number_mock_filters):
	q_filt_flux[:][j] = qtestfits[1].data[mock_filters[j]+'_fnu']#*1e-23*1e-9

full_all_IDs = np.append(full_sf_IDs, full_q_IDs)
full_all_redshifts = np.append(full_sf_redshifts, full_q_redshifts)
full_all_logmass = np.append(full_sf_logmass, full_q_logmass)
full_all_re_maj = np.append(full_sf_re_maj, full_q_re_maj)
all_filt_flux = np.zeros([number_mock_filters, n_sf_objects+n_q_objects])
for j in range(0, number_mock_filters):
	all_filt_flux[:][j] = np.append(sf_filt_flux[:][j], q_filt_flux[:][j])


print "Opening up Noisy Photometry"
# Open up the input photometry 
# args.input_photometry
if (args.noisy_input.endswith('.fits')):
	fitsinput = fits.open(args.noisy_input)
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
	full_input_file = np.loadtxt(args.noisy_input)
	ID_values = full_input_file[:,0]
	redshifts = full_input_file[:,1]
	number_objects = len(ID_values)
	flux_value_cat = np.zeros([number_objects, number_jades_filters])
	flux_value_err_cat = np.zeros([number_objects, number_jades_filters])

	print number_objects, number_jades_filters
	for j in range(0, number_objects):
		for y in range(0, number_jades_filters):
			flux_value_cat[j,y] = full_input_file[j,(y+1)*2]
			flux_value_err_cat[j,y] = full_input_file[j,((y+1)*2)+1]



# Look at the output from the Photo-Z programs. 

# EAZY
if (args.eazy_output_folder):

	output_zs = np.loadtxt(args.eazy_output_folder+'photz.zout') 
	ids = output_zs[:,0]
	z_spec = output_zs[:,1]
	z_phot = output_zs[:,13]

	print number_input_objects

	for z in range(0, number_input_objects):
		print "EAZY Analysis on object: "+str(int(ID_numbers[z]))
		version = 'EAZY'
		output_folder = optional_output_folder+'EAZY_analysis/'
	
		title_for_plot = 'EAZY Results'
		
		tempfilt, z_grid, obs_sed, templam, temp_sed, pz, chi2fit = generate_sed_arrays(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY=args.eazy_output_folder, CACHE_FILE='Same')
		idx_value = np.where(ids == ID_numbers[z])[0]
		mock_idx_value = np.where(full_all_IDs == ID_numbers[z])[0][0]
#		if not idx_value:
#		  sys.exit("Specified ID number ("+str(ID_numbers[z])+") not valid! ")
#		else:
		idx = idx_value[0]
		plt.errorbar(tempfilt['lc'], tempfilt['fnu'][:,idx]/1e-23/1e-9,
			tempfilt['efnu'][:,idx]/1e-23/1e-9, color='green', markersize='60', alpha = 0.8, 
			linestyle='None', zorder=1000)
			
		# GETTING THE TRUE SPECTRUM
		start_time = time.time()
		print "Getting true spectrum, OBJ ID "+str(ID_values[idx])
		true_wavelength, spectrum_obj = return_true_spectrum(int(ID_values[idx]))
		end_time = time.time()
		delta_time = end_time - start_time
		print "Done! It took: "+str(delta_time)
	

		#
		# START PLOTTING
		#
		# PLOT THE EAZY TEMPLATE PHOTOMETRY and SPECTRUM
		plt.clf() 
		plt.figure(figsize=(11,5))
		plt.subplot(1, 2, 1)
		plt.scatter(tempfilt['lc']/10000, obs_sed[:,idx]/1e-23/1e-9, color='green', edgecolors='none', s = 50, label='Template Photometry', zorder = 10)
		plt.plot(templam*(1+z_grid[idx])/10000, temp_sed[:,idx]/1e-23/1e-9, color='g', alpha=0.5, label='Template Spectrum')
		plt.errorbar(noisy_jades_wave[filters_to_use_indices], flux_value_cat[idx,filters_to_use_indices]/1e-23/1e-9, flux_value_err_cat[idx,filters_to_use_indices]/1e-23/1e-9, ecolor = 'red', fmt='none')
		plt.scatter(noisy_jades_wave[filters_to_use_indices], flux_value_cat[idx,filters_to_use_indices]/1e-23/1e-9, color = 'red', edgecolors='none', s = 60, label='Noisy Input Flux', zorder = 9)

		# PLOT THE TRUE SPECTRUM
		sf_spec_array_fnu = (spectrum_obj/(1.0+z_spec[idx]) * ((1.0+z_spec[idx]) * true_wavelength)**2 / c) 
		sf_wave_z = np.array((1.0+z_spec[idx]) * true_wavelength)
		sf_flux_nJy = np.array(sf_spec_array_fnu/1e-23/1e-9)
		plt.plot(sf_wave_z/10000, sf_flux_nJy, color = 'blue', label = 'Input Spectrum', alpha = 0.5)
		#plt.plot(sf_wave_z/10000, sf_flux_nJy/1.5, color = 'blue', label = 'Aperture Loss Spectrum', alpha = 0.3)
		#plt.plot(sf_wave_z/10000, sf_flux_nJy/2, color = 'blue', alpha = 0.2)
		#plt.plot(sf_wave_z/10000, sf_flux_nJy/2.5, color = 'blue', alpha = 0.1)
		#plt.plot(sf_wave_z/10000, sf_flux_nJy/3, color = 'blue', alpha = 0.05)
		
		plt.scatter(mock_filter_wavelength, all_filt_flux[:,mock_idx_value], color = 'blue', edgecolors='none', s = 60, label='Mock Catalog Flux', zorder = 9, alpha = 0.4)
		
		# GET THE MIN/MAX VALUES
		min_flux = np.min(flux_value_cat[idx,:]/1e-23/1e-9)
		if (min_flux < 1e-2):
			min_flux = 1e-2
		max_flux = np.max(flux_value_cat[idx,:]/1e-23/1e-9)
		y_axis_difference = max_flux - min_flux
		
		plt.semilogy()
		plt.title('ID = '+str(int(ID_values[idx]))+', log(Mass) = '+str(np.round(full_all_logmass[mock_idx_value],2))+', Re_maj = '+str(full_all_re_maj[mock_idx_value])+' \'\'')
		plt.xlim(3000/10000,6.e4/10000)
		plt.ylim(0.5*min_flux, max_flux+2*max_flux)
		plt.legend(loc='lower right')
#		plt.xlabel(r'Observed wavelength, $\mathrm{\AA}$')
		plt.xlabel(r'Observed wavelength, Microns')
		plt.ylabel(r'Observed flux, $f_\nu$ (nJy)')
		#plt.tight_layout()
		plt.subplot(1, 2, 2)
		#plt.title(z$_{spec}$ = '+str(z_spec[idx])+', z$_{phot}$ = '+str(z_phot[idx]))
		plt.plot(tempfilt['zgrid'], chi2fit[:,idx])
		plt.xlabel(r'z$_{phot}$')
		plt.ylabel(r'$\chi^2$')
		#plt.axvspan(l68[index_value[0]], u68[index_value[0]], color = 'green', alpha = 0.05)
		#plt.axvspan(l95[index_value[0]], u95[index_value[0]], color = 'green', alpha = 0.1)
		#plt.axvspan(l99[index_value[0]], u99[index_value[0]], color = 'green', alpha = 0.15)
		plt.axvspan(z_spec[idx]-0.05, z_spec[idx]+0.05, color = 'blue')
		plt.axvspan(z_phot[idx]-0.05, z_phot[idx]+0.05, color = 'green')
		plt.text(z_spec[idx]+0.15, 0.8*max(chi2fit[:,idx]), 'z$_{spec}$ = '+str(z_spec[idx]), color = 'blue', rotation = 90, alpha = 0.5)
		plt.text(z_phot[idx]+0.15, 0.8*max(chi2fit[:,idx]), 'z$_{phot}$ = '+str(z_phot[idx]), color = 'green', rotation = 90, alpha = 0.5)
		
		plt.savefig(args.eazy_output_folder+str(int(ID_numbers[z]))+'_SED.png', dpi = 300)


	
# BPZ 
if (args.bpz_output_folder):

	bpz_template_path = '/Users/knh/Desktop/NIRCam/photometric_redshifts/BPZ/bpz-1.99.3/SED/'
	bpz_template_file = 'new_templates.list'
	bpz_filters = ['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 'HST_F850LP', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W']
	bpz_wave = np.array([0.4297, 0.5907, 0.7764, 0.8057, 0.9033, 0.902, 1.154, 1.501, 1.989, 2.762, 3.362, 3.568, 4.082, 4.408])
	#bpz_filters = ['NRC_F070W' 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W']
	#bpz_wave = np.array([0.7088, 0.902, 1.154, 1.501, 1.989, 2.762, 3.362, 3.568, 4.082, 4.408])
	n_bpz_filters = len(bpz_wave)
	
	# Open up the templates 
	full_list_of_templates = np.loadtxt(bpz_template_path+bpz_template_file, dtype = str)
	template_name = full_list_of_templates[:]
	number_templates = len(template_name)

	output_zs = np.loadtxt(args.bpz_output_folder+'.bpz') 
	ids = output_zs[:,0]
	z_phot = output_zs[:,1]
	template_number = output_zs[:,4]
	z_spec = output_zs[:,9]

	# Open up the pz file
	pz_values_raw = np.loadtxt(args.bpz_output_folder+'.probs') 
	pz_redshifts = np.arange(0.2000,15.0100,0.0100)

	for z in range(0, number_input_objects):
		print "BPZ Analysis on object: "+str(int(ID_numbers[z]))
		version = 'BPZ'

		mock_idx_value = np.where(full_all_IDs == ID_numbers[z])[0][0]

		ID_index = np.where(ids == ID_numbers[z])[0][0]
		idx = ID_index
		print ID_index

		# Getting the probabilities
		pz_values = pz_values_raw[ID_index,:]
		pz_values = pz_values[1:len(pz_values)]
	
		# GETTING THE TRUE SPECTRUM
		#print "Getting true spectrum, OBJ ID "+str(ID_values[idx])
		true_wavelength, spectrum_obj = return_true_spectrum(int(ID_numbers[z]))
	
		# Open up flux comparison. 
		# BPZ_input_5_8_18.flux_comparison
		flux_comparison = np.loadtxt(args.bpz_output_folder+'.flux_comparison') 
		id_fluxcomp = int(flux_comparison[ID_index,0])
		template_fluxes = np.zeros(n_bpz_filters)
		for j in range(0, n_bpz_filters):
			template_ab_magnitude = -2.5*np.log10(flux_comparison[ID_index,j+5])
			template_fluxes[j] = ABMagtoFlux(template_ab_magnitude)/1e-23/1e-9

		# Get the correct template
		if ((int(template_number[ID_index]) - template_number[ID_index]) == 0.0):
			template_index = int(template_number[ID_index])
			print "One template fit the best: "+bpz_template_path+template_name[template_index-1]
			template_file = np.loadtxt(bpz_template_path+template_name[template_index-1]) 
			temp_wave_final = template_file[:,0]
			temp_sed = template_file[:,1]
			temp_sed_fnu = flambda_to_fnu(temp_wave_final,temp_sed)
			temp_AB_mag = FluxtoABMag(temp_sed_fnu) 
			temp_AB_mag_z = np.interp(temp_wave_final/(1+z_phot[ID_index]), temp_wave_final, temp_AB_mag)
			temp_fnu_final = ABMagtoFlux(temp_AB_mag_z)
		
			temp_fnu_final_bpz_filters = np.interp(bpz_wave*10000, temp_wave_final, temp_fnu_final)
			temp_fnu_final_bpz_filters_norm = temp_fnu_final_bpz_filters / template_fluxes
			temp_fnu_final_normalized = temp_fnu_final / np.mean(temp_fnu_final_bpz_filters_norm)
		
		else:
			template_one_number = int(template_number[ID_index])
			template_two_number = template_one_number + 1
			template_one_scale = template_two_number*1.0 - template_number[ID_index]
			template_two_scale = template_number[ID_index] - template_one_number#1 - template_one_scale
	
			print "Two templates fit the best: "+bpz_template_path+template_name[template_one_number-1]+" and "+bpz_template_path+template_name[template_two_number-1]
			template_one_file = np.loadtxt(bpz_template_path+template_name[template_one_number-1]) 
			temp_wave_final = template_one_file[:,0]
			temp_one_sed = template_one_file[:,1]
			temp_one_fnu = flambda_to_fnu(temp_wave_final,temp_one_sed)
			temp_one_AB_mag = FluxtoABMag(temp_one_fnu) 
			temp_one_AB_mag_z = np.interp(temp_wave_final/(1+z_phot[ID_index]), temp_wave_final, temp_one_AB_mag)
			
			template_two_file = np.loadtxt(bpz_template_path+template_name[template_two_number-1]) 
			temp_two_wave = template_two_file[:,0]
			temp_two_sed = template_two_file[:,1]
			temp_two_fnu = flambda_to_fnu(temp_two_wave,temp_two_sed)
			temp_two_AB_mag = FluxtoABMag(temp_two_fnu) 
			temp_two_AB_mag_z = np.interp(temp_two_wave/(1+z_phot[ID_index]), temp_two_wave, temp_two_AB_mag)
			
			temp_two_AB_mag_z_interp = np.interp(temp_wave_final,temp_two_wave,temp_two_AB_mag_z)
			temp_AB_final = (temp_one_AB_mag_z * template_one_scale) + (temp_two_AB_mag_z_interp * template_two_scale)
			temp_fnu_final = ABMagtoFlux(temp_AB_final)
			
			temp_fnu_final_bpz_filters = np.interp(bpz_wave*10000, temp_wave_final, temp_fnu_final)
			temp_fnu_final_bpz_filters_norm = temp_fnu_final_bpz_filters / template_fluxes
			temp_fnu_final_normalized = temp_fnu_final / np.mean(temp_fnu_final_bpz_filters_norm)
		
		#print temp_fnu_final_normalized

		
		#
		# START PLOTTING
		#
		# PLOT THE BPZ TEMPLATE PHOTOMETRY and SPECTRUM
		plt.figure(figsize=(11,5))
		plt.subplot(1, 2, 1)
		plt.scatter(bpz_wave, template_fluxes, color='green', edgecolors='none', s = 50, label='Template Photometry', zorder = 10)
		plt.plot(temp_wave_final/10000, temp_fnu_final_normalized, color='g', alpha=0.5, label='Template Spectrum')
		plt.errorbar(noisy_jades_wave[filters_to_use_indices], flux_value_cat[idx,filters_to_use_indices], flux_value_err_cat[idx,filters_to_use_indices], ecolor = 'red', fmt='none')
		plt.scatter(noisy_jades_wave[filters_to_use_indices], flux_value_cat[idx,filters_to_use_indices], color = 'red', edgecolors='none', s = 60, label='Noisy Input Flux', zorder = 9)
	
		# PLOT THE TRUE SPECTRUM	
		sf_spec_array_fnu = (spectrum_obj/(1.0+z_spec[idx]) * ((1.0+z_spec[idx]) * true_wavelength)**2 / c) 
		sf_wave_z = np.array((1.0+z_spec[idx]) * true_wavelength)
		sf_flux_nJy = np.array(sf_spec_array_fnu/1e-23/1e-9)
		plt.plot(sf_wave_z/10000, sf_flux_nJy, color = 'blue', label = 'Input Spectrum', alpha = 0.5)

		plt.scatter(mock_filter_wavelength, all_filt_flux[:,mock_idx_value], color = 'blue', edgecolors='none', s = 60, label='Mock Catalog Flux', zorder = 9, alpha = 0.4)
			
		# GET THE MIN/MAX VALUES
		min_flux = np.min(flux_value_cat[idx,:])
		if (min_flux < 1e-2):
			min_flux = 1e-2
		max_flux = np.max(flux_value_cat[idx,:])
		y_axis_difference = max_flux - min_flux
		
		#plt.text(3200, min_flux + 1.8*y_axis_difference, 'ID = '+str(int(ID_values[idx])))
		#plt.text(3200, min_flux + 1.2*y_axis_difference, 'z$_{spec}$ = '+str(z_spec[idx]))
		#plt.text(3200, min_flux + 0.8*y_axis_difference, 'z$_{phot}$ = '+str(z_phot[idx]))
		plt.semilogy()
		plt.title('ID = '+str(int(ID_values[idx]))+', log(Mass) = '+str(np.round(full_all_logmass[mock_idx_value],2))+', Re_maj = '+str(full_all_re_maj[mock_idx_value])+' \'\'')
		plt.xlim(3000/10000,6.e4/10000)
		plt.ylim(0.5*min_flux, max_flux+5*max_flux)
		plt.legend(loc='lower right')
		plt.xlabel(r'Observed wavelength, Microns')#$\mathrm{\AA}$')
		plt.ylabel(r'Observed flux, $f_\nu$ (nJy)')
	
		plt.subplot(1, 2, 2)
		plt.plot(pz_redshifts, pz_values)
		plt.xlabel(r'z$_{phot}$')
		plt.ylabel(r'p(z)')
		plt.axvspan(z_spec[idx]-0.05, z_spec[idx]+0.05, color = 'green', alpha = 0.2)
		plt.axvspan(z_phot[idx]-0.05, z_phot[idx]+0.05, color = 'blue', alpha = 0.2)
		plt.text(z_spec[idx]+0.15, 0.8*max(pz_values), 'z$_{spec}$ = '+str(z_spec[idx]), color = 'green', rotation = 90, alpha = 0.4)
		plt.text(z_phot[idx]+0.15, 0.8*max(pz_values), 'z$_{phot}$ = '+str(z_phot[idx]), color = 'blue', rotation = 90, alpha = 0.4)
	
		plt.savefig(args.bpz_output_folder+'_'+str(int(ID_numbers[z]))+'_SED.png', dpi = 300)
