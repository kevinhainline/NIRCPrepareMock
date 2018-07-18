import os
import sys
import math
import time
import argparse
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
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

def calculate_mean_median_offset(wave_array, values_array, bins):
	final_average_offsets = np.zeros(len(bins)-1)
	final_median_offsets = np.zeros(len(bins)-1)
	for w in range(0, len(bins)-1):
		final_average_offsets[w] = np.mean(values_array[np.where((wave_array > bins[w]) & (wave_array < bins[w+1]))[0]])
		final_median_offsets[w] = np.median(values_array[np.where((wave_array > bins[w]) & (wave_array < bins[w+1]))[0]])
	return final_average_offsets, final_median_offsets

# speed of light in angstroms per second
c = 2.99792458e18

noisy_jades_filters = np.array(['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 'HST_F850LP', 'NRC_F070W', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W'], dtype = 'str')
noisy_jades_wave = np.array([0.4297, 0.5907, 0.7764, 0.8057, 0.9033, 0.704, 0.902, 1.154, 1.501, 1.989, 2.762, 3.362, 3.568, 4.082, 4.408])
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

# ID list
parser.add_argument(
  '-idlist','--id_number_list',
  help="List of ID Numbers?",
  action="store",
  type=str,
  dest="id_number_list",
  required=False
)

# Ratio (instead of Flux Difference)
parser.add_argument(
  '-ra','--ratio',
  help="Calculate/Plot Flux Ratio instead of Flux Difference?",
  action="store_true",
  dest="ratio",
  required=False
)

# Logarithm of X Axis For Plots
parser.add_argument(
  '-logx','--logx',
  help="Use Logarithmic Scale for the X Axis on the Plots?",
  action="store_true",
  dest="logx",
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

if (args.id_number_list):
	ID_input_file = np.loadtxt(args.id_number_list)
	ID_numbers = ID_input_file[:,0]
	number_input_objects = len(ID_numbers)


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
	version = 'EAZY'
	title_for_plot = 'EAZY Results'

	output_zs = np.loadtxt(args.eazy_output_folder+'photz.zout') 
	ids = output_zs[:,0]
	z_spec = output_zs[:,1]
	z_phot = output_zs[:,13]
	total_number_of_objects = len(ids)

	photometric_difference = np.zeros([len(filters_to_use_indices),number_input_objects])
	spectroscopic_redshift = np.zeros(number_input_objects)
	wave_z = np.zeros([len(filters_to_use_indices),number_input_objects])
	tempfilt, z_grid, obs_sed, templam, temp_sed, pz, chi2fit = generate_sed_arrays(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY=args.eazy_output_folder, CACHE_FILE='Same')
	noisy_fluxes = flux_value_cat[idx,filters_to_use_indices]

	print "Going through all of the objects to calculate photometric offsets"
	for z in range(0, number_input_objects):
		
		idx_value = np.where(ids == ID_numbers[z])[0]
		idx = idx_value[0]
		
		template_photometry = obs_sed[:,idx]/1e-23/1e-9
		spectroscopic_redshift[z] = z_spec[idx]
		#photometric_difference[:,z] = (template_photometry - noisy_fluxes) / noisy_fluxes
		if (args.ratio):
			photometric_difference[:,z] = (template_photometry/noisy_fluxes)
		else:
			photometric_difference[:,z] = (template_photometry - noisy_fluxes) / noisy_fluxes
		wave_z[:,z] = noisy_jades_wave[filters_to_use_indices] / (1+spectroscopic_redshift[z])


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
	total_number_of_objects = len(ids)

	# Open up flux comparison. 
	flux_comparison = np.loadtxt(args.bpz_output_folder+'.flux_comparison') 
	#id_fluxcomp = int(flux_comparison[ID_index,0])

	photometric_difference = np.zeros([len(filters_to_use_indices),number_input_objects])
	spectroscopic_redshift = np.zeros(number_input_objects)
	wave_z = np.zeros([len(filters_to_use_indices),number_input_objects])
	for z in range(0, number_input_objects):

		ID_index = np.where(ids == ID_numbers[z])[0][0]
		idx = ID_index

		spectroscopic_redshift[z] = z_spec[idx]
		noisy_fluxes = flux_value_cat[idx,filters_to_use_indices]
		template_photometry = np.zeros(n_bpz_filters)
		for j in range(0, n_bpz_filters):
			template_ab_magnitude = -2.5*np.log10(flux_comparison[idx,j+5])
			template_photometry[j] = ABMagtoFlux(template_ab_magnitude)/1e-23/1e-9
		
		#print template_photometry
		#print noisy_fluxes
		#photometric_difference[:,z] = (template_photometry - noisy_fluxes) / noisy_fluxes
		if (args.ratio):
			photometric_difference[:,z] = (template_photometry/noisy_fluxes)
		else:
			photometric_difference[:,z] = (template_photometry - noisy_fluxes) / noisy_fluxes
		wave_z[:,z] = noisy_jades_wave[filters_to_use_indices] / (1+spectroscopic_redshift[z])




# NOW THIS	
median_values = np.zeros(len(filters_to_use_indices))
mean_values = np.zeros(len(filters_to_use_indices))
median_flux_waves = np.zeros(len(filters_to_use_indices))
mean_flux_waves = np.zeros(len(filters_to_use_indices))
for p in range(0, len(filters_to_use_indices)):
	print "Exploring Filter: "+noisy_jades_filters[filters_to_use_indices][p]
	lyman_break_wavelength = 0.1150 # 0.0912
	lyman_break_z = (noisy_jades_wave[filters_to_use_indices][p] / lyman_break_wavelength) - 1.0
	#print "Lyman break redshift = "+str(lyman_break_z)
	lyman_break_indices = np.where(spectroscopic_redshift > lyman_break_z)[0]
	good_indices = np.where(spectroscopic_redshift < lyman_break_z)[0]
	
	# Calculate the running averages
	#print np.min(wave_z[p,:][good_indices]), np.max(wave_z[p,:][good_indices])
	wave_bins = np.arange(np.min(wave_z[p,:][good_indices]), np.max(wave_z[p,:][good_indices]),0.05)
	wave_bins_centers = wave_bins[0:len(wave_bins)-1]+0.025
	final_average_offsets, final_median_offsets = calculate_mean_median_offset(wave_z[p,:][good_indices], photometric_difference[p,:][good_indices], wave_bins)

	median_values[p] = np.median(photometric_difference[p,:][good_indices])
	mean_values[p] = np.mean(photometric_difference[p,:][good_indices])
	plt.clf()
	fig = plt.figure(figsize=(8,6))
	ax1 = fig.add_subplot(111)

	cax = ax1.scatter(wave_z[p,:][good_indices], photometric_difference[p,:][good_indices], alpha = 0.3, edgecolors='none', c = spectroscopic_redshift[good_indices])
	cbar = fig.colorbar(cax, label = 'Redshift', orientation="vertical")
	ax1.scatter(wave_z[p,:][lyman_break_indices], photometric_difference[p,:][lyman_break_indices], s = 3, alpha = 0.1, edgecolors='none', c = 'black')
	ax1.scatter(wave_bins_centers, final_median_offsets, s = 20, color = 'green', label = 'Median Values')
	ax1.scatter(wave_bins_centers, final_average_offsets, s = 15, color = 'red', label = 'Mean Values')
	median_flux_waves[p] = np.median(wave_z[p,:]) 
	mean_flux_waves[p] = np.mean(wave_z[p,:]) 
	ax1.plot([0.0, 4.0], [median_values[p], median_values[p]], color = 'green', alpha = 0.5)
	ax1.plot([0.0, 4.0], [mean_values[p], mean_values[p]], color = 'red', alpha = 0.5)
	#plt.scatter(median_flux_waves, median_values, color = 'black', label = 'Median Values')
	#plt.scatter(mean_flux_waves, mean_values, color = 'red', label = 'Mean Values')
	plt.xlabel('Wavelength')
	maxx = np.max(wave_z[p,:]) + ((np.max(wave_z[p,:]) - 0.065)*0.1)

	if (args.logx):
		plt.semilogx()
	if (args.ratio):
		plt.plot([0.0, 4.0], [1.0, 1.0], color = 'black')
		plt.ylabel('template_photometry / true_photometry')
		plt.axis([0.065, maxx, 0.0, 2.0])
	else:
		plt.plot([0.0, 4.0], [0.0, 0.0], color = 'black')
		plt.ylabel('(template_photometry - true_photometry) / true_photometry')
		plt.axis([0.065, maxx, -2.0, 2.0])
	#print np.min(wave_z[p,:])-0.01
	#plt.axis([np.min(wave_z[p,:])-0.01, np.max(wave_z[p,:])+0.2, 0.5, 1.5])
	plt.legend()
	plt.title('Filter: '+noisy_jades_filters[filters_to_use_indices][p])
	plt.savefig(optional_output_folder+'Photometric_Offsets_Plot_'+noisy_jades_filters[filters_to_use_indices][p]+'.png', dpi=200)	

if (args.ratio):
	print "For template_photometry / true_photometry:"
else:
	print "For (template_photometry - true_photometry) / true_photometry:"
print "Median: "
print median_values
print "Mean: "
print mean_values
	
