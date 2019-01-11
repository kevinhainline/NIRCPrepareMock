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

# Which version of the mock are you using
JAGUAR_version = 'r1_v1.2'
JAGUAR_path = '/Users/knh/Desktop/NIRCam/mock_catalog_for_public/'

def FluxtoABMag(flux):
	return (-5.0 / 2.0) * np.log10(flux) - 48.60

region_width = 12.0#11.0 # In arcminutes
region_height = 4.84#11.0 # In arcminutes
region_area = region_width * region_height # In square arcminutes

######################
# Required Arguments #
######################

parser = argparse.ArgumentParser()

# Input file
parser.add_argument(
  '-in','--inputfile',
  help="Input file with noisy fluxes",
  action="store",
  type=str,
  dest="input_file",
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

# Filter 3
parser.add_argument(
  '-f3','--filter3',
  help="Filter 3",
  action="store",
  type=str,
  dest="filter3",
  required=True
)

# Color Limit 
parser.add_argument(
  '-clim','--climit',
  help="Color Limit? ",
  action="store",
  type=float,
  dest="climit",
  required=True
)

######################
# Optional Arguments #
######################

# SNR Limit? 
parser.add_argument(
  '-snr','--snr',
  help="Filter SNR (default = 0.0)",
  action="store",
  type=float,
  dest="snr",
  required=False
)

# Magnitude Limit?
parser.add_argument(
  '-maglim','--maglim',
  help="Magnitude Limit (will supersede SNR limits)",
  action="store",
  type=float,
  dest="maglim",
  required=False
)

# Blue Confirmation Filter
parser.add_argument(
  '-bf','--bfilter',
  help="X-Axis Filter 1",
  action="store",
  type=str,
  dest="bfilter",
  required=False
)

# blue SNR Limit? 
parser.add_argument(
  '-bsnr','--bsnr',
  help="Blue Filter SNR (default = 0.0)",
  action="store",
  type=float,
  dest="bsnr",
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

# Redshift Limit 
parser.add_argument(
  '-zlim','--zlimit',
  help="Redshift Limit? ",
  action="store",
  type=float,
  dest="zlimit",
  required=True
)



# Plot to Screen? 
parser.add_argument(
  '-ps','--plot_to_screen',
  help="Display Plot on Screen (Don't Save)?",
  action="store_true",
  dest="plot_to_screen",
  required=False
)

# Verbose? 
parser.add_argument(
  '-verb','--verbose',
  help="Verbose Mode?",
  action="store_true",
  dest="verbose",
  required=False
)

# Print Output Objects? 
parser.add_argument(
  '-prout','--proutput',
  help="Print the output objects?",
  action="store_true",
  dest="proutput",
  required=False
)

# Make Histogram of Galaxy Parameter
parser.add_argument(
  '-hp','--hparam',
  help="Galaxy Parameter for Histogram?",
  action="store",
  type=str,
  dest="histparam",
  required=False
)

# Upper limit on the filter2 - filter3 color
parser.add_argument(
  '-xclim','--xclimit',
  help="Upper limit on filter2 - filter3 color?",
  action="store",
  type=float,
  dest="xclimit",
  required=False
)

args=parser.parse_args()

# OPEN THE INPUT DATA FILE
if (args.input_file.endswith('.fits')):
	fitsinput = fits.open(args.input_file)
	ID_values = fitsinput[1].data['ID']
	redshifts = fitsinput[1].data['redshift']
	number_objects = ID_values.size

	number_input_objects = len(ID_values)

# Get the x and y filters.
filter_flux_one = fitsinput[1].data[args.filter1]
filter_flux_one_err = fitsinput[1].data[args.filter1+'_err']
filter_flux_two = fitsinput[1].data[args.filter2]
filter_flux_two_err = fitsinput[1].data[args.filter2+'_err']
filter_flux_three = fitsinput[1].data[args.filter3]
filter_flux_three_err = fitsinput[1].data[args.filter3+'_err']

if (args.bfilter):
	bfilter_flux_one = fitsinput[1].data[args.bfilter]
	bfilter_flux_one_err = fitsinput[1].data[args.bfilter+'_err']

# Have to do a SNR check. 
filter_flux_one_snr = filter_flux_one / filter_flux_one_err
filter_flux_two_snr = filter_flux_two / filter_flux_two_err
filter_flux_three_snr = filter_flux_three / filter_flux_three_err
if (args.bfilter):
	bfilter_flux_snr = bfilter_flux_one / bfilter_flux_one_err

if (args.snr):
	snr_limit = args.snr
else:
	snr_limit = 0.0

if (args.bsnr):
	bsnr_limit = args.bsnr
else:
	bsnr_limit = 0.0


if (args.maglim):
	filter_two_mag = FluxtoABMag(fitsinput[1].data[args.filter2]*1e-23*1e-9)
	filter_three_mag = FluxtoABMag(fitsinput[1].data[args.filter3]*1e-23*1e-9)
	SNR_check_objects = np.where((filter_two_mag < args.maglim) & (filter_three_mag < args.maglim))[0]
else:
	# Make sure that you have a detection in two filters red of the Lyman break. 
	SNR_check_objects = np.where((filter_flux_two_snr > snr_limit) & (filter_flux_three_snr > snr_limit))[0]

# Now, select the y_filter magnitudes above the SNR check
filter_one_mag_SNRcheck = FluxtoABMag(fitsinput[1].data[args.filter1][SNR_check_objects]*1e-23*1e-9)
filter_one_mag_SNRcheck[np.where(filter_flux_one[SNR_check_objects] < 0)[0]] = 99.00
filter_two_mag_SNRcheck = FluxtoABMag(fitsinput[1].data[args.filter2][SNR_check_objects]*1e-23*1e-9)
filter_three_mag_SNRcheck = FluxtoABMag(fitsinput[1].data[args.filter3][SNR_check_objects]*1e-23*1e-9)

# Get the redshifts above the SNR limit
redshifts_SNRcheck = redshifts[SNR_check_objects]

# And now the x and y axes. 
color_SNRcheck = filter_one_mag_SNRcheck - filter_two_mag_SNRcheck
xcolor_SNRcheck = filter_two_mag_SNRcheck - filter_three_mag_SNRcheck

# Now check whether or not objects are above the color limit(s). 
if (args.xclimit):
	objects_above_cut = np.where((color_SNRcheck > args.climit) & (xcolor_SNRcheck < args.xclimit))[0]
else:
	objects_above_cut = np.where(color_SNRcheck > args.climit)[0]

maglimprint = 0
snrlimprint = 0
xlimitprint = 0
#print "   "
if (args.maglim):
	maglimprint = 1
	#print "For "+args.filter2+" and "+args.filter3+" mag < "+str(args.maglim)
else:
	snrlimprint = 1
	#print "For "+args.filter2+" and "+args.filter3+" detections at SNR > "+str(snr_limit)
#print "And a cut of "+args.filter1+" - "+args.filter2+" > "+str(args.climit)
if (args.xclimit):
	xlimitprint = 1
	#print "And a cut of "+args.filter2+" - "+args.filter3+" < "+str(args.xclimit)

if (args.zlimit):
	objects_with_redshift_above_cut = np.where(redshifts_SNRcheck[objects_above_cut] >= args.zlimit)[0]
	objects_with_redshift_below_cut = np.where(redshifts_SNRcheck[objects_above_cut] < args.zlimit)[0]
	number_above_cut = len(objects_with_redshift_above_cut)/region_area
	number_below_cut = len(objects_with_redshift_below_cut)/region_area
	#print "There are "+str(round(number_above_cut,3))+" objects/square arcminute selected with z > "+str(args.zlimit)
	#print "There are "+str(round(number_below_cut,3))+" objects/square arcminute selected with z < "+str(args.zlimit)

if ((maglimprint == 1) & (xlimitprint == 0)):
	print "# maglim ycolorcut zcut highzdensity lowzdensity"
	print args.maglim, args.climit, args.zlimit, str(round(number_above_cut,3)), str(round(number_below_cut,3))
if ((maglimprint == 1) & (xlimitprint == 1)):
	print "# maglim ycolorcut xcolorcut zcut highzdensity lowzdensity"
	print args.maglim, args.climit, args.xclimit, args.zlimit, str(round(number_above_cut,3)), str(round(number_below_cut,3))
if ((snrlimprint == 1) & (xlimitprint == 0)):
	print "# SNRlim ycolorcut zcut highzdensity lowzdensity"
	print snr_limit, args.climit, args.zlimit, str(round(number_above_cut,3)), str(round(number_below_cut,3))
if ((snrlimprint == 1) & (xlimitprint == 1)):
	print "# SNRlim ycolorcut xcolorcut zcut highzdensity lowzdensity"
	print snr_limit, args.climit, args.xclimit, args.zlimit, str(round(number_above_cut,3)), str(round(number_below_cut,3))

# If you have a blueward filter, make sure that there's a nondetection
if (args.bfilter):
	if (args.maglim):
		SNR_check_objects_blue = np.where((filter_two_mag_snr < args.maglim) & (filter_three_mag < args.maglim) & (bfilter_flux_snr < bsnr_limit))[0]
	else:
		SNR_check_objects_blue = np.where((filter_flux_two_snr > snr_limit) & (filter_flux_three_snr > snr_limit) & (bfilter_flux_snr < bsnr_limit))[0]

	# Now, select the y_filter magnitudes above the SNR check
	filter_one_mag_SNRcheck_blue = FluxtoABMag(fitsinput[1].data[args.filter1][SNR_check_objects_blue]*1e-23*1e-9)
	filter_one_mag_SNRcheck_blue[np.where(filter_flux_one[SNR_check_objects_blue] < 0)[0]] = 99.00
	filter_two_mag_SNRcheck_blue = FluxtoABMag(fitsinput[1].data[args.filter2][SNR_check_objects_blue]*1e-23*1e-9)
	filter_three_mag_SNRcheck_blue = FluxtoABMag(fitsinput[1].data[args.filter3][SNR_check_objects_blue]*1e-23*1e-9)

	# Get the redshifts above the SNR limit
	redshifts_SNRcheck_blue = redshifts[SNR_check_objects_blue]
	
	# And now the x and y axes. 
	color_SNRcheck_blue = filter_one_mag_SNRcheck_blue - filter_two_mag_SNRcheck_blue
	xcolor_SNRcheck_blue = filter_two_mag_SNRcheck_blue - filter_three_mag_SNRcheck_blue
	
	# Now check whether or not objects are above the color limit(s). 
	if (args.xclimit):
		objects_above_cut_blue = np.where((color_SNRcheck_blue > args.climit) & (xcolor_SNRcheck_blue < args.xclimit))[0]
	else:
		objects_above_cut_blue = np.where(color_SNRcheck_blue > args.climit)[0]

	if (args.zlimit):
		objects_with_redshift_above_cut_blue = np.where(redshifts_SNRcheck_blue[objects_above_cut_blue] >= args.zlimit)[0]
		objects_with_redshift_below_cut_blue = np.where(redshifts_SNRcheck_blue[objects_above_cut_blue] < args.zlimit)[0]
		number_above_cut_blue = len(objects_with_redshift_above_cut_blue)/region_area
		number_below_cut_blue = len(objects_with_redshift_below_cut_blue)/region_area
		#print "    ...and only "+str(round(number_above_cut_blue,3))+" objects/square arcminute selected with z > "+str(args.zlimit)+" (with "+args.bfilter+" SNR < "+str(bsnr_limit)+")"
		#print "    ...and only "+str(round(number_below_cut_blue,3))+" objects/square arcminute selected with z < "+str(args.zlimit)+" (with "+args.bfilter+" SNR < "+str(bsnr_limit)+")"

redshift_deltaz = 0.2

if (args.histparam):

	# Let's open up the full SF catalog
	full_sf_mock_file = JAGUAR_path+'JADES_SF_mock_'+JAGUAR_version+'.fits'
	sftestfits = fits.open(full_sf_mock_file)
	
	# Get redshifts and galaxy properties
	catalog_sf_IDs = sftestfits[1].data['ID']
	catalog_sf_redshifts = sftestfits[1].data['redshift']
	catalog_sf_parameter = sftestfits[1].data[args.histparam]
	n_sf_objects = len(catalog_sf_IDs)
	
	# Let's open up the full Q catalog
	full_q_mock_file = JAGUAR_path+'JADES_Q_mock_'+JAGUAR_version+'.fits'
	qtestfits = fits.open(full_q_mock_file)
	
	# Get redshifts and galaxy properties
	catalog_q_IDs = qtestfits[1].data['ID']
	catalog_q_redshifts = qtestfits[1].data['redshift']
	catalog_q_parameter = qtestfits[1].data[args.histparam]
	n_q_objects = len(catalog_q_IDs)
	
	catalog_all_IDs = np.append(catalog_sf_IDs, catalog_q_IDs)
	catalog_all_redshifts = np.append(catalog_sf_redshifts, catalog_q_redshifts)
	catalog_all_param = np.append(catalog_sf_parameter, catalog_q_parameter)
	n_all_objects = len(catalog_all_IDs)
	
	# Get the indices for all of the possible objects
	# Objects above the color cut
	objects_above_cut_indices = np.zeros(len(objects_above_cut), dtype = int)
	for x in range(0, len(objects_above_cut)):
		objects_above_cut_indices[x] = np.where(catalog_all_IDs == ID_values[SNR_check_objects][objects_above_cut][x])[0][0]

	# Objects above the color cut, with redshifts *above* the redshift cut
	objects_with_redshift_above_cut_indices = np.zeros(len(objects_with_redshift_above_cut), dtype = int)
	for x in range(0, len(objects_with_redshift_above_cut)):
		objects_with_redshift_above_cut_indices[x] = np.where(catalog_all_IDs == ID_values[SNR_check_objects][objects_above_cut][objects_with_redshift_above_cut][x])[0][0]

	# Objects above the color cut, with redshifts *below* the redshift cut
	objects_with_redshift_below_cut_indices = np.zeros(len(objects_with_redshift_below_cut), dtype = int)
	for x in range(0, len(objects_with_redshift_below_cut)):
		objects_with_redshift_below_cut_indices[x] = np.where(catalog_all_IDs == ID_values[SNR_check_objects][objects_above_cut][objects_with_redshift_below_cut][x])[0][0]

	# And for those objects with the blue filter 
	if (args.bfilter):

		# Objects above the color cut with blue filter
		objects_above_cut_blue_indices = np.zeros(len(objects_above_cut_blue), dtype = int)
		for x in range(0, len(objects_above_cut_blue)):
			objects_above_cut_blue_indices[x] = np.where(catalog_all_IDs == ID_values[SNR_check_objects_blue][objects_above_cut_blue][x])[0][0]
	
		# Objects above the color cut, with redshifts *above* the redshift cut with blue filter
		objects_with_redshift_above_cut_blue_indices = np.zeros(len(objects_with_redshift_above_cut_blue), dtype = int)
		for x in range(0, len(objects_with_redshift_above_cut_blue)):
			objects_with_redshift_above_cut_blue_indices[x] = np.where(catalog_all_IDs == ID_values[SNR_check_objects_blue][objects_above_cut_blue][objects_with_redshift_above_cut_blue][x])[0][0]
	
		# Objects above the color cut, with redshifts *below* the redshift cut with blue filter
		objects_with_redshift_below_cut_blue_indices = np.zeros(len(objects_with_redshift_below_cut_blue), dtype = int)
		for x in range(0, len(objects_with_redshift_below_cut_blue)):
			objects_with_redshift_below_cut_blue_indices[x] = np.where(catalog_all_IDs == ID_values[SNR_check_objects_blue][objects_above_cut_blue][objects_with_redshift_below_cut_blue][x])[0][0]
	
	#print catalog_all_IDs[objects_with_redshift_above_cut_blue_indices]

# Make the two plots. 
#plt.clf()
if (args.histparam):
	plt.figure(figsize=(8, 11))
	plt.subplot(211)
else:
	plt.figure(figsize=(8,6))
redshift_array = np.arange(0,15.0+redshift_deltaz,redshift_deltaz)
redshift_hist = np.histogram(redshifts_SNRcheck[objects_above_cut], bins = redshift_array)
#plt.scatter(redshift_array[0:-1], redshift_hist[0] / region_area, color = 'red', label = 'All objects with '+args.filter1+' - '+args.filter2+' > '+str(args.climit))
plt.bar(redshift_array[0:-1], redshift_hist[0] / region_area, width = redshift_deltaz, color = 'red', label = 'All objects with '+args.filter1+' - '+args.filter2+' > '+str(args.climit))
if (args.bfilter):
	redshift_hist_blue = np.histogram(redshifts_SNRcheck_blue[objects_above_cut_blue], bins = redshift_array)
	#plt.scatter(redshift_array[0:-1], redshift_hist_blue[0] / region_area, color = 'blue', label = '...and '+args.bfilter+' SNR < '+str(bsnr_limit))
	plt.bar(redshift_array[0:-1], redshift_hist_blue[0] / region_area, width = redshift_deltaz, color = 'blue', label = '...and '+args.bfilter+' SNR < '+str(bsnr_limit))

plt.xlabel('Spectroscopic Redshift')
plt.ylabel('N/Area (square arcmin) in bins of delta_z = '+str(redshift_deltaz))
if (args.maglim):
	plt.title('All objects with '+args.filter2+' and '+args.filter3+' mag < '+str(args.maglim))
else:
	plt.title('All objects with '+args.filter2+' and '+args.filter3+' SNR > '+str(snr_limit))
plt.legend()


if (args.histparam):
	# And now for the second plot, the mass histogram. 
	plt.subplot(212)
	redshift_deltaparam = 0.1
	parameter_array = np.arange(6,13,redshift_deltaparam)
	# All objects
	parameter_hist = np.histogram(catalog_all_param[objects_above_cut_indices], bins = parameter_array)
	plt.bar(parameter_array[0:-1], parameter_hist[0] / region_area, width = redshift_deltaparam, color = 'red', label = 'All objects with '+args.filter1+' - '+args.filter2+' > '+str(args.climit))
	# Objects Above the Redshift Cut
	parameter_hist = np.histogram(catalog_all_param[objects_with_redshift_above_cut_indices], bins = parameter_array)
	plt.bar(parameter_array[0:-1], parameter_hist[0] / region_area, width = redshift_deltaparam, color = 'darkred', label = 'z > '+str(args.zlimit))
	# Objects Below the Redshift Cut
	parameter_hist = np.histogram(catalog_all_param[objects_with_redshift_below_cut_indices], bins = parameter_array)
	plt.bar(parameter_array[0:-1], parameter_hist[0] / region_area, width = redshift_deltaparam, color = 'lightsalmon', label = 'z < '+str(args.zlimit))

	if (args.bfilter):
		# All objects with blue filter
		parameter_hist = np.histogram(catalog_all_param[objects_above_cut_blue_indices], bins = parameter_array)
		plt.bar(parameter_array[0:-1], parameter_hist[0] / region_area, width = redshift_deltaparam, color = 'blue', label = 'All objects with '+args.filter1+' - '+args.filter2+' > '+str(args.climit))
		# Objects Above the Redshift Cut with blue filter
		parameter_hist = np.histogram(catalog_all_param[objects_with_redshift_above_cut_blue_indices], bins = parameter_array)
		plt.bar(parameter_array[0:-1], parameter_hist[0] / region_area, width = redshift_deltaparam, color = 'midnightblue', label = 'z > '+str(args.zlimit))
		# Objects Below the Redshift Cut with blue filter
		parameter_hist = np.histogram(catalog_all_param[objects_with_redshift_below_cut_blue_indices], bins = parameter_array)
		plt.bar(parameter_array[0:-1], parameter_hist[0] / region_area, width = redshift_deltaparam, color = 'skyblue', label = 'z < '+str(args.zlimit))

	plt.xlabel('log(Mass)')
	plt.ylabel('N/Area (square arcmin) in bins of delta_'+args.histparam+' = '+str(redshift_deltaparam))
	plt.legend()
plt.show()

if (args.proutput):
	print "    "
	print "Here are the objects above the color cut: "+args.filter1+" - "+args.filter2+" > "+str(args.climit)
	spacer = "  "
	print "ID  spec_z  "+args.filter1+"_fnu "+args.filter1+"_fnu_err "+args.filter2+"_fnu "+args.filter2+"_fnu_err "+args.filter3+"_fnu "+args.filter3+"_fnu_err "
	for q in range(0,len(objects_above_cut_indices)):
		object_ID = np.where(ID_values == catalog_all_IDs[objects_above_cut_indices[q]])[0][0]
		print catalog_all_IDs[objects_above_cut_indices[q]], redshifts[object_ID], spacer, np.round(filter_flux_one[object_ID],4), np.round(filter_flux_one_err[object_ID],4), spacer, np.round(filter_flux_two[object_ID],4), np.round(filter_flux_two_err[object_ID],4), spacer, np.round(filter_flux_three[object_ID],4), np.round(filter_flux_three_err[object_ID],4) 

	if (args.bfilter):
		print "    "
		print "And here are the objects above the color cut, with "+args.bfilter+" SNR < "+str(bsnr_limit)
		print "ID  spec_z  "+args.bfilter+"_fnu "+args.bfilter+"_fnu_err "+args.filter1+"_fnu "+args.filter1+"_fnu_err "+args.filter2+"_fnu "+args.filter2+"_fnu_err "+args.filter3+"_fnu "+args.filter3+"_fnu_err "
		for q in range(0,len(objects_above_cut_blue_indices)):
			object_ID = np.where(ID_values == catalog_all_IDs[objects_above_cut_blue_indices[q]])[0][0]
			print catalog_all_IDs[objects_above_cut_blue_indices[q]], redshifts[object_ID], spacer, np.round(bfilter_flux_one[object_ID],4), np.round(bfilter_flux_one_err[object_ID],4), spacer, np.round(filter_flux_one[object_ID],4), np.round(filter_flux_one_err[object_ID],4), spacer, np.round(filter_flux_two[object_ID],4), np.round(filter_flux_two_err[object_ID],4), spacer, np.round(filter_flux_three[object_ID],4), np.round(filter_flux_three_err[object_ID],4) 
		
		

