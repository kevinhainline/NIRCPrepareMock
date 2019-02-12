import os
import sys
import math
import astropy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
from astropy.table import Table


def FluxtoABMag(flux):
	return (-5.0 / 2.0) * np.log10(flux) - 48.60


# Define the cuts. 
lineone_isvertical = 0
lineslope_1 = 0
#lineyintercept_1 = 0.75
lineyintercept_1 = 1.3
xvalue_1 = -9999

linetwo_isvertical = 0
lineslope_2 = 0.5
#lineyintercept_2 = 0.85
lineyintercept_2 = 1.3
xvalue_2 = -9999

linethree_isvertical = 1
lineslope_3 = -9999
lineyintercept_3 = -9999
xvalue_3 = 0.2


# Define your color cut selection in terms of the x and y colors of the first plot
# and the second plot. 
def ColorCuts(colorx, colory):
	if ((lineone_isvertical == 0) & (linetwo_isvertical == 0) & (linethree_isvertical == 0)):
		return np.where(
			(colory > ((lineslope_1 * colorx) + lineyintercept_1))
			& (colory > ((lineslope_2 * colorx) + lineyintercept_2))
			& (colory > ((lineslope_3 * colorx) + lineyintercept_3))
			)[0]
	if ((lineone_isvertical == 0) & (linetwo_isvertical == 0) & (linethree_isvertical == 1)):
		return np.where(
			(colory > ((lineslope_1 * colorx) + lineyintercept_1))
			& (colory > ((lineslope_2 * colorx) + lineyintercept_2))
			& (colorx < xvalue_3)
			)[0]	
	if ((lineone_isvertical == 1) & (linetwo_isvertical == 0) & (linethree_isvertical == 1)):
		return np.where(
			(colorx > xvalue_1)
			& (colory > ((lineslope_2 * colorx) + lineyintercept_2))
			& (colorx < xvalue_3)
			)[0]	

# If you want to plot the selection function on the color-color plots, you can
# specify an array of X and Y values connecting the selection function:
if ((lineone_isvertical == 0) & (linetwo_isvertical == 0) & (linethree_isvertical == 0)):
	point_1_x = -10.0
	point_1_y = lineslope_1 * point_1_x + lineyintercept_1

	point_2_x = (lineyintercept_2 - lineyintercept_1) / (lineslope_1 - lineslope_2)
	point_2_y = lineslope_2 * point_2_x + lineyintercept_2

	point_3_x = (lineyintercept_3 - lineyintercept_2) / (lineslope_2 - lineslope_3)
	point_3_y = lineslope_2 * point_3_x + lineyintercept_2

	point_4_x = 10.0
	point_4_y = lineslope_3 * point_4_x + lineyintercept_3

	selection_points_plot1 = np.array([[point_1_x, point_1_y], [point_2_x, point_2_y], [point_3_x, point_3_y], [point_4_x, point_4_y]])

if ((lineone_isvertical == 0) & (linetwo_isvertical == 0) & (linethree_isvertical == 1)):
	point_1_x = -10.0
	point_1_y = lineslope_1 * -10.0 + lineyintercept_1

	point_2_x = (lineyintercept_2 - lineyintercept_1) / (lineslope_1 - lineslope_2)
	point_2_y = lineslope_2 * point_2_x + lineyintercept_2

	point_3_x = xvalue_3
	point_3_y = lineslope_2 * point_3_x + lineyintercept_2

	point_4_x = xvalue_3
	point_4_y = 10.0

	selection_points_plot1 = np.array([[point_1_x, point_1_y], [point_2_x, point_2_y], [point_3_x, point_3_y], [point_4_x, point_4_y]])

if ((lineone_isvertical == 1) & (linetwo_isvertical == 0) & (linethree_isvertical == 1)):
	point_1_x = xvalue_1
	point_1_y = 10.0

	point_2_x = xvalue_1
	point_2_y = lineslope_2 * xvalue_1 + lineyintercept_2

	point_3_x = xvalue_3
	point_3_y = lineslope_2 * xvalue_3 + lineyintercept_2

	point_4_x = xvalue_3
	point_4_y = 10.0

	selection_points_plot1 = np.array([[point_1_x, point_1_y], [point_2_x, point_2_y], [point_3_x, point_3_y], [point_4_x, point_4_y]])

# Which version of the mock are you using
JAGUAR_version = 'r1_v1.2'

# If you want to create an output plot, here's where you specify the title
output_filename = 'Color_Color_vs_Line_EW_HighZ.png'


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

# Filter 3
parser.add_argument(
  '-f3','--filter3',
  help="Filter 3",
  action="store",
  type=str,
  dest="filter3",
  required=True
)

# Filter 4
parser.add_argument(
  '-f4','--filter4',
  help="Filter 4",
  action="store",
  type=str,
  dest="filter4",
  required=True
)


######################
# Optional Arguments #
######################

# Comp Parameter
parser.add_argument(
  '-comp_param','--comp_param',
  help="Comparison Parameter?",
  action="store",
  type=str,
  dest="comp_param",
  required=False
)

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

# Blue Filter
parser.add_argument(
  '-bf','--bluefilter',
  help="Blue Rejection Filter?",
  action="store",
  type=str,
  dest="bluefilter",
  required=False
)

# SNR Limit 
parser.add_argument(
  '-snr','--snr',
  help="Filter 1, 2, 3, 4, 5 AND 6 SNR (default = 0.0)",
  action="store",
  type=float,
  dest="snrfilter",
  required=False
)

# SNR Limit, Filter One? 
parser.add_argument(
  '-snrf1','--snrf1',
  help="Filter 1 SNR (default = 0.0)",
  action="store",
  type=float,
  dest="snrfilter1",
  required=False
)

# SNR Limit, Filter Two? 
parser.add_argument(
  '-snrf2','--snrf2',
  help="Filter 2 SNR (default = 0.0)",
  action="store",
  type=float,
  dest="snrfilter2",
  required=False
)

# SNR Limit, Filter Three? 
parser.add_argument(
  '-snrf3','--snrf3',
  help="Filter 3 SNR (default = 0.0)",
  action="store",
  type=float,
  dest="snrfilter3",
  required=False
)

# SNR Limit, Filter Four? 
parser.add_argument(
  '-snrf4','--snrf4',
  help="Filter 4 SNR (default = 0.0)",
  action="store",
  type=float,
  dest="snrfilter4",
  required=False
)

# SNR Limit, Blue Filter? 
parser.add_argument(
  '-snrbf','--snrbf',
  help="Blue Rejection Filter SNR (default = 0.0)",
  action="store",
  type=float,
  dest="snrbluefilter",
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

# Redshift Minimum 
parser.add_argument(
  '-clean','--clean',
  help="Plot Clean Points instead of Noisy Points? ",
  action="store_true",
  dest="clean",
  required=False
)

# Set up the parsing from above.
args=parser.parse_args()

# Set the JAGUAR Path
JAGUAR_Path = args.input_folder

# Set the filters to explore
filter1 = args.filter1
filter2 = args.filter2
filter3 = args.filter3
filter4 = args.filter4


if (args.snrfilter):
	snr_limit_one = args.snrfilter
	snr_limit_two = args.snrfilter
	snr_limit_three = args.snrfilter
	snr_limit_four = args.snrfilter

else:
	snr_limit_one = 0.0
	snr_limit_two = 0.0
	snr_limit_three = 0.0
	snr_limit_four = 0.0


# Set the SNR limits on the filters
if (args.snrfilter1):
	snr_limit_one = args.snrfilter1
	
if (args.snrfilter2):
	snr_limit_two = args.snrfilter2
	
if (args.snrfilter3):
	snr_limit_three = args.snrfilter3

if (args.snrfilter4):
	snr_limit_four = args.snrfilter4

	
if (args.clean is False):

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
	filter_flux_three = fitsinput[1].data[filter3]
	filter_flux_three_err = fitsinput[1].data[filter3+'_err']
	filter_flux_four = fitsinput[1].data[filter4]
	filter_flux_four_err = fitsinput[1].data[filter4+'_err']
	
	# Let's do a SNR check. 
	filter_flux_one_snr = filter_flux_one / filter_flux_one_err
	filter_flux_two_snr = filter_flux_two / filter_flux_two_err
	filter_flux_three_snr = filter_flux_three / filter_flux_three_err
	filter_flux_four_snr = filter_flux_four / filter_flux_four_err

	SNR_check_objects = np.where((filter_flux_one_snr > snr_limit_one) & (filter_flux_two_snr > snr_limit_two) 
		& (filter_flux_three_snr > snr_limit_three) & (filter_flux_four_snr > snr_limit_four))[0]

	filter_one = FluxtoABMag(filter_flux_one[SNR_check_objects]*1e-23*1e-9)
	filter_one[np.where(filter_flux_one[SNR_check_objects] < 0)[0]] = 99.00
	filter_two = FluxtoABMag(filter_flux_two[SNR_check_objects]*1e-23*1e-9)
	filter_two[np.where(filter_flux_two[SNR_check_objects] < 0)[0]] = 99.00
	filter_three = FluxtoABMag(filter_flux_three[SNR_check_objects]*1e-23*1e-9)
	filter_three[np.where(filter_flux_three[SNR_check_objects] < 0)[0]] = 99.00
	filter_four = FluxtoABMag(filter_flux_four[SNR_check_objects]*1e-23*1e-9)
	filter_four[np.where(filter_flux_four[SNR_check_objects] < 0)[0]] = 99.00
	
	ID_values_SNR = ID_values[SNR_check_objects]	
	
	filter_color_y = filter_one - filter_two
	filter_color_x = filter_three - filter_four


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
if (args.comp_param):
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
sf_clean_filter_flux_three = sftestfits[1].data[filter3+'_fnu']
sf_clean_filter_flux_four = sftestfits[1].data[filter4+'_fnu']
sf_clean_filter_one = FluxtoABMag(sf_clean_filter_flux_one*1e-23*1e-9)
sf_clean_filter_one[np.where(sf_clean_filter_flux_one < 0)[0]] = 99.00
sf_clean_filter_two = FluxtoABMag(sf_clean_filter_flux_two*1e-23*1e-9)
sf_clean_filter_two[np.where(sf_clean_filter_flux_two < 0)[0]] = 99.00
sf_clean_filter_three = FluxtoABMag(sf_clean_filter_flux_three*1e-23*1e-9)
sf_clean_filter_three[np.where(sf_clean_filter_flux_three < 0)[0]] = 99.00
sf_clean_filter_four = FluxtoABMag(sf_clean_filter_flux_four*1e-23*1e-9)
sf_clean_filter_four[np.where(sf_clean_filter_flux_four < 0)[0]] = 99.00


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
if (args.comp_param):
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
q_clean_filter_flux_three = qtestfits[1].data[filter3+'_fnu']
q_clean_filter_flux_four = qtestfits[1].data[filter4+'_fnu']
q_clean_filter_one = FluxtoABMag(q_clean_filter_flux_one*1e-23*1e-9)
q_clean_filter_one[np.where(q_clean_filter_flux_one < 0)[0]] = 99.00
q_clean_filter_two = FluxtoABMag(q_clean_filter_flux_two*1e-23*1e-9)
q_clean_filter_two[np.where(q_clean_filter_flux_two < 0)[0]] = 99.00
q_clean_filter_three = FluxtoABMag(q_clean_filter_flux_three*1e-23*1e-9)
q_clean_filter_three[np.where(q_clean_filter_flux_three < 0)[0]] = 99.00
q_clean_filter_four = FluxtoABMag(q_clean_filter_flux_four*1e-23*1e-9)
q_clean_filter_four[np.where(q_clean_filter_flux_four < 0)[0]] = 99.00


# Combine the IDs, redshifts, masses, etc, from the two catalogs.
catalog_all_IDs = np.append(catalog_sf_IDs, catalog_q_IDs)
catalog_all_redshifts = np.append(catalog_sf_redshifts, catalog_q_redshifts)
catalog_all_logmass = np.append(catalog_sf_logmass, catalog_q_logmass)
catalog_all_sfr_100 = np.append(catalog_sf_sfr_100, catalog_q_sfr_100)
if (args.comp_param):
	catalog_all_comparison_parameter = np.append(catalog_sf_comparison_parameter, catalog_q_comparison_parameter)
catalog_all_ssfr = np.append(catalog_sf_ssfr, catalog_q_ssfr)
catalog_all_max_stellar_age = np.append(catalog_sf_max_stellar_age, catalog_q_max_stellar_age)

all_clean_filter_one = np.append(sf_clean_filter_one, q_clean_filter_one)
all_clean_filter_two = np.append(sf_clean_filter_two, q_clean_filter_two)
all_clean_filter_three = np.append(sf_clean_filter_three, q_clean_filter_three)
all_clean_filter_four = np.append(sf_clean_filter_four, q_clean_filter_four)

n_all_objects = len(catalog_all_IDs)

# What is the clean catalog color (for comparison)
clean_filter_color_y = all_clean_filter_one - all_clean_filter_two
clean_filter_color_x = all_clean_filter_three - all_clean_filter_four

if (args.clean is False):
	# Match the noisy catalog to the full catalog
	all_object_indices = np.zeros(number_objects, dtype = int)
	for x in range(0, number_objects):
		all_object_indices[x] = np.where(catalog_all_IDs == ID_values[x])[0][0]
	
# Make the Color Selections
if (args.clean is True):
	clean_color_selection = ColorCuts(clean_filter_color_x, clean_filter_color_y)
if (args.clean is False):
	noisy_color_selection = ColorCuts(filter_color_x, filter_color_y)


fig = plt.figure(figsize=(15, 8))
plt.set_cmap('cubehelix')
ax1 = fig.add_subplot(121)
if (args.clean is True):
	cax = ax1.scatter(clean_filter_color_x, clean_filter_color_y, c = catalog_all_redshifts, s = 3, label = 'All Objects')
	ax1.scatter(clean_filter_color_x[clean_color_selection], clean_filter_color_y[clean_color_selection], color = 'red', s = 10, label = 'Satisfies Color Cuts')
if (args.clean is False):
	cax = ax1.scatter(filter_color_x, filter_color_y, c = redshifts[SNR_check_objects], s = 3, label = 'All Objects')
	ax1.scatter(filter_color_x[noisy_color_selection], filter_color_y[noisy_color_selection], color = 'red', s = 10, label = 'Satisfies Color Cuts')
cbar = fig.colorbar(cax, label = 'Redshift', orientation="horizontal")
# Plot the Selection Region #
if (len(selection_points_plot1) > 1):
	for q in range(0, len(selection_points_plot1)-1):
		ax1.plot([selection_points_plot1[q][0], selection_points_plot1[q+1][0]], [selection_points_plot1[q][1], selection_points_plot1[q+1][1]], color = 'blue')
# # # # # # # # # # # # # # #
ax1.set_xlabel(filter3+' - '+filter4)
ax1.set_ylabel(filter1+' - '+filter2)
ax1.axis([-1.0,1.0,-0.8,2.2])
ax1.legend()

if (args.comp_param):
	ax3 = plt.subplot(122)
	parameter_name = args.comp_param
	if (args.clean is True):
		ax3.scatter(catalog_all_redshifts, catalog_all_comparison_parameter, s=2, c = 'black', alpha = 0.1, label = 'All Objects')
		ax3.scatter(catalog_all_redshifts[clean_color_selection], catalog_all_comparison_parameter[clean_color_selection], s=8, c = 'red', label = 'Satisfies Color Cuts')
	if (args.clean is False):
		ax3.scatter(redshifts[SNR_check_objects], catalog_all_comparison_parameter[all_object_indices][SNR_check_objects], s=2, c = 'black', alpha = 0.1, label = 'All Objects')
		ax3.scatter(redshifts[SNR_check_objects][noisy_color_selection], catalog_all_comparison_parameter[all_object_indices][SNR_check_objects][noisy_color_selection], s=8, c = 'red', label = 'Satisfies Color Cuts')
	ax3.set_xlabel('Redshift')
	ax3.set_ylabel(parameter_name)
else:
	ax3 = plt.subplot(122)
	redshift_bins = np.arange(0,15,0.2)
	if (args.clean is True):
		ax3.hist(catalog_all_redshifts, bins = redshift_bins, color = 'black', alpha = 0.1, label = 'All Objects')
		ax3.hist(catalog_all_redshifts[clean_color_selection], bins = redshift_bins, color = 'red', label = 'Satisfies Color Cuts')
	if (args.clean is False):
		ax3.hist(redshifts[SNR_check_objects], bins = redshift_bins, color = 'black', alpha = 0.1, label = 'All Objects')
		ax3.hist(redshifts[SNR_check_objects][noisy_color_selection], bins = redshift_bins, color = 'red', label = 'Satisfies Color Cuts')
	ax3.set_xlabel('Redshift')
	ax3.set_ylabel('Number of Objects')
	

# Plot to screen or Make a File? 
if (args.plot_to_screen):
	plt.show()
else:
	plt.savefig(output_filename, dpi = 200)

