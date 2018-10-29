import os
import sys
import math
import argparse
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
import astropy.stats
import scipy.stats
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table

# Which version of the mock are you using
JAGUAR_version = 'r1_v1.1'

# Parser for figuring out the various input parameters
parser = argparse.ArgumentParser()

######################
# Required Arguments #
######################

# NIRCam SNR Limit
parser.add_argument(
  '-snrl','--snr_limit',
  help="SNR Limit for analysis?",
  action="store",
  type=float,
  dest="snr_limit",
  required=True
)

# Original JAGUAR Catalog Path
parser.add_argument(
  '-jaguar','--jaguar_path',
  help="Path to JAGUAR Catalogs?",
  action="store",
  type=str,
  dest="jaguar_path",
  required=True
)

# Galaxy Parameter for Subplot
parser.add_argument(
  '-gpar','--galaxy_param',
  help="Galaxy Parameter For Subplot",
  action="store",
  type=str,
  dest="galaxy_param",
  required=True
)

######################
# Optional Arguments #
######################

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

# Second Galaxy Parameter for Subplot
parser.add_argument(
  '-sgpar','--second_galaxy_param',
  help="Galaxy Parameter For Subplot",
  action="store",
  type=str,
  dest="second_galaxy_param",
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

args=parser.parse_args()


class SelectFromCollection(object):
    def __init__(self, ax, collection, alpha_other=0.2):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, self.Npts).reshape(self.Npts, -1)
        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero([path.contains_point(xy) for xy in self.xys])[0]
        # This is an older code that would reset the colors and the alpha values for the
        # selected objects
        #self.fc[:, -1] = self.alpha_other
        #self.fc[self.ind, -1] = 1
        #print self.fc[self.ind]
        #self.collection.set_facecolors(self.fc)
        # However currently, we plot the points that are selected with black dots.
        n_points = len(self.ind)
        ax.scatter(outlier_z_spec, outlier_z_phot, s=30, c = outlier_snr_colors, edgecolor = 'None')
        zspec, zphot = make_into_points(self.xys[self.ind])
        ax.scatter(zspec, zphot, s = 5, color = 'black')
		# And now we plot the inset. 
        plot_point_properties(self.xys[self.ind])
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

# This function takes in an array of points and turns them into two separate arrays. 
def make_into_points(point_array):
	number_points = len(point_array)
	xvalues = np.zeros(number_points)
	yvalues = np.zeros(number_points)
	for j in range(0, number_points):
		xvalues[j] = point_array[j][0]
		yvalues[j] = point_array[j][1]
	return xvalues, yvalues

# This function plots the inset, which is either a histogram of one property (as
# compared to all objects in the same redshift range) or a scatter plot comparing two
# different properties. 
def plot_point_properties(indxys):
	# First, we get the zspec and zphot values of the selected outliers.
	zspec, zphot = make_into_points(indxys)
	# And now we get the indices for these objects 
	outlier_indices = np.zeros(len(zspec), dtype = 'int')
	for j in range(0, len(zspec)):
		outlier_indices[j] = get_index_from_zs(zspec[j], zphot[j])
	# We print the IDs out to the screen
	print "\n"
	print "The Outlier ID's you selected"
	print "z_spec = "+str(np.min(zspec))+" - "+str(np.max(zspec))
	print "z_phot = "+str(np.min(zphot))+" - "+str(np.max(zphot))
	print catalog_all_IDs[outlier_indices]
	# And now we plot them as a subsample
	ax2 = fig.add_subplot(122)
	a = ax2.set_position([0.5, 0.2, .3, .3])
	ax2.clear()
	# This is if you have a second parameter and want to make a scatter plot
	if (args.second_galaxy_param):
		outlier_parameter_values = catalog_all_parameter[outlier_indices]
		outlier_second_parameter_values = second_catalog_all_parameter[outlier_indices]
		all_object_parameter_values = catalog_all_parameter[get_all_object_parameter_indices(np.min(zspec), np.max(zspec))]
		all_object_second_parameter_values = second_catalog_all_parameter[get_all_object_parameter_indices(np.min(zspec), np.max(zspec))]
		ax2.scatter(all_object_parameter_values, all_object_second_parameter_values, color = 'grey', alpha = 0.5)
		ax2.scatter(outlier_parameter_values, outlier_second_parameter_values, color = 'red')
		ax2.set_xlim([param_min, param_max])
		ax2.set_ylim([second_param_min, second_param_max])
		ax2.set_xlabel(param_label)
		ax2.set_ylabel(second_param_label)
	# This is if you have only one parameter and want to make a histogram
	else:
		outlier_parameter_values = catalog_all_parameter[outlier_indices]
		all_object_parameter_values = catalog_all_parameter[get_all_object_parameter_indices(np.min(zspec), np.max(zspec))]
		bin_delta = 0.2
		param_bins = np.arange(np.min(all_object_parameter_values), np.max(all_object_parameter_values), bin_delta)
		weights = np.ones_like(all_object_parameter_values)/float(len(all_object_parameter_values))
		ax2.hist(all_object_parameter_values, weights=weights, bins = param_bins, alpha = 0.5, color = 'grey')
		weights = np.ones_like(outlier_parameter_values)/float(len(outlier_parameter_values))
		ax2.hist(outlier_parameter_values, weights=weights, bins = param_bins, alpha = 1.0, color = 'red')
		ax2.set_xlim([param_min, param_max])
		ax2.set_xlabel(param_label)
	# And now we print out the number of galaxies in the subsample vs. the comparison sample
	ax2.set_title('SUBSAMPLE GALAXIES')
	ax2.text(0.75, 0.85, 'z = '+str(round(np.min(zspec),2))+' - '+str(round(np.max(zspec),2)), horizontalalignment='right', verticalalignment='top', color = 'black', transform=ax2.transAxes)
	ax2.text(0.75, 0.75, str(len(zspec))+' Galaxies', horizontalalignment='right', verticalalignment='top', color = 'red', transform=ax2.transAxes)
	ax2.text(0.75, 0.65, str(len(all_object_parameter_values))+' Galaxies', horizontalalignment='right', verticalalignment='top', color = 'grey', transform=ax2.transAxes)

# This function gets the indices in the full entire catalog from the photo-z's and spec-z's
def get_index_from_zs(zspec, zphot):
	outlier_group_ID_index = np.where((outlier_z_spec == zspec) & (outlier_z_phot == zphot))[0][0]
	outlier_group_ALL_index = np.where(catalog_all_IDs == outlier_z_IDs[outlier_group_ID_index])[0][0]
	return outlier_group_ALL_index

# This function gets the parameter indices for all of the plotted objects in a given
# redshift range.
def get_all_object_parameter_indices(minzspec, maxzspec):
#	all_object_indices = np.where((z_spec[high_SNR_high_prob_indices] >= minzspec) & (z_spec[high_SNR_high_prob_indices] <= maxzspec))[0]
	if (args.eazy_output_file):
		all_object_indices = np.where((z_spec[high_SNR_high_prob_lowq_z_indices] >= minzspec) & (z_spec[high_SNR_high_prob_lowq_z_indices] <= maxzspec))[0]
	if (args.bpz_output_file):
		all_object_indices = np.where((z_spec[high_SNR_high_prob_lowchisq2_indices] >= minzspec) & (z_spec[high_SNR_high_prob_lowchisq2_indices] <= maxzspec))[0]
	number_of_all_objects = len(all_object_indices)
	full_catalog_indices = np.zeros(number_of_all_objects, dtype = 'int')
	for b in range(0, number_of_all_objects):
		#full_catalog_indices[b] = np.where(catalog_all_IDs == ids[high_SNR_high_prob_indices][all_object_indices][b])[0][0]
		if (args.eazy_output_file):
			full_catalog_indices[b] = np.where(catalog_all_IDs == ids[high_SNR_high_prob_lowq_z_indices][all_object_indices][b])[0][0]
		if (args.bpz_output_file):
			full_catalog_indices[b] = np.where(catalog_all_IDs == ids[high_SNR_high_prob_lowchisq2_indices][all_object_indices][b])[0][0]
	return full_catalog_indices

# This function gets the catastrophic outliers from a group of objects given their z_spec
# and z_phot values, along with an outlier_fraction_value. 
def find_catastrophic_outliers(IDs, z_spec, z_phot, outlier_faction_value): 
	residual = abs((z_spec - z_phot) / (1 + z_spec))
	catastrophic_indices = np.where(residual > outlier_faction_value)[0]
	catastrophic_outlier_IDs = IDs[catastrophic_indices]
	catastrophic_outlier_z_spec = z_spec[catastrophic_indices]
	catastrophic_outlier_z_phot = z_phot[catastrophic_indices]
	return catastrophic_outlier_IDs, catastrophic_outlier_z_spec, catastrophic_outlier_z_phot


#filter_file_name = 'NoisySubsample_to_PhotoZInput_filters.dat'
#noisy_jades_filters = ['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 'HST_F850LP', 'NRC_F070W', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W']
#number_jades_filters = len(noisy_jades_filters)


# Set the minimum and maximum redshifts
if (args.minimum_z):
	min_redshift = args.minimum_z
else:
	min_redshift = 0
	
if (args.maximum_z):
	max_redshift = args.maximum_z
else:
	max_redshift = 15
	
# Setting the limit on the SNR 
SNR_limit = args.snr_limit
#log_SNR_limit = math.log10(SNR_limit)

# Let's open up the full SF catalog
full_sf_mock_file = args.jaguar_path+'JADES_SF_mock_'+JAGUAR_version+'.fits'
sftestfits = fits.open(full_sf_mock_file)
# Get redshifts and galaxy properties
catalog_sf_IDs = sftestfits[1].data['ID']
catalog_sf_redshifts = sftestfits[1].data['redshift']
n_sf_objects = len(catalog_sf_IDs)

# Let's open up the full Q catalog
full_q_mock_file = args.jaguar_path+'JADES_Q_mock_'+JAGUAR_version+'.fits'
qtestfits = fits.open(full_q_mock_file)
# Get redshifts and galaxy properties
catalog_q_IDs = qtestfits[1].data['ID']
catalog_q_redshifts = qtestfits[1].data['redshift']
n_q_objects = len(catalog_q_IDs)

# Set the specified galaxy parameters. 
if (args.galaxy_param == 'mStar'):
	catalog_sf_param = sftestfits[1].data[args.galaxy_param]
	catalog_q_param = qtestfits[1].data[args.galaxy_param]
	param_min = 5.5
	param_max = 13.0
	param_label = 'log(Mass)'
elif (args.galaxy_param == 'sSFR'):
	catalog_sf_param = sftestfits[1].data[args.galaxy_param]
	catalog_q_param = qtestfits[1].data[args.galaxy_param]
	param_min = -20
	param_max = 0.5
	param_label = 'sSFR'
elif (args.galaxy_param == 'tau'):
	catalog_sf_param = sftestfits[1].data[args.galaxy_param]
	catalog_q_param = qtestfits[1].data[args.galaxy_param]
	param_min = 0
	param_max = 5
	param_label = 'tau'
elif (args.galaxy_param == 'max_stellar_age'):
	catalog_sf_param = sftestfits[1].data[args.galaxy_param]
	catalog_q_param = qtestfits[1].data[args.galaxy_param]
	param_min = 6
	param_max = 10
	param_label = 'max_stellar_age'
elif (args.galaxy_param == 'tauV_eff'):
	catalog_sf_param = sftestfits[1].data[args.galaxy_param]
	catalog_q_param = (qtestfits[1].data['ID']*0)
	param_min = 0
	param_max = 3.0
	param_label = 'tauV_eff'

else:
	sys.exit("The Specified Galaxy Parameter ("+args.galaxy_parameter+") Is Currently Not Supported")

# Append the IDs, redshifts, and parameters between the SF and Q galaxies.
catalog_all_IDs = np.append(catalog_sf_IDs, catalog_q_IDs)
catalog_all_redshifts = np.append(catalog_sf_redshifts, catalog_q_redshifts)
catalog_all_parameter = np.append(catalog_sf_param, catalog_q_param)

# In case there's a second parameter, do the same for that. 
if (args.second_galaxy_param):
	if (args.second_galaxy_param == 'mStar'):
		second_catalog_sf_param = sftestfits[1].data[args.second_galaxy_param]
		second_catalog_q_param = qtestfits[1].data[args.second_galaxy_param]
		second_param_min = 5.5
		second_param_max = 13.0
		second_param_label = 'log(Mass)'
	elif (args.second_galaxy_param == 'sSFR'):
		second_catalog_sf_param = sftestfits[1].data[args.second_galaxy_param]
		second_catalog_q_param = qtestfits[1].data[args.second_galaxy_param]
		second_param_min = -20
		second_param_max = 0.5
		second_param_label = 'sSFR'
	elif (args.second_galaxy_param == 'tau'):
		second_catalog_sf_param = sftestfits[1].data[args.second_galaxy_param]
		second_catalog_q_param = qtestfits[1].data[args.second_galaxy_param]
		second_param_min = 0
		second_param_max = 5
		second_param_label = 'tau'
	elif (args.second_galaxy_param == 'max_stellar_age'):
		second_catalog_sf_param = sftestfits[1].data[args.second_galaxy_param]
		second_catalog_q_param = qtestfits[1].data[args.second_galaxy_param]
		second_param_min = 6
		second_param_max = 10
		second_param_label = 'max_stellar_age'
	elif (args.second_galaxy_param == 'tauV_eff'):
		second_catalog_sf_param = sftestfits[1].data[args.second_galaxy_param]
		second_catalog_q_param = (qtestfits[1].data['ID']*0)
		second_param_min = 0
		second_param_max = 3.0
		second_param_label = 'tauV_eff'
	else:
		sys.exit("The Specified Second Galaxy Parameter ("+args.second_galaxy_parameter+") Is Currently Not Supported")

	second_catalog_all_parameter = np.append(second_catalog_sf_param, second_catalog_q_param)


# Now, we open up the summary files. 
# EAZY
if (args.eazy_output_file):
	version = 'EAZY'
	prob_limit = 0.9
	q_z_limit = 3.0
	#ID  z_spec  z_phot  z_prob  q_z NRC_F200W_SNR 

	output_zs = np.loadtxt(args.eazy_output_file) 
	ids = output_zs[:,0]
	z_spec = output_zs[:,1]
	z_phot = output_zs[:,2]
	z_prob = output_zs[:,3]
	q_z = output_zs[:,4]
	NIRc_SNR = output_zs[:,5]
	title_for_plot = 'EAZY Results'
	logNIRc_SNR = NIRc_SNR
	logNIRc_SNR[np.where(NIRc_SNR <= 0)[0]] = 1e-5
	logNIRc_SNR = np.log10(logNIRc_SNR)
	
# Le Phare
if (args.lephare_output_file):
	version = 'LePhare'

	output_zs = np.loadtxt(args.lephare_output_file) 
	ids = output_zs[:,0]
	z_spec = output_zs[:,1]
	z_phot = output_zs[:,2]
	z_prob = output_zs[:,3]
	NIRc_SNR = output_zs[:,4]
	logNIRc_SNR = NIRc_SNR
	logNIRc_SNR[np.where(NIRc_SNR <= 0)[0]] = 1e-5
	logNIRc_SNR = np.log10(logNIRc_SNR)

	title_for_plot = 'Le Phare Results'

# BPZ
if (args.bpz_output_file):
	version = 'BPZ'
	prob_limit = 0.95
	chisq2_limit = 1.0
	#ID  z_spec  z_phot  ODDS  chisq2 NRC_F200W_SNR 

	output_zs = np.loadtxt(args.bpz_output_file) 
	ids = output_zs[:,0]
	z_spec = output_zs[:,1]
	z_phot = output_zs[:,2]
	z_prob = output_zs[:,3]
	chisq2 = output_zs[:,4]
	NIRc_SNR = output_zs[:,5]
	logNIRc_SNR = NIRc_SNR
	logNIRc_SNR[np.where(NIRc_SNR <= 0)[0]] = 1e-5
	logNIRc_SNR = np.log10(logNIRc_SNR)

	title_for_plot = 'BPZ Results'

# Get High SNR Objects
high_SNR_indices = np.where(NIRc_SNR >= SNR_limit)[0]
high_SNR_high_prob_indices = np.where((NIRc_SNR >= SNR_limit) & (z_prob > prob_limit))[0]
if (args.eazy_output_file):
	high_SNR_high_prob_lowq_z_indices = np.where((NIRc_SNR >= SNR_limit) & (z_prob > prob_limit) & (q_z <= q_z_limit))[0]
if (args.bpz_output_file):
	high_SNR_high_prob_lowchisq2_indices = np.where((NIRc_SNR >= SNR_limit) & (z_prob > prob_limit) & (chisq2 <= chisq2_limit))[0]

# High SNR 
high_SNR_outlier_IDs, high_SNR_outlier_zspec, high_SNR_outlier_zphot = find_catastrophic_outliers(ids[high_SNR_indices], z_spec[high_SNR_indices], z_phot[high_SNR_indices], 0.15)
high_SNR_outlier_indices = np.zeros(len(high_SNR_outlier_IDs), dtype = 'int')
for k in range(0,len(high_SNR_outlier_IDs)):
	high_SNR_outlier_indices[k] = np.where(ids == high_SNR_outlier_IDs[k])[0][0]
high_SNR_outlier_NIRc_SNR = NIRc_SNR[high_SNR_outlier_indices]

# High SNR, High Probability
high_SNR_high_prob_outlier_IDs, high_SNR_high_prob_outlier_zspec, high_SNR_high_prob_outlier_zphot = find_catastrophic_outliers(ids[high_SNR_high_prob_indices], z_spec[high_SNR_high_prob_indices], z_phot[high_SNR_high_prob_indices], 0.15)
high_SNR_high_prob_outlier_indices = np.zeros(len(high_SNR_high_prob_outlier_IDs), dtype = 'int')
for k in range(0,len(high_SNR_high_prob_outlier_IDs)):
	high_SNR_high_prob_outlier_indices[k] = np.where(ids == high_SNR_high_prob_outlier_IDs[k])[0][0]
high_SNR_high_prob_outlier_NIRc_SNR = NIRc_SNR[high_SNR_high_prob_outlier_indices]
high_SNR_high_prob_outlier_IDs = ids[high_SNR_high_prob_outlier_indices]

if (args.eazy_output_file):
	# High SNR, High Probability, Low Q_z
	high_SNR_high_prob_lowq_z_outlier_IDs, high_SNR_high_prob_lowq_z_outlier_zspec, high_SNR_high_prob_lowq_z_outlier_zphot = find_catastrophic_outliers(ids[high_SNR_high_prob_lowq_z_indices], z_spec[high_SNR_high_prob_lowq_z_indices], z_phot[high_SNR_high_prob_lowq_z_indices], 0.15)
	high_SNR_high_prob_lowq_z_outlier_indices = np.zeros(len(high_SNR_high_prob_lowq_z_outlier_IDs), dtype = 'int')
	for k in range(0,len(high_SNR_high_prob_lowq_z_outlier_IDs)):
		high_SNR_high_prob_lowq_z_outlier_indices[k] = np.where(ids == high_SNR_high_prob_lowq_z_outlier_IDs[k])[0][0]
	high_SNR_high_prob_lowq_z_outlier_logNIRc_SNR = logNIRc_SNR[high_SNR_high_prob_lowq_z_outlier_indices]
	high_SNR_high_prob_lowq_z_outlier_IDs = ids[high_SNR_high_prob_lowq_z_outlier_indices]
	
if (args.bpz_output_file):
	# High SNR, High Probability, Low Chisq2
	high_SNR_high_prob_lowchisq2_outlier_IDs, high_SNR_high_prob_lowchisq2_outlier_zspec, high_SNR_high_prob_lowchisq2_outlier_zphot = find_catastrophic_outliers(ids[high_SNR_high_prob_lowchisq2_indices], z_spec[high_SNR_high_prob_lowchisq2_indices], z_phot[high_SNR_high_prob_lowchisq2_indices], 0.15)
	high_SNR_high_prob_lowchisq2_outlier_indices = np.zeros(len(high_SNR_high_prob_lowchisq2_outlier_IDs), dtype = 'int')
	for k in range(0,len(high_SNR_high_prob_lowchisq2_outlier_IDs)):
		high_SNR_high_prob_lowchisq2_outlier_indices[k] = np.where(ids == high_SNR_high_prob_lowchisq2_outlier_IDs[k])[0][0]
	high_SNR_high_prob_lowchisq2_outlier_logNIRc_SNR = logNIRc_SNR[high_SNR_high_prob_lowchisq2_outlier_indices]
	high_SNR_high_prob_lowchisq2_outlier_IDs = ids[high_SNR_high_prob_lowchisq2_outlier_indices]


# I may have to update this at some point. Right now, the code looks at the outliers
# that are high SNR and high probability
if (args.eazy_output_file):
	outlier_z_IDs = high_SNR_high_prob_lowq_z_outlier_IDs#high_SNR_high_prob_outlier_IDs
	outlier_z_spec = high_SNR_high_prob_lowq_z_outlier_zspec#high_SNR_high_prob_outlier_zspec
	outlier_z_phot = high_SNR_high_prob_lowq_z_outlier_zphot#high_SNR_high_prob_outlier_zphot
	outlier_snr_colors = high_SNR_high_prob_lowq_z_outlier_logNIRc_SNR#high_SNR_high_prob_outlier_logNIRc_SNR
if (args.bpz_output_file):
	outlier_z_IDs = high_SNR_high_prob_lowchisq2_outlier_IDs#high_SNR_high_prob_outlier_IDs
	outlier_z_spec = high_SNR_high_prob_lowchisq2_outlier_zspec#high_SNR_high_prob_outlier_zspec
	outlier_z_phot = high_SNR_high_prob_lowchisq2_outlier_zphot#high_SNR_high_prob_outlier_zphot
	outlier_snr_colors = high_SNR_high_prob_lowchisq2_outlier_logNIRc_SNR#high_SNR_high_prob_outlier_logNIRc_SNR


# Figure Size
# Get current size
fig_size = plt.rcParams["figure.figsize"]
 
# Prints: [8.0, 6.0]
#print "Current size:", fig_size
fig_size[0] = 15
fig_size[1] = 12
plt.rcParams["figure.figsize"] = fig_size

# Plotting the Outliers
plt.ion()
subplot_kw = dict(xlim=(0, 15), ylim=(0, 15), autoscale_on=False)
fig, ax = plt.subplots(subplot_kw=subplot_kw)
if (args.eazy_output_file):
	ax.scatter(z_spec[high_SNR_high_prob_lowq_z_indices], z_phot[high_SNR_high_prob_lowq_z_indices], edgecolor = 'None', alpha = 0.2, s = 8, color = 'grey')
if (args.bpz_output_file):
	ax.scatter(z_spec[high_SNR_high_prob_lowchisq2_indices], z_phot[high_SNR_high_prob_lowchisq2_indices], edgecolor = 'None', alpha = 0.2, s = 8, color = 'grey')
#ax.scatter(z_spec[high_SNR_high_prob_indices], z_phot[high_SNR_high_prob_indices], edgecolor = 'None', alpha = 0.2, s = 8, color = 'grey')
#pts = ax.scatter(outlier_z_spec, outlier_z_phot, s=8, color = 'red', edgecolor = 'None')
#colorbar = ax.scatter(outlier_z_spec, outlier_z_phot, s=50, c = outlier_snr_colors, edgecolor = 'None')
pts = ax.scatter(outlier_z_spec, outlier_z_phot, s=30, c = outlier_snr_colors, edgecolor = 'None')
plt.plot([0, 15],[0, 12.6], '--', color = 'grey')
plt.plot([0, 15],[0, 17.4], '--', color = 'grey')
ax.set_xlabel('z$_{phot}$')
ax.set_ylabel('z$_{spec}$')
ax.set_title('LASSO SUBSAMPLES (PRESS ENTER TO QUIT)')
cbar = fig.colorbar(pts, label = 'log(SNR)')
selector = SelectFromCollection(ax, pts)

plt.draw()
raw_input('Press Enter to Quit')
selected_points = selector.xys[selector.ind]
selector.disconnect()


