import os
import sys
import math
import argparse
import astropy
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def get_SF_SED_filename_v11(SF_ID_number):

	if (SF_ID_number <= 50000):
		SF_filename = 'JADES_SF_mock_r1_v1.1_spec_5A_ID_1_50000.fits'
	if ((SF_ID_number > 50000) & (SF_ID_number <= 100000)):
		SF_filename = 'JADES_SF_mock_r1_v1.1_spec_5A_ID_50001_100000.fits'
	if ((SF_ID_number > 100000) & (SF_ID_number <= 150000)):
		SF_filename = 'JADES_SF_mock_r1_v1.1_spec_5A_ID_100001_150000.fits'
	if ((SF_ID_number > 150000) & (SF_ID_number <= 200000)):
		SF_filename = 'JADES_SF_mock_r1_v1.1_spec_5A_ID_150001_200000.fits'
	if ((SF_ID_number > 200000) & (SF_ID_number <= 250000)):
		SF_filename = 'JADES_SF_mock_r1_v1.1_spec_5A_ID_200001_250000.fits'
	if ((SF_ID_number > 250000) & (SF_ID_number <= 300000)):
		SF_filename = 'JADES_SF_mock_r1_v1.1_spec_5A_ID_250001_300000.fits'
	if ((SF_ID_number > 300000) & (SF_ID_number <= 302515)):
		SF_filename = 'JADES_SF_mock_r1_v1.1_spec_5A_ID_300001_302515.fits'
		
	return SF_filename

def get_SF_SED_filename_v12(SF_ID_number):

	if (SF_ID_number <= 44829):
		SF_filename = 'JADES_SF_mock_r1_v1.2_spec_5A_30um_z_0p2_1.fits.gz'
	if ((SF_ID_number > 44829) & (SF_ID_number <= 81779)):
		SF_filename = 'JADES_SF_mock_r1_v1.2_spec_5A_30um_z_1_1p5.fits.gz'
	if ((SF_ID_number > 81779) & (SF_ID_number <= 115980)):
		SF_filename = 'JADES_SF_mock_r1_v1.2_spec_5A_30um_z_1p5_2.fits.gz'
	if ((SF_ID_number > 115980) & (SF_ID_number <= 172661)):
		SF_filename = 'JADES_SF_mock_r1_v1.2_spec_5A_30um_z_2_3.fits.gz'
	if ((SF_ID_number > 172661) & (SF_ID_number <= 216240)):
		SF_filename = 'JADES_SF_mock_r1_v1.2_spec_5A_30um_z_3_4.fits.gz'
	if ((SF_ID_number > 216248) & (SF_ID_number <= 248165)):
		SF_filename = 'JADES_SF_mock_r1_v1.2_spec_5A_30um_z_4_5.fits.gz'
	if ((SF_ID_number > 248165) & (SF_ID_number <= 302514)):
		SF_filename = 'JADES_SF_mock_r1_v1.2_spec_5A_30um_z_5_15.fits.gz'
		
	return SF_filename

JAGUAR_version = 'r1_v1.1'
JAGUAR_noLyA_filename = 'JADES_SF_mock_r1_v1.1_fluxes_noLya.fits'

def FluxtoABMag(flux):
	return (-5.0 / 2.0) * np.log10(flux) - 48.60

# speed of light in angstroms per second
c = 2.99792458e18

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

# Input folder
parser.add_argument(
  '-sedin','--sedinputfolder',
  help="Input Folder With Mock Catalog SED Files",
  action="store",
  type=str,
  dest="SED_input_folder",
  required=True
)

######################
# Optional Arguments #
######################

# SF Object ID?
parser.add_argument(
  '-SFID','--SFID',
  help="SF Object ID",
  action="store",
  type=int,
  dest="SF_ID",
  required=False
)

# Q Object ID?
parser.add_argument(
  '-QID','--QID',
  help="Q Object ID",
  action="store",
  type=int,
  dest="Q_ID",
  required=False
)

# Plot the Filters on Top of the Plot?
parser.add_argument(
  '-filt','--plot_filters',
  help="Plot the Filters on Top of the Plot?",
  action="store_true",
  dest="plot_filters",
  required=False
)

# Use the non-Lyman alpha photometry?
parser.add_argument(
  '-noLyA','--noLyA',
  help="Use the non-Lyman alpha photometry?",
  action="store_true",
  dest="noLyA",
  required=False
)

# Plot in AB Magnitude Units?
parser.add_argument(
  '-abmag','--abmag',
  help="Plot in AB Magnitude Units?",
  action="store_true",
  dest="abmag",
  required=False
)

# Save plot to png?
parser.add_argument(
  '-sp','--save_plot',
  help="Save plot to png?",
  action="store_true",
  dest="save_plot",
  required=False
)


args=parser.parse_args()

SF_ID = args.SF_ID#269865
quiescent_ID = args.Q_ID#309977


# The centroids and errors for the photometry
apparent_centroids_hst = np.array([0.4297, 0.5907, 0.7764, 0.8333, 0.9445])
apparent_error_hst = np.array([0.1038, 0.2342, 0.1528, 0.2511, 0.1229])/2.0

apparent_centroids_nircam = np.array([0.902, 1.154, 1.501, 1.989, 2.762, 3.362, 3.568, 4.082, 4.408])
apparent_error_nircam = np.array([0.132, 0.225, 0.318, 0.457, 0.683, 0.352, 0.781, 0.438, 1.029])/2.0

apparent_centroids_IRAC = np.array([3.550, 4.493])

# Quiescent Mock File
if (args.Q_ID):
	full_q_mock_file = args.input_folder+'JADES_Q_mock_'+JAGUAR_version+'.fits'
	testfits = fits.open(full_q_mock_file)
	full_q_IDs = testfits[1].data['ID']
	full_q_redshifts = testfits[1].data['redshift']
	full_q_logmass = testfits[1].data['mStar']
	full_q_F435W = testfits[1].data['HST_F435W_fnu']
	full_q_F606W = testfits[1].data['HST_F606W_fnu']
	full_q_F775W = testfits[1].data['HST_F775W_fnu']
	full_q_F814W = testfits[1].data['HST_F814W_fnu']
	full_q_F850LP = testfits[1].data['HST_F850LP_fnu']
	full_q_F070W = testfits[1].data['NRC_F070W_fnu']
	full_q_F090W = testfits[1].data['NRC_F090W_fnu']
	full_q_F115W = testfits[1].data['NRC_F115W_fnu']
	full_q_F150W = testfits[1].data['NRC_F150W_fnu']
	full_q_F200W = testfits[1].data['NRC_F200W_fnu']
	full_q_F277W = testfits[1].data['NRC_F277W_fnu']
	full_q_F335M = testfits[1].data['NRC_F335M_fnu']
	full_q_F356W = testfits[1].data['NRC_F356W_fnu']
	#full_q_fake_mag = testfits[1].data['NRC_F410M_F444W_fnu']
	full_q_F410M = testfits[1].data['NRC_F410M_fnu']
	full_q_F444W = testfits[1].data['NRC_F444W_fnu']
	full_q_sSFR = 9.0 + testfits[1].data['sSFR']
	
	full_q_IRAC_3p6 = testfits[1].data['IRAC_3p6_fnu']
	full_q_IRAC_4p5 = testfits[1].data['IRAC_4p5_fnu']

	n_objects_q = len(full_q_F435W)
	q_ID = np.where(full_q_IDs == quiescent_ID)[0][0]
	q_redshift = full_q_redshifts[q_ID]
	quiescent_mass = full_q_logmass[q_ID]
	
	# Quiescent SEDs
	full_q_SED = args.SED_input_folder + 'JADES_Q_mock_'+JAGUAR_version+'_spec_5A.fits'
	testfits = fits.open(full_q_SED)
	gal_properties = testfits[3].data
	wavelength = testfits[2].data
	spectrum = testfits[1].data
	
	quiescent_SED_ID = np.where(q_redshift == gal_properties['redshift'])[0][0]
	spectrum_obj = spectrum[quiescent_SED_ID] # in erg/s/cm^2/Angstrom
	q_spec_array_fnu = (spectrum_obj/(1.0+q_redshift) * ((1.0+q_redshift) * wavelength)**2 / c) 
	q_wave_z = np.array((1.0+q_redshift) * wavelength)/10000
	q_flux_nJy = np.array(q_spec_array_fnu/1e-23/1e-9)
	
	if (args.abmag):
		q_apparent_ABmag_hst = np.empty(5)
		q_apparent_ABmag_hst[0] = FluxtoABMag(full_q_F435W[q_ID]*1e-23*1e-9)
		q_apparent_ABmag_hst[1] = FluxtoABMag(full_q_F606W[q_ID]*1e-23*1e-9)
		q_apparent_ABmag_hst[2] = FluxtoABMag(full_q_F775W[q_ID]*1e-23*1e-9)
		q_apparent_ABmag_hst[3] = FluxtoABMag(full_q_F814W[q_ID]*1e-23*1e-9)
		q_apparent_ABmag_hst[4] = FluxtoABMag(full_q_F850LP[q_ID]*1e-23*1e-9)
		
		q_apparent_ABmag_nircam = np.empty(9)
		q_apparent_ABmag_nircam[0] = FluxtoABMag(full_q_F090W[q_ID]*1e-23*1e-9)
		q_apparent_ABmag_nircam[1] = FluxtoABMag(full_q_F115W[q_ID]*1e-23*1e-9)
		q_apparent_ABmag_nircam[2] = FluxtoABMag(full_q_F150W[q_ID]*1e-23*1e-9)
		q_apparent_ABmag_nircam[3] = FluxtoABMag(full_q_F200W[q_ID]*1e-23*1e-9)
		q_apparent_ABmag_nircam[4] = FluxtoABMag(full_q_F277W[q_ID]*1e-23*1e-9)
		q_apparent_ABmag_nircam[5] = FluxtoABMag(full_q_F335M[q_ID]*1e-23*1e-9)
		q_apparent_ABmag_nircam[6] = FluxtoABMag(full_q_F356W[q_ID]*1e-23*1e-9)
		q_apparent_ABmag_nircam[7] = FluxtoABMag(full_q_F410M[q_ID]*1e-23*1e-9)
		q_apparent_ABmag_nircam[8] = FluxtoABMag(full_q_F444W[q_ID]*1e-23*1e-9)
		
	else:
		q_apparent_fnu_nJy_hst = np.empty(5)
		q_apparent_fnu_nJy_hst[0] = full_q_F435W[q_ID]
		q_apparent_fnu_nJy_hst[1] = full_q_F606W[q_ID]
		q_apparent_fnu_nJy_hst[2] = full_q_F775W[q_ID]
		q_apparent_fnu_nJy_hst[3] = full_q_F814W[q_ID]
		q_apparent_fnu_nJy_hst[4] = full_q_F850LP[q_ID]
		
		q_apparent_fnu_nJy_nircam = np.empty(9)
		q_apparent_fnu_nJy_nircam[0] = full_q_F090W[q_ID]
		q_apparent_fnu_nJy_nircam[1] = full_q_F115W[q_ID]
		q_apparent_fnu_nJy_nircam[2] = full_q_F150W[q_ID]
		q_apparent_fnu_nJy_nircam[3] = full_q_F200W[q_ID]
		q_apparent_fnu_nJy_nircam[4] = full_q_F277W[q_ID]
		q_apparent_fnu_nJy_nircam[5] = full_q_F335M[q_ID]
		q_apparent_fnu_nJy_nircam[6] = full_q_F356W[q_ID]
		q_apparent_fnu_nJy_nircam[7] = full_q_F410M[q_ID]
		q_apparent_fnu_nJy_nircam[8] = full_q_F444W[q_ID]

	

# SF Mock File
if (args.SF_ID):
	full_sf_mock_file = args.input_folder+'JADES_SF_mock_'+JAGUAR_version+'.fits'
	testfits = fits.open(full_sf_mock_file)
	full_sf_IDs = testfits[1].data['ID']
	full_sf_redshifts = testfits[1].data['redshift']
	full_sf_logmass = testfits[1].data['mStar']
	full_sf_MUV = testfits[1].data['MUV']
	full_sf_sSFR = 9.0 + testfits[1].data['sSFR']
	if (args.noLyA):
		full_sf_mock_file = args.input_folder+JAGUAR_noLyA_filename
		testfits = fits.open(full_sf_mock_file)
	full_sf_F435W = testfits[1].data['HST_F435W_fnu']
	full_sf_F606W = testfits[1].data['HST_F606W_fnu']
	full_sf_F775W = testfits[1].data['HST_F775W_fnu']
	full_sf_F814W = testfits[1].data['HST_F814W_fnu']
	full_sf_F850LP = testfits[1].data['HST_F850LP_fnu']
	full_sf_F070W = testfits[1].data['NRC_F070W_fnu']
	full_sf_F090W = testfits[1].data['NRC_F090W_fnu']
	full_sf_F115W = testfits[1].data['NRC_F115W_fnu']
	full_sf_F150W = testfits[1].data['NRC_F150W_fnu']
	full_sf_F200W = testfits[1].data['NRC_F200W_fnu']
	full_sf_F277W = testfits[1].data['NRC_F277W_fnu']
	full_sf_F335M = testfits[1].data['NRC_F335M_fnu']
	full_sf_F356W = testfits[1].data['NRC_F356W_fnu']
	full_sf_fake_mag = testfits[1].data['NRC_F410M_F444W_fnu']
	full_sf_F410M = testfits[1].data['NRC_F410M_fnu']
	full_sf_F444W = testfits[1].data['NRC_F444W_fnu']
	
	full_sf_IRAC_3p6 = testfits[1].data['IRAC_3p6_fnu']#['_F444W_APP']
	full_sf_IRAC_4p5 = testfits[1].data['IRAC_4p5_fnu']#['_F444W_APP']


	n_objects_sf = len(full_sf_F435W)
	sf_ID = np.where(full_sf_IDs == float(SF_ID))[0][0]
	sf_redshift = full_sf_redshifts[sf_ID]

	# SF SEDs 
	if (JAGUAR_version == 'r1_v1.1'):
		full_sf_SED = args.SED_input_folder + get_SF_SED_filename_v11(SF_ID)
	if (JAGUAR_version == 'r1_v1.2'):
		full_sf_SED = args.SED_input_folder + get_SF_SED_filename_v12(SF_ID)
	testfits = fits.open(full_sf_SED)
	gal_properties = testfits[3].data
	wavelength = testfits[2].data
	spectrum = testfits[1].data
	
	sf_SED_ID = np.where(gal_properties['ID'] == SF_ID)[0][0]
	spectrum_obj = spectrum[sf_SED_ID] # in erg/s/cm^2/Angstrom
	
	sf_spec_array_fnu = (spectrum_obj/(1.0+sf_redshift) * ((1.0+sf_redshift) * wavelength)**2 / c) 
	sf_wave_z = np.array((1.0+sf_redshift) * wavelength)/10000
	sf_flux_nJy = np.array(sf_spec_array_fnu/1e-23/1e-9)

	apparent_fnu_nircam_error = np.array([0.52, 0.41, 0.37, 0.36, 0.53, 1.00, 0.57, 0.83, 0.76])

	if (args.abmag):
		sf_apparent_ABmag_hst = np.empty(5)
		sf_apparent_ABmag_hst[0] = FluxtoABMag(full_sf_F435W[sf_ID]*1e-23*1e-9)
		sf_apparent_ABmag_hst[1] = FluxtoABMag(full_sf_F606W[sf_ID]*1e-23*1e-9)
		sf_apparent_ABmag_hst[2] = FluxtoABMag(full_sf_F775W[sf_ID]*1e-23*1e-9)
		sf_apparent_ABmag_hst[3] = FluxtoABMag(full_sf_F814W[sf_ID]*1e-23*1e-9)
		sf_apparent_ABmag_hst[4] = FluxtoABMag(full_sf_F850LP[sf_ID]*1e-23*1e-9)
	
		sf_apparent_ABmag_nircam = np.empty(9)
		sf_apparent_ABmag_nircam[0] = FluxtoABMag(full_sf_F090W[sf_ID]*1e-23*1e-9)
		sf_apparent_ABmag_nircam[1] = FluxtoABMag(full_sf_F115W[sf_ID]*1e-23*1e-9)
		sf_apparent_ABmag_nircam[2] = FluxtoABMag(full_sf_F150W[sf_ID]*1e-23*1e-9)
		sf_apparent_ABmag_nircam[3] = FluxtoABMag(full_sf_F200W[sf_ID]*1e-23*1e-9)
		sf_apparent_ABmag_nircam[4] = FluxtoABMag(full_sf_F277W[sf_ID]*1e-23*1e-9)
		sf_apparent_ABmag_nircam[5] = FluxtoABMag(full_sf_F335M[sf_ID]*1e-23*1e-9)
		sf_apparent_ABmag_nircam[6] = FluxtoABMag(full_sf_F356W[sf_ID]*1e-23*1e-9)
		sf_apparent_ABmag_nircam[7] = FluxtoABMag(full_sf_F410M[sf_ID]*1e-23*1e-9)
		sf_apparent_ABmag_nircam[8] = FluxtoABMag(full_sf_F444W[sf_ID]*1e-23*1e-9)

		sf_apparent_ABmag_IRAC = np.empty(2)
		sf_apparent_ABmag_IRAC[0] = FluxtoABMag(full_sf_IRAC_3p6[sf_ID]*1e-23*1e-9)
		sf_apparent_ABmag_IRAC[1] = FluxtoABMag(full_sf_IRAC_4p5[sf_ID]*1e-23*1e-9)

		sf_flux_ABmag = FluxtoABMag(sf_spec_array_fnu)
	else:
		sf_apparent_fnu_nJy_hst = np.empty(5)
		sf_apparent_fnu_nJy_hst[0] = full_sf_F435W[sf_ID]
		sf_apparent_fnu_nJy_hst[1] = full_sf_F606W[sf_ID]
		sf_apparent_fnu_nJy_hst[2] = full_sf_F775W[sf_ID]
		sf_apparent_fnu_nJy_hst[3] = full_sf_F814W[sf_ID]
		sf_apparent_fnu_nJy_hst[4] = full_sf_F850LP[sf_ID]
	
		sf_apparent_fnu_nJy_nircam = np.empty(9)
		sf_apparent_fnu_nJy_nircam[0] = full_sf_F090W[sf_ID]
		sf_apparent_fnu_nJy_nircam[1] = full_sf_F115W[sf_ID]
		sf_apparent_fnu_nJy_nircam[2] = full_sf_F150W[sf_ID]
		sf_apparent_fnu_nJy_nircam[3] = full_sf_F200W[sf_ID]
		sf_apparent_fnu_nJy_nircam[4] = full_sf_F277W[sf_ID]
		sf_apparent_fnu_nJy_nircam[5] = full_sf_F335M[sf_ID]
		sf_apparent_fnu_nJy_nircam[6] = full_sf_F356W[sf_ID]
		sf_apparent_fnu_nJy_nircam[7] = full_sf_F410M[sf_ID]
		sf_apparent_fnu_nJy_nircam[8] = full_sf_F444W[sf_ID]

		sf_apparent_fnu_nJy_IRAC = np.empty(2)
		sf_apparent_fnu_nJy_IRAC[0] = full_sf_IRAC_3p6[sf_ID]
		sf_apparent_fnu_nJy_IRAC[1] = full_sf_IRAC_4p5[sf_ID]


# Open the filter file:
filter_names = 'HST_NIRCam_filters.txt'
tab = ascii.read(filter_names, names=['filter_names','filter_link'], data_start=0)
filter_names = tab['filter_names'] 
filter_link = tab['filter_link'] 
n_filters = len(filter_names)
hst_color = 'red'
hst_alpha = 0.3
nircam_color = 'blue'
nircam_alpha = 0.5
filter_alphas = [0.8, 0.6, 0.4, 0.2, 0.1, nircam_alpha, nircam_alpha, nircam_alpha, nircam_alpha, nircam_alpha, nircam_alpha, nircam_alpha, nircam_alpha, nircam_alpha, nircam_alpha]
filter_colors = ['black','black','black','black','black','indigo','violet','blue','green','greenyellow','yellow','orange','red','firebrick']

hst_filter_colors = ['black','black','black','black','black']
nircam_filter_colors = ['indigo','violet','blue','green','greenyellow','yellow','orange','red','firebrick']


if ((args.SF_ID is not None) & (args.Q_ID is not None)):

	fig1 = plt.figure(figsize=(8,10))
	ax1 = fig1.add_subplot(211)

else:

	fig1 = plt.figure(figsize=(9,6))
	ax1 = fig1.add_subplot(111)

if ((args.SF_ID is not None)):

	if (args.abmag):
		# Star Forming Galaxy
		# HST
		ax1.scatter(apparent_centroids_hst, sf_apparent_ABmag_hst, linewidths = 0, color = hst_filter_colors, alpha = 0.3, s = 80, zorder = 10)
		
		# NIRCam 
		ax1.scatter(apparent_centroids_nircam, sf_apparent_ABmag_nircam, linewidths = 0, color = nircam_filter_colors, s = 80, zorder = 10)
		
		# IRAC 
		ax1.scatter(apparent_centroids_IRAC, sf_apparent_ABmag_IRAC, linewidths = 0, color = 'black', s = 80, zorder = 10)
		
		# SED
		ax1.plot(sf_wave_z, sf_flux_ABmag, '-', color='black', alpha = 0.5)
		
		flux_max = np.max(sf_apparent_ABmag_nircam)
		y_max = flux_max
		flux_min = np.min(sf_apparent_ABmag_nircam)
		y_min = flux_min
		
		if (y_max > 30):
			y_max = 30.0
		
		if (args.plot_filters):
			for filter in range(0, n_filters):
				ind_filter_file = filter_link[filter]
				ind_filter_file_full = np.loadtxt(ind_filter_file)
				bp_wave = ind_filter_file_full[:,0]/10000
				bp_throughput = (-0.07)*(ind_filter_file_full[:,1]/np.max(ind_filter_file_full[:,1])*flux_max/2.5)+(y_max*1.05)+0.02#31.02-2
				ax1.plot(bp_wave, bp_throughput, color = filter_colors[filter], alpha = filter_alphas[filter])
				ax1.fill(bp_wave, bp_throughput, facecolor=filter_colors[filter], alpha=0.2)
		
		ax1.set_xlabel('Observed Wavelength (Microns)')
		ax1.set_ylabel(r'AB Mag')
		ax1.set_xlim([0.2,5.2])
		#ax1.set_ylim([31-2,27-2])
		ax1.set_ylim([y_max*1.05,y_min/1.1])
		
		if (.1216*(1.0+sf_redshift) < 5.2):
			ax1.text(.1216*(1.0+sf_redshift)+0.05, 28.9-2,r'Ly$\alpha$', color='black', alpha = 0.5)
		if (.6563*(1.0+sf_redshift) < 5.2):
			ax1.text(.6563*(1.0+sf_redshift)+0.05, 28.9-2,r'H$\alpha$', color='black', alpha = 0.5)
		if (.5007*(1.0+sf_redshift) < 5.2):
			ax1.text(.5007*(1.0+sf_redshift)+0.05, 28.9-2,'[OIII]', color='black', alpha = 0.5)
				
	else:
		# Star Forming Galaxy
		# HST
		ax1.scatter(apparent_centroids_hst, sf_apparent_fnu_nJy_hst, linewidths = 0, color = hst_filter_colors, alpha = 0.3, s = 80, zorder = 10)
		
		# NIRCam 
		ax1.scatter(apparent_centroids_nircam, sf_apparent_fnu_nJy_nircam, linewidths = 0, color = nircam_filter_colors, s = 80, zorder = 10)
		
		# IRAC 
		ax1.scatter(apparent_centroids_IRAC, sf_apparent_fnu_nJy_IRAC, linewidths = 0, color = 'black', s = 80, zorder = 10)
		
		# SED
		ax1.plot(sf_wave_z, sf_flux_nJy, '-', color='black', alpha = 0.5)

		flux_min = np.min(sf_apparent_fnu_nJy_nircam)
		y_min = flux_min
		flux_max = np.max(sf_apparent_fnu_nJy_nircam)
		y_max = flux_max
		
		if (y_max < 1.0):
			y_max = 1.0

		if (y_min < 1e-1):
			y_min = 1e-1

		if (args.plot_filters):
			for filter in range(0, n_filters):
				ind_filter_file = filter_link[filter]
				ind_filter_file_full = np.loadtxt(ind_filter_file)
				bp_wave = ind_filter_file_full[:,0]/10000
				bp_throughput = (0.2)*(ind_filter_file_full[:,1]/np.max(ind_filter_file_full[:,1])*flux_max/2.5)+(y_min/2)+0.02#31.02-2
				ax1.plot(bp_wave, bp_throughput, color = filter_colors[filter], alpha = filter_alphas[filter])
				ax1.fill(bp_wave, bp_throughput, facecolor=filter_colors[filter], alpha=0.2)

		ax1.set_xlabel('Observed Wavelength (Microns)')
		ax1.set_ylabel(r'F$_{\nu}$ (nJy)')
		ax1.set_xlim([0.2,5.2])
		ax1.semilogy()
		#ax1.set_ylim([31-2,27-2])
		ax1.set_ylim([y_min/2,y_max*2])

		if (.1216*(1.0+sf_redshift) < 5.2):
			ax1.text(.1216*(1.0+sf_redshift)+0.05, 60,r'Ly$\alpha$', color='black', alpha = 0.5)
		if (.6563*(1.0+sf_redshift) < 5.2):
			ax1.text(.6563*(1.0+sf_redshift)+0.05, 60,r'H$\alpha$', color='black', alpha = 0.5)
		if (.5007*(1.0+sf_redshift) < 5.2):
			ax1.text(.5007*(1.0+sf_redshift)+0.05, 60,'[OIII]', color='black', alpha = 0.5)
			
	ax1.text(0.05, 0.9,'Star-Forming Galaxy '+str(int(SF_ID)),
		horizontalalignment='left',
		verticalalignment='center',
		transform = ax1.transAxes)
	ax1.text(0.05, 0.85,r'$z$ = '+str(round(full_sf_redshifts[sf_ID],3)),
		horizontalalignment='left',
		verticalalignment='center',
		transform = ax1.transAxes)
	ax1.text(0.05, 0.8,r'log(sSFR/Gyr$^{-1}$) = '+str(round(full_sf_sSFR[sf_ID], 3)),
		horizontalalignment='left',
		verticalalignment='center',
		transform = ax1.transAxes)
		  			
	
	if (args.plot_filters):
		ax1.text(0.1, 0.25,r'HST',
		     horizontalalignment='right',
		     verticalalignment='center',
		     color='red',
		     alpha = 0.3,
		     transform = ax1.transAxes)
		ax1.text(0.95, 0.25,r'NIRCam',
		     horizontalalignment='right',
		     verticalalignment='center',
		     color='blue',
		     transform = ax1.transAxes)
		
			
if ((args.SF_ID is not None) & (args.Q_ID is not None)):
			
	ax1 = fig1.add_subplot(212)
	
if ((args.Q_ID is not None)):
	# Quiescent Galaxy
	# HST
	if (args.abmag):
		ax1.scatter(apparent_centroids_hst, q_apparent_ABmag_hst, linewidths = 0, color = hst_filter_colors, alpha = 0.3, s = 80, zorder = 10)
		
		# NIRCam 
		ax1.scatter(apparent_centroids_nircam, q_apparent_ABmag_nircam, linewidths = 0, color = nircam_filter_colors, s = 80, zorder = 10)
		
		# SED
		ax1.plot(q_wave_z, q_flux_ABmag, '-', color='black', alpha = 0.5)
		
		flux_max = np.max(q_apparent_ABmag_nircam)
		y_max = flux_max
		flux_min = np.min(q_apparent_ABmag_nircam)
		y_min = flux_min
	
		if (args.plot_filters):
			for filter in range(0, n_filters):
				ind_filter_file = filter_link[filter]
				ind_filter_file_full = np.loadtxt(ind_filter_file)
				bp_wave = ind_filter_file_full[:,0]/10000
				bp_throughput = (-0.13)*(ind_filter_file_full[:,1]/np.max(ind_filter_file_full[:,1])*flux_max/2.5)+(y_max*1.05)+0.02#28.02
				ax1.plot(bp_wave, bp_throughput, color = filter_colors[filter], alpha = filter_alphas[filter])
				ax1.fill(bp_wave, bp_throughput, facecolor=filter_colors[filter], alpha=0.2)
					
		ax1.set_xlabel('Observed Wavelength (Microns)')
		ax1.set_ylabel(r'AB Mag')
		ax1.set_xlim([0.2,5.2])
		#ax1.set_ylim([28,22])
		ax1.set_ylim([y_max*1.05,y_min/1.1])
				
	else:
		ax1.scatter(apparent_centroids_hst, q_apparent_fnu_nJy_hst, linewidths = 0, color = hst_filter_colors, alpha = 0.3, s = 80, zorder = 10)
		
		# NIRCam 
		ax1.scatter(apparent_centroids_nircam, q_apparent_fnu_nJy_nircam, linewidths = 0, color = nircam_filter_colors, s = 80, zorder = 10)
		
		# SED
		ax1.plot(q_wave_z, q_flux_nJy, '-', color='black', alpha = 0.5)
		
		flux_min = np.min(q_apparent_fnu_nJy_nircam)
		y_min = flux_min
		flux_max = np.max(q_apparent_fnu_nJy_nircam)
		y_max = flux_max
	
		if (args.plot_filters):
			for filter in range(0, n_filters):
				ind_filter_file = filter_link[filter]
				ind_filter_file_full = np.loadtxt(ind_filter_file)
				bp_wave = ind_filter_file_full[:,0]/10000
				bp_throughput = (0.2)*(ind_filter_file_full[:,1]/np.max(ind_filter_file_full[:,1])*flux_max/2.5)+(y_min/2.0)+0.02#28.02
				ax1.plot(bp_wave, bp_throughput, color = filter_colors[filter], alpha = filter_alphas[filter])
				ax1.fill(bp_wave, bp_throughput, facecolor=filter_colors[filter], alpha=0.2)
					
		ax1.set_xlabel('Observed Wavelength (Microns)')
		ax1.set_ylabel(r'F$_{\nu}$ (nJy)')
		ax1.set_xlim([0.2,5.2])
		ax1.semilogy()
		#ax1.set_ylim([28,22])
		ax1.set_ylim([y_min/2,y_max*2.0])
		
	ax1.text(0.05, 0.9,'Quiescent Galaxy '+str(int(quiescent_ID)),
	     horizontalalignment='left',
	     verticalalignment='center',
	     transform = ax1.transAxes)
	ax1.text(0.05, 0.85,r'$z$ = '+str(round(full_q_redshifts[q_ID],3)),
	     horizontalalignment='left',
	     verticalalignment='center',
	     transform = ax1.transAxes)
	ax1.text(0.05, 0.8,r'log(M / M$_\odot$) = '+str(round(quiescent_mass, 3)),
	     horizontalalignment='left',
	     verticalalignment='center',
	     transform = ax1.transAxes)
	 			
	if (args.plot_filters):
		ax1.text(0.1, 0.25,r'HST',
		     horizontalalignment='right',
		     verticalalignment='center',
		     color='red',
		     alpha = 0.3,
		     transform = ax1.transAxes)
		ax1.text(0.95, 0.25,r'NIRCam',
		     horizontalalignment='right',
		     verticalalignment='center',
		     color='blue',
		     transform = ax1.transAxes)

	
plt.tight_layout()

if (args.save_plot):
	if ((args.SF_ID is not None) & (args.Q_ID is not None)):
		object_name = 'SF_'+str(SF_ID)+'_Q_'+str(quiescent_ID)+'_SED.png'
	elif ((args.SF_ID is None) & (args.Q_ID is not None)):
		object_name = 'Q_'+str(quiescent_ID)+'_SED.png'
	elif ((args.SF_ID is not None) & (args.Q_ID is None)):
		object_name = 'SF_'+str(SF_ID)+'_SED.png'
	plt.savefig(object_name, dpi=300)
else:
	plt.show()

