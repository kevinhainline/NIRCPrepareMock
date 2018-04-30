import sys
import argparse
import random
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import fits
from astropy.io import ascii

def JytoABMag(flux):
	return (-5.0 / 2.0) * np.log10(flux) - 48.60

randomize_IDs = 0

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


######################
# Optional Arguments #
######################

# SF Output Filename
parser.add_argument(
  '-sfo','--sfoutputfilename',
  help="SF Output Filename",
  action="store",
  type=str,
  dest="sf_output_filename",
  required=False
)

# Q Output Filename
parser.add_argument(
  '-qo','--qoutputfilename',
  help="Q Output Filename",
  action="store",
  type=str,
  dest="q_output_filename",
  required=False
)

# Randomize IDs?
parser.add_argument(
  '-rid','--randomizeids',
  help="Randomize the ID Numbers? (1 = Yes, 0 = No)",
  action="store",
  type=float,
  dest="randomize_IDs",
  required=False
)

# Number of Objects per Redshift Bin?
parser.add_argument(
  '-nz','--nobjectsperz',
  help="Number of Objects Per Redshift Bin)",
  action="store",
  type=float,
  dest="number_per_redshift_bin",
  required=False
)

# Filters File
parser.add_argument(
  '-fi','--filters_file',
  help="Filters File Name",
  action="store",
  type=str,
  dest="filters_file",
  required=False
)

args=parser.parse_args()

if (args.sf_output_filename):
	make_sf = 1
else:
	make_sf = 0

if (args.q_output_filename):
	make_q = 1
else:
	make_q = 0

	
if ((make_sf == 0) & (make_q == 0)):
	sys.exit("No output SF or Q filename specified! Use -sfo or -q to provide a filename")

# Read in filters file

if (args.filters_file):
	filter_file_name = args.filters_file
	filters_full = np.loadtxt(filter_file_name, dtype='str')
	filters = filters_full
else:
	filters = ['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 
		'HST_F850LP', 'NRC_F070W', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W',
		'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 
		'NRC_F410M_F444W', 'NRC_F410M', 'NRC_F444W']

number_filters = len(filters)

#sf_output_filename = 'sf_fluxes_4_13_18.dat'
#q_output_filename = 'q_fluxes_4_13_18.dat'

#make_sf = 0
#make_q = 1
#randomize_IDs = 0

redshift_minimum = 0.0
redshift_maximum = 15.0
redshift_range = np.arange(redshift_minimum, redshift_maximum+1.0, 1.0)
n_redshift_bins = len(redshift_range - 1)
if (args.number_per_redshift_bin):
	number_per_redshift_bin = args.number_per_redshift_bin#1000
else:
	number_per_redshift_bin = 1000
	
mass_minimum = 6.0
mass_maximum = 12.0
mass_range = np.arange(mass_minimum, mass_maximum+1.0, 1.0)
n_mass_bins = len(mass_range)-1
number_per_mass_bin = number_per_redshift_bin / (n_mass_bins)

print number_per_mass_bin

# Full star-forming catalogue, fits file
full_sf_mock_file = args.input_folder+'JADES_SF_mock_r1_v1.1.fits'
sftestfits = fits.open(full_sf_mock_file)

# Get redshifts and galaxy properties
full_sf_IDs = sftestfits[1].data['ID']
full_sf_redshifts = sftestfits[1].data['redshift']
full_sf_logmass = sftestfits[1].data['mStar']
n_sf_objects = len(full_sf_IDs)

# Get the star-forming fluxes. 
sf_filt_flux = np.zeros([number_filters, n_sf_objects])
for j in range(0, number_filters):
	sf_filt_flux[:][j] = sftestfits[1].data[filters[j]+'_fnu']*1e-23*1e-9


# Full quiescent catalogue, fits file
full_q_mock_file = args.input_folder+'JADES_Q_mock_r1_v1.1.fits'
qtestfits = fits.open(full_q_mock_file)

# Get redshifts and galaxy properties
full_q_IDs = qtestfits[1].data['ID']
full_q_redshifts = qtestfits[1].data['redshift']
full_q_logmass = qtestfits[1].data['mStar']
n_q_objects = len(full_q_IDs)

# Get the quiescent fluxes. 
q_filt_flux = np.zeros([number_filters, n_q_objects])
for j in range(0, number_filters):
	q_filt_flux[:][j] = qtestfits[1].data[filters[j]+'_fnu']*1e-23*1e-9


# Open up the output files. 
if (make_sf == 1):
	sffile = open(args.sf_output_filename, 'w')
if (make_q == 1):
	qfile = open(args.q_output_filename, 'w')

# This is to be used if the randomized fluxes are requested
sf_ID_value = 1
q_ID_value = 1

# Go through all of the redshift bins.
for z in range(0, n_redshift_bins-1):
	zbin_min = redshift_range[z]
	zbin_max = redshift_range[z+1]
	print "z = "+str(zbin_min)+" - "+str(zbin_max)
	
	# Find the indices for the objects in the redshift bin being examined.
	z_indices_sf = np.where((full_sf_redshifts > zbin_min) & (full_sf_redshifts <= zbin_max))[0]
	z_indices_q = np.where((full_q_redshifts > zbin_min) & (full_q_redshifts <= zbin_max))[0]

	# Go through the masses. 
	for m in range(0, n_mass_bins-1):
		massbin_min = mass_range[m]
		massbin_max = mass_range[m+1]
		print "   log(Mass) = "+str(massbin_min)+" - "+str(massbin_max)
		
		# Find the indices for the objects in the mass bin being examined. 
		mass_indices_sf = np.where((full_sf_logmass[z_indices_sf] > massbin_min) & (full_sf_logmass[z_indices_sf] <= massbin_max))[0]
		full_sf_IDs_z_mass = full_sf_IDs[z_indices_sf][mass_indices_sf]
		full_sf_redshifts_z_mass = full_sf_redshifts[z_indices_sf][mass_indices_sf]
		full_sf_logmass_z_mass = full_sf_logmass[z_indices_sf][mass_indices_sf]

		# Find the SF flux values for the objects in the mass and redshift bin being examined.
		sf_filt_flux_z_mass = np.zeros([number_filters, len(full_sf_IDs_z_mass)])
		for j in range(0, number_filters):
			flux_1D_array = sf_filt_flux[:][j]
			sf_filt_flux_z_mass[:][j] = flux_1D_array[z_indices_sf][mass_indices_sf]

		# Now, choose some objects randomly from within the bin (providing there are enough objects)
		if (len(full_sf_redshifts_z_mass) < number_per_mass_bin):
			ri = random.sample(xrange(0,len(full_sf_redshifts_z_mass)), len(full_sf_redshifts_z_mass))
		else:
			ri = random.sample(xrange(0,len(full_sf_redshifts_z_mass)), number_per_mass_bin)

		# Now, go and find the ID, redshift, mass, etc, and write it to a file
		if (make_sf == 1):
			print "    SF ri = "+str(len(ri))
			for x in range(0, len(ri)):
				if (randomize_IDs == 0):
					sf_ID_value = full_sf_IDs_z_mass[ri[x]]
				sffile.write(str(sf_ID_value)+' '+str(full_sf_redshifts_z_mass[ri[x]])+' '+str(full_sf_logmass_z_mass[ri[x]])+' ')
				for j in range(0, number_filters):
					flux_1D_array = sf_filt_flux_z_mass[:][j]
					sffile.write(str(flux_1D_array[ri[x]])+' ')
				sffile.write(' \n')
				if (randomize_IDs == 1):
					sf_ID_value = sf_ID_value + 1
			
		# Find the indices for the objects in the mass bin being examined. 
		mass_indices_q = np.where((full_q_logmass[z_indices_q] > massbin_min) & (full_q_logmass[z_indices_q] <= massbin_max))[0]
		full_q_IDs_z_mass = full_q_IDs[z_indices_q][mass_indices_q]
		full_q_redshifts_z_mass = full_q_redshifts[z_indices_q][mass_indices_q]
		full_q_logmass_z_mass = full_q_logmass[z_indices_q][mass_indices_q]

		# Find the Q flux values for the objects in the mass and redshift bin being examined.
		q_filt_flux_z_mass = np.zeros([number_filters, len(full_q_IDs_z_mass)])
		for j in range(0, number_filters):
			flux_1D_array = q_filt_flux[:][j]
			q_filt_flux_z_mass[:][j] = flux_1D_array[z_indices_q][mass_indices_q]
		
		# Now, choose some objects randomly from within the bin (providing there are enough objects)
		if (len(full_q_redshifts_z_mass) < number_per_mass_bin):
			ri = random.sample(xrange(0,len(full_q_redshifts_z_mass)), len(full_q_redshifts_z_mass))
				
		else:
			ri = random.sample(xrange(0,len(full_q_redshifts_z_mass)), number_per_mass_bin)

		# Now, go and find the ID, redshift, mass, etc, and write it to a file
		if (make_q == 1):
			print "    Q ri = "+str(len(ri))
			for x in range(0, len(ri)):
				if (randomize_IDs == 0):
					q_ID_value = full_q_IDs_z_mass[ri[x]]
				qfile.write(str(q_ID_value)+' '+str(full_q_redshifts_z_mass[ri[x]])+' '+str(full_q_logmass_z_mass[ri[x]])+' ')
				for j in range(0, number_filters):
					flux_1D_array = q_filt_flux_z_mass[:][j]
					qfile.write(str(flux_1D_array[ri[x]])+' ')
				qfile.write(' \n')
				if (randomize_IDs == 1):
					q_ID_value = q_ID_value + 1
		
if (make_sf == 1):
	sffile.close()
if (make_q == 1):
	qfile.close()
