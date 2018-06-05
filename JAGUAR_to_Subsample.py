import os
import sys
import argparse
import random
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table

def JytoABMag(flux):
	return (-5.0 / 2.0) * np.log10(flux) - 48.60

JAGUAR_version = 'r1_v1.1'


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

# Input ID Filename
parser.add_argument(
  '-iID','--inputIDs',
  help="Input ID Filename",
  action="store",
  type=str,
  dest="inputIDs",
  required=False
)

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

# Combine Files
parser.add_argument(
  '-co','--coutputfilename',
  help="Filename for combined file?",
  action="store",
  type=str,
  dest="combine_filename",
  required=False
)

# Generate Fits File
parser.add_argument(
  '-mf','--make_fits',
  help="Make Fits File?",
  action="store_true",
  dest="make_fits",
  required=False
)

args=parser.parse_args()

if (args.sf_output_filename):
	make_sf = 1
	if ((args.sf_output_filename.endswith('.txt')) or (args.sf_output_filename.endswith('.dat'))):
		sf_filenameroot = ('.').join(args.sf_output_filename.split('.')[:-1])
	else:
		sf_filenameroot = args.sf_output_filename
else:
	make_sf = 0
	
if (args.q_output_filename):
	make_q = 1
	if ((args.q_output_filename.endswith('.txt')) or (args.q_output_filename.endswith('.dat'))):
		q_filenameroot = ('.').join(args.q_output_filename.split('.')[:-1])
	else:
		q_filenameroot = args.q_output_filename
else:
	make_q = 0

if (args.combine_filename):
	make_combine = 1
	if ((args.combine_filename.endswith('.txt')) or (args.combine_filename.endswith('.dat'))):
		c_filenameroot = ('.').join(args.combine_filename.split('.')[:-1])
	else:
		c_filenameroot = args.combine_filename
else:
	make_combine = 0


if (((make_sf == 0) or (make_q == 0)) and (make_combine == 1)):
	sys.exit("I can't combine the SF and Q files, since you only request one of the two. ")
 
if (args.make_fits):
	make_fits_file = 1
else: 
	make_fits_file = 0
	
#if ((make_sf == 0) & (make_q == 0)):
#	sys.exit("No output SF or Q filename specified! Use -sfo or -q to provide a filename")

if (args.inputIDs):
	id_file = np.loadtxt(args.inputIDs, dtype = 'int')
	input_ID_numbers = id_file[:]
	n_input_ID_numbers = len(input_ID_numbers)
	
	id_filenameroot = 'ID_output_list'
	use_IDs = 1
	make_sf = 0
	make_q = 0
	make_combine = 0

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


redshift_minimum = 0.0
redshift_maximum = 15.0
redshift_range = np.arange(redshift_minimum, redshift_maximum+1.0, 1.0)
n_redshift_bins = len(redshift_range - 1)
if (args.number_per_redshift_bin):
	number_per_redshift_bin = int(args.number_per_redshift_bin)#1000
else:
	number_per_redshift_bin = 1000
	
mass_minimum = 6.0
mass_maximum = 12.0
mass_range = np.arange(mass_minimum, mass_maximum+1.0, 1.0)
n_mass_bins = len(mass_range)-1
number_per_mass_bin = number_per_redshift_bin / (n_mass_bins)

#print number_per_mass_bin

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
sf_filt_flux = np.zeros([number_filters, n_sf_objects])
for j in range(0, number_filters):
	sf_filt_flux[:][j] = sftestfits[1].data[filters[j]+'_fnu']#*1e-23*1e-9

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
q_filt_flux = np.zeros([number_filters, n_q_objects])
for j in range(0, number_filters):
	q_filt_flux[:][j] = qtestfits[1].data[filters[j]+'_fnu']#*1e-23*1e-9

# Open up the output files. 
if (make_sf == 1):
	sffile = open(sf_filenameroot+'.dat', 'w')
	sffile.write('#ID    Redshift    Log(Mass)    Re_Maj    Sersic_n     ')
	for j in range(0, number_filters):
		sffile.write(filters[j]+'    ')
	sffile.write(' \n')
if (make_q == 1):
	qfile = open(q_filenameroot+'.dat', 'w')
	qfile.write('#ID    Redshift    Log(Mass)    Re_Maj    Sersic_n     ')
	for j in range(0, number_filters):
		qfile.write(filters[j]+'    ')
	qfile.write(' \n')
if (use_IDs == 1):
	idfile = open(id_filenameroot+'.dat', 'w')
	idfile.write('#ID    Redshift    Log(Mass)    Re_Maj    Sersic_n     ')
	for j in range(0, number_filters):
		idfile.write(filters[j]+'    ')
	idfile.write(' \n')

	
# This is to be used if the randomized fluxes are requested
sf_ID_value = 1
q_ID_value = 1

if (use_IDs == 1):
	for j in range(0, n_input_ID_numbers):
		if (input_ID_numbers[j] <= np.max(full_sf_IDs)):
			object_ID = input_ID_numbers[j]
			object_index = np.where(full_sf_IDs == object_ID)[0][0]
			object_redshifts = full_sf_redshifts[object_index]
			object_logmass = full_sf_logmass[object_index]
			object_re_maj = full_sf_re_maj[object_index]
			object_sersic_n = full_sf_sersic_n[object_index]
			object_fluxes = np.zeros(number_filters)
			for j in range(0, number_filters):
				object_fluxes[j] = sf_filt_flux[object_index][j]
		else:
			object_ID = input_ID_numbers[j]
			object_index = np.where(full_q_IDs == object_ID)[0][0]
			object_redshifts = full_q_redshifts[object_index]
			object_logmass = full_q_logmass[object_index]
			object_re_maj = full_q_re_maj[object_index]
			object_sersic_n = full_q_sersic_n[object_index]
			for j in range(0, number_filters):
				object_fluxes[j] = q_filt_flux[object_index][j]

		idfile.write(str(object_ID)+' '+str(object_redshifts)+' '+str(object_logmass)
					+' '+str(object_re_maj)+' '+str(object_sersic_n)+' ')
		for j in range(0, number_filters):
			idfile.write(str(object_fluxes[j])+' ')
		idfile.write(' \n')
		
	idfile.close()

else:
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
			full_sf_re_maj_z_mass = full_sf_re_maj[z_indices_sf][mass_indices_sf]
			full_sf_sersic_n_z_mass = full_sf_sersic_n[z_indices_sf][mass_indices_sf]
	
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
					sffile.write(str(sf_ID_value)+' '+str(full_sf_redshifts_z_mass[ri[x]])+' '+str(full_sf_logmass_z_mass[ri[x]])
						+' '+str(full_sf_re_maj_z_mass[ri[x]])+' '+str(full_sf_sersic_n_z_mass[ri[x]])+' ')
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
			full_q_re_maj_z_mass = full_q_re_maj[z_indices_q][mass_indices_q]
			full_q_sersic_n_z_mass = full_q_sersic_n[z_indices_q][mass_indices_q]
	
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
					qfile.write(str(q_ID_value)+' '+str(full_q_redshifts_z_mass[ri[x]])+' '+str(full_q_logmass_z_mass[ri[x]])
						+' '+str(full_q_re_maj_z_mass[ri[x]])+' '+str(full_q_sersic_n_z_mass[ri[x]])+' ')
					for j in range(0, number_filters):
						flux_1D_array = q_filt_flux_z_mass[:][j]
						qfile.write(str(flux_1D_array[ri[x]])+' ')
					qfile.write(' \n')
					if (randomize_IDs == 1):
						q_ID_value = q_ID_value + 1
			
	if (make_sf == 1):
		sffile.close()
		os.system('sort -bn -k1 '+sf_filenameroot+'.dat > '+sf_filenameroot+'.sorted.dat')
		os.system('rm '+sf_filenameroot+'.dat')
		os.system('mv '+sf_filenameroot+'.sorted.dat '+sf_filenameroot+'.dat')
	if (make_q == 1):
		qfile.close()
		os.system('sort -bn -k1 '+q_filenameroot+'.dat > '+q_filenameroot+'.sorted.dat')
		os.system('rm '+q_filenameroot+'.dat')
		os.system('mv '+q_filenameroot+'.sorted.dat '+q_filenameroot+'.dat')
	if (args.combine_filename):
		os.system("sed '1d' "+q_filenameroot+".dat > tmp.dat")
		os.system('cat '+sf_filenameroot+'.dat tmp.dat > '+c_filenameroot+'.dat')
		os.system('rm tmp.dat')
	
if (args.make_fits):
	# First, let's make the  dtype and colnames arrays
	colnames = np.zeros(5+number_filters, dtype ='S20')
	dtype = np.zeros(5+number_filters, dtype ='str')
	colnames[0] = 'ID'
	colnames[1] = 'redshift'
	colnames[2] = 'logmass'
	colnames[3] = 're_major'
	colnames[4] = 'sersic_n'
		
	dtype[0] = 'I'
	dtype[1] = 'd'
	dtype[2] = 'd'
	dtype[3] = 'd'
	dtype[4] = 'd'
	for j in range(0, number_filters):
		colnames[j+5] = filters[j]
		dtype[j+5] = 'd'

	if (make_sf == 1):
		catalogue_file = sf_filenameroot+'.dat'
		cat_file_full = np.loadtxt(catalogue_file)
		ID_numbers = cat_file_full[:,0]
		redshifts = cat_file_full[:,1]
		logmasses = cat_file_full[:,2]
		re_major = cat_file_full[:,3]
		sersic_n = cat_file_full[:,4]
		n_objects = ID_numbers.size
	
		apparent_flux = np.zeros([number_filters, n_objects])
		for j in range(0, number_filters):
			apparent_flux[:][j] = cat_file_full[:,5+j]
		
		# And now let's assemble the data array
		output_data = np.zeros([n_objects, 5+number_filters])
		output_data[:,0] = ID_numbers
		output_data[:,1] = redshifts
		output_data[:,2] = logmasses
		output_data[:,3] = re_major
		output_data[:,4] = sersic_n
		for j in range(0, number_filters):
			output_data[:,j+5] = apparent_flux[:][j]
		
		# And finally, let's write out the output file.
		outtab = Table(output_data, names=colnames, dtype=dtype)
		outtab.write(sf_filenameroot+'.fits')

	if (make_q == 1):
		catalogue_file = q_filenameroot+'.dat'
		cat_file_full = np.loadtxt(catalogue_file)
		ID_numbers = cat_file_full[:,0]
		redshifts = cat_file_full[:,1]
		logmasses = cat_file_full[:,2]
		re_major = cat_file_full[:,3]
		sersic_n = cat_file_full[:,4]
		n_objects = ID_numbers.size
	
		apparent_flux = np.zeros([number_filters, n_objects])
		for j in range(0, number_filters):
			apparent_flux[:][j] = cat_file_full[:,5+j]
		
				
		# And now let's assemble the data array
		output_data = np.zeros([n_objects, 5+number_filters])
		output_data[:,0] = ID_numbers
		output_data[:,1] = redshifts
		output_data[:,2] = logmasses
		output_data[:,3] = re_major
		output_data[:,4] = sersic_n
		for j in range(0, number_filters):
			output_data[:,j+5] = apparent_flux[:][j]
		
		# And finally, let's write out the output file.
		outtab = Table(output_data, names=colnames, dtype=dtype)
		outtab.write(q_filenameroot+'.fits')

	if (make_combine == 1):
		catalogue_file = c_filenameroot+'.dat'
		cat_file_full = np.loadtxt(catalogue_file)
		ID_numbers = cat_file_full[:,0]
		redshifts = cat_file_full[:,1]
		logmasses = cat_file_full[:,2]
		re_major = cat_file_full[:,3]
		sersic_n = cat_file_full[:,4]
		n_objects = ID_numbers.size
	
		apparent_flux = np.zeros([number_filters, n_objects])
		for j in range(0, number_filters):
			apparent_flux[:][j] = cat_file_full[:,5+j]
		
		# And now let's assemble the data array
		output_data = np.zeros([n_objects, 5+number_filters])
		output_data[:,0] = ID_numbers
		output_data[:,1] = redshifts
		output_data[:,2] = logmasses
		output_data[:,3] = re_major
		output_data[:,4] = sersic_n
		for j in range(0, number_filters):
			output_data[:,j+5] = apparent_flux[:][j]
		
		# And finally, let's write out the output file.
		outtab = Table(output_data, names=colnames, dtype=dtype)
		outtab.write(c_filenameroot+'.fits')

	if (use_IDs == 1):
		catalogue_file = id_filenameroot+'.dat'
		cat_file_full = np.loadtxt(catalogue_file)
		ID_numbers = cat_file_full[:,0]
		redshifts = cat_file_full[:,1]
		logmasses = cat_file_full[:,2]
		re_major = cat_file_full[:,3]
		sersic_n = cat_file_full[:,4]
		n_objects = ID_numbers.size
	
		apparent_flux = np.zeros([number_filters, n_objects])
		for j in range(0, number_filters):
			apparent_flux[:][j] = cat_file_full[:,5+j]
		
		# And now let's assemble the data array
		output_data = np.zeros([n_objects, 5+number_filters])
		output_data[:,0] = ID_numbers
		output_data[:,1] = redshifts
		output_data[:,2] = logmasses
		output_data[:,3] = re_major
		output_data[:,4] = sersic_n
		for j in range(0, number_filters):
			output_data[:,j+5] = apparent_flux[:][j]
		
		# And finally, let's write out the output file.
		outtab = Table(output_data, names=colnames, dtype=dtype)
		outtab.write(id_filenameroot+'.fits')

