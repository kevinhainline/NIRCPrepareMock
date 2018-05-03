import os
import sys
import math
import astropy
import numpy as np
import argparse
from astropy.io import fits
from astropy.table import Table

filter_file_name = 'NoisySubsample_to_PhotoZInput_filters.dat'

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

######################
# Optional Arguments #
######################

# Output Filename 
parser.add_argument(
  '-out','--output_name',
  help="Output Filename Extension?",
  action="store",
  type=str,
  dest="output_name",
  required=False
)

# Generate BEAGLE input file
parser.add_argument(
  '-beagle','--beagle_input_file',
  help="Make BEAGLE file?",
  action="store_true",
  dest="beagle_input_file",
  required=False
)

# Generate EAZY input file
parser.add_argument(
  '-eazy','--eazy_input_file',
  help="Make EAZY file?",
  action="store_true",
  dest="eazy_input_file",
  required=False
)

# Generate BPZ input file
parser.add_argument(
  '-bpz','--bpz_input_file',
  help="Make BPZ file?",
  action="store_true",
  dest="bpz_input_file",
  required=False
)

# Generate Le Phare input file
parser.add_argument(
  '-lep','--lephare_input_file','--lephare',
  help="Make Le Phare file?",
  action="store_true",
  dest="lephare_input_file",
  required=False
)

# Generate Zebra input file
parser.add_argument(
  '-zeb','--zebra_input_file','--zebra',
  help="Make Zebra file?",
  action="store_true",
  dest="zebra_input_file",
  required=False
)

args=parser.parse_args()

if (args.output_name):
	BEAGLE_output_filename = 'BEAGLE_input_'+args.output_name+'.fits'
	EAZY_output_filename = 'EAZY_input_'+args.output_name+'.dat'
	bpz_output_filename = 'BPZ_input_'+args.output_name+'.cat'
	bpz_output_columns_filename = 'BPZ_input_'+args.output_name+'.columns'
	lephare_output_filename = 'LEPHARE_input_'+args.output_name+'.in'
	zebra_output_filename = 'ZEBRA_input_'+args.output_name+'.cat'
	zebra_output_definitions_filename = 'ZEBRA_input_'+args.output_name+'.cat.def'
else:
	filename_extension = ''
	BEAGLE_output_filename = 'BEAGLE_input.fits'
	EAZY_output_filename = 'EAZY_input.dat'
	bpz_output_filename = 'BPZ_input.cat'
	bpz_output_columns_filename = 'BPZ_input.columns'
	lephare_output_filename = 'LEPHARE_input.in'
	zebra_output_filename = 'ZEBRA_input.cat'
	zebra_output_definitions_filename = 'ZEBRA_input.cat.def'


# Convert the fluxes (in erg/s/cm^2/Hz) and errors to magnitudes and magnitude errors for BPZ
def convert_fluxes_to_mags_bpz(flux_value,flux_error):

	c = 2.998e+18
	if (flux_value < 0) or (flux_value / flux_error < 1.0):
		mag = 99.0
		magerror = (-5.0/2.0) * math.log10(flux_error) - 48.6
	else:
		mag = (-5.0/2.0) * math.log10(flux_value) - 48.6
		magerror = (5.0/2.0) * (0.434) * abs(flux_error / flux_value)	
		
	return mag, magerror

# Convert the fluxes (in erg/s/cm^2/Hz) and errors to magnitudes and magnitude errors for ZEBRA
def convert_fluxes_to_mags_zebra(flux_value,flux_error):
	
	c = 2.998e+18
	if (flux_value < 0) or (flux_value / flux_error < 1.0):
		mag = 99.0
		magerror = (-5.0/2.0) * math.log10(flux_error) - 48.6
	else:
		mag = (-5.0/2.0) * math.log10(flux_value) - 48.6
		magerror = (5.0/2.0) * (0.434) * abs(flux_error / flux_value)	
		
	return mag, magerror


# OPEN THE FILTERS FILE
filters_full = np.loadtxt(filter_file_name, dtype='str')
header_filters = filters_full[:,0]
filters = filters_full[:,1]
bpz_filters_files = filters_full[:,2]
zebra_filters_files = filters_full[:,3]
number_filters = len(filters)

# OPEN THE INPUT DATA FILE
if (args.input_file.endswith('.fits')):
	fitsinput = fits.open(args.input_file)
	ID_values = fitsinput[1].data['ID']
	redshifts = fitsinput[1].data['redshift']
	number_objects = ID_values.size
	fluxes_all = np.zeros([number_objects, number_filters])
	errors_all = np.zeros([number_objects, number_filters])

	for j in range(0, number_filters):
		fluxes_all[:,j] = fitsinput[1].data[header_filters[j]]
		errors_all[:,j] = fitsinput[1].data[header_filters[j]+'_err']

else:
	full_input_file = np.loadtxt(args.input_file)
	number_objects = len(full_input_file)
	
	ID_values = np.zeros(number_objects)
	redshifts = np.zeros(number_objects)
	fluxes_all = np.zeros([number_objects,number_filters])
	errors_all = np.zeros([number_objects,number_filters])
	
	for n in range(0, number_objects):
		ID_values[n] = full_input_file[n][0]
		redshifts[n] = full_input_file[n][1]
		fluxes_all[n,:] = full_input_file[n][::2][1:number_filters+1]
		errors_all[n,:] = full_input_file[n][1::2][1:number_filters+1]
	
# Convert the fluxes and errors to ergs/s/cm^2/Hz, since they're in nJy (I assume)
fluxes_all_ergs = fluxes_all * 1e-23 * 1e-9
errors_all_ergs = errors_all * 1e-23 * 1e-9

if (args.beagle_input_file):
	print "Making BEAGLE file!"

	# ID redshift F070W F070W_err F090W F090W_err F115W F115W_err F140M F140M_err F150W F150W_err F162M F162M_err F182M F182M_err F200W F200W_err F210M F210M_err F250M F250M_err F277W F277W_err F300M F300M_err F335M F335M_err F356W F356W_err F360M F360M_err F410M F410M_err F430M F430M_err F444W F444W_err F460M F460M_err F480M F480M_err 
	# I d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d d

	# First, let's make the dtype and colnames arrays
	colnames = np.zeros(number_filters*2 + 2, dtype ='S10')
	dtype = np.zeros(number_filters*2 + 2, dtype ='str')
	colnames[0] = 'ID'
	colnames[1] = 'redshift'
	dtype[0] = 'I'
	dtype[1] = 'd'
	filter_number = 0
	for n in range(2, number_filters*2 + 2):
		if (n%2 == 1):
			colnames[n] = filters[filter_number]+'_err'
			filter_number = filter_number + 1
		else: 
			colnames[n] = filters[filter_number]
		dtype[n] = 'd'

	# And now let's assemble the data array
	beagle_data = np.zeros([number_objects, number_filters*2 + 2])
	beagle_data[:,0] = ID_values
	beagle_data[:,1] = redshifts
	column_number = 2
	for n in range(0, number_filters):
		beagle_data[:,column_number] = fluxes_all[:,n]
		column_number = column_number+1
		beagle_data[:,column_number] = errors_all[:,n]
		column_number = column_number+1
	
	# And finally, let's write out the output file.
	outtab = Table(beagle_data, names=colnames, dtype=dtype)
	outtab.write(BEAGLE_output_filename)

if (args.eazy_input_file):
	print "Making EAZY file!"

	# # id z_spec f_F435W e_F435W f_F606W e_F606W f_F775W e_F775W f_F814W e_F814W f_F850LP e_F850LP f_F090W e_F090W f_F115W e_F115W f_F150W e_F150W f_F200W e_F200W f_F277W e_F277W f_F335M e_F335M f_F356W e_F356W f_fake e_fake f_F410M e_F410M f_F444W e_F444W
	# # id z_spec F1 E1 F2 E2 F3 E3 F4 E4 F5 E5 F6 E6 F7 E7 F8 E8 F9 E9 F10 E10 F11 E11 F12 E12 F13 E13 F14 E14 F15 E15
	#   1.0 0.200166 3.00182640838e-30 2.40452886923e-32 8.12124135018e-30 1.51715515006e-32 1.18825750443e-29 3.16978638492e-32 1.23124056011e-29 2.63651347711e-32 1.35826456148e-29 2.89087954149e-32 1.35171226454e-29 5.70164272281e-33 1.48656227646e-29 4.20726628384e-33 1.44055231374e-29 4.0926065973e-33 1.43908346615e-29 4.1686938347e-33 9.5075489726e-30 5.80764417521e-33 7.12009252788e-30 1.05681750921e-32 6.6450080072e-30 6.25172692776e-33 5.06226663417e-30 6.13762005165e-33 5.374062576e-30 9.20449571753e-33 4.83285627801e-30 8.16582371359e-33 

	# Write the first two lines of the file
	f = open(EAZY_output_filename, 'w')
	f.write('# id z_spec ')
	for n in range(0, number_filters):
		f.write('f_'+filters[n]+' e_'+filters[n]+' ')
	f.write(' \n')

	f.write('# id z_spec ')
	for n in range(0, number_filters):
		f.write('F'+str((n+1))+' E'+str((n+1))+' ')
	f.write(' \n')

	# Write out the data 
	for j in range(0, number_objects):
		f.write(str(int(ID_values[j]))+' '+str(redshifts[j])+' ')
		for n in range(0, number_filters):
			f.write(str(fluxes_all_ergs[j,n])+' '+str(errors_all_ergs[j,n])+' ')
		f.write(' \n')

	f.close() 

if (args.bpz_input_file):
	print "Making BPZ files!"

	m_0_column = 0

	# Convert the fluxes and errors to mags and mag errors
	bpz_mags_all = np.zeros([number_objects,number_filters])
	bpz_mag_errors_all = np.zeros([number_objects,number_filters])
	
	for j in range(0, number_objects):
		for n in range(0, number_filters):
			bpz_mags_all[j,n], bpz_mag_errors_all[j,n] = convert_fluxes_to_mags_bpz(fluxes_all_ergs[j,n],errors_all_ergs[j,n])

	# 1 99.000 30.447  30.042 0.471   99.000 30.147  29.814 0.664   29.440 0.516    30.066 0.186   30.115 0.188   30.334 0.247  29.986 0.151   29.845 0.178  29.698 0.150  30.046 0.200  29.978 0.285  30.597 0.447  4.0010 

	f = open(bpz_output_filename, 'w')
	f.write('#  1 id \n')
	column_number = 2
	for n in range(0, number_filters):
		f.write('#  '+str(column_number)+' '+filters[n]+' \n')
		column_number = column_number+1
		f.write('#  '+str(column_number)+' d'+filters[n]+' \n')
		column_number = column_number+1
	f.write('#  '+str(column_number)+' zspec \n')
	f.write('# \n')
	
	f.write('# id ')
	for n in range(0, number_filters):
		f.write('f_'+filters[n]+' e_'+filters[n]+' ')
	f.write('zspec \n')

	for j in range(0, number_objects):
		f.write(str(int(ID_values[j]))+' ')
		for n in range(0, number_filters):
			f.write(str(round(bpz_mags_all[j,n],3))+' '+str(round(bpz_mag_errors_all[j,n],3))+'  ')
		f.write(str(redshifts[j])+' \n')
	
	f.close() 

	f = open(bpz_output_columns_filename, 'w')
	f.write('# Filter            columns  AB/Vega  zp_error  zp_offset \n')

	column_numbers = 2
	for n in range(0, number_filters):
		if (filters[n] == 'F200W'):
			m_0_column = column_numbers
		f.write(bpz_filters_files[n]+'    '+str(column_numbers)+', '+str(column_numbers+1)+' AB     0.01     0.0 \n')
		column_numbers = column_numbers + 2
	f.write('Z_S             '+str(column_numbers)+'\n')
	f.write('M_0             '+str(m_0_column)+'\n')
	f.write('ID              1 \n')
	
	f.close()

if (args.lephare_input_file):
	print "Making LEPHARE files!"

	f = open(lephare_output_filename, 'w')
	for j in range(0, number_objects):
		f.write(str(int(ID_values[j]))+' ')
		for n in range(0, number_filters):
			f.write(str(fluxes_all_ergs[j,n])+' '+str(errors_all_ergs[j,n])+'  ')
		f.write(' \n')
	f.close()

if (args.zebra_input_file):
	print "Making ZEBRA files!"

	# Convert the fluxes and errors to mags and mag errors
	zebra_mags_all = np.zeros([number_objects,number_filters])
	zebra_mag_errors_all = np.zeros([number_objects,number_filters])
	
	for j in range(0, number_objects):
		for n in range(0, number_filters):
			zebra_mags_all[j,n], zebra_mag_errors_all[j,n] = convert_fluxes_to_mags_zebra(fluxes_all_ergs[j,n],errors_all_ergs[j,n])

	# Write the first line of the file
	f = open(zebra_output_filename, 'w')
	f.write('# ')
	for n in range(0, number_filters):
		f.write('f_'+filters[n]+' e_'+filters[n]+' ')
	f.write('z_spec id')
	f.write(' \n')

	# Write out the data 
	for j in range(0, number_objects):
		for n in range(0, number_filters):
			f.write(str(round(zebra_mags_all[j,n],3))+' '+str(round(zebra_mag_errors_all[j,n],3))+'  ')
		f.write(str(redshifts[j])+' '+str(int(ID_values[j]))+' ')
		f.write(' \n')

	f.close() 
	
	# And now the configurations file
	f = open(zebra_output_definitions_filename, 'w')
	f.write('# column_name      type(F,dF,Z)     column_number \n')
	flux_column_number = 1
	for n in range(0, number_filters):
		f.write(bpz_filters_files[n]+'     F      '+str(flux_column_number)+'\n')
		flux_column_number = flux_column_number + 2
	flux_error_column_number = 2
	for n in range(0, number_filters):
		f.write(bpz_filters_files[n]+'     dF      '+str(flux_error_column_number)+'\n')
		flux_error_column_number = flux_error_column_number + 2
	
	f.write('z_cat     Z      '+str(flux_error_column_number-1)+'\n')
	f.close()
