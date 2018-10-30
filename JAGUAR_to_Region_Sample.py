import os
import sys
import math
import argparse
import random
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
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

# Width of Subsample
parser.add_argument(
  '-width',
  help="Width of subsample (in arcminutes, max 11')",
  action="store",
  type=float,
  dest="width",
  required=True
)

# Height of Subsample
parser.add_argument(
  '-height',
  help="Height of subsample (in arcminutes, max 11')",
  action="store",
  type=float,
  dest="height",
  required=True
)

######################
# Optional Arguments #
######################


args=parser.parse_args()

# Set survey size 
# width of region, in arcminutes
width = args.width #6
height = args.height #12

total_jaguar_width = 11.0
total_jaguar_height = 11.0

ra_center = 53.1625
dec_center = -27.78141667
cosdec_center = math.cos(dec_center * 3.141593 / 180.0)

if ((width > total_jaguar_width) | (height > total_jaguar_height)):
	sys.exit("Width or height is larger than the 11x11 JAGUAR survey size!")
else:
	new_mock_center_ra = np.random.uniform(ra_center - (((total_jaguar_width - width)/2.0)/60.0 / cosdec_center), ra_center + (((total_jaguar_width - width)/2.0)/60.0 / cosdec_center))
	new_mock_center_dec = np.random.uniform(dec_center - ((total_jaguar_height - height)/2.0)/60.0, dec_center + ((total_jaguar_height - height)/2.0)/60.0)
	mock_positions_ra_max = new_mock_center_ra + ((width/2.0)/60.0 / cosdec_center)
	mock_positions_ra_min = new_mock_center_ra - ((width/2.0)/60.0 / cosdec_center)
	mock_positions_dec_max = new_mock_center_dec + ((height/2.0)/60.0)
	mock_positions_dec_min = new_mock_center_dec - ((height/2.0)/60.0)

# Full star-forming catalogue, fits file
full_sf_mock_file = args.input_folder+'JADES_SF_mock_'+JAGUAR_version+'.fits'
sftestfits = fits.open(full_sf_mock_file)

# Get redshifts and galaxy properties
full_sf_IDs = sftestfits[1].data['ID']
full_sf_RAs = sftestfits[1].data['RA']
full_sf_DECs = sftestfits[1].data['DEC']
n_sf_objects = len(full_sf_IDs)

# Full quiescent catalogue, fits file
full_q_mock_file = args.input_folder+'JADES_Q_mock_'+JAGUAR_version+'.fits'
qtestfits = fits.open(full_q_mock_file)

# Get redshifts and galaxy properties
full_q_IDs = qtestfits[1].data['ID']
full_q_RAs = qtestfits[1].data['RA']
full_q_DECs = qtestfits[1].data['DEC']
n_q_objects = len(full_q_IDs)

for j in range(0, n_sf_objects):
	if ((full_sf_RAs[j] >= mock_positions_ra_min) & (full_sf_RAs[j] <= mock_positions_ra_max) & (full_sf_DECs[j] >= mock_positions_dec_min) & (full_sf_DECs[j] <= mock_positions_dec_max)):
		print full_sf_IDs[j]#, full_sf_RAs[j], full_sf_DECs[j]

for j in range(0, n_q_objects):
	if ((full_q_RAs[j] >= mock_positions_ra_min) & (full_q_RAs[j] <= mock_positions_ra_max) & (full_q_DECs[j] >= mock_positions_dec_min) & (full_q_DECs[j] <= mock_positions_dec_max)):
		print full_q_IDs[j]#, full_q_RAs[j], full_q_DECs[j]

