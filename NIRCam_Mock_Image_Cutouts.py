import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import make_lupton_rgb
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

JAGUAR_version = 'r1_v1.1'

# Parser for figuring out the various input parameters
parser = argparse.ArgumentParser()

# Filter file names for R, G, and B filters
r_name = 'goods_s_F200W_2018_08_29.fits'
g_name = 'goods_s_F115W_2018_08_29.fits'
b_name = 'goods_s_F090W_2018_08_29.fits'

# Lupton image stretch values, you should change
# these if you're not super happy with the way the
# cutout looks. 
# See: http://adsabs.harvard.edu/abs/2004PASP..116..133L
# for more details
Qval = 8 
stretchval = 0.03
minimum = 0.0001

######################
# Required Arguments #
######################

# Input folder
parser.add_argument(
  '-in','--inputfolder',
  help="Input Folder With Mock Catalog Files (and images)",
  action="store",
  type=str,
  dest="input_folder",
  required=True
)

######################
# Optional Arguments #
######################

# Image Folder
parser.add_argument(
  '-image','--imagefolder',
  help="(Optional) Alternate Image Folder",
  action="store",
  type=str,
  dest="image_folder",
  required=True
)

# Object ID
parser.add_argument(
  '-objid','--objid',
  help="Object ID",
  action="store",
  type=str,
  dest="objid",
  required=False
)

# Center RA 
parser.add_argument(
  '-ra','--ra',
  help="Center RA of Box",
  action="store",
  type=str,
  dest="center_ra",
  required=False
)

# Center DEC
parser.add_argument(
  '-dec','--dec',
  help="Center DEC of Box",
  action="store",
  type=str,
  dest="center_dec",
  required=False
)

# Height of Cutout
parser.add_argument(
  '-rasize','--rasize',
  help="Width of Cutout in Arcseconds (Default = 5)",
  action="store",
  type=str,
  dest="rasize",
  required=False
)

# Width of Cutout
parser.add_argument(
  '-decsize','--decsize',
  help="Height of Cutout in Arcseconds (Default = 5)",
  action="store",
  type=str,
  dest="decsize",
  required=False
)

# RGB Scaling
parser.add_argument(
  '-rgb','--rgblist',
  help="RGB Color Scaling List, e.g. 1,1,1",
  action="store",
  type=str,
  dest="rgblist",
  required=False
)
# Show crosshairs on the image? 
parser.add_argument(
  '-cross','--cross',
  help="Show crosshairs?",
  action="store_true",
  dest="show_crosshairs",
  required=False
)

# Show scale on the image? 
#parser.add_argument(
#  '-scale','--scale',
#  help="Show scale bar?",
#  action="store_true",
#  dest="show_scale",
#  required=False
#)

# Show info on image?
parser.add_argument(
  '-info','--info',
  help="Show info on image?",
  action="store_true",
  dest="show_info",
  required=False
)

# Show info on image?
parser.add_argument(
  '-filters','--filters',
  help="Show filter names on image?",
  action="store_true",
  dest="show_filters",
  required=False
)

args=parser.parse_args()

# Folder with the JAGUAR files
jaguar_folder = args.input_folder

# Image folder and names. 
if (args.image_folder):
	image_folder = args.image_folder
else:
	image_folder = args.input_folder

r_name = image_folder+r_name
g_name = image_folder+g_name
b_name = image_folder+b_name

if (args.rasize):
	rasize_value = np.float(args.rasize)
else:
	rasize_value = 5.0
if (args.decsize):
	decsize_value = np.float(args.decsize)
else:
	decsize_value = 5.0

# Full star-forming catalogue, fits file
full_sf_mock_file = jaguar_folder+'JADES_SF_mock_'+JAGUAR_version+'.fits'
sftestfits = fits.open(full_sf_mock_file)
full_sf_IDs = sftestfits[1].data['ID']
full_sf_RA = sftestfits[1].data['RA']
full_sf_DEC = sftestfits[1].data['DEC']
full_sf_redshifts = sftestfits[1].data['redshift']
n_sf_objects = len(full_sf_redshifts)
max_SF_ID = np.max(full_sf_IDs)

# Full star-forming catalogue, fits file
full_q_mock_file = jaguar_folder+'JADES_Q_mock_'+JAGUAR_version+'.fits'
qtestfits = fits.open(full_q_mock_file)
full_q_IDs = qtestfits[1].data['ID']
full_q_RA = qtestfits[1].data['RA']
full_q_DEC = qtestfits[1].data['DEC']
full_q_redshifts = qtestfits[1].data['redshift']
n_q_objects = len(full_q_redshifts)

# If the user wants to make a cutout of a specific object
if (args.objid):
	object_ID = np.int(args.objid)
	
	if (object_ID > max_SF_ID):
		#print "Q Object!"
		ID_index = np.where(full_q_IDs == object_ID)[0]
		objRA = full_q_RA[ID_index][0]
		objDEC = full_q_DEC[ID_index][0]
		objredshift = full_q_redshifts[ID_index][0]
	else:
		#print "SF Object!"
		ID_index = np.where(full_sf_IDs == object_ID)[0]
		objRA = full_sf_RA[ID_index][0]
		objDEC = full_sf_DEC[ID_index][0]
		objredshift = full_sf_redshifts[ID_index][0]
else:
	if (args.center_ra is None):
		sys.exit("You must specify Object ID or RA/DEC!")
	if (args.center_dec is None):
		sys.exit("You must specify Object ID or RA/DEC!")
	if (args.center_ra):
		objRA = np.float(args.center_ra)
		objDEC = np.float(args.center_dec)

cosdec_center = math.cos(objDEC * 3.141593 / 180.0)


# Open the R band image
r = fits.open(r_name)[0].data
rhdu = fits.open(r_name)[0]
rwcs = WCS(rhdu.header)

# Open the g band image
g = fits.open(g_name)[0].data
ghdu = fits.open(b_name)[0]
gwcs = WCS(ghdu.header)

# Open the b band image
b = fits.open(b_name)[0].data
bhdu = fits.open(g_name)[0]
bwcs = WCS(bhdu.header)

if (args.rgblist):
	rgb_list = [np.float(item) for item in args.rgblist.split(',')]
	print "Scaling by: r:"+str(rgb_list[0])+", g:"+str(rgb_list[1])+", b:"+str(rgb_list[2])
	r = r*rgb_list[0]
	g = g*rgb_list[1]
	b = b*rgb_list[2]

# Set the position of the object
position = SkyCoord(str(objRA)+'d '+str(objDEC)+'d', frame='fk5')

# Make the cutout
size = u.Quantity((decsize_value, rasize_value), u.arcsec)
rcutout = Cutout2D(r, position, size, wcs=rwcs)
gcutout = Cutout2D(g, position, size, wcs=gwcs)
bcutout = Cutout2D(b, position, size, wcs=bwcs)

# Use the method outlined in Lupton et al. (2004) to generate the image
# http://docs.astropy.org/en/stable/visualization/lupton_rgb.html
rgbcutout = make_lupton_rgb(rcutout.data, gcutout.data, bcutout.data, Q=Qval, stretch=stretchval, minimum = minimum)

fig = plt.figure()
ax = plt.subplot(projection=rcutout.wcs)
plt.imshow(rgbcutout, origin = 'lower')

# Add a crosshair
if (args.show_crosshairs):
	closeness_to_object = 0.05
	farposition_from_object = 0.12
	ax.plot([objRA, objRA],[objDEC + (closeness_to_object)*decsize_value/3600.0 , objDEC + (farposition_from_object)*decsize_value/3600.0], transform=ax.get_transform('fk5'), color='white')
	ax.plot([objRA, objRA],[objDEC - (closeness_to_object)*decsize_value/3600.0 , objDEC - (farposition_from_object)*decsize_value/3600.0], transform=ax.get_transform('fk5'), color='white')
	ax.plot([objRA + (closeness_to_object) * rasize_value/3600.0, objRA + (farposition_from_object) * rasize_value/3600.0],[objDEC, objDEC], transform=ax.get_transform('fk5'), color='white')
	ax.plot([objRA - (closeness_to_object) * rasize_value/3600.0, objRA - (farposition_from_object) * rasize_value/3600.0],[objDEC, objDEC], transform=ax.get_transform('fk5'), color='white')

# Add the scale bar 
#if (rasize_value < 8):
#	scale_bar_length = 1.0
#if (rasize_value >= 8):
#	scale_bar_length = 3.0
#if (rasize_value >= 60):
#	scale_bar_length = 30.0
#if (args.show_scale):
#	leftside_RA = (objRA+scale_bar_length/(2*3600)/cosdec_center)+(0.4)*(rasize_value/3600.0)
#	rightside_RA = (objRA-scale_bar_length/(2*3600)/cosdec_center)+(0.4)*(rasize_value/3600.0)
#	bar_average_RA = (leftside_RA + rightside_RA)/2.0
#	ax.plot([leftside_RA, rightside_RA],[objDEC + (0.4)*decsize_value/3600.0, objDEC + (0.4)*decsize_value/3600.0], transform=ax.get_transform('fk5'), color='white')
#	ax.text(bar_average_RA, objDEC + (0.43)*decsize_value/3600.0,str(scale_bar_length)+'\'\'', transform=ax.get_transform('fk5'), color='white', horizontalalignment='center')

if (args.objid):

	# Add the Information about the object
	if (args.show_info):
		rightside_RA = (objRA-scale_bar_length/(2*3600))-(0.3)*rasize_value/3600.0*cosdec_center
		ax.text(rightside_RA, objDEC + (0.43)*decsize_value/3600.0,'Obj '+str(object_ID), transform=ax.get_transform('fk5'), color='white', horizontalalignment='right')
		ax.text(rightside_RA, objDEC + (0.38)*decsize_value/3600.0,'z = '+str(objredshift), transform=ax.get_transform('fk5'), color='white', horizontalalignment='right')
	
# Add the filters in the RGB colors they correspond to
if (args.show_filters):
	rightside_RA = (objRA-scale_bar_length/(2*3600))-(0.3)*rasize_value/3600.0*cosdec_center
	ax.text(rightside_RA, objDEC - (0.43)*decsize_value/3600.0, r_filtername, transform=ax.get_transform('fk5'), color = 'red', horizontalalignment='right')
	ax.text(rightside_RA, objDEC - (0.38)*decsize_value/3600.0, g_filtername, transform=ax.get_transform('fk5'), color = 'green', horizontalalignment='right')
	ax.text(rightside_RA, objDEC - (0.33)*decsize_value/3600.0, b_filtername, transform=ax.get_transform('fk5'), color = 'blue', horizontalalignment='right')

# Now turn off the axis labeling
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_ticks_visible(False)
lon.set_ticklabel_visible(False)
lat.set_ticks_visible(False)
lat.set_ticklabel_visible(False)
lon.set_axislabel('')
lat.set_axislabel('')

# And save the figure
if (args.objid):
	plt.savefig(str(object_ID)+'_cutout_'+str(rasize_value)+'-'+str(decsize_value)+'.png', bbox_inches='tight', dpi = 300)
if (args.center_ra):
	plt.savefig(str(objRA)+'_'+str(objDEC)+'_cutout_'+str(rasize_value)+'-'+str(decsize_value)+'.png', bbox_inches='tight', dpi = 300)
