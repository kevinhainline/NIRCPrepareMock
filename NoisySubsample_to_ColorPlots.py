import os
import sys
import math
import astropy
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import argparse
from astropy.io import fits
from astropy.table import Table

noisy_jades_filters = ['HST_F435W', 'HST_F606W', 'HST_F775W', 'HST_F814W', 'HST_F850LP', 'NRC_F070W', 'NRC_F090W', 'NRC_F115W', 'NRC_F150W', 'NRC_F200W', 'NRC_F277W', 'NRC_F335M', 'NRC_F356W', 'NRC_F410M', 'NRC_F444W']
number_jades_filters = len(noisy_jades_filters)

def FluxtoABMag(flux):
	return (-5.0 / 2.0) * np.log10(flux) - 48.60

# Finds the intersection of two curves with the same xvalues
# crossover_value = find_intersection(selection_limit_highres, selection_fraction_highres, total_selection_fraction_highres)
def find_intersection(xvalues, yvalues1, yvalues2):
	yvalue_return = 0
	xvalue_return = 0
	for a in range(0, len(xvalues)-1):
		#for b in range(0, len(xvalues)-1):
			#if ((xvalues[a] > 1.24) & (xvalues[a] < 1.56)):
				#print xvalues[a], yvalues2[a], yvalues1[a], yvalues2[a+1]
			#if ((yvalues2[a+1] <= yvalues1[a]) & (yvalues2[a] >= yvalues1[a])):
			if ((yvalues2[a+1] <= yvalues1[a+1]) & (yvalues2[a] >= yvalues1[a])):
				#print xvalues[a], yvalues2[a], yvalues1[a], yvalues2[a+1]
				slope_line_one = (yvalues1[a+1] - yvalues1[a]) / (xvalues[a+1] - xvalues[a])
				y_intercept_line_one = yvalues1[a] - (slope_line_one * xvalues[a])
				slope_line_two = (yvalues2[a+1] - yvalues2[a]) / (xvalues[a+1] - xvalues[a])
				y_intercept_line_two = yvalues2[a] - (slope_line_two * xvalues[a])
				#print slope_line_one, y_intercept_line_one, slope_line_two, y_intercept_line_two
				xvalue_intersect = (-1.0)*(y_intercept_line_one - y_intercept_line_two)/(slope_line_one - slope_line_two)
				#print "X_value intersect = "+str(xvalue_intersect)
				xvalue_return = xvalue_intersect#xvalues[a]
				yvalue_return = (slope_line_one * xvalue_intersect) + y_intercept_line_one
	return xvalue_return, yvalue_return

	


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


# Y Axis, Filter 1
parser.add_argument(
  '-yf1','--yfilter1',
  help="Y-Axis Filter 1",
  action="store",
  type=str,
  dest="yfilter1",
  required=True
)

# Y Axis, Filter 2
parser.add_argument(
  '-yf2','--yfilter2',
  help="Y-Axis Filter 2",
  action="store",
  type=str,
  dest="yfilter2",
  required=True
)


# X Axis, Filter 1
parser.add_argument(
  '-xf1','--xfilter1',
  help="X-Axis Filter 1",
  action="store",
  type=str,
  dest="xfilter1",
  required=True
)

# X Axis, Filter 2
parser.add_argument(
  '-xf2','--xfilter2',
  help="X-Axis Filter 2",
  action="store",
  type=str,
  dest="xfilter2",
  required=True
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

######################
# Optional Arguments #
######################

# Optional ID list
parser.add_argument(
  '-idlist','--id_number_list',
  help="List of ID Numbers?",
  action="store",
  type=str,
  dest="id_number_list",
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

# X-Axis Limit 
parser.add_argument(
  '-xlim','--xlimit',
  help="Limit on X-Axis Color? ",
  action="store",
  type=float,
  dest="xlimit",
  required=False
)

# Positive X-Axis Limit 
parser.add_argument(
  '-xplim','--xplimit',
  help="Positive Limit on X-Axis Color? ",
  action="store",
  type=float,
  dest="xplimit",
  required=False
)

# Negative X-Axis Limit 
parser.add_argument(
  '-xnlim','--xnlimit',
  help="Negative Limit on X-Axis Color? ",
  action="store",
  type=float,
  dest="xnlimit",
  required=False
)

# Y Slope Value
parser.add_argument(
  '-yslope','--yslope',
  help="Slope on the Selection Line?",
  action="store",
  type=float,
  dest="yslope",
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

# Verbose? 
parser.add_argument(
  '-verb','--verbose',
  help="Verbose Mode?",
  action="store_true",
  dest="verbose",
  required=False
)



args=parser.parse_args()

if (args.id_number_list):
	# Get the IDs for the objects to plot. 
	ID_input_file = np.loadtxt(args.id_number_list)
	ID_numbers = ID_input_file[:,0]
	number_input_objects = len(ID_numbers)

# OPEN THE INPUT DATA FILE
if (args.input_file.endswith('.fits')):
	fitsinput = fits.open(args.input_file)
	ID_values = fitsinput[1].data['ID']
	redshifts = fitsinput[1].data['redshift']
	number_objects = ID_values.size



if (args.id_number_list):
	correct_object_indices = np.zeros(number_input_objects, dtype = int)
	for z in range(0, number_input_objects):
		correct_object_indices[z] = np.where(ID_values == ID_numbers[z])[0][0]
else:
	number_input_objects = len(ID_values)
	correct_object_indices = np.arange(0,number_input_objects)

y_filter_flux_one = fitsinput[1].data[args.yfilter1][correct_object_indices]
y_filter_flux_one_err = fitsinput[1].data[args.yfilter1+'_err'][correct_object_indices]
y_filter_flux_two = fitsinput[1].data[args.yfilter2][correct_object_indices]
y_filter_flux_two_err = fitsinput[1].data[args.yfilter2+'_err'][correct_object_indices]

x_filter_flux_one = fitsinput[1].data[args.xfilter1][correct_object_indices]
x_filter_flux_one_err = fitsinput[1].data[args.xfilter1+'_err'][correct_object_indices]
x_filter_flux_two = fitsinput[1].data[args.xfilter2][correct_object_indices]
x_filter_flux_two_err = fitsinput[1].data[args.xfilter2+'_err'][correct_object_indices]

# Have to do a SNR check. 
y_filter_flux_one_snr = y_filter_flux_one / y_filter_flux_one_err
y_filter_flux_two_snr = y_filter_flux_two / y_filter_flux_two_err
x_filter_flux_one_snr = x_filter_flux_one / x_filter_flux_one_err
x_filter_flux_two_snr = x_filter_flux_two / x_filter_flux_two_err

snr_limit = 5.0
SNR_check_objects = np.where((y_filter_flux_two_snr > snr_limit) & (x_filter_flux_one_snr > snr_limit) & (x_filter_flux_two_snr > snr_limit))[0]

#print np.sort(y_filter_flux_one_snr)
#print len(y_filter_flux_two_snr)
#print len(x_filter_flux_one_snr)
#print len(x_filter_flux_two_snr)


y_filter_one = FluxtoABMag(fitsinput[1].data[args.yfilter1][correct_object_indices][SNR_check_objects]*1e-23*1e-9)
y_filter_one[np.where(y_filter_flux_one[SNR_check_objects] < 0)[0]] = 99.00
#y_filter_one_err = fitsinput[1].data[args.yfilter1+'_err'][correct_object_indices]
y_filter_two = FluxtoABMag(fitsinput[1].data[args.yfilter2][correct_object_indices][SNR_check_objects]*1e-23*1e-9)
#y_filter_two[np.where(y_filter_flux_two < 0)[0]] = 99.00
#y_filter_two_err = fitsinput[1].data[args.yfilter2+'_err'][correct_object_indices]

x_filter_one = FluxtoABMag(fitsinput[1].data[args.xfilter1][correct_object_indices][SNR_check_objects]*1e-23*1e-9)
#x_filter_one[np.where(x_filter_flux_one < 0)[0]] = 99.00
#x_filter_one_err = fitsinput[1].data[args.xfilter1+'_err'][correct_object_indices]
x_filter_two = FluxtoABMag(fitsinput[1].data[args.xfilter2][correct_object_indices][SNR_check_objects]*1e-23*1e-9)
#x_filter_two[np.where(x_filter_flux_two < 0)[0]] = 99.00
#x_filter_two_err = fitsinput[1].data[args.xfilter2+'_err'][correct_object_indices]

correct_redshifts = redshifts[correct_object_indices][SNR_check_objects]

x_axis_color = x_filter_one - x_filter_two
y_axis_color = y_filter_one - y_filter_two

if (args.xlimit):
	x_axis_limit = args.xlimit
if (args.xplimit):
	x_axis_plimit = args.xplimit
if (args.xnlimit):
	x_axis_nlimit = args.xnlimit


# Find all high-z objects
high_z_objects_i = np.where(redshifts[correct_object_indices][SNR_check_objects] > args.zlimit)[0]
print 'There are '+str(len(high_z_objects_i))+' objects at z > '+str(args.zlimit)

# Selection Lines
selection_limit = np.arange(0.1, 2.8, 0.2)#np.arange(0.5, 3.2, 0.2)
n_selection_limit = len(selection_limit)
selection_fraction = np.zeros(n_selection_limit)
selection_fraction_err = np.zeros(n_selection_limit)
total_selection_fraction = np.zeros(n_selection_limit)
total_selection_fraction_err = np.zeros(n_selection_limit)
for z in range(0, n_selection_limit):
#for z in range(0, 1):
	if (args.verbose):
		print "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
	# Did you set an X color limit?
	if ((args.xplimit) and not (args.xnlimit)):
		# Did you set a slope value? 
		if (args.yslope):
			selected_objects_i = np.where((y_axis_color > ((args.yslope*x_axis_color) + selection_limit[z])) & (x_axis_color < x_axis_plimit))[0]
		else:
			selected_objects_i = np.where((y_axis_color > selection_limit[z]) & (x_axis_color < x_axis_plimit))[0]
	elif ((args.xnlimit) and not (args.xplimit)):
		# Did you set a slope value? 
		if (args.yslope):
			selected_objects_i = np.where((y_axis_color > ((args.yslope*x_axis_color) + selection_limit[z])) & (x_axis_color > x_axis_nlimit))[0]
		else:
			selected_objects_i = np.where((y_axis_color > selection_limit[z]) & (x_axis_color > x_axis_nlimit))[0]
	elif ((args.xplimit) and (args.xnlimit)):
		# Did you set a slope value? 
		if (args.yslope):
			selected_objects_i = np.where((y_axis_color > ((args.yslope*x_axis_color) + selection_limit[z])) & (x_axis_color < x_axis_plimit) & (x_axis_color > x_axis_nlimit))[0]
		else:
			selected_objects_i = np.where((y_axis_color > selection_limit[z]) & (x_axis_color < x_axis_nlimit) & (x_axis_color > x_axis_nlimit))[0]
	else:
		if (args.yslope):
			selected_objects_i = np.where((y_axis_color > ((args.yslope*x_axis_color) + selection_limit[z])))[0]
		else:
			selected_objects_i = np.where((y_axis_color > selection_limit[z]))[0]
#	if (args.xlimit):
#		# Did you set a slope value? 
#		if (args.yslope):
#			selected_objects_i = np.where((y_axis_color > ((args.yslope*x_axis_color) + selection_limit[z])) | (x_axis_color > x_axis_limit))[0]
#		else:
#			selected_objects_i = np.where((y_axis_color > selection_limit[z]) | (x_axis_color > x_axis_limit))[0]
	n_z_above = 0
	for w in range(0, len(selected_objects_i)):
		if (correct_redshifts[selected_objects_i][w] > args.zlimit):
			#print "Yes!"
			n_z_above = n_z_above + 1
	selection_fraction[z] = n_z_above / (1.0*len(selected_objects_i))
	selection_fraction_err[z] = np.sqrt(n_z_above) / (1.0*len(selected_objects_i))
	total_selection_fraction[z] = n_z_above / (1.0*len(high_z_objects_i))
	total_selection_fraction_err[z] = np.sqrt(n_z_above) / (1.0*len(high_z_objects_i))
	if (args.verbose):
		if (args.yslope):
			print 'For yslope = '+str(args.yslope)+', and a y-intercept of '+args.yfilter1+' - '+args.yfilter2+' = '+str(selection_limit[z])+':'
		else:
			print 'For '+args.yfilter1+' - '+args.yfilter2+' > '+str(selection_limit[z])+':'
		print '    There are '+str(len(selected_objects_i))+' total objects with these colors at all redshifts.'
		print '    Of these objects, there are '+str(n_z_above)+' objects at z > '+str(args.zlimit)
	#print n_z_above, len(selected_objects_i), selection_fraction[z]	
selection_colors = cm.rainbow(np.linspace(0, 1, n_selection_limit))

if (args.verbose):
	print "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"

# Now, let's make the higher resolution version for the bottom plot. 
selection_limit_highres = np.arange(0.1, 2.8, 0.2)#np.arange(0.5, 3.2, 0.05)
n_selection_limit_highres = len(selection_limit_highres)
selection_fraction_highres = np.zeros(n_selection_limit_highres)
total_selection_fraction_highres = np.zeros(n_selection_limit_highres)
for z in range(0, n_selection_limit_highres):
	# Did you set an X color limit?
#	if (args.xlimit):
#		# Did you set a slope value? 
#		if (args.yslope):
#			selected_objects_highres_i = np.where((y_axis_color > ((args.yslope*x_axis_color) + selection_limit_highres[z])) | (x_axis_color > x_axis_limit))[0]
#		else:
#			selected_objects_highres_i = np.where((y_axis_color > selection_limit_highres[z]) | (x_axis_color > x_axis_limit))[0]
	if ((args.xplimit) and not (args.xnlimit)):
		# Did you set a slope value? 
		if (args.yslope):
			selected_objects_highres_i = np.where((y_axis_color > ((args.yslope*x_axis_color) + selection_limit_highres[z])) & (x_axis_color < x_axis_plimit))[0]
		else:
			selected_objects_highres_i = np.where((y_axis_color > selection_limit_highres[z]) & (x_axis_color < x_axis_plimit))[0]
	elif ((args.xnlimit) and not (args.xplimit)):
		if (args.yslope):
			selected_objects_highres_i = np.where((y_axis_color > ((args.yslope*x_axis_color) + selection_limit_highres[z])) & (x_axis_color > x_axis_nlimit))[0]
		else:
			selected_objects_highres_i = np.where((y_axis_color > selection_limit_highres[z]) & (x_axis_color > x_axis_nlimit))[0]
	elif ((args.xplimit) and (args.xnlimit)):
		if (args.yslope):
			selected_objects_highres_i = np.where((y_axis_color > ((args.yslope*x_axis_color) + selection_limit_highres[z])) & (x_axis_color < x_axis_plimit) & (x_axis_color > x_axis_nlimit))[0]
		else:
			selected_objects_highres_i = np.where((y_axis_color > selection_limit_highres[z]) & (x_axis_color < x_axis_nlimit) & (x_axis_color > x_axis_nlimit))[0]
	else:
		if (args.yslope):
			selected_objects_highres_i = np.where((y_axis_color > ((args.yslope*x_axis_color) + selection_limit_highres[z])))[0]
		else:
			selected_objects_highres_i = np.where((y_axis_color > selection_limit_highres[z]))[0]
	n_z_above_highres = 0
	for w in range(0, len(selected_objects_highres_i)):
		if (correct_redshifts[selected_objects_highres_i][w] > args.zlimit):
			n_z_above_highres = n_z_above_highres + 1
	selection_fraction_highres[z] = n_z_above_highres / (1.0*len(selected_objects_highres_i))
	total_selection_fraction_highres[z] = n_z_above_highres / (1.0*len(high_z_objects_i))
	#print n_z_above, len(selected_objects_i), selection_fraction[z]	


#plt.clf()
# Color Color Plot
fig = plt.figure(figsize=(10,12))
ax1 = fig.add_subplot(311)
#ax1.scatter(x_axis_color[selected_objects_i], y_axis_color[selected_objects_i], s = 20, c = 'red')
ax1.scatter(x_axis_color[high_z_objects_i], y_axis_color[high_z_objects_i], s = 10, c = 'black')
cax = ax1.scatter(x_axis_color, y_axis_color, s = 5, c = correct_redshifts, edgecolors='none')
cbar = fig.colorbar(cax, label = 'Redshift', orientation="vertical")
for z in range(0, n_selection_limit):
#for z in range(0, 1):
	if ((args.xplimit) and not (args.xnlimit)):
		if (args.yslope):
			ax1.plot([-1, x_axis_plimit],[(args.yslope * (-1))+selection_limit[z], (args.yslope * (x_axis_plimit))+selection_limit[z]], color = selection_colors[z], alpha = 0.2)
		else:
			ax1.plot([-1, x_axis_plimit],[selection_limit[z], selection_limit[z]], color = selection_colors[z], alpha = 0.2)
	elif ((args.xnlimit) and not (args.xplimit)):
		if (args.yslope):
			ax1.plot([args.xnlimit, 4],[(args.yslope * (x_axis_nlimit))+selection_limit[z], (args.yslope * (4))+selection_limit[z]], color = selection_colors[z], alpha = 0.2)
		else:
			ax1.plot([args.xnlimit, 4],[selection_limit[z], selection_limit[z]], color = selection_colors[z], alpha = 0.2)
	elif ((args.xplimit) and (args.xnlimit)):
		if (args.yslope):
			ax1.plot([args.xnlimit, args.xplimit],[(args.yslope * (x_axis_nlimit))+selection_limit[z], (args.yslope * (x_axis_plimit))+selection_limit[z]], color = selection_colors[z], alpha = 0.2)
		else:
			ax1.plot([args.xnlimit, args.xplimit],[selection_limit[z], selection_limit[z]], color = selection_colors[z], alpha = 0.2)
	else:
		if (args.yslope):
			ax1.plot([-1, 4],[(args.yslope * (-1))+selection_limit[z], (args.yslope * (4))+selection_limit[z]], color = selection_colors[z], alpha = 0.2)
		else:
			ax1.plot([-1, 4],[selection_limit[z], selection_limit[z]], color = selection_colors[z], alpha = 0.2)

if (args.xplimit):
	ax1.plot([x_axis_plimit, x_axis_plimit],[(args.yslope * (x_axis_plimit))+selection_limit[0], 6], color = 'black')
if (args.xnlimit):
	ax1.plot([x_axis_nlimit, x_axis_nlimit],[(args.yslope * (x_axis_nlimit))+selection_limit[0], 6], color = 'black')
#ax1.plot([x_axis_limit, x_axis_limit],[-1, 6], color = 'magenta', alpha = 0.2)
#if (args.xlimit):
#	ax1.axvspan(x_axis_limit, 4, facecolor='black', alpha=0.2)
ax1.set_xlabel(args.xfilter1+' - '+args.xfilter2)
ax1.set_ylabel(args.yfilter1+' - '+args.yfilter2)
if (args.yslope):
	plt.title('Slope = '+str(args.yslope))
plt.axis([-1,2,-1,6])

# Y Color vs. Redshift 
ax2 = fig.add_subplot(312)
if (args.yslope):
	ax2.plot([args.zlimit, args.zlimit],[-1,6])
	ax2.scatter(correct_redshifts, y_axis_color - (args.yslope * x_axis_color), color = 'black', s = 5, edgecolors='none')
	ax2.set_ylabel('('+args.yfilter1+' - '+args.yfilter2+') - '+str(args.yslope)+' * ('+args.xfilter1+' - '+args.xfilter2+')', fontsize=7)
else:
	ax2.plot([args.zlimit, args.zlimit],[-1,6])
	ax2.scatter(correct_redshifts, y_axis_color, color = 'black', s = 5, edgecolors='none')
	ax2.set_ylabel(args.yfilter1+' - '+args.yfilter2)

for z in range(0, n_selection_limit):
	if (args.yslope):
		# Is this actually the way to go here? 
		ax2.plot([0, 15],[selection_limit[z], selection_limit[z]], color = selection_colors[z], alpha = 0.2)
	else:
		ax2.plot([0, 15],[selection_limit[z], selection_limit[z]], color = selection_colors[z], alpha = 0.2)

ax2.set_xlabel('Redshift')
plt.axis([0,15,-1,6])


# Selection Fraction
ax3 = fig.add_subplot(313)
#ax3.plot(selection_limit, selection_fraction, color = 'blue', label = 'Fraction of Selected Objects At z > '+str(args.zlimit))
ax3.plot(selection_limit_highres, selection_fraction_highres, color = 'blue', label = 'Fraction of Selected Objects At z > '+str(args.zlimit))
ax3.scatter(selection_limit, selection_fraction, color = 'blue', edgecolors='none')
ax3.errorbar(selection_limit, selection_fraction, selection_fraction_err, color = 'blue', fmt='', ls='none')
#ax3.plot(selection_limit, total_selection_fraction, color = 'red', label = 'Fraction of All z > '+str(args.zlimit)+' Objects')
ax3.plot(selection_limit_highres, total_selection_fraction_highres, color = 'red', label = 'Fraction of All z > '+str(args.zlimit)+' Objects')
ax3.scatter(selection_limit, total_selection_fraction, color = 'red', edgecolors='none')
ax3.errorbar(selection_limit, total_selection_fraction, total_selection_fraction_err, color = 'red', fmt='', ls='none')
crossover_value, crossover_percentage = find_intersection(selection_limit_highres, selection_fraction_highres, total_selection_fraction_highres)

print "\n\n"
if (args.yslope):
	ax3.set_xlabel(args.yfilter1+' - '+args.yfilter2+' (Y-Intercept)')
	print "Crossover at "+args.yfilter1+" - "+args.yfilter2+" (Y-Intercept) = "+str(round(crossover_value,2))
	ax3.text(0.45, 0.4, "Crossover at "+args.yfilter1+" - "+args.yfilter2+" (Y-Intercept) = "+str(round(crossover_value,2)), transform=ax3.transAxes)

else:
	ax3.set_xlabel(args.yfilter1+' - '+args.yfilter2)
	print "Crossover at "+args.yfilter1+" - "+args.yfilter2+" = "+str(round(crossover_value,2))
	ax3.text(0.45, 0.4, "Crossover at "+args.yfilter1+" - "+args.yfilter2+" = "+str(round(crossover_value,2)), transform=ax3.transAxes)

print " Completeness: "+str(round(crossover_percentage,4)*100.0)+"%"
print " Contamination: "+str((1.0 - round(crossover_percentage,4))*100.0)+"%"
ax3.text(0.5, 0.3, "Completeness: "+str(round(crossover_percentage,4)*100.0)+"%", transform=ax3.transAxes)
ax3.text(0.5, 0.2, "Contamination: "+str((1.0 - round(crossover_percentage,4))*100.0)+"%", transform=ax3.transAxes)

ax3.plot([crossover_value, crossover_value], [np.min(selection_limit)-0.2, np.max(selection_limit)+0.2], color = 'black', linestyle = '--')

ax3.set_ylabel('Fraction')
ax3.plot([np.min(selection_limit), np.max(selection_limit)],[0.9,0.9], color = 'black', alpha = 0.4)
ax3.legend(loc = 3)
plt.axis([np.min(selection_limit)-0.2,np.max(selection_limit)+0.2,0.0,1.1])

# Go back and plot on the first plot the crossover point

if ((args.xplimit) and not (args.xnlimit)):
	if (args.yslope):
		ax1.plot([-1, x_axis_plimit],[(args.yslope * (-1))+crossover_value, (args.yslope * (x_axis_plimit))+crossover_value], color = 'black', linestyle = '--')
	else:
		ax1.plot([-1, x_axis_plimit],[crossover_value, crossover_value], color = 'black', linestyle = '--')
elif ((args.xnlimit) and not (args.xplimit)):
	if (args.yslope):
		ax1.plot([args.xnlimit, 4],[(args.yslope * (x_axis_nlimit))+crossover_value, (args.yslope * (4))+crossover_value], color = 'black', linestyle = '--')
	else:
		ax1.plot([args.xnlimit, 4],[crossover_value, crossover_value], color = 'black', linestyle = '--')
elif ((args.xplimit) and (args.xnlimit)):
	if (args.yslope):
		ax1.plot([args.xnlimit, args.xplimit],[(args.yslope * (x_axis_nlimit))+crossover_value, (args.yslope * (x_axis_plimit))+crossover_value], color = 'black', linestyle = '--')
	else:
		ax1.plot([args.xnlimit, args.xplimit],[crossover_value, crossover_value], color = 'black', linestyle = '--')
else:
	if (args.yslope):
		ax1.plot([-1, 4],[(args.yslope * (-1))+crossover_value, (args.yslope * (4))+crossover_value], color = 'black', linestyle = '--')
	else:
		ax1.plot([-1, 4],[crossover_value, crossover_value], color = 'black', linestyle = '--')

# Do the same on the second plot
ax2.plot([0, 15],[crossover_value, crossover_value], color = 'black', linestyle = '--')



if (args.plot_to_screen):
	plt.show()
else:
	if (args.yslope):
		plt.savefig('Color_Color_'+args.yfilter1+'-'+args.yfilter2+'_vs_'+args.xfilter1+'-'+args.xfilter2+'_slope_'+str(args.yslope)+'.png', dpi=300)
	else:
		plt.savefig('Color_Color_'+args.yfilter1+'-'+args.yfilter2+'_vs_'+args.xfilter1+'-'+args.xfilter2+'_slope_0.0.png', dpi=300)
