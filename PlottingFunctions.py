import os
import sys
import math
import argparse
import numpy as np
#import seaborn as sns
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
#import matplotlib
import matplotlib.pyplot as plt
import astropy.stats
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table

#SMALL_SIZE = 8
#MEDIUM_SIZE = 10
#BIGGER_SIZE = 12

#plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
#plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
#plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
#plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
#plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
#plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
#plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def residuals_68(residuals):
	sorted_residuals = np.sort(abs(residuals))
	number_objects = len(sorted_residuals)
	return sorted_residuals[int(number_objects * 0.68)-1]

def find_catastrophic_outliers(IDs, z_spec, z_phot, outlier_faction_value): 
	residual = abs((z_spec - z_phot) / (1 + z_spec))
	catastrophic_indices = np.where(residual > outlier_faction_value)[0]
	catastrophic_outlier_IDs = IDs[catastrophic_indices]
	catastrophic_outlier_z_spec = z_spec[catastrophic_indices]
	catastrophic_outlier_z_phot = z_phot[catastrophic_indices]
	return catastrophic_outlier_IDs, catastrophic_outlier_z_spec, catastrophic_outlier_z_phot
	

def plotty(data1, data2, color_data, name, filtername, title_for_plot, min_redshift, max_redshift, colorlabel):
	fig, ax = plt.subplots()

	b = np.arange(0,max_redshift+1)

	cax = ax.scatter(data1, data2, s=5.0, c=color_data, edgecolors='none')
	
	ax.plot(b, b, c='black', linewidth = 1.5)

	ax.set_xlim([0.2, max_redshift])
	ax.set_ylim([-0.1, max_redshift])

	plt.plot([0, 15],[0, 12.6], '--', color = 'grey', label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = 0.15')
	plt.plot([0, 15],[0, 17.4], '--', color = 'grey')

	ax.set_xlabel('$z_{spec}$', fontsize=20)
	ax.set_ylabel('$z_{phot}$', fontsize=20)
	cbar = fig.colorbar(cax, label = colorlabel)
	ax.set_title(title_for_plot)
	#plt.show()
	plt.savefig(name, format='png', dpi=600)

def offset_vs_specz(z_spec, z_phot, title, filename, min_redshift, max_redshift):

	# Offset vs. Spec-z
	norm_diff = (z_spec - z_phot) / (1 + z_spec)
	
	plt.clf()
	plt.figure(figsize=(10,10))
	plt.scatter(z_spec, norm_diff, color = 'red', s = 5)
	plt.title(title)
	plt.xlabel('z$_{spec}$')
	plt.ylabel('(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$)')
	plt.plot([min_redshift, max_redshift],[0, 0], color = 'black')
	plt.plot([min_redshift, max_redshift],[-0.15, -0.15], '-', color = 'grey', alpha = 0.2, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = +/- 0.15')
	plt.plot([min_redshift, max_redshift],[0.15, 0.15], '-', color = 'grey', alpha = 0.2)
	plt.plot([min_redshift, max_redshift],[-0.10, -0.10], '-', color = 'grey', alpha = 0.5, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = +/- 0.10')
	plt.plot([min_redshift, max_redshift],[0.10, 0.10], '-', color = 'grey', alpha = 0.5)
	plt.plot([min_redshift, max_redshift],[-0.05, -0.05], '-', color = 'grey', alpha = 0.9, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = +/- 0.05')
	plt.plot([min_redshift, max_redshift],[0.05, 0.05], '-', color = 'grey', alpha = 0.9)
	plt.axis([min_redshift, max_redshift, -5, 1])
	plt.legend()
	plt.savefig(filename, dpi=200)

def photz_specz_offset(z_spec, z_phot, color_data, name, filtername, title_for_plot, min_redshift, max_redshift, colorlabel):
	fig = plt.figure(figsize=(6,8))
	ax1 = fig.add_subplot(211)

	b = np.arange(0,max_redshift+1)

	ax1.scatter(z_spec, z_phot, s=5.0, c=color_data, edgecolors='none')
	
	ax1.plot(b, b, c='black', linewidth = 1.5)

	ax1.set_xlim([0.2, max_redshift])
	ax1.set_ylim([-0.1, max_redshift])

	zspec_vals = np.arange(0,15,0.1)
	zphot_vals_015_pos = zspec_vals + (0.15 * (1+zspec_vals))
	zphot_vals_015_neg = zspec_vals - (0.15 * (1+zspec_vals))
	zphot_vals_050_pos = zspec_vals + (0.50 * (1+zspec_vals))
	zphot_vals_050_neg = zspec_vals - (0.50 * (1+zspec_vals))

	ax1.plot(zspec_vals, zphot_vals_015_pos, '-', color = 'grey', alpha = 0.7, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = 0.15')
	ax1.plot(zspec_vals, zphot_vals_015_neg, '-', color = 'grey', alpha = 0.7, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = 0.15')

	ax1.plot(zspec_vals, zphot_vals_050_pos, '-', color = 'grey', alpha = 0.2, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = 0.50')
	ax1.plot(zspec_vals, zphot_vals_050_neg, '-', color = 'grey', alpha = 0.2, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = 0.50')
	
	ax1.text(0.05, 0.9, 'N$_{objects}$ = '+str(len(z_spec)), transform=ax1.transAxes)
	ax1.set_xlabel('$z_{spec}$')
	ax1.set_ylabel('$z_{phot}$')
	#cbar = fig.colorbar(cax, label = colorlabel, orientation="horizontal")
	ax1.set_title(title_for_plot)
	#plt.show()

	ax2 = fig.add_subplot(212)
	# Offset vs. Spec-z
	norm_diff = (z_spec - z_phot) / (1 + z_spec)
	
	cax = ax2.scatter(z_spec, norm_diff, c = color_data, s = 5, edgecolors='none')
	#ax2.title(title)
	ax2.set_xlabel('z$_{spec}$')
	ax2.set_ylabel('(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$)')
	ax2.plot([min_redshift, max_redshift],[0, 0], color = 'black')
	ax2.plot([min_redshift, max_redshift],[-0.15, -0.15], '-', color = 'grey', alpha = 0.7, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = +/- 0.15')
	ax2.plot([min_redshift, max_redshift],[0.15, 0.15], '-', color = 'grey', alpha = 0.7)
	ax2.plot([min_redshift, max_redshift],[-0.5, -0.5], '-', color = 'grey', alpha = 0.2, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = +/- 0.5')
	ax2.plot([min_redshift, max_redshift],[0.5, 0.5], '-', color = 'grey', alpha = 0.2)
	ax2.axis([min_redshift, max_redshift, -1, 1])
	ax2.legend()
	cbar = fig.colorbar(cax, label = colorlabel, orientation="horizontal")

	plt.savefig(name, format='png', dpi=600)

def photz_specz_density(z_spec, z_phot, name, filtername, title_for_plot, min_redshift, max_redshift, colorlabel):
	fig = plt.figure(figsize=(7,6))
	ax1 = fig.add_subplot(111)

	b = np.arange(0,max_redshift+1)

	# Hexbin plot 
	#good_objects = np.where(z_phot > 0)
	#plt.hexbin(z_spec[good_objects], z_phot[good_objects], gridsize=300, norm=LogNorm(), cmap=plt.cm.viridis)
	
	# Hist2D
	#plt.hist2d(z_spec, z_phot, bins=(100, 100))
	histbins = 100
	good_objects = np.where(z_phot > 0)
	plt.hist2d(z_spec[good_objects], z_phot[good_objects], bins=histbins, norm=LogNorm(), cmap='viridis')
	plt.colorbar(label = colorlabel)

	# SNS scatterplot
	#sns.set_style("white")
	#sns.kdeplot(z_spec, z_phot, n_levels=30, gridsize=300)
	#ax1.scatter(z_spec, z_phot, s=5.0, alpha = 0.1)
	
	# joint plot
	#g = (sns.jointplot(x=z_spec, y=z_phot, xlim = [0.2, max_redshift], ylim = [-0.1, max_redshift], s=5, edgecolor="None"))
	
	ax1.plot(b, b, c='black', linewidth = 1.5, alpha = 0.1)

	ax1.set_xlim([0.2, max_redshift])
	ax1.set_ylim([-0.1, max_redshift])

	ax1.text(0.05, 0.9, 'N$_{objects}$ = '+str(len(z_spec)), transform=ax1.transAxes)
	ax1.set_xlabel('$z_{spec}$')
	ax1.set_ylabel('$z_{phot}$')
	ax1.set_title(title_for_plot)

	plt.savefig(name, format='png', dpi=600)

def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]

def mag_error_histogram(z_spec, z_phot, mag_min, mag_max, name):
	
	sigma_outlier_value = 2.0

	fig = plt.figure(figsize=(6,5))
	ax1 = fig.add_subplot(111)
	
	norm_diff = (z_spec - z_phot) / (1 + z_spec)
	#print "Average = "+str(np.mean(norm_diff))
	#print "Median = "+str(np.median(norm_diff))
	#print "standard deviation = "+str(np.std(norm_diff))
	nooutliers = reject_outliers(norm_diff, m = sigma_outlier_value)
	#print "Average (no outliers) = "+str(np.mean(nooutliers))
	#print "Median (no outliers) = "+str(np.median(nooutliers))
	#print "standard deviation (no outliers) = "+str(np.std(nooutliers))
	mag_error_bins = np.arange(-2.0, 2.0, 0.01)
	plt.hist(norm_diff, bins = mag_error_bins)
	ax1.set_xlabel('(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$)')
	ax1.set_ylabel('Number of Objects Per Bin')
	ax1.set_xlim([-0.3, 0.9])
	plt.title('m(AB,F200W) = '+str(mag_min)+' - '+str(mag_max))

	plt.axvline(x = np.mean(norm_diff), color = 'grey', alpha = 0.7)
	plt.axvline(x=np.mean(norm_diff)+np.std(norm_diff), color = 'grey', alpha = 0.3)
	plt.axvline(x=np.mean(norm_diff)-np.std(norm_diff), color = 'grey', alpha = 0.3)

	plt.axvline(x = np.mean(nooutliers), color = 'red', alpha = 0.7)
	plt.axvline(x=np.mean(nooutliers)+np.std(nooutliers), color = 'red', alpha = 0.3)
	plt.axvline(x=np.mean(nooutliers)-np.std(nooutliers), color = 'red', alpha = 0.3)

	plt.text(0.57, 0.95,'<(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$)> Statistics', horizontalalignment='left',  verticalalignment='center', size = 5, transform = ax1.transAxes)
	plt.text(0.57, 0.90,'All Objects', horizontalalignment='left',  verticalalignment='center', size = 5, transform = ax1.transAxes)
	plt.text(0.6, 0.85,'Average (Median) = '+str(round(np.mean(norm_diff),4))+' ('+str(round(np.median(norm_diff),4))+') +/- '+str(round(np.std(norm_diff),4)), horizontalalignment='left',  verticalalignment='center', size = 5, transform = ax1.transAxes)
	plt.text(0.57, 0.80,'Without '+str(sigma_outlier_value)+' Sigma Outliers', horizontalalignment='left',  verticalalignment='center', size = 5, color = 'red', transform = ax1.transAxes)
	plt.text(0.6, 0.75,'Average (Median) = '+str(round(np.mean(nooutliers),4))+' ('+str(round(np.median(nooutliers),4))+') +/- '+str(round(np.std(nooutliers),4)), horizontalalignment='left',  verticalalignment='center', size = 5, color = 'red', transform = ax1.transAxes)
 
 	plt.savefig(name, format='png', dpi=600)
	

def catastrophic_outlier_plots(z_spec, z_phot, catastrophic_z_spec, catastrophic_z_phot, name, title_for_plot, min_redshift, max_redshift):
	fig = plt.figure(figsize=(6,8))
	ax1 = fig.add_subplot(211)

	b = np.arange(0,max_redshift+1)

#	cax = ax1.scatter(z_spec, z_phot, s=5.0, c=color_data, edgecolors='none')
	ax1.scatter(z_spec, z_phot, s=5.0, color = 'grey', edgecolors='none')
	ax1.scatter(catastrophic_z_spec, catastrophic_z_phot, s=5.0, color = 'red', edgecolors='none')
	
	ax1.plot(b, b, c='black', linewidth = 1.5)

	ax1.set_xlim([0.2, max_redshift])
	ax1.set_ylim([-0.1, max_redshift])

	zspec_vals = np.arange(0,15,0.1)
	zphot_vals_015_pos = zspec_vals + (0.15 * (1+zspec_vals))
	zphot_vals_015_neg = zspec_vals - (0.15 * (1+zspec_vals))
	zphot_vals_050_pos = zspec_vals + (0.50 * (1+zspec_vals))
	zphot_vals_050_neg = zspec_vals - (0.50 * (1+zspec_vals))

	ax1.plot(zspec_vals, zphot_vals_015_pos, '-', color = 'grey', alpha = 0.7, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = 0.15')
	ax1.plot(zspec_vals, zphot_vals_015_neg, '-', color = 'grey', alpha = 0.7, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = 0.15')

	ax1.plot(zspec_vals, zphot_vals_050_pos, '-', color = 'grey', alpha = 0.2, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = 0.50')
	ax1.plot(zspec_vals, zphot_vals_050_neg, '-', color = 'grey', alpha = 0.2, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = 0.50')


	ax1.set_xlabel('$z_{spec}$')
	ax1.set_ylabel('$z_{phot}$')
	ax1.set_title(title_for_plot)

	ax2 = fig.add_subplot(212)
	# Offset vs. Spec-z
	norm_diff = (z_spec - z_phot) / (1 + z_spec)
	norm_diff_catastrophic = (catastrophic_z_spec - catastrophic_z_phot) / (1 + catastrophic_z_spec)
	
	cax = ax2.scatter(z_spec, norm_diff, color = 'grey', s = 5, edgecolors='none')
	cax = ax2.scatter(catastrophic_z_spec, norm_diff_catastrophic, color = 'red', s = 5, edgecolors='none')
	#ax2.title(title)
	ax2.set_xlabel('z$_{spec}$')
	ax2.set_ylabel('(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$)')
	ax2.plot([min_redshift, max_redshift],[0, 0], color = 'black')
	ax2.plot([min_redshift, max_redshift],[-0.15, -0.15], '-', color = 'grey', alpha = 0.7, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = +/- 0.15')
	ax2.plot([min_redshift, max_redshift],[0.15, 0.15], '-', color = 'grey', alpha = 0.7)
	ax2.plot([min_redshift, max_redshift],[-0.5, -0.5], '-', color = 'grey', alpha = 0.2, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) = +/- 0.5')
	ax2.plot([min_redshift, max_redshift],[0.5, 0.5], '-', color = 'grey', alpha = 0.2)
	ax2.axis([min_redshift, max_redshift, -1, 1])
	ax2.legend()

	plt.savefig(name, format='png', dpi=600)
	


def plot_deltazvsz(data1, data2, NIRcSNR, name, filtername):
	fig = plt.figure(figsize=(7,6))

	#f, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(2, 3, sharex=True)
	
	zspec_4_to_6 = data1[np.where((data1 > 4.0) & (data1 < 6.0))]
	zphot_4_to_6 = data2[np.where((data1 > 4.0) & (data1 < 6.0))]
	NIRcSNR_4_to_6 = NIRcSNR[np.where((data1 > 4.0) & (data1 < 6.0))]
	zdiff_4_to_6 = zspec_4_to_6 - zphot_4_to_6
	zfracdiff_4_to_6 = zdiff_4_to_6 / (1 + zspec_4_to_6)
	
	ax1 = fig.add_subplot(321)
	ax1.set_xlim([-0.2,max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0])
	ax1.set_ylim([-0.2, 1.0])
	ax1.set_xlabel('log(SNR$_{'+filtername+'}$)')
	ax1.set_ylabel('($z_{s}$ - $z_{p})$ / (1 + $z_{s}$)')
	ax1.scatter(NIRcSNR_4_to_6, zfracdiff_4_to_6, edgecolor='none', s = 1.8)
	ax1.plot([-0.2, max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0], [0, 0], color='grey', ls='--', alpha = 0.6)
	ax1.plot([math.log10(3.0), math.log10(3.0)], [-10, 10], color='grey', alpha = 0.6)
	ax1.set_title('$z_{spec}$ = 4 - 6')
	
	zspec_6_to_8 = data1[np.where((data1 > 6.0) & (data1 < 8.0))]
	zphot_6_to_8 = data2[np.where((data1 > 6.0) & (data1 < 8.0))]
	NIRcSNR_6_to_8 = NIRcSNR[np.where((data1 > 6.0) & (data1 < 8.0))]
	zdiff_6_to_8 = zspec_6_to_8 - zphot_6_to_8
	zfracdiff_6_to_8 = zdiff_6_to_8 / (1 + zspec_6_to_8)
	
	ax1 = fig.add_subplot(322)
	ax1.set_xlim([-0.2,max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0])
	ax1.set_ylim([-0.2, 1.0])
	ax1.set_xlabel('log(SNR$_{'+filtername+'}$)')
	ax1.set_ylabel('($z_{s}$ - $z_{p})$ / (1 + $z_{s}$)')
	ax1.scatter(NIRcSNR_6_to_8, zfracdiff_6_to_8, edgecolors='none', s = 1.8)
	ax1.plot([-0.2, max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0], [0, 0], color='grey', ls='--', alpha = 0.6)
	ax1.plot([math.log10(3.0), math.log10(3.0)], [-10, 10], color='grey', alpha = 0.6)
	ax1.set_title('$z_{spec}$ = 6 - 8')
	
	zspec_8_to_10 = data1[np.where((data1 > 8.0) & (data1 < 10.0))]
	zphot_8_to_10 = data2[np.where((data1 > 8.0) & (data1 < 10.0))]
	NIRcSNR_8_to_10 = NIRcSNR[np.where((data1 > 8.0) & (data1 < 10.0))]
	zdiff_8_to_10 = zspec_8_to_10 - zphot_8_to_10
	zfracdiff_8_to_10 = zdiff_8_to_10 / (1 + zspec_8_to_10)
	
	ax1 = fig.add_subplot(323)
	ax1.set_xlim([-0.2,max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0])
	ax1.set_ylim([-0.2, 1.0])
	ax1.set_xlabel('log(SNR$_{'+filtername+'}$)')
	ax1.set_ylabel('($z_{s}$ - $z_{p})$ / (1 + $z_{s}$)')
	ax1.scatter(NIRcSNR_8_to_10, zfracdiff_8_to_10, edgecolors='none', s = 1.8)
	ax1.plot([-0.2, max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0], [0, 0], color='grey', ls='--', alpha = 0.6)
	ax1.plot([math.log10(3.0), math.log10(3.0)], [-10, 10], color='grey', alpha = 0.6)
	ax1.set_title('$z_{spec}$ = 8 - 10')

	zspec_10_to_12 = data1[np.where((data1 > 10.0) & (data1 < 12.0))]
	zphot_10_to_12 = data2[np.where((data1 > 10.0) & (data1 < 12.0))]
	NIRcSNR_10_to_12 = NIRcSNR[np.where((data1 > 10.0) & (data1 < 12.0))]
	zdiff_10_to_12 = zspec_10_to_12 - zphot_10_to_12
	zfracdiff_10_to_12 = zdiff_10_to_12 / (1 + zspec_10_to_12)
	
	ax1 = fig.add_subplot(324)
	ax1.set_xlim([-0.2,max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0])
	ax1.set_ylim([-0.2, 1.0])
	ax1.set_xlabel('log(SNR$_{'+filtername+'}$)')
	ax1.set_ylabel('($z_{s}$ - $z_{p})$ / (1 + $z_{s}$)')
	ax1.scatter(NIRcSNR_10_to_12, zfracdiff_10_to_12, edgecolors='none', s = 1.8)
	ax1.plot([-0.2, max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0], [0, 0], color='grey', ls='--', alpha = 0.6)
	ax1.plot([math.log10(3.0), math.log10(3.0)], [-10, 10], color='grey', alpha = 0.6)
	ax1.set_title('$z_{spec}$ = 10 - 12')

	zspec_12_to_14 = data1[np.where((data1 > 12.0) & (data1 < 14.0))]
	zphot_12_to_14 = data2[np.where((data1 > 12.0) & (data1 < 14.0))]
	NIRcSNR_12_to_14 = NIRcSNR[np.where((data1 > 12.0) & (data1 < 14.0))]
	zdiff_12_to_14 = zspec_12_to_14 - zphot_12_to_14
	zfracdiff_12_to_14 = zdiff_12_to_14 / (1 + zspec_12_to_14)
	
	ax1 = fig.add_subplot(325)
	ax1.set_xlim([-0.2,max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0])
	ax1.set_ylim([-0.2, 1.0])
	ax1.set_xlabel('log(SNR$_{'+filtername+'}$)')
	ax1.set_ylabel('($z_{s}$ - $z_{p})$ / (1 + $z_{s}$)')
	ax1.scatter(NIRcSNR_12_to_14, zfracdiff_12_to_14, edgecolors='none', s = 1.8)
	ax1.plot([-0.2, max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0], [0, 0], color='grey', ls='--', alpha = 0.6)
	ax1.plot([math.log10(3.0), math.log10(3.0)], [-10, 10], color='grey', alpha = 0.6)
	ax1.set_title('$z_{spec}$ = 12 - 14')

	zspec_14_to_16 = data1[np.where((data1 > 14.0) & (data1 < 16.0))]
	zphot_14_to_16 = data2[np.where((data1 > 14.0) & (data1 < 16.0))]
	NIRcSNR_14_to_16 = NIRcSNR[np.where((data1 > 14.0) & (data1 < 16.0))]
	
	if (len(zspec_14_to_16) > 0):
		zdiff_14_to_16 = zspec_14_to_16 - zphot_14_to_16
		zfracdiff_14_to_16 = zdiff_14_to_16 / (1 + zspec_14_to_16)
		
		ax1 = fig.add_subplot(326)
		ax1.set_xlim([-0.2,max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0])
		ax1.set_ylim([-0.2, 1.0])
		ax1.set_xlabel('log(SNR$_{'+filtername+'}$)')
		ax1.set_ylabel('($z_{s}$ - $z_{p})$ / (1 + $z_{s}$)')
		ax1.scatter(NIRcSNR_14_to_16, zfracdiff_14_to_16, edgecolors='none', s = 1.8)
		ax1.plot([-0.2, max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0], [0, 0], color='grey', ls='--', alpha = 0.6)
		ax1.plot([math.log10(3.0), math.log10(3.0)], [-10, 10], color='grey', alpha = 0.6)
		ax1.set_title('$z_{spec}$ = 14 - 16')
	
	plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
	
	plt.savefig(name, format='png', dpi = 600)


def plot_histograms(data1, data2, NIRcSNR, name, filtername, title_for_plot):
	fig = plt.figure(figsize=(7,6))

	nbins = 10
	
	xmin = -0.4
	xmax = 4.5
	
	NIRcSNR_4_to_6 = NIRcSNR[np.where((data1 > 4.0) & (data1 < 6.0))]
	NIRcSNR_4_to_6_hist = NIRcSNR_4_to_6[np.where((NIRcSNR_4_to_6 > xmin) & (NIRcSNR_4_to_6 < xmax))]
	ax1 = fig.add_subplot(321)
	ax1.set_xlim([xmin,xmax])
	ax1.hist(NIRcSNR_4_to_6_hist, bins = nbins)
	ax1.set_title('$z_{spec}$ = 4 - 6')
	#ax1.plot([math.log10(3.0), math.log10(3.0)], [0, 10], color='k')
	ax1.set_xlabel('log(SNR$_{'+filtername+'}$)')
	
	NIRcSNR_6_to_8 = NIRcSNR[np.where((data1 > 6.0) & (data1 < 8.0))]
	NIRcSNR_6_to_8_hist = NIRcSNR_6_to_8[np.where((NIRcSNR_6_to_8 > xmin) & (NIRcSNR_6_to_8 < xmax))]
	ax2 = fig.add_subplot(322)
	ax2.set_xlim([xmin,xmax])
	ax2.hist(NIRcSNR_6_to_8_hist, bins = nbins)
	ax2.set_title('$z_{spec}$ = 6 - 8')
	#ax2.plot([math.log10(3.0), math.log10(3.0)], [0, 10], color='k')
	ax2.set_xlabel('log(SNR$_{'+filtername+'}$)')
	
	NIRcSNR_8_to_10 = NIRcSNR[np.where((data1 > 8.0) & (data1 < 10.0))]
	NIRcSNR_8_to_10_hist = NIRcSNR_8_to_10[np.where((NIRcSNR_8_to_10 > xmin) & (NIRcSNR_8_to_10 < xmax))]
	ax3 = fig.add_subplot(323)
	ax3.set_xlim([xmin,xmax])
	ax3.hist(NIRcSNR_8_to_10_hist, bins = nbins)
	ax3.set_title('$z_{spec}$ = 8 - 10')
	#ax3.plot([math.log10(3.0), math.log10(3.0)], [0, 10], color='k')
	ax3.set_xlabel('log(SNR$_{'+filtername+'}$)')

	NIRcSNR_10_to_12 = NIRcSNR[np.where((data1 > 10.0) & (data1 < 12.0))]
	NIRcSNR_10_to_12_hist = NIRcSNR_10_to_12[np.where((NIRcSNR_10_to_12 > xmin) & (NIRcSNR_10_to_12 < xmax))]
	ax4 = fig.add_subplot(324)
	ax4.set_xlim([xmin,xmax])
	ax4.hist(NIRcSNR_10_to_12_hist, bins = nbins)
	ax4.set_title('$z_{spec}$ = 10 - 12')
	#ax4.plot([math.log10(3.0), math.log10(3.0)], [0, 10], color='k')
	ax4.set_xlabel('log(SNR$_{'+filtername+'}$)')

	NIRcSNR_12_to_14 = NIRcSNR[np.where((data1 > 12.0) & (data1 < 14.0))]
	NIRcSNR_12_to_14_hist = NIRcSNR_12_to_14[np.where((NIRcSNR_12_to_14 > xmin) & (NIRcSNR_12_to_14 < xmax))]
	ax5 = fig.add_subplot(325)
	ax5.set_xlim([xmin,xmax])
	ax5.hist(NIRcSNR_12_to_14_hist, bins = nbins)
	ax5.set_title('$z_{spec}$ = 12 - 14')
	#ax5.plot([math.log10(3.0), math.log10(3.0)], [0, 10], color='k')
	ax5.set_xlabel('log(SNR$_{'+filtername+'}$)')

	NIRcSNR_14_to_16 = NIRcSNR[np.where((data1 > 14.0) & (data1 < 16.0))]
	NIRcSNR_14_to_16_hist = NIRcSNR_14_to_16[np.where((NIRcSNR_14_to_16 > xmin) & (NIRcSNR_14_to_16 < xmax))]
	
	if (len(NIRcSNR_14_to_16_hist) > 0):
		ax6 = fig.add_subplot(326)
		ax6.set_xlim([xmin,xmax])
		ax6.hist(NIRcSNR_14_to_16_hist, bins = nbins)
		ax6.set_title('$z_{spec} > 14$')
		#ax6.plot([math.log10(3.0), math.log10(3.0)], [0, 10], color='k')
		ax6.set_xlabel('log(SNR$_{'+filtername+'}$)')

	plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
	
	#fig, ax = plt.subplots()

	#cax = ax.scatter(data1, data2, s=20.0, c=color_data, edgecolors='none')
	
	#ax.plot(b, b, c='black', linewidth = 1.5)

	#ax.set_xlim([3.5, max_redshift])
	#ax.set_ylim([-0.1, max_redshift])

	#ax.set_xlabel('$z_{spec}$', fontsize=20)
	#ax.set_ylabel('$z_{phot}$', fontsize=20)
	#cbar = fig.colorbar(cax, label = 'log(SNR_F444W)')
	#plt.show()
	#ax.set_title(title_for_plot)
	plt.savefig(name, format='png', dpi = 600)

def outlier_fraction_vs_mag(z_spec, z_phot, magnitude_values, name, filtername, title_for_plot):
	# Outlier fraction vs. F200W Magnitude
	norm_diff = (z_spec - z_phot) / (1 + z_spec)
	outlier_value = 0.15
	catastrophic_value = 0.5
	magnitude_min = 20
	magnitude_max = 35
	magnitude_binsize = 0.5
	magnitude_bins = np.arange(magnitude_min,magnitude_max + magnitude_binsize, magnitude_binsize)
	
	outlier_fraction_mag = np.zeros(len(magnitude_bins))
	outlier_fraction_mag_error = np.zeros(len(magnitude_bins))
	catastrophic_outlier_fraction_mag = np.zeros(len(magnitude_bins))
	catastrophic_outlier_fraction_mag_error = np.zeros(len(magnitude_bins))
		
	for mag in range(0, len(magnitude_bins)-1):
		mag_min = magnitude_bins[mag]
		mag_max = magnitude_bins[mag+1]
	
		norm_diff_in_bin = norm_diff[np.where((magnitude_values > mag_min) & (magnitude_values < mag_max))[0]]
		outliers = np.where(norm_diff_in_bin > outlier_value)[0]

		catastrophic_outliers = np.where(norm_diff_in_bin > catastrophic_value)[0]
		if (len(norm_diff_in_bin) > 0):
			outlier_fraction_mag[mag] = len(outliers)*1.0 / len(norm_diff_in_bin)*1.0
			outlier_fraction_mag_error[mag] = np.sqrt(len(outliers)) / len(norm_diff_in_bin)*1.0
			catastrophic_outlier_fraction_mag[mag] = len(catastrophic_outliers)*1.0 / len(norm_diff_in_bin)*1.0
			catastrophic_outlier_fraction_mag_error[mag] = np.sqrt(len(catastrophic_outliers)) / len(norm_diff_in_bin)*1.0
		if (outlier_fraction_mag[mag] == 0):
			outlier_fraction_mag[mag] = -9999
		if (catastrophic_outlier_fraction_mag[mag] == 0):
			catastrophic_outlier_fraction_mag[mag] = -9999	

	#outlier_fraction_mag_error = np.sqrt(outlier_fraction_mag)
	#catastrophic_outlier_fraction_mag_error = np.sqrt(catastrophic_outlier_fraction_mag)
	 
	plt.clf()
	plt.figure(figsize=(10,6))
	plt.scatter(magnitude_bins, outlier_fraction_mag, color = 'black', s = 15, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) $> '+str(outlier_value)+'$')
	plt.errorbar(magnitude_bins, outlier_fraction_mag, yerr=outlier_fraction_mag_error, fmt='none', color='black')
	plt.scatter(magnitude_bins, catastrophic_outlier_fraction_mag, color = 'red', s = 15, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) $> '+str(catastrophic_value)+'$')
	plt.errorbar(magnitude_bins, catastrophic_outlier_fraction_mag, yerr=catastrophic_outlier_fraction_mag_error, fmt='none', color='red')
	plt.legend()
	plt.title(title_for_plot+', Outlier Fraction')
	plt.xlabel(filtername+' magnitude')
	plt.ylabel('f$_{outliers}$')
	plt.plot([24, 32],[0, 0], color = 'black')
	plt.axis([24, 32, -0.05, 0.7])
	plt.savefig(name, dpi = 300)


def outlier_fraction_vs_redshift(z_spec, z_phot, name, title_for_plot):
	# Outlier fraction vs. F200W Magnitude
	norm_diff = (z_spec - z_phot) / (1 + z_spec)
	outlier_value = 0.15
	catastrophic_value = 0.5
	z_min = 0
	z_max = 15
	z_binsize = 0.5
	z_bins = np.arange(z_min,z_max + z_binsize, z_binsize)
	
	outlier_fraction_z = np.zeros(len(z_bins))
	outlier_fraction_z_error = np.zeros(len(z_bins))
	catastrophic_outlier_fraction_z = np.zeros(len(z_bins))
	catastrophic_outlier_fraction_z_error = np.zeros(len(z_bins))
		
	for z in range(0, len(z_bins)-1):
		redshift_min = z_bins[z]
		redshift_max = z_bins[z+1]
	
		norm_diff_in_bin = norm_diff[np.where((z_spec > redshift_min) & (z_spec < redshift_max))[0]]
		outliers = np.where(norm_diff_in_bin > outlier_value)[0]

		catastrophic_outliers = np.where(norm_diff_in_bin > catastrophic_value)[0]
		if (len(norm_diff_in_bin) > 0):
			outlier_fraction_z[z] = len(outliers)*1.0 / len(norm_diff_in_bin)*1.0
			outlier_fraction_z_error[z] = np.sqrt(len(outliers)) / len(norm_diff_in_bin)*1.0
			catastrophic_outlier_fraction_z[z] = len(catastrophic_outliers)*1.0 / len(norm_diff_in_bin)*1.0
			catastrophic_outlier_fraction_z_error[z] = np.sqrt(len(catastrophic_outliers)) / len(norm_diff_in_bin)*1.0
		if (len(norm_diff_in_bin) == 0):
			outlier_fraction_z[z] = -9999
			catastrophic_outlier_fraction_z[z] = -9999	
#		if (outlier_fraction_z[z] == 0):
#			outlier_fraction_z[z] = 0#-9999
#		if (catastrophic_outlier_fraction_z[z] == 0):
#			catastrophic_outlier_fraction_z[z] = 0#-9999	

	plt.clf()
	plt.figure(figsize=(10,6))
	plt.scatter(z_bins, outlier_fraction_z, color = 'black', s = 15, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) $> '+str(outlier_value)+'$')
	plt.errorbar(z_bins, outlier_fraction_z, yerr=outlier_fraction_z_error, fmt='none', color='black')
	plt.scatter(z_bins, catastrophic_outlier_fraction_z, color = 'red', s = 15, label = '(z$_{spec}$ - z$_{phot}$) / (1 + z$_{spec}$) $> '+str(catastrophic_value)+'$')
	plt.errorbar(z_bins, catastrophic_outlier_fraction_z, yerr=catastrophic_outlier_fraction_z_error, fmt='none', color='red')
	plt.legend()
	plt.title(title_for_plot+', Outlier Fraction')
	plt.xlabel('redshift')
	plt.ylabel('f$_{outliers}$')
	plt.plot([0, 15],[0, 0], color = 'black')
	max_y = np.max(catastrophic_outlier_fraction_z)
	#plt.axis([0, 15, -0.05, max_y + (2)*max_y])
	plt.axis([0, 15, -0.05, 0.3])
	plt.savefig(name, dpi = 300)

def q_z_analysis_plots(z_spec, z_phot, q_z, limit_value, name, title_for_plot):
	q_z_limit_values = np.arange(0.0,20.2,0.2)
	residual_all = (z_spec - z_phot) / (1.0 + z_spec)
	
	n_objects_values = np.zeros(len(q_z_limit_values)-1)
	residual_median_values = np.zeros(len(q_z_limit_values)-1)
	residual_mean_values = np.zeros(len(q_z_limit_values)-1)
	residual_stdev_values = np.zeros(len(q_z_limit_values)-1)
	residual_68_values = np.zeros(len(q_z_limit_values)-1)
	NMAD_values = np.zeros(len(q_z_limit_values)-1)
	fraction_gt_15_values = np.zeros(len(q_z_limit_values)-1)
	
	for z in range(0, len(q_z_limit_values)-1):
		good_objects = np.where(q_z < q_z_limit_values[z+1])[0]
		n_objects_values[z] = len(good_objects)
		residuals = (z_spec[good_objects] - z_phot[good_objects]) / (1.0 + z_spec[good_objects])
		residual_68_values[z] = residuals_68(residuals)
		NMAD_values[z] = 1.48 * astropy.stats.median_absolute_deviation(residuals)
		fraction_gt_15_values[z] = len(np.where(abs(residuals) > 0.15)[0])*1.0/len(abs(residuals))*1.0
	
	for z in range(0, len(q_z_limit_values)-1):
		good_objects = np.where((q_z >= q_z_limit_values[z]) & (q_z < q_z_limit_values[z+1]))[0]
		residuals = (z_spec[good_objects] - z_phot[good_objects]) / (1.0 + z_spec[good_objects])
		residual_median_values[z] = np.median(residuals)
		residual_mean_values[z] = np.mean(residuals)
		residual_stdev_values[z] = np.std(residuals)
	
	fig = plt.figure(figsize=(20, 10))
	
	ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=2, rowspan=2) #fig.add_subplot(221)
	ax1.scatter(q_z, residual_all, s = 1, color = 'black', edgecolor = None, alpha = 0.1, label = 'All Objects')
	test_objects = np.where(q_z < limit_value)[0]
	ax1.scatter(q_z[test_objects], residual_all[test_objects], s = 1, color = 'green', edgecolor = None, label = 'Objects with $Q_{z}$ < '+str(limit_value) )
	#ax1.scatter(q_z_limit_values[1:len(q_z_limit_values)]-0.5, residual_mean_values, s = 8, color = 'orange', edgecolor = None)
	ax1.scatter(q_z_limit_values[1:len(q_z_limit_values)]-0.5, residual_median_values, s = 8, color = 'red', edgecolor = None, label = 'Median Values')
	ax1.errorbar(q_z_limit_values[1:len(q_z_limit_values)]-0.5, residual_mean_values, yerr = residual_stdev_values, markersize = 8, color = 'blue', edgecolor = None, alpha = 0.5, label = 'Average Values')
	ax1.set_xlabel('$Q_{z}$', fontsize=15)
	ax1.set_ylabel('<$z_{spec}-z_{phot}$)/(1+$z_{spec}$)>', fontsize=15)
	ax1.set_xlim(0,20)
	ax1.set_ylim(-5,1)
	ax1.plot([0,20], [0,0], color = 'black')
	ax1.text(0.95, 0.05, title_for_plot,
	     horizontalalignment='right',
	     verticalalignment='center',
	     transform = ax1.transAxes)
	ax1.legend()

	ax2 = plt.subplot2grid((2, 4), (0, 2)) #fig.add_subplot(221)
	ax2.scatter(q_z_limit_values[1:len(q_z_limit_values)], n_objects_values, s = 8, color = 'red', edgecolor = None, zorder=10)
	ax2.set_xlabel('Maximum $Q_{z}$', fontsize=10)
	ax2.set_ylabel('Number of Objects', fontsize=10)
	ax2.axvline(x=limit_value, alpha = 0.5, color = 'green')
	f2 = interp1d(q_z_limit_values[1:len(q_z_limit_values)], n_objects_values, kind='cubic')
	ax2.axhline(y=f2(limit_value), alpha = 0.5, color = 'green')
	ax2.text(0.95, 0.05,'For $Q_{z}$ < '+str(limit_value)+', there are '+str(round(f2(limit_value),0))+' Objects',
	     horizontalalignment='right',
	     verticalalignment='center',
	     transform = ax2.transAxes)

	ax3 = plt.subplot2grid((2, 4), (0, 3)) #fig.add_subplot(221)
	ax3.scatter(q_z_limit_values[1:len(q_z_limit_values)], residual_68_values, s = 8, color = 'red', edgecolor = None, zorder=10)
	ax3.set_xlabel('Maximum $Q_{z}$', fontsize=10)
	ax3.set_ylabel('$\sigma_{68}$', fontsize=10)
	ax3.axvline(x=limit_value, alpha = 0.5, color = 'green')
	f2 = interp1d(q_z_limit_values[1:len(q_z_limit_values)], residual_68_values, kind='cubic')
	ax3.axhline(y=f2(limit_value), alpha = 0.5, color = 'green')
	ax3.text(0.95, 0.05,'For $Q_{z}$ < '+str(limit_value)+', $\sigma_{68}$ = '+str(round(f2(limit_value),4)),
	     horizontalalignment='right',
	     verticalalignment='center',
	     transform = ax3.transAxes)

	ax4 = plt.subplot2grid((2, 4), (1, 2)) #fig.add_subplot(221)
	ax4.scatter(q_z_limit_values[1:len(q_z_limit_values)], NMAD_values, s = 8, color = 'red', edgecolor = None, zorder=10)
	ax4.set_xlabel('Maximum $Q_{z}$', fontsize=10)
	ax4.set_ylabel('NMAD', fontsize=10)
	ax4.axvline(x=limit_value, alpha = 0.5, color = 'green')
	f2 = interp1d(q_z_limit_values[1:len(q_z_limit_values)], NMAD_values, kind='cubic')
	ax4.axhline(y=f2(limit_value), alpha = 0.5, color = 'green')
	ax4.text(0.95, 0.05,'For $Q_{z}$ < '+str(limit_value)+', NMAD = '+str(round(f2(limit_value),4)),
	     horizontalalignment='right',
	     verticalalignment='center',
	     transform = ax4.transAxes)

	ax5 = plt.subplot2grid((2, 4), (1, 3)) #fig.add_subplot(221)
	ax5.scatter(q_z_limit_values[1:len(q_z_limit_values)], fraction_gt_15_values, s = 8, color = 'red', edgecolor = None, zorder=10)
	ax5.set_xlabel('Maximum $Q_{z}$', fontsize=10)
	ax5.set_ylabel('fraction > 0.15', fontsize=10)
	ax5.axvline(x=limit_value, alpha = 0.5, color = 'green')
	f2 = interp1d(q_z_limit_values[1:len(q_z_limit_values)], fraction_gt_15_values, kind='cubic')
	ax5.axhline(y=f2(limit_value), alpha = 0.5, color = 'green')
	ax5.text(0.95, 0.05,'For $Q_{z}$ < '+str(limit_value)+', fraction > 0.15 = '+str(round(f2(limit_value),4)),
	     horizontalalignment='right',
	     verticalalignment='center',
	     transform = ax5.transAxes)
	
	plt.tight_layout()
	plt.savefig(name, dpi = 300)
