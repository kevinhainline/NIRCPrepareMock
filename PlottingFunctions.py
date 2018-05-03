import os
import sys
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table


def plotty(data1, data2, color_data, name, filtername, title_for_plot, min_redshift, max_redshift, colorlabel):
	fig, ax = plt.subplots()

	b = np.arange(0,max_redshift+1)

	cax = ax.scatter(data1, data2, s=5.0, c=color_data, edgecolors='none')
	
	ax.plot(b, b, c='black', linewidth = 1.5)

	ax.set_xlim([0.2, max_redshift])
	ax.set_ylim([-0.1, max_redshift])

	ax.set_xlabel('$z_{spec}$', fontsize=20)
	ax.set_ylabel('$z_{phot}$', fontsize=20)
	cbar = fig.colorbar(cax, label = colorlabel)
	ax.set_title(title_for_plot)
	#plt.show()
	plt.savefig(name, format='png', dpi=600)

def plot_deltazvsz(data1, data2, NIRcSNR, name, filtername):
	fig = plt.figure()

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
	ax1.plot([-0.2, max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0], [0, 0], color='k', ls='--')
	ax1.plot([math.log10(3.0), math.log10(3.0)], [-10, 10], color='k')
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
	ax1.plot([-0.2, max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0], [0, 0], color='k', ls='--')
	ax1.plot([math.log10(3.0), math.log10(3.0)], [-10, 10], color='k')
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
	ax1.plot([-0.2, max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0], [0, 0], color='k', ls='--')
	ax1.plot([math.log10(3.0), math.log10(3.0)], [-10, 10], color='k')
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
	ax1.plot([-0.2, max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0], [0, 0], color='k', ls='--')
	ax1.plot([math.log10(3.0), math.log10(3.0)], [-10, 10], color='k')
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
	ax1.plot([-0.2, max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0], [0, 0], color='k', ls='--')
	ax1.plot([math.log10(3.0), math.log10(3.0)], [-10, 10], color='k')
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
		ax1.plot([-0.2, max(NIRcSNR_4_to_6) + max(NIRcSNR_4_to_6)/10.0], [0, 0], color='k', ls='-')
		ax1.plot([math.log10(3.0), math.log10(3.0)], [-10, 10], color='k')
		ax1.set_title('$z_{spec}$ = 14 - 16')
	
	plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
	
	plt.savefig(name, format='png', dpi = 600)


def plot_histograms(data1, data2, NIRcSNR, name, filtername, title_for_plot):
	fig = plt.figure()

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
	plt.savefig(name, format='png')