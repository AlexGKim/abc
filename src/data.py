#!/usr/bin/env python

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pystan
import numpy.random
from astropy.cosmology import FlatwCDM
import scipy
import sncosmo
import pickle
import astropy.units as u
import abc
import copy
import Sim
		

def data():
	survey = Sim.Survey()

	host  = Sim.HostGalaxy
	cosmology = Sim.Cosmology
	rates = Sim.RelativeRates
	global_throughput=  Sim.GlobalThroughput

	sim   = Sim.SurveyModel(cosmology, rates, host, global_throughput)
	ans = sim.realize(survey)

	for k in iter(ans):
		print 'SN ',k
		print ans[k]['Spec Type']
		print ans[k]['Phot Type']
		print ans[k]['Spec z']
		for pkey in ans[k]['Photometry'].keys():
			print pkey.name
			for mkey in ans[k]['Photometry'][pkey]:
				print mkey, ans[k]['Photometry'][pkey][mkey]

def main():
	ans = data()

if __name__ == "__main__":
	main()
