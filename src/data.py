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

def data2(N_sn):

	mjds=numpy.arange(53010, 53210,6)
	out=dict()


	#	gThroughput = GlobalThroughput()
	for i in xrange(N_sn):
		host  = Sim.HostGalaxy()
	#	print host.host_z

	#	model_pars = {'x0':1, 'x1':0, 'c':0 } #,'ebv':0, 'r_v':3.1}
	#	print snIa.luminosity(0.,4400.,**model_pars)



		modeltype = Sim.ModelType(Sim.RelativeRates, host)
		# print ttype.ttype, ttype.subtype

		distance = Sim.Distance(Sim.Cosmology, host)
	#	print distance.luminosity_distance()

		luminosity = Sim.Luminosity(modeltype,host)

		flux = Sim.Flux(luminosity,host,distance)


		#Data
		specType = Sim.SpecType(modeltype)
		print Sim.RelativeRates.sources[specType.type_o].__name__

		photRedshift = Sim.PhotRedshift(host)
		print photRedshift.redshift_phot_o

		specRedshift = Sim.SpecRedshift(host)
		print specRedshift.redshift_spec_o

		throughput  = Sim.Throughput(Sim.GlobalThroughput)


		throughputs = dict()
		for b in Sim.GlobalThroughput.filters():
			throughputs[b] = dict()
			for mjd in mjds:
				throughputs[b][mjd] = throughput

		photometry = Sim.Photometry(flux, throughputs)

		phot = photometry.photometry(mjds, Sim.GlobalThroughput.filters())

		for pkey in phot.keys():
			print pkey.name
			for mkey in phot[pkey]:
				print mkey, phot[pkey][mkey]

def main():
	ans = data()

if __name__ == "__main__":
	main()
