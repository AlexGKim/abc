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

class Source(object):

	__metaclass__ = abc.ABCMeta
	__input__ = None
	
	#must have rates rule and luminosity 
	# @abc.abstractmethod
	# def luminosity(self):
	# 	return

	# @abc.abstractmethod
	# def relativeRates(self):
	# 	return

class Cosmology(object):
	"""docstring for Cosmology"""

	__input__ = None
	__par_names__ = ['H0', 'Omega_M' 'w']

	H0 = 72.
	Omega_M = 0.28
	w =-1

	cosmology = FlatwCDM(H0,Omega_M,w)


class Throughput(object):
	"""docstring for Throughput"""
	__input__ = None
	__par_names__ = ['zp']

	filters = [sncosmo.get_bandpass('desg'),sncosmo.get_bandpass('desr'),sncosmo.get_bandpass('desi'),
		sncosmo.get_bandpass('desz')]

	zp = numpy.array([22.,22,22,22])


class HostGalaxy(object):
	"""docstring for HostGalaxy"""

	__input__ = None
	__par_names__ = ['host_z']

	zmin = 0.1
	zmax = 1.4

	def __init__(self):
		super(HostGalaxy, self).__init__()
		self.host_z = numpy.random.uniform(HostGalaxy.zmin**3, HostGalaxy.zmax**3)**(1./3)

		


class SNIa(Source):
	"""docstring for SNIa"""
	__input__ = None
	__par_names__ = ['alpha', 'beta', 'x1_sigma', 'c_sigma', 'ebv_sigma','r_v_sigma']


	model_pars = ['x1', 'c','hostebv', 'hostr_v']

	alpha= 0.14
	beta= 0.31
	x1_sigma = 0.1
	c_sigma = 0.1
	ebv_sigma = 0.1
	r_v_sigma = 0.1

	sphere_area_cm2 = 4 * numpy.pi * (10. * u.pc.to(u.cm))**2

	@staticmethod
	def luminosity(**kwargs):
		#still have to do dust
		# self.
		# self.model = sncosmo.Model(source=self.source, effects=[self.dust],effect_names=['host'], effect_frames=['rest'])

		source = sncosmo.get_source('salt2')

		source.set_peakmag(-19.5-SNIa.alpha*kwargs['x1']+SNIa.beta*kwargs['c'], 'bessellb', 'ab') # magnitude at 10pc

		dust = sncosmo.CCM89Dust()
		model = sncosmo.Model(source=source,effects=[dust], effect_names=['host'], effect_frames=['rest'])
		model.set(**kwargs)
		return model

#		wave = numpy.linspace(source.minwave(), source.maxwave(), 1000)
#		flux = source.flux(phase, wave) # flux in erg/s/cm^2/A at phase=0 at 10pc 
#		flux = source.flux(phase, wave) # flux in erg/s/cm^2/A at phase=0 at 10pc 
		# now integrate flux over wavelength and sphere:
		
		#luminosity = numpy.sum(flux * numpy.gradient(wave)) * sphere_area_cm2 # result is 1.349e43 erg/s
#		return flux * SNIa.sphere_area_cm2

		# return self.model.flux(phase, 4000.)

	@staticmethod
	def realize_model_params(host):
		params= numpy.random.normal([0,0,0,3.1],[SNIa.x1_sigma,SNIa.c_sigma,SNIa.ebv_sigma,SNIa.r_v_sigma])
		return dict(zip(SNIa.model_pars,params))

class NonIa(Source):
	"""docstring for NonIa"""
	__input__=None
	__model_pars__ = None

	@staticmethod
	def luminosity(**kwargs):
		source = sncosmo.get_source("nugent-sn1bc")	
		source.set_peakmag(-20.5, 'bessellb', 'ab') # magnitude at 10pc	
		return sncosmo.Model(source=source)
	@staticmethod
	def realize_model_params(host):
		return dict()

class RelativeRates(object):
	"""docstring for RelativeRates"""

	__input__ = None
	__model_pars_ = ['iarate_zmin','iarate_zmax']


	populations = [SNIa, NonIa]

	zmin=0.
	zmax=1.5

	iarate_zmin = 0.85
	iarate_zmax = 0.2

	#for debugging
	# iarate_zmin = 1
	# iarate_zmax = 1
	slope = (iarate_zmax-iarate_zmin)/(zmax-zmin)

	@staticmethod
	def rates(z):
		iarate = RelativeRates.iarate_zmin + RelativeRates.slope * (z-RelativeRates.zmin)
		return [iarate,1-iarate]

class TType(object):
	"""docstring for TType"""

	__input__=[RelativeRates, HostGalaxy]

	def __init__(self, rates, host):
		super(TType, self).__init__()
		self.rates = rates
		self.host = host
		self.realize()

	def realize(self):
		rate =  self.rates.rates(self.host.host_z)
		draw = numpy.random.uniform()

		self.ttype = int(draw > rate[0])
		self.population = self.rates.populations[self.ttype]
		self.subtype = self.population.realize_model_params(self.host)
		self.t0 = numpy.random.uniform(53000, 53000+1000)

		
class Distance(object):
	"""docstring for Distance"""

	__input__=[Cosmology, HostGalaxy]
	def __init__(self, cosmology, host):
		super(Distance, self).__init__()
		self.cosmology = cosmology
		self.host = host

		self.luminosity_distance = self.cosmology.cosmology.luminosity_distance(self.host.host_z).value

class Luminosity(object):
	"""docstring for Luminosity"""

	__input__=[TType, HostGalaxy, Source]
	def __init__(self, ttype_, host):
		super(Luminosity, self).__init__()
		self.ttype = ttype_
		self.host = host

		self.model = self.ttype.population.luminosity(**self.ttype.subtype)


class Flux(object):
	"""docstring for Flux"""

	__input__=[Luminosity,HostGalaxy,Distance]
	def __init__(self, luminosity, host, distance):
		super(Flux, self).__init__()
		self.luminosity = luminosity
		self.host = host
		self.distance = distance

		if self.luminosity.ttype.ttype ==0:
			index = self.luminosity.model.param_names =='x0'
			x0=self.luminosity.model.parameters[index]
			self.luminosity.model.set(z=host.host_z,x0=x0/4/numpy.pi/distance.luminosity_distance**2)
		else:
			index =self.luminosity.model.param_names =='amplitude'
			x0=self.luminosity.model.parameters[index]
			self.luminosity.model.set(z=host.host_z,amplitude=x0/4/numpy.pi/distance.luminosity_distance**2)
		self.luminosity.model.set(t0=self.luminosity.ttype.t0)
		self.model = self.luminosity.model

class SpecType(object):
	"""docstring for SpecType"""

	__input__ = [TType]

	def __init__(self, ttype):
		super(SpecType, self).__init__()
		self.ttype = ttype.ttype


class SpecRedshift(object):
	"""docstring for ORedshift"""

	__input__ = [HostGalaxy]
	def __init__(self, host):
		super(SpecRedshift, self).__init__()
		self.host = host
		self.redshift = host.host_z

class PhotRedshift(object):
	"""docstring for PhotRedshift"""

	__input__ = [HostGalaxy]
	probability = [0.98,0.02]
	zmin =0.1
	zmax=2
	def __init__(self, host):
		super(PhotRedshift, self).__init__()
		self.host = host
		ran = numpy.random.uniform()
		if ran<PhotRedshift.probability[0]:
			self.redshift=[host.host_z, numpy.random.uniform((PhotRedshift.zmin)**3, (PhotRedshift.zmax)**3)**(1./3)]
		else:
			self.redshift=[numpy.random.uniform((PhotRedshift.zmin)**3, (PhotRedshift.zmax)**3)**(1./3),host.host_z]



						
class Photometry(object):
	"""docstring for OPhotometry"""
	def __init__(self, arg):
		super(Photometry, self).__init__()
		self.arg = arg
		



def main():

	N_sn=5
	cosmology=Cosmology()
#	print cosmology.cosmology.luminosity_distance(0.5)

	throughput = Throughput()
#	print throughput.zp

	rates = RelativeRates()
#	print rates.rates(0.5)

	snIa = SNIa()

	for i in xrange(N_sn):
		host  = HostGalaxy()
	#	print host.host_z

	#	model_pars = {'x0':1, 'x1':0, 'c':0 } #,'ebv':0, 'r_v':3.1}
	#	print snIa.luminosity(0.,4400.,**model_pars)



		ttype = TType(rates, host)
		# print ttype.ttype, ttype.subtype

		distance = Distance(cosmology, host)
	#	print distance.luminosity_distance()

		luminosity = Luminosity(ttype,host)

		flux = Flux(luminosity,host,distance)

		#Data
		specType = SpecType(ttype)
		print specType.ttype

		photRedshift = PhotRedshift(host)
		print photRedshift.redshift

		specRedshift = SpecRedshift(host)
		print specRedshift.redshift

if __name__ == "__main__":
	main()
