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

"""
Each bubble of the PGM is described by a class.  The has internal variables __input__ and __par_names__
to specify the anteceding bubbles and the PGM parameters that are realized in the bubble.  The values
of the parameters are stored in class/object variables.  The parameters are either a class variable or
created by the constructor.


There is one object per bubble.  For example, there is only one Cosmology, therefore the Cosmology
class is a singleton.  There are many HostGalaxy's, and a distinct object is used
to represent each one.

The top-level Universal bubbles (Cosmolgy, Rates, Source, Throughput) have parameters that are fixed.
For simulating data, they are fixed so we can precisely control the top-level configuration.
We may want to change this.  In the analysis these parameters are probabilisitic.

The Luminosity and Flux information are stored as an sncosmo.Model.

All the observed parameters are distinguished with _o.

"""

class Cosmology(object):
	"""docstring for Cosmology"""

	__input__ = None
	__par_names__ = ['H0', 'Omega_M' 'w']

	H0 = 72.
	Omega_M = 0.28
	w =-1

	cosmology = FlatwCDM(H0,Omega_M,w)


class GlobalThroughput(object):
	"""docstring for GlobalThroughput"""
	__input__ = None
	__par_names__ = ['zeropoints'] 

	zeropoints = dict(zip([sncosmo.get_bandpass('desg'),sncosmo.get_bandpass('desr'),sncosmo.get_bandpass('desi'),
		sncosmo.get_bandpass('desz')],numpy.array([20.,20,20,20])))

	@staticmethod
	def filters():
		return GlobalThroughput.zeropoints.keys()


class Source(object):
	"""docstring for Source"""
	__metaclass__ = abc.ABCMeta
	__input__ = None

	model_pars =None
	
	@abc.abstractmethod
	def luminosity(**kwargs):
		return

	@abc.abstractmethod
	def realize_model_pars(**kwargs):
		return	


class SNIa(Source):
	"""docstring for SNIa"""
	__input__ = None
	__par_names__ = ['alpha', 'beta', 'x0_sigma', 'x1_sigma', 'c_sigma', 'ebv_sigma','r_v_sigma']

	alpha= 0.14
	beta= 0.31
	x0_sigma = 0.1 	# mag dispersion
	x1_sigma = 0.1
	c_sigma = 0.1
	ebv_sigma = 0.1
	r_v_sigma = 0.1

	model_pars = ['x0','x1', 'c','hostebv', 'hostr_v']
	@staticmethod
	def luminosity(**kwargs):
		source = sncosmo.get_source('salt2')
		source.set_peakmag(-19.5-SNIa.alpha*kwargs['x1']+SNIa.beta*kwargs['c']+kwargs['x0'], 'desg', 'ab') # magnitude at 10pc
		dust = sncosmo.CCM89Dust()
		model = sncosmo.Model(source=source,effects=[dust], effect_names=['host'], effect_frames=['rest'])
		kwargs2 = dict(kwargs)
		del kwargs2['x0']

		model.set(**kwargs2)
		return model

	@staticmethod
	def realize_model_pars(host):
		params= numpy.random.normal([0,0,0,0,3.1],[SNIa.x0_sigma, SNIa.x1_sigma,SNIa.c_sigma,SNIa.ebv_sigma,SNIa.r_v_sigma])
		return dict(zip(SNIa.model_pars,params))

class NonIa(Source):
	"""docstring for NonIa"""
	__input__=None
	__par_names__ = ['x0_sigma','ebv_sigma','r_v_sigma']

	x0_sigma = 0.6
	ebv_sigma = 0.1
	r_v_sigma = 0.1
	model_pars = ['x0', 'hostebv', 'hostr_v']

	@staticmethod
	def luminosity(**kwargs):
		source = sncosmo.get_source("nugent-sn1bc")
		source.set_peakmag(-17.5+kwargs['x0'], 'desg', 'ab') # magnitude at 10pc	
		dust = sncosmo.CCM89Dust()
		model = sncosmo.Model(source=source,effects=[dust], effect_names=['host'], effect_frames=['rest'])
		kwargs2 = dict(kwargs)
		del kwargs2['x0']

		model.set(**kwargs2)
		return model

	@staticmethod
	def realize_model_pars(host):
		params= numpy.random.normal([0,0,3.1],[NonIa.x0_sigma, NonIa.ebv_sigma,NonIa.r_v_sigma])
		return dict(zip(NonIa.model_pars,params))

class RelativeRates(object):
	"""docstring for RelativeRates"""

	__input__ = None
	__model_pars_ = ['iarate_zmin','iarate_zmax']


	sources = [SNIa, NonIa]

	zmin=0.
	zmax=1.5

	iarate_zmin = 0.85
	iarate_zmax = 0.2

	slope = (iarate_zmax-iarate_zmin)/(zmax-zmin)

	@staticmethod
	def rates(z, ia_only = False):
		if ia_only:
			return [1.,0.]
		else:
			iarate = RelativeRates.iarate_zmin + RelativeRates.slope * (z-RelativeRates.zmin)
			return [iarate,1-iarate]

class HostGalaxy(object):
	"""docstring for HostGalaxy"""

	__input__ = None
	__par_names__ = ['host_z']

	zmin = 0.1
	zmax = 1.4

	def __init__(self):
		super(HostGalaxy, self).__init__()
		self.host_z = numpy.random.uniform(HostGalaxy.zmin**3, HostGalaxy.zmax**3)**(1./3)		


class Throughput(object):
	"""docstring for Throughput"""

	__input__=[GlobalThroughput]

	def __init__(self, gthroughput):
		super(Throughput, self).__init__()
		self.global_throughput = gthroughput
		self.zeropoints = gthroughput.filters

	def filters(self):
		return self.global_throughput.filters()
		

class ModelType(object):
	"""docstring for ModelType"""

	__input__=[RelativeRates, HostGalaxy]
	__par_names__ = ['type','subtype','t0']

	mjd_start = 53000.
	mjd_end = mjd_start + 60

	def __init__(self, rates, host):
		super(ModelType, self).__init__()
		self.rates = rates
		self.host = host
		self.realize()

	def realize(self):
		rate =  self.rates.rates(self.host.host_z)
		draw = numpy.random.uniform()

		self.type = int(draw > rate[0])
		self.source = self.rates.sources[self.type]
		self.subtype = self.source.realize_model_pars(self.host)
		self.t0 = numpy.random.uniform(ModelType.mjd_start, ModelType.mjd_end)

		
class Distance(object):
	"""docstring for Distance"""

	__input__=[Cosmology, HostGalaxy]
#	__par_names__ = ['luminosity_distance']

	def __init__(self, cosmology, host):
		super(Distance, self).__init__()
		self.cosmology = cosmology
		self.host = host
		self.luminosity_distance = cosmology.cosmology.luminosity_distance(host.host_z).value

class Luminosity(object):
	"""docstring for Luminosity"""

	__input__=[ModelType, HostGalaxy, Source]
	__par_names__=['model']

	def __init__(self, modeltype, host):
		super(Luminosity, self).__init__()
		self.type = modeltype
		self.host = host
		self.source = modeltype.source

		self.model = self.source.luminosity(**self.type.subtype)


class Flux(object):
	"""docstring for Flux"""

	__input__=[Luminosity,HostGalaxy,Distance]
	__par_names__ = ['model']

	def __init__(self, luminosity, host, distance):
		super(Flux, self).__init__()
		self.luminosity = luminosity
		self.host = host
		self.distance = distance

		self.model = copy.deepcopy(self.luminosity.model)

		if self.luminosity.type.type ==0:
			x0=self.model.get('x0')
			self.model.set(z=host.host_z,x0=x0/distance.luminosity_distance**2/1e10)
		else:
			x0=self.model.get('amplitude')
			self.model.set(z=host.host_z,amplitude=x0/distance.luminosity_distance**2/1e10)

		self.model.set(t0=self.luminosity.type.t0)
		

class SpecType(object):
	"""docstring for SpecType"""

	__input__ = [ModelType]
	__par_names__ = ['type_o']

	def __init__(self, modeltype):
		super(SpecType, self).__init__()
		self.modeltype = modeltype
		self.type_o = modeltype.type


class SpecRedshift(object):
	"""docstring for ORedshift"""

	__input__ = [HostGalaxy]
	__par_names__ = 'redshift_spec_o'
	def __init__(self, host):
		super(SpecRedshift, self).__init__()
		self.host = host
		self.redshift_spec_o = host.host_z

class PhotRedshift(object):
	"""docstring for PhotRedshift"""

	__input__ = [HostGalaxy]
	__par_names__ = 'redshift_phot_o'

	probability = [0.98,0.02]
	zmin =0.1
	zmax=2
	def __init__(self, host):
		super(PhotRedshift, self).__init__()
		self.host = host

		ran = numpy.random.uniform()
		otherz = numpy.random.uniform((PhotRedshift.zmin)**3, (PhotRedshift.zmax)**3)**(1./3)
		if ran<PhotRedshift.probability[0]:
			self.redshift_phot_o=[host.host_z, otherz]
		else:
			self.redshift_phot_o=[otherz, host.host_z]



						
class Photometry(object):
	"""docstring for Photometry"""

	__input__ = [Flux, Throughput]
	__par_names__ = ["Photometry_o"]

	sky_sigma = 1.
	correlated_noise = .01
	def __init__(self, flux, throughput):
		super(Photometry, self).__init__()
		self.flux = flux
		self.throughput = throughput

	def photometry(self,mjds,bands):
		ans=dict()
		for b in bands:
			ans[b] = dict()
			corr = numpy.random.normal(0,Photometry.correlated_noise)
			for mjd in mjds:
				try:
					mn=self.flux.model.bandflux(b, mjd, zp = self.throughput.zeropoints[b], zpsys='ab')
				except:
					mn=0.

				mn = numpy.random.normal(mn,Photometry.sky_sigma)
				ans[b][mjd]=mn+corr
		return ans


def data(N_sn):

	out=dict()


	gThroughput = GlobalThroughput()
	for i in xrange(N_sn):
		host  = HostGalaxy()
	#	print host.host_z

	#	model_pars = {'x0':1, 'x1':0, 'c':0 } #,'ebv':0, 'r_v':3.1}
	#	print snIa.luminosity(0.,4400.,**model_pars)



		modeltype = ModelType(RelativeRates, host)
		# print ttype.ttype, ttype.subtype

		distance = Distance(Cosmology, host)
	#	print distance.luminosity_distance()

		luminosity = Luminosity(modeltype,host)

		flux = Flux(luminosity,host,distance)

		photometry = Photometry(flux, gThroughput)

		#Data
		specType = SpecType(modeltype)
		print RelativeRates.sources[specType.type_o].__name__

		photRedshift = PhotRedshift(host)
		print photRedshift.redshift_phot_o

		specRedshift = SpecRedshift(host)
		print specRedshift.redshift_spec_o


		phot = photometry.photometry(numpy.arange(53010, 53210,6), gThroughput.filters())

		for pkey in phot.keys():
			print pkey.name
			for mkey in phot[pkey]:
				print mkey, phot[pkey][mkey]

def main():
	ans = data(5)

if __name__ == "__main__":
	main()
