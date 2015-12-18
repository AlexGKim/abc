#!/usr/bin/env python

import numpy.random
from astropy.cosmology import FlatwCDM
from collections import OrderedDict
import scipy
import sncosmo
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
	__par_names__ = ['H0', 'Om0' ,'w0']
	
	def __init__(self):
		super(Cosmology, self).__init__()
		self.__par_values__ = dict(zip(Cosmology.__par_names__,numpy.array([72,0.28,-1])))

	def distance(self, z):
		cosmology = FlatwCDM(**self.__par_values__)
		return cosmology.luminosity_distance(z).value


class GlobalThroughput(object):
	"""docstring for GlobalThroughput"""
	__input__ = None
	__par_names__ = ['globalzero_des'+x for x in ['g','r','i','z']]

	__filters__ = [sncosmo.get_bandpass('desg'),sncosmo.get_bandpass('desr'),sncosmo.get_bandpass('desi'),
		sncosmo.get_bandpass('desz')]

	def __init__(self):
		super(GlobalThroughput, self).__init__()
		self.__par_values__ = dict(zip(GlobalThroughput.__par_names__,numpy.array([20.,20,20,20])))

	@staticmethod
	def filters():
		return __filters__

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

	model_pars = ['x0','x1', 'c','hostebv', 'hostr_v']
	model_zeros= numpy.array([0,0,0,0,3.1])

	def __init__(self):
		super(SNIa, self).__init__()
		self.__par_values__ = OrderedDict(zip(SNIa.__par_names__,numpy.array([0.14,.31,.1,.1,.1,.1,.1])))

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

	def realize_model_pars(self,host):
		sigs = numpy.array([self.__par_values__['x0_sigma'],
			self.__par_values__['x1_sigma'], self.__par_values__['c_sigma'],
			self.__par_values__['ebv_sigma'], self.__par_values__['r_v_sigma']])
		params= numpy.random.normal(SNIa.model_zeros,sigs)
		return dict(zip(SNIa.model_pars,params))

	def lnlike(self, y):
		sigs = numpy.array([self.__par_values__['x0_sigma'],
			self.__par_values__['x1_sigma'], self.__par_values__['c_sigma'],
			self.__par_values__['ebv_sigma'], self.__par_values__['r_v_sigma']])
		return -0.5*(numpy.sum(((y-SNIa.model_zeros)/sigs)**2)) -0.5 *numpy.log(numpy.sum(sigs**2))


class NonIa(Source):
	"""docstring for NonIa"""
	__input__=None
	__par_names__ = ['x0_sigma','ebv_sigma','r_v_sigma']

	model_pars = ['x0', 'hostebv', 'hostr_v']
	model_zeros = numpy.array([0,0,3.1])

	def __init__(self):
		super(NonIa, self).__init__()
		self.__par_values__ = dict(zip(NonIa.__par_names__,numpy.array([0.6,0.1,0.1])))

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

	def realize_model_pars(self, host):
		sigs=numpy.array([self.__par_values__['x0_sigma'], self.__par_values__['ebv_sigma'],
			self.__par_values__['r_v_sigma']])
		params= numpy.random.normal(NonIa.model_zeros,sigs)
		return dict(zip(NonIa.model_pars,params))

	def lnlike(self, y):
		sigs=numpy.array([self.__par_values__['x0_sigma'], self.__par_values__['ebv_sigma'],
			self.__par_values__['r_v_sigma']])
		return -0.5*(numpy.sum(((y-NonIa.model_zeros)/sigs)**2)) -0.5 *numpy.log(numpy.sum(sigs**2))

class RelativeRates(object):
	"""docstring for RelativeRates"""

	__input__ = None
	__model_pars__ = ['iarate_zmin','iarate_zmax']

	sources = [SNIa, NonIa]

	zmin=0.
	zmax=1.5

	def __init__(self):
		super(RelativeRates, self).__init__()
		self.__par_values__ = dict(zip(RelativeRates.__model_pars__, numpy.array([0.85, 0.2])))

	def rates(self, z, ia_only = False):
		if ia_only:
			return [1.,0.]
		else:
			slope = (self.__par_values__['iarate_zmax']-self.__par_values__['iarate_zmin'])/(RelativeRates.zmax-RelativeRates.zmin)
			iarate = self.__par_values__['iarate_zmin'] + slope * (z-RelativeRates.zmin)
			return [iarate,1-iarate]

class HostGalaxy(object):
	"""docstring for HostGalaxy"""

	__input__ = None
	__par_names__ = ['host_z']

	zmin = 0.1
	zmax = 1.4

	def __init__(self):
		super(HostGalaxy, self).__init__()
		self.__par_values__ = dict(zip(HostGalaxy.__par_names__,
			[numpy.random.uniform(HostGalaxy.zmin**3, HostGalaxy.zmax**3)**(1./3)]))

class Throughput(object):
	"""docstring for Throughput"""

	__input__=[GlobalThroughput]
	__par_names__ = ['zero_des'+x for x in ['g','r','i','x']]

	sigma = 0.01

	def __init__(self, gthroughput):
		super(Throughput, self).__init__()
		self.globalThroughput = gthroughput
		self.__par_values__ = OrderedDict()
		for f in ['g','r','i','z']:
			self.__par_values__['zero+des'+f] = self.globalThroughput.__par_values__['globalzero_des'+f]

	def lnlike(self, y):
		return -0.5*(numpy.sum(((y-self.globalThroughput.__par_values__.values())/Throughput.sigma)**2))

	def realize(self):
		self.__par_values__ = dict(zip(Throughput.__par_names__, numpy.random.normal(loc=self.globalThroughput.__par_values__.values(), scale=Throughput.sigma,
			size=len(Throughput.__par_names__))))


class TypeSubtype(object):
	"""docstring for TypeSubtype"""

	__input__=[RelativeRates, HostGalaxy]
	__par_names0__ = ['type','t0']

	mjd_start = 53000.
	mjd_end = mjd_start + 60

	def __init__(self, rates, host):
		super(TypeSubtype, self).__init__()
		self.rates = rates
		self.host = host

	def realize(self):
		rate =  self.rates.rates(self.host.__par_values__['host_z'])
		draw = numpy.random.uniform()
		self.type = int(draw > rate[0])
		self.source = self.rates.sources[self.type]()
		self.__par_names__=TypeSubtype.__par_names0__ + self.source.model_pars

		self.subtype = self.source.realize_model_pars(self.host)
		self.t0 = numpy.random.uniform(TypeSubtype.mjd_start, TypeSubtype.mjd_end)
		self.__par_values__=dict(zip(self.__par_names__,[self.type, self.t0]))
		self.__par_values__.update(self.subtype)

	def lnlike(self, y):
		pass
		
class Distance(object):
	"""docstring for Distance"""

	__input__=[Cosmology, HostGalaxy]
	__par_names__ = ['luminosity_distance']

	def __init__(self, cosmology, host):
		super(Distance, self).__init__()
		self.cosmology = cosmology
		self.host = host
		self.__par_values__=dict(zip(Distance.__par_names__, [self.cosmology.distance(host.__par_values__['host_z'])]))

class Luminosity(object):
	"""docstring for Luminosity"""

	__input__=[TypeSubtype, HostGalaxy, Source]
	__par_names__=['Luminosity model']

	def __init__(self, modeltype, host):
		super(Luminosity, self).__init__()
		self.type = modeltype
		self.host = host
		self.source = modeltype.source

		self.__par_values__ = dict(zip(Luminosity.__par_names__),[self.source.luminosity(**self.type.subtype)])


class Flux(object):
	"""docstring for Flux"""

	__input__=[Luminosity,HostGalaxy,Distance]
	__par_names__ = ['Flux model']

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

	__input__ = [TypeSubtype]
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
	def __init__(self, flux, throughputs):
		super(Photometry, self).__init__()
		self.flux = flux
		self.throughputs = throughputs

	def photometry(self,mjds,bands):
		ans=dict()
		for b in bands:
			ans[b] = dict()
			corr = numpy.random.normal(0,Photometry.correlated_noise)
			for mjd in mjds:
				try:
					mn=self.flux.model.bandflux(b, mjd, zp = self.throughputs[b][mjd].zeropoints[b], zpsys='ab')
				except:
					mn=0.
				mn = numpy.random.normal(mn,Photometry.sky_sigma)
				ans[b][mjd]=mn+corr
		return ans

class ExtractedPhotometry(object):
	"""docstring for MeasuredPhotometry"""

	__input__ = [Photometry]
	__par_names__ = ["ExtractedPhotometry_o"]
	def __init__(self, arg):
		super(ExtractedPhotometry, self).__init__()
		self.arg = arg

class Survey(object):
	"""docstring for Survey"""
	def __init__(self):
		super(Survey, self).__init__()
		self.mjds = mjds=numpy.arange(53010, 53210,6)
		self.n_sn = 1
		
# class Likelihood(object):
# 	"""docstring for Likelihood"""
# 	def __init__(self, arg):
# 		super(Likelihood, self).__init__()
# 		self.arg = arg

# 	@staticmethod
# 	def spec_redshift0_host0(redhift0_host0, host, cov):

# 	@staticmethod
# 	def phot_redshift0_host0(redhift0_host0, host, cov):

# 	@staticmethod
# 	def type(type, rates, host):

# 	@staticmethod
# 	def type0(type0, type, cov):

# 	@staticmethod
# 	def luminosity(luminosity,type, host, population):

# 	@staticmethod
# 	def photometry0(photometry0, luminosity, cosmology, host, throughput, cov):
		

class SurveyModel(object):
	"""docstring for SurveyModel"""
	def __init__(self, cosmology, rates, host, global_throughput):
		super(SurveyModel, self).__init__()
		self.cosmology = cosmology
		self.rates=rates
		self.host=host
		self.global_throughput=global_throughput

	def realize(self, survey):

		ans = dict()

		for i in xrange(survey.n_sn):
			ans[i] = dict()
			host  = self.host()

			modeltype = ModelType(self.rates, host)
			# print ttype.ttype, ttype.subtype

			distance = Distance(self.cosmology, host)
		#	print distance.luminosity_distance()

			luminosity = Luminosity(modeltype,host)

			flux = Flux(luminosity,host,distance)


			#Data
			specType = SpecType(modeltype)
			ans[i]['Spec Type']= RelativeRates.sources[specType.type_o].__name__

			photRedshift = PhotRedshift(host)
			ans[i]['Phot Type']= photRedshift.redshift_phot_o

			specRedshift = SpecRedshift(host)
			ans[i]['Spec z']= specRedshift.redshift_spec_o

			throughput  = Throughput(self.global_throughput)


			throughputs = dict()
			for b in self.global_throughput.filters():
				throughputs[b] = dict()
				for mjd in survey.mjds:
					throughputs[b][mjd] = throughput

			ans[i]['Throughputs'] = throughputs

			photometry = Photometry(flux, throughputs)

			phot = photometry.photometry(survey.mjds, GlobalThroughput.filters())

			ans[i]['Photometry'] = phot
		return ans
		
"""
A generator of all combinations of observation success given a minimum number of successes. 
"""		
import itertools
# def sampleCombinations(obs_list, ntrue):
# 	num = ntrue
# 	while num <= len(obs_list):
# 		for a in itertools.combinations(obs_list,num):
# 			yield a
# 		num+=1

def sampleCombinations(nlist, ntrue):
	num = ntrue
	indeces=xrange(nlist)
	while num <= nlist:
		for a in itertools.combinations(indeces,num):
			a = numpy.array(a)
			ans = numpy.zeros(nlist,dtype='bool')
			if len(a) != 0:
				ans[a]=True
			yield ans
		num+=1

def sampleProbability(obsnode, conditions, ntrue):
	nlist = len(obsnode)
	ans=0.
	for combo in sampleCombinations(nlist, ntrue):
		combo = numpy.array(combo)
		ans +=numpy.prod(obsnode[combo])* numpy.prod(1-obsnode[numpy.logical_not(combo)])
			# ans +=numpy.prod(obsnode[combo].rcdf(conditions[combo]))* numpy.prod(
			# obsnode[not combo].cdf(conditions[not combo]))
	return ans

print sampleProbability(numpy.arange(.1,.61,.1),xrange(5),5)
