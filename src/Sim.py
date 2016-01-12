import numpy.random
from astropy.cosmology import FlatwCDM
import scipy
import sncosmo
import abc
import copy

""" 

This is a sandbox for the development of code to describe models based on PGMs
for the SN Hubble cosmology analysis.  It contains implementations of
supernova-specific models.  However, it would be prefered to have an
abstraction beyond what is currently implemented : a framework for describing
models based on PGMs.

The framework is meant to describe models for use in the

1) Simulation of data (in conjunction with a description of the experiment).

2) Fitting a model (in conjunction with data from an experiment).

Therefore the framework needs to make available the likeihood both as a pdf
for data simulation and as a function for model fitting. Maybe the framework
will encompass the simulation and fitting. For the moment we have not settled
on what code will be used:  emcee is an example of code that for a set of
parameters returns a scalar likelihood. real likelihood(array parameters);
there is no candidate code at the moment that realizes data, we may have to
write our own.

The model is assembled through a probabilistic graphical model (PGM).  A PGM
is composed of nodes and nodes are connected by edges.  Each node has its set
of model parameters (names and values).  There are two kinds of parameters:

The probabilistic are described by a pdf

pdf(A|B)

while the fixed are described by a function

A=f(B),

where there is an edge that connects node B to A.

The likelihood for a model is described as follows.  Lets call B the set of
nodes that have no inputs, A the nodes without outputs (i.e. the A are data),
and C are all other nodes.  The pdfs of interest are pdf(AC|B) and
pdf(A|B)=Sum_C pdf(AC|B).  The former is used in fitting, the latter in
simulation.

For practical purposes it is useful to describe the model only in terms of
nodes with continuous parameters, that is to say with no discrete parameters.
If the PGM does have nodes with discrete parameters, a "supernode" that
margninalizes over the discrete parameters.

Each node has

	parameters and associated unique parameter names
	knowledge of other nodes on which it depends (edges)
	pdf(A|B) of A=f(B), where A are the node parameters and B are the
	parameters of what it depends on.  Practically for supernovae we use the
	sncosmo package which has an object from which A=f(B) is a method.

STAN has sections for "data", "parameters", and "transformed parameters" where
the parameters and corresponding pdf's are provided.  STAN knows how to
connect these together.

emcee requires a likelihood.  Our framework could provide a way to take a
bunch of nodes and from them derive the total likelihood.  Danny in his code
just does this by hand.  I think for now doing it by hand is OK. If we find
that users would benefit by having the framework take care of this we can
expand the framework.

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
		self.zeropoints = gthroughput.zeropoints

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
		