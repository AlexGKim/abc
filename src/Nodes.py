#!/usr/bin/env python

import numpy
from pymc3 import NUTS, Model, Normal, Lognormal, Flat, Bernoulli, Uniform
from astropy.cosmology import FlatwCDM
from pymc3.distributions import Continuous
from pymc3.distributions.dist_math import bound, std_cdf

from astropy import constants as const
from astropy import units as u

import matplotlib.pyplot as plt

import theano.tensor as T

class LuminosityMarginalizedOverType(Continuous):
    r"""The distribution for the luminosity marginalized over two kinds
    of astronomical sources:

	.. math::
    	pdf(L|X) = \sum_i pdf(L|T_ii,X) pdf(T_i|X).

    This class should be generalized to handle multiple types

		
    Parameters
    -----------
    pdf1 : Continuous
    	pdf(L|T1, X)
    pdf2 : Continuous
    	pdf(L|T2, X)
    p : theano.tensor
    	pdf(T1|X), as a consequence pdf(T2|X) = 1-p
    """

    def __init__(self, pdf1, pdf2, p, *args, **kwargs):
    	super(LuminosityMarginalizedOverType, self).__init__(*args, **kwargs)
    	self.pdf1 = pdf1
    	self.pdf2 = pdf2
    	self.p = p

	def logp(self, value):
		"""
		Implementation of Sum_i pdf(L|Ti,X) pdf(Ti|X).
		"""
		return T.log(self.p*numpy.exp(self.pdf1.logp(value)) +
    		(1-self.p)*T.exp(self.pdf2.logp(value)))

class LogLuminosityGivenSpectype(Normal):
    r"""The distribution for the joint spectype and log-luminosity

    It is the product of the probability of the type times the pdf of the luminosity
    which in this case is Lognormal (the Parent class).  Do templates exist in Python?

	.. math::
		p(spectype, luminosity|rate_II_r,L_Ia,L_II) =
						p(spectype,luminosity| ttype=snIa,rate_II_r,L_Ia,L_II)p(ttype=snIa|rate_II_r) +
						p(spectype,luminosity| ttype=snII,rate_II_r,L_Ia,L_II)p(ttype=snII|rate_II_r)
						= p(luminosity | ttype=spectype,L_Ia,L_II) p((ttype=spectype|rate_II_r))

    This class should be generalized to handle multiple types

		
    Parameters
    -----------
    p : Theano.Tensor
    	pdf(T|X)
    """

    def __init__(self, p=1, *args, **kwargs):
    	super(LogLuminosityGivenSpectype, self).__init__(*args, **kwargs)
    	self.p = p

	def logp(self, value):
		"""
		Implementation of Sum_i pdf(L|Ti,X) pdf(Ti|X).
		""" 
		return bound(T.log(p)+ super(LuminosityGivenSpectype, self).logp(value), p > 0, p<=1)

class CountsWithThreshold(Normal):
    r"""The distribution for the joint spectype and log-luminosity

    pdf of counts given a threshold

	.. math::

    This class should be generalized to handle multiple types

		
    Parameters
    -----------
    threshold : theano.tensor
    	minimum counts in order to be discovered

    """

    def __init__(self, threshold=0, *args, **kwargs):
    	super(CountsWithThreshold, self).__init__(*args, **kwargs)
    	self.threshold = threshold
    	self.normalization  = T.log(1-std_cdf((self.threshold-self.mu)/self.sd))

	def logp(self, value):
		"""
		Implementation of Sum_i pdf(L|Ti,X) pdf(Ti|X).
		""" 
		return bound(-self.normalization+ super(CountsWithThreshold, self).logp(value), value >= self.threshold)

def pgm():
	from daft import PGM, Node, Plate
	from matplotlib import rc
	rc("font", family="serif", size=8)
	rc("text", usetex=True)

	pgm = PGM([9.5, 8.5], origin=[0., 0.2], observed_style='inner')

	#pgm.add_node(Node('dispersion',r"\center{$\sigma_{Ia}$ \newline $\sigma_{!Ia}$}", 1,6,scale=1.2,aspect=1.8))
	pgm.add_node(Node('Rate_Ia',r"{SNIa Rate}", 1,8, fixed=1))
	pgm.add_node(Node('Rate_II',r"{SNII Rate}", 2,8,scale=1.6,aspect=1.2))
	pgm.add_node(Node('L_Ia',r"{SNIa L, $\sigma_L$}", 3,8,scale=1.6,aspect=1.2))
	pgm.add_node(Node('L_II',r"{SNII L, $\sigma_L$}", 4,8,scale=1.6,aspect=1.2))
	pgm.add_node(Node('Cosmology',r"Cosmology", 7,8, scale=1.6,aspect=1.2))
	pgm.add_node(Node('Calibration',r"Calibration", 8, 8, scale=1.6,aspect=1.2))

	pgm.add_node(Node('Neighbors',r"\centering{Neighbor \newline Redshift}", 5,7, scale=1.6,aspect=1.2))
	pgm.add_node(Node('Redshift',r"{Redshift}", 6,7, scale=1.6,aspect=1.2))

	pgm.add_node(Node('Type_prob',r"Type prob", 1,6, fixed=1,offset=(20,-10)))
	pgm.add_node(Node('Distance',r"$L_D$", 7,6, fixed=1,offset=(10,10)))

	pgm.add_node(Node('Type',r"Type", 1, 5, scale=1.6,aspect=1.2))

	pgm.add_node(Node('Luminosity',r"Luminosity", 4, 4, scale=1.6,aspect=1.2))
	pgm.add_node(Node('Flux',r"Flux", 7, 3, scale=1.2,fixed=True,offset=(-20,-20)))


	pgm.add_node(Node('Spec_type',r"Spec type", 1, 1, scale=1.6,aspect=1.2,observed=1))
	pgm.add_node(Node('Galaxy_redshift',r"\centering{Galaxy \newline redshift}", 5, 1, scale=1.6,aspect=1.2,observed=1))
	pgm.add_node(Node('Spec_redshift',r"Spec redshift", 6, 1, scale=1.6,aspect=1.2,observed=1))
	pgm.add_node(Node('Counts',r"Counts", 7, 1, scale=1.2,observed=1))


	pgm.add_edge("Rate_Ia","Type_prob")
	pgm.add_edge("Rate_II","Type_prob")

	pgm.add_edge("Cosmology","Distance")
	pgm.add_edge("Redshift","Distance")

	pgm.add_edge("Type_prob", "Type")

	pgm.add_edge("Type","Luminosity")
	pgm.add_edge("L_Ia", "Luminosity")
	pgm.add_edge("L_II", "Luminosity")

	pgm.add_edge("Luminosity","Flux")
	pgm.add_edge("Redshift","Flux")
	pgm.add_edge("Distance","Flux")

	pgm.add_edge("Type","Spec_type")
	pgm.add_edge("Neighbors","Galaxy_redshift")
	pgm.add_edge("Redshift","Galaxy_redshift")
	pgm.add_edge("Redshift","Spec_redshift")

	pgm.add_edge("Flux","Counts")
	pgm.add_edge("Calibration","Counts")

	# Big Plate: Galaxy
	pgm.add_plate(Plate([0.4, 0.5, 8.2, 7.],
	                    label=r"SNe $i = 1, \cdots, N_{SN}$",
	                    shift=-0.2,label_offset=[20,2]))

	pgm.add_plate(Plate([0.5, 3.5, 4., 2.],
	                    label=r"Type $\in \{Ia, II\}$",
	                    shift=-0.2,label_offset=[20,2]))
	# Render and save.

	pgm.render()

	# pgm.figure.text(0.01,0.9,r'\underline{UNIVERSAL}',size='large')
	# pgm.figure.text(0.01,0.55,r'{\centering \underline{INDIVIDUAL} \newline \underline{SN}}',size='large')
	# pgm.figure.text(0.01,0.2,r'\underline{OBSERVATORY}',size='large')
	# pgm.figure.text(0.01,0.1,r'\underline{DATA}',size='large')


	pgm.figure.savefig("nodes_pgm.pdf")

pgm()
wedwe

# the number of transients
nTrans = 15

# set the state of the random number generator
seed=0
numpy.random.seed(seed)

# simulated data in the dictionary observation, including photometry at peak,
# spectroscopic redshift, and spectroscopic type.
# the convention is SNIa are '0', SNII are '1'
# the current implementation is barebones

observation=dict()
observation['specz'] = numpy.random.uniform(low=0.1, high=0.8, size=nTrans)
spectype = numpy.random.uniform(low=0, high=1, size=nTrans)
observation['spectype'] = spectype.round().astype('int')
luminosity = (1.-observation['spectype'])*10**(numpy.random.normal(0, 0.1/2.5, size=nTrans)) \
	+ observation['spectype']*.5**10**(numpy.random.normal(0, 0.4/2.5,size=nTrans))
cosmo = FlatwCDM(H0=72, Om0=0.28, w0=-1)
ld = cosmo.luminosity_distance(observation['specz']).value
h0 = (const.c/cosmo.H0).to(u.Mpc).value


observation['counts'] = luminosity / 4/numpy.pi/ld/ld*10**(0.02/2.5)

count_lim = .4e-8

found  = observation['counts'] >= count_lim
nTrans = found.sum()
print nTrans
observation['specz'] = observation['specz'][found]
observation['spectype'] =observation['spectype'][found]
observation['counts'] =observation['counts'][found]

# plt.scatter(observation['specz'],observation['counts'])
# plt.ylim((0, 2e-8))
# plt.show()


# Create the pymc3 model and fill it with the distributions and parameters
# of the model
basic_model = Model()

with basic_model:

	r"""
	Cosmology Node.  The FlatwCDM cosmology.  We need the flexibility
	to switch in and out different models

	The luminosity distance is specific to the model: it should be the
	job of the class that describes this to provide this.

	Presumably there will be a numerical integration function that inherits
	from theano.Op with gradient implmented so that HMC can be run.  The
	class can be specified by the integrand, which is the gradient.

	Parameters
	----------
	Om0:	Omega_M
	w0:		constant equation of state w
	"""

	Om0 = Lognormal('Om0', mu=numpy.log(0.28), tau=1/.2/.2)
	w0 = Normal('w0', mu=-1, sd=0.2)

	"""
	Calibration Node.  These are the global zeropoints for each band.
	For now we consider one band.

	This class needs to provide the transmission function of each band.

	More complicated parameterizations of calibration are expected.

	Parameters
	-----------
	Z:	zeropoint (in mag) for the bands

	"""
	n_bands = 1
	zeropoints = Normal('zeropoints', mu=0, sd=.02, shape = n_bands)

	"""
	SN Ia Rate Node.  For SN cosmology the relative rates between different
	populations are sufficient.  We will do rates relative to the
	type Ia rate:

	Parameters
	-----------
	rate_Ia_r =1 	: the relative rates are relative to type Ia. Fixed.

	"""

	rate_Ia_r = 1.


	"""
	SN II Rate Node.  For the moment a two-population model is considered;
	generally we will want to consider more populations.

	Parameters
	----------

	z0_snII_r	: relative rate of SNe II compared to SNe Ia.  In the future we
		want the model to be a function of host-galaxy parameters (including
		redshift)

	"""
	rate_II_r = Uniform('rate_II_r', lower=0.1, upper=10)

	"""
	SN Ia luminosity Node.  For the moment consider the SN to be phase-indepemdent
	with no internal parameters.  Eventually this will represent time-evolving
	SED.

	Actually work in log-luminosity

	Parameters
	----------

	logL_snIa 	:		SN Ia mean luminosity
	sigma_snIa :	intrinsic luminosity dispersion (mag)

	"""
	logL_snIa = Normal('logL_snIa', mu=numpy.log(1), sd = 0.02)
	sigma_snIa = Lognormal('sigma_snIa', mu=numpy.log(0.02), tau=1./0.1/0.1)

	"""
	SN II luminosity Node.  For the moment consider the SN to be time-indepemdent
	with no internal parameters.  Eventually this will represent time-evolving
	SED.

	Actually work in log-luminosity

	Parameters
	----------

	L_snII 	: 		SN II mean luminosity
	sigma_snII : 	intrinsic luminosity dispersion (mag)

	"""
	logL_snII = Normal('logL_snII', mu=numpy.log(0.5), sd=0.02)
	sigma_snII = Lognormal('sigma_snII', mu=numpy.log(0.2), tau=1./0.1/0.1)

	# Loop through parameters that are object-specific.
	# Distribution names have the transient index
	for i in xrange(nTrans):

		"""
		Type probability Node.  Probabilities of being a type of object.  For now only SN Ia, and SN II.

		Dependencies
		-------------

		rate_Ia_r	:	Type Ia rate
		rate_II_r	:	Type II rate
		host galaxy :	Not implemented now but eventually depends on host properties


		Parameters
		----------

		prob :			probability of the object being a type Ia.  Fixed.
		"""

		prob = rate_Ia_r/(rate_Ia_r+rate_II_r)


		"""
		Type Node.

		In this implementation the type node is not explicitly considered.

		There are two possibilities:

		1. If there is a spectroscopic type it perfectly determines the type.
		2. Otherwise, type is marginalized in the input model.

		Dependencies
		------------

		prob :

		Parameters
		----------

		ttype :	the type, SN Ia=0, SNII=1

		"""

		if observation['spectype'][i] is not None:

			"""
			Luminosity Node in the case of spectroscopic classification.

			spectype gives perfect spectroscopic typing sp

			p(spectype, luminosity|X) = p(luminosity| type_i, X) p(ttype=spectype| X)

			The distribution for luminosity is handled by the class LuminosityGivenSpectype

			Dependencies
			-------------

			prob
			L_Ia, L_II, sigma_L_snIa, sigma_L_snII

			Parameters
			-----------

			luminosity
			"""

			if observation['spectype'][i] == 0:
				logluminosity = LogLuminosityGivenSpectype('logluminosity'+str(i),
					mu=logL_snIa,sd=sigma_snIa, p=prob)
			else:
				logluminosity = LogLuminosityGivenSpectype('logluminosity'+str(i),
					mu=logL_snII,sd=sigma_snII, p=1-prob)

			luminosity = T.exp(logluminosity)

		else:
			"""
			luminosity Node for the case where the type is not known.  The type Node
			is marginalized internally.

			 pdf(L|X) = \sum_i pdf(L|T_i,X) pdf(T_i|X).

			Dependencies
			------------

			prob 	:	probability of types
			L_snIa, sigma_L_snIa : SN Ia parameters
			L_snII, sigma_L_snII : SN II parameters

			Parameters
			-----------

			luminosity : intrinsic luminosity marginalized over the types
			"""

			logluminosity = LuminosityMarginalizedOverType(
				Normal('dum1'+str(i), mu=logL_snIa,sd=sigma_snIa),
				Normal('dum2'+str(i), mu=logL_snII,sd=sigma_snII), prob)
			luminosity = T.exp(logluminosity)

		if observation['spectype'][i] is not None:

			"""
			Redshift Node.  Is equal to the observed redshift

			See below.
			"""

			redshift = observation['specz']

			# luminosity_distance = Om0 +  w0* observation['specz']
			# flux = luminosity/4/numpy.pi/luminosity_distance/luminosity_distance

		else:

			"""
			Host Redshift Node.  Eventually we will want to model the host,
			for which redshift is just one parameter.  This is assumed
			to be measured perfectly if there is a spectrum.

			There are therefore two scenarios:

			1. If there is a spectroscopic redshift this node is not explicitly considered.
			2. If there is no spectroscopic redshift it is considered.

			Parameters
			----------

			redshift :		Redshift of host galaxy (and hence of transient)
			"""

			redshift = Uniform('redshift'+str(i),lower =0.01, upper =5)

		"""
		Luminosity distance, Flux nodes.  Fixed.

		If there is a spectroscopic redshift

		p(specz, luminosity_distance, flux | redshift, X) = p(luminosity_distance, flux | redshift=specz, X)

		Otherwise

		p(luminosity_distance, flux | redshift, X) = p(luminosity_distance, flux | redshift, X)

		current implementation does not have redshift-dependence of flux but eventually it will.

		Dependencies
		------------

		luminosity  :	luminosity
		redshift  	:	host redshift
		cosmology 	:   cosmology

		Parameters
		-----------

		luminosity_distance : fixed
		flux 				: fixed
		"""

		"""
		poorman luminosity distance
		"""

		# Having problems with theano levels of nesting in finer numerical integration 
		luminosity_distance = (1+redshift)*(0.5 +0.5/T.sqrt(Om0*(1+redshift)**3 + (1-Om0)*(1+redshift)**(3*(1+w0))) \
			+1./T.sqrt(Om0*(1+.25*redshift)**3 + (1-Om0)*(1+.25*redshift)**(3*(1+w0))) \
			+1./T.sqrt(Om0*(1+.50*redshift)**3 + (1-Om0)*(1+.50*redshift)**(3*(1+w0))) \
			+1./T.sqrt(Om0*(1+.75*redshift)**3 + (1-Om0)*(1+.75*redshift)**(3*(1+w0))) \
			)/h0*redshift/4

		flux = luminosity/4/numpy.pi/luminosity_distance**2		


		"""
		Per observation calibration

		This is not implemented for the moment
		"""

		"""
		counts Node.  This is observed.

		eventually will have selection criteria

		Dependencies
		-------------

		flux
		calibration

		Parameters
		----------

		counts : measured counts

		"""
		counts_mu = flux*10**(-zeropoints/2.5)
		counts = CountsWithThreshold('counts'+str(i), threshold = count_lim, mu= counts_mu, sd=0.02*counts_mu, observed=observation['counts'][i])


		"""
		spectral type Node.

		the spectroscopic typing is assumed to be perfect.  It is therefore handled upstream

		Dependencies
		------------
		ttype

		Parameters
		----------
		spectype

		"""

		"""
		spectral redshift Node.

		the spectroscopic redshift is assumed to be perfect.  It is therefore handled upstream

		Dependencies
		------------
		redshift

		Parameters
		----------
		specredshift

		"""

from pymc3 import find_MAP, NUTS, sample
from scipy import optimize
with basic_model:

    # obtain starting values via MAP
    start = find_MAP(fmin=optimize.fmin_powell)

    # draw 2000 posterior samples
    trace = sample(2000, start=start)

from pymc3 import traceplot

traceplot(trace);