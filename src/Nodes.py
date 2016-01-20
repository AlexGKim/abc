#!/usr/bin/env python

import numpy
from pymc3 import NUTS, Model, Normal, Lognormal, Flat, Bernoulli, Uniform
from astropy.cosmology import FlatwCDM
from pymc3.distributions import Continuous
from pymc3.distributions.dist_math import bound
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
    p : Theano.Tensor
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

class LuminosityGivenSpectype(Lognormal):
    r"""The distribution for the joint spectype and luminosity

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
    	super(LuminosityGivenSpectype, self).__init__(*args, **kwargs)
    	self.p = p

	def logp(self, value):
		"""
		Implementation of Sum_i pdf(L|Ti,X) pdf(Ti|X).
		""" 
		mu = self.mu
		tau = self.tau
		return bound(T.log(p) + super(LuminosityGivenSpectype, self).logp(value), tau > 0)

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
	pgm.add_node(Node('Cosmology',r"Cosmology", 6,8, scale=1.6,aspect=1.2))
	pgm.add_node(Node('Calibration',r"Calibration", 7, 8, scale=1.6,aspect=1.2))

	pgm.add_node(Node('Redshift',r"{Redshift}", 5,7, scale=1.6,aspect=1.2))

	pgm.add_node(Node('Type_prob',r"Type prob", 1,6, fixed=1,offset=(20,-10)))
	pgm.add_node(Node('Distance',r"$L_D$", 6,6, fixed=1,offset=(10,10)))

	pgm.add_node(Node('Type',r"Type", 1, 5, scale=1.6,aspect=1.2))

	pgm.add_node(Node('Luminosity',r"Luminosity", 4, 4, scale=1.6,aspect=1.2))
	pgm.add_node(Node('Flux',r"Flux", 6, 3, scale=1.2,fixed=True,offset=(-20,-20)))


	pgm.add_node(Node('Spec_type',r"Spec type", 1, 1, scale=1.6,aspect=1.2,observed=1))
	pgm.add_node(Node('Spec_redshift',r"Spec redshift", 5, 1, scale=1.6,aspect=1.2,observed=1))
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
	pgm.add_edge("Redshift","Spec_redshift")

	pgm.add_edge("Flux","Counts")
	pgm.add_edge("Calibration","Counts")

	# Big Plate: Galaxy
	pgm.add_plate(Plate([0.4, 0.5, 7.2, 7.],
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

# pgm()

# wefwe

# the number of transients
nTrans = 20

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
luminosity = (1.-observation['spectype']) + observation['spectype']*.5
ld = FlatwCDM(H0=72, Om0=0.28, w0=-1).luminosity_distance(observation['specz']).value
observation['photometry'] = luminosity / 4/numpy.pi/ld/ld

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

	Om0 = Uniform('Om0',lower=0.0, upper=1)
	w0 = Uniform('w0', lower=-2, upper=2)

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
	SN Ia rates.  For SN cosmology the relative rates between different
	populations are sufficient.  We will do rates relative to the
	type Ia rate:

	Parameters
	-----------
	rate_Ia_r =1 	: the relative rates are relative to type Ia. Fixed.

	"""

	rate_Ia_r = 1.


	"""
	SN II rates.  For the moment a two-population model is considered;
	generally we will want to consider more populations.

	Parameters
	----------

	z0_snII_r	: relative rate of SNe II compared to SNe Ia.  In the future we
		want the model to be a function of host-galaxy parameters (including
		redshift)

	"""
	rate_II_r = Uniform('rate_II_r', lower=0.1, upper=10)

	"""
	SN Ia luminosity.  For the moment consider the SN to be phase-indepemdent
	with no internal parameters.  Eventually this will represent time-evolving
	SED.


	Parameters
	----------

	L_snIa 	:		SN Ia mean luminosity
	sigma_L_snIa :	intrinsic luminosity dispersion (mag)

	"""
	L_snIa = Uniform('L_snIa', lower=0.01, upper=10)
	sigma_L_snIa = Uniform('sigma_L_snIa', lower=0.01, upper=.2)
	logL_snIa = T.log(L_snIa)
	tau_snIa = 1/sigma_L_snIa/sigma_L_snIa


	"""
	SN II luminosity.  For the moment consider the SN to be time-indepemdent
	with no internal parameters.  Eventually this will represent time-evolving
	SED.

	Parameters
	----------

	L_snII 	: 		SN II mean luminosity
	sigma_L_snII : 	intrinsic luminosity dispersion (mag)

	"""
	L_snII = Uniform('L_snII', lower=-10, upper=10)
	sigma_L_snII = Uniform('sigma_L_snII', lower=0.01, upper=1)
	logL_snII = T.log(L_snII)
	tau_snII = 1/sigma_L_snII/sigma_L_snII


	# Loop through parameters that are object-specific.
	# Distribution names have the transient index
	for i in xrange(nTrans):

		"""
		Probabilities of being a type of object.  For now only SN Ia, and SN II.

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
		Type of the object. The type is a discrete parameter for the type
		of the object.  To simplify the model we do the following:

		1. Assume that spectroscopic classification is perfect.
		2. When spectroscopic classification is not available, marginalize 
		over this parameter.

		Otherwise, we would have to represent

		pdf(type, luminosity | X),

		which is challenging.  I am not ever sure that mcpy3 can handle a
		mixed discrete and continuous pdf.

		Therefore in this code type not considered explicitly

		Dependencies
		------------

		prob :

		Parameters
		----------

		ttype :	the type, SN Ia=0, SNII=1

		"""

		if observation['spectype'][i] is not None:

			"""
			luminosity and spectype when spectype gives perfect spectroscopic typing.

			p(spectype, luminosity|X) = p(luminosity| type_i, X) p(ttype=spectype| X)

			which is handled by the class LuminosityGivenSpectype

			Dependencies
			-------------

			prob
			L_Ia, L_II, sigma_L_snIa, sigma_L_snII

			Parameters
			-----------

			spectype
			luminosity
			"""

			if observation['spectype'][i] == 0:
				luminosity = LuminosityGivenSpectype('luminosity'+str(i),
					mu=logL_snIa,tau=tau_snIa, p=prob)
			else:
				luminosity = LuminosityGivenSpectype('luminosity'+str(i),
					mu=logL_snII,tau=tau_snII, p=1-prob)

		else:
			"""
			luminosity for the case where the type is not known

			 pdf(L|X) = \sum_i pdf(L|T_ii,X) pdf(T_i|X).

			Dependencies
			------------

			prob 	:	probability of types
			L_snIa, sigma_L_snIa : SN Ia parameters
			L_snII, sigma_L_snII : SN II parameters

			Parameters
			-----------

			luminosity : intrinsic luminosity marginalized over the types
			"""
			luminosity = LuminosityMarginalizedOverType(
				Lognormal('dum1'+str(i), mu=logL_snIa,tau=tau_snIa),
				Lognormal('dum2'+str(i), mu=logL_snII,tau=tau_snII), prob)

		"""
		Host Redshift.  Eventually we will want to model the host,
		for which redshift is just one parameter.  This is assumed
		to be measured perfectly if there is a spectrum.


		Parameters
		----------

		redshift :		Redshift of host galaxy (and hence of transient)
		"""
		if observation['spectype'][i] is not None:
			redshift = observation['spectype'][i]
		else:
			redshift = Uniform('redshift'+str(i),lower =0.01, upper =5)

		"""
		luminosity distance.  This is fixed given the cosmology.

		for the moment the described by a linear model.  Need to implement
		the real solution that needs an integration implemented in theano

		Dependencies
		------------

		redshift  	:	host redshift
		cosmology 	:   cosmology

		Parameters
		-----------

		luminosity_distance : fixed
		"""

		luminosity_distance = Om0 +  w0* redshift




		"""
		flux that arrives to the telescope is deterministic

		Dependencies
		------------

		luminosity
		redshift
		luminosity_distance

		Parameters
		-----------

		flux: 		flux of object that arrives to observatory. Fixed.
		"""

		flux = luminosity/4/numpy.pi/luminosity_distance/luminosity_distance


		"""
		Per observation calibration

		This is not implemented for the moment
		"""

		"""
		counts.  This is observed.

		Dependencies
		-------------

		flux
		calibration

		Parameters
		----------

		counts : measured counts

		"""
		counts_mu = flux*10**(-zeropoints/2.5)
		counts = Normal('counts'+str(i), mu= counts_mu, sd=0.02*counts_mu, observed=observation['photometry'][i])



from pymc3 import find_MAP, NUTS, sample
from scipy import optimize
with basic_model:

    # obtain starting values via MAP
    start = find_MAP(fmin=optimize.fmin_powell)

    # draw 2000 posterior samples
    trace = sample(2000, start=start)

from pymc3 import traceplot

traceplot(trace);