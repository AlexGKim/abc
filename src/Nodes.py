#!/usr/bin/env python

import numpy
from pymc3 import NUTS, Model, Normal, Lognormal, Flat, Bernoulli, Uniform
from astropy.cosmology import FlatwCDM
from pymc3.distributions import Continuous

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
		return numpy.log(self.p*numpy.exp(self.pdf1.logp(value)) +
    		(1-self.p)*numpy.exp(self.pdf2.logp(value)))

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
	rate_Ia_r =1 	: the relative rates are relative to type Ia

	"""

	rate_Ia_r = 1.


	"""
	SN II rates.  For the moment a two-population model is considered;
	generally we will want to consider more populations.

	Parameters
	----------

	z0_snII_r	: float (>=0)
		relative rate of SNe II compared to type Ia.  In the future we
		want the model to be a function of host-galaxy parameters (including
			redshift)

	"""
	rate_II_r = Uniform('rate_II_r', lower=0, upper=10)

	"""
	SN Ia luminosity.  For the moment consider the SN to be time-indepemdent
	with no internal parameters.

	Parameters
	----------

	L_snIa 	: float (>0)
		SN Ia mean luminosity
	sigma_L_snIa : float (>0)
		intrinsic luminosity dispersion (mag)

	"""
	L_snIa = Uniform('L_snIa', lower=0.01, upper=10)
	sigma_L_snIa = Uniform('sigma_L_snIa', lower=0.01, upper=1)


	"""
	SN II luminosity.  For the moment consider the SN to be time-indepemdent
	with no internal parameters.

	Parameters
	----------

	L_snII 	: float (>0)
		SN II mean luminosity
	sigma_L_snII : float (>0)
		intrinsic luminosity dispersion (mag)

	"""
	L_snII = Uniform('L_snII', lower=0.01, upper=10)
	sigma_L_snII = Uniform('sigma_L_snII', lower=0.01, upper=1)


	# Loop through parameters that are object-specific.  Note that
	# the distribution names have the transient index
	for i in xrange(nTrans):


		"""
		Host Redshift.  Eventually we will want to model the host,
		for which redshift is just one parameter.  This is assumed
		to be measured perfectly if there is a spectrum.

		Parameters
		----------

		redshift : float
		"""
		if observation['spectype'][i] is not None:
			redshift = Uniform('redshift'+str(i),lower =0.01, upper =5, observed=observation['spectype'][i])
		else:
			redshift = Uniform('redshift'+str(i),lower =0.01, upper =5)

		"""
		luminosity distance.  This is fixed given the cosmology.

		for the moment the described by a linear model.  Need to implement
		the real solution that needs an integration implemented in theano

		Parameters
		-----------

		luminosity_distance : float
		"""

		luminosity_distance = Om0 +  w0* redshift

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

		Parameters
		----------

		prob : float
			probability of the object being a type Ia
		"""

		prob = rate_Ia_r/(rate_Ia_r+rate_II_r)

		if observation['spectype'][i] is not None:
			"""
			luminosity for the case where there is a spectrum and the type
			is known

			Parameters
			----------

			luminosity 	: float
				intrinsic luminosity of the object
			"""
			if observation['spectype'][i] == 0:
				luminosity = Lognormal('luminosity'+str(i), mu=L_snIa,tau=sigma_L_snIa)
			else:
				luminosity = Lognormal('luminosity'+str(i), mu=L_snII,tau=sigma_L_snII)

			"""
			The observed spectroscopic type

			Note that for bernoulli p is the probability of success and 0 (our
				lablel for SNe Ia) is a failure.

			Parameters
			----------

			ttype :	int
				the type of the object
			"""
				
			ttype = Bernoulli('type '+str(i), 1-prob, observed=observation['spectype'][i])


		else:
			"""
			luminosity for the case where the type is not known

			Parameters
			----------

			luminosity : float
				intrinsic luminosity marginalized over the types
			"""
			luminosity = LuminosityMarginalizedOverType(
				Lognormal('dum1'+str(i), mu=L_snIa,tau=sigma_L_snIa),
				Lognormal('dum2'+str(i), mu=L_snII,tau=sigma_L_snII), prob)

		"""
		flux that arrives to the telescope is deterministic

		Parameters
		-----------

		flux: float
			flux of object that arrives to observatory
		"""

		flux = luminosity/4/numpy.pi/luminosity_distance/luminosity_distance


		"""
		Per observation calibration

		This is not implemented for the moment
		"""

		"""
		counts that are measured

		Parameters
		----------

		counts : float

		"""

		counts = Lognormal('counts'+str(i), mu=flux*10**(-zeropoints/2.5) , tau=0.05, observed=observation['photometry'][i])



from pymc3 import find_MAP, NUTS, sample
from scipy import optimize
with basic_model:

    # obtain starting values via MAP
    start = find_MAP(fmin=optimize.fmin_powell)

    # draw 2000 posterior samples
    trace = sample(2000, start=start)

from pymc3 import traceplot

traceplot(trace);