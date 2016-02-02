#!/usr/bin/env python
import numpy
from pymc3 import NUTS, Model, Normal, Lognormal, Flat, Bernoulli, Uniform
from astropy.cosmology import FlatwCDM
from pymc3.distributions import Continuous
from pymc3.distributions.dist_math import bound, std_cdf
from pymc3.backends import SQLite

from astropy import constants as const
from astropy import units as u

import matplotlib.pyplot as plt

import theano
import theano.tensor as T
from theano import pp
from theano.compile.ops import as_op

cosmo = FlatwCDM(H0=72, Om0=0.28, w0=-1)
h0 = (const.c/cosmo.H0).to(u.Mpc).value

def hinv(z, Om0, w0):
    return T.pow((1+z)**3 * (Om0 + (1-Om0)*(1+z)**(3*w0)),-0.5)

def luminosity_distance(z, Om0, w0):
    return 0.5/h0*(z+z**2)*(1+ hinv(z, Om0, w0))

def normalization_integrand(lnLoversigma, lnL, sigma, Z, threshold, ld):
    luminosity = T.exp(lnLoversigma*sigma-lnL)
    flux = (1./4/numpy.pi)*luminosity/ld**2
    counts = flux*numpy.exp(Counts.ln10/2.5)*T.exp(Z)
    ccdf  = 1-std_cdf((threshold-counts)/1e-9)
    return numpy.exp(-lnLoversigma**2/2.)*ccdf

def normalization_integral(lnL, sigma, Z, threshold, ld):
    # integral with coordinates in lnL/sigma units
    # coarse trapezoidal integral over L
    # This ugliness is due to nesting problems in theano
    #for index in xrange(-3,4):
    return 0.5*normalization_integrand(-3, lnL, sigma, Z, threshold, ld) +\
        normalization_integrand(-2, lnL, sigma, Z, threshold, ld) + \
        normalization_integrand(-1, lnL, sigma, Z, threshold, ld) + \
        normalization_integrand(0, lnL, sigma, Z, threshold, ld) + \
        normalization_integrand(1, lnL, sigma, Z, threshold, ld) + \
        normalization_integrand(2, lnL, sigma, Z, threshold, ld) + \
        0.5 * normalization_integrand(3, lnL, sigma, Z, threshold, ld)

#   Maybe want to define a theano function for luminosity distance.  This is a first crack
#   that doesn't quire work
#
# class LuminosityDistance(theano.Op):
#     """
#     This creates an Op that takes x to a*x+b.
#     """
#     __props__ = ("z")

#     def __init__(self, z):
#         super(LuminosityDistance, self).__init__()
#         self.z = z

#     def make_node(self, Om0, w0):
#         # check that the theano version has support for __props__.
#         assert hasattr(self, '_props'), "Your version of theano is too old to support __props__."
#         Om0 = theano.tensor.as_tensor_variable(Om0)
#         w0 = theano.tensor.as_tensor_variable(w0)

#         return theano.Apply(self, [Om0, w0], [theano.tensor.lscalar()])

#     def h2(self, Om0, w0):
#         # h2(z=0) = 1
#         return (1+self.z)**3 * (Om0 + (1-Om0)*(1+self.z)**(3*w0))

#     def dh2dOm0 (self, Om0, w0):
#         # dh2dOm0 h2(z=0) = 0
#         return (1+self.z)**3 * (1 - (1+self.z)**(3*w0))

#     def dh2dw0(self, Om0, w0):
#         # dh2dw0 = 0
#         return (1+self.z)**3 * T.log(1+self.z) * (1+self.z)**(3*w0)

#     def perform(self, node, inputs, output_storage):

#         """
#         Poorman luminosity distance.

#         Presumably there will be a numerical integration function that inherits from theano.Op with
#         gradient implmented so that HMC can be run.  The class can be specified by the integrand,
#         which is the gradient.

#         Inputs
#         ------
#         z :         Theano.lscalar
#             redshift
#         Om0 :       Theano.lscalar
#             Omega_M
#         w :         Theano.lscalar
#             w0

#         Output
#         ------

#         luminosity distance : Theano.lscaoar
#             luminosity distance in Mp units
#         """

#         Om0 = inputs[0]
#         w0  = inputs[1]

#         z = output_storage[0]

#         # luminosity_distance = (1+z)*(0.5 \
#         #     +0.5/T.sqrt(self.Om0*(1+z)**3 + (1-self.Om0)*(1+z)**(3*(1+self.w0))) \
#         #     +1./T.sqrt(self.Om0*(1+.25*z)**3 + (1-self.Om0)*(1+.25*z)**(3*(1+self.w0))) \
#         #     +1./T.sqrt(self.Om0*(1+.50*z)**3 + (1-self.Om0)*(1+.50*z)**(3*(1+self.w0))) \
#         #     +1./T.sqrt(self.Om0*(1+.75*z)**3 + (1-self.Om0)*(1+.75*z)**(3*(1+self.w0))) \
#         #     )/CountsWithThreshold.h0*z/4
#         z[0] = 0.5/h0*(self.z+T.sqr(self.z))*(1+ 1//T.sqrt(self.h2(Om0, w0)))


#     def infer_shape(self, node, i0_shapes):
#         return i0_shapes

#     def grad(self, inputs, output_grads):
#         Om0 = inputs[0]
#         w0  = inputs[1]
#         return [-0.5 * 0.5/h0*(self.z+T.sqr(self.z))* \
#             T.pow(h2(self.z, Om0, w0),-1.5) * self.dh2dOm0(Om0, w0)* output_grads[0],  \
#             + -0.5 * 0.5/h0*(self.z+T.sqr(self.z))* \
#             T.pow(h2(self.z, Om0, w0),-1.5) * self.dh2dw0(Om0, w0)* output_grads[1]]

class LogLuminosityMarginalizedOverType(Continuous):
    r"""The distribution for the luminosity marginalized over two kinds
    of astronomical sources:

    .. math::
        pdf(Luminosity | Type prob, logL_snIa, logL_snII)
            = sum_i pdf(Luminosity | Type_i, logL_snIa, logL_snII) *
                pdf(Type_i | Type prob)


    This class should be generalized to handle multiple types each with its own model

    Parameters
    -----------
    mus : array
        the logL_X
    sds : array
        the sigma_X
    p : theano.tensor
        The probability of the first case
        pdf(T1|X), as a consequence pdf(T2|X) = 1-p
    """

    def __init__(self, mus=numpy.zeros(2), sds=1.+numpy.zeros(2), p=0.5, *args, **kwargs):
        super(LogLuminosityMarginalizedOverType, self).__init__(*args, **kwargs)
        self.mus =mus
        self.sds = sds
        self.p = p

    def logp(self, value):
        """
        Implementation of Sum_i pdf(L|Ti,X) pdf(Ti|X).
        """
        taus = T.pow(self.sds,-2)
        return T.log(self.p) + (-taus[0] * (value - self.mus[0])**2 + T.log(taus[0] / numpy.pi / 2.)) / 2. + \
            T.log(1-self.p) + (-taus[1] * (value - self.mus[1])**2 + T.log(taus[1] / numpy.pi / 2.)) / 2.

class LogLuminosityGivenSpectype(Normal):
    r"""The distribution for the joint spectype and log-luminosity.

    The input variable ls log-luminosity, not luminosity.  This therefore subclasses
    Normal.

    It is the product of the probability of the type times the pdf of the luminosity
    which in this case is Lognormal (the Parent class).  Do templates exist in Python?

    .. math::
        pdf(Obs type, Luminosity | Type prob, logL_snIa, logL_snII)
            = sum_i pdf(Obs type| Type_i) *
                pdf(Luminosity | Type_i, logL_snIa, logL_snII) *
                pdf(Type_i | Type prob)
            = pdf(Luminosity | Type=Obs type, logL_snIa, logL_snII) *
                pdf(Type=Obs type | Type prob)

    This class should be generalized to handle multiple types

        
    Parameters
    -----------
    p : Theano.Tensor
        pdf(T|X)
    """

    def __init__(self, p=1, *args, **kwargs):
        super(LogLuminosityGivenSpectype, self).__init__(*args, **kwargs)
        self.p = p
#        self.logp = T.log(self.p)  #For unknown reasons this causes ths code to crash

    def logp(self, value):
        return T.log(self.p) + super(LogLuminosityGivenSpectype, self).logp(value)

class Counts(Continuous):
    r"""The distribution for the joint spectype and log-luminosity

    pdf of counts given a threshold

    .. math::
        pdf(observed redshift, Counts | Luminosity, Redshift, Cosmology, Calibration)
            = sum_i p_i pdf(Counts | Flux_i, Redshift=observer_redshift_i, Calibration)
        
    Parameters
    -----------
    threshold : theano.tensor
        minimum counts in order to be discovered

    """


    ln10 = numpy.log(10.)

    def __init__(self, fluxes=None, pzs = None, Z=None, *args, **kwargs):
        super(Counts, self).__init__(*args, **kwargs)
        self.fluxes = fluxes
        self.pzs=pzs
        self.Z = Z

    def logp(self, value):
        ans=0.
        for index in xrange(len(self.pzs)):
            flux = self.fluxes[index]
            counts  = flux*T.exp(Counts.ln10/2.5*self.Z)
            tau = 1/1e-9/1e-9
            ans = ans + T.log(self.pzs[index]) + (-tau * T.sqr(value - counts) + T.log(tau / numpy.pi / 2.)) / 2.

        #need to add threshold

        return ans
        #self.normalization  = T.log(1-std_cdf((self.threshold-self.mu)/self.sd))
        #return bound(-self.normalization+ super(CountsWithThreshold, self).logp(value), value >= self.threshold)

class SampleRenormalization(Continuous):
    r""" Renormalization factor from sample selection for a typed supernova with fixed
    redshift.

    .. math:
        P(S_c, S_T| z=z_o, X) =
            sum_i  P(T_i|z=z_o, X)
                \int dL p(L|T_i, z=z_o, X)  \left[\int_{c_T}^{\infty} dc_o  p(c_o | T=T_i, L, z=z_o, X)\right]


    """
    def __init__(self, threshold = 0, logL_snIa=0, sigma_snIa=1, logL_snII=0, sigma_snII=1,
             luminosity_distances=None, pzs=0, prob = None, Z=None, *args, **kwargs):
        super(SampleRenormalization, self).__init__()
        self.threshold = threshold
        self.logL_snIa=logL_snIa
        self.sigma_snIa = sigma_snIa
        self.logL_snII=logL_snII
        self.sigma_snII = sigma_snII
        self.pzs = pzs
        self.prob = prob
        self.luminosity_distances=luminosity_distances
        self.Z = Z

    def logp(self, value):
        # this is independent of value

        #sum over redshift
        for ld, pz, indz in zip(self.luminosity_distances,self.pzs,xrange(len(self.pzs))):
            #sum over type
            tsum = self.prob*normalization_integral(self.logL_snIa, self.sigma_snIa, self.Z, self.threshold, ld) + \
                (1-self.prob)*normalization_integral(self.logL_snII, self.sigma_snII, self.Z, self.threshold, ld)

            if indz == 0:
                ans = pz*tsum
            else: 
                ans = ans + pz*tsum
        return -T.log(ans)

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

 #   pgm.add_node(Node('Neighbors',r"\centering{Neighbor \newline Redshifts}", 5,7, scale=1.6,aspect=1.2))
    pgm.add_node(Node('Redshift',r"{Redshift}", 6,7, scale=1.6,aspect=1.2))

    pgm.add_node(Node('Type_prob',r"Type prob", 1,6, fixed=1,offset=(20,-10)))
    pgm.add_node(Node('Distance',r"$L_D$", 7,6, fixed=1,offset=(10,10)))

    pgm.add_node(Node('Type',r"Type", 1, 5, scale=1.6,aspect=1.2))

    pgm.add_node(Node('Luminosity',r"Luminosity", 4, 4, scale=1.6,aspect=1.2))
    pgm.add_node(Node('Flux',r"Flux", 7, 3, scale=1.2,fixed=True,offset=(-20,-20)))


    pgm.add_node(Node('Obs_Type',r"Obs type", 1, 1, scale=1.6,aspect=1.2,observed=1))
    pgm.add_node(Node('Obs_Redshift',r"Obs redshift", 6, 1, scale=1.6,aspect=1.2,observed=1))
    pgm.add_node(Node('Counts',r"Counts", 8, 1, scale=1.2,observed=1))


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

    pgm.add_edge("Type","Obs_Type")
#    pgm.add_edge("Neighbors","Obs_Redshift")
    pgm.add_edge("Redshift","Obs_Redshift")

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


    pgm.figure.savefig("../results/nodes_pgm.pdf")


def simulateData():
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
    observation['zprob'] = numpy.zeros(nTrans)+1.
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
    nTrans =  found.sum()
    observation['specz'] = numpy.reshape(observation['specz'][found],(nTrans,1))
    observation['zprob'] = numpy.reshape(observation['zprob'][found],(nTrans,1))
    observation['spectype'] =observation['spectype'][found]
    observation['spectype'][0] = -1   # case of no spectral type
    observation['counts'] =observation['counts'][found]
    return observation



# plt.scatter(observation['specz'],observation['counts'])
# plt.ylim((0, 2e-8))
# plt.show()

def runModel():

    observation = simulateData()
    nTrans = len(observation['spectype'])

    # Create the pymc3 model and fill it with the distributions and parameters
    # of the model
    basic_model = Model()

    with basic_model:

        r"""
        Cosmology Node.

        The FlatwCDM cosmology.  

        pdf(Om0, w0)

        We need the flexibility to switch in and out different cosmological models.  The function
        that describes luminosity distance is specific to the model: the parameters and function
        should be packaged together.

        Parameters
        ----------
        Om0:    Omega_M
        w0:     constant equation of state w
        """

        Om0 = Lognormal('Om0', mu=numpy.log(0.28), tau=1/.1/.1)
        w0 = Normal('w0', mu=-1, sd=0.05)

        """
        Calibration Node.

        Global zeropoints for each band.

        pdf(Z)

        The transmission function of the bands will be used later.  The transmission and zeropoints
        should be packaged together. More complicated parameterizations of calibration are expected.

        Parameters
        -----------
        Z:  zeropoint (in mag) for the bands

        """
        n_bands = 1
        zeropoints = Normal('zeropoints', mu=0, sd=.02, shape = n_bands)

        """
        SN Ia Rate Node.  

        rate_Ia_r = constant

        For SN cosmology the relative rates between different populations are sufficient.  Rates of
        all types are relative the snIa rate, so snIa rate is taken to be 1.

        Parameters
        -----------
        rate_Ia_r =1    : the relative rates are relative to type Ia. Fixed.

        """

        rate_Ia_r = 1.


        """
        SN II Rate Node.

        The rate of SNe II realtiave SNIa.

        pdf(rate_II_r)

        Along with the rate parameters is a rate model.

        There should be equivalent nodes for all other transient types being modeled.

        Parameters
        ----------

        rate_II_r     : relative rate of SNe II compared to SNe Ia. 

        """

        rate_II_r = Uniform('rate_II_r', lower=0.25, upper=4)

        """
        SN Ia luminosity Node.  (actually working in log-L)

        pdf(logL_snIa, sigma_snIa)

        For the moment consider the SN to be phase-indepemdent with no internal parameters.  Eventually
        this will represent time-evolving SED, e.g. SALT2.


        Parameters
        ----------

        logL_snIa   :       SN Ia mean log-luminosity
        sigma_snIa :        intrinsic dispersion (mag)

        """
        logL_snIa = Normal('logL_snIa', mu=numpy.log(1), sd = 0.02)
        sigma_snIa = Lognormal('sigma_snIa', mu=numpy.log(0.1), tau=1./0.1/0.1)

        """
        SN Ia luminosity Node.  (actually working in log-L)

        pdf(logL_snII, sigma_snIa)

        Parameters
        ----------

        logL_snII   :       SN II mean log-luminosity
        sigma_snII :        intrinsic dispersion (mag)

        """
        logL_snII = Normal('logL_snII', mu=numpy.log(0.5), sd=0.02)
        sigma_snII = Lognormal('sigma_snII', mu=numpy.log(0.4), tau=1./0.1/0.1)

        """
        Enter the plate that considers one supernova at a time
        """

        for i in xrange(nTrans):

            """
            Type Probability Node.

            Probabilities of being a type of object.  For now only SN Ia, and SN II.

            Dependencies
            -------------

            rate_Ia_r   :   Type Ia rate
            rate_II_r   :   Type II rate
            host galaxy :   Not implemented now but eventually depends on host properties

            Parameters
            ----------

            prob :          probability of the object being a type Ia.  Fixed.
            """

            prob = rate_Ia_r/(rate_Ia_r+rate_II_r)


            """
            Type Node.

            Not explicitly considered in our model.
            """

            """
            Observed Type Node and Luminosity Node.

            pdf(Obs type, Luminosity | Type prob, logL_snIa, logL_snII)

            There are two possibilities:

            1. There is an observed type assumed to be perfect.

                pdf(Obs type | Type) = delta(Obs type - Type)

                then 
                
                pdf(Obs type, Luminosity | Type prob, logL_snIa, logL_snII)
                    = sum_i pdf(Obs type| Type_i) *
                        pdf(Luminosity | Type_i, logL_snIa, logL_snII) *
                        pdf(Type_i | Type prob)
                    = pdf(Luminosity | Type=Obs type, logL_snIa, logL_snII) *
                        pdf(Type=Obs type | Type prob)

                The class LogLuminosityGivenSpectype is responsible for providing this pdf

            2. There is no observed type.

                pdf(Luminosity | Type prob, logL_snIa, logL_snII)
                    = sum_i pdf(Luminosity | Type_i, logL_snIa, logL_snII) *
                        pdf(Type_i | Type prob)

                The class LuminosityMarginalizedOverType is responsible for providing this pdf

            Dependencies
            ------------

            prob        :
            logL_snIa   :
            sigma_snIa  :
            logL_snII   :
            sigma_snII  :

            Parameters
            ----------

            obstype         :   observed type, SN Ia=0, SNII=1 Marginalized over
            Luminosity      :

            """
            if observation['spectype'][i] == -1 :
                logluminosity = LogLuminosityMarginalizedOverType('logluminosity'+str(i), 
                    mus=[logL_snIa, logL_snII], \
                    sds = [numpy.log(10)/2.5*sigma_snIa,numpy.log(10)/2.5*sigma_snII], p=prob, \
                    testval = 1.)
            else:
                if observation['spectype'][i] == 0:
                    usemu = logL_snIa
                    usesd = numpy.log(10)/2.5*sigma_snIa
                    usep = prob
                else:
                    usemu = logL_snII
                    usesd = numpy.log(10)/2.5*sigma_snII
                    usep = 1-prob

                logluminosity = LogLuminosityGivenSpectype('logluminosity'+str(i), \
                        mu=usemu,sd=usesd, p=usep)
                
            luminosity = T.exp(logluminosity)

            """
            Redshift Node.

            Not considered explicitly in our model.

            """

            """
            Observed Redshift, Counts Node.

            pdf(observed redshift, Counts | Luminosity, Redshift, Cosmology, Calibration)
                = pdf(observed redshift| Redshift) *
                    pdf(Counts | Luminosity, Redshift, Cosmology, Calibration)

            The pdf of the observed redshift is assumed to be a sum of delta functions, perfectly
            measured redshift of the supernova or redshifts of potential galaxy hosts.

            pdf(observed redshift | Redshift) = sum_i p_i delta(observer redshift_i - Redshift)

            where p_i is the probability of observer redshift_i being the correct redshift.

            so

            pdf(observed redshift, Counts | Luminosity, Redshift, Cosmology, Calibration)
                = sum_i p_i pdf(Counts | Luminosity, Redshift=observer_redshift_i, Cosmology, Calibration)

            The class CountsWithThreshold handles this pdf

            Dependencies
            ------------

            luminosity  :   luminosity
            redshift    :   host redshift
            cosmology   :   cosmology
            Calibration :   calibration

            Parameters
            -----------

            observed_redshift   Marginalized over
            counts

            """

            lds=[]
            fluxes=[]
            for z_ in observation['specz'][i]:
                # ld = 0.5/h0*(z_+T.sqr(z_))* \
                #     (1+ 1//T.sqrt((1+z_)**3 * (Om0 + (1-Om0)*(1+z_)**(3*w0))))
                ld = luminosity_distance(z_, Om0, w0)
                lds.append(ld)
                fluxes.append(luminosity/4/numpy.pi/ld**2)

            counts = Counts('counts'+str(i),fluxes =fluxes,  \
                pzs = observation['zprob'][i], Z=zeropoints, observed=observation['counts'][i])

            if observation['spectype'][i] == -1 :
                pass
            else:
                normalization=SampleRenormalization('normalization'+str(i), threshold = 1e-9, 
                    logL_snIa=logL_snIa, sigma_snIa=sigma_snIa, logL_snII=logL_snII, sigma_snII=sigma_snII,
                    luminosity_distances=lds, Z=zeropoints, pzs=observation['zprob'][i], prob=prob, observed=1)

    from pymc3 import find_MAP, NUTS, sample, summary
    from scipy import optimize
    with basic_model:

        backend = SQLite('trace.sqlite')

        # obtain starting values via MAP
        start = find_MAP(fmin=optimize.fmin_bfgs, disp=True)

        # draw 2000 posterior samples
        trace = sample(500, start=start, trace=backend)

        summary(trace)
runModel()
#pgm()