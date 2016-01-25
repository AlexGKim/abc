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
from theano import pp

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
        return bound(T.log(self.p)+ super(LogLuminosityGivenSpectype, self).logp(value), self.p > 0, \
            self.p<=1)

class CountsWithThreshold(Continuous):
    r"""The distribution for the joint spectype and log-luminosity

    pdf of counts given a threshold

    .. math::

    This class should be generalized to handle multiple types

        
    Parameters
    -----------
    threshold : theano.tensor
        minimum counts in order to be discovered

    """

    cosmo = FlatwCDM(H0=72, Om0=0.28, w0=-1)
    h0 = (const.c/cosmo.H0).to(u.Mpc).value

    def __init__(self, threshold=0, luminosity=None, zs=None, pzs = None, Om0= None, w0=None, \
            Z=None, *args, **kwargs):
        super(CountsWithThreshold, self).__init__(*args, **kwargs)
        self.threshold = threshold
        self.zs=zs
        self.pzs=pzs
        self.luminosity = luminosity
        self.Om0 = Om0
        self.w0 = w0
        self.Z = Z
        
    def h_inv(self,z):
        # return 1./T.sqrt(self.Om0*(1+z)**3 + (1-self.Om0)*(1+z)**(3*(1+self.w0)))
        return (1+z)**(-1.5)/T.sqrt(self.Om0 + (1-self.Om0)*(1+z)**(3*self.w0))
        
    def luminosity_distance(self, z):
        """
        Poorman luminosity distance.

        Presumably there will be a numerical integration function that inherits from theano.Op with
        gradient implmented so that HMC can be run.  The class can be specified by the integrand,
        which is the gradient.
        """

        # luminosity_distance = (1+z)*(0.5 \
        #     +0.5/T.sqrt(self.Om0*(1+z)**3 + (1-self.Om0)*(1+z)**(3*(1+self.w0))) \
        #     +1./T.sqrt(self.Om0*(1+.25*z)**3 + (1-self.Om0)*(1+.25*z)**(3*(1+self.w0))) \
        #     +1./T.sqrt(self.Om0*(1+.50*z)**3 + (1-self.Om0)*(1+.50*z)**(3*(1+self.w0))) \
        #     +1./T.sqrt(self.Om0*(1+.75*z)**3 + (1-self.Om0)*(1+.75*z)**(3*(1+self.w0))) \
        #     )/CountsWithThreshold.h0*z/4
        luminosity_distance = 0.5/CountsWithThreshold.h0*(z+z*2)*(1+ self.h_inv(z))
        return luminosity_distance

    def logp(self, value):
        
        for index in xrange(len(self.zs)):
            ld = self.luminosity_distance(self.zs[index])
            mu = self.luminosity/(4*numpy.pi)/(ld**2)*10**(-self.Z/2.5)    #average counts
            tau = (0.02*mu)**(-2)
            if index == 0:
                ans = bound(T.log(self.pzs[index]) + (-tau * (value - mu)**2 + T.log(tau / numpy.pi / 2.)) / 2., \
                    tau > 0,self.pzs[index] >0,self.pzs[index] <=1 )
            else:
                ans = bound(ans + T.log(self.pzs[index]) + (-tau * (value - mu)**2 + T.log(tau / numpy.pi / 2.)) / 2., \
                    tau > 0,self.pzs[index] >0,self.pzs[index] <=1 )
        #need to add threshold

        return ans
        #self.normalization  = T.log(1-std_cdf((self.threshold-self.mu)/self.sd))
        #return bound(-self.normalization+ super(CountsWithThreshold, self).logp(value), value >= self.threshold)


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

#   pgm.add_node(Node('Neighbors',r"\centering{Neighbor \newline Redshift}", 5,7, scale=1.6,aspect=1.2))
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


    pgm.figure.savefig("nodes_pgm.pdf")


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

        Om0 = Lognormal('Om0', mu=numpy.log(0.28), tau=1/.2/.2)
        w0 = Normal('w0', mu=-1, sd=0.2)

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

        rate_snII_r     : relative rate of SNe II compared to SNe Ia. 

        """
        rate_II_r = Uniform('rate_II_r', lower=0.1, upper=10)

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
        sigma_snIa = Lognormal('sigma_snIa', mu=numpy.log(0.02), tau=1./0.1/0.1)

        """
        SN Ia luminosity Node.  (actually working in log-L)

        pdf(logL_snII, sigma_snIa)

        Parameters
        ----------

        logL_snII   :       SN II mean log-luminosity
        sigma_snII :        intrinsic dispersion (mag)

        """
        logL_snII = Normal('logL_snII', mu=numpy.log(0.5), sd=0.02)
        sigma_snII = Lognormal('sigma_snII', mu=numpy.log(0.2), tau=1./0.4/0.4)

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

            obstype         :   observed type, SN Ia=0, SNII=1
            Luminosity      :

            """
            if observation['spectype'][i] is not None:

                if observation['spectype'][i] == 0:
                    logluminosity = LogLuminosityGivenSpectype('logluminosity'+str(i),
                        mu=logL_snIa,sd=sigma_snIa, p=prob)
                else:
                    logluminosity = LogLuminosityGivenSpectype('logluminosity'+str(i),
                        mu=logL_snII,sd=sigma_snII, p=1-prob)

                luminosity = T.exp(logluminosity)

            else:

                logluminosity = LuminosityMarginalizedOverType(
                    Normal('dum1'+str(i), mu=logL_snIa,sd=sigma_snIa),
                    Normal('dum2'+str(i), mu=logL_snII,sd=sigma_snII), prob)
                luminosity = T.exp(logluminosity)

            """
            Redshift Node.

            Not considered explicitly in our model.

            """

            """
            Observed Redshift, Counts Node.

            pdf(observed redshift, Counts | Luminosity, Redshift, Cosmology, Calibration)
                = pdf(observed redshift| Redshift) *
                    pdf(Counts | Redshift, Cosmology, Calibration)

            The pdf of the observed redshift is assumed to be a sum of delta functions, prefectly
            measured redshift of the supernova or redshifts of potential galaxy hosts.

            pdf(observed redshift | Redshift) = sum_i p_i delta(observer redshift_i - Redshift)

            where p_i is the probability of observer redshift_i being the correct redshift.

            so

            pdf(observed redshift, Counts | Luminosity, Redshift, Cosmology, Calibration)
                = sum_i p_i pdf(Counts | Redshift=observer_redshift_i, Cosmology, Calibration)

            The class CountsWithThreshold handles this pdf


            Dependencies
            ------------

            luminosity  :   luminosity
            redshift    :   host redshift
            cosmology   :   cosmology
            Calibration :   calibration

            Parameters
            -----------

            observed_redshift   
            counts

            """

            counts = CountsWithThreshold('counts'+str(i),luminosity=luminosity, \
                zs=observation['specz'][i],  pzs = observation['zprob'][i], Om0= Om0, \
                w0=w0, Z=zeropoints, observed=observation['counts'][i])


    from pymc3 import find_MAP, NUTS, sample
    from scipy import optimize
    with basic_model:

        # obtain starting values via MAP
        start = find_MAP(fmin=optimize.fmin_powell)

        # draw 2000 posterior samples
        trace = sample(500, start=start)

runModel()
#pgm()