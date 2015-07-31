#!/usr/bin/env python

import pystan
import matplotlib.pyplot as plt
import numpy.random
from astropy.cosmology import WMAP9 as cosmo

N_sn=100
N_s_obs=90
zs = numpy.sort(numpy.random.uniform(0.1,1,N_sn))
s_obs = numpy.sort(numpy.random.choice(N_sn, N_s_obs,replace=False))
s_mis = numpy.zeros(N_sn,dtype=bool)
s_mis[s_obs] = True
s_mis = numpy.where(numpy.logical_not(s_mis))[0]

sigma_snIa=0.1
sigma_nonIa=.5

alpha_snIa=2.
alpha_nonIa=1.
frac_Ia=0.8

snIa = numpy.random.binomial(1, frac_Ia, size=N_sn)

wsnIa = numpy.where(snIa)[0]
wnonIa = numpy.where(numpy.logical_not(snIa))[0]

adu = 1/(cosmo.luminosity_distance(zs).value/cosmo.hubble_distance.value)**2
adu[wsnIa]*= alpha_snIa*adu[wsnIa]*10**(numpy.random.normal(loc=0, scale=sigma_snIa,size=len(wsnIa))/2.5)
adu[wnonIa]*= alpha_nonIa*adu[wnonIa]*10**(numpy.random.normal(loc=0, scale=sigma_nonIa,size=len(wnonIa))/2.5)


data = {'N_sn':N_sn,
		'N_obs':N_s_obs,
		'N_mis':N_sn-N_s_obs,

		'N_adu_max':1,

		'zmin':0.1,
		'zmax':2.0,

		'adu_obs': adu[s_obs,None],
		'adu_mis': adu[s_mis,None],

		'trans_zs_obs': zs[s_obs],
		'snIa_obs': snIa[s_obs],

		'host_zs_obs': zs[s_obs],
		'host_zs_mis': zs[s_mis]
		}

fit1 = pystan.stan(file='des.stan', data=data, iter=1000, chains=4)

print fit1
fit1.plot()
plt.show()
