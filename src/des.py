#!/usr/bin/env python

import pystan
import matplotlib.pyplot as plt
import numpy.random
from astropy.cosmology import FlatLambdaCDM

omega_M=0.28
cosmo=FlatLambdaCDM(70,omega_M)

N_sn=50
N_s_obs=50
zmin=0.1
zmax=1.
zs = numpy.sort(numpy.random.uniform(zmin,zmax,N_sn))
s_obs = numpy.sort(numpy.random.choice(N_sn, N_s_obs,replace=False))
s_mis = numpy.zeros(N_sn,dtype=bool)
s_mis[s_obs] = True
s_mis = numpy.where(numpy.logical_not(s_mis))[0]

sigma_snIa=0.1
sigma_nonIa=.3

alpha_snIa=2.
alpha_nonIa=1.
frac_Ia=.8

snIa = numpy.random.binomial(1, frac_Ia, size=N_sn)

wsnIa = numpy.where(snIa)[0]
wnonIa = numpy.where(numpy.logical_not(snIa))[0]

adu = 1/(cosmo.luminosity_distance(zs).value/cosmo.hubble_distance.value)**2
adu[wsnIa] = alpha_snIa*adu[wsnIa]*10**(numpy.random.normal(loc=0, scale=sigma_snIa/2.5,size=len(wsnIa)))
adu[wnonIa] = alpha_nonIa*adu[wnonIa]*10**(numpy.random.normal(loc=0, scale=sigma_nonIa/2.5,size=len(wnonIa)))

data = {'N_sn':N_sn,
		'N_obs':N_s_obs,

		'N_adu_max':1,

		'zmin':zmin,
		'zmax':zmax,

		'adu_obs': adu[s_obs,None],
		'adu_mis': adu[s_mis,None],

		'trans_zs_obs': zs[s_obs],
		'snIa_obs': snIa[s_obs],

		'host_zs_obs': zs[s_obs],
		'host_zs_mis': zs[s_mis]
		}

init = {
	'Omega_M':omega_M,
	'Omega_L':1-omega_M,
	'zs_true_obs': zs[s_obs],
  	'zs_true_mis': zs[s_mis],
  	'alpha_Ia': alpha_snIa,
 	'alpha_nonIa': alpha_nonIa,
  	'sigma_Ia': sigma_snIa,
  	'sigma_nonIa':sigma_nonIa,
  	'snIa_rate':frac_Ia
	}

fit1 = pystan.stan(file='des.stan', data=data, iter=1000, chains=4, init=[init,init,init,init])

print fit1
fit1.plot()
plt.show()
