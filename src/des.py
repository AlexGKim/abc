#!/usr/bin/env python

import pystan
import matplotlib.pyplot as plt
import numpy.random

N_sn=50
N_s_obs=45
zs = numpy.sort(numpy.random.uniform(0.1,1,N_sn))
s_obs = numpy.sort(numpy.random.choice(N_sn, N_s_obs,replace=False))
s_mis = numpy.zeros(N_sn,dtype=bool)
s_mis[s_obs] = True
s_mis = numpy.where(numpy.logical_not(s_mis))[0]

sigma=0.02

alpha_snIa=0.
alpha_nonIa=1.
frac_Ia=0.9

snIa = numpy.random.binomial(1, frac_Ia, size=N_sn)

wsnIa = numpy.where(snIa)[0]
adu = zs + alpha_nonIa + numpy.random.normal(loc=0, scale=sigma,size=N_sn)
adu[wsnIa]+= (alpha_snIa-alpha_nonIa)


data = {'N_sn':N_sn,
		'N_s_obs':N_s_obs,
		'N_s_mis':N_sn-N_s_obs,
		'N_adu_max':1,
		'adu_s_obs': adu[s_obs,None],
		'adu_s_mis': adu[s_mis,None],
		'zs_obs': zs[s_obs],
		'snIa_obs': snIa[s_obs]
		}

fit1 = pystan.stan(file='des.stan', data=data, iter=1000, chains=4)

print fit1
fit1.plot()
plt.show()
