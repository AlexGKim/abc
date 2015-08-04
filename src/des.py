#!/usr/bin/env python

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pystan
import numpy.random
from astropy.cosmology import FlatLambdaCDM

import pickle


def genData(N_sn, N_s_obs, Ninit, seed):

	numpy.random.seed(seed)

	omega_M=0.28
	cosmo=FlatLambdaCDM(70,omega_M)

	zmin=0.1
	zmax=1.

	sigma_snIa=0.1
	sigma_nonIa=1

	alpha_snIa=2.
	alpha_nonIa=1.
	frac_Ia=.8

	zs = numpy.sort(numpy.random.uniform(zmin,zmax,N_sn))

	snIa = numpy.random.binomial(1, frac_Ia, size=N_sn)

	wsnIa = numpy.where(snIa)[0]
	wnonIa = numpy.where(numpy.logical_not(snIa))[0]

	adu = 1/(cosmo.luminosity_distance(zs).value/cosmo.hubble_distance.value)**2

	adu_random = numpy.random.normal(size=N_sn)

	adu[wsnIa] = alpha_snIa*adu[wsnIa]*10**(adu_random[wsnIa]*sigma_snIa/2.5)
	adu[wnonIa] = alpha_nonIa*adu[wnonIa]*10**(adu_random[wnonIa]*sigma_nonIa/2.5)

	host_zs_random = zs*numpy.random.binomial(1,0.98,size=N_sn)
	wrong_host_zs = host_zs_random == 0
	host_zs_random[wrong_host_zs] = 	numpy.random.uniform(zmin/1.5,zmax*1.5,size=wrong_host_zs.sum())

	numpy.random.seed()

	s_obs = numpy.sort(numpy.random.choice(N_sn, N_s_obs,replace=False))
	s_mis = numpy.zeros(N_sn,dtype=bool)
	s_mis[s_obs] = True
	s_mis = numpy.where(numpy.logical_not(s_mis))[0]

	host_zs_mis=host_zs_random[s_mis]

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
			'host_zs_mis': host_zs_mis
			}

	init=[]
	for i in xrange(Ninit):
		host_zs_mis = zs[s_mis]*numpy.random.binomial(1,0.98,size=len(s_mis))
		wrong_host_zs = host_zs_mis == 0
		host_zs_mis[wrong_host_zs] = numpy.random.uniform(zmin/1.5,zmax*1.5,size=wrong_host_zs.sum())
		init.append ( {
			'Omega_M':omega_M,
			'Omega_L':1-omega_M,
			'w': -1,
			'zs_true_obs': zs[s_obs],
		  	'zs_true_mis': host_zs_mis,
		  	'alpha_Ia': alpha_snIa,
		 	'alpha_nonIa': alpha_nonIa,
		  	'sigma_Ia': sigma_snIa,
		  	'sigma_nonIa':sigma_nonIa,
		  	'snIa_rate':frac_Ia
			} )

	return data, init

Nchains=16

sm = pystan.StanModel(file='des.stan')

nspec = numpy.arange(500,1001,100)
mn=[]
std=[]
for ns in nspec:
	data, init = genData(1000,ns,Nchains,1)

	fit = sm.sampling(data=data, iter=1000, chains=Nchains, init=init)
	print fit
	samples = fit.extract(['Omega_M','zs_true_mis','w'])
	mn.append(samples['w'].mean())
	std.append(samples['w'].std())

	with open('model'+str(ns)+'.pkl', 'wb') as f:
		pickle.dump(fit.extract(), f)

# plt.errorbar(nspec,mn,yerr=std,fmt='o')
# plt.savefig('plot.pdf')


# fit1 = pystan.stan(file='des.stan', data=data, iter=1000, chains=4, init=[init,init,init,init])

# print fit1
# fit1.plot()
# plt.show()
