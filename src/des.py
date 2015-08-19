#!/usr/bin/env python

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pystan
import numpy.random
from astropy.cosmology import FlatLambdaCDM

import pickle


def genData(N_sn, N_s_obs, Ninit, seed, ia_only=False):

	numpy.random.seed(seed)

	omega_M=0.28
	cosmo=FlatLambdaCDM(70,omega_M)

	zmin=0.1
	zmax=1.

	sigma_snIa=0.1
	sigma_nonIa=1

	alpha_snIa=2.
	alpha_nonIa=alpha_snIa*10**(-2./2.5)
	frac_Ia=.8

	# zs : the true redshifts
	zs = numpy.sort(numpy.random.uniform(zmin,zmax,N_sn))

	# snIa : the true types
	snIa = numpy.random.binomial(1, frac_Ia, size=N_sn)

	wsnIa = numpy.where(snIa)[0]
	wnonIa = numpy.where(numpy.logical_not(snIa))[0]

	# adu : the observed counts
	adu = 1/(cosmo.luminosity_distance(zs).value/cosmo.hubble_distance.value)**2
	adu_random = numpy.random.normal(size=N_sn)
	adu[wsnIa] = alpha_snIa*adu[wsnIa]*10**(adu_random[wsnIa]*sigma_snIa/2.5)
	adu[wnonIa] = alpha_nonIa*adu[wnonIa]*10**(adu_random[wnonIa]*sigma_nonIa/2.5)

	# host_choice : 1 if correct host is chosen to have highest probability, 0 otherwise
	host_choice = numpy.random.binomial(1,0.98,size=N_sn)
	neighbor_zs=numpy.random.uniform((1+zmin/1.5)**3,(1+zmax*1.5)**3,size=N_sn)**(1./3)-1
	host_zs_random = zs*host_choice + neighbor_zs*(1-host_choice)
	neighbor_zs_random = zs*(1-host_choice) + neighbor_zs*host_choice
	# wrong_host_zs : indeces if incorrect host galaxy association
	#wrong_host_zs = host_zs_random_ == 0
	#host_zs_random_[wrong_host_zs] = 


	# s_obs_random : order in which supernovae get a spectrum
	s_obs_random = numpy.arange(N_sn,dtype='int')
	numpy.random.shuffle(s_obs_random)

	s_obs=numpy.sort(s_obs_random[:N_s_obs])
	s_mis= numpy.sort(s_obs_random[N_s_obs:])



	if ia_only:
		s_snIa = snIa[s_obs] ==1 
		s_obs = s_obs[s_snIa]
		s_snIa = snIa[s_mis] ==1 
		s_mis = s_mis[s_snIa]
		N_sn=snIa.sum()
		N_s_obs = len(s_obs)


	s_obs = s_obs.tolist()
	s_mis = s_mis.tolist()

	data = {'N_sn':N_sn,
			'N_obs':N_s_obs,

			'N_adu_max':1,

			'zmin':zmin,
			'zmax':zmax,

			'adu_obs': adu[s_obs,None],
			'adu_mis': adu[s_mis,None],

			'trans_ainv_obs': 1+zs[s_obs],
			'snIa_obs': snIa[s_obs],

			'host_zs_obs': host_zs_random[s_obs], #zs[s_obs],
			'host_zs_mis': host_zs_random[s_mis],
			'host2_zs_mis': neighbor_zs_random[s_mis],

			'n_int': 25,
			'N_SNIa': snIa[s_obs].sum()
			}

	init=[]
	for i in xrange(Ninit):
		# host_zs_mis_init : initial conditions of host redshift _mis for MCMC
		#host_choice = numpy.random.binomial(1,0.98,size=N_sn)
		#host_zs_mis_init = (host_zs_random_*host_choice + neighbor_zs_random_*(1-host_choice))[s_mis]
		# neighbor_zs_random__ = host_zs*_random_(1-host_choice) + neighbor_zs_random_*host_choice
		# host_zs_mis_init=[]
		# for j in xrange(len(s_mis)):
		# 	z_cand = numpy.array([host_zs_random[s_mis[j]][0][0],host_zs_random[s_mis[j]][1][0],host_zs_random[s_mis[j]][2][0]])
		# 	ind = numpy.random.multinomial(1,[host_zs_random[s_mis[j]][0][1],host_zs_random[s_mis[j]][1][1],host_zs_random[s_mis[j]][2][1]])
		# 	ind = numpy.where(ind ==1)[0][0]
		# 	host_zs_mis_init.append(z_cand[ind])
		# host_zs_mis_init = numpy.array(host_zs_mis_init)

#		host_zs_mis_init = host_zs_random[s_mis]*numpy.random.binomial(1,0.98,size=len(s_mis))
#		wrong_host_zs_init = host_zs_mis_init == 0
#		host_zs_mis_init[wrong_host_zs_init] = numpy.random.uniform(zmin/1.5,zmax*1.5,size=wrong_host_zs_init.sum())
		rate = numpy.random.lognormal(numpy.log(frac_Ia),0.1)
		init.append ( {
			'Omega_M':numpy.random.normal(omega_M,0.1),
			'Omega_L':1-omega_M,
			'w': numpy.random.normal(-0.9,0.1),
			'ainv_true_obs': 1+zs[s_obs],
		  	'ainv_true_mis': (i % 2)*(1+host_zs_random[s_mis]) + (1-(i%2))*(1+zs[s_mis]), #host_zs_mis_init,
		  	'alpha_Ia': numpy.random.normal(alpha_snIa,0.1),
		 	'alpha_nonIa': numpy.random.normal(alpha_nonIa,0.1),
		  	'sigma_Ia': numpy.random.normal(sigma_snIa,0.05),
		  	'sigma_nonIa':numpy.random.normal(sigma_nonIa,0.5),
		  	'snIa_rate': [rate,1-rate]
			} )

	info = {
			'zs' : zs,							# the true redshifts
			'snIa' : snIa,						# the true types
			'host_zs_random' : host_zs_random, 	# redshift of associated host galaxy IF there were no spectra
			'host_choice' : host_choice,	# indeces if incorrect host galaxy association
			's_obs' : s_obs,
			's_mis' : s_mis
		}

	return data, init, info

def main():
	Nchains=4
	N_sn=500
	ia_only=False
	sm = pystan.StanModel(file='des.stan')

	nspec = numpy.arange(300,N_sn+1,100)
	mn=[]
	std=[]
	for ns in nspec:
		data, init, info = genData(N_sn,ns,Nchains,1, ia_only=ia_only)

		fit = sm.sampling(data=data, iter=50000, thin=5, chains=Nchains, init=init)
		#print fit
		#samples = fit.extract(['Omega_M','ainv_true_mis','w'])

		logposterior = fit.get_logposterior()

		app=''
		if ia_only:
			app+='.ia_only.'
		with open('../results/model'+app+str(ns)+'.pkl', 'wb') as f:
			pickle.dump([fit.extract(),info, logposterior], f)


if __name__ == "__main__":
    main()
