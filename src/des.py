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
	host_zs_random_ = zs*host_choice + neighbor_zs*(1-host_choice)
	neighbor_zs_random_ = zs*(1-host_choice) + neighbor_zs*host_choice
	# wrong_host_zs : indeces if incorrect host galaxy association
	#wrong_host_zs = host_zs_random_ == 0
	#host_zs_random_[wrong_host_zs] = 

	# host_zs_random : redshifts of potential host galaxies with the probability of each
	host_zs_random=[]
	for i in xrange(N_sn):
		ans=[]
		ans.append(numpy.array([host_zs_random_[i],numpy.log(0.98/(1-0.98)]))
		ans.append(numpy.array([neighbor_zs_random_[i],numpy.log(0.02/(1-0.02)]))
		host_zs_random.append(ans)

	# s_obs_random : order in which supernovae get a spectrum
	s_obs_random = numpy.arange(N_sn,dtype='int')
	numpy.random.shuffle(s_obs_random)
	s_obs = numpy.sort(s_obs_random[:N_s_obs]).tolist()
	s_mis = numpy.sort(s_obs_random[N_s_obs:]).tolist()

	host_zs_obs=[]
	for i in xrange(len(s_obs)):
		host_zs_obs.append(host_zs_random[s_obs[i]])
	host_zs_mis=[]
	for i in xrange(len(s_mis)):
		host_zs_mis.append(host_zs_random[s_mis[i]])

	host_zs_mis=numpy.array(host_zs_mis)

	data = {'N_sn':N_sn,
			'N_obs':N_s_obs,

			'N_adu_max':1,

			'zmin':zmin,
			'zmax':zmax,

			'adu_obs': adu[s_obs,None],
			'adu_mis': adu[s_mis,None],

			'trans_zs_obs': zs[s_obs],
			'snIa_obs': snIa[s_obs],

			'host_zs_obs': host_zs_obs, #zs[s_obs],
			'host_zs_mis_': host_zs_mis.flatten(),

			'n_int': 50
			}

	init=[]
	for i in xrange(Ninit):
		# host_zs_mis_init : initial conditions of host redshift _mis for MCMC
		host_choice = numpy.random.binomial(1,0.98,size=N_sn)
		host_zs_mis_init = (host_zs_random_*host_choice + neighbor_zs_random_*(1-host_choice))[s_mis]
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
		init.append ( {
			'Omega_M':numpy.random.normal(omega_M,0.1),
			'Omega_L':1-omega_M,
			'w': numpy.random.normal(-0.9,0.1),
			'ainv_true_obs': 1+zs[s_obs],
		  	'ainv_true_mis': numpy.random.uniform((1+zmin/1.5)**3,(1+zmax*1.5)**3,size=N_sn-N_s_obs)**(1./3), #host_zs_mis_init,
		  	'alpha_Ia': numpy.random.normal(alpha_snIa,0.1),
		 	'alpha_nonIa': numpy.random.normal(alpha_nonIa,0.1),
		  	'sigma_Ia': numpy.random.normal(sigma_snIa,0.05),
		  	'sigma_nonIa':numpy.random.normal(sigma_nonIa,0.5),
		  	'snIa_rate':numpy.random.lognormal(numpy.log(frac_Ia),0.1)
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
	N_sn=200

	sm = pystan.StanModel(file='des.stan')

	nspec = numpy.arange(100,N_sn+1,250)
	mn=[]
	std=[]
	for ns in nspec:
		data, init, info = genData(N_sn,ns,Nchains,1)

		fit = sm.sampling(data=data, iter=1000,  chains=Nchains, init=init)
		samples = fit.extract(['Omega_M','ainv_true_mis','w'])

		logposterior = fit.get_logposterior()

		with open('model'+str(ns)+'.pkl', 'wb') as f:
			pickle.dump([fit.extract(),info, logposterior], f)


if __name__ == "__main__":
    main()