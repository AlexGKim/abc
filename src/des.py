#!/usr/bin/env python

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pystan
import numpy.random
from astropy.cosmology import FlatLambdaCDM
import scipy

import pickle

class Data(object):
	"""docstring for Data"""
	def __init__(self, N_sn, seed, ia_only=False):
		super(Data, self).__init__()
		self.N_sn = N_sn
		self.seed = seed
		self.ia_only = ia_only

		self.omega_M=0.28
		self.cosmo=FlatLambdaCDM(70,self.omega_M)

		self.zmin=0.1
		self.zmax=1.

		self.sigma_snIa=0.1
		self.sigma_nonIa=1

		self.alpha_snIa=2.
		self.alpha_nonIa=self.alpha_snIa*10**(-2./2.5)
		self.frac_Ia_0=.5
		self.frac_Ia_1=.1
		numpy.random.seed(seed)
		self.initialize_()

	def initialize_(self):
		self.redshifts_()
		self.types_()
		self.adus_()
		self.hosts_()

	def redshifts_(self):
		# zs : the true redshifts
		volume_zero = self.cosmo.comoving_volume(self.zmin).value
		volume_norm = self.cosmo.comoving_volume(self.zmax).value- volume_zero

		self.zs = numpy.sort(numpy.random.uniform(size=self.N_sn))
		for i in xrange(len(self.zs)):
			self.zs[i]=scipy.optimize.newton(lambda z: (self.cosmo.comoving_volume(z).value - volume_zero)/volume_norm-self.zs[i], 0.5)		

	def types_(self):
		self.snIa = numpy.random.uniform(size=self.N_sn)
		self.snIa = numpy.less_equal(self.snIa, self.frac_Ia_0 + (self.frac_Ia_1-self.frac_Ia_0)*self.zs/self.zmax/1.5).astype(int)

	def adus_(self):
		adu = 1/(self.cosmo.luminosity_distance(self.zs).value/self.cosmo.hubble_distance.value)**2
		adu_random = numpy.random.normal(size=self.N_sn)
		wsnIa = numpy.where(self.snIa)[0]
		wnonIa = numpy.where(numpy.logical_not(self.snIa))[0]
		adu[wsnIa] = self.alpha_snIa*adu[wsnIa]*10**(adu_random[wsnIa]*self.sigma_snIa/2.5)
		adu[wnonIa] = self.alpha_nonIa*adu[wnonIa]*10**(adu_random[wnonIa]*self.sigma_nonIa/2.5)
		self.adu = adu

	def hosts_(self):
		host_choice = numpy.random.binomial(1,0.98,size=self.N_sn)
		self.neighbor_zs=numpy.random.uniform((1+self.zmin/1.5)**3,(1+self.zmax*1.5)**3,size=self.N_sn)**(1./3)-1
		self.host_zs_random = self.zs*host_choice + self.neighbor_zs*(1-host_choice)
		self.neighbor_zs_random = self.zs*(1-host_choice) + self.neighbor_zs*host_choice

	def found(self, ADU0):
		self.ADU0 =ADU0
		self.found = self.adu>=ADU0
		self.s_obs_random = numpy.where(self.found)[0]
		numpy.random.shuffle(self.s_obs_random)

	def spectrum(self, frac_obs):
		self.N_s_obs = frac_obs*self.found.sum()
		s_obs=numpy.sort(self.s_obs_random[:self.N_s_obs])
		s_mis= numpy.sort(self.s_obs_random[self.N_s_obs:])
		self.s_obs = s_obs.tolist()
		self.s_mis = s_mis.tolist()

	def observe(self, ADU0, frac_obs):
		data.found(ADU0)
		data.spectrum(frac_obs)

	def dict(self):
		return 	{'N_sn': self.found.sum(),
			'N_obs': len(self.s_obs),

			'N_adu_max':1,

			'zmin':self.zmin,
			'zmax':self.zmax,

			'adu_obs': self.adu[self.s_obs,None],
			'adu_mis': self.adu[self.s_mis,None],

			'trans_ainv_obs': 1+self.zs[self.s_obs],
			'snIa_obs': self.snIa[self.s_obs],

			'host_zs_obs': self.host_zs_random[self.s_obs], #zs[s_obs],
			'host_zs_mis': self.host_zs_random[self.s_mis],
			'host2_zs_mis': self.neighbor_zs_random[self.s_mis],

			'ADU0': self.ADU0,
#			'n_int': 25,
			'N_SNIa': self.snIa[self.s_obs].sum()
			}

	def init(self, n):
		ans=[]
		for i in xrange(n):
			rate = numpy.random.lognormal(numpy.log(self.frac_Ia_0),0.1)
			rate1 = numpy.random.lognormal(numpy.log(self.frac_Ia_1),0.1)
			ans.append( {
				'Omega_M':numpy.random.normal(self.omega_M,0.1),
				'Omega_L':1-self.omega_M,
				'w': numpy.random.normal(-0.9,0.1),
				# 'ainv_true_obs': 1+zs[s_obs],
			 #  	'ainv_true_mis': (i % 2)*(1+host_zs_random[s_mis]) + (1-(i%2))*(1+zs[s_mis]), #host_zs_mis_init,
			  	'alpha_Ia': numpy.random.normal(self.alpha_snIa,0.1),
			 	'alpha_nonIa': numpy.random.normal(self.alpha_nonIa,0.1),
			  	'sigma_Ia': numpy.random.lognormal(numpy.log(self.sigma_snIa),0.05),
			  	'sigma_nonIa':numpy.random.lognormal(numpy.log(self.sigma_nonIa),0.5),
			  	'snIa_rate_0': [rate,1-rate],
			  	'snIa_rate_1': [rate1,1-rate1]
				})
		return ans

	def plot(self):
		with PdfPages('foo.pdf') as pdf:
		     # As many times as you like, create a figure fig and save it:
		     fig = plt.figure()
		     f = plt.hist([self.zs[numpy.logical_and(self.found, self.snIa ==1)], self.zs[numpy.logical_and(self.found, self.snIa ==0)]],label=['SN Ia','non-Ia'])
		     plt.legend()
		     pdf.savefig(fig)
		     fig = plt.figure()
		     plt.scatter(self.zs[numpy.logical_and(self.found, self.snIa ==1)], -2.5*numpy.log10(self.adu[numpy.logical_and(self.found, self.snIa ==1)]),label='SN Ia',color='b')
		     plt.scatter(self.zs[numpy.logical_and(self.found, self.snIa ==0)], -2.5*numpy.log10(self.adu[numpy.logical_and(self.found, self.snIa ==0)]),label='non-Ia',color='r')
		     plt.scatter(self.zs[self.snIa ==1], -2.5*numpy.log10(self.adu[self.snIa ==1]),alpha=0.1,color='b')
		     plt.scatter(self.zs[self.snIa ==0], -2.5*numpy.log10(self.adu[self.snIa ==0]),alpha=0.1,color='r')
		     plt.scatter(self.zs[self.s_obs_random[:self.N_s_obs]], -2.5*numpy.log10(self.adu[self.s_obs_random[:self.N_s_obs]]),marker='+',color='k',label='Spectrum',s=40)
		     ax=plt.gca()
		     ax.invert_yaxis()
		     plt.legend()
		     pdf.savefig(fig)



# def genData(N_sn, frac_obs, Ninit, seed, ia_only=False, ADU0=None):

# 	numpy.random.seed(seed)

# 	omega_M=0.28
# 	cosmo=FlatLambdaCDM(70,omega_M)

# 	zmin=0.1
# 	zmax=1.

# 	sigma_snIa=0.1
# 	sigma_nonIa=1

# 	alpha_snIa=2.
# 	alpha_nonIa=alpha_snIa*10**(-2./2.5)
# 	frac_Ia_0=.5
# 	frac_Ia_1=.1

# 	# zs : the true redshifts
# 	volume_zero = cosmo.comoving_volume(zmin).value
# 	volume_norm = cosmo.comoving_volume(zmax).value- volume_zero

# 	cosmo.differential_comoving_volume
# 	zs = numpy.sort(numpy.random.uniform(size=N_sn))
# 	for i in xrange(len(zs)):
# 		zs[i]=scipy.optimize.newton(lambda z: (cosmo.comoving_volume(z).value - volume_zero)/volume_norm-zs[i], 0.5)

# 	# snIa : the true types
# #	snIa = numpy.random.binomial(1, frac_Ia, size=N_sn)	
# 	snIa = numpy.random.uniform(size=N_sn)
# 	snIa = numpy.less_equal(snIa, frac_Ia_0 + (frac_Ia_1-frac_Ia_0)*zs/zmax/1.5).astype(int)

# 	wsnIa = numpy.where(snIa)[0]
# 	wnonIa = numpy.where(numpy.logical_not(snIa))[0]

# 	# adu : the observed counts
# 	adu = 1/(cosmo.luminosity_distance(zs).value/cosmo.hubble_distance.value)**2
# 	adu_random = numpy.random.normal(size=N_sn)
# 	adu[wsnIa] = alpha_snIa*adu[wsnIa]*10**(adu_random[wsnIa]*sigma_snIa/2.5)
# 	adu[wnonIa] = alpha_nonIa*adu[wnonIa]*10**(adu_random[wnonIa]*sigma_nonIa/2.5)

# 	# host_choice : 1 if correct host is chosen to have highest probability, 0 otherwise
# 	host_choice = numpy.random.binomial(1,0.98,size=N_sn)
# 	neighbor_zs=numpy.random.uniform((1+zmin/1.5)**3,(1+zmax*1.5)**3,size=N_sn)**(1./3)-1
# 	host_zs_random = zs*host_choice + neighbor_zs*(1-host_choice)
# 	neighbor_zs_random = zs*(1-host_choice) + neighbor_zs*host_choice
# 	# wrong_host_zs : indeces if incorrect host galaxy association
# 	#wrong_host_zs = host_zs_random_ == 0
# 	#host_zs_random_[wrong_host_zs] = 



# 	if ADU0 is None:
# 		found = numpy.zeros(len(adu),dtype='bool')
# 		found[:]=True
# 		ADU0=0
# 	else:
# 		found  = adu>=ADU0


# 	N_s_obs = frac_obs*found.sum()

# 	if ia_only:
# 		s_obs_random = numpy.where(numpy.logical_and(found, snIa ==1))[0]
# 	else:
# 		s_obs_random = numpy.where(found)[0]

# 	# s_obs_random : order in which supernovae get a spectrum
# 	numpy.random.shuffle(s_obs_random)
# 	s_obs=numpy.sort(s_obs_random[:N_s_obs])
# 	s_mis= numpy.sort(s_obs_random[N_s_obs:])
# 	N_sn = len(s_obs)+len(s_mis)
# 	N_s_obs=len(s_obs)
# 	s_obs = s_obs.tolist()
# 	s_mis = s_mis.tolist()

# 	data = {'N_sn':N_sn,
# 			'N_obs':N_s_obs,

# 			'N_adu_max':1,

# 			'zmin':zmin,
# 			'zmax':zmax,

# 			'adu_obs': adu[s_obs,None],
# 			'adu_mis': adu[s_mis,None],

# 			'trans_ainv_obs': 1+zs[s_obs],
# 			'snIa_obs': snIa[s_obs],

# 			'host_zs_obs': host_zs_random[s_obs], #zs[s_obs],
# 			'host_zs_mis': host_zs_random[s_mis],
# 			'host2_zs_mis': neighbor_zs_random[s_mis],

# 			'ADU0': ADU0,
# #			'n_int': 25,
# 			'N_SNIa': snIa[s_obs].sum()
# 			}

# 	init=[]
# 	for i in xrange(Ninit):
# 		# host_zs_mis_init : initial conditions of host redshift _mis for MCMC
# 		#host_choice = numpy.random.binomial(1,0.98,size=N_sn)
# 		#host_zs_mis_init = (host_zs_random_*host_choice + neighbor_zs_random_*(1-host_choice))[s_mis]
# 		# neighbor_zs_random__ = host_zs*_random_(1-host_choice) + neighbor_zs_random_*host_choice
# 		# host_zs_mis_init=[]
# 		# for j in xrange(len(s_mis)):
# 		# 	z_cand = numpy.array([host_zs_random[s_mis[j]][0][0],host_zs_random[s_mis[j]][1][0],host_zs_random[s_mis[j]][2][0]])
# 		# 	ind = numpy.random.multinomial(1,[host_zs_random[s_mis[j]][0][1],host_zs_random[s_mis[j]][1][1],host_zs_random[s_mis[j]][2][1]])
# 		# 	ind = numpy.where(ind ==1)[0][0]
# 		# 	host_zs_mis_init.append(z_cand[ind])
# 		# host_zs_mis_init = numpy.array(host_zs_mis_init)

# #		host_zs_mis_init = host_zs_random[s_mis]*numpy.random.binomial(1,0.98,size=len(s_mis))
# #		wrong_host_zs_init = host_zs_mis_init == 0
# #		host_zs_mis_init[wrong_host_zs_init] = numpy.random.uniform(zmin/1.5,zmax*1.5,size=wrong_host_zs_init.sum())
# 		rate = numpy.random.lognormal(numpy.log(frac_Ia_0),0.1)
# 		rate1 = numpy.random.lognormal(numpy.log(frac_Ia_1),0.1)
# 		init.append ( {
# 			'Omega_M':numpy.random.normal(omega_M,0.1),
# 			'Omega_L':1-omega_M,
# 			'w': numpy.random.normal(-0.9,0.1),
# 			# 'ainv_true_obs': 1+zs[s_obs],
# 		 #  	'ainv_true_mis': (i % 2)*(1+host_zs_random[s_mis]) + (1-(i%2))*(1+zs[s_mis]), #host_zs_mis_init,
# 		  	'alpha_Ia': numpy.random.normal(alpha_snIa,0.1),
# 		 	'alpha_nonIa': numpy.random.normal(alpha_nonIa,0.1),
# 		  	'sigma_Ia': numpy.random.lognormal(numpy.log(sigma_snIa),0.05),
# 		  	'sigma_nonIa':numpy.random.lognormal(numpy.log(sigma_nonIa),0.5),
# 		  	'snIa_rate_0': [rate,1-rate],
# 		  	'snIa_rate_1': [rate1,1-rate1]
# 			} )

# 	info = {
# 			'zs' : zs,							# the true redshifts
# 			'snIa' : snIa,						# the true types
# 			'host_zs_random' : host_zs_random, 	# redshift of associated host galaxy IF there were no spectra
# 			'host_choice' : host_choice,	# indeces if incorrect host galaxy association
# 			's_obs' : s_obs,
# 			's_mis' : s_mis
# 		}

# 	return data, init, info

def main():
	Nchains=4
	N_sn=2000
	ADU0=1.
	ia_only=False

	data= Data(N_sn, 1, ia_only=ia_only)
	data.found(ADU0)

	sm = pystan.StanModel(file='des.stan')

	fracspec = numpy.arange(.2,1.01,.4)
	for ns in fracspec:
		data.spectrum(ns)

#		data, init, info = genData(N_sn,ns,Nchains,1, ia_only=ia_only, ADU0=ADU0)

		fit = sm.sampling(data=data.dict(), iter=1000, thin=1, n_jobs=-1, chains=Nchains, init=data.init(Nchains))
		#print fit
		#samples = fit.extract(['Omega_M','ainv_true_mis','w'])

		logposterior = fit.get_logposterior()

		app=''
		if ia_only:
			app+='.ia_only.'
		with open('../results/temp/model'+app+str(ns)+'.pkl', 'wb') as f:
			pickle.dump([fit.extract(), logposterior], f)

if __name__ == "__main__":
    main()
