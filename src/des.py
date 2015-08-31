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
		self.zmax=1.4

		self.sigma_snIa=0.1
		self.sigma_nonIa=1

		self.alpha_snIa=2.
		self.alpha_nonIa=self.alpha_snIa*10**(-2./2.5)
		self.frac_Ia_0=.95
		self.frac_Ia_1=.2
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
		self.host_choice = numpy.random.binomial(1,0.98,size=self.N_sn)

		# zs : the true redshifts
#		volume_zero = self.cosmo.comoving_volume(self.zmin/1.1).value
#		volume_norm = self.cosmo.comoving_volume(self.zmax*1.1).value- volume_zero
#		self.neighbor_zs = numpy.sort(numpy.random.uniform(size=self.N_sn))
#		for i in xrange(len(self.zs)):
#			self.neighbor_zs[i]=scipy.optimize.newton(lambda z: (self.cosmo.comoving_volume(z).value - volume_zero)/volume_norm-self.zs[i], 0.5)		
		self.neighbor_zs=numpy.random.uniform((1+self.zmin/1.1)**3,(1+self.zmax*1.1)**3,size=self.N_sn)**(1./3)-1
		self.host_zs_random = self.zs*self.host_choice + self.neighbor_zs*(1-self.host_choice)
		self.neighbor_zs_random = self.zs*(1-self.host_choice) + self.neighbor_zs*self.host_choice

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
			'N_SNIa': self.snIa[self.s_obs].sum(),
			'N_adu_max':1,

			'zmin':self.zmin,
			'zmax':self.zmax,

			'adu_obs': self.adu[self.s_obs],
			'adu_mis': self.adu[self.s_mis],

			'trans_ainv_obs': 1+self.zs[self.s_obs],
			'snIa_obs': self.snIa[self.s_obs],

			'host_zs_obs': self.host_zs_random[self.s_obs], #zs[s_obs],
			'host_zs_mis': self.host_zs_random[self.s_mis],
			'host2_zs_mis': self.neighbor_zs_random[self.s_mis],

			'ADU0': self.ADU0
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
			  	'sigma_nonIa':numpy.random.lognormal(numpy.log(self.sigma_nonIa),0.1),
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
		     plt.scatter(self.zs[numpy.logical_and(self.found, self.snIa ==1)], -2.5*numpy.log10(self.adu[numpy.logical_and(self.found, self.snIa ==1)]),label='SN Ia',color='b',facecolors='none')
		     plt.scatter(self.zs[numpy.logical_and(self.found, self.snIa ==0)], -2.5*numpy.log10(self.adu[numpy.logical_and(self.found, self.snIa ==0)]),label='non-Ia',color='r',facecolors='none')
		     plt.scatter(self.zs[self.snIa ==1], -2.5*numpy.log10(self.adu[self.snIa ==1]),alpha=0.1,color='b')
		     plt.scatter(self.zs[self.snIa ==0], -2.5*numpy.log10(self.adu[self.snIa ==0]),alpha=0.1,color='r')
		     plt.scatter(self.zs[self.s_obs_random[:self.N_s_obs]], -2.5*numpy.log10(self.adu[self.s_obs_random[:self.N_s_obs]]),marker='.',color='k',label='Spectrum',s=20)
		     w  = self.s_obs_random[self.N_s_obs:]
		     w2 = numpy.logical_and(numpy.logical_not(self.host_choice[w]), self.found[w])
		     plt.scatter(self.host_zs_random[w[w2]], -2.5*numpy.log10(self.adu[w[w2]]),marker='x',color='g',label='Wrong Host',s=40)
		     plt.ylim((-6,2))
		     plt.xlim((0,1.5))
		     ax=plt.gca()
#		     ax.invert_yaxis()
		     plt.legend(loc=4)
		     pdf.savefig(fig)


def dataPlot():
	N_sn=2000
	ADU0=.5
	ia_only=False

	data= Data(N_sn, 1, ia_only=ia_only)
	data.found(ADU0)
	data.spectrum(0.2)
	data.plot()

def main():
	Nchains=4
	N_sn=2000
	ADU0=0.2 #0.2
	ia_only=False

	data= Data(N_sn, 1, ia_only=ia_only)
	data.found(ADU0)


	sm = pystan.StanModel(file='des.stan')

	fracspec = numpy.arange(.2,1.01,.4)
	for ns in fracspec:
		data.spectrum(ns)

		fit = sm.sampling(data=data.dict(), iter=2000, thin=1, n_jobs=-1, chains=Nchains, init=data.init(Nchains))

		logposterior = fit.get_logposterior()

		app=''
		if ia_only:
			app+='.ia_only.'
		with open('../results/temp/model'+app+str(ns)+'.pkl', 'wb') as f:
			pickle.dump([fit.extract(), logposterior], f)

if __name__ == "__main__":
#	dataPlot()
    main()
