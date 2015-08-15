#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pystan
import numpy.random
from astropy.cosmology import FlatLambdaCDM

import pickle
import des

import triangle

def individual(N_s):
	[extract, info, logposterior] = pickle.load(file('../results/model'+str(N_s)+'.pkl','rb'))

	n_samples = len(extract[extract.keys()[0]])
	samples=numpy.zeros((n_samples,2))
	samples[:,0] = extract['Omega_M']
	samples[:,1] = extract['w']
	figure = triangle.corner(samples, labels=[r"$\Omega_M$", r"$w$"],truths=[0.28,-1],
		quantiles=[0.16, 0.5, 0.84], plot_datapoints=False,plot_density=True)
	plt.savefig('../results/contour.'+str(N_s)+'.pdf')
	plt.clf()

	matplotlib.rcParams.update(matplotlib.rcParamsDefault)

	plt.plot(logposterior[0][500:],label=0)
	plt.plot(logposterior[1][500:],label=1)
	plt.plot(logposterior[2][500:],label=2)
	plt.plot(logposterior[3][500:],label=3)
	plt.ylabel('log-posterior')
	plt.xlabel('link')
	plt.legend(loc=3,prop={'size':9})
	plt.savefig('../results/posterior.'+str(N_s)+'.pdf')
	plt.clf()

	# plt.plot(extract['w'])
	# plt.show()

	plt.hist(extract['w'],normed=True)
	plt.xlabel(r'$w$')
	plt.legend(loc=2)
	plt.tight_layout()
	plt.savefig('../results/w.'+str(N_s)+'.pdf')

	plt.clf()
	plt.hist([numpy.append(extract['w'][:1000],extract['w'][1500:]),extract['w'][1000:1500]],normed=True,label=['Chains 0,1,3','Chain 2'])
	plt.legend(loc=2)
	plt.xlabel(r'$w$')
	# plt.hist(extract['w'][1000:1500],label='Chain 2')
	plt.tight_layout()
	plt.savefig('../results/wsub.'+str(N_s)+'.pdf')
	plt.clf()
	missed_wrong = info['host_choice'][info['s_mis']] == 0
	print info['zs'][info['s_mis']][missed_wrong]
	print info['host_zs_random'][info['s_mis']][missed_wrong]
	label=[]
	for i in xrange(missed_wrong.sum()):
		ans= r"$T="+info['snIa'][info['s_mis']][missed_wrong][i].astype(str)+"$, $z_{true}="+info['zs'][info['s_mis']][missed_wrong][i].astype('S4')+"$, "+"$z_{host}="+info['host_zs_random'][info['s_mis']][missed_wrong][i].astype('S4')+"$"
		label.append(ans)
	if missed_wrong.sum() !=0:
		lineobjects=plt.plot(extract['ainv_true_mis'][:,missed_wrong] - 1)
	plt.legend(lineobjects,label,loc=2,prop={'size':9})
	plt.tight_layout()
	plt.savefig('../results/missed_wrong.'+str(N_s)+'.pdf')
	plt.clf()
	# plt.plot(extract['zs_true_mis'][:,missed_wrong] - info['host_zs_random'][info['s_mis']][missed_wrong])
	# plt.show()
	print info['snIa'][info['s_mis']][missed_wrong]



def group(nspec):

	for n in nspec:
		[extract, info] = pickle.load(file('model'+int(n)+'.pkl','rb'))
		plt.hist(extract['w'])
		plt.show()


def main(n_s):
	individual(n_s)

if __name__ == "__main__":
    main(100)