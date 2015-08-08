#!/usr/bin/env python

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pystan
import numpy.random
from astropy.cosmology import FlatLambdaCDM

import pickle
import des

def individual(N_sn, N_s):
	[extract, info, logposterior] = pickle.load(file('model'+str(N_s)+'.pkl','rb'))

	plt.plot(logposterior[0][200:],label=0)
	plt.plot(logposterior[1][200:],label=1)
	plt.plot(logposterior[2][200:],label=2)
	plt.plot(logposterior[3][200:],label=3)
	plt.legend()
	plt.show()
	# plt.plot(extract['w'])
	# plt.show()
	# plt.hist(extract['w'])
	# plt.show()

	plt.plot(extract['w'])
	plt.show()

	missed_wrong = info['host_choice'][info['s_mis']] == 0
	print missed_wrong.sum()
	if missed_wrong.sum() !=0:
		plt.plot(extract['ainv_true_mis'][:,missed_wrong] - info['zs'][info['s_mis']][missed_wrong])
		plt.show()
	# plt.plot(extract['zs_true_mis'][:,missed_wrong] - info['host_zs_random'][info['s_mis']][missed_wrong])
	# plt.show()
	print info['snIa'][info['s_mis']][missed_wrong]



def group(N_sn, nspec):

	for n in nspec:
		[extract, info] = pickle.load(file('model'+int(n)+'.pkl','rb'))
		plt.hist(extract['w'])
		plt.show()


def main(n_s):
	individual(200,n_s)

if __name__ == "__main__":
    main(200)