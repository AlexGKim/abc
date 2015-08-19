#!/usr/bin/env python

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pystan
import numpy.random
from astropy.cosmology import FlatLambdaCDM

import pickle
import des

import triangle

levels = numpy.array([0.025, 0.05, 0.16, 1-0.16, 1-0.05, 1-.025])

def individual(N_s, ia_only = False):
	app = ''
	if ia_only:
		app+='.ia_only.'
	[extract, info, logposterior] = pickle.load(file('../results/model'+app+str(N_s)+'.pkl','rb'))

	n_samples = len(extract[extract.keys()[0]])
	samples=numpy.zeros((n_samples,2))
	samples[:,0] = extract['Omega_M']
	samples[:,1] = extract['w']
	figure = triangle.corner(samples, labels=[r"$\Omega_M$", r"$w$"],truths=[0.28,-1],
		quantiles=[0.16, 0.5, 0.84], plot_datapoints=False,plot_density=True)
	plt.savefig('../results/contour.'+app+str(N_s)+'.pdf')
	plt.clf()

	matplotlib.rcParams.update(matplotlib.rcParamsDefault)

	plt.plot(logposterior[0][500:],label=0)
	plt.plot(logposterior[1][500:],label=1)
	plt.plot(logposterior[2][500:],label=2)
	plt.plot(logposterior[3][500:],label=3)
	plt.ylabel('log-posterior')
	plt.xlabel('link')
	plt.legend(loc=3,prop={'size':9})
	plt.tight_layout()
	plt.savefig('../results/posterior.'+app+str(N_s)+'.pdf')
	plt.clf()

	# plt.plot(extract['w'])
	# plt.show()

	plt.hist(extract['w'],normed=True)
	plt.xlabel(r'$w$')
	plt.legend(loc=2)
	plt.tight_layout()
	plt.savefig('../results/w.'+app+str(N_s)+'.pdf')

	plt.clf()

	ans = numpy.zeros((2,len(levels)))
	deltas = numpy.zeros((2,len(levels)/2))
	threechain = numpy.append(extract['w'][:1000],extract['w'][1500:])
	wsort = numpy.sort(threechain)
	ans[0,:] = wsort[numpy.round(levels*len(wsort)).astype(int)]
	wsort = numpy.sort(extract['w'][1000:1500])
	ans[1,:] = wsort[numpy.round(levels*len(wsort)).astype(int)]
	for i in xrange(len(levels)/2):
		deltas[:,i] = ans[:,-1-i]-ans[:,i]

	label=['013', '2']
	for j in xrange(2):
		for i in xrange(len(levels)/2-1,-1,-1):
			print '{}& ${:4.2f}$ & $[{:5.3f}, {:5.3f}]$ & ${:5.3f}$ \\\\'.format(label[j],1-levels[i]*2, ans[j,i],ans[j,-1-i], deltas[j,i])

	ans = numpy.zeros(len(levels))
	deltas = numpy.zeros(len(levels)/2)
	wsort = numpy.sort(extract['w'])
	ans = wsort[numpy.round(levels*len(wsort)).astype(int)]
	for i in xrange(len(levels)/2):
		deltas[i] = ans[-1-i]-ans[i] 
	for i in xrange(len(levels)/2-1,-1,-1):
		print '${:4.2f}$ & $[{:5.3f}, {:5.3f}]$ & ${:5.3f}$ \\\\'.format(1-levels[i]*2, ans[i],ans[-1-i], deltas[i])



	# for j in xrange(len(nspec)):
	# 	for i in xrange(len(levels)/2-1,-1,-1):
	# 		print '{:3.0f}& ${:4.2f}$ & $[{:5.3f}, {:5.3f}]$ & ${:5.3f}$ \\\\'.format(nspec[j]/2.,1-levels[i]*2, ans[j,i],ans[j,-1-i], deltas[j,i])


	plt.hist([threechain,extract['w'][1000:1500]],normed=True,label=['Chains 0,1,3','Chain 2'])
	plt.legend(loc=2)
	plt.xlabel(r'$w$')
	# plt.hist(extract['w'][1000:1500],label='Chain 2')
	plt.tight_layout()
	plt.savefig('../results/wsub.'+app+str(N_s)+'.pdf')
	plt.clf()

	missed_wrong = info['host_choice'][info['s_mis']] == 0
	if missed_wrong.sum() > 0:
		print info['zs'][info['s_mis']][missed_wrong]
		print info['host_zs_random'][info['s_mis']][missed_wrong]
		label=[]
		for i in xrange(missed_wrong.sum()):
			ans= r"$T="+info['snIa'][info['s_mis']][missed_wrong][i].astype(str)+"$, $z_{true}="+info['zs'][info['s_mis']][missed_wrong][i].astype('S4')+"$, "+"$z_{host}="+info['host_zs_random'][info['s_mis']][missed_wrong][i].astype('S4')+"$"
			label.append(ans)
		lineobjects=plt.plot(extract['ainv_true_mis'][:,missed_wrong] - 1)
		plt.legend(lineobjects,label,loc=2,prop={'size':9})
		plt.ylabel(r'$z$')
		plt.xlabel('link')
		plt.tight_layout()
		plt.savefig('../results/missed_wrong.'+app+str(N_s)+'.pdf')
		plt.clf()
		# plt.plot(extract['zs_true_mis'][:,missed_wrong] - info['host_zs_random'][info['s_mis']][missed_wrong])
		# plt.show()
		print info['snIa'][info['s_mis']][missed_wrong]



def group(nspec):

	ans = numpy.zeros((len(nspec), len(levels)))
	deltas = numpy.zeros((len(nspec),len(levels)/2))
	ind=0
	for n in nspec:
		[extract, info, logposterior] = pickle.load(file('../results/model'+str(n)+'.pkl','rb'))
		wsort = numpy.sort(extract['w'])
		ans[ind,:] = wsort[numpy.round(levels*len(wsort)).astype(int)]
		for i in xrange(len(levels)/2):
			deltas[ind,i] = ans[ind,-1-i]-ans[ind,i] 
		ind +=1


	for j in xrange(len(nspec)):
		for i in xrange(len(levels)/2-1,-1,-1):
			print '{:3.0f}& ${:4.2f}$ & $[{:5.3f}, {:5.3f}]$ & ${:5.3f}$ \\\\'.format(nspec[j]/2.,1-levels[i]*2, ans[j,i],ans[j,-1-i], deltas[j,i])


def main():
	individual(100,ia_only=False)
	wefew
	group([100,150,200])

if __name__ == "__main__":
    main()