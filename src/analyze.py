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

def individual(N_s, ia_only = False, ADU0=None,N_sn=None,dir='temp'):
	app=''
	if N_sn is not None:
		app='.'+str(N_sn)+'.'
	if ia_only:
		app+='ia_only.'
	if ADU0 == 0.:
			app+='noMalm.'
	[extract,logposterior] = pickle.load(file('../results/'+dir+'/model'+app+str(N_s)+'.pkl','rb'))
#	print numpy.abs((extract['w'] < -1).sum()*1. / len(extract['w'])-0.5)*2
	n_samples = len(extract[extract.keys()[0]])

	if ia_only:
		samples=numpy.zeros((n_samples,4))
		samples[:,0] = extract['Omega_M']
		samples[:,1] = extract['w']
		samples[:,2] = extract['alpha_Ia']
		samples[:,3] = extract['sigma_Ia']
		figure = triangle.corner(samples, labels=[r"$\Omega_M$", r"$w$", r"$\alpha_{Ia}$", r"$\sigma_{Ia}$"],truths=[0.28,-1, 2.,0.1],
			quantiles=[0.16, 0.5, 0.84], plot_datapoints=False,plot_density=True)
		plt.savefig('../results/'+dir+'/contour'+app+str(N_s)+'.pdf')
		plt.clf()
	else:
		samples=numpy.zeros((n_samples,8))
		samples[:,0] = extract['Omega_M']
		samples[:,1] = extract['w']
		samples[:,2] = extract['alpha_Ia']
		samples[:,3] = extract['sigma_Ia']
		samples[:,4] = extract['alpha_nonIa']
		samples[:,5] = extract['sigma_nonIa']
		samples[:,6] = extract['snIa_rate_0'][:]
		samples[:,7] = extract['snIa_rate_1'][:]

		kwargs = {'levels':1.0 - numpy.exp(-0.5 * numpy.arange(1, 3.1, 1) ** 2)}
		figure = triangle.corner(samples, labels=[r"$\Omega_M$", r"$w$", r"$\alpha_{Ia}$", r"$\sigma_{Ia}$", r"$\alpha_{non-Ia}$", r"$\sigma_{non-Ia}$", r"SN Ia Rate 0", r"SN Ia Rate1"],
			truths=[0.28,-1, 2.,0.1,2*10**(-2./2.5),1,0.95,0.2],
			quantiles=[0.16, 0.5, 0.84], plot_datapoints=False,plot_density=True, **kwargs)
		plt.savefig('../results/'+dir+'/contour'+app+str(N_s)+'.analyze.pdf')
		plt.clf()


	plt.plot(logposterior[0][500:],label=0)
	plt.plot(logposterior[1][500:],label=1)
	plt.plot(logposterior[2][500:],label=2)
	plt.plot(logposterior[3][500:],label=3)
	plt.ylabel('log-posterior')
	plt.xlabel('link')
	plt.legend(loc=3,prop={'size':9})
	plt.tight_layout()
	plt.savefig('../results/'+dir+'/posterior'+app+str(N_s)+'.pdf')
	plt.clf()


	plt.plot(extract['snIa_rate_0'])
	plt.ylabel(r'snIa_rate_0')
	plt.xlabel('link')
	plt.tight_layout()
	plt.savefig('../results/'+dir+'/snIa_rate_0'+app+str(N_s)+'.pdf')
	plt.clf()

	plt.plot(extract['snIa_rate_1'])
	plt.ylabel(r'snIa_rate_1')
	plt.xlabel('link')
	plt.tight_layout()
	plt.savefig('../results/'+dir+'/snIa_rate_1'+app+str(N_s)+'.pdf')
	plt.clf()

	plt.plot(extract['w'])
	plt.ylabel(r'$w$')
	plt.xlabel('link')
	plt.tight_layout()
	plt.savefig('../results/'+dir+'/w'+app+str(N_s)+'.pdf')
	plt.clf()
	
	# plt.hist(extract['w'],normed=True)
	# plt.xlabel(r'$w$')
	# plt.legend(loc=2)
	# plt.tight_layout()
	# plt.savefig('../results/'+dir+'/w.'+app+str(N_s)+'.pdf')

	# plt.clf()

	ans = numpy.zeros((2,len(levels)))
	deltas = numpy.zeros((2,len(levels)/2))
	threechain = numpy.append(extract['w'][:1000],extract['w'][1500:])
	wsort = numpy.sort(threechain)
	ans[0,:] = wsort[numpy.round(levels*len(wsort)).astype(int)]
	wsort = numpy.sort(extract['w'][1000:1500])
	ans[1,:] = wsort[numpy.round(levels*len(wsort)).astype(int)]
	for i in xrange(len(levels)/2):
		deltas[:,i] = ans[:,-1-i]-ans[:,i]

	# label=['013', '2']
	# for j in xrange(2):
	# 	for i in xrange(len(levels)/2-1,-1,-1):
	# 		print '{}& ${:4.2f}$ & $[{:5.3f}, {:5.3f}]$ & ${:5.3f}$ \\\\'.format(label[j],1-levels[i]*2, ans[j,i],ans[j,-1-i], deltas[j,i])

	n, b = numpy.histogram(samples[:,1], bins=20)

	ans = numpy.zeros(len(levels))
	deltas = numpy.zeros(len(levels)/2)
	wsort = numpy.sort(extract['w'])
	ans = wsort[numpy.round(levels*len(wsort)).astype(int)]
	median  = wsort[len(wsort)/2]
	mean = numpy.mean(ans)
	std = numpy.std(ans)
	for i in xrange(len(levels)/2):
		deltas[i] = ans[-1-i]-ans[i] 

	index = 0

	for i in xrange(len(levels)/2-1,-1,-1):
#		print '1.0& ${:4.2f}$ &${:4.2f}$ & $[{:5.3f}, {:5.3f}]$ & ${:5.3f}$ \\\\'.format(1-levels[i]*2, b[numpy.argmax(n)],  ans[i],ans[-1-i], deltas[i])
		if index ==0:
			print '$1.00$& ${:4.2f} \pm {:4.2f}$ &${:4.2f}$ &${:4.2f}$ & $[{:5.3f}, {:5.3f}]$ & ${:5.3f}$ \\\\'.format(mean, std/numpy.sqrt(len(ans)), median, 1-levels[i]*2,  ans[i],ans[-1-i], deltas[i])
		else:
			print '&  &  & ${:4.2f}$ & $[{:5.3f}, {:5.3f}]$ & ${:5.3f}$ \\\\'.format(1-levels[i]*2,  ans[i],ans[-1-i], deltas[i])

		index +=1


	# for j in xrange(len(nspec)):
	# 	for i in xrange(len(levels)/2-1,-1,-1):
	# 		print '{:3.0f}& ${:4.2f}$ & $[{:5.3f}, {:5.3f}]$ & ${:5.3f}$ \\\\'.format(nspec[j]/2.,1-levels[i]*2, ans[j,i],ans[j,-1-i], deltas[j,i])


	# plt.hist([threechain,extract['w'][1000:1500]],normed=True,label=['Chains 0,1,3','Chain 2'])
	# plt.legend(loc=2)
	# plt.xlabel(r'$w$')
	# # plt.hist(extract['w'][1000:1500],label='Chain 2')
	# plt.tight_layout()
	# plt.savefig('../results/temp/wsub.'+app+str(N_s)+'.pdf')
	# plt.clf()

	# missed_wrong = info['host_choice'][info['s_mis']] == 0
	# if missed_wrong.sum() > 0:
	# 	plt.plot(extract[2000,missed_wrong,0])

	# 	print info['zs'][info['s_mis']][missed_wrong]
	# 	print info['host_zs_random'][info['s_mis']][missed_wrong]
	# 	label=[]
	# 	for i in xrange(missed_wrong.sum()):
	# 		ans= r"$T="+info['snIa'][info['s_mis']][missed_wrong][i].astype(str)+"$, $z_{true}="+info['zs'][info['s_mis']][missed_wrong][i].astype('S4')+"$, "+"$z_{host}="+info['host_zs_random'][info['s_mis']][missed_wrong][i].astype('S4')+"$"
	# 		label.append(ans)
	# 	lineobjects=plt.plot(extract['ainv_true_mis'][:,missed_wrong] - 1)
	# 	plt.legend(lineobjects,label,loc=2,prop={'size':9})
	# 	plt.ylabel(r'$z$')
	# 	plt.xlabel('link')
	# 	plt.tight_layout()
	# 	plt.savefig('../results/temp/missed_wrong.'+app+str(N_s)+'.pdf')
	# 	plt.clf()
	# 	# plt.plot(extract['zs_true_mis'][:,missed_wrong] - info['host_zs_random'][info['s_mis']][missed_wrong])
	# 	# plt.show()
	# 	print info['snIa'][info['s_mis']][missed_wrong]



def group(nspec):

	ans = numpy.zeros((len(nspec), len(levels)))
	deltas = numpy.zeros((len(nspec),len(levels)/2))
	ind=0
	for n in nspec:
		[extract, logposterior] = pickle.load(file('../results/temp/model'+str(n)+'.pkl','rb'))
		wsort = numpy.sort(extract['w'])
		ans[ind,:] = wsort[numpy.round(levels*len(wsort)).astype(int)]
		for i in xrange(len(levels)/2):
			deltas[ind,i] = ans[ind,-1-i]-ans[ind,i]

		print numpy.abs((wsort < -1).sum()*1. / len(wsort)-0.5)*2
		ind +=1


	for j in xrange(len(nspec)):
		for i in xrange(len(levels)/2-1,-1,-1):
			print '{:3.0f}& ${:4.2f}$ & $[{:5.3f}, {:5.3f}]$ & ${:5.3f}$ \\\\'.format(nspec[j],1-levels[i]*2, ans[j,i],ans[j,-1-i], deltas[j,i])


def main():
	individual(0.0,ia_only=False,ADU0=0.75,N_sn=2000,dir='seed2_pop2')
#	individual(500,ia_only=False)
	wefew
	group([200,350,500])

if __name__ == "__main__":
    main()
