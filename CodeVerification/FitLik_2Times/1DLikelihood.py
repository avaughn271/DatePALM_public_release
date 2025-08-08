import numpy as np
from hmm_utils import backward_algorithm, _nstep_log_trans_prob
from scipy.special import logsumexp
from scipy.optimize import minimize, minimize_scalar
import argparse
import os
from scipy.stats import chi2, norm, beta
from scipy.stats import multivariate_normal
from inference import load_data, likelihood_wrapper, likelihood_wrapper_scalar

def parse_args():
	"""Define the Arguments"""
	parser = argparse.ArgumentParser()
	parser.add_argument('--times',type=str, default=None)
	parser.add_argument('--popFreq',type=float,default=None)

	parser.add_argument('--ancientSamps',type=str,default=None)
	parser.add_argument('--ancientHaps',type=str,default=None)
	parser.add_argument('--out',type=str,default=None)

	parser.add_argument('--N',type=float,default=None)
	parser.add_argument('--coal',type=str,default=None,help='path to Relate .coal file. Negates --N option.')

	parser.add_argument('--tCutoff',type=float,default=None)
	parser.add_argument('--timeBins',type=str,default=None, nargs='+')
	parser.add_argument('--sMax',type=float,default=0.1)
	parser.add_argument('--df',type=int,default=450)
	parser.add_argument('--noAlleleTraj', default=False, action='store_true', help='whether to compute the posterior allele frequency trajectory or not.')
	parser.add_argument('--integration_points', type=int, default = -1)
	parser.add_argument('--h', type=float, default = 0.5)

	return parser.parse_args()

def likelihood(theta, args):
	Xvals = args[0]
	Yvals = args[1]
	numdimensions = round((np.sqrt( len(theta) * 8 + 1)  - 1)/2.0)
	if numdimensions == 1: #1d optimization
		standarddev = theta[0]
		scalarr  = 1/norm.pdf( args[2][0],  args[2][0], standarddev)
		if standarddev <=0:
			return(10000000000.0)
		FUNC = 0
		for i in range(len(Xvals)):
			FUNC = FUNC + (Yvals[i] - scalarr * norm.pdf(Xvals[i], loc = args[2][0], scale = standarddev))**2
		return(FUNC)
	else:		
		cholfactor = np.zeros((numdimensions, numdimensions))
		elementindex = 0

		for difference in range(numdimensions):
			for col in range(numdimensions - difference):
				cholfactor[col + difference, col] = theta[elementindex]
				elementindex = elementindex + 1
		for row in range(numdimensions):
			for col in range(numdimensions):
				if cholfactor[row, col] == 0.0:
					cholfactor[row, col] = cholfactor[col, row]
		if not np.all(np.linalg.eigvals(cholfactor) > 0):
			return(10000000000.0)
		try:
			scalarr  = 1/multivariate_normal.pdf(args[2], mean = args[2], cov = cholfactor)
		except:
			return(10000000000.0)

		FUNC = 0
		for i in range(len(Xvals)):
			try:
				FUNC = FUNC + (Yvals[i] - scalarr * multivariate_normal.pdf(Xvals[i], mean = args[2], cov = cholfactor))**2
			except:
				return(10000000000.0)
				
		return(FUNC)

if __name__ == "__main__":
	args = parse_args()
	if args.times == None and args.ancientSamps == None and args.ancientHaps == None:
		print('You need to supply coalescence times (--times) and/or ancient samples (--ancientSamps)')
	
	# load data and set up model
	#should delete the previous verion
	if os.path.exists(args.out+"_tempfile.txt"):
		os.remove(args.out+"_tempfile.txt")
	functionvals = open(args.out+"_tempfile.txt", "a")
	sMax = args.sMax
	timeBins,times,epochs,Ne,freqs,ancientGLs,ancientHapGLs,noCoals,currFreq,logfreqs,log1minusfreqs,derSampledTimes,ancSampledTimes,h = load_data(args)
	# read in global Phi(z) lookups
	linearspacing = np.linspace(0.0, 1.0, 2000)
	linearspacing[0] = 1e-10
	linearspacing[len(linearspacing) - 1] = 1 - 1e-10

	z_bins = norm.ppf(linearspacing)
	z_logcdf = norm.cdf(z_bins)
	z_logsf = z_logcdf # can delete this later.

	Ne *= 1/2
	noCoals = int(noCoals)

	# optimize over selection parameters
	T = len(timeBins)
	S0 = 0.0 * np.ones(T-1)
	opts = {}

	if T == 2:
		Simplex = np.reshape(np.array([-0.05,0.05]),(2,1))
	elif T > 2:
		Simplex = np.zeros((T,T-1))
		for i in range(Simplex.shape[1]):
			Simplex[i,:] = 0.0
			Simplex[i,i] = 0.01
		Simplex[-1,:] = 0.0
	else:
		raise ValueError
	opts['disp']=False
	opts['maxfev'] = (T - 1) * 20
	opts['initial_simplex']=Simplex

	ImpSamp = False
	if times.shape[2] > 1:
		print('\t(Importance sampling with M = %d samples)'%(times.shape[2]))
		print()
		ImpSamp = True
	if not ImpSamp: # to account for whether we return the likelihood or the log likelihood

		logL0 = likelihood_wrapper(S0,timeBins,Ne,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,functionvals)

		minargs = (timeBins,Ne,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,functionvals)

		if len(S0) == 1:
			for iiiiii in np.arange(-0.01,0.03, 0.0002):
				print(iiiiii/2.0, -likelihood_wrapper([iiiiii],timeBins,Ne,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,functionvals,[]))
		else:
			for iiiiii in np.arange(0.0,0.08, 0.002):
				for jj in np.arange(-0.02,0.02, 0.0002):
					print(iiiiii/2.0, jj/2.0, -likelihood_wrapper([iiiiii, jj],timeBins,Ne,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,functionvals,[]))

	functionvals.close()

	if os.path.exists(args.out+"_tempfile.txt"):
		os.remove(args.out+"_tempfile.txt")