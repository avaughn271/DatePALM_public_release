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
	parser.add_argument('--sMax',type=float,default=1.1)
	parser.add_argument('--df',type=int,default=450)
	parser.add_argument('--noAlleleTraj', default=False, action='store_true', help='whether to compute the posterior allele frequency trajectory or not.')
	parser.add_argument('--integration_points', type=int, default = -1)
	parser.add_argument('--h', type=float, default = 0.5)

	return parser.parse_args()

def likelihood(theta, args):
	Xvals = args[0]
	Yvals = args[1]
	numdimensions = int(round((np.sqrt( len(theta) * 8 + 1)  - 1)/2.0))
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
			try:
				res = (minimize_scalar(likelihood_wrapper_scalar, bracket = [0.9,1.0,1.1],args=minargs, method = "Brent", tol = 1e-4))
				S = [res.x - 1.0] # adjusted wrapper to work on selection + 1, so that tolerance makes more sense.
				L = res.fun

				if times.shape[1] < 31: # with small number of samples, result function can be multimodal, so do multiple optims and take best.
					try:
						res = (minimize_scalar(likelihood_wrapper_scalar, bracket = [1.0001,1.00011,1.02],args=minargs, method = "Brent", tol = 1e-4))
						Stemp = [res.x - 1.0] # adjusted wrapper to work on selection + 1, so that tolerance makes more sense.
						Ltemp = res.fun
						if Ltemp < L:
							S = Stemp
							L = Ltemp
					except ValueError:
						pass
			except ValueError:
				try:
					print("Selection MLE not found in [-0.1,0.1], possibly due to noninformative data. Expanding search to [-1,1].")
					res = (minimize_scalar(likelihood_wrapper_scalar, bracket = [0.00000001,1.0,2.0],args=minargs, method = "Brent", tol = 1e-4))
					S = [res.x - 1.0] # adjusted wrapper to work on selection + 1, so that tolerance makes more sense.
					L = res.fun
				except ValueError:
					print("Selection MLE not found in [-1,1], possibly due to noninformative data. SNP likelihood fit will not be saved.")
					if os.path.exists(args.out+"_tempfile.txt"):
						os.remove(args.out+"_tempfile.txt")
					exit(0)

		else:
			res = minimize(likelihood_wrapper, S0, args=minargs, options=opts, method='Nelder-Mead')
			S = res.x
			L = res.fun

		toprint = '%.4f'%(-L+logL0)
		numericloglik = -L+logL0

		Weights = []
	else:
		M = times.shape[2]
		Weights = np.zeros(M)
		precompute = 0
		tranmatrix = _nstep_log_trans_prob(Ne[0],0.0,freqs,z_bins,z_logcdf,z_logsf,h) # this only handles the precompute case, just use initial values
		if len(np.unique(Ne)) == 1:
			precompute = 1
		for i in range(M):
			betaMatl0 = backward_algorithm(np.zeros(len(Ne)),times[:,:,i],derSampledTimes,ancSampledTimes,epochs,Ne,h,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,tranmatrix,noCoals=noCoals,precomputematrixboolean=precompute,currFreq=currFreq)
			Weights[i] = logsumexp(betaMatl0[-2,:])

		minargs = (timeBins,Ne,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,functionvals,Weights)

		if len(S0) == 1:
			try:
				res = (minimize_scalar(likelihood_wrapper_scalar, bracket = [0.9,1.0,1.1],args=minargs, method = "Brent", tol = 1e-4))
				S = [res.x - 1.0] # adjusted wrapper to work on selection + 1, so that tolerance makes more sense.
				L = res.fun
				if times.shape[1] < 31: # with small number of samples, result function can be multimodal, so do multiple optims and take best.
					try:
						res = (minimize_scalar(likelihood_wrapper_scalar,  bracket = [1.0001,1.00011,1.02] ,args=minargs, method = "Brent", tol = 1e-4))
						Stemp = [res.x - 1.0] # adjusted wrapper to work on selection + 1, so that tolerance makes more sense.
						Ltemp = res.fun
						if Ltemp < L:
							S = Stemp
							L = Ltemp
					except ValueError:
						pass
			except ValueError:
				try:
					print("Selection MLE not found in [-0.1,0.1], possibly due to noninformative data. Expanding search to [-1,1].")
					res = (minimize_scalar(likelihood_wrapper_scalar, bracket = [0.00000001,1.0,2.0],args=minargs, method = "Brent", tol = 1e-4))
					S = [res.x - 1.0] # adjusted wrapper to work on selection + 1, so that tolerance makes more sense.
					L = res.fun
				except ValueError:
					print("Selection MLE not found in [-1,1], possibly due to noninformative data. SNP likelihood fit will not be saved.")
					if os.path.exists(args.out+"_tempfile.txt"):
						os.remove(args.out+"_tempfile.txt")
					exit(0)
		else:
			res = minimize(likelihood_wrapper, S0, args=minargs, options=opts, method='Nelder-Mead')
			S = res.x
			L = res.fun
		numericloglik = -L
		toprint = '%.4f'%(-L)

	FirstLine = "logLR" + "\t" + "-log10(p-value)"
	epochnum = 1
	degreesoffreedom = len(S)

	toprint = toprint + "\t" + '%.2f'%(-(chi2.logsf(numericloglik + numericloglik, degreesoffreedom ) ) / np.log(10) )
	for s,t,u in zip(S,timeBins[:-1],timeBins[1:]):
		toprint = toprint + "\t" + '%d'%(t)
		if u <= args.tCutoff:
			toprint = toprint + "\t" + '%d'%(u)
		else:
			toprint = toprint + "\t" + '%d'%(args.tCutoff)
		toprint = toprint + "\t" + '%.5f'%(s)
		FirstLine = FirstLine + "\t" + "Epoch" + str(epochnum) + "_start" + "\t" +  "Epoch" + str(epochnum)   + "_end"  + "\t" +  "SelectionMLE" + str(epochnum)
		epochnum = epochnum + 1
	functionvals.close()

	if True:
		file1 = open(args.out+"_tempfile.txt", 'r')
		Lines = file1.readlines()
		Xvals = []
		Yvals = []
		indexxx = 0
		for i in Lines:
			Xvals.append([])
			DerivedSampleTimes = i.split(",")
			for j in range(len(DerivedSampleTimes) - 1 ):
				(Xvals[indexxx]).append(float(DerivedSampleTimes[j]) / 2.0 ) #HERE IS WHERE WE CONVERT BACK TO THE 1, 1+S,1+2S FORMULATION!!!!!
			Yvals.append(float(DerivedSampleTimes[len(DerivedSampleTimes) - 1 ]))
			indexxx = indexxx + 1
		Yvals = np.exp(np.subtract(Yvals, max(Yvals)))
		muu = Xvals[(list(Yvals)).index(max(Yvals))]

		if len(Xvals[0]) == 1:
			S0 =[1.0]
			fullresult =  minimize(likelihood, S0, args=[Xvals, Yvals, muu], method='Nelder-Mead', options={"maxfev":1000, "fatol":1e-20, "xatol":1e-20})
			res = fullresult.x
			print(fullresult.fun)
			standard_dev = res[0]
			toprint = ('%.16f'%(muu[0]) + "\n" +  '%.16f'%(standard_dev**2) + "\n")
			if fullresult.fun > 0.05:
				print("Poor fit of normal distribution to data. SNP likelihood fit will not be saved. Try increasing the df parameter.")
				if os.path.exists(args.out+"_tempfile.txt"):
					os.remove(args.out+"_tempfile.txt")
				exit(0)

			else:
				outputfile = open(args.out+".txt", "w+")
				outputfile.write(toprint)
				outputfile.close()
			#####result is mu=muu, var=standard_dev**2
		else:
			S0 =[0.0] * (round((len(Xvals[0])*len(Xvals[0])+len(Xvals[0]))/2)  )
			for innn in range(len(Xvals[0]) ):
				S0[innn] = 1.0
			res = minimize(likelihood, S0, args=[Xvals, Yvals, muu], method='Nelder-Mead', options={"maxfev":1000, "fatol":1e-40, "xatol":1e-40})
			for ifi in range(1,9):
				S0 =[0.0] * (round((len(Xvals[0])*len(Xvals[0])+len(Xvals[0]))/2)  )
				for innn in range(len(Xvals[0])  ):
					S0[innn] = 10**(-ifi)
				res1 = minimize(likelihood, S0, args=[Xvals, Yvals, muu], method='Nelder-Mead', options={"maxfev":1000, "fatol":1e-40, "xatol":1e-40})

				if res1.fun < res.fun:
					res = res1
			print("least-squares residual:", res.fun)
			res = res.x
			print(res)
			for iggi in res:
				if iggi > 0.2 or iggi < -0.2:
					print("Poor fit of normal distribution to data. Unreliable results follow.")

			#print("mu2: ", muu)
			#print("sd2: ", res)
			standard_dev = res
			numdimensions=  len(muu)

			covarmat = np.zeros((numdimensions, numdimensions))
			elementindex = 0

			for difference in range(numdimensions):
				for col in range(numdimensions - difference):
					covarmat[col + difference, col] = res[elementindex]
					elementindex = elementindex + 1
			for row in range(numdimensions):
				for col in range(numdimensions):
					if covarmat[row, col] == 0.0:
						covarmat[row, col] = covarmat[col, row]
			#result is mu=muu, cov=covarmat
			toprint = ""
			for i in muu:
				toprint = toprint + '%.16f'%(i) + "\n"
			for i in range(len(covarmat)):
				for j in range(len(covarmat[i])):
					toprint = toprint + '%.16f'%(covarmat[i][j]) + " "
				toprint = toprint[:-1]
				toprint = toprint + "\n"
			outputfile = open(args.out+".txt", "w+")
			outputfile.write(toprint)
			outputfile.close()
	if os.path.exists(args.out+"_tempfile.txt"):
		os.remove(args.out+"_tempfile.txt")