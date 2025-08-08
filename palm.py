import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import norm, multivariate_normal
from sys import exit

def _args(parser):
	# mandatory inputs:
	required = parser.add_argument_group('required arguments')
	required.add_argument('--metadata',type=str,help='A dataframe holding attributes for each SNP')
	required.add_argument('--snpDir',type=str,help='A directory holding the info of each snp')
	parser.add_argument('--makePlot', default=False, action='store_true')
	parser.add_argument('--skipMissingSNPs', default=False, action='store_true')

	parser.add_argument('--out',type=str,default=None) #output prefix.
	parser.add_argument('--B',type=int,default=1000) #number of bootstrap replicates
	parser.add_argument('--maxp',type=float,default=1) #only snps with p vals below this threshold should be included

	return parser

def _parse_loci_stats(args): # have checked for case of a single trait
	
	coeffs = []
	betas = []
	mults = []
	LDblocks = []
	standarderrors = []
	df = pd.read_csv(args.metadata,sep='\t',index_col=(0,1),header=0)
	#metadata, the first column are the linkage blocks
	traitNames = []
	for i in range(len(df.columns)):
		if '@' in df.columns[i]:
			TraitName = df.columns[i].split('@')[1]
			if TraitName not in traitNames:
				traitNames.append(TraitName)
	traitNames.sort()
		
	if len(traitNames) == 0: # only analyze one trait.
		betaColumns = ['beta']
		pColumns = ['pval']
		traitNames = ['']
		seColumns = ['se']

	else:
		betaColumns = []
		pColumns = []
		seColumns = []
		for trait in traitNames:
			betaColumns.append('beta@' + trait)
			pColumns.append('pval@' + trait)
			seColumns.append('se@' + trait)

	dfFiltered = df
	idxs = dfFiltered.index.values

	RelevantSNPIndices = []
	for (k,idx) in enumerate(idxs):
		dfRow = df.loc[idx]

		loc_pvals = np.array(dfRow[pColumns],dtype=float)

		if np.any(np.isnan(list(np.array(dfRow[betaColumns])))): #remove missing data rows!!!!!!!!!!!!!
			print("Skipping row %d due to missing data"%(k))
			continue
		if np.all(loc_pvals >= args.maxp): #if all p-vals for the traits are above threshold, continue.
			continue
		RelevantSNPIndices.append(k)

	df = df.iloc[RelevantSNPIndices]
	idxs = df.index.values

	for (k,idx) in enumerate(idxs):
		dfRow = df.loc[idx]
		variant = idx[1]
		cols = variant.split(':')
		if len(cols) != 4:
			print("Improperly formatted variant name: %s"%(variant))
			
		if dfRow.derived_allele != cols[-1]:
			flipper = -1.0
		else:
			flipper = 1.0

		loc_pvals = np.array(dfRow[pColumns],dtype=float)
		LISTOFCOEFS = []
		if not args.skipMissingSNPs:
			with open(args.snpDir + '/%s.txt'%(cols[0] + "_" + cols[1])) as f:
				LISTOFCOEFS =  [line for line in f]
		else:
			try:
				with open(args.snpDir + '/%s.txt'%(cols[0] + "_" + cols[1])) as f:
					LISTOFCOEFS =  [line for line in f]
			except:
				print("Skipping SNP ", cols[0] + "_" + cols[1] + " as it is missing.")
				continue
		coeffstr = []
		for i in range(len(LISTOFCOEFS)):
			coeffstr.extend(LISTOFCOEFS[i].split())
		coeff = []
		for i in range(len(coeffstr)):
			coeff.append(float(coeffstr[i].rstrip('\n')))
		coeffs.append(coeff)
		print(cols, coeff)
		betas.append(list(flipper * np.array(dfRow[betaColumns])))
		standarderrors.append(list(np.array(dfRow[seColumns])))
		mults.append(1/len(df.loc[idx[0]].index))
		LDblocks.append(idx[0])
	if len(betas) == 0:
		print("No SNPs left to analyze!")
		exit(0)

	betas = np.array(betas)
	standarderrors = np.array(standarderrors)
	mults = np.array(mults)
	coeffs = np.array(coeffs)
	LDblocks = np.array(LDblocks)
	# better way to read in the coefficients
	#for single test, coeffs has #snps rows and 3 columns, each column is a coefficient. Same for joint test.
	#for single test, betas is array of betas if derived=minor, else -betas. dimension number of traits by number of snps.
	#for single test, mults is all 1.0. For joint, it is 1/K where K is the number of snps that are in that block.
	#for single test, traitNames will be ['']
	return coeffs,betas,mults,LDblocks,traitNames,standarderrors

def calculatelogLR(omega, stats):
	
	coeffs,betas,mults,LDblocks,betaerrors = stats
	numberofepochs = int(np.floor(np.sqrt(len(coeffs[0]))))
	L = betas.shape[0]
	SUM = 0

	if numberofepochs == 1 and len(betas[0]) == 1:
		for l in range(L):
			SUM += mults[l] * norm.logpdf(betas[l,0] * omega[0], coeffs[l,0],  (coeffs[l,1] + betaerrors[l,0]*betaerrors[l,0]*omega[0] * omega[0]) ** 0.5)
			
	elif numberofepochs == 1:
		for l in range(L):
			SUM += mults[l]  * norm.logpdf(betas[l,0] * omega[0] + betas[l,1] * omega[1],
								   coeffs[l,0], 
									 (coeffs[l,1]  + betaerrors[l,0]*betaerrors[l,0]*omega[0] * omega[0]+ betaerrors[l,1]*betaerrors[l,1]*omega[1] * omega[1]) ** 0.5)
			
	elif len(betas[0]) == 1: #seems to work but did not rigorously check this.
		for l in range(L):
			xval = np.multiply(betas[l,0], omega)
			fittedmu = coeffs[l,0:numberofepochs]
			fittedcovarmatrix = np.reshape(coeffs[l,numberofepochs:], (numberofepochs, numberofepochs))
			covarresult = fittedcovarmatrix + np.outer(omega, omega) * betaerrors[l,0] * betaerrors[l,0]
			try:
				SUM += mults[l] * multivariate_normal.logpdf(xval, fittedmu, covarresult)
			except:
				for temporaryvar in [1e+00,2e+00,5e+00,1e+01,2e+01,5e+01,1e+02,2e+02,5e+02,1e+03,2e+03,5e+03,1e+04,2e+04,5e+04,1e+05,2e+05,5e+05,1e+06,2e+06,5e+06,1e+07,2e+07,5e+07, 1e+08,2e+08,5e+08,1e+09,2e+09,5e+09,1e+10,2e+10,5e+10]:
					try:
						SUM += mults[l] * multivariate_normal.logpdf(xval, fittedmu, covarresult + temporaryvar * 1e-9 * np.eye(covarresult.shape[0]))
						break
					except:
						pass
	else:
		for l in range(L):
			fittedmu = coeffs[l,0:numberofepochs]
			fittedcovarmatrix = np.reshape(coeffs[l,numberofepochs:], (numberofepochs, numberofepochs))
			xval = [0] * numberofepochs
			for i in range(len(omega)):
				xval[i % numberofepochs] = xval[i % numberofepochs] + betas[l, int(np.floor(i/numberofepochs))] * omega[i]
			covarresult = fittedcovarmatrix +  \
				np.outer(omega[0:numberofepochs], omega[0:numberofepochs]) * betaerrors[l,0] * betaerrors[l,0] +  \
				np.outer(omega[numberofepochs:], omega[numberofepochs:]) * betaerrors[l,1] * betaerrors[l,1]
			try:
				SUM += mults[l] * multivariate_normal.logpdf(xval, fittedmu, covarresult)
			except:
				for temporaryvar in [1e+00,2e+00,5e+00,1e+01,2e+01,5e+01,1e+02,2e+02,5e+02,1e+03,2e+03,5e+03,1e+04,2e+04,5e+04,1e+05,2e+05,5e+05,1e+06,2e+06,5e+06,1e+07,2e+07,5e+07, 1e+08,2e+08,5e+08,1e+09,2e+09,5e+09,1e+10,2e+10,5e+10]:
					try:
						SUM += mults[l] * multivariate_normal.logpdf(xval, fittedmu, covarresult + temporaryvar * 1e-9 * np.eye(covarresult.shape[0]))
						break
					except:
						pass

	return SUM

def negativelogLR(omega, coeffs,betas,mults,LDblocks,betaerrors):
	return(-calculatelogLR(omega, (coeffs,betas,mults,LDblocks,betaerrors) ))

def _opt_omega(stats, numberofepochs):
	coeffs,betas,mults,LDblocks,betaerrors = stats
	J = betas.shape[1]
	statstemp = coeffs,betas,mults,LDblocks,betaerrors
	initialomega = [0.0] * J * numberofepochs
	return(minimize(negativelogLR, initialomega, statstemp, method='BFGS', options = {'maxiter': 10, "gtol" : 1e-200, "xrtol" : 0.000005}).x)

def _inference(statistics,args):
	numberofepochs = int(np.floor(np.sqrt(len(statistics[0][0]))))

	omega = _opt_omega(statistics, numberofepochs)

	print('Analyzing %d loci...'%(len(statistics[2])))
	B = args.B
	
	omegaJK = np.zeros((numberofepochs,B))
	for b in range(B):
		coeffs,betas,mults,LDblocks,betaerrors = statistics
		UNIQUEBLOCKS = np.unique(LDblocks)
		I = np.random.choice(UNIQUEBLOCKS,len(UNIQUEBLOCKS),replace=True)
		INDICES = []
		for i in range(len(I)):
			INDICES.extend(np.where(LDblocks == I[i])[0])
		statsDK =  coeffs[INDICES,:],betas[INDICES,:],mults[INDICES],LDblocks[INDICES],betaerrors[INDICES,:]
		omegaJK_b = _opt_omega(statsDK, numberofepochs)
		omegaJK[:,b] = omegaJK_b
	ses = np.std(omegaJK,axis=1)

	EpochDiffs = []
	EpochDiffStandardErrors = []
	for i in range(numberofepochs - 1):
		EpochDiffs.append(omega[i] - omega[i + 1])
		EpochDiffStandardErrors.append(np.std(omegaJK[i, :] - omegaJK[i + 1, :]))
	EpochDiffStandardErrors.append(-1)
	EpochDiffs.append(-1)

	return omega,ses,EpochDiffs,EpochDiffStandardErrors

def _T_inference(statistics,args):
	numberofepochs = int(np.floor(np.sqrt(len(statistics[0][0]))))
	#Calculate the joint fit omega.
	coeffs,betas,mults,LDblocks,betaerrors = statistics
	L = len(statistics[2])
	print('Analyzing %d loci...'%(L))
	omega = _opt_omega(statistics,numberofepochs)
	
	J = betas.shape[1]
	margOmega = np.zeros(J * numberofepochs)
	for j in range(J):
		mbetas = np.reshape(betas[:,j], (L,1))
		mbetaerrors = np.reshape(betaerrors[:,j], (L,1))
		mstats = coeffs,mbetas,mults,LDblocks,mbetaerrors
		margOmega[(j*numberofepochs):(j * numberofepochs + numberofepochs)] = _opt_omega(mstats,numberofepochs)

	#Calculate R
	R_statistic = omega - margOmega

	B = args.B
	###calculate the joint standard errors:
	omega_joint = np.zeros((J * numberofepochs,B))
	omega_marginal = np.zeros((J * numberofepochs,B))
	
	for b in range(B):
		UNIQUEBLOCKS = np.unique(LDblocks)
		I = np.random.choice(UNIQUEBLOCKS,len(UNIQUEBLOCKS),replace=True)
		INDICES = []
		for i in range(len(I)):
			INDICES.extend(np.where(LDblocks == I[i])[0])
		
		#This is the joint matrix
		statsjoint =  coeffs[INDICES,:],betas[INDICES,:],mults[INDICES],LDblocks[INDICES],betaerrors[INDICES,:]
		omega_joint[:,b] = _opt_omega(statsjoint,numberofepochs)

		#this is the marginal matrix
		for j in range(J):
			mbetas = np.reshape(betas[:,j], (L,1))
			mbetaerrors = np.reshape(betaerrors[:,j], (L,1))
			stats_marginal =  coeffs[INDICES,:],mbetas[INDICES,:],mults[INDICES],LDblocks[INDICES],mbetaerrors[INDICES,:]
			omega_marginal[(j*numberofepochs):(j * numberofepochs + numberofepochs), b] = _opt_omega(stats_marginal,numberofepochs)
	JointSE = np.std(omega_joint,axis=1)
	margSE = np.std(omega_marginal,axis=1)

	Rstandarderrors = np.std(omega_joint - omega_marginal, axis = 1)

	Delta_jointSE = np.zeros((J * numberofepochs,B))
	for traitindex in range(J):
		for epochindex in range(numberofepochs):
			if (numberofepochs * traitindex + epochindex + 1) < (J * numberofepochs):
				Delta_jointSE[numberofepochs * traitindex + epochindex,:] = omega_joint[(numberofepochs * traitindex + epochindex+1) , :] - omega_joint[(numberofepochs * traitindex + epochindex) , :]
	Delta_joint = omega[0:(J * numberofepochs - 1)] - omega[1:]

	np.append(Delta_joint , -1)

	Delta_RSE = np.zeros((J * numberofepochs,B))
	for traitindex in range(J):
		for epochindex in range(numberofepochs):
			if (numberofepochs * traitindex + epochindex + 1) < (J * numberofepochs):
				Delta_RSE[numberofepochs * traitindex + epochindex,:] = (omega_joint - omega_marginal)[(numberofepochs * traitindex + epochindex+1) , :] - (omega_joint - omega_marginal)[(numberofepochs * traitindex + epochindex) , :]
	Delta_R = R_statistic[0:(J * numberofepochs - 1)] - R_statistic[1:]
	np.append(Delta_R , -1)

	return omega,JointSE,margOmega,margSE,R_statistic,Rstandarderrors,Delta_joint,np.std(Delta_jointSE, axis = 1),Delta_R,np.std(Delta_RSE, axis = 1)

def prettyformat(numberr):
	stringnumber = "{:10.4f}".format(numberr)
	while stringnumber[0] == " ":
		stringnumber = stringnumber[1:]
	return(stringnumber)

def determineoptimalplotdimensions(margOmega, statistics):
	Lowerbound = 0.1 #CHANGED FROM 0.1 TO 1.0
	Upperbound = 200.0 #CHANGED FROM 200.0 TO 20.0
	while Upperbound - Lowerbound > 0.25:
		xsize = (Lowerbound + Upperbound)/2.0
		xlower = margOmega[0] - xsize
		ylower = margOmega[1] - xsize
		xupper =  margOmega[0] + xsize
		yupper = margOmega[1] + xsize

		delta =  (xupper - xlower) / 40.0
		x = np.arange(xlower, xupper, delta)
		y = np.arange(ylower, yupper, delta)
		X, Y = np.meshgrid(x, y)
		Z = np.zeros((len(x),len(y)))
		for i in range(len(x)):
			for j in range(len(y)):
				Z[i,j] = calculatelogLR([X[i,j], Y[i,j]], statistics)

		##Plot marginal
		Z_marginal = np.zeros((len(x),len(y)))
		coeffs,betas,mults,LDblocks,betaerrors = statistics
		L = len(statistics[2])

		statstrait1 = coeffs,np.reshape(betas[:,0], (L,1)),mults,LDblocks,np.reshape(betaerrors[:,0], (L,1))
		statstrait2 = coeffs,np.reshape(betas[:,1], (L,1)),mults,LDblocks,np.reshape(betaerrors[:,1], (L,1))
		for i in range(len(x)):
			for j in range(len(y)):
				Z_marginal[i,j] = calculatelogLR( [X[i,j]], statstrait1) + calculatelogLR( [Y[i,j]],statstrait2)
		Z = np.exp(Z - np.max(Z))
		Z_marginal = np.exp(Z_marginal - np.max(Z_marginal))
		max1 = np.max([np.max(Z_marginal[0,:]),  np.max(Z_marginal[:,0]), np.max(Z_marginal[(len(x)-1),:]), np.max(Z_marginal[:,(len(x)-1)])])
		max2 = np.max([np.max(Z[0,:]),  np.max(Z[:,0]), np.max(Z[(len(x)-1),:]), np.max(Z[:,(len(x)-1)])])
		if (max1 > 0.005) | (max2 > 0.005):
			Lowerbound = xsize
		else:
			Upperbound = xsize
	return((Lowerbound + Upperbound)/2)

def determineoptimalplotdimensionsmultiple(margOmega, omega, statistics, plotnumber, numberofepochs):
	Lowerbound = 0.1 #CHANGED FROM 0.1 TO 1.0
	Upperbound = 200.0 #CHANGED FROM 200.0 TO 20.0
	while Upperbound - Lowerbound > 0.25:
		xsize = (Lowerbound + Upperbound)/2.0
		xlower = margOmega[plotnumber] - xsize
		ylower = margOmega[plotnumber + numberofepochs] - xsize
		xupper =  margOmega[plotnumber] + xsize
		yupper = margOmega[plotnumber + numberofepochs] + xsize

		delta =  (xupper - xlower) / 40.0
		x = np.arange(xlower, xupper, delta)
		y = np.arange(ylower, yupper, delta)
		X, Y = np.meshgrid(x, y)
		Z = np.zeros((len(x),len(y)))
		for i in range(len(x)):
			for j in range(len(y)):
				fullomega = np.copy(omega)
				fullomega[plotnumber] = X[i,j]
				fullomega[plotnumber + numberofepochs] = Y[i,j]
				Z[i,j] = calculatelogLR(fullomega,statistics) # maybe unnecessary if statement

		##Plot marginal
		Z_marginal = np.zeros((len(x),len(y)))
		coeffs,betas,mults,LDblocks,betaerrors = statistics
		L = len(statistics[2])

		statstrait1 = coeffs,np.reshape(betas[:,0], (L,1)),mults,LDblocks,np.reshape(betaerrors[:,0], (L,1))
		statstrait2 = coeffs,np.reshape(betas[:,1], (L,1)),mults,LDblocks,np.reshape(betaerrors[:,1], (L,1))
		for i in range(len(x)):
			for j in range(len(y)):
				fullomega1 = np.copy(margOmega[0:numberofepochs])
				fullomega1[plotnumber] = X[i,j]
				fullomega2 = np.copy(margOmega[numberofepochs:])
				fullomega2[plotnumber] = Y[i,j]
				Z_marginal[i,j] = calculatelogLR(fullomega1,statstrait1) + calculatelogLR(fullomega2,statstrait2) # maybe unnecessary statement.
		Z = np.exp(Z - np.max(Z))
		Z_marginal = np.exp(Z_marginal - np.max(Z_marginal))
		max1 = np.max([np.max(Z_marginal[0,:]),  np.max(Z_marginal[:,0]), np.max(Z_marginal[(len(x)-1),:]), np.max(Z_marginal[:,(len(x)-1)])])
		max2 = np.max([np.max(Z[0,:]),  np.max(Z[:,0]), np.max(Z[(len(x)-1),:]), np.max(Z[:,(len(x)-1)])])
		if (max1 > 0.005) | (max2 > 0.005):
			Lowerbound = xsize
		else:
			Upperbound = xsize
	return((Lowerbound + Upperbound)/2)

def findquantile(MATRIX,quantile):
	lowerbound = 0
	upperbound = 1
	MatrixSum = np.sum(MATRIX)
	while upperbound - lowerbound > 0.001:
		threshold = (upperbound + lowerbound)/2
		if np.sum(MATRIX[MATRIX > threshold]) > quantile * MatrixSum:
			lowerbound = threshold
		else:
			upperbound = threshold
	return((upperbound + lowerbound)/2)

def _main(args):
	coeffs,betas,mults,LDblocks,traitNames,betaerrors = _parse_loci_stats(args)
	# L is total number of loci being analyzed
	# L_byTrait is a list of number of loci retained for each trait??

	statistics = coeffs,betas,mults,LDblocks,betaerrors
	IsJointTest = (betas.shape[1] > 1)

	if IsJointTest:
		# run test jointly
		omega,JointSE,margOmega,margSE,R_statistic,Rstandarderrors,deltastat,deltastatse,deltarstat,deltarstatse = _T_inference(statistics,args)
	else:
		omega,ses,delta_stat,delta_se = _inference(statistics,args)

	if args.out != None:
		if IsJointTest:
			if len(omega) == 2:
				toprint = "Trait" + "\t" + \
						"jgrad_mle" + "\t" + "jgrad_se" + "\t" + "jgrad_z" + "\t" \
						+ "mgrad_mle" + "\t" + "mgrad_se" + "\t" + "mgrad_z" + "\t" \
						+ "R_mle" + "\t"+ "R_se" + "\t" + "R_z" + "\n"
				j = 0
				for trait in traitNames:
					toprint = toprint + trait + "\t" + \
					prettyformat(omega[j]) + "\t" + prettyformat(JointSE[j]) + "\t" + prettyformat(omega[j]/JointSE[j]) + "\t"  + \
					prettyformat(margOmega[j]) + "\t" + prettyformat(margSE[j]) + "\t" + prettyformat(margOmega[j]/margSE[j]) + "\t" + \
					prettyformat(R_statistic[j]) + "\t" + prettyformat(Rstandarderrors[j]) + "\t" + prettyformat(R_statistic[j]/Rstandarderrors[j])  + "\n"
					j = j + 1
			else:
				toprint = "Trait_epoch" + "\t" + \
						"jgrad_mle" + "\t" + "jgrad_se" + "\t" + "jgrad_z" + "\t" \
						+ "mgrad_mle" + "\t" + "mgrad_se" + "\t" + "mgrad_z" + "\t" \
						+ "R_mle" + "\t"+ "R_se" + "\t" + "R_z" +  "\t"  \
						+ "delta_mle" + "\t"+ "delta_se" + "\t" + "delta_z" +  "\t"  \
							+ "deltaR_mle" + "\t"+ "deltaR_se" + "\t" + "deltaR_z" +  "\n"
				j = 0
				for trait in traitNames:
					for timeindex in range(int(len(omega)/len(traitNames))):
						if timeindex != (int(len(omega)/len(traitNames)) - 1):
							toprint = toprint + trait + "_" + str(timeindex + 1) + "\t" + \
							prettyformat(omega[j]) + "\t" + prettyformat(JointSE[j]) + "\t" + prettyformat(omega[j]/JointSE[j]) + "\t"  + \
							prettyformat(margOmega[j]) + "\t" + prettyformat(margSE[j]) + "\t" + prettyformat(margOmega[j]/margSE[j]) + "\t" + \
							prettyformat(R_statistic[j]) + "\t" + prettyformat(Rstandarderrors[j]) + "\t" + prettyformat(R_statistic[j]/Rstandarderrors[j])  + "\t" + \
							prettyformat(deltastat[j]) + "\t" + prettyformat(deltastatse[j]) + "\t" + prettyformat(deltastat[j]/deltastatse[j])  + "\t" + \
							prettyformat(deltarstat[j]) + "\t" + prettyformat(deltarstatse[j]) + "\t" + prettyformat(deltarstat[j]/deltarstatse[j])  + "\n"
						else:
							toprint = toprint + trait + "_" + str(timeindex + 1) + "\t" + \
							prettyformat(omega[j]) + "\t" + prettyformat(JointSE[j]) + "\t" + prettyformat(omega[j]/JointSE[j]) + "\t"  + \
							prettyformat(margOmega[j]) + "\t" + prettyformat(margSE[j]) + "\t" + prettyformat(margOmega[j]/margSE[j]) + "\t" + \
							prettyformat(R_statistic[j]) + "\t" + prettyformat(Rstandarderrors[j]) + "\t" + prettyformat(R_statistic[j]/Rstandarderrors[j])  + "\t" + \
							"-1"+ "\t" + "-1" + "\t" + "-1"  + "\t" + \
							"-1" + "\t" + "-1" + "\t" + "-1"  + "\n"
						j = j + 1
			f = open(args.out + ".txt", "w+")
			f.writelines(toprint)
			f.close()
		else:
			if len(omega) == 1:
				toprint = "sgrad_mle" + "\t" + "sgrad_se" + "\t" + "sgrad_z" + "\n"
				toprint = toprint + prettyformat(omega[0]) +  "\t" + prettyformat(ses[0]) + "\t" +  prettyformat(omega[0]/ses[0]) + "\n"
			else:
				toprint = "sgrad_mle" + "\t" + "sgrad_se" + "\t" + "sgrad_z" + "\t" + "delta_mle" + "\t" + "delta_se" + "\t" + "delta_z" + "\n"
				for numberofomega in range(len(omega)):
					if numberofomega != (len(omega) - 1):
						toprint = toprint + prettyformat(omega[numberofomega]) +  "\t" + prettyformat(ses[numberofomega]) + "\t" +  prettyformat(omega[numberofomega]/ses[numberofomega]) + "\t" + prettyformat(delta_stat[numberofomega]) +  "\t" + prettyformat(delta_se[numberofomega]) + "\t" +  prettyformat(delta_stat[numberofomega]/delta_se[numberofomega]) + "\n"
					else:
						toprint = toprint + prettyformat(omega[numberofomega]) +  "\t" + prettyformat(ses[numberofomega]) + "\t" +  prettyformat(omega[numberofomega]/ses[numberofomega]) + "\t" + "-1" +  "\t" + "-1" + "\t" +  "-1" + "\n"
			f = open(args.out + ".txt", "w+")
			f.writelines(toprint)
			f.close()
	if IsJointTest and args.makePlot:
		print("A")
		numberofepochs = int(np.floor(np.sqrt(len(statistics[0][0]))))
		##should still rigorusly check that this does what we want.
		for plotnumber in range(numberofepochs):
			if numberofepochs == 1:
				xsize = determineoptimalplotdimensions(margOmega, statistics) #LOCATION1
			else:
				xsize = determineoptimalplotdimensionsmultiple(margOmega, omega, statistics, plotnumber, numberofepochs) #LOCATION1

			QUANTILES = [1-((i+1)/(8.0)) for i in range(7)]   #7
			print("B")
			Grayscalecolors = []
			LightestGray = 0.85
			DarkestGray = 0.3
			for i in range(len(QUANTILES)):
				Grayscalecolors.append(str(LightestGray - i/(len(QUANTILES)-1) * (len(QUANTILES))*(LightestGray - DarkestGray)/len(QUANTILES)))

			JOINTContours = []
			MARGINALContours = []
			#we evaluate on a square grid.
			xlower = margOmega[plotnumber] - xsize #LOCATION2
			ylower = margOmega[plotnumber + numberofepochs] - xsize
			xupper =  margOmega[plotnumber] + xsize
			yupper =  margOmega[plotnumber + numberofepochs] + xsize

			delta =  (xupper - xlower) / 200.0
			x = np.arange(xlower, xupper, delta)
			y = np.arange(ylower, yupper, delta)
			X, Y = np.meshgrid(x, y)
			Z = np.zeros((len(x),len(y)))
			for i in range(len(x)):
				for j in range(len(y)):
					if numberofepochs == 1:
						Z[i,j] = calculatelogLR( [X[i,j], Y[i,j]],statistics )  #LOCATION3
					else:
						fullomega = np.copy(omega)
						fullomega[plotnumber] = X[i,j]
						fullomega[plotnumber + numberofepochs] = Y[i,j]
						Z[i,j] = calculatelogLR(fullomega,statistics ) # maybe unnecessary if statement

			Z =  np.exp(Z - np.max(Z))
			print("C")
			
			##Plot marginal
			Z_marginal = np.zeros((len(x),len(y)))
			coeffs,betas,mults,LDblocks,betaerrors = statistics
			L = len(statistics[2])

			statstrait1 = coeffs,np.reshape(betas[:,0], (L,1)),mults,LDblocks,np.reshape(betaerrors[:,0], (L,1))
			statstrait2 = coeffs,np.reshape(betas[:,1], (L,1)),mults,LDblocks,np.reshape(betaerrors[:,1], (L,1))
			for i in range(len(x)):
				for j in range(len(y)):
					if numberofepochs == 1:
						Z_marginal[i,j] = calculatelogLR( [X[i,j]],statstrait1 ) + calculatelogLR( [Y[i,j]],statstrait2 )  #LOCATION4
					else:
						fullomega1 = np.copy(margOmega[0:numberofepochs])
						fullomega1[plotnumber] = X[i,j]
						fullomega2 = np.copy(margOmega[numberofepochs:])
						fullomega2[plotnumber] = Y[i,j]
						Z_marginal[i,j] = calculatelogLR(fullomega1,statstrait1 ) + calculatelogLR(fullomega2,statstrait2 ) # maybe unnecessary statement.
			Z_marginal = np.exp(Z_marginal - np.max(Z_marginal))
			#We have now calculated Z and Z_marginal, so we can find the appropriate quantiles.
			print("D")
			for i in range(len(QUANTILES)):
				JOINTContours.append(findquantile(Z, QUANTILES[i]))
			for i in range(len(QUANTILES)):
				MARGINALContours.append(findquantile(Z_marginal, QUANTILES[i]))
			MARGINALContours.append(0.99999)
			plt.contour(X, Y, Z,JOINTContours, cmap = 'plasma') #cividis is bad, magma is bad, inferno is bad
			plt.contourf(X, Y, Z_marginal,MARGINALContours ,colors=Grayscalecolors)
			plt.suptitle('Likelihood Function of Selection Gradient')
			plt.xlabel(traitNames[0])
			plt.ylabel(traitNames[1])
			plt.plot(omega[plotnumber], omega[plotnumber + numberofepochs],  'o', markersize=4, markerfacecolor='#ffff66',
				markeredgewidth = 0.5, markeredgecolor = "k")

			plt.plot(margOmega[plotnumber], margOmega[plotnumber + numberofepochs],  "Hk", markersize=4)
			#We then readjust the bounds of the plot to fully contain the contours.
			bounds = 0
			for i in range(X.shape[0]):
				for j in range(X.shape[1]):
					if Z[i,j] >= JOINTContours[0] and abs(X[i,j] - margOmega[plotnumber]) > bounds:
						bounds =  abs(X[i,j] - margOmega[plotnumber])
					if Z_marginal[i,j] >= MARGINALContours[0] and abs(X[i,j] - margOmega[plotnumber]) > bounds:
						bounds =  abs(X[i,j] - margOmega[plotnumber])
					if Z[i,j] >= JOINTContours[0] and abs(Y[i,j] - margOmega[plotnumber + numberofepochs]) > bounds:
						bounds =  abs(Y[i,j] - margOmega[plotnumber + numberofepochs])
					if Z_marginal[i,j] >= MARGINALContours[0] and abs(Y[i,j] - margOmega[plotnumber + numberofepochs]) > bounds:
						bounds =  abs(Y[i,j] - margOmega[plotnumber + numberofepochs])
			plt.xlim(margOmega[plotnumber] - bounds*1.1, margOmega[plotnumber] + bounds*1.1)
			plt.ylim(margOmega[plotnumber + numberofepochs] - bounds*1.1, margOmega[plotnumber + numberofepochs] + bounds*1.1)
			plt.gca().set_aspect(1)
			plt.savefig(args.out + '_' + str(plotnumber + 1) + '.pdf',format='pdf', bbox_inches='tight')
			plt.clf()

super_parser = argparse.ArgumentParser()
parser = _args(super_parser)
args = parser.parse_args()
_main(args)