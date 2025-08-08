import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal

mu = [0.0203558504581451 , -0.0040019392967224] # hard-coded this for now!!!!!!!!!!!!!!
covar = [[0.0001885077209431, -0.0001442197279983],[-0.0001442197279983, 0.0001150900507139]]

A  = pd.read_csv("truevalues2d.txt", sep = " ", header = None)
A = A.to_numpy()
X = A[:,0]
Y = A[:,1]
Z = A[:,2]
Z = Z  - np.max(Z)
Z = np.exp(Z)
nx = len(np.unique(X))
ny = len(np.unique(Y))

X  = np.reshape(X, (nx,ny))
Y = np.reshape(Y, (nx,ny))

QUANTILES = [0.1*i for i in range(1,11)]

Grayscalecolors = []
LightestGray = 0.85
DarkestGray = 0.3
for i in range(len(QUANTILES)):
	Grayscalecolors.append(str(LightestGray - i/(len(QUANTILES)-1) * (len(QUANTILES))*(LightestGray - DarkestGray)/len(QUANTILES)))

ZFitted = np.zeros((nx,ny))
for i in range(nx):
	for j in range(ny):
		ZFitted[i,j] = multivariate_normal.logpdf([X[i,j], Y[i,j]], mu, covar)
ZFitted = ZFitted  - np.max(ZFitted)
ZFitted = np.exp(ZFitted)
plt.contour(X, Y, ZFitted, QUANTILES, cmap = 'plasma') #cividis is bad, magma is bad, inferno is bad
plt.contourf( X ,  Y, np.reshape(Z, (nx,ny)), QUANTILES ,colors=Grayscalecolors)
plt.suptitle('Likelihood Function of Selection Gradient')

plt.gca().set_aspect(1)
plt.savefig('Output.pdf',format='pdf', bbox_inches='tight')
plt.clf()