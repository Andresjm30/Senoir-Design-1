import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
from sympy import symbols, Eq, solve, diff

A = np.array([[1,-1,0],[-2,-4,1],[0,1,1]])
eigvals, eigvecs = la.eig(A)
print("Eigenvalues: " + str(eigvals))
eigvals = eigvals.real
print("Real Eigenvalues: " + str(eigvals))
print("Eigenvectors: \n" + str(eigvecs))

maxval = abs(eigvals[0])
col = 0
for i in range(0, len(eigvals)):
     if(abs(eigvals[i]) > maxval):    
       maxval = eigvals[i];
       col = i;
maxvec = eigvecs[:,col].reshape(3,1)
print("Dominant Eiganvalue: " + str(maxval))
print("Dominant Eiganvector: \n" + str(maxvec))