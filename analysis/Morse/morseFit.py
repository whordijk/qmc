import numpy as np
from scipy import optimize,exp
import matplotlib.pyplot as plt

def errorQuadratic(p,s,E,weight):
    k,se,c = p
    fit = 0.5*k*(s-se)**2 + c
    error = fit-E
    error = error*weight
    return error

def errorMorse(p,s,E,weight):
    De,a,c,se = p
    fit = De*(1-exp(-a*(s-se)))**2 + c
    error = fit-E
    error = error*weight
    return error

def fitData(s,E,dE):
    weight = exp(-100**2*(E-min(E))**2)

    pQuad = [1,1.4, -1]
    pQuad,success = optimize.leastsq(errorQuadratic,pQuad,args=(s,E,weight))

    pMorse = [10,10,1.4,1]
    pMorse, success = optimize.leastsq(errorMorse,pMorse,args=(s,E,weight))
    return pQuad,pMorse

def makePlot(s,E,dE,p1,p2):
    k,se,c = p1
    fit1 = 0.5*k*(s-se)**2 + c
    De,a,c,se = p2
    fit2 = De*(1-exp(-a*(s-se)))**2 + c
    plt.plot(s,E,'o',label="Data")
    plt.plot(s,fit1,label="Quadratic")
    plt.plot(s,fit2,label="Morse")
    plt.ylim(min(E),max(E))
    plt.legend()
    plt.show()

#Filename for input data to minimize
fName = 'data.dat'

#Loading data
s,E,dE = np.loadtxt(fName,unpack=True)

#Fit data
pQuad,pMorse=fitData(s,E,dE)

#Print summary
print "Parameters for quadratic fit:"
print "     1/2*k*(s - se)**2"
print "     k = {0}".format(pQuad[0])    
print "     se = {0}".format(pQuad[1])
print "     c = {0}".format(pQuad[2])    

print ""

print "Parameters for Morse fit:"
print "     De*(1-exp(-a(s-se)))^2 + c"
print "     De = {0}".format(pMorse[0])
print "     a = {0}".format(pMorse[1])
print "     se = {0}".format(pMorse[2])
print "     c = {0}".format(pMorse[3])

#Plot results
makePlot(s,E,dE,pQuad,pMorse)
