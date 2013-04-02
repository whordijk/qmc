import numpy as np
from scipy import optimize,exp
import matplotlib.pyplot as plt

def errorQuadratic(p,s,E,weight):
    a,b,c=p
    fit = a*s**2+b*s+c
    error = fit-E
    error = error*weight
    return error

def errorMorse(p,s,E,weight):
    De,a,se = p
    fit = De*(1-exp(-a*(s-se)))**2
    error = fit-E
    error = error*weight
    return error

def fitData(s,E,dE):
    weight = exp(-100**2*(E-min(E))**2)

    pQuad = [1,-2.8,-1.1]
    pQuad,success = optimize.leastsq(errorQuadratic,pQuad,args=(s,E,weight))

    pMorse = [10,10,1.4]
    pMorse, success = optimize.leastsq(errorMorse,pMorse,args=(s,E,weight))
    return pQuad,pMorse

def makePlot(s,E,dE,p1,p2):
    a,b,c=p1
    fit1 = a*s**2+b*s+c
    De,a,se = p2
    fit2 = De*(1-exp(-a*(s-se)))**2
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
print "     a*x^2 + b*x + c"
print "     a = {0}".format(pQuad[0])    
print "     b = {0}".format(pQuad[1])    
print "     c = {0}".format(pQuad[2])    

print ""

print "Parameters for Morse fit:"
print "     De*(1-exp(-a(s-se)))^2"
print "     De = {0}".format(pMorse[0])
print "     a = {0}".format(pMorse[1])
print "     se = {0}".format(pMorse[2])

#Plot results
makePlot(s,E,dE,pQuad,pMorse)
