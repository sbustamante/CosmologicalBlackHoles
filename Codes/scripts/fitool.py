#!/usr/bin/env python
import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import leastsq
random=np.random.random
mag=np.linalg.norm

def residuals(pars,x,y,func):
    err=y-func(x,pars)
    return err

def fit(xdata,ydata,func,pguess,verbose=True,logfile=""):
    """
    Description:
      Fit data to a func
      pbest=fit(xdata,ydata,func,pguess)
      
    Inputs:
    xdata: array -- measured abcisas
    ydata: array -- measured ordinates

    func: function -- function to be fitted.  
    Function prototype must be func(x,pars)

    pguess: array -- Guess values of parameters
    
    Outputs:
    pbest: array -- fitted values of parameters
    """
    info=leastsq(residuals,pguess,args=(xdata,ydata,func),full_output=1)
    
    pbest=info[0]
    covm=info[1]
    nfev=info[2]['nfev']
    fvec=info[2]['fvec']
    msg=info[3]
    ier=info[4]
    
    if ier in [1,2,3,4]:
        #print "Successful fit."
        #print msg
        None
    else:
        #print "Fit has failed" 
        #print msg
        None
        
    if verbose is True:
        output="""
The best parameters: %s
Covariance matrix:
%s
Number of function calls: %g
Average residual at the best fit: %g
"""%(str(pbest),str(covm),nfev,mag(fvec))
        #print output
        if logfile is not "":
            f=open(logfile,'w')
            print >>f,output
            f.close()
            #print "Information saved in ",logfile
    
    return pbest
