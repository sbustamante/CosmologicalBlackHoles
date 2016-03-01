#****************************************************************************************
#	DYNAMICAL FRICTION
#	Sebastian Bustamante (macsebas33@gmail.com)
#****************************************************************************************

#========================================================================================
#		IMPORTS
#========================================================================================
import numpy as np


#========================================================================================
#		GLOBAL VARIABLES
#========================================================================================
#Cavendish constant [SI]
GC	=	6.67384e-11


#========================================================================================
#		POTENTIALS
#========================================================================================
class hernquist(object):
    """Hernquist sphere
    
    Class of Hernquist sphere
    
    Attributes
    ----------
    M : total mass of the system
    a : length scale of the system
    units: dictionary with units to be used (transformed from SI)
    """
    
    
    def __init__( self, 
		  M, 
		  a, 
		  units = {'M':1, 'L':1, 'T':1} ):
	self.M = M
	self.a = a
	self.units = units
	self.GC = GC*units['M']*units['T']**2/units['L']**3
	
	
    def density( self, 
		 r ):
	"""
	Name: density
	Function: calculates the density of a hernquist sphere
	Arguments:
	    self: hernquist object
	    r: radial coordinate to be evaluated
	"""
	return self.M/( 2*np.pi )*( self.a/( r*(r + self.a)**3 ) )
	
	
    def potential( self, 
		   r ):
	"""
	Name: potential
	Function: calculates the potential of a hernquist sphere
	Arguments:
	    self: hernquist object
	    r: radial coordinate to be evaluated
	"""
	return -self.GC*self.M/( r + self.a )
      
      
    def v_circular( self, 
		    r ):
	"""
	Name: v_circular velocity ar radius r
	Function: calculates the circular velocity of a hernquist sphere
	Arguments:
	    self: hernquist object
	    r: radial coordinate to be evaluated
	"""
	return np.sqrt(self.GC*self.M*r)/( r + self.a )
