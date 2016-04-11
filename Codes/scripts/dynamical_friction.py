#****************************************************************************************
#	DYNAMICAL FRICTION
#	Sebastian Bustamante (macsebas33@gmail.com)
#****************************************************************************************

#========================================================================================
#		IMPORTS
#========================================================================================
import numpy as np
from scipy import integrate as integ
import h5py


#========================================================================================
#		GLOBAL VARIABLES
#========================================================================================
#Cavendish constant [SI]
GC	=	6.67384e-11
#Standard units
unitsstd = {'M':1.989e40, 'L':3.08568e19, 'T':3.15569e16, 'V':1000. }


#========================================================================================
#		FUNCTIONS
#========================================================================================
#RK4 integrator
def RK4_step( f, r, v, t, dt ):
    #Creating solutions
    K0r, K0v = f( r, v, t )
    K1r, K1v = f( r + 0.5*dt*K0r, v + 0.5*dt*K0v, t + 0.5*dt )
    K2r, K2v = f( r + 0.5*dt*K1r, v + 0.5*dt*K1v, t + 0.5*dt )
    K3r, K3v = f( r + dt*K2r, v + dt*K2v, t + dt )
    rf = np.array(r + dt/6.0*( K0r + 2*K1r + 2*K2r + K3r ))
    vf = np.array(v + dt/6.0*( K0v + 2*K1v + 2*K2v + K3v ))

    #Returning solution
    return np.array([rf, vf])
    

#========================================================================================
#		INTERACTIONS
#========================================================================================
class hernquist_sphere(object):
    """Hernquist sphere
    
    Class of Hernquist sphere
    
    Attributes
    ----------
    M : total mass of the system
    a : length scale of the system
    units : dictionary with units to be used (transformed from SI). 
	   Default: M [1e10 Msun]  L [kpc]  T [Gyr]  V [km/s]
    kargs : extra arguments
    """
    def __init__( self, M, a, units = unitsstd, kargs={} ):
	self.M = M
	self.a = a
	self.units = units
	self.kargs = kargs
	
	#Units
	self.GC = GC*self.units['M']*self.units['T']**2/self.units['L']**3
	self.Vf = self.units['V']*self.units['T']/self.units['L']
	
	
    def density( self, r ):
	return self.M/( 2*np.pi )*( self.a/( r*(r + self.a)**3 ) )
	

    def potential( self, r ):
	return -self.GC*self.M/( r + self.a )
      

    def v_circular( self, r ):
	return np.sqrt(self.GC*self.M*r)/( r + self.a )/self.Vf


    def v_escape( self, r ):
	return np.sqrt(-2*self.potential(r) )/self.Vf


    def distribution_function( self, r, v ):
	#Velocity conversion
	v *= self.Vf
	#Specific energy
	E = 0.5*v**2 + self.potential( r )
	#q factor
	q = np.sqrt( -self.a*E/(self.GC*self.M) )
	#characteristic velocity
	vg = np.sqrt( self.GC*self.M/self.a )
	#distribution function
	f = self.M/( 8*np.sqrt(2)*np.pi**3*self.a**3*vg**3 )*1/( 1-q*q )**2.5*\
	( 3*np.arcsin(q) + q*np.sqrt(1 - q*q)*(1 - 2*q*q)*(8*q**4 - 8*q**2 - 3) )
	return f


    def force_hernquist( self, r, v=None, t=None ):
	#Norm of input vectors
	rm = np.linalg.norm( r )
	#Force
	force = -self.GC*self.M/( rm + self.a )**2*np.array(r)/rm

	return force


    def chandrasekhar_friction( self, r, v, t=None ):
	"""
	Chandrasekhar dynamical friction formula
	
	Required extra arguments
	------------------------
	Mb : mass of the body that experiences dynamical friction
	bmin : 90deg deflection radius
	bmax : maximum deflection radius
	"""
	#Extracting arguments
	self.Mb = self.kargs['Mb']
	#Coulomb Logarithm
	self.LogL = lambda bmin, bmax:  np.log( bmax/bmin )
	#Calculating minimum and maximum impact parameters
	self.bmin = lambda v: np.max( [ 2*self.GC*self.Mb/(3.0e5*self.Vf)**2, self.GC*self.Mb/(v*self.Vf)**2 ] ) #Tremmel 2015
	self.bmax = lambda r: r*(r + self.a)/( 4*r + self.a ) #Just, et al. 2011
	
	#Norm of input vectors
	rm = np.linalg.norm( r )
	vm = np.linalg.norm( v )
	#Contribution of slower particles
	if vm <= self.v_escape( rm )*self.Vf:
	    rho_slow = integ.quad( lambda vl: vl**2*self.distribution_function( rm, vl ), 0, vm )[0]
	else:
	    rho_slow = integ.quad( lambda vl: vl**2*self.distribution_function( rm, self.v_escape( rm )*self.Vf ), 0, vm )[0]
	#Dynamical friction
	a_dyn = -16*self.GC**2*np.pi**2*self.Mb*self.LogL( self.bmin(vm), self.bmax(rm) )*rho_slow*np.array(v)/vm**3

	return a_dyn


    def drag_friction( self, r, v, t=None ):
      	"""
	Drag dynamical friction formula
	
	Required extra arguments
	------------------------
	coef : friction coefficient
	"""
	self.coef = self.kargs['coef']
	
	#Norm of input vectors
	rm = np.linalg.norm( r )
	vm = np.linalg.norm( v )
	#Calculating dynamical time
	t_dyn = 1/np.sqrt( 4*np.pi*self.GC*self.density( rm ) )
	#Dynamical friction
	a_dyn = -self.coef*np.array(v)/t_dyn

	return a_dyn
	
      
#========================================================================================
#		BH
#========================================================================================      
class black_hole(object):
    """black hole
    
    Class with functions to evolve a BH
    
    Attributes
    ----------
    M : mass of the black hole
    r : position of the BH
    v : velocity of the BH
    acc : accretion mode
    force : forces acting on the BH
    units : dictionary with units to be used (transformed from SI). 
	   Default: M [1e10 Msun]  L [kpc]  T [Gyr]
    """  
    def __init__( self, M, r, v, acc, force, units = unitsstd ):
	self.M = M
	self.r = np.array(r)
	self.t = 0
	self.acc = acc
	self.force = force
	self.units = units
	
	#Units
	self.GC = GC*self.units['M']*self.units['T']**2/self.units['L']**3
	self.Vf = self.units['V']*self.units['T']/self.units['L']
	self.v = np.array(v)/self.Vf

	#Creating trayectory of the BH
	self.rt = [self.r]
	self.vt = [self.v]
	self.time = [self.t]
	
	
    def f_evol( self, r, v, t ):
	return np.array([np.array(v), np.array(self.force( r, v, t ))])
	
	
    def time_step( self, dt, scheme='rk4' ):
	if scheme == 'rk4':
	    rf, vf = RK4_step( self.f_evol, self.r, self.v, self.t, dt )

	#Updating variables
	self.r = rf
	self.v = vf
	self.t += dt
	#Updating trajectories
	self.rt.append( self.r )
	self.vt.append( self.v )
	self.time.append( self.t )
	
	
#========================================================================================
#		BH_SIMULATED
#========================================================================================
class black_hole_sim(object):
    """black hole sim
    
    Class of simulated black holes
    
    Attributes
    ----------
    simulation : name of simulation
    snapbase : snap base 
    n_snap : number of snapshots
    datafolder : folder where simulations are stored
    resultsfolder : folder where trajectories will be stored
    center: center of the simulation
    units : dictionary with units to be used (transformed from SI). 
	   Default: M [1e10 Msun]  L [kpc]  T [Gyr]  V [km/s]
    kargs : extra arguments
    """
    def __init__( self, simulation, snapbase, n_snap, datafolder, resultsfolder, center = [0,0,0], units = unitsstd, kargs={} ):
	self.simulation = simulation
	self.snapbase = snapbase
	self.n_snap = n_snap
	self.datafolder = datafolder
	self.resultsfolder = resultsfolder
	
	self.center = np.array(center)	
	self.units = units
	self.kargs = kargs
	
	#Units
	self.GC = GC*self.units['M']*self.units['T']**2/self.units['L']**3
	self.Vf = self.units['V']*self.units['T']/self.units['L']
	
    
    def read_trajectory( self, total_energy = False ):
	#Snap function
	filename = lambda snap : '%s/%s/output/%s_%03d.hdf5'%( self.datafolder, self.simulation, self.snapbase, snap )
	#Initializing position and velocity
	self.r = []
	self.v = []
	self.r_mp = []
	self.v_mp = []
	self.t = []
	for i in xrange( self.n_snap ):
	    datafile = h5py.File(filename(i), 'r')
	    #Reading information of minimum of potential
	    i_min = np.argsort( datafile['PartType1']['Potential'] + total_energy*0.5*np.linalg.norm(datafile['PartType1']['Velocities'], axis=1) )[0]
	    self.r_mp.append( datafile['PartType1']['Coordinates'][i_min] - self.center )
	    self.v_mp.append( datafile['PartType1']['Velocities'][i_min] )
	    #Reading information of black hole trajectory
	    self.r.append( datafile['PartType5']['Coordinates'][0] - self.center )
	    self.v.append( datafile['PartType5']['Velocities'][0] )
	    self.t.append( datafile['Header'].attrs['Time'] )

	self.r = np.array( self.r )
	self.rm = np.linalg.norm( self.r, axis=1 )
	self.v = np.array( self.v )
	self.vm = np.linalg.norm( self.v, axis=1 )
	
	self.r_mp = np.array( self.r_mp )
	self.rm_mp = np.linalg.norm( self.r_mp, axis=1 )
	self.v_mp = np.array( self.v_mp )
	self.vm_mp = np.linalg.norm( self.v_mp, axis=1 )
	
	self.t = np.array(self.t)
	#Centralizing trajectory
	self.centralization( )
	#Saving trajectory
	self.save_trajectory( )
	
	
    def centralization( self ):
	#This function centralizes the trajectory of the BH with respect to the minimum of the potential
	self.r_c = self.r - self.r_mp
	self.rm_c = np.linalg.norm( self.r_c, axis=1 )
	self.v_c = self.v - self.v_mp
	self.vm_c = np.linalg.norm( self.v_c, axis=1 )
	

    def save_trajectory( self ):
	#Saving data from BH and minimum of potential
	data = np.transpose( [ self.t,
	self.r[:,0], self.r[:,1], self.r[:,2], self.rm,
	self.v[:,0], self.v[:,1], self.v[:,2], self.vm,
	self.r_mp[:,0], self.r_mp[:,1], self.r_mp[:,2], self.rm_mp,
	self.v_mp[:,0], self.v_mp[:,1], self.v_mp[:,2], self.vm_mp ] )
	
	np.savetxt( "%s/BH_%s.dat"%(self.resultsfolder, self.simulation), data )
	
	
    def load_trajectory( self ):
	data = np.loadtxt( "%s/BH_%s.dat"%(self.resultsfolder, self.simulation) )

	self.t = data[:,0]
	
	self.r = data[:,[1,2,3]]
	self.rm = data[:,4]
	self.v = data[:,[5,6,7]]
	self.vm = data[:,8]
	
	self.r_mp = data[:,[9,10,11]]
	self.rm_mp = data[:,12]
	self.v_mp = data[:,[13,14,15]]
	self.vm_mp = data[:,16]
	
	#Centralizing trajectory
	self.centralization( )
	
	
    def noise_filter( self ):
	#Function that filters noisy points using numerical estimative
	