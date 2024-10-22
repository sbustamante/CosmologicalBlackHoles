#****************************************************************************************
#	DYNAMICAL FRICTION
#	Sebastian Bustamante (macsebas33@gmail.com)
#****************************************************************************************

#========================================================================================
#		IMPORTS
#========================================================================================
import numpy as np
from scipy import integrate as integ
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from astropy.cosmology import WMAP9
import h5py

try:
    execfile( './scripts/fitool.py' )
except:
    None

#========================================================================================
#		GLOBAL VARIABLES
#========================================================================================
#Cavendish constant [SI]
GC	=	6.67384e-11
#Light speed
CLIGHT	=	2.99792458e8
#Proton Mass
PROTONMASS =    1.67262178e-27
#Thomson constant
THOMPSON =      6.65245873e-29
#1e8 Solar units
E8MSUN = 	1.989e38
#Seconds per year
SECPERYEAR =    3.15569e7
#Solar mass
SOLARMASS =     1.989e30


#Standard units
unitsstd = {'M':1.989e40, 'L':3.08568e19, 'T':3.15569e16, 'V':1000. }

#Creating time for cosmological setups
atmp = np.arange( 0.001, 1.02, 0.001)
t_tmp = WMAP9.lookback_time(1.0/atmp[0] - 1) - WMAP9.lookback_time( 1.0/atmp - 1 )
def physical_time( a ):
    time = np.interp( a, atmp, t_tmp)
    return time

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
  
  
#Numerical derivative
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
  
  
#Function to load coordinates and potential from a snapshot
def snapshot_reader( datafolder, simulation, snap, snapbase = 'snapshot', parttype = 1, BH_offset = False, units = unitsstd ):
  
    filename = lambda snap : '%s/%s/output/%s_%03d.hdf5'%( datafolder, simulation, snapbase, snap )
    datafile = h5py.File(filename(snap), 'r')
    #Calculating potential vector due to central BH
    coordinates = datafile['PartType%d'%(parttype)]['Coordinates']
    potential = datafile['PartType%d'%(parttype)]['Potential']

    r_mp = []
    if BH_offset:
        #Units
	GCl = GC*units['M']*units['T']**2/units['L']**3
	#Potential of BH
	dist = np.linalg.norm( coordinates - datafile['PartType5']['Coordinates'][0], axis=1 )
	potBH = -GCl*datafile['PartType5']['Masses'][0]/dist
	#Reading information of minimum of potential
	i_min = np.argsort( potential - potBH )[0]
	r_mp = coordinates[i_min]
	potential = np.array(potential)-potBH
      
    #Finding ID of most bounded particle
    id_bound = np.argsort( potential )[0]
    
    try:
	return np.array(coordinates), np.array(potential), id_bound, r_mp, datafile['PartType5']['Coordinates'][0]
    except:
	return np.array(coordinates), np.array(potential), id_bound, r_mp, []
  
  
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
    def __init__( self, simulation, snapbase, n_snap, datafolder, resultsfolder, center = [0,0,0], comoving = False, units = unitsstd, kargs={} ):
	self.simulation = simulation
	self.snapbase = snapbase
	self.n_snap = n_snap
	self.datafolder = datafolder
	self.resultsfolder = resultsfolder
	
	self.comoving = comoving
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
	self.a = []
	self.adir = []
	self.mbh = []
	self.mbh_dot = []
	
	#Reference Index
        datafile0 = h5py.File(filename(0), 'r')
        IDs = datafile0['PartType5']['ParticleIDs']
        
	for i in xrange( self.n_snap ):
	    datafile = h5py.File(filename(i), 'r')
	    
	    #Time
	    if self.comoving:
		self.t.append( physical_time(datafile['Header'].attrs['Time']) )
	    else:
		self.t.append( datafile['Header'].attrs['Time'] )
		
	    if len(datafile['PartType5']['ParticleIDs'])>1:
		mask_ID1 = np.array(datafile['PartType5']['ParticleIDs']) == IDs[1]               
		_, mask_ID1_3D = np.broadcast_arrays(datafile['PartType5']['Coordinates'], mask_ID1[...,None])
		if sum(mask_ID1)==0:
		    mask_ID1 = 0  
		    mask_ID1_3D = 0
	    else:
		mask_ID1 = 0
		mask_ID1_3D = 0
		
	    #Calculating potential vector due to central BH
	    dist = np.linalg.norm( datafile['PartType1']['Coordinates'] - datafile['PartType5']['Coordinates'][mask_ID1_3D], axis=1 )
	    potBH = -self.GC*datafile['PartType5']['Masses'][mask_ID1]/dist
	    
	    #Reading information of minimum of potential
	    i_min = np.argsort( datafile['PartType1']['Potential'] + total_energy*0.5*np.linalg.norm(datafile['PartType1']['Velocities'], axis=1) - potBH )[0]
	    self.r_mp.append( datafile['PartType1']['Coordinates'][i_min] - self.center )
	    self.v_mp.append( datafile['PartType1']['Velocities'][i_min] )
	    #Reading information of black hole trajectory
	    self.r.append( datafile['PartType5']['Coordinates'][mask_ID1_3D] - self.center )
	    self.v.append( datafile['PartType5']['Velocities'][mask_ID1_3D] )
	    
	    self.mbh.append( datafile['PartType5']['BH_Mass'][mask_ID1] )
	    self.mbh_dot.append( datafile['PartType5']['BH_Mdot'][mask_ID1] )

	    try:
		self.a.append( datafile['PartType5']['BH_SpinParameter'][mask_ID1] )
		self.adir.append( np.array(datafile['PartType5']['BH_SpinOrientation'][mask_ID1_3D]) )
	    except:
		self.a.append( 0 )
		self.adir.append( [0,0,0] )

	self.r = np.array( self.r )
	self.rm = np.linalg.norm( self.r, axis=1 )
	self.v = np.array( self.v )
	self.vm = np.linalg.norm( self.v, axis=1 )
	
	self.r_mp = np.array( self.r_mp )
	self.rm_mp = np.linalg.norm( self.r_mp, axis=1 )
	self.v_mp = np.array( self.v_mp )
	self.vm_mp = np.linalg.norm( self.v_mp, axis=1 )
	#Reading spin	
	self.a = np.array( self.a )
	self.adir = np.array( self.adir )
	self.mbh = np.array( self.mbh )
	self.mbh_dot = np.array( self.mbh_dot )
	
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

	
    def centralization_fit( self ):
	#This function centralizes the trajectory of the BH with respect to the fitted minimum of the potential
	self.r_cf = self.r - self.r_mp_fit
	self.rm_cf = np.linalg.norm( self.r_cf, axis=1 )


    def save_trajectory( self ):
	#Saving data from BH and minimum of potential
	data = np.transpose( [ self.t,
	self.r[:,0], self.r[:,1], self.r[:,2], self.rm,
	self.v[:,0], self.v[:,1], self.v[:,2], self.vm,
	self.r_mp[:,0], self.r_mp[:,1], self.r_mp[:,2], self.rm_mp,
	self.v_mp[:,0], self.v_mp[:,1], self.v_mp[:,2], self.vm_mp, 
	self.a, self.mbh, self.mbh_dot,
	self.adir[:,0], self.adir[:,1], self.adir[:,2]] )
	
	np.savetxt( "%s/BH_%s.dat"%(self.resultsfolder, self.simulation), data )
	
	
    def load_trajectory( self, Distance = 3 ):
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
	
	try:
	    self.a = data[:,17]
	    self.mbh = data[:,18]
	    self.mbh_dot = data[:,19]
	    self.adir = data[:,[20,21,22]]
	except:
	    None
	
	try:
	    self.r_mp_fit = np.loadtxt( "%s/BH_%s_fit_%1.3f.dat"%(self.resultsfolder, self.simulation, Distance) )
	except:
	    None
	    
	try:
	    self.r_mp_fit_no = np.loadtxt( "%s/BH_%s_fit_no_%1.3f.dat"%(self.resultsfolder, self.simulation, Distance) )
	except:
	    None
	    
	try:
	    self.r_mp_3Dfit = np.loadtxt( "%s/BH_%s_3Dfit.dat"%(self.resultsfolder, self.simulation) )
	except:
	    None
	
	
    def fitted_potential_minimum( self, BH_offset = False, Nbins = 30, Distance = 3 ):
	#Function to track orbit of potential minimum
    
    	filename = lambda snap : '%s/%s/output/%s_%03d.hdf5'%( self.datafolder, self.simulation, self.snapbase, snap )
    
	#Fitting Function
	def function( x, args ):
	    A = args[0]
	    B = args[1]
	    C = args[2]
	    return A*x**2 + B*x + C
	
	#Looping throughout all snapshots
	self.contours = [] 
	self.r_mp_fit = []
	self.opt_args = []
	self.t = []

	for snap in xrange(self.n_snap):
	    CoordinatesC, PotentialC, id_boundC, r_mpC, r_bhC = snapshot_reader( self.datafolder, self.simulation, snap, BH_offset=BH_offset )

	    datafile = h5py.File(filename(snap), 'r')
	    self.t.append( datafile['Header'].attrs['Time'] )
	    #Filtering particles closer to the box's center
	    CenterC = CoordinatesC[id_boundC]
	    #Distance from most bounded particle
	    Distance = 3
	    #Properties
	    mask = np.linalg.norm( CoordinatesC - 1*CenterC, axis=1 ) <= Distance
	    CoordClsC = CoordinatesC[mask]
	    PotClsC = PotentialC[mask]
	    #Finding most bounded particles delimiting potential in every projection
	    RpsC = []
	    PpsC = []
	    for ip in xrange(3):
		#Building grid
		R = []
		P = []
		for ig in xrange(Nbins):
		    coords_shift = CoordClsC[:,ip]-(CenterC[ip]-Distance)
		    mask = (2.0*ig*Distance/Nbins <= coords_shift)*(coords_shift < 2.0*(ig+1)*Distance/Nbins)
		    try:
			id_min = np.argsort(PotClsC[mask])[0]
			R.append( CoordClsC[mask][id_min,ip] )
			P.append( PotClsC[mask].min() )
		    except:
			None
		R = np.array( R )
		P = np.array( P )
		RpsC.append( R )
		PpsC.append( P )
	    RpsC = np.array( RpsC )
	    PpsC = np.array( PpsC )
	    #Fitting every projection
	    r_mpOpt = []
	    OptArgs = []
	    for ip in xrange(3):
		OptArg = fit(RpsC[ip], PpsC[ip],function,[1076.,-538000,67000000]) 
		OptArgs.append( OptArg )
		r_mpOpt.append(-OptArg[1]/(2*OptArg[0]))
	    r_mpOpt = np.array(r_mpOpt) - self.center
	    self.contours.append( [RpsC.T  - self.center,PpsC] )
	    self.r_mp_fit.append( r_mpOpt )
	    self.opt_args.append( OptArgs )
	self.contours = np.array( self.contours )
	self.r_mp_fit = np.array( self.r_mp_fit )
	self.opt_args = np.array( self.opt_args )
	self.t = np.array(self.t)

	if BH_offset:
	    np.savetxt( "%s/BH_%s_fit_%1.3f.dat"%(self.resultsfolder, self.simulation, Distance ), self.r_mp_fit )
	else:
	    np.savetxt( "%s/BH_%s_fit_no_%1.3f.dat"%(self.resultsfolder, self.simulation, Distance ), self.r_mp_fit )
	
	
    def fitted_potential_minimum3D( self, BH_offset = False, Distance = 3 ):
	#Function to track orbit of potential minimum
    
	#Fitting Function
	def fitFunc(x, a, x0, x1, x2, y0):
	    return a*( (x[0]-x0)**2 + (x[1]-x1)**2 + (x[2]-x2)**2 ) + y0
	
	#Looping throughout all snapshots
	self.r_mp_fit = []
	for snap in xrange(self.n_snap):
	    CoordinatesC, PotentialC, id_boundC, r_mpC, r_bhC = snapshot_reader( self.datafolder, self.simulation, snap, BH_offset=BH_offset )
	    #Filtering particles closer to the box's center
	    CenterC = CoordinatesC[id_boundC]
	    #Distance from most bounded particle
	    Distance = 3
	    #Properties
	    mask = np.linalg.norm( CoordinatesC - 1*CenterC, axis=1 ) <= Distance
	    CoordClsC = CoordinatesC[mask]
	    PotClsC = PotentialC[mask]
	    fitParams, fitCovariances = curve_fit(fitFunc, CoordClsC.T, PotClsC, [ 1076., 250., 250., 250., -210000 ])

	    self.r_mp_fit.append( fitParams[1:4] - self.center )
	    
	self.r_mp_fit = np.array( self.r_mp_fit )
	np.savetxt( "%s/BH_%s_3Dfit.dat"%(self.resultsfolder, self.simulation), self.r_mp_fit )

#........................................................................................
#		SPIN_EVOLUTION
#........................................................................................
    def m_dot_eddignton( self, Mbh, rad_eff = 0.1 ):
	return (4*np.pi*GC*CLIGHT*PROTONMASS/(rad_eff*CLIGHT**2*THOMPSON))*Mbh*self.units['T']
      
    def spin_parameter_evolution( self, Mbh0, Mbhf, a0 ):
	Z1 = 1 + (1-a0**2)**(1/3.0)*( (1+a0)**(1/3.0) + (1-a0)**(1/3.0) )
	Z2 = ( 3*a0**2 + Z1**2 )**(0.5)
	Rlso = 3 + Z2 - abs(a0)/a0*( (3-Z1)*(3+Z1+2*Z2) )**(0.5)
      
	if Mbhf/Mbh0 > Rlso**0.5:
	    return 0.998
	else:
	    return (1/3.0)*Rlso**0.5*(Mbh0/Mbhf)*( 4 - (3*Rlso*(Mbh0/Mbhf)**2 - 2)**0.5 )
      
      
    def spin_evolution( self, mode='continuous', rad_eff = 0.1, alpha = 0.1, a0 = 0.5 ):
	"""
	Spin evolution function
	
	Required arguments
	------------------------
	mode: mode for spin evolution. It can be either 'continuous' or 'chaotic'.
	rad_eff: radiative efficiency of the disk
	alpha: Shakura-Sunyaev parameter
	a0: initial spin parameter
	"""
	#Disk properties
	self.rad_eff = rad_eff
	self.alpha = alpha
	
	#Interpolating functions of mass and m_dot
	self.Mbh = interp1d( self.t, self.mbh )
	self.Mbh_dot = interp1d( self.t, self.mbh_dot )
	Mbh_inv = interp1d( self.mbh, self.t )
	
	#Iterating in time
	aps = []
	tms = []
	dms = []
	t = self.t[0]
	a = a0
	i = 0
	dm = 0
	while t < self.t[-1] and i<50000:
	  #Properties at this time
	  lambd = self.Mbh_dot(t)/self.m_dot_eddignton( self.Mbh(t), self.rad_eff )
	  mass = self.Mbh(t)*self.units['M']/E8MSUN
	  time_acc = 3e6*abs(a)**(7/8.)*mass**(11/8.)*lambd**(-3/4.)*self.alpha**(1/4.)*SECPERYEAR/self.units['T']
	  Mself = 2.13e5*(self.rad_eff/lambd)**(-5/27.)*mass**(23/27.)*self.alpha**(-2/17.)*SOLARMASS/self.units['M']
	  Rself = 1.5e3*(self.rad_eff/lambd)**(8/27.)*mass**(-26/27.)*self.alpha**(14/27.)
	  #Storing properties
	  aps.append( a )
	  tms.append( t )
	  dms.append( dm )
	  i += 1
	  if lambd >= 0.001:
	      if mode == 'continuous':
		  try:
		      a = self.spin_parameter_evolution( self.Mbh(t), self.Mbh(t + time_acc), a )
		      dm = self.Mbh(t + time_acc)-self.Mbh(t)
		      if self.Mbh(t) != 0:
			  t = t + time_acc
		      else:
			  t = t+0.001*self.t[-1]
		  except:
		      break
	      if mode == 'chaotic':
		  try:
		      Jb = a*mass**2
		      Jd = 3.6*abs(a)**(19/16.)*mass**(55/16.)*lambd**(1/8.)*self.alpha**(5/8.)
		      Jd = mass*Mself*Rself**0.5
		      dm = Mself
		      
		      if np.random.random()<(0.5*(1 - 0.5*Jd/Jb)):#np.min( (abs(0.5*(1 - 0.5*Jd/Jb)),0.5) ):
			  a = -self.spin_parameter_evolution( self.Mbh(t), self.Mbh(t) + Mself, -a )
		      else:
			  a = self.spin_parameter_evolution( self.Mbh(t), self.Mbh(t) + Mself, a )
			  
		      if self.Mbh(t) != 0:
			  t = Mbh_inv( self.Mbh(t) + Mself )
		      else:
			  t = t+0.001*self.t[-1]
		  except:
		      break
	  else:
	      t = t+0.001*self.t[-1]
	self.aps = np.array(aps)
	self.tms = np.array(tms)
	self.dms = np.array(dms)
	      