import numpy as np


# define default coefficients
v= 0.25    # drift term
lbg=0.1    # Specific algal maintenance respiration losses
mumax=1.2  # Maximum specific algal production rate
rhomax=0.2 # Maximum specific algal nutrient uptake rate
qmax=0.04  # Maximum algal nutrient quota
qmin=0.004 # Minimum algal nutrient quota 
m=1.5      # Half-saturation constant of algal nutrient uptake
h=120.0    # Half-saturation constant of light-dependent algal production
I0 = 300   # Light intensity at the surface 
kbg=0.4    # Background light-attenuation coefficient - change to 0.02 to see some oscillations? Originally 0.4
r = 0.02   # Specific mineralization rate of sedimented nutrients

def rnorm(mu,sigma,z):
    """ Use to create create normal hump around mu over z, of std. dev. sigma for initial conditions A_0,Rb_0,Rd_0."""
    return np.exp( -(z-mu)**2/(2*sigma**2) )/ np.sqrt(2*np.pi*sigma**2)

def I(z,A,I_0=I0, k = 0.0003, dz=0.1):
    """Function to plot I using array A[:,i], default k=0.0003, larger values of k make effect of A on I more apparent"""
    integral = np.zeros(len(z))
    integral[1:] = k*np.cumsum(A[1:])*dz
    return I_0 * np.exp( - integral - kbg*z)

def p(I,q):
    """Specific algal production (growth) rate"""
    return mumax * (1.0 - qmin/q) * (I/(h + I))

def rho(q, Rd):
    """Specific algal uptake ("feeding") rate"""
    return rhomax * (qmax-q)/(qmax-qmin) * ( Rd/(m + Rd) )

def next_step(z, A, Rb, Rd, Rs, dz, dt, d):
    """Calculates next step for input arrays of length zmax"""
    
    II = I(z,A)
    q = Rb[1:-1] / A[1:-1]
    pp = p(II[1:-1],q)
    rrho = rho(q,Rd[1:-1])
    
    A_next = np.zeros(len(A))
    Rb_next = np.zeros(len(Rb))
    Rd_next = np.zeros(len(Rd))
    
    A_drift = v * (A[2:]-A[:-2]) / (2  *dz)
    Rb_drift = v * (Rb[2:]-Rb[:-2]) / (2 * dz)
    
    A_diffusion = d * (A[2:]-2*A[1:-1] + A[:-2]) / (dz**2)
    Rb_diffusion = d * (Rb[2:]-2*Rb[1:-1] + Rb[:-2]) / (dz**2)
    Rd_diffusion = d * (Rd[2:]-2*Rd[1:-1] + Rd[:-2]) / (dz**2)
    
    A_next[1:-1] = A[1:-1] + dt * ( pp*A[1:-1] -lbg*A[1:-1] - A_drift + A_diffusion )
    A_next[0] = d *(4*A_next[1] - A_next[2] ) / (2*v*dz + 3*d)
    A_next[-1] = (4*A_next[-2] - A_next[-3]) / 3
    
    Rb_next[1:-1] = Rb[1:-1] + dt * (rrho * A[1:-1] -lbg * Rb[1:-1] - Rb_drift + Rb_diffusion )
    Rb_next[0] = d  * ( 4 * Rb_next[1] - Rb_next[2] ) / (2*v*dz + 3*d)
    Rb_next[-1] = (4*Rb_next[-2] - Rb_next[-3]) / 3
    
    Rs_next = Rs + dt*(v*Rb[-1] - r*Rs)
    
    Rd_next[1:-1] = Rd[1:-1] + dt*( lbg*Rb[1:-1] - rrho*A[1:-1] + Rd_diffusion )
    Rd_next[0] = (4 * Rd_next[1] - Rd_next[2]) / 3 
    Rd_next[-1] = (2*dz*r*Rs_next + d*(4*Rd_next[-2] - Rd_next[-3]))/(3*d)
    
    return A_next, Rb_next, Rd_next, Rs_next
    
    
def get_stationary(zmax=30, tmax=100, d=1.0, I0=300.0, dz=0.1):
    """ Obtains stationary distrbution via finite difference methods. \nA, Rb, Rd, Rs = equations.get_stationary(zmax=30, tmax=100, d=1.0, I0=300.0, dz=0.1) """
    
    dt = dz/1000
    Nz = int(zmax/dz)
    Nt = int(tmax/ (1000*dt) )

    z_grid = np.arange(0,zmax,dz)
    time_steps = np.arange(0,tmax,dt)

    A_0 = 100 
    Rb_0 = 2.2
    Rd_0 = 30 
    Rs_0 = 0
    
    # RESULTS MATRICES: rows: deeper z-values, cols: time steps forward
    A = np.zeros((Nz,Nt)) 
    Rb = np.zeros((Nz,Nt))
    Rd = np.zeros((Nz,Nt))
    
    # create homogenous initial conditions
    A_next = A_0 * np.ones(Nz) 
    Rb_next = Rb_0 * np.ones(Nz)
    Rd_next = Rd_0 * np.ones(Nz)
    Rs_next = 0
    # iterate next solution over time until tmax
    for t in time_steps:
        A_next, Rb_next, Rd_next, Rs_next = next_step(z_grid, A_next, Rb_next, Rd_next, Rs_next, dz, dt, d)
        
    # return final values
    return A_next, Rb_next, Rd_next,Rs_next

def get_time_evolution(zmax=10.0, tmax=10.0, d=1.0, I0=300.0, dz=0.1):
    """Return matrices showing time evolution of A, Rb, Rd and Rs. \nA, Rb, Rd, Rs, z_grid, time_steps, Nz, Nt = my_equations.get_time_evolution(zmax=10.0, tmax=10.0, d=1.0, I0=300.0, dz=0.1) \nNote: this function is quite slow."""
    
    dt = dz/1000
    Nz = int(zmax/dz)
    Nt = int(tmax/ (1000*dt) )

    z_grid = np.arange(0,zmax,dz)
    time_steps = np.arange(0,tmax,dt)

    A_0 = 100 
    Rb_0 = 2.2
    Rd_0 = 30 
    Rs_0 = 0
    
    # RESULTS MATRICES: rows: deeper z-values, cols: time steps forward
    A = np.zeros((Nz,Nt))
    Rb = np.zeros((Nz,Nt))
    Rd = np.zeros((Nz,Nt))
    Rs = np.zeros(Nt)
    
    # Create homogenous initial conditions with values as specified in Jaeger et al. 2010
    A[:,0] = A_0 * np.ones(Nz) 
    Rb[:,0] = Rb_0 * np.ones(Nz)
    Rd[:,0] = Rd_0 * np.ones(Nz)
    
    # begin loop
    A_next,Rb_next,Rd_next, Rs_next = A[:,0], Rb[:,0], Rd[:,0], Rs[0]
    i = 1
    counter = 1
    
    for t in time_steps[:-1]:
        A_next, Rb_next, Rd_next, Rs_next = next_step(z_grid, A_next, Rb_next, Rd_next, Rs_next, dz ,dt ,d)
        
        # Record every 1000th value
        if ( i % 1000 == 0):
            A[:,counter],Rb[:,counter],Rd[:,counter],Rs[counter] = A_next,Rb_next,Rd_next,Rs_next
            counter += 1
            
        i += 1
    
    # return result matrices and grid values to be used for plotting
    return A, Rb, Rd, Rs, z_grid, time_steps, Nz, Nt

def get_R(Rd,Rb):
    # iterate through rows
    R = np.zeros(Rd.shape)
    for i in range(len(Rd)):
       # iterate through columns
       for j in range(len(Rd[0])):
            R[i,j] = Rd[i,j] + Rb[i,j]
    return R
    