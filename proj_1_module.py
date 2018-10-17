import numpy as np

def dynamics_solve(f, D = 1, t_0 = 0.0, s_0 = 1, h = 0.1, N = 100, method = "Euler"):
    
    """ Solves for dynamics of a given dynamical system
    
    - User must specify dimension D of phase space.
    - Includes Euler, RK2, RK4, that user can choose from using the keyword "method"
    
    Args:
        f: A python function f(t, s) that assigns a float to each time and state representing
        the time derivative of the state at that time.
        
    Kwargs:
        D: Phase space dimension (int) set to 1 as default
        t_0: Initial time (float) set to 0.0 as default
        s_0: Initial state (float for D=1, ndarray for D>1) set to 1.0 as default
        h: Step size (float) set to 0.1 as default
        N: Number of steps (int) set to 100 as default
        method: Numerical method (string), can be "Euler", "RK2", "RK4"
    
    Returns:
        T: Numpy array of times
        S: Numpy array of states at the times given in T
    """
    
    T = np.array([t_0 + n * h for n in range(N + 1)])
    
    if D == 1:
        S = np.zeros(N + 1)
    
    if D > 1:
        S = np.zeros((N + 1, D))
        
    S[0] = s_0
    
    if method == 'Euler':
        for n in range(N):
            S[n + 1] = S[n] + h * f(T[n], S[n])
    
    if method == 'RK2':
        for n in range(N):
            k1 = 
            k2 = 
            S[n + 1] = 
    
    if method == 'RK4':
        for n in range(N):
            k1 = 
            k2 = 
            k3 = 
            k4 = 
            S[n + 1] = 
            
    return T, S

def hamiltonian_solve(d_qH, d_pH, d = 1, t_0 = 0.0, q_0 = 0.0, p_0 = 1.0, h = 0.1, N = 100, method = "Euler",):
    
    """ Solves for dynamics of Hamiltonian system
    
    - User must specify dimension d of configuration space.
    - Includes Euler, RK2, RK4, Symplectic Euler (SE) and Stormer Verlet (SV) 
      that user can choose from using the keyword "method"
    
    Args:
        d_qH: Partial derivative of the Hamiltonian with respect to coordinates (float for d=1, ndarray for d>1)
        d_pH: Partial derivative of the Hamiltonian with respect to momenta (float for d=1, ndarray for d>1)
        
    Kwargs:
        d: Spatial dimension (int) set to 1 as default
        t_0: Initial time (float) set to 0.0 as default
        q_0: Initial position (float for d=1, ndarray for d>1) set to 0.0 as default
        p_0: Initial momentum (float for d=1, ndarray for d>1) set to 1.0 as default
        h: Step size (float) set to 0.1 as default
        N: Number of steps (int) set to 100 as default
        method: Numerical method (string), can be "Euler", "RK2", "RK4", "SE", "SV"
    
    Returns:
        T: Numpy array of times
        Q: Numpy array of positions at the times given in T
        P: Numpy array of momenta at the times given in T
    """
    T = 
    
    if d == 1:
        P = 
        Q = 
        
        Q[0] = 
        P[0] =
    
    if d > 1:
        
    
    if method == 'Euler':
        for n in range(N):
            Q[n + 1] = Q[n] + h * d_pH(Q[n], P[n])
            P[n + 1] = 
    
    if method == 'RK2':
        for n in range(N):
            k1_Q = h * d_pH(Q[n], P[n])
            k1_P = h * (- d_qH(Q[n], P[n]))
            
            k2_Q = 
            k2_P = 
            
            Q[n + 1] = 
            P[n + 1] = 
        
    if method == 'RK4':
        for n in range(N): 
            
        
    if method == 'SE':
        
    
    if method == "SV":
        
        
    return T, Q, P
