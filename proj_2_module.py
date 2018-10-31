import numpy as np
from numpy.linalg import norm

def hermitian_eigensystem(H, tolerance):
    
    """ Solves for the eigenvalues and eigenvectors of a hermitian matrix
    
    Args:
        H: Hermitian matrix for which we want to compute eigenvalues and eigenvectors
        
        tolerance: A number that sets the tolerance for the accuracy of the computation.  This number
        is multiplied by the norm of the matrix H to obtain a number delta.  The algorithm successively
        applies (via similarity transformation) Jacobi rotations to the matrix H until the sum of the
        squares of the off-diagonal elements are less than delta.
    
    
    
    Returns:
        d: Numpy array containing eigenvalues of H in non-decreasing order
        
        U: A 2d numpy array whose columns are the eigenvectors corresponding to the computed
        eigenvalues.
        
        
    Checks you might need to do:
        
        H * U[:,k] = d[k] *　U[:,k]      k=0,1,2,...,(n-1)
        
        d[0] <= d[1] <= ... <= d[n-1]     (where n is the dimension of H)
       
        np.transpose(U) * U = U * np.transpose(U) = np.eye(n)
        
    """
    
    return d, U
