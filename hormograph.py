"""
    hormograph.py 
"""

from numpy import *

def normalize(points):
    """"
        Normalize a collection of points in homogeneous coordinates so that last row = 1
    """ 
    for row in points:
        row /= points[-1]
        print row
    return points
    
def make_homog(points):
    """
        convert a set of points (dim*array) to homogenous coordinates.
    """
    print ones
    
    return vstack((points,ones((1,points.shape[1]))))

def H_from_points(fp,tp):
    """ Find hormography H, such that fp is mapped to tp using the linear DLT method. Points are conditioned automatically
    """
    
    if fp.shape != tp.shape:
        raise RuntimeError('Number of points do not match')
        
    # conditioned points (important for numerical reasons)
    # -- from points --
    m = mean(fp[:2], axis=1)
    maxstd = max(std(fp[:2], axis=1)) + 1e-9
    C1 = diag([1/maxstd, 1/maxstd, 1])
    C1[0][2] = -m[0]/maxstd
    C1[1][2] = -m[1]/maxstd
    
    fp = dot(C1,fp)
    
    # --- to points
    m = mean(tp[:2], axis=1)
    maxstd = max(std(tp[:2], axis=1)) + 1e-9
    C2 = diag([1/maxstd, 1/maxstd, 1])
    C2[0][2] = -m[0]/maxstd
    C2[1][2] = -m[1]/maxstd
    
    tp = dot(C2,tp)
    
    # create matrix for linear method, 2 rows for each correspondence pair
    nbr_correspondences = fp.shape[1]
    A = zeros((2*nbr_correspondences,9)) 
    for i in range(nbr_correspondences):
        A[2*i] = [-fp[0][i],-fp[1][i],-1,0,0,0, tp[0][i]*fp[0][i],tp[0][i]*fp[1][i],tp[0][i]]
        A[2*i+1] = [0,0,0,-fp[0][i],-fp[1][i],-1, tp[1][i]*fp[0][i],tp[1][i]*fp[1][i],tp[1][i]]
    
    U,S,V = linalg.svd(A)
    H = V[8].reshape((3,3))
    
    # decondition
    H = dot(linalg.inv(C2),dot(H,C1))
    
    # normalize and return
    return H / H[2,2]

#"""
def main():
    
    p1 = array([1, 4, 3, 1, 4, 3, 1, 4, 3])
    p2 = array([2, 5, 2, 1, 4, 3, 1, 4, 3])
    
    
    H_from_points(p1,p2)
    
if __name__ == "__main__":
    main()
#"""