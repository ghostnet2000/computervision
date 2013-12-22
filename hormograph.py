"""
    hormograph.py 
"""

def normalize(points):
    """"
        Normalize a collection of points in homogeneous coordinates so that last row = 1
    """ 
    for row in points:
        row /= points[-1]
    return points
    
def make_homog(points):
    """
        convert a set of points (dim*array) to homogenous coordinates.
    """
    return vstack((points,ones((1,points.shape[1]))))