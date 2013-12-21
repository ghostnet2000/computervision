from PIL import Image
from pylab import *
from scipy.ndimage import filters

def compute_harris_response(im, sigma=3):
    """
    Compute the Harris corner detector response function for each pixel in a graylevel image
    """
    #derivatives
    imx = zeros(im.shape)
    filters.gaussian_filter(im, (sigma,sigma), (0,1), imx)
    imy = zeros(im.shape)
    filters.gaussian_filter(im, (sigma,sigma), (1,0), imy)
    
    Wxx = filters.gaussian_filter(imx*imx,sigma)
    Wxy = filters.gaussian_filter(imx*imy,sigma)
    Wyy = filters.gaussian_filter(imy*imy,sigma)

    Wdet = Wxx*Wyy - Wxy**2
    Wtr = Wxx + Wyy
    
    return (Wdet / Wtr)
    
def get_harris_points(harrisim,min_dist=10,threshold=0.1):
    """
    """
    #find top corner candidates above a threshold
    corner_threshold = harrisim.max() * threshold
    
    harrisim_t = (harrisim > corner_threshold) * 1
    
    # get coordinates of candidates
    coords = array(harrisim_t.nonzero()).T
    
    # ... and their values
    candidate_values = [harrisim[c[0], c[1]] for c in coords]
    print candidate_values
    
    #sort candidates
    index = argsort(candidate_values)
    
    #stored allowed point location in array
    allowed_locations = zeros(harrisim.shape)
    allowed_locations[min_dist:-min_dist, min_dist:-min_dist ] = 1
    
    # select the best points talking min_distance into account
    filtered_coords = []
    for i in index:
        if allowed_locations[coords[i, 0], coords[i, 1]] == 1:
            filtered_coords.append(coords[i])
            allowed_locations[(coords[i,0]-min_dist):(coords[i,0]+min_dist),(coords[i,1]-min_dist):(coords[i,1]+min_dist)] = 0

    return filtered_coords


def get_descriptors(image,filtered_coords, wid=5):
    """ For each point return, pixel values around the point
        using a  neighbourhood of width 2*wid+1
        Assume points are extracted with min_distance > wid
    
    """  
    desc = []
    for coords in filtered_coords:
        patch = image[coords[0]-wid:coords[0]+wid+1,coords[1]-wid:coords[1]+wid+1].flatten()
        desc.append(patch)
        
    return desc
    
def match(desc1, desc2,threshold=0.5):
    """
    For each corner point descriptor in the first image, select its match to the second image using normalized cross-correlation.
        
    """
    n = len(desc1[0])
    
    # pair-wise distance
    d = -ones((len(desc1),len(desc2)))
    for i in range(len(desc1)):
        for j in range(len(desc2)):
            d1 = (desc1[i] - mean(desc1[i])) / std(desc1[i])
            d2 = (desc2[j] - mean(desc2[j])) / std(desc2[j])
            ncc_value = sum(d1 * d2) / (n-1)
            if ncc_value > threshold:
                d[i,j] = ncc_value
    ndx = argsort(-d)
    matchscores = ndx[:,0]
    
    return matchscores
    
def match_twosided(desc1,desc2,threshold=0.5)
    """ 
    Two sided symmetric version of match(). 
    """
    matches_12 = match(desc1,desc2,threshold)
    matches_21 = match(desc2,desc1,threshold)
    
    ndx_12 = where(matches_12 >= 0)[0]
    
    #remove matches that are not symmetric
    for n in ndx_12:
        if matches_21[matches_12[n]] !n n:
            matches_12[n] = -1
    
    return matches_12[n]
    
def appendimages(im1,im2):

    """ Return a new image that appends the two images side-by-side. """
    
    # select the image with the fewest rows and fill in enough empty rows
    row1 = im1.shape[0]
    row2 = im2.shape[0]
    
    if row1 < row2:
        im1 = concatenate((im1,zeros((rows2-rows1,im1.shape[1]))), axis=0)
    elif rows1 > rows2:
        im1 = concatenate((im2,zeros((rows2-rows1,im2.shape[1]))), axis=0)
    #if none of these case they are equal
    
    return concatenate((im1,im2), axis=1)
            
            
    
def plot_harris_points(image, filtered_coords):
    """ Plots corners found in image. """

    figure()
    gray()
    imshow(image)
    plot([p[1] for p in filtered_coords], [p[0] for p in filtered_coords], '*')
    axis('off')
    show()
       
    
    
def main():
    im = array(Image.open('nelson.jpg').convert('L'))
    harrisim = compute_harris_response(im)
    filtered_coords = get_harris_points(harrisim,6)
    plot_harris_points(im, filtered_coords)

if __name__ == "__main__":
    main()

