# This function is used to offset a datacube
#
# INPUT:     - img_file
#            - offset x, offset y (in pixels)
#            - offset can be determined using find_center.py
#
# OUTPUT:    - shifted datacube for combining 
#
# INPUT EXAMPLE:
#   
#   img, ih = fits.getdata("kb171021_00077_icuber.fits", header=True)
#   img_shift= cube_shift(img,shift_x,shift_y)
#
# EXAMPLE 2: 
#  
#   img1, ih1 = fits.getdata("kb171021_00077_icuber.fits", header=True)
#   img2, ih2 = fits.getdata("kb171021_00082_icuber.fits", header=True)
#
#   center1=find_center(img1)
#   center2=find_center(img2)
#   shift_x= center2[0]-center1[0]
#   shift_y= center2[1]-center1[1]
#
#   img2_shift= cube_shift(img2,shift_x,shift_y)
#   center_new= find_center(img2_shift)
#'''One can see new center for img2 is the same as that of img1. Ready for combining'''
#
#
# MODIFICATION HISTORY:
#   2017-11-02 Initial version (ZC)
#

def cube_shift(image,shift_x,shift_y):
    
    image_sh= image  # construct a new datacube called image_sh (image_shift)
    for i in np.arange(0, len(image[:,1,1]), 1):  # loop layers 
        # offset each layer using shift_x, shift_y 
        image_sh[i,:,:]= np.roll(image[i,:,:], -int(round(shift_x)), axis=0)
        image_sh[i,:,:]= np.roll(image[i,:,:], -int(round(shift_y)), axis=1)

    return image_sh
 
