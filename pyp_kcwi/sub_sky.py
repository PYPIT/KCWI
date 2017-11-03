""" Functions to subtract sky:
"""

def median(sky_file, img_file, mask_limit, cut_ch=300):
    """ This function is used to subtract sky median

    Parameters
    ----------
    sky_file
    img_file
    mask_limit : mask limit
    cut_ch     : cut bad channel number

    Eaxmple:
    median('kb171021_00083_icuber.fits', 'kb171021_00082_icuber.fits', 40, cut_ch=300)

    Returns
    -------
    mask.fits  : masked image (mask the bad pixels and qso)
    img.fits   : sky median subtracted image
    """
    from astropy.io import fits
    import numpy as np
    from astropy.stats import sigma_clip
    sky_input, sh = fits.getdata(sky_file, header=True)
    img_input, ih = fits.getdata(img_file, header=True)

    # goodwave layers, default clip 300 layers at the beginning and the end
    img=img_input[cut_ch:(len(img_input[:,1,1])-cut_ch),:,:]
    sky=sky_input[cut_ch:(len(img_input[:,1,1])-cut_ch),:,:]

    ### begin to mask qso on sky
    # define final mask image
    mask_sum = np.zeros(sky[1,:,:].shape)

    ## looping over the layers
    for j in np.arange(0, len(sky[:,1,1]), 1):
        # iterated 10 times, rejecting points by > 3 sigma
        filtered_data = sigma_clip(sky[j,:,:], sigma=3, iters=10)
        sky_value = np.median(filtered_data.data[filtered_data.mask == False])
        # subtract the sky median value
        sky[j,:,:]-=sky_value
        img[j,:,:]-=sky_value
        # record the masked counts
        mask = np.zeros(sky[1,:,:].shape)
        mask[filtered_data.mask == True] = 1
        mask_sum+=mask

    # use mask_limit to avoid masking the good pixels (noise)
    mask_sum[mask_sum<mask_limit]=0
    mask_sum[mask_sum>=mask_limit]=1

    fits.writeto('mask.fits',mask_sum)
    fits.writeto('img.fits',img, ih)








def psf(img_file, dpx, dpy, scale):
    """ This function is used to subtract PSF for bright point source

    Parameters
    ----------
    img_file  :   default 'img.fits'
    dpx       :   cut PSF region in x
    dpy       :   cut PSF region in y
    scale     :   sclae method: 'peak', 'p', 'integrated', 'i'

    Eaxmple:
    psf('img.fits', dpx=8, dpy=10, scale='integrated')

    Returns
    -------
    img_psf.fits  : psf subtracted image
    """
    from astropy.io import fits
    import numpy as np
    img, ih = fits.getdata(img_file, header=True)

    ### PSF reduced
    img_median= np.median(img, axis=0)

    # find peak location
    loc = np.where(img_median > 0.9*np.amax(img_median))
    px= int(round(np.mean(loc[1])))
    py= int(round(np.mean(loc[0])))

    a = np.zeros(img[:,1,1].shape)
    psf = np.zeros(img.shape)

    ## scale to each channel
    if scale not in ['peak', 'p', 'integrated', 'i']:
        raise IOError('Need input the scale method! scale = peak / integrated')
    # peak value to scale
    if scale in ['peak', 'p']:
        peak_value = np.max(img_median)
        a=img[:,py,px]/peak_value

    # cut the region
    for j in np.arange(0, len(img[:,1,1]), 1):
        # integrated to scale
        if scale in ['integrated', 'i']:
            integrated_value = np.sum(img[j,(py-dpy):(py+dpy),(px-dpx):(px+dpx)])
            a[j]=integrated_value/np.sum(img_median[(py-dpy):(py+dpy),(px-dpx):(px+dpx)])
        psf[j,(py-dpy):(py+dpy),(px-dpx):(px+dpx)]=img_median[(py-dpy):(py+dpy),(px-dpx):(px+dpx)]*a[j]

    # reduce psf
    img-=psf

    fits.writeto('img_psf.fits',img, ih)
