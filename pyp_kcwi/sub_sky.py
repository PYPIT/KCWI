# This function is used to subtract sky median
#
# INPUT:     - sky_cube
#            - img_cube
#            - mask_limit   mask limit (e.g., 10 means the pixel in the same position of each layer is masked 10 times throughout the datacube)
#            - cut_ch       cut layer number (the first and last few layers are bad)
#
# INPUT EXAMPLE:
#   median('kb171021_00083_icuber.fits', 'kb171021_00082_icuber.fits', 8, cut_ch=300)
#
# MODIFICATION HISTORY:
#   Written by: Qiong Li
#   Modified by Zheng Cai
#   2017-11-01 Initial version
#   2017-11-02 Second version

def median(sky_file, img_file, mask_limit, cut_ch=300):
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    sky_input, sh = fits.getdata(sky_file, header=True)
    img_input, ih = fits.getdata(img_file, header=True)

    img=img_input[cut_ch:(len(img_input[:,1,1])-cut_ch),:,:]
    sky=sky_input[cut_ch:(len(img_input[:,1,1])-cut_ch),:,:]

    ### begin to mask qso on sky
    # define final mask image
    mask_sum = np.zeros(sky[1,:,:].shape)

    ## sky[w,x,y]; for every wavelength: len(sky[:,1,1])
    for j in np.arange(0, len(sky[:,1,1]), 1):

        # sky image for each channel
        s=sky[j,:,:]

        # define mask image for each channel
        mask = np.zeros(sky[1,:,:].shape)

    # begin 10 loops
        for i in np.arange(0, 10, 1):

            # select: s>0(remove egde) & mask<1(loop new sky)
            need=(mask<1)&(s!=0)

            # mask out of 'mean + [-3, 3] sigma' as 1, sky image: mask == 0
            max_3 = np.median(s[need])+3*np.std(s[need])
            min_3 = np.median(s[need])-3*np.std(s[need])
            mask[(s>max_3) | (s<min_3) ]=1

        sky_value = np.median(s[mask == 0])

        sky[j,:,:]-=sky_value
        img[j,:,:]-=sky_value
        mask_sum+=mask


    ################ check and select limit:
    # transfer to 1D
    #mask_sum_1D=mask_sum[mask_sum>(mask_sum-1)]
    #hist=plt.hist(mask_sum_1D,np.arange(-10, 3000,1),histtype='step')
    #plt.xlabel('counts of mask')
    #plt.ylabel('counts of pixel')
    #plt.text(1100,1300,'e.g. 1400 pixels are masked 1000 times(channels)')
    #plt.show()
    ################
    mask_sum[mask_sum<mask_limit]=0     
    mask_sum[mask_sum>=mask_limit]=1   

    ################ check mask image:
    #face = plt.imshow(mask_sum, aspect='equal', interpolation='None',origin='lower',vmin=0,vmax=1.0,cmap = mpl.cm.gray)
    #cbar = plt.colorbar(face)
    #plt.text(55,95,'mask')
    #plt.text(55,-2,'sky')
    #plt.show()
    ################
    fits.writeto('mask.fits',mask_sum)
    fits.writeto('img.fits',img, ih)







# This function is used to subtract PSF for bright point source
#
# INPUT:     - img_file      default 'img.fits'
#            - dpx           cut PSF region in x
#            - dpy           cut PSF region in y
#
# INPUT EXAMPLE:
#   psf('img.fits', dpx=8, dpy=10)
#
# MODIFICATION HISTORY:
#   Written by: Qiong Li 
#   Modified by Zheng Cai
#   2017-11-01 Initial version

def psf(img_file, dpx, dpy):
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    img, ih = fits.getdata(img_file, header=True)

    ### PSF reduced
    #img_median = np.zeros(img[1,:,:].shape)
    #for x in np.arange(0, len(img[1,:,1]), 1):
    #    for y in np.arange(0, len(img[1,1,:]), 1):
    #        img_median[x,y] = np.median(img[:,x,y])
    img_median= np.median(img, axis=0)
    
    # find peak location
    #peak = np.max(img_median)
    #loc = np.where(img_median == peak)
    #px=loc[1][0]  # add[0]: Error only integer scalar arrays can be converted to a scalar index
    #py=loc[0][0]
    
    loc = np.where(img_median > 0.95*np.amax(img_median))    
    px= int(round(np.mean(indices[0])))
    py= int(round(np.mean(indices[1])))
    

    #################### test peak
    #face = plt.imshow(img_median, aspect='equal', interpolation='None',origin='lower',cmap = mpl.cm.gray)
    #cbar = plt.colorbar(face)
    #plt.scatter(px, py,color='r', marker='+', s=60,linewidths=1.5)
    #plt.show()
    ####################

    a = np.zeros(img[:,1,1].shape)
    psf = np.zeros(img.shape)

    # scale to each channel
    peak_value = np.max(img_median)
    a=img[:,py,px]/peak_value


    # cut the region
    for j in np.arange(0, len(img[:,1,1]), 1):
        psf[j,(py-dpy):(py+dpy),(px-dpx):(px+dpx)]=img_median[(py-dpy):(py+dpy),(px-dpx):(px+dpx)]*a[j]

    # reduce psf
    img-=psf

    fits.writeto('img_psf.fits',img, ih)
