def find_center(image): 
    '''
    Find the center of QSOs for each target,
    we will return the center of the QSO; 
    and later, we will use the center to
    calculate the offset.  (Z. Cai)

    input: datacube 
    '''
    # construct a 2D broadband image using median stacking along the wavelength
    sci_med=np.median(image, axis=0)

    # interpolate to a finer grid:  
    x_ori= np.linspace(0, len(sci_med[:,0])-1, len(sci_med[:,0])) # 0--23 (24 pixels)
    y_ori= np.linspace(0, len(sci_med[0,:])-1, len(sci_med[0,:])) # 0--69 (70 pixels)
    x_new= np.linspace(0, len(sci_med[:,0])-1, 3.*len(sci_med[:,0])-2)  # new grid, 3x pixel number
    y_new= np.linspace(0, len(sci_med[0,:])-1, 3.*len(sci_med[0,:])-2)  # new grid, 3x pixel number

    f = interpolate.interp2d(y_ori, x_ori, sci_med,kind='cubic')
    sci_interp = f(y_new,x_new)
    
    # find center using pixels > 0.93* the peak value: 
    indices = np.where(sci_med > 0.93*np.amax(sci_med))
    indices_interp = np.where(sci_interp > 0.93*np.amax(sci_interp))

    # center using the original pixel grid      
    x_mean0= np.mean(indices[0])  
    y_mean0= np.mean(indices[1])

    # center using the new finer pixel grid (1/3 x original pixel scale) 
    x_mean_interp= np.mean(indices_interp[0])*1/3.
    y_mean_interp= np.mean(indices_interp[1])*1/3.

    # return center determined using both the original and new pixel grid.  
    center= [x_mean0, y_mean0, x_mean_interp,y_mean_interp]
    return center   


img1, ih1 = fits.getdata("kb171021_00077_icuber.fits", header=True)
center1=find_center(img1)

