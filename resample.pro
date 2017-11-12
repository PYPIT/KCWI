; the simplest example to understand resampling image

pro resample

; Create the display window
WINDOW, /FREE, XSIZE=20000, YSIZE=10000

; Set up the arrays of original points
yo = INDGEN(100)
xo = fix(sqrt(yo))

; Create and display the original image:
img = FLTARR(100,100)
xt1 = fix(sqrt(yo))+30
xt2 = fix(sqrt(yo))+60
img[yo,xo]=1
img[yo,xt1]=1
img[yo,xt2]=1

; X and Y coordinates to be fit
xi = fix(FLTARR(100) + median(xo))
yi = yo

; Run POLYWARP to obtain a Kx and Ky:
polywarp, xi, yi, xo, yo, 3, kx, ky

; Create a warped image based on Kx and Ky with POLY_2D:
wimg = poly_2d(img, kx, ky, 2,cubic=-0.5)

TVSCL, img, 0
TVSCL, wimg, 1

end