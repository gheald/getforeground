#!/usr/bin/env python

import h5py
import healpy
from numpy import *
from pylab import hist,show
import argparse
from astropy import units as u
from astropy.coordinates import SkyCoord

def main(args):
	# This script requires the Galactic foreground map from:
	# https://wwwmpa.mpa-garching.mpg.de/ift/faraday/2014/index.html
	f = h5py.File('foreground.hdf5','r')

	rm = f['maps']['phi'].value
	rmunc = f['uncertainty']['phi'].value

	if args.coords:
		c = SkyCoord(ra=args.lon*u.degree,dec=args.lat*u.degree,frame='icrs')
		lon = c.galactic.l.degree
		lat = c.galactic.b.degree
		print 'Interpreting input latitude and longitude as equatorial (RA,DEC)'
		print 'Converted to Galactic lon,lat = (%f,%f) degrees'%(lon,lat)
	else:
		print 'Interpreting input latitude and longitude as being specified in Galactic coordinates'
		lon = args.lon
		lat = args.lat

	pix = healpy.pixelfunc.ang2pix(128,(90.-lat)*pi/180.,lon*pi/180.)

	print 'RM +/ uncertainty is %f +/- %f rad/m2'%(rm[pix],rmunc[pix])

	nxy = asarray(args.nsamples.split(','),dtype=float)
	rmvals = []
	rmuncvals = []
	for i in arange(nxy[0])-floor(nxy[0]/2.):
		for j in arange(nxy[1])-floor(nxy[1]/2.):
			lat = args.lat + j*args.dang
			lon = args.lon + i*args.dang
			pix = healpy.pixelfunc.ang2pix(128,(90.-lat)*pi/180.,lon*pi/180.)
			rmvals.append(rm[pix])
			rmuncvals.append(rmunc[pix])
	rmvals = array(rmvals)
	rmuncvals = array(rmuncvals)
	medval = median(rmvals)
	rmsval = sqrt(sum(rmuncvals**2)/len(rmuncvals))
	print 'Median RM in the region is %f +/- %f rad/m2'%(medval,rmsval)
	print 'Size of region is %f deg x %f deg (lat,lon)'%((nxy[1]-1.)*args.dang,(nxy[0]-1.)*args.dang)

	if args.showmap:
		healpy.mollview(rm,title='Oppermann RM sky',min=rm[pix]-args.mult*rmunc[pix],max=rm[pix]+args.mult*rmunc[pix])
		healpy.graticule()
		show()
		hist(rmvals)
		show()

ap = argparse.ArgumentParser()
ap.add_argument('lat',help='Galactic latitude (degrees)',type=float)
ap.add_argument('lon',help='Galactic longitude (degrees)',type=float)
ap.add_argument('--coords','-c',help='Convert coordinates from RA/DEC to Galactic l/b? If so, provide DEC as lat, and RA and lon. [default False] NOTE: if declination is negative, precede it with "--" e.g. "./getforeground.py -c -- -34. 128."',default=False,action='store_true')
ap.add_argument('--showmap','-s',help='Show map and histogram? [default False]',action='store_true',default=False)
ap.add_argument('--mult','-m',help='Multiplier to set colour range [default 3]',default=3.,type=float)
ap.add_argument('--nsamples','-n',help='Number of samples to use for developing a wide view and region statistics (x,y) [default 1,1 i.e. no broader view/region]',default='1,1')
ap.add_argument('--dang','-d',help='Spacing in angular size to use for the samples (degrees) [default 0.4581, which is correct for NSIDE=128 in HealPIX maps]',default=0.4581,type=float)
args = ap.parse_args()
main(args)

