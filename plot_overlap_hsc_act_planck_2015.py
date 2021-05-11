import numpy as np
import matplotlib.pyplot as plt

# for polygons
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

#import moby2
import pyfits
import healpy as hp

# to use flipper
import liteMap
import astLib.astWCS
import fftTools

###############################################################################

# takes an angle in deg,
# and returns it in [cen-180, cen+180],
# with cen in degrees
def angle(x, cen=110.):
   x *= np.pi/180.
   cen *= np.pi/180.
   c = np.cos(x-cen)
   s = np.sin(x-cen)
   t = np.arctan2(s, c) + cen
   t *= 180./np.pi
   return t

# center we adopt for ra and dec
# in degrees, such that
# ra, dec in [center-180, center+180]
center = 110.

###############################################################################

# compute area of polygon in sq. deg.
# corners = [[ra1,dec1], [ra2,dec2], ...] in deg.
# assumes flat sky
def polygonArea(corners):
   n = len(corners) # of corners
   area = 0.0
   for i in range(n):
      j = (i + 1) % n
      area += corners[i][0] * corners[j][1]
      area -= corners[j][0] * corners[i][1]
   area = abs(area) / 2.0
   return area

###############################################################################
###############################################################################
# HSC 1st year fields

# XMM
ra = np.array([29.5, 38., 41.5, 29.5])
dec = np.array([-1.5, -1.5, -7., -7.])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
xmm = angle(radec, cen=center)

# VVDS
ra = np.array([330., 344., 344., 330.])
dec = np.array([3.5, 3.5, -1.5, -1.5])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
vvds = angle(radec, cen=center)

# GAMA09
ra = np.array([128., 142.5, 142.5, 128])
dec = np.array([3., 3., -2., -2.])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
gama09 = angle(radec, cen=center)

# GAMA15
ra = np.array([213., 223., 223., 213.])
dec = np.array([2., 2., -2.5, -2.5])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
gama15 = angle(radec, cen=center)

# HECTOMAP
ra = np.array([243., 248., 248., 243.])
dec = np.array([45., 45., 42., 42.])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
hectomap = angle(radec, cen=center)

# WIDE12
# split into _a and _b, because crosses 180 deg
ra = np.array([175.5, 183., 183., 175.5])
dec = np.array([2., 2., -2.5, -2.5])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
wide12 = angle(radec, cen=center)
#
ra = np.array([175.5, 179.99, 179.99, 175.5])
dec = np.array([2., 2., -2.5, -2.5])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
wide12_a = angle(radec, cen=center)
#
ra = np.array([183, 180.01, 180.01, 183])
dec = np.array([2., 2., -2.5, -2.5])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
wide12_b = angle(radec, cen=center)


###############################################################################
# HSC full fields
# here the names are chosen randomly (by me...)

# XMM FULL
ra = np.array([-29.5, 39.5, 39.5, 29., 29., -29.5])
dec = np.array([6.5, 6.5, -6.5, -6.5, -1.5, -1.5])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
xmmfull = angle(radec, cen=center)

# GAMA FULL
# split into _a and _b, because crosses 180 deg
ra = np.array([128., 226., 226., 128.])
dec = np.array([5., 5., -2.5, -2.5])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
gamafull = angle(radec, cen=center)
#
ra = np.array([128., 179.99, 179.99, 128.])
dec = np.array([5., 5., -2.5, -2.5])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
gamafull_a = angle(radec, cen=center)
#
ra = np.array([180.01, 226., 226., 180.01])
dec = np.array([5., 5., -2.5, -2.5])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
gamafull_b = angle(radec, cen=center)


# HECTOMAP FULL
ra = np.array([200., 250., 250., 200.])
dec = np.array([44.5, 44.5, 42.5, 42.5])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
hectomapfull = angle(radec, cen=center)


###############################################################################
###############################################################################
# ACTPol

# D56
ra = np.array([-9., 41., 41., -9.])
dec = np.array([5., 5., -8., -8.])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
d56 = angle(radec, cen=center)

# BOSSNorth
# split into _a and _b, because crosses 180 deg
ra = np.array([147., 229.5, 229.5, 147.])
dec = np.array([19.5, 19.5, -4., -4.])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
bossnorth = angle(radec, cen=center)
#
ra = np.array([147., 179.99, 179.99, 147.])
dec = np.array([19.5, 19.5, -4., -4.])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
bossnorth_a = angle(radec, cen=center)
#
ra = np.array([180.01, 229.5, 229.5, 180.01])
dec = np.array([19.5, 19.5, -4., -4.])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
bossnorth_b = angle(radec, cen=center)


###############################################################################
###############################################################################
# AdvACT

ra = np.array([center-179.999, center+179.999, center+179.999, center-179.999])
ra = angle(ra, cen=center)
dec = np.array([28., 28., -75., -75.])
radec = np.array([(ra[i], dec[i]) for i in range(len(ra))])
advact = radec #= angle(radec, cen=center)

###############################################################################
###############################################################################
# CMASS

cmasssouth = angle(np.genfromtxt("./input/boss/cmass_south_dr12.txt"), cen=center)
cmassnorth = angle(np.genfromtxt("./input/boss/cmass_north_dr12.txt"), cen=center)
cmass = np.concatenate((cmassnorth, cmasssouth))


###############################################################################
###############################################################################
# plot overlap
'''
fig = plt.figure(0, figsize=(12,8))
ax=fig.add_subplot(111)
#
ax.axhline(0., c='k', ls='--')
ax.axvline(0., c='k', ls='--')

# AdvACT
patches = []
#
polygon = Polygon(advact, closed=True)
patches.append(polygon)
#
p = PatchCollection(patches, facecolor='deepskyblue', edgecolor='', alpha=0.1)
ax.add_collection(p)

# CMASS
plt.plot(cmass[:,0], cmass[:,1], 'y+', alpha=0.02)

# ACTPol
patches = []
# d56
polygon = Polygon(d56, closed=True)
patches.append(polygon)
# bossnorth
polygon = Polygon(bossnorth_a, closed=True)
patches.append(polygon)
polygon = Polygon(bossnorth_b, closed=True)
patches.append(polygon)
#
p = PatchCollection(patches, facecolor='c', edgecolor='', alpha=0.5)
ax.add_collection(p)

# HSC full
patches = []
# hectomap full
polygon = Polygon(xmmfull, closed=True)
patches.append(polygon)
#
# gama full
polygon = Polygon(gamafull_a, closed=True)
patches.append(polygon)
polygon = Polygon(gamafull_b, closed=True)
patches.append(polygon)
#
# hectomap full
polygon = Polygon(hectomapfull, closed=True)
patches.append(polygon)
#
p = PatchCollection(patches, facecolor='m', edgecolor='', alpha=0.4)
ax.add_collection(p)

# HSC first year
patches = []
# xmm
polygon = Polygon(xmm, closed=True)
patches.append(polygon)
# vvds
polygon = Polygon(vvds, closed=True)
patches.append(polygon)
# gama09
polygon = Polygon(gama09, closed=True)
patches.append(polygon)
# gama15
polygon = Polygon(gama15, closed=True)
patches.append(polygon)
# hectomap
polygon = Polygon(hectomap, closed=True)
patches.append(polygon)
# wide12
polygon = Polygon(wide12_a, closed=True)
patches.append(polygon)
polygon = Polygon(wide12_b, closed=True)
patches.append(polygon)
#
p = PatchCollection(patches, facecolor='m', edgecolor='', alpha=1.)
ax.add_collection(p)

ax.plot(np.mean(xmm[:,0]), np.mean(xmm[:,1]), alpha=0.)  # annoying: needed to make polygons happy
ax.set_xlim((center-180., center+180.))
ax.set_ylim((-20., 70.))
#ax.set_aspect('equal', adjustable='box')
ax.invert_xaxis()
ax.grid()
ax.set_xlabel('RA [deg]')
ax.set_ylabel('Dec [deg]')
#
#fig.savefig("./figures/overlap.png", bbox_inches='tight', dpi=400)

plt.show()
'''
###############################################################################
# compute approximate areas


# HSC 1st year
print ""
print "HSC 1st year"
print "XMM:", polygonArea(xmm), "sq.deg."
print "VVDS:", polygonArea(vvds), "sq.deg."
print "GAMA15:", polygonArea(gama15), "sq.deg."
print "WIDE12:", polygonArea(wide12), "sq.deg."
print "GAMA09:", polygonArea(gama09), "sq.deg."
print "HECTOMAP:", polygonArea(hectomap), "sq.deg."
print "total:", polygonArea(xmm)+polygonArea(vvds)+polygonArea(gama15)+polygonArea(wide12)+polygonArea(gama09)+polygonArea(hectomap), "sq.deg."

# HSC full
print ""
print "HSC full"
print "XMM:", polygonArea(xmmfull), "sq.deg."
print "GAMA:", polygonArea(gamafull), "sq.deg."
print "HECTOMAP:", polygonArea(hectomapfull), "sq.deg."
print "total:", polygonArea(xmmfull)+polygonArea(gamafull)+polygonArea(hectomapfull), "sq.deg."

# ACTPol
print ""
print "ACTPol"
print "D56:", polygonArea(d56), "sq.deg."
print "BOSS North:", polygonArea(bossnorth), "sq.deg."
print "total:", polygonArea(d56)+polygonArea(bossnorth), "sq.deg."








