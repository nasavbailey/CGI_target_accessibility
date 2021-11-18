import datetime
import numpy as np

from astropy.io import fits, ascii

import astropy.units as u
import astropy.coordinates

from pandas import DataFrame

from os import mkdir

# rot? returns rotation matrix for rotation by inp radians about axis xyz
# usage is np.dot(rot?(ang), vec) rotates x,y,z vector vec by ang about ?

def rotx(inp): return np.array([[1., 0., 0.], [0., np.cos(inp), -np.sin(inp)], [0., np.sin(inp), np.cos(inp)]])

def roty(inp): return np.array([[np.cos(inp), 0., np.sin(inp)], [0., 1., 0.], [-np.sin(inp), 0., np.cos(inp)]])

def rotz(inp): return np.array([[np.cos(inp), -np.sin(inp), 0.], [np.sin(inp), np.cos(inp), 0.], [0., 0., 1.]])

# nrm returns inp normalized to unity length

def nrm(inp): return inp / np.sqrt((inp**2).sum() + (not inp.any()))

def rottox(ax):
    """rottox returns rotation matrix that will rotate input ax to x-axis"""
    axn = nrm(ax)
    az = np.arctan2(axn[1], axn[0])
    alt = np.arcsin(axn[2])
    return np.dot(roty(alt), rotz(-az))

def rotfromx(ax):
    """rotfromx returns rotation matrix that will rotate x-axis to input ax"""
    axn = nrm(ax)
    az = np.arctan2(axn[1], axn[0])
    alt = np.arcsin(axn[2])
    return np.dot(rotz(az), roty(-alt))

def rotabout(ax, inp):
    """rotabout returns rotation matrix, right-handed about direction ax by inp radians"""
    axn = nrm(ax)
    az = np.arctan2(axn[1], axn[0])
    alt = np.arcsin(axn[2])
    return np.linalg.multi_dot([rotz(az), roty(-alt), rotx(inp), roty(alt), rotz(-az)])

def radectoxyzq(radec):
    return np.array([np.cos(radec[0])*np.cos(radec[1]), np.sin(radec[0])*np.cos(radec[1]), np.sin(radec[1])])

def xyzqtoradec(xyzq):
    xyzqn = nrm(xyzq)
    return np.array([np.arctan2(xyzqn[1], xyzn[0]), np.arcsin(xyzqn[2])])

# q = equatorial
# l = ecliptic
# b = observatory bccs

# rd = ra,dec [radians], always equatorial
# xyz = direction cosines [dimensionless]

# r = rotation matrix

rdsun = fits.getdata('sunradec2026.fits') # 2x365, ICRS equatorial position of sun per day of 2026
xyzqsun = np.array([radectoxyzq(rd) for rd in rdsun.T]).T # 3x365, position of sun on unit equatorial sphere (direction cosines) per day of 2026

nlp = astropy.coordinates.SkyCoord(0*u.degree, 90*u.degree, \
    equinox=astropy.time.Time(datetime.datetime(2026,1,1,0,0,0) + datetime.timedelta(days=365//2), scale='utc'), \
    frame=astropy.coordinates.GeocentricTrueEcliptic) # north ecliptic pole in middle of year (assume negligible motion)
rdnlp = np.array([nlp.icrs.ra.radian, nlp.icrs.dec.radian]) # ra, dec of north ecliptic pole
xyzqnlp = radectoxyzq(rdnlp) # equatorial direction cosines of north ecliptic pole


def yp_xyzq(xyzqinp, iday):
    """given (xyz)q and iday return (yaw, pitch)"""
    rsuntoz = np.dot(roty(-np.pi/2.), rottox(xyzqsun[:,iday])) # rotation matrix that would take xyzqsun to z-axis (BCCS); mutiply by this to align sun with solar array normal
    xyzsunznlp = np.dot(rsuntoz, xyzqnlp) # north ecliptic pole after alignment of sun to bccs z-axis (intermediate coordinate system)
    rznlp = rotz(-np.arctan2(xyzsunznlp[1], xyzsunznlp[0])) # rotation about BCCS z that puts xyznpsunz on x-axis
    rb = np.dot(rznlp, rsuntoz) # rotation matrix from equatorial xyz that puts sun on x-axis and ecliptic north pole on x-axis (BCCS)

    xyzbinp = np.dot(rb, xyzqinp) # input in bccs
    yaw = np.arctan2(xyzbinp[1], xyzbinp[0]) # rotation of observatory about z will put observatory x-axis at same azimuth as target
    pitch = -np.arcsin(xyzbinp[2]) # rotation of observatory about new y-axis will put obs x-axis at same altitude as target
    return np.array([yaw,pitch])

def yp_rd(rdinp, iday):
    """given (ra,dec) [radians] and iday return (yaw, pitch)"""
    return yp_xyzq(radectoxyzq(rdinp), iday)

def is_observable(rdinp, iday, max_dgr=34):
    """given (ra,dec) [radians] and iday return (yaw, pitch).
    Observable when abs(pitch) <=34 degrees.
    Returns True/False observability."""

    # require units on input ra and dec
    if type(rdinp) is astropy.coordinates.sky_coordinate.SkyCoord:
        rdinp = np.array([rdinp.icrs.ra.radian, rdinp.icrs.dec.radian]) * u.radian
    elif not all(isinstance(n, astropy.units.quantity.Quantity) for n in rdinp):
        raise Exception('rdinp must be either astropy sky coord type OR a tuple of with each entry having an astropy unit')

    tmp = yp_rd(rdinp,iday) * 180./np.pi
    return np.abs(tmp[1]) <= max_dgr



