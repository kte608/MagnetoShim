#!/usr/bin/env python
from math import pi,sqrt,atan2,pow
from KStandard import inclusive_range
#from SH_Definitions import *
from KGeometry import Vector

def declination_measurement_points(bandwidth):
    """This function returns a list of declination angles in
    radians where measurements should be made if the
    transform has a certain bandwidth"""
    thetas=[pi*(2*j+1)/(4*bandwidth) for j in inclusive_range(0,2*bandwidth-1)]
    return thetas

def azimuthal_measurement_points(bandwidth):
    """This function returns a list of azimuthal angles in
    radians where measurements should be made if the
    transform has a certain bandwidth"""
    phis=[2*pi*k/(2*bandwidth) for k in inclusive_range(0,2*bandwidth-1)]
    return phis

def azimuthal_delta_phi(bandwidth):
    return 2*pi/(2*bandwidth)

def declination_measurement_points_deg(bandwidth):
    start=declination_measurement_points(bandwidth)
    result=[]
    for i in start:
        result.append(i*180/pi)
    return result

def azimuthal_measurement_points_deg(bandwidth):
    start=azimuthal_measurement_points(bandwidth)
    result=[]
    for i in start:
        result.append(i*180/pi)
    return result

def ExperimentOrder2SHTOrder(points,data):
    """
    When I perform experiments I add an extra point
    at the top and bottom of the sphere and I alternate
    directions of travel in the azimuthal direction. 
    
    This function takes a list of points and a list
    of data and re-orders them to be of the form
    presented by measurement_points (you must explicitly
    exclude the the top and bottom cap points before
    you pass the data to this function)
    """
    BW=BandwidthFromNumPoints(len(points))
    #print "BW:",BW
    numaz=len(azimuthal_measurement_points(BW))
    numdec=len(declination_measurement_points(BW))
    reordered_points=[]
    reordered_data=[]
    #pnum=0
    
    for d in range(numdec):
        #print d
        for a in range(numaz):                
            if d%2==0:
                if a==0:
                    index=d*numaz+a
                    #print "not reordering:",index," d:",d," a:",a
                else:
                    index=(d+1)*(numaz)-a
                    #print "reordering:",index," d:",d," a:",a
                
                reordered_points.append(points[index])
                reordered_data.append(data[index])
            else:   
                index=None
                if a==0:
                    index=(d+1)*numaz-1
                    #print "reordering:",index," d:",d," a:",a
                if index==None:
                    index=d*numaz+a-1
                    #print "not reordering:",index," d:",d," a:",a
                
                reordered_points.append(points[index])
                reordered_data.append(data[index])
                

    #reordered_points=points
    #reordered_data=data
    return BW,reordered_points,reordered_data,numaz,numdec
    

def pprint_measurement_points(bandwidth):
    """Prints the azimuthal and declination angles in a nice format"""
    dcp=declination_measurement_points(bandwidth)
    azp=azimuthal_measurement_points(bandwidth)
    if(len(dcp)!=len(azp)):
        print "This is all wrong"

    print "declination\tazimuthal"
    print "______________________"

    for d,a in zip(dcp,azp):
        print d*180/pi,
        print "\t\t",
        print a*180/pi

    

def measurement_points(bandwidth):
    """This function returns a list of ordered pairs
    (theta,phi) of the declination and azimuthal angle of
    points that should be measured if the transform has a
    certain bandwidth"""
    result=[]
    phis=azimuthal_measurement_points(bandwidth)
    for theta in declination_measurement_points(bandwidth):
        for phi in phis:
            result.append((theta,phi))

    return result

def BandwidthFromNumPoints(numpoints):
    BWnumpoints={0:0,
                 4:1,
                 16:2,
                 36:3,
                 64:4,
                 100:5,
                 144:6,
                 196:7,
                 256:8}
    if BWnumpoints.has_key(numpoints):
        result=BWnumpoints[numpoints]
    else:
        result=None
    return result
                 


def measurement_points_deg(bandwidth):
    """This function does the same thing as measurement_points but returns
    the values in degrees"""
    temp=measurement_points(bandwidth)
    result=[]
    for pair in temp:
        result.append((pair[0]*180/pi,pair[1]*180/pi))

    return result
            

def genindicies(bandwidth):
    """This function generates a list of ordered pairs (n,m) that
    describe the order of the harmonics coming out of SHT"""
    result=genindicies_top(bandwidth)
    result.extend(genindicies_bottom(bandwidth))
    return result


def genindicies_bottom(bandwidth):
    result=[]
    for m in range(-(bandwidth-1),0):
        for n in range(-m,bandwidth):# This is correct and agrees with page 8 of s2kit_fx.pdf
            result.append((n,m))
    return result

def genindicies_top(bandwidth):
    result=[]
    for m in range(bandwidth):
        for n in range(m,bandwidth):
            result.append((n,m))
    return result


def genindicies_topplus(bandwidth):
    """
    Just like genindicies_top except this one also includes m-1,m
    """
    result=[]
    for m in range(bandwidth):
        for n in range(m-1,bandwidth):
            result.append((n,m))
    return result


def HarmIndex(n,m,BW):
    """
    returns the index into an array generated by genindicies
    that stores the harmonic n,m
    """
    if n>=BW or abs(m)>= BW or n<0 or m>n or m<0:
        return None
    am=abs(m)
    def beforeposfinder(absm,BW):
        return absm*(2*BW-absm+1)/2

    if m<0:
        posbefore=beforeposfinder(BW-1,BW)+1
        #print "PosBefore: ",posbefore
        negbefore=(BW-am)*(BW-1-am)/2
 
        #print "NegBefore: ",negbefore
        before=posbefore+negbefore
    else:
        before=beforeposfinder(am,BW)
    index=before+(n-am)
    #print "Before: ",before
    #print "Index: ",index
    return index


def HarmOrderWalker(nmax,BW):
    """ given a pre-established dataset this will give a list
    of indicies that will walk through the data set in the following order:

    0,0
    1,0
    1,1
    2,0
    2,1
    2,2
    etc..

    """
    result=[]
    for n in range(nmax+1):
        for m in range(n+1):
            #print n," ",m
            result.append(HarmIndex(n,m,BW))

    return result

def HarmPairs(BW,
              supersectorals=False #include terms where m=n+1
              ):
    """ return a list of harmonic orders and degrees
        in order for a certain BW """
    result=[]
    if supersectorals:
        inc=2
    else:
        inc=1
    for n in range(BW):
        for m in range(n+inc):
            result.append((n,m))
    return result


if __name__ == "__main__":
    print genindicies(6)
    #print measurement_points(6)
    print len(measurement_points(6))
    print HarmOrderWalker(5,21)
