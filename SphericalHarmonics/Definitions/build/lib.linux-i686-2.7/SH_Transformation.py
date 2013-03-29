#!/usr/bin/env python
from copy import deepcopy
from SH_Definitions import SphereHarmonic
import SHT_Points
from KGeometry import Vector
import KGeometry
from SHT import SHT
from math import pi,cos,sin
import numpy as np
def spherical_harmonic_translation_and_rotation(harmonics,
                                                new_zvect=Vector((0,0,1)), # The direction of the new z-axis.
                                                new_xvect=Vector((1,0,0)),
                                                new_origin=Vector((0,0,0)),    # The location of the new coordinate origin.
                                                phirot=0,             # We can numerically rotate about phi or do so in a closed form 
                                                                      # with the azimuthal rotation functions below. 
                                                BW=6,                 # Harmonic bandwidth
                                                harmtype='_{n,m}',    # Choose a spherical harmonic definition convention.
                                                verbose=True,
                                                show_points=False       # Graphical Representation of measurement points.
                                                ):
    """
    In this function we assume that we have a harmonic description around the z-vector (0,0,z) and we want a new description about
    a new vector. We also rotate the description azimuthally about this new vector by the amount phirot in radians.

    The principle is entirely numerical. We start by producing a single function of (x,y,z) from the given spherical harmonics.
    Then we define some measurement points about the new origin and oriented about the new direction. The function is then evaluated
    at these measurement points and -- provided we have chosen the points correctly -- the spherical harmonic transform of the computed
    values is used to produce a set of spherical harmonic coefficients valid for the new frame of reference.

    This problem can also be solved in closed form but I have not yet succeeded in doing so. See SphericalHarmonics.lyx for my attempts
    at the closed form solution.
    """
    from KGeometry import find_vector_zx_rotation_transform
    rotationTransformation=find_vector_zx_rotation_transform(new_zvect,new_xvect)

    # Setup something to evaluate the total contribution of the harmonics. 
    sharms=[]
    for (n,m),(ccoeff,scoeff) in harmonics:
        sharms.append(SphereHarmonic(n,m,ccoeff,scoeff,mode="COS_SIN",harmtype=harmtype))

    # First we make some measurement points for computing the spherical harmonic transform
    measpoints=[Vector((1.0,theta,phi),vtype="SPHERICAL") for theta,phi in SHT_Points.measurement_points(BW)]
    # Then we rotate and shift them appropriately for the new frame.
    newmeaspoints=[m.rotate(rotationTransformation)+new_origin for m in measpoints]
    
    if show_points:
        from KGeomVisualization import Sphere,Colours,Visualization
        from KGeomVisualization import CoordinateAxes
        objs=[]
        objs.extend([Sphere(Vector(m.xyz()),
                            0.04,
                            color=Colours["red"]) for m in measpoints])
        objs.extend([Sphere(Vector(m.xyz()),
                            0.04,
                            color=Colours["blue"]) for m in newmeaspoints])
        coordsOrig=CoordinateAxes(colors=[Colours["red"] for i in range(4)])
        rotmat=KGeometry.find_vector_zx_rotation_transform(new_zvect,
                                                           new_xvect)
        coordsNew=CoordinateAxes(origin=new_origin,
                                 rot_matrix=rotmat,
                                 colors=[Colours["blue"] for i in range(4)])
        objs.append(coordsOrig)
        objs.append(coordsNew)
        Visualization(objs)

    # Now we compute the contribution of the harmonics on the new data points.
    newdata=[]
    for m in newmeaspoints:
        tmp=0
        for sh in sharms:
            tmp=tmp+sh.val_point(m) # Just find the contribution of the harmonic to the field point.
        newdata.append(tmp)
    
    # Now we find the Spherical Harmonic Transform of the newdata
    result=SHT(BW,realdata=newdata,mode="COS_SIN",harmtype=harmtype,Rref=1.0)
    
    # The result is a list of spherical harmonics in the new coordinate frame.
    return result

def spherical_harmonic_rot_azimuth_kernel(m,alpha,Qa,Qb):
    cma=cos(m*alpha)
    sma=sin(m*alpha)
    Qap=Qa*cma-Qb*sma
    Qbp=Qb*cma+Qa*sma
    return (Qap,Qbp)
    
def spherical_harmonic_rot_azimuth(harmonics,alpha):
    """
    See the lyx document on spherical harmonics to see how this works.
    
    alpha is the azimuthal rotation in radians

    Essentially if Qa(n,m) and Qb(n,m) are the coefficients of T(n,m) and Tp(n,m)
    in the unrotated coordinate frame then Qap(n,m) and Qbp(n,m) are the coefficients
    in the rotated frame.

    Qap(n,m)=Qa(n,m) cos(m alpha)-Qb(n,m) sin(m alpha)
    Qbp(n,m)=Qb(n,m) cos(m alpha)+Qa(n,m) sin(m alpha)
    """
    oldharmonics=deepcopy(harmonics)
    newharmonics=[]

    for (n,m),(Qa,Qb) in oldharmonics:
        Qap,Qbp=spherical_harmonic_rot_azimuth_kernel(m,alpha,Qa,Qb)
        
        newharmonics.append(((n,m),(Qap,Qbp)))
    return newharmonics

def spherical_harmonic_rotation(harmonics,theta,phi):
    """
    rotate the harmonics by theta in declination and phi in azimuth
    """ 
    harmonics=spherical_harmonic_rot_azimuth(harmonics,phi)
    harmonics=spherical_harmonic_rot_declination(harmonics,theta)
    return harmonics

def _testspherical_harmonic_translation_and_rotation(BW=6):
    from pprint import pprint
    from KGeometry import Vector
    harmonics=[[(1,1),(0.0,-1.0)],[(5,1),(0.0,0.1)]]
    #pprint(harmonics)
    new_zvect=Vector((0,0,1)) # The direction of the new z-axis.
    new_xvect=Vector((1,0,0)) # The direction of the new x-axis.
    new_origin=Vector((0.0,0,0.5))    # The new origin
    forward=spherical_harmonic_translation_and_rotation(harmonics,
                                                        new_zvect=new_zvect,
                                                        new_xvect=new_xvect,
                                                        new_origin=new_origin,
                                                        BW=BW)
    #pprint(forward)
    new_zvect=Vector((0,0,1)) # The direction of the new z-axis.
    new_xvect=Vector((0,0,0)) # The direction of the new x-axis.
    new_origin=Vector((0.0,0,-0.5))    # The new origin
    backagain=spherical_harmonic_translation_and_rotation(forward,
                                                          new_zvect=new_zvect,
                                                          new_xvect=new_xvect,
                                                          new_origin=new_origin,
                                                          BW=BW)
    #pprint(backagain)
    
    

    # Let's plot the results as a bar graph
    # First we need to fill in the zero values in the
    # list called 'harmonics'
    indicies=SHT_Points.genindicies_top(BW)
    filledHarmonics=[]
    for n,m in indicies:
        tmp=[(n,m),(0,0)]
        for (N,M),(BA,BB) in harmonics:
            if N==n and M==m:
                tmp=[(N,M),(BA,BB)]
                break
        filledHarmonics.append(tmp)
    import matplotlib.pyplot as plt
    # Now we make some lists we can give to the 
    # bar graph routines.
    BaOriginals=[]
    BbOriginals=[]
    for (n,m),(Ba,Bb) in filledHarmonics:
        BaOriginals.append(Ba)
        BbOriginals.append(Bb)
    BaForwards=[]
    BbForwards=[]
    for (n,m),(Ba,Bb) in forward:
        BaForwards.append(Ba)
        BbForwards.append(Bb)
    BaBackAgain=[]
    BbBackAgain=[]
    for (n,m),(Ba,Bb) in backagain:
        BaBackAgain.append(Ba)
        BbBackAgain.append(Bb)
    #pprint(BbOriginals)
    #pprint(BbBackAgain)
    # Now we do the bar graph.
    N = len(filledHarmonics)

    ind = np.arange(N)  # the x locations for the groups
    width = 0.1       # the width of the bars

    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    rects1 = ax.bar(ind,         BaOriginals, width, color='r')#, yerr=menStd)
    rects2 = ax.bar(ind+width,   BbOriginals, width, color='b')
    rects3 = ax.bar(ind+width*2, BaForwards,  width, color='y')#, yerr=menStd)
    rects4 = ax.bar(ind+width*3, BbForwards,  width, color='g')
    rects5 = ax.bar(ind+width*4, BaBackAgain, width, color='r')#, yerr=menStd)
    rects6 = ax.bar(ind+width*5, BbBackAgain, width, color='b')

    # add some labels
    ax.set_ylabel('Harmonic')
    ax.set_title('Harmonic Order and Degree')
    ax.set_xticks(ind+width)
    labels=[]
    for n,m in indicies:
        labels.append(str(n)+","+str(m))
    ax.set_xticklabels( labels )

    ax.legend( (rects1[0], rects2[0],rects3[0],rects4[0],rects5[0],rects6[0]), ('Ba Original', 'Bb Original','Ba Forward', 'Bb Forward','Ba backagain','Bb backagain') )

    def autolabel(rects):
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            print height
            ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                    ha='center', va='bottom')

    #autolabel(rects1)
    #autolabel(rects2)

    plt.show()


if __name__ == "__main__":
    from pprint import pprint
    _testspherical_harmonic_translation_and_rotation()
