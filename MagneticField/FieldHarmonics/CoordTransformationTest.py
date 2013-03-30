#!/usr/bin/env python
from KStandard import attributesFromDict
from numpy import *
from KGeometry import *
from SHT import *
from FieldHarmonics import *
from pprint import pprint
import random
from SHT_Points import HarmPairs
from KColours import ColourEnumerator
# So we can do +/- 20 mm about the origin.
# What I want is a z-shift of 20 mm and 
# a rotation of both 10 degrees in each
# direction and 54.74 degress in each direction
# too.
# Use a harmonic bandwidith of 6.

def PointGenerator(radius,
                   alpha,beta,gamma,
                   x,y,z,BW=6):
    new_origin=Vector((x,y,z))
    measpoints=[Vector((radius,
                        theta,
                        phi),vtype="SPHERICAL")
                for theta,phi in measurement_points(BW)]
    rot_matrix=find_rotation_transform((alpha,beta,gamma))
    orient=[rot_matrix,new_origin]
    # Then we rotate and shift them appropriately for the new frame.
    newmeaspoints=[m.transform(orient) for m in measpoints]
    
    return newmeaspoints



if __name__ == '__main__': 
    radius=8 # millimeters # Maximum is 10 
    AllPoints=[]
    measpoints=PointGenerator(radius,
                              0,0,0,
                              0,0,0)
    AllPoints.append(["OrigPoints",measpoints])
    ##################################################
    measpoints=PointGenerator(radius,
                              0,10*pi/180.0,0,
                              0,0,0)
    #AllPoints.append(["SmallRotA",measpoints])
    measpoints=PointGenerator(radius,
                               pi/2.0,10*pi/180.0,0,
                               0,0,0)
    #AllPoints.append(["SmallRotB",measpoints])

    #################################################
    measpoints=PointGenerator(radius,
                              0,54.74*pi/180.0,0,
                              0,0,0)
    AllPoints.append(["MagicRotA",measpoints])
    measpoints=PointGenerator(radius,
                               pi/2.0,54.74*pi/180.0,0,
                               0,0,0)
    AllPoints.append(["MagicRotB",measpoints])
    
    #################################################
    measpoints=PointGenerator(radius,
                              0,0,0,
                              0,0,radius)
    AllPoints.append(["TransZ",measpoints])
    measpoints=PointGenerator(radius,
                              0,10*pi/180.0,0,
                              0,0,radius)
    #AllPoints.append(["SmallRotATransZ",measpoints])
    measpoints=PointGenerator(radius,
                               pi/2.0,10*pi/180.0,0,
                               0,0,radius)
    #AllPoints.append(["SmallRotBTransZ",measpoints])
    
    measpoints=PointGenerator(radius,
                              0,54.74*pi/180.0,0,
                              0,0,radius)
    AllPoints.append(["MagicRotATransZ",measpoints])
    measpoints=PointGenerator(radius,
                               pi/2.0,54.74*pi/180.0,0,
                               0,0,radius)
    AllPoints.append(["MagicRotBTransZ",measpoints])


    from KGeomVisualization import Sphere,Colours,Visualization
    from KGeomVisualization import CoordinateAxes
    
    objs=[]
    fout=file("DataPoints.txt","w+")
    fout.write("# Field point coordinates in millimeters X,Y,Z\n")
    pointnum=0
    for i,(name,measpoints) in enumerate(AllPoints):
        #objs.extend([Sphere(Vector(m.xyz()),
        #                    0.4,
        #                    color=ColourEnumerator(i,Mode="1.0")) 
        #             for m in measpoints])
        fout.write("# "+name+"\n")
        for m in measpoints:
            objs.append(Sphere(Vector(m.xyz()),
                            0.4,
                            color=ColourEnumerator(i,Mode="1.0")))
            fout.write(str(pointnum)+" : "+repr(m.x())+" "+repr(m.y())+" "+repr(m.z())+"\n")
            pointnum=pointnum+1
    fout.close()
    Visualization(objs)
    
