#!/usr/bin/env python
from KStandard import attributesFromDict
from numpy import *
from KGeometry import *
from SHT import *
from FieldHarmonics import *
from pprint import pprint
import random
from SHT_Points import HarmPairs
#
# Here we do a coordinate transformation of the magnetic field harmonics.
#

def FieldHarmonicCoordTransform(harmonics,
                                orient=[Vector((0,0,0)),
                                        [[1,0,0],[0,1,0],[0,0,1]]],
                                #new_zvect=Vector((0,0,1)),
                                #new_xvect=Vector((1,0,0)),
                                #phi_rot=0.0,
                                BW=6,
                                show_points=False,
                                compute_pedantic=False):
 
    """
    orient has the form [RotationMatrix,TranslationVector] or
    [TranslationVector,RotationMatrix] depending on the order
    in the list the order of operations is switched.
    """
    fharms=[]
    for k,(A,B) in harmonics.iteritems():
        ns,ms=k.split("_")
        n=int(ns)
        m=int(ms)
        fharms.append(FieldHarmonic(n,m,A,B))

 
    #rotationTransformation=find_vector_zx_rotation_transform(new_zvect,new_xvect)
    #rot_matrix=find_vector_z_rotation_transform(new_zvect,phi_rot)
    measpoints=[Vector((1.0,theta,phi),vtype="SPHERICAL") for theta,phi in measurement_points(BW)]
    
    # Then we rotate and shift them appropriately for the new frame.
    newmeaspoints=[m.transform(orient) for m in measpoints]
    if hasattr(orient[0],'__len__'):
        # The first one is the rotation matrix
        rot_matrix=orient[0]    
    else:
        # The first one is the translation vector
        rot_matrix=orient[1]
    inv_rot_matrix=inv(rot_matrix)
    

    # Now we compute the contribution of the harmonics on the new data points.
    newdataBz=[]
    newdataBx=[]
    newdataBy=[]
    for m in newmeaspoints:
        tmpBx=0
        tmpBy=0
        tmpBz=0
        for fh in fharms:
            B=fh.B_val_point(m) # Just find the contribution of the harmonic to the field point.
            B=B.rotate(inv_rot_matrix)
            tmpBx=tmpBx+B.x()
            tmpBy=tmpBy+B.y()
            tmpBz=tmpBz+B.z()
        newdataBx.append(tmpBx)
        newdataBy.append(tmpBy)
        newdataBz.append(tmpBz)
    
    

    # Now we find the Spherical Harmonic Transform of the Bz data
    resultBz=SHT(BW,
                 realdata=newdataBz,
                 mode="COS_SIN",
                 harmtype="_{n,m}",
                 Rref=1.0)
    # But this misses out the super sectoral harmonics so we need more data
    resultBx=SHT(BW,
                 realdata=newdataBx,
                 mode="COS_SIN",
                 harmtype="_{n,m}",
                 Rref=1.0)
    resultBy=SHT(BW,
                 realdata=newdataBy,
                 mode="COS_SIN",
                 harmtype="_{n,m}",
                 Rref=1.0)
    # Combine the Bx and By data to figure out the super-sectoral harmonics.
    # See page 96 of the thesis.
    #pprint(resultBz)
    result={}
    for ((nx,mx),(Bxa,Bxb)),((ny,my),(Bya,Byb)),((nz,mz),(Bza,Bzb)) in zip(resultBx,resultBy,resultBz):
        if nx!=ny!=nz or mx!=my!=mz:
            raise Exception
        result[str(nx)+"_"+str(mx)]=(Bza,Bzb)
        if nx==mx:
            # We can compute a super sectoral here
            Bza_m_mp1=((1+delta(mx,0))/float(2.0*mx+1.0))*(Bxa-Byb)
            Bzb_m_mp1=((1+delta(mx,0))/float(2.0*mx+1.0))*(Bxb+Bya)
            result[str(mx)+"_"+str(mx+1)]=(Bza_m_mp1,Bzb_m_mp1)


   
    if show_points:
        from KGeomVisualization import Sphere,Colours,Visualization
        from KGeomVisualization import CoordinateAxes
        #pprint(Colours.keys())
        objs=[]
        objs.extend([Sphere(Vector(m.xyz()),
                            0.04,
                            color=Colours["orchid"]) for m in measpoints])
        objs.extend([Sphere(Vector(m.xyz()),
                            0.04,
                            color=Colours["AliceBlue"]) for m in newmeaspoints])
        coordsOrig=CoordinateAxes()
        coordsNew=CoordinateAxes(origin=new_origin,
                                 rot_matrix=rot_matrix)
                                 #new_zvect=new_zvect,
                                 #phi_rot=phi_rot)#new_xvect=new_xvect)
        objs.append(coordsOrig)
        objs.append(coordsNew)
        Visualization(objs)

    
    
    if compute_pedantic:
        # Now we compute the contribution of the harmonics on the old data points.
        olddataBxyz=[]
        for m in measpoints:
            tmpBx=0
            tmpBy=0
            tmpBz=0
            for fh in fharms:
                B=fh.B_val_point(m) # Just find the contribution of the harmonic to the field point.
                tmpBx=tmpBx+B.x()
                tmpBy=tmpBy+B.y()
                tmpBz=tmpBz+B.z()
            olddataBxyz.append((tmpBx,tmpBy,tmpBz))

        print "############"
        print "# measpoints"
        for m,Bxyz in zip(measpoints,olddataBxyz): 
            print "mp:",m.xyz(),"\nB:",Bxyz,"\n\n"


        print "################"
        print "# newmeaspoints"
        for m,Bx,By,Bz in zip(newmeaspoints,newdataBx,newdataBy,newdataBz): 
            print "mp:",m.xyz(),"\nB:",(Bx,By,Bz),"\n\n"


    # The result is a list of spherical harmonics in the new coordinate frame.
    return result
    

def _testRandom(DataBW=2,TotalBW=6,verbose=False):
    harmonics={} 
    if TotalBW<DataBW:
        raise Exception
    for n,m in HarmPairs(DataBW,supersectorals=True):
        Ba=random.uniform(-10,10)+100
        if m!=0:
            Bb=random.uniform(-10,10)+100
        else:
            Bb=0.0
        harmonics[str(n)+"_"+str(m)]=(Ba,Bb)


    new_origin=Vector((random.uniform(-10.0,10.0),
                       random.uniform(-10.0,10.0),
                       random.uniform(-10.0,10.0)))
                       
    alpha=random.uniform(0.0,2*pi)
    beta=random.uniform(0.0,2*pi)
    gamma=random.uniform(0.0,2*pi)

    rot_matrix=find_rotation_transform((alpha,beta,gamma))
    inv_rot_matrix=inv(rot_matrix)
    
    newharms=FieldHarmonicCoordTransform(harmonics,
                                         orient=[new_origin,rot_matrix],
                                         BW=TotalBW)
    mirror_origin=new_origin.negate()
    recovharms=FieldHarmonicCoordTransform(newharms,
                                           orient=[inv_rot_matrix,mirror_origin],
                                           BW=TotalBW)
    
    if verbose: print "harm:: original : new : recovered"
    testPassed=True
    maxDeviant=0
    for k,v in newharms.iteritems():
        if harmonics.has_key(k):
            origvals=harmonics[k]
        else:
            origvals=(0.0,0.0)
        if verbose: print k,"::",origvals,":",v,":",recovharms[k]
        for a,b in zip(origvals,recovharms[k]):
            diff = abs(a-b)
            if diff>maxDeviant:
                maxDeviant=diff
            Z=basicallyzero(diff,epsilon=1e-4)
            if verbose: print diff,Z
            if not Z:
                testPassed=False
    if verbose: print "TestPassed?:",testPassed
    return testPassed,maxDeviant
            


def _testBasic():
    harmonics={"0_0":(2,0),
               "0_1":(1,-1),
               "1_0":(4.0,0.0),
               "1_1":(-5.0,5.0),
               "1_2":(8.0,6.0),
               #"2_0":(2.0,0.0),
               #"2_1":(2.0,5.0),
               #"2_2":(4.0,6.6),
               #"2_3":(-9,-13)
               }
    BW=2
    new_origin=Vector((1,1,1))
    alpha=pi/3.2
    beta=pi/3.2 #pi/2.0
    gamma=pi/12.0

    rot_matrix=find_rotation_transform((alpha,beta,gamma))
    
    inv_rot_matrix=inv(rot_matrix)
    
    newharms=FieldHarmonicCoordTransform(harmonics,
                                         orient=[new_origin,rot_matrix],
                                         BW=BW)
    mirror_origin=new_origin.negate()
    recovharms=FieldHarmonicCoordTransform(newharms,
                                           orient=[inv_rot_matrix,mirror_origin],
                                           BW=BW)
    
    print "harm:: original : new : recovered"
    for k,v in newharms.iteritems():
        if harmonics.has_key(k):
            origvals=harmonics[k]
        else:
            origvals=(0.0,0.0)
        print k,"::",origvals,":",v,":",recovharms[k]


if __name__ == '__main__':
    for i in range(10):
        print _testRandom(DataBW=4,
                          TotalBW=8,
                          verbose=False)

        
