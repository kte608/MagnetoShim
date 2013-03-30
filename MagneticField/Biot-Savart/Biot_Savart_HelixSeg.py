#!/usr/bin/env python
#
# What we are doing is computing the equations given in
# SHID/Documents/MagneticField/BiotSavart_Helix_Seg.lyx
#

from RingArrayFieldViewUtils import fieldplot
from scipy.integrate.quadpack import quad as quadrature 
from scipy import *
from KGeometry import HelixSeg,Cylinder,Vector
from copy import deepcopy
from constants import mu_0
import SHT
import cPickle

def fieldfind(R,HelixSegList,filename):
    harmonicBW=12
    RSphereList=[R,0.95*R,0.5*R,0.25*R]
    Big_rvect,Big_zvect,Big_fieldBox=Boxfieldfind(R,HelixSegList)
    Small_rvect,Small_zvect,Small_fieldBox=Boxfieldfind(R,HelixSegList,BoxSize=R)
    theta_phi_points=SHT.measurement_points(harmonicBW)
    BFieldonSpheres=[Sphericalfieldfind(R,HelixSegList,Rsphere,theta_phi_points) for Rsphere in RSphereList] 

    BFieldSHTs=[]
    for BField in BFieldonSpheres:
        fields=transpose(BField)
        Bx=fields[0]
        By=fields[1]
        Bz=fields[2]
        Bs=fields[3]
        tmp=[SHT.SHT(harmonicBW,realdata=f,mode="COS_SIN") for f in fields]
        BFieldSHTs.append(tmp)

    q=["['Units_r','Units_z','Units_B'],Big_rvect,Big_zvect,Big_fieldBox(Bx,By,Bz,fieldsquared),Small_rvect,Small_zvect,Small_fieldBox(Bx,By,Bz,fieldsquared),harmonicBW,RSphereList,BFieldonSpheres,BFieldSHTs(one for each RSphere and (Bx,By,Bz,B2) for each one)",
       ['meters','meters','teslas'],Big_rvect,Big_zvect,Big_fieldBox,Small_rvect,Small_zvect,Small_fieldBox,harmonicBW,RSphereList,BFieldonSpheres,BFieldSHTs]
    f=file(filename,"w+")
    cPickle.dump(q,f)
    f.close()


def Boxfieldfind(R,HelixSegList,BoxSize="AUTO",dyz="AUTO"):
    # All dimensions should be in meters to work correctly
    # The field results are in Teslas


    if BoxSize=="AUTO":
        plotradial=3*R #meters (a box this size will be plotted in the middle)
        plotaxial=3*R  #meters
    else:
        plotradial=BoxSize
        plotaxial=BoxSize
    if dyz=="AUTO":
        dyz=0.5*R/3.0       #meters (plot interval)


    rvect=arange(-plotradial,plotradial+dyz,dyz)
    zvect=arange(-plotaxial,plotaxial+dyz,dyz)

    result=[]
    y=0
    for x in rvect:
        print "x:",x
        tmp=[MagFieldHelixSegArray(HelixSegList,Vector((x,y,z),type="CARTESIAN")) for z in zvect]
        result.append(tmp)
    return rvect,zvect,result

def Sphericalfieldfind(R,HelixSegList,Rsphere,theta_phi_points):
    """ Find all components of the magnetic field on points on a sphere that are condusive to computing the
    spherical harmonic transform """
    BField=[MagFieldHelixSegArray(HelixSegList,Vector((Rsphere,theta,phi),type="SPHERICAL")) for theta,phi in theta_phi_points]
    return BField

def MagFieldHelixSegArray(seglist,fieldpoint,returnasVectors=False,returnErrors=False):
    """
    

    """
    
    Btmp=Vector((0,0,0),type="CARTESIAN")
    dBtmp=Vector((0,0,0),type="CARTESIAN")
    for s in seglist:
        B,dB=_Biot_Savart_SingleHelixSeg(s,fieldpoint,returnasVectors=True)
        Btmp=Btmp.add(B)
        dBtmp=dBtmp.add(dB)
    if returnasVectors:
        if returnErrors:
            return ((deepcopy(Btmp),deepcopy(dBtmp)))
        else:
            return deepcopy(Btmp)
    else:
        if returnErrors:
            return ((Btmp.x,Btmp.y,Btmp.z,(Btmp.magnitude)*(Btmp.magnitude)),
                    (dBtmp.x,dBtmp.y,dBtmp.z,(dBtmp.magnitude)*(dBtmp.magnitude)))
        else:
            return (Btmp.x,Btmp.y,Btmp.z,(Btmp.magnitude)*(Btmp.magnitude))
    

def _Biot_Savart_SingleHelixSeg(helixseg,fieldpoint_vector,returnasVectors=False):
    """
    Given a helix segment on the surface of a cylinder
    and a vector indicating the field point, this function
    computes Bx,By,Bz,and B2. Errors are also given. The results
    are Teslas per amp of current so You can multiply the result
    by the actual current to get the field.
    """
    R=helixseg.cylind.radius
    const=mu_0/float(4*pi*R)
 
    M=helixseg.M/float(R)
    phiS=helixseg.cylind_startvect.phi # this is guaranteed to be zero
    if phiS != 0.0:
        print "Someone has changed the definition of HelixSeg"
        raise Exception
    phiE=helixseg.cylind_endvect.phi
    x=fieldpoint_vector.x/float(R)
    y=fieldpoint_vector.y/float(R)
    fieldpoint_cylind=helixseg.transform_to_helixcoords(fieldpoint_vector) # We need some transformation to flip into coordinates referenced to the cylinder.
    
    z=(helixseg.cylind_startvect.z-fieldpoint_cylind.z)/float(R)

    #print "x,y,z,M,phiS,phiE:",x," ",y," ",z," ",M," ",phiS," ",phiE

    Bx,dBx=quadrature(Afunc,phiS,phiE,args=(x,y,z,M,phiS))
    By,dBy=quadrature(Bfunc,phiS,phiE,args=(x,y,z,M,phiS))
    Bz,dBz=quadrature(Cfunc,phiS,phiE,args=(x,y,z,M,phiS))
    
    Bvect=Vector((const*Bx,const*By,const*Bz),type="CARTESIAN")
    dBvect=Vector((const*dBx,const*dBy,const*dBz),type="CARTESIAN")

    Bvect=helixseg.transform_from_helixcoords(Bvect)
    dBvect=helixseg.transform_from_helixcoords(dBvect)

    if returnasVectors==True:
        return (Bvect,dBvect)
    else:
        return ((Bvect.x,Bvect.y,Bvect.z,Bvect.magnitude*Bvect.magnitude),
                (dBvect.x,dBvect.y,dBvect.z,(dBvect.magnitude)*(dBvect.magnitude)))

def denominator(theta,x,y,z,M,deltatheta):
    Q2=x*x+y*y+z*z+1
    t=deltatheta
    return 1.0*(Q2-2*(x+y)*cos(theta)+2*z*M*t+M*M*t*t)

def Afunc(theta,x,y,z,M,thetaS):
    t=theta-thetaS
    #print "Theta S is:",thetaS
    #print "Theta is:",theta
    #print "(x,y,z,M):",x," ",y," ",z," ",M
    return 1.0*(-M*y-(z+M*t)*cos(theta)+M*sin(theta))/denominator(theta,x,y,z,M,t)


def Bfunc(theta,x,y,z,M,thetaS):
    t=theta-thetaS
    return 1.0*(M*x-(z+M*t)*sin(theta)-M*cos(theta))/denominator(theta,x,y,z,M,t)

def Cfunc(theta,x,y,z,M,thetaS):
    t=theta-thetaS
    return 1.0*(1-x*cos(theta)-y*sin(theta))/denominator(theta,x,y,z,M,t)

def KernelTest():
    print Cfunc(0,0,0,0,0,0)
    print quadrature(Cfunc,0,pi,args=(0,0,0,0,0))

def test_simplering():
    """
    Define a ring of helix segments and see what the field looks like.
    """
    R=0.10 #meters
    cylind=Cylinder(Vector((0,0,0),type="CARTESIAN"),Vector((0,0,1),type="CARTESIAN"),R)
    SegList=[HelixSeg(Vector((R,0,0),type="CYLINDRICAL"),
                      Vector((R,pi,0),type="CYLINDRICAL"),cylind),
             HelixSeg(Vector((R,pi,0),type="CYLINDRICAL"),
                      Vector((R,2*pi,0),type="CYLINDRICAL"),cylind)]
    Fieldpoint=Vector((0,0,0),type="CARTESIAN")
    print MagFieldHelixSegArray(SegList,Fieldpoint)


def RingMaker(zpos,direction,numsegs,cylind):
    """
    Returns a list of helix segments
    """
    dphi=direction*2*pi/numsegs
    result=[]
    phi=0
    R=cylind.radius
    if direction != 1 and direction != -1:
        raise Exception

    for i in range(numsegs):
        result.append(HelixSeg(Vector((R,phi,zpos),type="CYLINDRICAL"),
                               Vector((R,phi+dphi,zpos),type="CYLINDRICAL"),
                               cylind))
        phi=phi+dphi
    return result

def Maxwell():
    current=1.0
    Rring=1.0
    # [turns,z,radius]
    cylind=Cylinder(Vector((0,0,-10),type="CARTESIAN"),
                    Vector((0,0,10),type="CARTESIAN"),Rring)
    HelixSegList=RingMaker(sqrt(3)/2.0,1,16,cylind)
    HelixSegList.extend(RingMaker(-sqrt(3)/2.0,-1,16,cylind))
    fieldfile="Maxwellfield.txt"
    fieldfind(Rring,HelixSegList,fieldfile)
    fieldplot(fieldfile)

def Helmholtz():
    current=1.0    #amps
    Rring=1.0
    RingList=[[current,Rring/2.0,Rring],
              [current,-Rring/2.0,Rring]]
    fieldfile="Helmfield.txt"
    fieldfind(Rring,RingList,fieldfile)
    fieldplot(fieldfile)

if __name__ == "__main__":
    #KernelTest()
    #test_simplering()
    Maxwell()
