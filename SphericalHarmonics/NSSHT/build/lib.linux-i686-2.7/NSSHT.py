#!/usr/bin/python
# Non-Standard Spherical Harmonic Transform using the Pseudo-Inverse
# The idea here is to construct a spherical harmonic transform on a sausage
# shape instead of a sphere.
import SHT
import SH_Definitions 
from pprint import pprint
#import CollectionViewer

from scipy.linalg import pinv as p_inv
from scipy import dot,pi,sin,cos
from KStandard import is_odd
from KGeometry import Vector

#class point(CollectionViewer.DataPoint):
#    def __init__(self,coords,type="NONE",color=(0,1,0)):
#        self._pointvals=SH_Definitions.point(coords,type=type)
#        self.x=self._pointvals.x
#        self.y=self._pointvals.y
#        self.z=self._pointvals.z
#        self.rsphere=self._pointvals.rsphere
#        self.theta=self._pointvals.theta
#        self.phi=self._pointvals.phi
#        self.rcylind=self._pointvals.rcylind
#        CollectionViewer.DataPoint.__init__(self,self.x,self.y,self.z,color=color)

def matrixbuilder(harmoniclist,pointlist):
    """
    N is a vector of harmonic amplitudes T,Tprime,T,Tprime,...
    M is a matrix relating those amplitudes to the field value
    at the points described in the pointlist
    """
    result=[]
    numharmonics=len(harmoniclist)
    for point in pointlist:
        tmprow=[]
        for i,harm in enumerate(harmoniclist):
            #print "i,numharm,end:",i,numharmonics,numharmonics-i-1
            tmprow.append(harm.val_point(point))            
        result.append(tmprow)
    return result

def SHTmeasurement_points(bw,r):#,color=(0,1,0)):
    """Return the SHT.measurement_points in terms of the point construct"""
    meas_points=SHT.measurement_points(bw)
    result=[Vector((r,theta,phi),type="SPHERICAL") for theta,phi in meas_points]#,color=color)
    return result

def UnityWeight_COSSIN_Harmonics(indicies,
                                 harmtype='_{n}^{m}',
                                 RSphere=1.0):
    """
    result forms a chain T,Tprime,T,Tprime,...
    """
    result=[]
    for n,m in indicies:
        result.append(SH_Definitions.SphereHarmonic(n,m,1.0,0.0,
                                                    Rref=RSphere,
                                                    harmtype=harmtype,mode="COS_SIN"))
        result.append(SH_Definitions.SphereHarmonic(n,m,0.0,1.0,
                                                    Rref=RSphere,
                                                    harmtype=harmtype,mode="COS_SIN"))
    return result

def NSSHT(bw,points,
          rdata,
          mode="COS_SIN",
          RSphere=1.0, #This is the sphere that the points would
                       #would be on if they were not constrained
                       #to a cylinder
          harmtype='_{n}^{m}',
          returnIntermediates=False):
    """ given a bandwidth, and some real data collected at a number of data points this should return the
    Spherical Harmonic Transform using the Pseudo-Inverse"""
    indicies=SHT.genindicies_top(bw)
    harms=UnityWeight_COSSIN_Harmonics(indicies,
                                       harmtype=harmtype,
                                       RSphere=RSphere)
    M=matrixbuilder(harms,points)
    Minv=p_inv(M)

    tmp=dot(Minv,rdata)
    # Since the result is in a form 
    # T,Tprime,T,Tprime we reorganize it
    # to be like (T,Tprime),(T,Tprime),...
    Q=[]
    for i in range(len(tmp)):
        if is_odd(i):
            Q.append((tmp[i-1],tmp[i]))

    cos_sinresult=zip(indicies,Q)
    if mode=="COS_SIN":
        if returnIntermediates:
            return cos_sinresult,tmp,M,Minv
        else:
            return cos_sinresult
    else:
        print "I have not thought of this yet, maybe later"
        raise Exception

def getrdata(meas_points,harms):
    rdata=[]
    for p in meas_points:
        tmp=0
        for h in harms:
            tmp=tmp+h.val_point(p)
        rdata.append(tmp)
    return rdata


def Cylinder_measpoints(bw,Rcylind,Rsphere):#,color=(0,1,0)
    """Slide point that are outside the cylinder down their radial component
    until they are on the surface of the cylinder. The cylinder has infinite
    axial extent"""
    spherepoints=SHTmeasurement_points(bw,Rsphere)#,color=color
    result=[]
    for p in spherepoints:
        r=p.rsphere()
        theta=p.theta()
        delta=r-(Rcylind/sin(theta))
        
        if delta > 0:
            newr=r-delta
            result.append(Vector((newr,theta,p.phi()),type="SPHERICAL"))#,color=color)
        else:
            result.append(p)
    return result
        


def testlevel3():
    """Now it is time for an altered set of measurement points on the surface of a cylinder rather than the surface of a sphere"""
    import KGeomVisualization as KGV
    n=1
    m=0
    cosmag=1
    sinmag=0
    bw=4
    r=20.0
    mode="COS_SIN"
    colorpoints=(1,0,0)
    colorpointsinner=(0,1,0)
    harms=[SH_Definitions.SphereHarmonic(n,m,cosmag,sinmag,mode=mode)]#,Rref=r)]
    meas_points=SHTmeasurement_points(bw,r)#,color=(1,0,0))
    meas_points_inner=Cylinder_measpoints(bw,r*0.15,r)#,color=(0,1,0))
    Q=[KGV.Sphere(p,r/20,color=colorpoints) for p in meas_points]#(0,1,meas_points
    Q.extend([KGV.Sphere(p,r/20,color=colorpointsinner) for p in meas_points_inner])#measpoints_inner)
    V=KGV.Visualization(Q)
    rdata_inner=getrdata(meas_points_inner,harms)
    result=NSSHT(bw,meas_points_inner,rdata_inner,Rref=r)
    pprint(result)

def testlevel2():
    """Using a sphere that is just a little smaller (for the measurement points) determine the harmonic transform on a larger sphere"""
    import KGeomVisualization as KGV
    n=1
    m=0
    cosmag=1
    sinmag=0
    bw=4
    r=20.0
    mode="COS_SIN"
    colorpoints=(1,0,0)
    colorpointsinner=(0,1,0)
    harms=[SH_Definitions.SphereHarmonic(n,m,cosmag,sinmag,mode=mode,Rref=r)]
    meas_points=SHTmeasurement_points(bw,r)#,color=(1,0,0))
    meas_points_inner=SHTmeasurement_points(bw,r*0.75)#,color=(0,1,0))
    Q=[CollectionViewer.DisplayPoint(p.x,p.y,p.z,color=colorpoints) for p in meas_points]#(0,1,meas_points
    Q.extend([CollectionViewer.DisplayPoint(p.x,p.y,p.z,color=colorpointsinner) for p in meas_points_inner])#measpoints_inner)
    #Q=meas_points
    #Q.extend(meas_points_inner)
    CollectionViewer.collectionviewer(SphereList=Q)
    rdata_inner=getrdata(meas_points_inner,harms)
    result=NSSHT(bw,meas_points_inner,rdata_inner,Rref=r)
    pprint(result)
    
    
        
def testlevel1():
    """Using exactly the same data and points as the SHT can we get the same
    result?"""
    #import KGeomVisualization as KGV
    n=1
    m=0
    cosmag=1
    sinmag=0
    bw=4
    r=1.0
    mode="COS_SIN"
    harms=[SH_Definitions.SphereHarmonic(n,m,cosmag,sinmag,mode=mode)]
    meas_points=SHTmeasurement_points(bw,r)#SHT.measurement_points(bw)
    rdata=getrdata(meas_points,harms)
        
    resultSHT=SHT.SHT(bw,realdata=rdata,mode=mode,Rref=r)
    resultNSSHT=NSSHT(bw,meas_points,rdata)

    
    #pprint(S)
    print "rdata:",len(rdata)
    #resultNSSHT=dot(Minv,rdata)
    print "NSSHT:",len(resultNSSHT)
    pprint(resultNSSHT)
    #pprint(M)
    pprint(resultSHT)


if __name__ == "__main__":
    """ I need to understand the meaning of negative harmonics!
    and even and odd!"""
    testlevel1()
    #testlevel2()
    #testlevel3()
    
