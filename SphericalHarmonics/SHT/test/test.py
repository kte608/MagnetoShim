#!/usr/bin/env python
import SHT
from pprint import pprint
from KGeometry import Vector
from numpy import *


def multitester():
    """
     This function generates some data
     and then when the SHT is taken it should
     have the same harmonic components as were
     origininally generated
    """
    bw=8;
    #harmtype='_{n,m}'
    harmtype='_{n}^{m}'
    Rsphere=0.08
    meas_points=[]
    for theta,phi in SHT.measurement_points(bw):
        meas_points.append(Vector([Rsphere,theta,phi],type="SPHERICAL"))
    
    mode = "COS_SIN"
    harms=[((0,0),(1.0,0.0)),
           ((1,0),(1.0,0.0)),
           ((2,0),(1.0,0.0)),
           ((3,0),(1.0,0.0)),
           ((4,0),(2.0,0.0)),
           ((1,1),(1.0,1.0)),
           ((2,1),(1.0,1.0)),
           ((2,2),(0.0,1.0)),
           ((7,6),(0.0,20))]
    #harms=[((1,1),(1.0,0.0))]
           
    d,rdata=SHT.harmonizer(meas_points,harms,mode=mode,pointType="VECTOR",harmtype=harmtype)
  


    SHT.harmplotterMAT(bw,harms,Rsphere,mode=mode,harmtype=harmtype,type="3D")
    print "I am about to SHT"
    l=SHT.SHT(bw,realdata=rdata,mode=mode,Rref=Rsphere,harmtype=harmtype)
    print "I SHTed"
    pprint(l)

    return



def harmplotterGNU(bw,harms,mode="NONE"):
    """GNUPlot some harmonics"""
    import Gnuplot, Gnuplot.funcutils    
    g=Gnuplot.Gnuplot()
    g.title('Input Data')
    y=array(SHT.azimuthal_measurement_points(bw))
    x=array(SHT.declination_measurement_points(bw))
    g('set parametric')
    g('set xrange [0:pi]')
    g('set yrange [0:2*pi]')
    g('set data style lines')
    g('set hidden3d')
  #  g('set pm3d')
    g('set contour base')
    g.xlabel('theta (declination)')
    g.ylabel('phi (azimuth)')
    class funcboy:
        def __init__(self,harms,mode="NONE"):
            self.harms=harms
            self.mode=mode
        def val(self,theta,phi):
            d,rdata=SHT.harmonizer([(theta,phi)],self.harms,mode=self.mode)
            result=float(rdata[0])
            return result

      
    f=funcboy(harms,mode=mode)
    g.splot(Gnuplot.funcutils.compute_GridData(x,y,f.val))#GridData(m,x,y))
    raw_input('Please press return to continue...\n')

def harmonizer_tester(n,m,harms,bw,Rsphere,mode="NONE"):
    """
    This function shows the data produced by the harmonizer function
    at a number of data points given a certain set of harmonic components

    see 'anothertester'
    """

    meas_points=[]#SHT.measurement_points(bw)
    for theta,phi in SHT.measurement_points(bw):
        meas_points.append(Vector([Rsphere,theta,phi],type="SPHERICAL"))
    d,rdata=SHT.harmonizer(meas_points,harms,mode=mode,pointType="VECTOR");

    harmplotterMAT(bw,harms,mode=mode)
 
    l=SHT.SHT(bw,realdata=rdata,mode=mode,Rref=Rsphere)
    pprint(l)
    return

def tester():
    bw=4;
    num=len(SHT.measurement_points(bw))
    rdata=[1 for i in range(num)];
    idata=[0 for i in range(num)];
    pprint(bw)
    pprint(rdata)
    pprint(idata)
    #pprint(len(SHT.measurement_points(bw)))
    l=SHT.SHT(bw,rdata,idata)
    pprint(l)
    return

def anothertester():
    """Just cycle the values and see what you get!"""
    n=4
    m=1
    bw=5
    Rsphere=0.08
    mode = "COS_SIN"#"MAG_PHASE"#
    if mode == "COS_SIN":
        cosmag=0
        sinmag=1
        harms=[(n,m,cosmag,sinmag)]
    if mode == "MAG_PHASE":
        mag=1
        phase=0.501#pi/2.0
        harms=[(n,m,mag,phase)]

    
    harmonizer_tester(n,m,harms,bw,Rsphere,mode=mode);
    return


if __name__=="__main__":
    #tester()
    multitester()
    #anothertester()
