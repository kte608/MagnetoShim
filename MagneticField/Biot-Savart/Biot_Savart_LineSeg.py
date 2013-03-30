#!/usr/bin/env python
from RingArrayFieldViewUtils import fieldplot
from scipy.integrate.quadpack import quad as quadrature 
from scipy import *
from KGeometry import *
from copy import deepcopy
from constants import mu_0
import SHT
import cPickle
import BiotSavartLineSeg

print "I think you should use BiotSavartLineSeg instead... see _Biot_Savart_LineSeg"
raise Exception

def SimpleSeg():
    from CollectionViewer import collviewer_natural
    current=1.0
    L=1.0
    LineSegList=[LineSeg(Vector((0,0,0),type="CARTESIAN"),Vector((0,0,L),type="CARTESIAN"))]
    #collviewer_natural(linelist=[[l,(0,0,0)] for l in LineSegList])
    fieldfile="SingSegfield.txt"
    fieldfind(L,LineSegList,fieldfile)
    fieldplot(fieldfile)
    
def MultiSimpleSeg():
    from CollectionViewer import collviewer_natural
    from scipy import arange
    current=1.0
    numsegs=10
    L=1.0
    dL=L/float(numsegs)
    starts=arange(0,L,dL)
    LineSegList=[LineSeg(Vector((0,0,a),type="CARTESIAN"),Vector((0,0,a+dL),type="CARTESIAN")) for a in starts]
    #collviewer_natural(linelist=[[l,(0,0,0)] for l in LineSegList])
    fieldfile="MultiSingSegfield.txt"
    fieldfind(L,LineSegList,fieldfile)
    fieldplot(fieldfile)
    
    
    
if __name__ == "__main__":
    #KernelTest()
    #test_simplering()
    #Maxwell()
    #Helmholtz()
    Golay()
    #MultiSimpleSeg()
    #SimpleSeg()
