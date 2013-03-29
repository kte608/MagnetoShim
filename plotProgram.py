#!/usr/bin/env python 
import cPickle
from matplotlib import pyplot
from KGeometry import *
from numpy import * 
from LDpsA3A4 import *

#To see T(n,m) take file Tnm.cPickle

def plotfunction(fname="out.cPickle",paper="A4"): #fname depends on the name of input file
    """
    it plots InkShimMatrix and launches the function that writes this matrix in a postscript.
    """
    
    fin=file(fname,"rb")
    [patchlocations,
     resultVector,
     ShimMatrix,
     ShimFieldHarmonics, #The desired harmonics
     ProducedHarmonics,
     numZ,
     numPhi,
     Z,
     R,
     magLimitPerPixel,
     totalmagnetization]=cPickle.load(fin)
    #print patchlocations
    #Magnetization Amounts Matrix
    inkShimMatrix=resultVector.reshape(numZ,numPhi)
    #pyplot.matshow(inkShimMatrix)
    #pyplot.show()
    
    ### GrayScale definition ###   
    GrayScale=1-inkShimMatrix/magLimitPerPixel
    pyplot.matshow(GrayScale)
    pyplot.show()
    ### Paper format: A3 or A4 ###
    if paper=="prompt":
        paper=raw_input("which paper format do you need: A3 or A4 ?\n")
    ### call the function that write the postscript     ###
    ### and allows to print it in the good paper format ###
    if paper=="A3":
        writepsA3("InkShimMatrixT22A3.ps",Z,R,numZ,numPhi,GrayScale)
    elif paper=="A4":
        writepsA4("InkShimMatrixT22A4.ps",Z,R,numZ,numPhi,GrayScale)
    else:
        print("unknown paper type\n")
        paper=raw_input("which paper format do you need: A3 or A4 ?\n")
        
    fin.close()

def plottest(paper="A3"):
    numZ=2
    numPhi=2
    Z=0.29
    R=0.07
    dx=0.05
    dy=0.02
    GrayScale=1-array([[0,0.1],[0,0.2]])/0.2
    if paper=="A4":
        writepsA4("printps.ps",Z,R,numZ,numPhi,GrayScale)
    elif paper=="A3":
        writepsA3("printps.ps",Z,R,numZ,numPhi,GrayScale)
    else:
        print("unknown paper type\n")

if __name__=="__main__":
    import sys
    plotfunction(fname=sys.argv[1])
    #plottest()    
