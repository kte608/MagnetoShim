#!/usr/bin/env python
from KStandard import attributesFromDict
import cPickle
from pprint import pprint
from matplotlib import pyplot
from math import *
import SHT_Points
from SHT import SHT
from KGeometry import Vector
from FieldHarmonics import  ScalarPot2BzField
#
# I want to numerically check the design that
# has be produced by the linear programming routines.
#
# This means:
# - Load the design from the cPickle file produced by
#   InkShim.py
# - Break the design into really little magnetization
#   patches
# - Calculate the scalar potential from the patches
#   at many locations. Specifically those which can
#   be used to perform a spherical harmonic transform.
# - Determine the Bz field (maybe Bx and By too) from
#   the scalar potential harmonics.
# - Make a graph comparing what the design was intended
#   to produce to what the simulations say it will produce.
# - Try shifting the and rotating the origin of the harmonics
#   and see what happens.
#

def ScalarPotential(magnetization,sourcepoint,fieldpoint):
    """
    Given an infinitesimal magnetization at a source point,
    what is the scalar potential at the field point?
    """
    v=fieldpoint-sourcepoint # the vector from sourcepoint to fieldpoint
    absv=v.magnitude()
    result=magnetization.dot(v)/(4*pi*absv*absv*absv)
    return result

def ShimAnalysis(fname,options):
    # First we load all the data
    fin=file(fname,"rb")
    [patchlocations,
     resultVector,
     ShimMatrix,
     DesiredHarmonics,
     ProducedHarmonicsFromDesign,
     numZ,
     numPhi,
     Z,
     R,
     magLimitPerPixel,
     totalMagnetization]=cPickle.load(fin)
    
    if options.verbose:
        inkShimMatrix=resultVector.reshape(numZ,numPhi)
        pyplot.matshow(inkShimMatrix)
        pyplot.show()
    
    print "Z:",Z," numZ:",numZ," R:",R," numPhi:",numPhi
    pprint(DesiredHarmonics)
    pprint(ProducedHarmonicsFromDesign)
    patchHeight=Z/float(numZ)
    patchWidth =R*(2*pi)/float(numPhi)
    print "patchHeight:",patchHeight
    print "patchWidth:",patchWidth
    # Now we need to break the design into very small bits for
    # analysis. Forget that for now. Let's just use the patches
    # as they are
    print "Computing Scalar Potential"
    harmR=R*options.rSphereFact
    ThetaPhis=SHT_Points.measurement_points(options.HarmBW)
    ScalarPotVals=[]
    for theta,phi in ThetaPhis:
        fieldpoint=Vector((harmR,theta,phi),vtype="SPHERICAL")
        scpotval=0
        # Probably should use KFarm here to split the computation
        # load across many cores
        for sourcepoint,pixelmag in zip(patchlocations,resultVector):
            if(pixelmag<-magLimitPerPixel*0.01 or pixelmag>magLimitPerPixel*1.01):
                print "pixelmag:",pixelmag
                print "magLimitPerPixel:",magLimitPerPixel
                raise Exception
            magnetization=Vector((0,0,1)).scale(pixelmag)
            scpotval=scpotval+ScalarPotential(magnetization,sourcepoint,fieldpoint)
        ScalarPotVals.append(scpotval)
    ScalarPotHarms=SHT(options.HarmBW,harmtype="_{n,m}",mode="COS_SIN",
                       realdata=ScalarPotVals,
                       Rref=harmR)
    # Now we need to convert this into harmonics of Bz.
    BzHarms={}
    for (n,m),(A,B) in ScalarPotHarms:
        (n,m),(A,B)=ScalarPot2BzField(n,m,A,B)
        if (n!=None):
            BzHarms[(n,m)]=(A,B)
    # print the result
    print "Here are the Harmonics of the Bz field"
    pprint(BzHarms)

    # and write it to a file
    fnameout=fname.split(".")[0]+"Analysis.cPickle"
    fout=file(fnameout,"wb+")
    cPickle.dump([BzHarms,fname],fout)
    fout.close()

if __name__ == '__main__': 
    # This is where we make the commandline utility
    import sys
    from optparse import OptionParser

    usage="usage: %prog DesignFile.cPickle [options]"
    parser=OptionParser(usage)
    parser.add_option('-o','--output',
                      action="store",type="string",dest="outfname",
                      default=None,
                      help="The output filename. By default we make a filename using the inputfilename.")
    parser.add_option('-v','--verbose',
                      action='store_true',dest='verbose',
                      default=False,help="Write extra stuff to stdout.")
    parser.add_option('-a','--analysis_patch_size_mm',
                      action="store",type="float",
                      dest="AnaylsisPatchSize",
                      default=1.0,
                      help="The maximum patch width or height (in mm) to use for analysis.")
    parser.add_option('--HarmBW',
                      action="store",type="int",dest="HarmBW",
                      default=12,
                      help="The harmonic bandwidth of the analysis.")
    parser.add_option('--r_sphere_factor',
                      action="store",type="float",dest="rSphereFact",
                      default=0.5,
                      help="We multiply the coil radius by this factor to determine the radius of the sphere for finding spherical harmonics.")
    (options,args)=parser.parse_args()
    if len(args) < 1:
        parser.print_help()
        parser.error("Incorrect number of arguments")
        
    import os.path
    DesignFname  =os.path.abspath(args[0])
    
    # Do the deed
    ShimAnalysis(DesignFname,options)
