#!/usr/bin/env python
from KStandard import attributesFromDict
from numpy import *
from KGeometry import *
from SHT import *
from FieldHarmonics import *
from pprint import pprint
import random
from SHT_Points import HarmPairs
import os.path
import fileaccess

#
# This code has something to do with rotating and translating spherical harmonics like I would want to do with Dimitris.
#

def SHTFieldPoints(args):
    from optparse import OptionParser
    usage="usage: %prog SHTfieldpoints HarmBW [options]"
    parser=OptionParser(usage)
    parser.add_option('--VectorFormat',
                      action="store",type="string",dest="vectorformat",
                      default="CARTESIAN",
                      help="Valid vector formats are CARTESIAN, SPHERICAL, and CYLINDRICAL"),
    parser.add_option("--radius",action="store",type="float",dest="radius",default=1.0,
                      help="The desired measurement radius"),
    parser.add_option('-o','--output',
                      action="store",type="string",dest="outfname",
                      default="fieldpoints.txt",
                      help="The output filename.")
    parser.add_option('-v','--verbose',
                      action='store_true',dest='verbose',
                      default=False,help="Write extra stuff to stdout.")
    (options,args)=parser.parse_args(args)
    if len(args) < 1:
        parser.print_help()
        parser.error("Incorrect number of arguments")
            
    BW=int(args[0])
    # Do the deed
    measpoints=[Vector((options.radius,theta,phi),vtype="SPHERICAL") for theta,phi in measurement_points(BW)]

    fout=file(os.path.abspath(os.path.expanduser(options.outfname)),"w+")
    fout.write("# Field points for Spherical Harmonic Transform\n")
    fout.write("# Harmonic Bandwidth="+str(BW)+"\n")
    fout.write("# Radius="+repr(options.radius)+"\n")
    fout.write("# VectorFormat="+options.vectorformat+"\n")
    for i, m in enumerate(measpoints):
        fout.write(str(i)+" -> "+m.__str__(Mode=options.vectorformat)+"\n")
    fout.close()

def Evaluate(args):
    from optparse import OptionParser
    usage="usage: %prog Evaluate HarmMagnitudes FieldPoints [options]"
    parser=OptionParser(usage)
    parser.add_option('-o','--output',
                      action="store",type="string",dest="outfname",
                      default="fieldvalues.txt",
                      help="The output filename.")
    parser.add_option('-v','--verbose',
                      action='store_true',dest='verbose',
                      default=False,help="Write extra stuff to stdout.")
    (options,args)=parser.parse_args(args)
    if len(args) < 2:
        parser.print_help()
        parser.error("Incorrect number of arguments")

    ##########################################################
    # First we read in the Harms text file. This file is of the form:
    #
    # BW=5
    # (1,2)=(0.325,6.625)
    # (3,2)=(0.251,1.215)
    lines=fileacess.getlineslist(args[0])
    # The first line is supposed to be the BW
    #########################################################
    tmp=l.split("=")        
    if tmp[0].strip() != "BW":
        raise Exception
    BW=int(tmp[1])
    harmonics={}
    for n,m in HarmPairs(DataBW,supersectorals=True):
        harmonics[str(n)+"_"+str(m)]=(0.0,0.0)
    for l in lines[1:]:
        tmp=split("=")
        n,m=tmp[0].split(",")
        n=n.replace("(","").strip()
        m=m.replace(")","").strip()
        key=n+"_"+m
        Ba,Bb=tmp[1].split(",")
        Ba=float(Ba.replace("(","").strip())
        Bb=float(Bb.replace(")","").strip())
        if int(m)==0 and Bb!=0.0:
            print "Warning: ("+n+","+m+") is a zonal harmonic but Bb is not zero!"
            raise Exception
        if not harmonics.has_key(key):
            raise Exception
        else:
            harmonics[key]=(Ba,Bb)

    #############################################
    # Now we read the file with the field point
    # locations. The field point locations
    # are in the form:
    # 0 -> XYZ:(0.2,0.1,-0.25)
    # 1 -> RTP:(0.2,0.12,1.75)
    # 2 -> RZP:(0.1,0.75,3.0)
    lines=fileacess.getlineslist(args[1])
    fieldpoints=[]
    for l in lines:
        v=split("->")
        i=v[0].strip()
        tmp=v[1].split(":")
        a,b,c=tmp[1].replace("(","").replace(")").split(",")
        a=float(a)
        b=float(b)
        c=float(c)
        vtype={"XYZ":"CARTESIAN",
               "RTP":"SPHERICAL",
               "RZP"."CYLINDRICAL"}[tmp[0].strip()]
        fieldpoints.append([i,Vector((a,b,c),vtype=vtype)])

    # Do the deed
    fieldvalues=[]
    for i,v in fieldpoints:
        tmpBx=0
        tmpBy=0
        tmpBz=0
        for fh in harmonics:
            B=fh.B_val_point(m) # Just find the contribution of the harmonic to the field point.
            tmpBx=tmpBx+B.x()
            tmpBy=tmpBy+B.y()
            tmpBz=tmpBz+B.z()
        fieldvalues.append(Vector((tmpBx,tmpBy,tmpBz)))

    # Write the result
    fout=file(os.path.abspath(os.path.expanduser(options.outfname)),"w+")
    fout.write("# Field points for Spherical Harmonic Transform\n")
    fout.write("# Harmonic Bandwidth="+str(BW)+"\n")
    fout.write("# Radius="+repr(options.radius)+"\n")
    fout.write("# VectorFormat="+options.vectorformat+"\n")
    for i, m in enumerate(measpoints):
        fout.write(str(i)+" -> "+m.__str__(Mode=options.vectorformat)+"\n")
    fout.close()

def FindHarmonics(args):
    pass

def TransformHarmonics(args):
    pass

commandDict={"SHTfieldpoints":SHTFieldPoints, # Find the field points required for the SHT
             "Evaluate":Evaluate,   # Evaluate the field and/or potentials at given points from harmonic magnitudes
             "FindHarmonics":FindHarmonics, # Given some field points and the field values at those points we develop a list of spherical harmonics. May include the NSSHT.
             "TransformHarmonics":TransformHarmonics, # Given a harmonic representation and new orientation information (rotation and translation) we find a new set of harmonic magnitudes in the new coordinate frame
             
                 }
def Help():
    print "Invalid command type. Valid commands:"
    for k in commandDict.keys():
        print "\t",k
    sys.exit()

if __name__ == '__main__': 
    import sys
    

    if len(sys.argv)<2:
        Help()
        
    commandType=sys.argv[1]
    if commandDict.has_key(commandType):
        commandDict[commandType](sys.argv[2:])
    else:
        Help()
