#!/usr/bin/env python
from KStandard import attributesFromDict
from numpy import *
from KGeometry import *
from SH_Definitions import legendrePnm_theta as Pnmcos
from SH_Definitions import factorialQuotient,epsilonm
# You will need:
# sudo apt-get install coinor-csdp
# sudo apt-get install coinor-libcbc-dev
# sudo apt-get install coinor-libcgl-dev
# sudo apt-get install coinor-libclp-dev
# sudo apt-get install coinor-coinutils-dev
# sudo apt-get install coinor-libdylp-dev
# sudo apt-get install coinor-libflopc++-dev
# sudo apt-get install coinor-libipopt-dev
# sudo apt-get install coinor-libosi-dev
# sudo apt-get install coinor-libsymphony-dev
# sudo apt-get install coinor-libvol-dev
# To install easy_install use: sudo apt-get install python-pip
# To install PuLP use: sudo easy_install -U pulp
from pulp import *
from SHT_Points import *
from FieldHarmonics import *
from numpy.linalg import pinv
import cPickle
import os


class PatchShimmerDirect:
    """
    Given a bunch of locations where infinitesimal
    amounts of magnetization can be placed and some
    field points where a desired Bz field is indicated,
    this code attempts to determine the magnetizations.
    This is to be done directly from the magnetic
    scalar potential without any word about spherical
    harmonics.
    """
    def __init__(self,patchPositions,fieldPointsAndVals):
        attributesFromDict(locals())
        
    def solve(self,verbose=True):
        from ShimAnalysis import ScalarPotential
        prob=LpProblem("Shimming",LpMinimize)
        
        variables={} 
        patchlocations=[]

        # Make the variables
        patchlabels=[str(i) for i in range(len(self.patchPositions))]
        magnetizations=LpVariable.dicts("Mz_",patchlabels,0)
        
        # Set the objective function
        prob+=lpSum([magnetizations[i] for i in patchlabels]),"Minimize the Magnetization"

        def PatchContribution(fp,i):
            sp=self.patchPositions[int(i)]
            mag=Vector((0,0,1.0))
            return ScalarPotential(mag,sp,fp)

        # Set the field constraints.
        for j,(fp,val) in enumerate(self.fieldPointsAndVals):
            prob+=lpSum([PatchContribution(fp,i)*magnetizations[i] for i in patchlabels]),"Constrain field point"+str(i)
                
        # Write the problem to a file
        prob.writeLP("MyProblem.lp")
    
        print "Problem written to file ... Beginning Solver"
        prob.solve()
        print "Status:", LpStatus[prob.status]

        resultVector=array([v.varValue for v in prob.variables()])
        #from pprint import pprint
        #pprint(resultVector)
        totalmagnetization=value(prob.objective)
        return resultVector,totalmagnetization
        

class CylindricalShimmer:
    """
    Given some array of magnetization points, large static magnetic
    field, and a small set of shim harmonics to produce, this thing
    computes the required magnetizations at the points.
    """
        
    def __init__(self,
                 R=0.162,       #0.16 #0.165 too big #radius meters
                 Z=0.26,        #length meters
                 numZ=50,#64,    #T20 50   #number of segmentation
                 numPhi=36,#64,  #T20 36   #number of segmentation
                 HarmBW=4,      #shm spherical harmonic max
                 magLimitPerArea=0.3, #T20 magLimit =0.05 #magnetization max
                 fname="fout", #depends on the name of spherical harmonics
                 ShimFieldHarmonics={(2,2):(6,0),
                                     (2,0):(6,0)}):

        fname=fname+".cPickle"
        self.dz=Z/(numZ-1)
        self.dphi=2*pi/(numPhi) # No -1 because the topology is circular
        PixelArea=self.dz*R*self.dphi
        magLimitPerPixel=magLimitPerArea*PixelArea
        attributesFromDict(locals())

    def solve(self,verbose=True):
        """
        Here we actually try to solve the problem.
        """
        
        prob=LpProblem("Shimming",LpMinimize)
        
        dz=self.Z/(self.numZ)   # note that the patch position is at the patch center and not on the edge. Computational geometery...
        dphi=2*pi/(self.numPhi) # circular topology
        variablenames=[] # we create void lists
        patchlocations=[]

        

        ###########################################################
        # Generate the variables (the point magnetizations)
        ###########################################################
        for iZ in range(self.numZ):
            for iPhi in range(self.numPhi):
                name="Mz_"+str(iZ)+"_"+str(iPhi)
                variablenames.append(name)
                patchlocations.append(Vector((self.R,
                                              iPhi*dphi,
                                              -self.Z*0.5+iZ*dz+dz*0.5),
                                             vtype="CYLINDRICAL"))
                exec(name+"=LpVariable(name,0,"+repr(self.magLimitPerPixel)+",LpContinuous)")
        print "Variables Set"
        ############################################################
        # Set the objective function which represents 
        # the total amount of magnetization
        ############################################################
        objString=variablenames[0]
        for vn in variablenames[1:]:
            objString=objString+" + "+vn
        #print objString
        exec('prob += '+objString+', "Minimize the Magnetization"')

        print "Objective Function Set"

        ############################################################
        # Describe the problem with a matrix where the rows 
        # correspond to the contributions of the magnetization
        # patches to each harmonic and each columns corresponds
        # to a patch. 
        ############################################################
        ShimMatrix=[]
        #print HarmPairs(self.HarmBW)
        #raw_input()
        for [n,m] in HarmPairs(self.HarmBW):
            coeffG=-(epsilonm(m)/(4*pi))
            coeffG=coeffG*factorialQuotient(n-m+2,n+m)
            rowCoeffs=[]
            for pl in patchlocations:
                f=pl.rsphere()
                alpha=pl.theta()
                psi=pl.phi()
                coeff=coeffG*Pnmcos(n+2,m,alpha)*(1.0/f**(n+3))
                rowCoeffs.append((coeff,(f,alpha,psi)))
            # First a row for Ba
            row=[]    
            for coeff,(f,alpha,psi) in rowCoeffs:
                row.append(coeff*cos(m*psi))
            ShimMatrix.append(row)
            
            # Then a row for Bb
            row=[]
            for coeff,(f,alpha,psi) in rowCoeffs:
                row.append(coeff*sin(m*psi))
            ShimMatrix.append(row)
        print "Shim Matrix Online"
        
        ##############################################################
        # Now we constrain the problem such that the ShimMatrix is
        # true. That is we want the total field to add up correctly.
        # We have two constraints for every n,m pair because there
        # are two field terms (Ba and Bb).
        ##############################################################
        print "Generating constraints"
        for i,[n,m] in enumerate(HarmPairs(self.HarmBW)):
            if self.ShimFieldHarmonics.has_key((n,m)):
                va,vb=self.ShimFieldHarmonics[n,m]
            else:
                va=vb=0
            # First for Ba
            if not (n==0 and m==0):
                constraintString=repr(ShimMatrix[2*i+0][0])+"*"+variablenames[0]
                for vname,coeff in zip(variablenames[1:],ShimMatrix[2*i+0][1:]):
                    constraintString=constraintString+" + "+repr(coeff)+"*"+vname
                exec('prob += '+constraintString+'=='+repr(va)+' ,"Ba_'+str(n)+'_'+str(m)+'"')
            # Then for Bb
            if m!=0:
                constraintString=repr(ShimMatrix[2*i+1][0])+"*"+variablenames[0]
                for vname,coeff in zip(variablenames[1:],ShimMatrix[2*i+1][1:]):
                    constraintString=constraintString+" + "+repr(coeff)+"*"+vname
                exec('prob += '+constraintString+'=='+repr(vb)+' ,"Bb_'+str(n)+'_'+str(m)+'"')
            
        print "Constraints Constrained"
        
        # Write the problem to a file
        prob.writeLP("MyProblem.lp")
    
        print "Problem written to file ... Beginning Solver"
        prob.solve()

        resultVector=array([v.varValue for v in prob.variables()])

        Harmonics=dot(array(ShimMatrix),resultVector)

        #if verbose:
        #    # Each of the variables is printed with it's resolved optimum value
        #    for v in prob.variables():
        #        print v.name, "=", v.varValue

        # The status of the solution is printed to the screen
        ProducedHarmonics={} #assuming infinitesimal patches
        print "Status:", LpStatus[prob.status]

        for i,[n,m] in enumerate(HarmPairs(self.HarmBW)):
            j=2*i
            tmpA=Harmonics[j]
            print "Ba_"+str(n)+"_"+str(m),"=",tmpA
            
            j=j+1
            tmpB=Harmonics[j]
            print "Bb_"+str(n)+"_"+str(m)+"=",tmpB
            ProducedHarmonics[(n,m)]=(tmpA,tmpB)

        # The optimised objective function value is printed to the screen
        totalmagnetization=value(prob.objective)
        print "The total magnetization = ", totalmagnetization
        #write values on a file
        fout=file(self.fname,"wb+")
        cPickle.dump([patchlocations,
                      resultVector,
                      ShimMatrix,
                      self.ShimFieldHarmonics, #The desired harmonics
                      ProducedHarmonics,
                      self.numZ,
                      self.numPhi,
                      self.Z,
                      self.R,
                      self.magLimitPerPixel,
                      totalmagnetization],fout)
        fout.close()
        
        return self.fname

def OriginalInkShim():
    from optparse import OptionParser
    usage="usage: %prog [options] '(n,m):(Ba,Bb)' '(n,m):(Ba,Bb)' ..."
    parser=OptionParser(usage)

    parser.add_option('--radius',"-R",
                      action="store",type="float",dest="R",
                      default=0.162,
                      help="The cylindrical former's radius.")
    parser.add_option('--length',"-Z",
                      action="store",type="float",dest="Z",
                      default=0.26,
                      help="The cylindrical former's length.")
    parser.add_option('--numZ',
                      action="store",type="int",dest="numZ",
                      default=50,
                      help="The number of segments along the length.")
    parser.add_option('--numPhi',
                      action="store",type="int",dest="numPhi",
                      default=50,
                      help="The number of segments around the circumference.")
    parser.add_option('--HarmBW',
                      action="store",type="int",dest="HarmBW",
                      default=4,
                      help="The harmonic bandwidth.")
    parser.add_option('--magLimitPerArea',
                      action="store",type="float",dest="magLimitPerArea",
                      default=0.3,
                      help="Just like it sounds.")
    parser.add_option('--fname',"-o",
                      action="store",type="string",dest="fname",
                      default="out",
                      help="The outputfilename which will have '.cPickle' appended.")
    (options,args)=parser.parse_args()

    
    if len(args) < 1:
        parser.print_help()
        parser.error("Incorrect number of arguments. You need to specify at least one harmonic to produce!")
    FieldHarmonics={}
    for arg in args:
        nm,vals=arg.split(":")
        vals=vals.replace("(","").replace(")","")
        nm=nm.replace("(","").replace(")","")
        n,m=nm.split(",")
        Ba,Bb=vals.split(",")
        FieldHarmonics[(int(n),int(m))]=(float(Ba),float(Bb))
    
    C=CylindricalShimmer(R=options.R,Z=options.Z,
                         numZ=options.numZ,numPhi=options.numPhi,HarmBW=options.HarmBW,
                         magLimitPerArea=options.magLimitPerArea,fname=options.fname,
                         ShimFieldHarmonics=FieldHarmonics)
                         
    fname=C.solve()
    print fname

class DirectShimmer():
    def __init__(self,
                 R=0.162,       #0.16 #0.165 too big #radius meters
                 Z=0.26,        #length meters
                 numZ=50,#64,    #T20 50   #number of segmentation
                 numPhi=36,#64,  #T20 36   #number of segmentation
                 HarmBW=4,      #shm spherical harmonic max
                 magLimitPerArea=0.3, #T20 magLimit =0.05 #magnetization max
                 fname="fout", #depends on the name of spherical harmonics
                 ShimFieldHarmonics={(2,2):(6,0),
                                     (2,0):(6,0)}):

        self.ShimField=FieldHarmonicEnsemble(ShimFieldHarmonics)
        fname=fname+".cPickle"
        self.dz=Z/(numZ-1)
        self.dphi=2*pi/(numPhi) # No -1 because the topology is circular
        PixelArea=self.dz*R*self.dphi
        magLimitPerPixel=magLimitPerArea*PixelArea
        attributesFromDict(locals())
    def shim(self):
        fieldPointsAndVals=[]
        for t,p in measurement_points(self.HarmBW):
            fp=Vector((self.R*0.75,t,p),
                      vtype="CYLINDRICAL")
            fv=self.ShimField.B_val_point(fp)
            fieldPointsAndVals.append([fp,fv.z()])

        dz=self.Z/(self.numZ)   # note that the patch position is at the patch center and not on the edge. Computational geometery...
        dphi=2*pi/(self.numPhi) # circular topology
    
        patchlocations=[]
        for iZ in range(self.numZ):
            for iPhi in range(self.numPhi):
                patchlocations.append(Vector((self.R,
                                              iPhi*dphi,
                                              -self.Z*0.5+iZ*dz+dz*0.5),
                                             vtype="CYLINDRICAL"))
    
        ps=PatchShimmerDirect(patchlocations,fieldPointsAndVals)
        resultVector,totalmagnetization=ps.solve()
        
        fout=file(self.fname,"wb+")
        cPickle.dump([patchlocations,
                      resultVector,
                      [],#ShimMatrix
                      self.ShimFieldHarmonics, #The desired harmonics
                      [],#ProducedHarmonics,
                      self.numZ,
                      self.numPhi,
                      self.Z,
                      self.R,
                      self.magLimitPerPixel,
                      totalmagnetization],fout)
        fout.close()

if __name__ == '__main__':
    #OriginalInkShim()
    ds=DirectShimmer()
    ds.shim()
